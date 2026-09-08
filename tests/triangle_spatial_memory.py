import itertools
import json
import random
import sys
import unittest
from collections import Counter, defaultdict
from fractions import Fraction as F
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'tools'))
from triangle_charge_response import ChargeMarginal, coefficients
from triangle_elastic_scattering import ElasticExperiment
from triangle_reference import AxialReference, subgroup
from triangle_relational_memory import conditional_marginals
from triangle_spatial_memory import SpatialObserver, Moments, conditioned_coefficients, conditioned_variances, trial


class SpatialMemoryTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.bank = json.loads((ROOT/'data/triangle-feedback.json').read_text())['bank']
        cls.elastic = json.loads((ROOT/'data/triangle-elastic-scattering.json').read_text())['pair_census']['elastic_tables']
        cls.x = ElasticExperiment(3, cls.bank, cls.elastic)
        cls.o = SpatialObserver(cls.x)
        cls.g, cls.q = cls.x.e.geometry.group, cls.x.e.charges

    def test_conditioning_on_every_small_complete_charge_field_removes_charge_information(self):
        g, q = self.g, self.q
        fields = defaultdict(lambda: [0, 0, 0])
        for values in itertools.product(range(6), repeat=5):
            field = tuple(q[h] for h in values)
            if field.count(1) < 2 or field.count(1) % 2:
                continue
            p = 0
            for h in values:
                p = g.mul[p][h]
            for a, b in itertools.product(range(6), repeat=2):
                commutator = g.mul[g.inv[b]][g.mul[g.inv[a]][g.mul[b][a]]]
                if commutator == p and len(subgroup(g, (a, b, *values))) == 6:
                    row = fields[field]
                    row[0] += 1
                    row[1] += field[:2] == (1, 1) and values[0] != values[1]
                    row[2] += field[:3] == (1, 1, 1) and ((values[0] == values[1]) != (values[1] == values[2]))
        self.assertGreater(len(fields), 40)
        for field, (count, pairs, fans) in fields.items():
            beta, alpha = conditioned_coefficients(field.count(1), field.count(2))
            self.assertEqual(F(pairs, count)-beta*int(field[:2] == (1, 1)), 0)
            self.assertEqual(F(fans, count)-alpha*int(field[:3] == (1, 1, 1)), 0)
        self.assertEqual(conditioned_coefficients(2, 0), (F(9, 16), F()))
        self.assertEqual(conditioned_coefficients(4, 0), (F(27, 40), F(9, 20)))
        for bad in ((0, 0), (3, 0), (4, -1), (True, 1)):
            with self.assertRaises(ValueError):
                conditioned_coefficients(*bad)

    def test_group_convolution_independently_recovers_the_character_formula(self):
        g, q = self.g, self.q
        comm = Counter(g.mul[g.inv[b]][g.mul[g.inv[a]][g.mul[b][a]]] for a, b in itertools.product(range(6), repeat=2))
        for n1, n2 in itertools.product((2, 4, 6), range(4)):
            beta, alpha = conditioned_coefficients(n1, n2)
            for k, expected in ((2, beta), (3, alpha)):
                if n1 < k:
                    continue
                products = Counter({0: 1})
                for charge in [1]*(n1-k)+[2]*n2:
                    updated = Counter()
                    for a, count in products.items():
                        for h in range(6):
                            if q[h] == charge:
                                updated[g.mul[a][h]] += count
                    products = updated
                all_weight, active_weight = 0, 0
                for values in itertools.product((1, 2, 5), repeat=k):
                    p = 0
                    for h in values:
                        p = g.mul[p][h]
                    weight = sum(count*comm[g.mul[p][h]] for h, count in products.items())
                    if len(set(values)) == 1 and n2 == 0:
                        weight -= 4
                    all_weight += weight
                    active = values[0] != values[1] if k == 2 else (values[0] == values[1]) != (values[1] == values[2])
                    active_weight += weight*active
                self.assertEqual(F(active_weight, all_weight), expected)

    def test_full_field_conditioning_reduces_or_equals_the_old_residual_variance(self):
        for faces in (18, 72, 288):
            variance = conditioned_variances(faces, faces)
            old = conditional_marginals(self.g, self.q, faces, faces)
            for index, name in ((1, 'distinct_reflection_pair'), (2, 'active_reflection_fan')):
                self.assertGreater(variance[index], 0)
                self.assertLessEqual(variance[index], F(*old['observables'][name]['residual_variance']))

    def test_translation_channels_cover_real_supports_and_commute_with_raw_updates(self):
        x, o, rng = self.x, self.o, random.Random(895614)
        translation = o.translations
        self.assertEqual([layout.shape for layout in translation.layouts], [(2, 3, 3), (6, 3, 3), (6, 3, 3)])
        for layout in translation.layouts:
            self.assertEqual(sorted(layout.flat), list(range(layout.size)))
        initial = [rng.randrange(6) for _ in x.e.geometry.edges]
        # Use a valid reflection-present charge sector for the observer.
        initial = AxialReference(3, self.bank, 18).sample(rng, 'nonabelian_reflections')['links']
        original = o.measure(initial)
        for dy, dx in itertools.product(range(3), repeat=2):
            shifted = translation.links(initial, dx, dy)
            measured = o.measure(shifted)
            for a, b in zip(original['fields'], measured['fields']):
                np.testing.assert_allclose(np.roll(a, (dy, dx), axis=(1, 2)), b, atol=1e-12)
            for rule in range(15):
                kind = 1 if rule < 3 else 2
                layout = translation.layouts[kind]
                patch = int(layout[0, 0, 0])
                target_patch = int(layout[0, dy, dx])
                a, b = initial[:], shifted[:]
                x.apply(a, rule*x.supports+patch)
                x.apply(b, rule*x.supports+target_patch)
                self.assertEqual(translation.links(a, dx, dy), b)

    def test_dual_laplacian_fourier_blocks_match_the_actual_graph(self):
        x, o, rng = self.x, self.o, random.Random(195042)
        raw = np.array([rng.randrange(-4, 5) for _ in range(18)])
        layout = o.translations.layouts[0]
        actual = (x.geometry.l1@raw)[layout]
        fft = np.fft.fft2(raw[layout])
        transformed = np.einsum('yxab,byx->ayx', o.blocks, fft)
        np.testing.assert_allclose(np.fft.ifft2(transformed).real, actual, atol=1e-12)
        np.testing.assert_allclose(o.eigenvalues.sum(axis=-1), 6, atol=1e-12)
        self.assertAlmostEqual(o.eigenvalues[0, 0, 0], 0)
        self.assertAlmostEqual(o.eigenvalues[0, 0, 1], 6)

    def test_fft_correlation_matches_every_direct_displacement_and_finds_a_moving_zero_sum_pattern(self):
        o, s = self.o, self.o.translations.side
        fields = [np.zeros(layout.shape) for layout in o.translations.layouts]
        for field in fields:
            field[0, 0, 0] = 1; field[0, 0, 1] = -1
        pack = lambda values: {'fields': values, 'fourier': [np.fft.fft2(v) for v in values]}
        shifted = [np.roll(v, (1, 2), axis=(1, 2)) for v in fields]
        first, later = pack(fields), pack(shifted)
        maps, _, imaginary = o.correlate(first, later)
        for index, (a, b, variance) in enumerate(zip(fields, shifted, o.variances)):
            for dy, dx in itertools.product(range(s), repeat=2):
                direct = np.sum(a*np.roll(b, (-dy, -dx), axis=(1, 2)))/(a.size*float(variance))
                self.assertAlmostEqual(maps[index, dy, dx], direct)
            self.assertEqual(np.unravel_index(np.argmax(maps[index]), (s, s)), (1, 2))
            self.assertAlmostEqual(maps[index].sum(), 0)
        self.assertTrue(any(abs(v) > 1e-8 for v in imaginary))

    def test_even_time_spectral_constraints_and_jensen_bound_include_negative_one_step_eigenvalues(self):
        # An explicit reversible stochastic matrix with a negative eigenvalue.
        p = np.array([[F(1, 4), F(3, 4)], [F(3, 4), F(1, 4)]], dtype=object)
        h = np.array([F(1), F(-1)], dtype=object)
        c = [sum(h*(np.linalg.matrix_power(p, 2*n)@h))/2 for n in range(12)]
        self.assertEqual(c, [F(1, 4)**n for n in range(12)])
        for difference in range(6):
            values = c[:]
            for _ in range(difference):
                values = [values[i]-values[i+1] for i in range(len(values)-1)]
            self.assertTrue(all(v >= 0 for v in values))
        for n in range(10):
            self.assertGreaterEqual(c[n]*c[n+2], c[n+1]**2)
        # Multiple modes: Jensen uses the measured one-step moment, not a
        # guessed exponential fit to even-time samples.
        poles, weights = (F(-1, 2), F(4, 5)), (F(1, 3), F(2, 3))
        mean = sum(p*w for p, w in zip(poles, weights))
        for n in range(1, 10):
            self.assertGreaterEqual(sum(w*p**(2*n) for p, w in zip(poles, weights)), mean**(2*n))

    def test_running_uncertainty_uses_independent_replicates(self):
        data = np.arange(120).reshape(10, 3, 4).astype(float)
        moment = Moments()
        for row in data:
            moment.add(row)
        result = moment.result()
        np.testing.assert_allclose(result['mean'], data.mean(axis=0))
        np.testing.assert_allclose(result['standard_error'], data.std(axis=0, ddof=1)/10**0.5)
        with self.assertRaises(ValueError):
            Moments().result()

    def test_every_saved_summary_counts_independent_pairs(self):
        saved = json.loads((ROOT/'data/triangle-spatial-memory.json').read_text())
        for row in saved['whole_mesh']:
            self.assertEqual(row['trials'], len(row['replicates']))
            for key, target in (('shell_correlations', 'shell_estimates'),
                                ('oriented_mode_imaginary_parts', 'imaginary_estimates')):
                arrays = {name: np.array([r['runs'][name][key] for r in row['replicates']])
                          for name in ('baseline_padded', 'elastic')}
                arrays['paired_elastic_minus_baseline'] = arrays['elastic']-arrays['baseline_padded']
                for name, a in arrays.items():
                    np.testing.assert_allclose(a.mean(axis=0), row[target][name]['mean'])
                    np.testing.assert_allclose(a.std(axis=0, ddof=1)/len(a)**.5,
                                               row[target][name]['standard_error'])

    def test_complete_observation_is_gauge_invariant(self):
        x, g = self.x, self.g
        rng = random.Random(881903)
        sampler = AxialReference(3, self.bank, x.geometry.faces)
        for _ in range(6):
            raw = sampler.sample(rng, 'nonabelian_reflections')['links']
            frames = [rng.randrange(g.n) for _ in range(x.e.geometry.size)]
            changed = [g.mul[frames[v]][g.mul[h][g.inv[frames[u]]]]
                       for (u, v), h in zip(x.e.geometry.edges, raw)]
            first, second = self.o.measure(raw), self.o.measure(changed)
            self.assertEqual(first['populations'], second['populations'])
            for a, b in zip(first['fields'], second['fields']):
                np.testing.assert_array_equal(a, b)

    def test_even_time_cross_channel_matrix_is_psd(self):
        # Independent complex features; retain off-diagonal orientation terms.
        p = np.array([[.25, .75], [.75, .25]])
        features = np.array([[1, 2j, 3+1j], [-2j, 1, -1+2j]])
        for n in range(7):
            kernel = features.conj().T@np.linalg.matrix_power(p, 2*n)@features/2
            np.testing.assert_allclose(kernel, kernel.conj().T)
            self.assertGreaterEqual(np.linalg.eigvalsh(kernel).min(), -1e-12)
            for i, j in itertools.product(range(3), repeat=2):
                self.assertLessEqual(abs(kernel[i, j])**2,
                                     (kernel[i, i]*kernel[j, j]).real+1e-12)

    def test_saved_wavevector_runs_and_analytic_bounds_reproduce(self):
        path = ROOT/'data/triangle-spatial-memory.json'
        if not path.exists():
            self.fail('the checked spatial dataset is required')
        saved = json.loads(path.read_text())
        for row in saved['whole_mesh']:
            x = ElasticExperiment(row['side'], self.bank, self.elastic)
            o = SpatialObserver(x)
            self.assertEqual(o.shells, row['shell_wavevectors_yx'])
            self.assertEqual([list((v.numerator, v.denominator)) for v in o.variances], row['field_variances'])
            control = ElasticExperiment(row['side'], self.bank, [list(range(36))]*2)
            sampler = AxialReference(row['side'], self.bank, row['faces'])
            fresh, maps = trial(o, control, sampler, 0, row['times_in_attempts_per_face'])
            old = row['replicates'][0]
            for name, run in fresh['runs'].items():
                for key in ('events', 'populations', 'final_links_sha256'):
                    self.assertEqual(run[key], old['runs'][name][key])
                for key in ('shell_correlations', 'oriented_mode_imaginary_parts'):
                    np.testing.assert_allclose(run[key], old['runs'][name][key], atol=1e-11, rtol=1e-11)
            _, _, _, k1, k2 = coefficients(ChargeMarginal(row['faces'], row['faces'], 'nonabelian_reflections'))
            rates = [float(k1)*v+float(k2)*(6*v-v*v) for v in [*row['acoustic_shell_laplacian_eigenvalues'], 6]]
            bounds = [[(1-rate/x.operator_count)**(t*row['faces']) for rate in rates] for t in row['times_in_attempts_per_face']]
            np.testing.assert_allclose(bounds, row['charge_jensen_lower_bounds'])


if __name__ == '__main__':
    unittest.main()
