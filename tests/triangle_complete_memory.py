import itertools
import json
import random
import sys
import unittest
from collections import Counter, defaultdict
from fractions import Fraction as F
from pathlib import Path

import numpy as np
from flint import fmpq_mat

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'tools'))
from triangle_complete_memory import CompleteReference, CompleteObserver, run_pair, rotation_gate
from triangle_elastic_scattering import ElasticExperiment
from triangle_patch_observer import PatchUnion
from triangle_reference import AxialReference, subgroup


class CompleteMemoryTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.bank = json.loads((ROOT/'data/triangle-feedback.json').read_text())['bank']
        cls.elastic = json.loads((ROOT/'data/triangle-elastic-scattering.json').read_text())['pair_census']['elastic_tables']
        cls.x = ElasticExperiment(3, cls.bank, cls.elastic)
        cls.o = CompleteObserver(cls.x)
        cls.r, cls.g, cls.q = cls.o.reference, cls.x.e.geometry.group, cls.x.e.charges

    def test_complete_conditional_distribution_matches_all_small_torus_presentations(self):
        g, q = self.g, self.q
        for faces, charge in ((4, 4), (5, 4), (5, 6)):
            ref = CompleteReference(g, q, faces, charge)
            counts = defaultdict(Counter)
            for values in itertools.product(range(6), repeat=faces):
                field = tuple(q[h] for h in values)
                if sum(field) != charge or 1 not in field:
                    continue
                p = 0
                for h in values:
                    p = g.mul[p][h]
                state = ref.ids[ref.words.canonical(values[:3])]
                for a, b in itertools.product(range(6), repeat=2):
                    comm = g.mul[g.inv[b]][g.mul[g.inv[a]][g.mul[b][a]]]
                    if comm == p and len(subgroup(g, (a, b, *values))) == 6:
                        counts[field][state] += 1
            total = sum(sum(c.values()) for c in counts.values())
            self.assertEqual(total, ref.counts['canonical_connections'])
            covariance = np.full((49, 49), F(), dtype=object)
            for field, observed in counts.items():
                pop = tuple(field.count(i) for i in range(3))
                ids = ref.buckets[field[:3]]
                expected = ref.distribution(pop, field[:3])
                n = sum(observed.values())
                self.assertEqual(tuple(F(observed[i], n) for i in ids), expected)
                # Direct residual outer products, independently from the
                # population-weighted diag(p)-p p^T implementation.
                for i in ids:
                    residual = [F(int(i == j))-p for j, p in zip(ids, expected)]
                    for a, j in enumerate(ids):
                        for b, k in enumerate(ids):
                            covariance[j, k] += F(observed[i], total)*residual[a]*residual[b]
            np.testing.assert_array_equal(covariance, ref.covariance)

    def test_exact_rank_and_charge_nullspace_then_whitening(self):
        r = self.r
        self.assertEqual((len(r.reps), len(r.buckets), r.dimension), (49, 27, 22))
        matrix = fmpq_mat([[str(v) for v in row] for row in r.covariance])
        self.assertEqual(matrix.rank(), 22)
        for ids in r.buckets.values():
            indicator = np.array([int(i in ids) for i in range(49)], dtype=object)
            self.assertTrue(all(v == 0 for v in r.covariance@indicator))
        np.testing.assert_allclose(r.whitener.T@np.array(r.covariance, float)@r.whitener,
                                   np.eye(22), atol=1e-12)

    def test_real_patch_orbits_have_no_missing_gauge_state(self):
        union = PatchUnion(self.x.e, (0,))
        states = []
        for rep in union.reps:
            code = self.x.e.factor.oracle.code(union.lift(rep), 1, 0)
            states.append(int(self.r.codes[code]))
        self.assertEqual(sorted(states), list(range(49)))

    def test_conditional_projection_loses_only_charge_functions(self):
        r = self.r
        for pop in ((8, 2, 8), (9, 4, 5), (0, 18, 0)):
            vectors, valid = r.lookup(pop)
            for pattern, ids in r.buckets.items():
                if not all(valid[ids]):
                    continue
                p = np.array(r.distribution(pop, pattern), float)
                np.testing.assert_allclose(p@vectors[ids], 0, atol=1e-12)
                # Within a charge bucket all genuinely distinct states remain
                # distinct; no nonlinear local invariant can be merged away.
                for a, b in itertools.combinations(ids, 2):
                    self.assertGreater(np.linalg.norm(vectors[a]-vectors[b]), 1e-5)
        # Rotation triples have identical charges and zero previous RR/RRR
        # flags, but multiple inequivalent relative alignments.
        ids = r.buckets[(2, 2, 2)]
        self.assertEqual(len(ids), 4)
        vectors, _ = r.lookup((8, 2, 8))
        self.assertGreater(np.linalg.norm(vectors[ids[0]]-vectors[ids[1]]), 1)

    def test_kernel_is_independent_of_the_contrast_basis(self):
        r = self.r
        c = np.array(r.covariance, float)
        # Independent reference via the full 49-category Moore-Penrose inverse.
        inverse = np.linalg.pinv(c, rcond=1e-12, hermitian=True)
        w = r.whitener
        for pop in ((8, 2, 8), (0, 18, 0)):
            residuals = []
            for pattern, ids in r.buckets.items():
                if any(pattern.count(q) > n for q, n in enumerate(pop)):
                    continue
                p = np.array(r.distribution(pop, pattern), float)
                for i in ids:
                    residual = np.zeros(49); residual[ids] -= p; residual[i] += 1
                    residuals.append(residual)
            a = np.array(residuals)
            np.testing.assert_allclose((a@w)@(a@w).T, a@inverse@a.T, atol=1e-10)

    def test_raw_gauge_transformations_and_all_translations_preserve_observation(self):
        x, o, g = self.x, self.o, self.g
        rng = random.Random(881150)
        raw = AxialReference(3, self.bank, 18).sample(rng, 'nonabelian_reflections')['links']
        first = o.measure(raw)
        frames = [rng.randrange(6) for _ in range(9)]
        changed = [g.mul[frames[v]][g.mul[h][g.inv[frames[u]]]]
                   for (u, v), h in zip(x.e.geometry.edges, raw)]
        np.testing.assert_array_equal(first['states'], o.measure(changed)['states'])
        np.testing.assert_array_equal(first['fields'], o.measure(changed)['fields'])
        for dy, dx in itertools.product(range(3), repeat=2):
            translated = o.measure(o.translations.links(raw, dx, dy))
            np.testing.assert_allclose(translated['fields'], np.roll(first['fields'], (dy, dx), axis=(2, 3)))

    def test_complete_fft_kernel_matches_direct_all_displacement_comparison(self):
        o = self.o
        rng = np.random.default_rng(818150)
        shape = (6, 22, 3, 3)
        a, b = rng.normal(size=shape), rng.normal(size=shape)
        _, spatial = o.correlate({'fourier': np.fft.fft2(a)}, {'fourier': np.fft.fft2(b)})
        for dy, dx in itertools.product(range(3), repeat=2):
            direct = np.sum(a*np.roll(b, (-dy, -dx), axis=(2, 3)))/a.size
            self.assertAlmostEqual(spatial[dy, dx], direct)

    def test_vacancy_free_rotation_gate_is_exhaustive_and_falsifiable(self):
        result = rotation_gate(self.x)
        self.assertEqual(sum(r['raw_inputs'] for r in result['primitive_rows']), 2700)
        self.assertGreater(sum(r['enabled_rotation_inputs'] for r in result['primitive_rows']), 0)
        self.assertEqual(result['primitive_rows'][1]['enabled_rotation_inputs'], 0)
        self.assertEqual(result['primitive_rows'][2]['enabled_rotation_inputs'], 0)
        for row in result['primitive_rows']:
            for p in row['enabled_rotation_charge_patterns']:
                self.assertEqual(sorted(p['pattern']), [0, 2] if row['arity'] == 2 else [0, 1, 2])
        # Inject a known forbidden change in a copy to make the test exercise
        # rejection, not merely restate that the current census has no failures.
        from types import SimpleNamespace
        tables = [row[:] for row in self.x.tables]
        tables[1][3*6+3] = 4*6+4
        with self.assertRaisesRegex(ValueError, 'lacks a vacancy'):
            rotation_gate(SimpleNamespace(e=self.x.e, tables=tables))

    def test_saved_statistics_and_one_pair_per_size_reproduce(self):
        saved = json.loads((ROOT/'data/triangle-complete-memory.json').read_text())
        for row in saved['whole_mesh']:
            x = ElasticExperiment(row['side'], self.bank, self.elastic)
            o = CompleteObserver(x)
            self.assertEqual(json.loads(json.dumps(o.reference.report())), row['reference'])
            self.assertEqual(fmpq_mat([[str(v) for v in a] for a in o.reference.gram]).rank(), 22)
            self.assertEqual(json.loads(json.dumps(rotation_gate(x))), row['rotation_gate'])
            control = ElasticExperiment(row['side'], self.bank, [list(range(36))]*2)
            sampler = AxialReference(row['side'], self.bank, row['faces'])
            fresh, _ = run_pair(o, control, sampler, 0, row['times_in_attempts_per_face'])
            old = row['replicates'][0]
            for name, run in fresh['runs'].items():
                for key in ('events', 'populations', 'final_links_sha256'):
                    self.assertEqual(run[key], old['runs'][name][key])
                for key in run['observations'][0]:
                    np.testing.assert_allclose([a[key] for a in run['observations']],
                                               [a[key] for a in old['runs'][name]['observations']], atol=1e-10)
            for key, estimates in row['estimates'].items():
                arrays = {name: np.array([[a[key] for a in r['runs'][name]['observations']] for r in row['replicates']])
                          for name in ('baseline_padded', 'elastic')}
                arrays['paired_elastic_minus_baseline'] = arrays['elastic']-arrays['baseline_padded']
                for name, values in arrays.items():
                    np.testing.assert_allclose(values.mean(axis=0), estimates[name]['mean'])
                    np.testing.assert_allclose(values.std(axis=0, ddof=1)/len(values)**.5, estimates[name]['standard_error'])


if __name__ == '__main__':
    unittest.main()
