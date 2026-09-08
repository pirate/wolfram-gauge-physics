import itertools
import json
import random
import sys
import unittest
from collections import Counter
from fractions import Fraction as F
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'tools'))
from triangle_elastic_scattering import ElasticExperiment
from triangle_reference import AxialReference, reference, subgroup
from triangle_relational_memory import (RelationalObserver, conditional_marginals, patch_calibration,
                                        refine_partition, replay_samples, whole_mesh)


class RelationalMemoryTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.bank = json.loads((ROOT/'data/triangle-feedback.json').read_text())['bank']
        cls.source = json.loads((ROOT/'data/triangle-elastic-scattering.json').read_text())
        cls.data = json.loads((ROOT/'data/triangle-relational-memory.json').read_text())
        cls.elastic = cls.source['pair_census']['elastic_tables']
        cls.x = ElasticExperiment(3, cls.bank, cls.elastic)
        cls.g, cls.q = cls.x.e.geometry.group, cls.x.e.charges

    def test_exact_patch_calibration_reproduces_and_charge_requires_one_extra_bit(self):
        p = patch_calibration(self.bank, self.elastic, self.source)
        self.assertEqual(p, self.data['patch_calibration'])
        levels = p['charge_partition_refinement']
        self.assertEqual([len(set(v)) for v in levels], [7, 8])
        desired = [(tuple(q), h > 0) for q, h in zip(p['charge_values'], p['relational_residual'])]
        for i, j in itertools.product(range(10), repeat=2):
            self.assertEqual(levels[-1][i] == levels[-1][j], desired[i] == desired[j])
        counts = p['transition_counts']['elastic']
        for i, j in itertools.product(range(10), repeat=2):
            if levels[-1][i] == levels[-1][j]:
                for block in set(levels[-1]):
                    self.assertEqual(sum(counts[i][k] for k in range(10) if levels[-1][k] == block),
                                     sum(counts[j][k] for k in range(10) if levels[-1][k] == block))
        self.assertEqual(refine_partition([[1, 1], [1, 1]], [0, 0]), [[0, 0]])

    def test_rational_spectral_measures_match_direct_integer_matrix_powers(self):
        p = self.data['patch_calibration']
        observables = {'charge': [q[0]-1 for q in p['charge_values']], 'relational_residual': p['relational_residual']}
        for name, counts in p['transition_counts'].items():
            counts = np.array(counts, dtype=object)
            for label, values in observables.items():
                vector = np.array(values, dtype=object)
                evolved = vector.copy()
                for tick in range(24):
                    direct = F(int(vector@evolved), 10*24**tick)
                    spectral = sum(F(*r['covariance_weight'])*F(*r['pole'])**tick for r in p['spectral_measures'][name][label])
                    self.assertEqual(direct, spectral)
                    evolved = counts@evolved

    def test_frozen_plateau_is_removed_without_changing_any_charge_autocorrelation(self):
        p = self.data['patch_calibration']
        weights = lambda mode, label: {F(*r['pole']): F(*r['covariance_weight']) for r in p['spectral_measures'][mode][label] if r['covariance_weight'][0]}
        self.assertEqual(weights('baseline_padded', 'charge'), weights('elastic', 'charge'))
        self.assertEqual(weights('baseline_padded', 'relational_residual'), {F(1, 3): F(3, 20), F(1): F(1, 4)})
        parts = self.source['patch_activation']['component_mergers'][0]['previous_components']
        h = dict(zip(p['states'], p['relational_residual']))
        plateau = sum(F(len(part), 10)*F(sum(h[i] for i in part), len(part))**2 for part in parts)
        self.assertEqual(plateau, F(1, 4))
        self.assertEqual(weights('elastic', 'relational_residual'), {F(1, 6): F(1, 4), F(5, 6): F(3, 20)})
        for label, expected in (('charge', 7), ('relational_residual', 5)):
            w = weights('elastic', label)
            self.assertEqual(1+2*sum(weight*pole/(1-pole) for pole, weight in w.items())/sum(w.values()), expected)

    def test_disk_marginals_match_exhaustive_small_torus_presentations(self):
        g, q = self.g, self.q
        for faces, total_charge in ((4, 4), (5, 4), (5, 6)):
            total, conditioning, active = 0, Counter(), Counter()
            for values in itertools.product(range(6), repeat=faces):
                if sum(q[h] for h in values) != total_charge or not any(q[h] == 1 for h in values):
                    continue
                product = 0
                for h in values:
                    product = g.mul[product][h]
                for a, b in itertools.product(range(6), repeat=2):
                    commutator = g.mul[g.inv[b]][g.mul[g.inv[a]][g.mul[b][a]]]
                    if commutator != product or len(subgroup(g, (a, b, *values))) != 6:
                        continue
                    total += 1
                    for k in (2, 3):
                        if all(q[h] == 1 for h in values[:k]):
                            conditioning[k] += 1
                            active[k] += (values[0] != values[1] if k == 2 else (values[0] == values[1]) != (values[1] == values[2]))
            m = conditional_marginals(g, q, faces, total_charge)
            self.assertEqual(total, m['reference_connections'])
            for k, name in ((2, 'distinct_reflection_pair'), (3, 'active_reflection_fan')):
                row = m['observables'][name]
                self.assertEqual(F(*row['conditioning_probability']), F(conditioning[k], total))
                self.assertEqual(F(*row['conditional_mean']), F(active[k], conditioning[k]))
                alpha = F(active[k], conditioning[k])
                variance = (active[k]*(1-alpha)**2+(conditioning[k]-active[k])*alpha**2)/total
                self.assertEqual(F(*row['residual_variance']), variance)

    def test_whole_mesh_marginals_reproduce_without_independent_reflection_approximation(self):
        for row in self.data['whole_mesh']:
            m = conditional_marginals(self.g, self.q, row['faces'], row['charge'])
            self.assertEqual(m, row['exact_marginals'])
        small = conditional_marginals(self.g, self.q, 4, 4)
        self.assertNotEqual(F(*small['observables']['distinct_reflection_pair']['conditional_mean']), F(2, 3))
        self.assertNotEqual(F(*small['observables']['active_reflection_fan']['conditional_mean']), F(4, 9))

    def test_observer_is_gauge_invariant_and_integer_moments_equal_direct_rationals(self):
        x, g, rng = self.x, self.g, random.Random(7118)
        m = conditional_marginals(g, self.q, 18, 18)
        observer = RelationalObserver(x, m)
        sampler = AxialReference(3, self.bank, 18)
        for _ in range(12):
            raw = sampler.sample(rng, 'nonabelian_reflections')['links']
            frames = [rng.randrange(6) for _ in range(9)]
            gauged = [g.mul[frames[v]][g.mul[z][g.inv[frames[u]]]] for (u, v), z in zip(x.e.geometry.edges, raw)]
            a, b = observer.measure(raw), observer.measure(gauged)
            for first, second in zip(a, b):
                np.testing.assert_array_equal(first, second)
            later = raw[:]
            for _ in range(25):
                x.apply(later, rng.randrange(x.operator_count))
            b = observer.measure(later)
            values, moments = observer.correlation(a, b)
            centered = lambda sample: [[F(int(v))-observer.mean_charge for v in sample[0]]]+[[F(int(v[0]))-alpha*int(v[1]) for v in row] for row, alpha in zip(sample[1:], (observer.beta, observer.alpha))]
            aa, bb = centered(a), centered(b)
            for i in range(3):
                direct = sum(u*v for u, v in zip(aa[i], bb[i]))/len(aa[i])/observer.variances[i]
                self.assertEqual(values[i], float(direct))
            self.assertEqual([len(row) for row in moments], [3, 4, 4])
            permuted = [rng.sample(list(v), len(v)) for v in b]
            permuted = [np.array(v, dtype=np.int64) for v in permuted]
            self.assertEqual(observer.integrated_correlation(a, b), observer.integrated_correlation(a, permuted))

    def test_sampled_replay_matches_every_available_cpp_raw_snapshot(self):
        x, rng = self.x, random.Random(80345)
        initial = [rng.randrange(6) for _ in x.e.geometry.edges]
        observer = RelationalObserver(x, conditional_marginals(self.g, self.q, 18, 18))
        schedule = [rng.randrange(x.operator_count) for _ in range(120)]
        samples, check = replay_samples(x, initial, schedule, list(range(121)), observer)
        full = x.compiled([initial], schedule)['runs'][1]
        events = {row[0]: row[-1] for row in full['events']}
        raw = initial
        for tick, measured in enumerate(samples):
            raw = events.get(tick, raw)
            for a, b in zip(measured, observer.measure(raw)):
                np.testing.assert_array_equal(a, b)
        self.assertTrue(check['cpp_all_event_targets_final_state_and_inverse_verified'])
        with self.assertRaises(ValueError):
            replay_samples(x, initial, schedule, [0, 120, 120], observer)

    def test_saved_correlations_use_independent_replicates_and_paired_uncertainty(self):
        for row in self.data['whole_mesh']:
            n = row['independent_trials']
            self.assertEqual(len(row['replicates']), n)
            self.assertEqual(len({r['initial_seed'] for r in row['replicates']}), n)
            marginal = row['exact_marginals']
            alpha = F(*marginal['observables']['active_reflection_fan']['conditional_mean'])
            beta = F(*marginal['observables']['distinct_reflection_pair']['conditional_mean'])
            variances = [F(*marginal['charge_variance'])]+[F(*marginal['observables'][name]['residual_variance'])
                          for name in ('distinct_reflection_pair', 'active_reflection_fan')]
            # Reconstruct EVERY reported covariance from saved integer counts,
            # not merely its aggregate from previously exported floats.
            for replicate in row['replicates']:
                for run in replicate['runs'].values():
                    initial = run['spatial_indicator_totals'][0]
                    for tick, (moments, totals) in enumerate(zip(run['integer_moment_counts'], run['spatial_indicator_totals'])):
                        local = [F(moments[0][0]-moments[0][1]-moments[0][2], row['faces'])+1]
                        integrated = [F((initial[0][0]-row['faces'])*(totals[0][0]-row['faces']), row['faces'])]
                        for index, coefficient in ((1, beta), (2, alpha)):
                            a, b, c, d = moments[index]
                            local.append((a-coefficient*(b+c)+coefficient**2*d)/(3*row['faces']))
                            first = initial[index][0]-coefficient*initial[index][1]
                            later = totals[index][0]-coefficient*totals[index][1]
                            integrated.append(first*later/(3*row['faces']))
                        self.assertEqual([float(v/variance) for v, variance in zip(local, variances)], run['normalized_autocovariances'][tick])
                        self.assertEqual([float(v/variance) for v, variance in zip(integrated, variances)], run['integrated_correlations'][tick])
            arrays = {name: np.array([r['runs'][name]['normalized_autocovariances'] for r in row['replicates']]) for name in ('baseline_padded', 'elastic')}
            np.testing.assert_array_equal(arrays['baseline_padded'][:, 0], arrays['elastic'][:, 0])
            arrays['paired_elastic_minus_baseline'] = arrays['elastic']-arrays['baseline_padded']
            for name, a in arrays.items():
                np.testing.assert_allclose(a.mean(axis=0), row['estimates'][name]['mean'], atol=1e-14, rtol=1e-14)
                np.testing.assert_allclose(a.std(axis=0, ddof=1)/n**0.5, row['estimates'][name]['standard_error'], atol=1e-14, rtol=1e-14)
            self.assertEqual(row['attempt_ticks'], [row['faces']*t for t in row['times_in_attempts_per_face']])
            self.assertTrue(all(r['initial_spatial_means'][0] == 0 for r in row['replicates']))
            integrated = {name: np.array([r['runs'][name]['integrated_correlations'] for r in row['replicates']]) for name in ('baseline_padded', 'elastic')}
            integrated['paired_elastic_minus_baseline'] = integrated['elastic']-integrated['baseline_padded']
            for name, a in integrated.items():
                np.testing.assert_array_equal(a[:, :, 0], np.zeros((n, len(row['attempt_ticks']))))
                np.testing.assert_allclose(a.mean(axis=0), row['integrated_estimates'][name]['mean'], atol=1e-14, rtol=1e-14)
                np.testing.assert_allclose(a.std(axis=0, ddof=1)/n**0.5, row['integrated_estimates'][name]['standard_error'], atol=1e-14, rtol=1e-14)

    def test_selected_saved_raw_trajectories_reproduce_on_each_mesh_size(self):
        for stored in self.data['whole_mesh']:
            # Reproduce the first two independent, paired trajectories in full;
            # compare integer sufficient statistics and raw hashes, not floats.
            fresh = whole_mesh(stored['side'], self.bank, self.elastic, 2, stored['times_in_attempts_per_face'])
            for a, b in zip(fresh['replicates'], stored['replicates'][:2]):
                for name in ('baseline_padded', 'elastic'):
                    for key in ('integer_moment_counts', 'spatial_indicator_totals', 'events', 'final_links_sha256'):
                        self.assertEqual(a['runs'][name][key], b['runs'][name][key])


if __name__ == '__main__':
    unittest.main()
