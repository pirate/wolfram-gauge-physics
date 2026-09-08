import itertools
import json
import random
import sys
import unittest
from fractions import Fraction as F
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'tools'))
from triangle_complete_memory import CompleteObserver
from triangle_charge_response import ChargeMarginal
from triangle_elastic_scattering import ElasticExperiment
from triangle_encounter_memory import EncounterObserver, run_pair
from triangle_reference import AxialReference
from triangle_relational_memory import replay_samples


class EncounterMemoryTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.bank = json.loads((ROOT/'data/triangle-feedback.json').read_text())['bank']
        cls.elastic = json.loads((ROOT/'data/triangle-elastic-scattering.json').read_text())['pair_census']['elastic_tables']
        cls.x = ElasticExperiment(3, cls.bank, cls.elastic)
        cls.complete = CompleteObserver(cls.x)
        cls.sampler = AxialReference(3, cls.bank, 18)

    def test_arbitrary_rotation_disks_have_uniform_relative_sign_orbits(self):
        ref = self.complete.reference
        g = self.x.e.geometry.group
        rotations = ref.words.rotations
        for k in range(1, 9):
            populations = (18-2-k, 2, k)
            counts = {}
            for values in itertools.product(rotations, repeat=k):
                orbit = ref.words.canonical(values)
                counts[orbit] = counts.get(orbit, 0)+ref.completion(values, populations)
            self.assertEqual(len(counts), 2**(k-1))
            self.assertEqual(len(set(counts.values())), 1)
            self.assertEqual(sum(counts.values()), 12*3**2*2**k)
        observer = EncounterObserver(self.x, self.complete)
        self.assertEqual(observer.rotation_probability,
                         ChargeMarginal(18, 18, 'nonabelian_reflections').probability((0, 0, 3)))
        # Independent exact four-state kernel, normalized by charge-222 chance.
        pop = (8, 2, 8)
        values, _ = ref.lookup(pop)
        ids = observer.rotation_states
        expected = (4*np.eye(4)-1)/float(observer.rotation_probability)
        np.testing.assert_allclose(values[ids]@values[ids].T, expected, atol=1e-11)

    def test_first_encounter_operator_identity_without_commuting_kernels(self):
        p = np.array([[F(1, 2), F(1, 3), F(1, 6), F()],
                      [F(1, 3), F(1, 2), F(), F(1, 6)],
                      [F(1, 6), F(), F(1, 3), F(1, 2)],
                      [F(), F(1, 6), F(1, 2), F(1, 3)]], dtype=object)
        labels = [0, 0, 1, 1]
        physical = np.array([[p[i, j]*int(labels[i] == labels[j]) for j in range(4)] for i in range(4)], dtype=object)
        raw = np.diag(np.diag(p))
        self.assertFalse(np.array_equal(p@physical, physical@p))
        f = np.array([F(1), F(-1), F(2), F(-2)], dtype=object)
        for n in range(5):
            pn, kn = np.linalg.matrix_power(p, n), np.linalg.matrix_power(physical, n)
            expansion = kn.copy()
            for j in range(n):
                expansion += np.linalg.matrix_power(physical, j)@(p-physical)@np.linalg.matrix_power(p, n-j-1)
            np.testing.assert_array_equal(expansion, pn)
            # Independently sum actual histories, including paths that return
            # to their starting state after crossing a physical-state boundary.
            sums = [F(), F(), F()]
            for path in itertools.product(range(4), repeat=n+1):
                weight = F(1, 4)
                for a, b in zip(path, path[1:]):
                    weight *= p[a, b]
                touched = any(a != b for a, b in zip(path, path[1:]))
                changed = any(labels[a] != labels[b] for a, b in zip(path, path[1:]))
                sums[2 if changed else int(touched)] += weight*f[path[0]]*f[path[-1]]
            rn = np.linalg.matrix_power(raw, n)
            expected = [f@rn@f/4, f@(kn-rn)@f/4, f@(pn-kn)@f/4]
            self.assertEqual(sums, expected)

    def test_incremental_visits_match_full_raw_and_orbit_history_at_every_attempt(self):
        x, c = self.x, self.complete
        rng = random.Random(918382)
        raw = self.sampler.sample(rng, 'nonabelian_reflections')['links']
        observer = EncounterObserver(x, c)
        observer.measure(raw)
        first_raw, first_orbit = np.zeros(x.supports, int), np.zeros(x.supports, int)
        changes = np.zeros(x.supports, int)
        for tick in range(1, 151):
            old = raw[:]
            before = c.measure(old)['states']
            op = rng.randrange(x.operator_count)
            code, target = x.apply(raw, op)
            after = c.measure(raw)['states']
            edges = {i for i, (a, b) in enumerate(zip(old, raw)) if a != b}
            for patch, reads in enumerate(x.e.factor.oracle.fan.reads):
                if edges & reads and not first_raw[patch]:
                    first_raw[patch] = tick
                if before[patch] != after[patch]:
                    if not first_orbit[patch]:
                        first_orbit[patch] = tick
                    changes[patch] += 1
            if code != target:
                observer.event(raw, tick, op, code, target)
            else:
                self.assertEqual(raw, old)
            observer.measure(raw)
            np.testing.assert_array_equal(observer.first_raw, first_raw)
            np.testing.assert_array_equal(observer.first_orbit, first_orbit)
            np.testing.assert_array_equal(observer.orbit_changes, changes)

    def test_actual_rotation_change_and_inverse_is_not_misclassified_as_untouched(self):
        x, c = self.x, self.complete
        rng = random.Random(81744)
        witness = None
        for _ in range(20):
            raw = self.sampler.sample(rng, 'nonabelian_reflections')['links']
            initial = c.measure(raw)
            rotation = [p for p, i in enumerate(initial['states']) if c.reference.patterns[i] == (2, 2, 2)]
            if not rotation:
                continue
            for op in range(x.operator_count):
                out = raw[:]
                code, target = x.apply(out, op)
                if code == target:
                    continue
                after = c.measure(out)
                changed = [p for p in rotation if initial['states'][p] != after['states'][p]]
                if changed:
                    witness = raw, op, changed
                    break
            if witness:
                break
        self.assertIsNotNone(witness)
        raw, op, changed = witness
        observer = EncounterObserver(x, c)
        samples, check = replay_samples(x, raw, [op, x.inverse_operator(op)], [0, 1, 2], observer, observer.event)
        np.testing.assert_array_equal(samples[0]['states'], samples[-1]['states'])
        np.testing.assert_allclose(samples[0]['fields'], samples[-1]['fields'])
        for p in changed:
            self.assertEqual(observer.first_raw[p], 1)
            self.assertEqual(observer.first_orbit[p], 1)
            self.assertEqual(observer.orbit_changes[p], 2)
        self.assertTrue(check['cpp_all_event_targets_final_state_and_inverse_verified'])
        self.assertGreater(samples[-1]['initial_rotation_bin_counts'][2], 0)
        self.assertEqual(samples[-1]['initial_rotation_bin_counts'][1], 0)

    def test_history_masks_sum_to_complete_detector_without_survivor_normalization(self):
        x, c = self.x, self.complete
        rng = random.Random(671882)
        raw = self.sampler.sample(rng, 'nonabelian_reflections')['links']
        schedule = [rng.randrange(x.operator_count) for _ in range(300)]
        observer = EncounterObserver(x, c)
        samples, _ = replay_samples(x, raw, schedule, [0, 30, 100, 300], observer, observer.event)
        for sample in samples:
            split, maps = observer.correlate(samples[0], sample)
            ordinary, spatial = c.correlate(samples[0], sample)
            np.testing.assert_allclose(maps[0].sum(axis=0), spatial, atol=1e-12)
            self.assertAlmostEqual(sum(split['local_contributions'][0]), ordinary['local_trace'])
            self.assertAlmostEqual(sum(split['lowest_shell_contributions'][0]), ordinary['lowest_shell_trace'])
            self.assertAlmostEqual(sum(split['history_bin_fractions']), 1)
            # The rotation-only charge block cannot have a raw-touch history
            # with its orbit unchanged: its first touch must leave charge 222.
            np.testing.assert_array_equal(maps[1, 1], 0)
            both = np.array(split['rotation_both_endpoint_counts'])
            same = np.array(split['rotation_same_endpoint_counts'])
            direct = (4*same-both)/(3*x.supports*float(observer.rotation_probability))
            np.testing.assert_allclose(direct, split['local_contributions'][1], atol=1e-11)
            self.assertEqual(int(both[0]), sample['initial_rotation_bin_counts'][0])
            self.assertEqual(int(same[0]), int(both[0]))

    def test_gauge_frames_preserve_actual_visitation_and_every_correlation_bin(self):
        x, c, g = self.x, self.complete, self.x.e.geometry.group
        rng = random.Random(71431)
        raw = self.sampler.sample(rng, 'nonabelian_reflections')['links']
        frames = [rng.randrange(6) for _ in range(9)]
        gauged = [g.mul[frames[v]][g.mul[h][g.inv[frames[u]]]] for (u, v), h in zip(x.e.geometry.edges, raw)]
        schedule = [rng.randrange(x.operator_count) for _ in range(300)]
        records = []
        for initial in (raw, gauged):
            observer = EncounterObserver(x, c)
            samples, _ = replay_samples(x, initial, schedule, [0, 30, 100, 300], observer, observer.event)
            records.append((observer.report(), [observer.correlate(samples[0], s) for s in samples]))
        self.assertEqual(records[0][0], records[1][0])
        for (a, am), (b, bm) in zip(records[0][1], records[1][1]):
            self.assertEqual(a, b)
            np.testing.assert_array_equal(am, bm)

    def test_saved_data_replay_first_pair_and_reconstruct_all_paired_uncertainties(self):
        data = json.loads((ROOT/'data/triangle-encounter-memory.json').read_text())
        for row in data['whole_mesh']:
            x = ElasticExperiment(row['side'], self.bank, self.elastic)
            c = CompleteObserver(x)
            control = ElasticExperiment(row['side'], self.bank, [list(range(36))]*2)
            sampler = AxialReference(row['side'], self.bank, row['faces'])
            fresh, _ = run_pair(c, control, sampler, 0, row['times_in_attempts_per_face'])
            for name, run in fresh['runs'].items():
                old = row['replicates'][0]['runs'][name]
                for key in ('events', 'populations', 'final_links_sha256', 'visit_history'):
                    self.assertEqual(run[key], old[key])
                for key in run['observations'][0]:
                    np.testing.assert_allclose([a[key] for a in run['observations']], [a[key] for a in old['observations']], atol=1e-10)
            for key, estimates in row['estimates'].items():
                arrays = {name: np.array([[a[key] for a in r['runs'][name]['observations']] for r in row['replicates']])
                          for name in ('baseline_padded', 'elastic')}
                arrays['paired_elastic_minus_baseline'] = arrays['elastic']-arrays['baseline_padded']
                for name, values in arrays.items():
                    np.testing.assert_allclose(values.mean(axis=0), estimates[name]['mean'])
                    np.testing.assert_allclose(values.std(axis=0, ddof=1)/len(values)**.5, estimates[name]['standard_error'])
                    if key in row['total_estimates']:
                        total = values.sum(axis=-1)
                        expected = row['total_estimates'][key][name]
                        np.testing.assert_allclose(total.mean(axis=0), expected['mean'], atol=1e-12)
                        # Summing paired differences versus differencing two
                        # sums can differ at roundoff when the true SE is zero.
                        np.testing.assert_allclose(total.std(axis=0, ddof=1)/len(values)**.5,
                                                   expected['standard_error'], atol=1e-12)

    def test_all_saved_integer_returns_and_first_visits_reconstruct_the_measurements(self):
        data = json.loads((ROOT/'data/triangle-encounter-memory.json').read_text())
        for row in data['whole_mesh']:
            n = 3*row['faces']
            probability = F(*row['rotation_222_probability'])
            self.assertEqual(probability, ChargeMarginal(row['faces'], row['faces'],
                             'nonabelian_reflections').probability((0, 0, 3)))
            q = self.x.e.charges
            rotation_ids = [i for i, values in enumerate(row['reference']['representatives'])
                            if all(q[h] == 2 for h in values)]
            for pair in row['replicates']:
                for run in pair['runs'].values():
                    history = run['visit_history']
                    raw = np.array(history['first_raw_change_attempt'])
                    orbit = np.array(history['first_orbit_change_attempt'])
                    initial_rotation = np.isin(history['initial_states'], rotation_ids)
                    for t, sample in zip(row['times_in_attempts_per_face'], run['observations']):
                        tick = t*row['faces']
                        labels = np.where((orbit > 0) & (orbit <= tick), 2,
                                          np.where((raw > 0) & (raw <= tick), 1, 0))
                        np.testing.assert_array_equal(np.bincount(labels, minlength=3)/n,
                                                      sample['history_bin_fractions'])
                        np.testing.assert_array_equal(np.bincount(labels[initial_rotation], minlength=3)/n,
                                                      sample['initial_rotation_bin_fractions_of_all_patches'])
                        both = np.array(sample['rotation_both_endpoint_counts'])
                        same = np.array(sample['rotation_same_endpoint_counts'])
                        self.assertTrue(np.all(same <= both))
                        self.assertEqual(int(both[0]), int(np.sum(initial_rotation & (labels == 0))))
                        self.assertEqual(int(same[0]), int(both[0]))
                        self.assertEqual(int(both[1]), 0)
                        direct = (4*same-both)/(3*n*float(probability))
                        np.testing.assert_allclose(direct, sample['local_contributions'][1], atol=1e-11)


if __name__ == '__main__':
    unittest.main()
