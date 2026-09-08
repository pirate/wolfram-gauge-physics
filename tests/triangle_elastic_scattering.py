import itertools
import json
import random
import subprocess
import sys
import unittest
from collections import Counter
from fractions import Fraction
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'tools'))
from triangle_elastic_scattering import ElasticExperiment, audit, homogeneous_state, inverse, order
from triangle_patch_observer import PatchUnion, WordObserver
from triangle_reference import AxialReference, constant_connection, subgroup
from triangle_response_memory import activity_mask


class TriangleElasticScatteringTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.bank = json.loads((ROOT/'data/triangle-feedback.json').read_text())['bank']
        cls.data = json.loads((ROOT/'data/triangle-elastic-scattering.json').read_text())
        cls.census = cls.data['pair_census']
        cls.x = ElasticExperiment(3, cls.bank, cls.census['elastic_tables'])
        cls.g, cls.q = cls.x.e.geometry.group, cls.x.e.charges

    def test_complete_census_kernel_raw_trajectories_and_response_reproduce(self):
        self.assertEqual(audit(), self.data)
        self.assertEqual(self.census['bijections'], 192)
        self.assertEqual(self.census['order_counts'], {'1': 1, '2': 47, '3': 2, '4': 16, '6': 94, '12': 32})
        self.assertEqual(sum(self.census['reversing_orientation_order_counts'].values()), 144)

    def test_every_law_and_every_product_frame_satisfies_all_constraints(self):
        g, q = self.g, self.q
        restrictions = {p: set() for p in set(g.sectors)}
        for t in self.census['tables']:
            self.assertEqual(sorted(t), list(range(36)))
            for a, b in itertools.product(range(6), repeat=2):
                u, v = divmod(t[6*a+b], 6)
                self.assertEqual(g.mul[a][b], g.mul[u][v])
                self.assertEqual(q[a]+q[b], q[u]+q[v])
                for row in g.conj:
                    self.assertEqual(t[6*row[a]+row[b]], 6*row[u]+row[v])
                if q[a] == q[b] == 1:
                    # Equivariance prevents conversion of two reflections to
                    # a vacancy/rotation pair, even without involutivity.
                    self.assertEqual((q[u], q[v]), (1, 1))
            for p in restrictions:
                restrictions[p].add(tuple(t[6*a+g.mul[g.inv[a]][p]]//6 for a in range(6)))
        # Exhaust ALL 6! choices on each representative product fiber,
        # independently of how the production census transported its frames.
        for p, found in restrictions.items():
            centralizer = [h for h in range(6) if g.conj[h][p] == p]
            cost = lambda a: q[a]+q[g.mul[g.inv[a]][p]]
            expected = {f for f in itertools.permutations(range(6))
                        if all(cost(a) == cost(f[a]) for a in range(6))
                        and all(f[g.conj[h][a]] == g.conj[h][f[a]] for h in centralizer for a in range(6))}
            self.assertEqual(found, expected)
        self.assertEqual([len(restrictions[p]) for p in sorted(restrictions)], [2, 16, 6])

    def test_two_unique_minimal_noninvolutions_reverse_rather_than_commute_with_orientation(self):
        candidates = self.census['elastic_tables']
        noninvolutions = [t for t in self.census['tables'] if order(t) > 2]
        support = lambda t: sum(i != x for i, x in enumerate(t))
        self.assertEqual(min(map(support, noninvolutions)), 6)
        self.assertEqual(sorted(t for t in noninvolutions if support(t) == 6), sorted(candidates))
        self.assertEqual(inverse(candidates[0]), candidates[1])
        self.assertEqual(self.census['elastic_braid_failures_on_all_triples'], [36, 36])
        j = [6*self.g.inv[b]+self.g.inv[a] for a, b in itertools.product(range(6), repeat=2)]
        for t in candidates:
            self.assertEqual(order(t), 3)
            self.assertEqual([j[t[j[x]]] for x in range(36)], inverse(t))
            self.assertNotEqual([t[j[x]] for x in range(36)], [j[t[x]] for x in range(36)])
            self.assertNotEqual([t[t[x]] for x in range(36)], list(range(36)))
            for code, target in enumerate(t):
                self.assertEqual(subgroup(self.g, divmod(code, 6)), subgroup(self.g, divmod(target, 6)))
        for bad in ([0, 0], [1], [-1, 0]):
            with self.assertRaises(ValueError):
                inverse(bad)

    def test_frozen_homogeneous_states_are_exactly_parallel_reflection_reductions(self):
        expected = [{'holonomy_group_order': 2, 'original_frozen': True, 'candidate_frozen': True, 'count': 12},
                    {'holonomy_group_order': 6, 'original_frozen': False, 'candidate_frozen': False, 'count': 90},
                    {'holonomy_group_order': 6, 'original_frozen': True, 'candidate_frozen': False, 'count': 6}]
        for record in self.data['homogeneous_frozen_reduction']['records']:
            self.assertEqual(record['all_charge_one_assignments'], 108)
            self.assertEqual(record['classification'], expected)
        for record in self.data['uniform_attempt_runs']:
            for check in record['nonabsorption_checks']:
                self.assertGreater(check['states_checked'], 1)
                self.assertTrue(check['full_S3_Q_equals_F_and_enabled_pair_throughout'])
        x, g, rng = self.x, self.g, random.Random(31829)
        reflections = [h for h, q in enumerate(self.q) if q == 1]
        for _ in range(100):
            # Any product of three reflections in S3 is a reflection, so these
            # are all-q=1 meshes without assuming a homogeneous link pattern.
            raw = [rng.choice(reflections) for _ in x.e.geometry.edges]
            row = homogeneous_state(x, raw)
            self.assertEqual(row['candidate_frozen'], row['holonomy_group_order'] == 2)
        raw = constant_connection(x.p, (1, 1, 1))
        for _ in range(50):
            frames = [rng.randrange(6) for _ in range(x.e.geometry.size)]
            moved = [g.mul[frames[v]][g.mul[z][g.inv[frames[u]]]] for (u, v), z in zip(x.e.geometry.edges, raw)]
            row = homogeneous_state(x, moved)
            self.assertTrue(row['candidate_frozen'])
            self.assertEqual(row['parallel_reflection_section'], [g.conj[h][1] for h in frames])
        with self.assertRaises(ValueError):
            homogeneous_state(x, [0]*len(raw))

    def test_collision_is_invisible_to_complete_pair_orbit_but_not_global_gauge(self):
        words = WordObserver(self.g, self.q)
        for t in self.census['elastic_tables']:
            for code, target in enumerate(t):
                self.assertEqual(words.canonical(divmod(code, 6)), words.canonical(divmod(target, 6)))
        w = self.data['encounter_witness']
        a, b = w['initial_links'], w['scattered_links']
        self.assertEqual(sum(x != y for x, y in zip(a, b)), 1)
        self.assertEqual(self.x.p.charge(a), self.x.p.charge(b))
        self.assertNotEqual(words.canonical(self.x.e.forest.based_loops(a)[0]),
                            words.canonical(self.x.e.forest.based_loops(b)[0]))
        self.assertNotEqual(activity_mask(self.x.geometry, a), activity_mask(self.x.geometry, b))
        self.assertEqual(Fraction(w['one_step_charge_law_l1_numerator'], 2*w['readout_operators']), Fraction(2, 117))

    def test_all_five_link_assignments_have_exact_inverse_and_unchanged_spectators(self):
        x, oracle = self.x, self.x.e.factor.oracle
        support = sorted({edge for path in oracle.pairs[0] for edge, _ in path})
        self.assertEqual(len(support), 5)
        write = oracle.pairs[0][0][0][0]
        for values in itertools.product(range(6), repeat=5):
            raw = [i % 6 for i in range(len(x.e.geometry.edges))]
            for edge, value in zip(support, values):
                raw[edge] = value
            before, q = raw[:], x.p.charge(raw)
            x.apply(raw, x.supports)
            self.assertEqual(x.p.charge(raw), q)
            self.assertTrue(all(a == b for i, (a, b) in enumerate(zip(raw, before)) if i != write))
            x.apply(raw, 2*x.supports)
            self.assertEqual(raw, before)

    def test_actual_trajectory_is_covariant_under_nonconstant_local_frames(self):
        x, g, rng = self.x, self.g, random.Random(60821)
        raw = [rng.randrange(6) for _ in x.e.geometry.edges]
        frames = [rng.randrange(6) for _ in range(x.e.geometry.size)]
        transform = lambda links: [g.mul[frames[v]][g.mul[z][g.inv[frames[u]]]]
                                   for (u, v), z in zip(x.e.geometry.edges, links)]
        gauged, initial = transform(raw), raw[:]
        schedule = [rng.randrange(x.operator_count) for _ in range(500)]
        for op in schedule:
            x.apply(raw, op); x.apply(gauged, op)
            self.assertEqual(transform(raw), gauged)
        for op in reversed(schedule):
            x.apply(raw, op, reverse=True)
        self.assertEqual(raw, initial)
        runs = x.compiled([initial, transform(initial)], schedule, stride=100)['runs']
        for a, b in zip(runs[:2], runs[2:]):
            self.assertEqual(transform(a['final_links']), b['final_links'])

    def test_two_charge_sector_still_has_the_exact_mean_heat_law(self):
        x, reference, rng = self.x, AxialReference(3, self.bank, 2), random.Random(89536)
        for _ in range(12):
            initial = reference.sample(rng, 'nonabelian_reflections')['links']
            q, drift = np.array(x.p.charge(initial)), np.zeros(x.geometry.faces, dtype=np.int64)
            for op in range(x.operator_count):
                raw = initial[:]; x.apply(raw, op)
                drift += np.array(x.p.charge(raw))-q
            np.testing.assert_array_equal(drift, -2*x.geometry.l1@q)

    def test_patch_kernel_unfreezes_two_states_with_exact_geometric_escape(self):
        k = self.data['patch_activation']
        self.assertEqual((k['gauge_states'], k['baseline_operators'], k['candidate_operators']), (49, 16, 24))
        self.assertEqual((len(k['baseline_components']), len(k['candidate_components'])), (30, 28))
        self.assertEqual(len(k['component_mergers']), 1)
        merge = k['component_mergers'][0]
        self.assertEqual(sorted(map(len, merge['previous_components'])), [1, 1, 8])
        self.assertEqual(Fraction(*merge['escape_probability']), Fraction(1, 6))
        self.assertEqual(Fraction(*merge['mean_first_escape_attempts']), 6)
        frozen, tables = merge['previously_frozen_states'], k['transition_tables']
        counts = np.array([[sum(t[i] == j for t in tables) for j in frozen] for i in frozen], dtype=object)
        survival = np.ones(len(frozen), dtype=object)
        for n in range(1, 9):
            survival = counts@survival
            self.assertEqual(list(survival), [20**n]*len(frozen))
        all_counts = np.array([[sum(t[i] == j for t in tables) for j in range(49)] for i in range(49)])
        np.testing.assert_array_equal(all_counts, all_counts.T)

    def test_exact_matrix_powers_check_matched_clock_response_identity_on_every_patch_state(self):
        k, union = self.data['patch_activation'], PatchUnion(self.x.e, (0,))
        tables = k['transition_tables']
        original = tables[:16]+[list(range(49)) for _ in range(8)]
        counts = lambda rows: np.array([[sum(t[i] == j for t in rows) for j in range(49)] for i in range(49)], dtype=object)
        a, b = counts(original), counts(tables)
        q = np.array([[self.x.p.charge(union.lift(rep))[f] for f in self.x.e.fan_faces[0]] for rep in union.reps], dtype=object)
        np.testing.assert_array_equal(a@q, b@q)
        np.testing.assert_array_equal(a@a@q, b@b@q)
        delta = b@b@b@q-a@a@a@q
        l0 = a-24*np.eye(49, dtype=int)
        np.testing.assert_array_equal(delta, (b-a)@l0@l0@q)
        # The isolated patch is a NEGATIVE control for mean response. Its
        # fluctuations still change; neighboring mesh dynamics enables the
        # nonzero whole-mesh third-step mean witness below.
        self.assertFalse(np.any(delta))
        values = q
        for _ in range(49):
            np.testing.assert_array_equal((b-a)@values, np.zeros((49, 3), dtype=int))
            values = a@values
        square_delta = b@b@(q*q)-a@a@(q*q)
        self.assertTrue(np.any(square_delta))
        np.testing.assert_array_equal(square_delta, k['two_step_square_charge_difference_numerator'])
        whole = self.data['automatic_matched_clock_response']
        self.assertEqual(whole['attempted_operators'], 810)
        self.assertEqual(whole['third_step_mean_denominator'], 810**3)
        self.assertEqual(sum(whole['third_step_mean_difference_numerator']), 0)
        self.assertTrue(any(whole['third_step_mean_difference_numerator']))

    def test_cpp_protocol_rejects_invalid_pair_banks_and_preserves_legacy_dispatch(self):
        cmd = ['build/wgphysics_mixed_bank_experiments', '--cycle', '3', '--pair-bank']
        for count in (0, 33, -1):
            run = subprocess.run(cmd, input=f'3 12 1 1 1 {count}\n', capture_output=True, text=True)
            self.assertNotEqual(run.returncode, 0)
            self.assertIn('pair bank size', run.stderr)
        # Legacy Python wrapper also independently verifies the old protocol,
        # original rule indices and exact reversal against its raw oracle.
        raw = [i % 6 for i in range(len(self.x.e.geometry.edges))]
        self.x.e.compiled_bank([raw], list(range(702)), 702, capture_links=True)
        for bad in (-1, self.x.operator_count, True, 1.5):
            with self.assertRaises(ValueError):
                self.x.apply(raw, bad)


if __name__ == '__main__':
    unittest.main()
