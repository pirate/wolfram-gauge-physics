import itertools
import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'tools'))
from face_energy_obstruction import derive_fiber_group
from triple_rule_search import algebra, assess, make_table
from three_face_reachability import (AbelianFaceFactor, compiled_replay, generated_subgroup,
                                    link_pair, random_schedule, replay, search, sector_count,
                                    verify_local_factor)


class ReachabilityTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.census = json.loads((ROOT/'data/d4-triple-rule-search.json').read_text())
        cls.g = derive_fiber_group()
        cls.collision = cls.census['selected_rule']['table']
        cls.table = cls.census['combined_rule']['table']
        cls.subgroup = generated_subgroup(cls.g, [1, 4])
        cls.factor = AbelianFaceFactor(3, cls.g, cls.table, cls.subgroup)

    def test_factor_exhausts_all_seven_link_assignments(self):
        self.assertEqual(self.subgroup, (0, 1, 4, 5))
        self.assertEqual(verify_local_factor(self.factor), 4**7)
        with self.assertRaisesRegex(ValueError, 'abelian'):
            AbelianFaceFactor(3, self.g, self.table, tuple(range(8)))
        with self.assertRaisesRegex(ValueError, 'closed subgroup'):
            AbelianFaceFactor(3, self.g, self.table, (0, 1, 4))
        links = link_pair(self.factor)
        links[0] = 3
        with self.assertRaisesRegex(ValueError, 'outside'):
            self.factor.project(links)

    def test_every_minimal_rule_preserves_exact_frame_stabilizers(self):
        states, _, action, _ = algebra(self.g)
        stabilizers = [tuple(frame for frame, row in enumerate(action) if row[x] == x)
                       for x in range(len(states))]
        for rule in self.census['minimal_rules']:
            table = make_table(rule['transpositions'])
            self.assertTrue(all(stabilizers[x] == stabilizers[y] for x, y in enumerate(table)))
        # An isolated reflection-pair-to-central fusion loses the frame stabilizer.
        self.assertEqual(stabilizers[12], (0, 1, 4, 5))  # (1,r,r')
        self.assertEqual(stabilizers[5], tuple(range(8)))  # (1,1,z)
        with self.assertRaisesRegex(ValueError, 'gauge covariance'):
            # Include the orientation-reversed pair so covariance is isolated.
            assess(self.g, make_table([(5, 12), (264, 320)]))
        # The selected catalyst retains that stabilizer on both sides.
        self.assertEqual(stabilizers[76], stabilizers[13])
        self.assertEqual(stabilizers[76], (0, 1, 4, 5))

    def test_two_exhaustive_components_and_exact_stationary_counts(self):
        f = self.factor
        initial = f.project(link_pair(f))
        counts = sector_count(f)
        result = search(f, initial, self.collision, exhaustive=True)
        self.assertTrue(result['complete'])
        self.assertEqual(result['states_discovered'], 6480)
        self.assertEqual(result['central_face_histogram'], counts['reflection_present_central_histogram'])
        self.assertEqual(result['directed_collision_events'], counts['reflection_present_directed_collision_events'])
        self.assertGreater(result['minimum_idle_patches_per_state'], 0)
        central_links = [0]*len(f.oracle.edges)
        central_links[0] = 5
        central = search(f, f.project(central_links), self.collision, exhaustive=True)
        self.assertEqual(central['states_discovered'], 81)
        self.assertIsNone(central['shortest_collision_path'])
        self.assertEqual(central['directed_collision_events'], 0)
        self.assertEqual(result['states_discovered']+central['states_discovered'], counts['all_charge_and_product_compatible_states'])
        truncated = search(f, initial, self.collision, max_states=1, exhaustive=True)
        self.assertFalse(truncated['complete'])
        self.assertEqual(truncated['stop_reason'], 'state budget')

    def test_shortest_encounter_replays_and_is_frame_covariant(self):
        f = self.factor
        initial = link_pair(f)
        path = search(f, f.project(initial), self.collision)['shortest_collision_path']
        self.assertEqual(len(path), 3)
        expected = replay(f, initial[:], path, self.collision)
        self.assertEqual([r['collision'] for r in expected['steps']], [False, False, True])
        self.assertTrue(compiled_replay(f, initial, path, expected))
        frames = [(v*3+1) % 8 for v in range(9)]
        def change(values):
            return [self.g.mul[frames[v]][self.g.mul[a][self.g.inv[frames[u]]]]
                    for a, (u, v) in zip(values, f.oracle.edges)]
        inputs = [3, *self.table, *change(initial), len(path), *path]
        process = subprocess.run([str(ROOT/'build/wgphysics_three_face_replay')],
                                 input=' '.join(map(str, inputs))+'\n', text=True, capture_output=True, check=True)
        actual = json.loads(process.stdout)
        self.assertEqual(actual['final_links'], change(expected['final_links']))
        self.assertTrue(actual['exact_link_inverse'])

    def test_dual_tree_lift_and_sparse_frontier(self):
        f = self.factor
        # Every formal state in this charge/product sector has a link preimage.
        choices = []
        for part in (0, 1):
            faces = [f for f, color in enumerate(self.factor.oracle.colors) if color == part]
            choices.append([[(a, 5)] for a in faces]
                           + [[(a, 1), (b, 4)] for a in faces for b in faces if a != b])
        checked = 0
        for left, right in itertools.product(*choices):
            state = bytearray(18)
            for a, value in left+right:
                state[a] = value
            self.assertEqual(f.project(f.lift_faces(state)), state)
            checked += 1
        self.assertEqual(checked, 6561)
        state = bytearray(18)
        state[0] = 1
        with self.assertRaisesRegex(ValueError, 'global product'):
            f.lift_faces(state)
        initial = f.project(link_pair(f))
        candidates = f.candidate_patches(initial)
        self.assertEqual([p for p in candidates if f.step(initial, p) != initial],
                         [p for p in range(len(f.patches)) if f.step(initial, p) != initial])

    def test_random_schedule_counts_attempts_and_replays(self):
        control = AbelianFaceFactor(3, self.g, self.census['transport_control']['table'], self.subgroup)
        run = random_schedule(self.factor, link_pair(self.factor), self.collision, 1000, 5920941, control)
        self.assertEqual(sum(run['central_face_residence_attempts'].values()), 1000)
        self.assertEqual(run['changing_events'], len(run['changing_patch_path']))
        self.assertTrue(run['compiled_raw_link_replay'])
        self.assertGreater(run['collision_events'], 0)
        self.assertTrue(run['transport_control']['compiled_raw_link_replay'])
        self.assertEqual(run, random_schedule(self.factor, link_pair(self.factor), self.collision, 1000, 5920941, control))


if __name__ == '__main__':
    unittest.main()
