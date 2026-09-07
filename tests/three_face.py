import copy
import itertools
import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'tools'))
from face_energy_obstruction import derive_fiber_group
from triple_rule_search import analyze as analyze_census, assess, make_table, minimal_closures
from run_three_face import analyze, LinkOracle


class ThreeFaceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.census = json.loads((ROOT/'data/d4-triple-rule-search.json').read_text())
        inputs = [str(x) for name in ('selected_rule', 'transport_control', 'combined_rule') for x in cls.census[name]['table']]
        process = subprocess.run([str(ROOT/'build/wgphysics_three_face_experiments'), '--side', '4', '--layers', '24', '--trials', '1'],
                                 input=' '.join(inputs)+'\n', capture_output=True, text=True, check=True)
        cls.raw = json.loads(process.stdout)

    def test_exhaustive_minimal_closure_coverage(self):
        group = derive_fiber_group()
        expected = minimal_closures(group)
        actual = [tuple(map(tuple, r['transpositions'])) for r in self.census['minimal_rules']]
        self.assertEqual(actual, expected)
        self.assertEqual(len(actual), 945)
        self.assertEqual(sum(r['same_classes_and_boundary_feedback'] is not None for r in self.census['minimal_rules']), 354)
        raw = {'group_order': 8, 'minimal_rules': [{'id': r['id'], 'transpositions': r['transpositions']} for r in self.census['minimal_rules']]}
        raw['minimal_rules'].pop()
        with self.assertRaisesRegex(ValueError, 'independent Python action'):
            analyze_census(raw)

    def test_combination_rederived_without_an_inserted_charge(self):
        group = derive_fiber_group()
        collision, transport, combined = [self.census[k]['table'] for k in ('selected_rule', 'transport_control', 'combined_rule')]
        self.assertEqual(combined, [transport[collision[i]] for i in range(512)])
        self.assertFalse(any(collision[i] != i and transport[i] != i for i in range(512)))
        derived = assess(group, combined)
        self.assertEqual(derived['moved_states'], 22)
        self.assertEqual(derived['charge_basis'], [[0, 1, 0, 0], [0, 0, 1, 0], [1, 0, 0, 2]])
        witness = derived['same_classes_and_boundary_feedback']
        self.assertIsNotNone(witness)
        self.assertNotEqual(*witness['output_charge_vectors'])
        with self.assertRaisesRegex(ValueError, 'overlap'):
            make_table([(1, 2), (2, 3)])

    def test_full_independent_link_replay_and_fixed_boundary_witness(self):
        result = analyze(copy.deepcopy(self.raw), self.census, True)
        a = result['analysis']
        self.assertEqual(a['independent_replayed_updates'], a['forward_updates'])
        self.assertGreater(a['forward_updates'], 0)
        self.assertEqual(len(a['bipartite_charge_basis']), 6)
        self.assertEqual(a['positive_weight_per_element'], [0, 1, 1, 1, 1, 2, 1, 1])
        witness = a['fixed_exterior_witness']
        self.assertFalse(witness['gauge_equivalent'])
        oracle = LinkOracle(result['side'], derive_fiber_group())
        self.assertEqual(oracle.sectors(witness['initial_left']), oracle.sectors(witness['initial_right']))
        self.assertNotEqual(oracle.sectors(witness['after_left']), oracle.sectors(witness['after_right']))
        self.assertTrue(all(x == y for i, (x, y) in enumerate(zip(witness['initial_left'], witness['initial_right']))
                            if i not in oracle.writes[witness['prefix_patch']]))

    def test_proven_no_reaction_subgroup_controls(self):
        collision = self.census['selected_rule']['table']
        combined = self.census['combined_rule']['table']
        for subgroup in ((0, 1), (0, 3, 5, 6)):
            for a, b, c in itertools.product(subgroup, repeat=3):
                i = (a*8+b)*8+c
                self.assertEqual(collision[i], i)
                j = combined[i]
                self.assertTrue(all(x in subgroup for x in (j//64, (j//8) % 8, j % 8)))
        lookup = {(r['condition'], r['mode']): r for r in self.raw['runs']}
        for condition in ('flat', 'reflection_link', 'rotation_link', 'reflection_pair'):
            a, b = lookup[condition, 'combined'], lookup[condition, 'transport']
            self.assertEqual(a['final_links'], b['final_links'])
            self.assertEqual(sum(row[1] for row in a['trajectory']), 0)

    def test_corrupt_controls_schedule_or_event_counts_rejected(self):
        raw = copy.deepcopy(self.raw)
        raw['schedules'][0][0].pop()
        with self.assertRaisesRegex(ValueError, 'drops or duplicates'):
            analyze(raw, self.census, True)
        raw = copy.deepcopy(self.raw)
        raw['runs'][0]['initial_links'][0] = 1
        with self.assertRaisesRegex(ValueError, 'share their initial'):
            analyze(raw, self.census, True)
        raw = copy.deepcopy(self.raw)
        raw['runs'][0]['trajectory'][1][1] = 1
        with self.assertRaisesRegex(ValueError, 'event classification'):
            analyze(raw, self.census, True)

    def test_gauge_oracle_accepts_actual_frame_copies(self):
        group = derive_fiber_group()
        oracle = LinkOracle(self.raw['side'], group)
        values = next(r['initial_links'] for r in self.raw['runs'] if r['condition'] == 'witness_left')
        frames = [(v*3+1) % 8 for v in range(self.raw['side']**2)]
        changed = [group.mul[frames[v]][group.mul[a][group.inv[frames[u]]]] for a, (u, v) in zip(values, oracle.edges)]
        self.assertTrue(oracle.gauge_equivalent(values, changed))


if __name__ == '__main__':
    unittest.main()
