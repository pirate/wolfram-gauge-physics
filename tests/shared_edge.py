import copy
import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'tools'))
from run_shared_edge import analyze, bipartition, class_memory_model
from face_energy_obstruction import derive_fiber_group


class SharedEdgeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        census = json.loads((ROOT/'data/d4-involution-braid-search.json').read_text())
        cls.table = census['solutions'][11229]['table']
        process = subprocess.run([str(ROOT/'build/wgphysics_mesh_experiments'), '--lift', 'shared-edge',
                                  '--side', '3', '--layers', '12', '--trials', '1'],
                                 input=' '.join(map(str, cls.table))+'\n', text=True,
                                 capture_output=True, check=True)
        cls.raw = json.loads(process.stdout)

    def test_complete_link_replay_and_derived_currents(self):
        result = analyze(copy.deepcopy(self.raw), self.table, full_replay=True)
        check = result['analysis']
        self.assertEqual(len(check['additive_charge_basis']), 3)
        self.assertEqual(check['forward_updates'], check['independent_replayed_updates'])
        self.assertGreater(check['forward_updates'], 0)
        memory = check['derived_class_memory']
        self.assertEqual(memory['vacancy_sectors'], [0, 5])
        self.assertEqual(memory['flip_species'], [1])
        self.assertEqual(memory['local_staggered_parity_checks'], 50)
        self.assertTrue(any(r['independent_transport']['charge_hops'] for r in result['runs']))
        self.assertTrue(any(r['independent_transport']['vacancy_flip_hops'] for r in result['runs']))
        for run in result['runs']:
            self.assertIn(run['independent_transport']['conserved_staggered_parity'], (0, 1))
            self.assertLessEqual(run['independent_transport']['vacancy_flip_hops'],
                                 run['independent_transport']['charge_hops'])

    def test_wrong_lift_or_missing_spectators_rejected(self):
        result = copy.deepcopy(self.raw)
        result['lift'] = 'exclusive-closing'
        with self.assertRaisesRegex(ValueError, 'explicitly selected'):
            analyze(result, self.table)
        result = copy.deepcopy(self.raw)
        result['local_census']['transitions'][0][6] = 1
        with self.assertRaisesRegex(ValueError, 'spectator'):
            analyze(result, self.table)

    def test_corrupted_schedule_and_links_rejected_by_independent_replay(self):
        result = copy.deepcopy(self.raw)
        result['runs'][0]['schedule'][0].pop()
        with self.assertRaisesRegex(ValueError, 'drops or duplicates'):
            analyze(result, self.table, True)
        result = copy.deepcopy(self.raw)
        result['runs'][0]['final_links'][0] = 1
        with self.assertRaisesRegex(ValueError, 'final links'):
            analyze(result, self.table, True)

    def test_charge_and_closed_surface_constraint_reject_corrupted_histogram(self):
        result = copy.deepcopy(self.raw)
        row = result['runs'][0]['trajectory'][1]
        row[5] -= 1; row[6] += 1; row[1] += 1
        with self.assertRaisesRegex(ValueError, 'neutrality|charge drifted'):
            analyze(result, self.table)

    def test_bipartite_scope_is_explicit(self):
        self.assertEqual(bipartition(4, [(0, 1), (1, 2), (2, 3), (3, 0)]), [0, 1, 0, 1])
        with self.assertRaisesRegex(ValueError, 'not bipartite'):
            bipartition(3, [(0, 1), (1, 2), (2, 0)])

    def test_memory_factor_is_not_assumed_for_arbitrary_tables(self):
        group = derive_fiber_group()
        charges = [[int(s == label) for s in group.sectors] for label in (1, 2, 3)]
        self.assertIsNone(class_memory_model(group, list(range(64)), charges))


if __name__ == '__main__':
    unittest.main()
