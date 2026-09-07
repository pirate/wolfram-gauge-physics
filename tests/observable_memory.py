import copy
import itertools
import json
import sys
import unittest
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'tools'))
from mesh_observable_audit import analyze_orbits, memory_audit, verify_link_records


class ObservableMemoryTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.mesh = json.loads((ROOT/'data/d4-shared-mesh.json').read_text())
        cls.data = json.loads((ROOT/'data/d4-mesh-observable-audit.json').read_text())
        census = json.loads((ROOT/'data/d4-involution-braid-search.json').read_text())
        cls.table = census['solutions'][cls.mesh['solution_id']]['table']

    def verify(self, data):
        return verify_link_records(data, self.table, self.mesh['group_permutations'])

    def test_independent_paths_gauge_quotient_and_transition_features(self):
        result = self.verify(self.data)
        self.assertEqual(result['tetrahedron_successors'], 8448)
        self.assertEqual(result['octahedron_feature_transitions'], 2048)
        self.assertTrue(result['central_type_frozen_gate_certificate'])

    def test_corrupted_records_are_rejected(self):
        mutations = [
            (lambda d: d['tetrahedron']['states'][0]['successors'].__setitem__(0, 1), 'successor'),
            (lambda d: d['tetrahedron']['states'][0]['face_sectors'].__setitem__(0, 1), 'face sectors'),
            (lambda d: d['tetrahedron']['states'][0]['successors'].pop(), 'generator count'),
            (lambda d: d['octahedron']['checks'][0]['after_features'].__setitem__(0, -1), 'feature vector'),
            (lambda d: d['octahedron']['checks'][1]['before_links'].__setitem__(0, -1), 'continuity'),
        ]
        for mutate, error in mutations:
            with self.subTest(error=error):
                data = copy.deepcopy(self.data)
                mutate(data)
                with self.assertRaisesRegex(ValueError, error):
                    self.verify(data)

    def test_exact_two_step_memory_error(self):
        result = memory_audit(self.mesh['local_census'])
        self.assertEqual(result['naive_markov_stationary_average_two_step_return'], [467, 512])
        self.assertEqual(result['naive_markov_minimum_two_step_return'], [1, 4])
        self.assertEqual(result['coarse_states_with_false_two_step_spreading'], 192)
        self.assertEqual(result['one_generator_memory_refinement_states'], 757)
        total = sum(row['rooted_multiplicity'] for row in result['rows'])
        weighted = sum(row['rooted_multiplicity']*Fraction(*row['markov_two_step_return'])
                       for row in result['rows'])/total
        self.assertEqual(weighted, Fraction(467, 512))

    def test_closed_coarse_control_has_no_false_spreading(self):
        # A nontrivial involution on four variables, with an autonomous coarse map.
        states = itertools.product(range(3), repeat=4)
        local = {'element_sectors': [0, 1, 1],
                 'transitions': [list(s)+[s[1], s[0], s[3], s[2]] for s in states]}
        result = memory_audit(local)
        self.assertEqual(result['naive_markov_minimum_two_step_return'], [1, 1])
        self.assertEqual(result['coarse_states_with_false_two_step_spreading'], 0)
        self.assertEqual(result['coarse_states'], result['one_generator_memory_refinement_states'])
        local['transitions'].append(local['transitions'][0])
        with self.assertRaisesRegex(ValueError, 'exhaustive involution'):
            memory_audit(local)

    def test_components_frozen_states_and_schedule_are_separate(self):
        data = copy.deepcopy(self.data)
        data['independent_verification'] = self.verify(data)
        result = analyze_orbits(data)['analysis']
        self.assertEqual(result['component_size_counts'], [(1, 56), (6, 2), (12, 1), (24, 2), (48, 1)])
        self.assertEqual(result['frozen_gauge_states'], {'all_central': 8, 'all_noncentral': 48})
        self.assertEqual(result['predictive_face_partition']['refinement_sizes'], [149])
        self.assertEqual(result['periodic_sweeps']['index_order']['cycle_counts'], [[1, 176]])
        self.assertNotEqual(result['periodic_sweeps']['permuted_order']['cycle_counts'], [[1, 176]])

    def test_pair_invariant_falsification_and_global_constant_certificate(self):
        result = analyze_orbits(copy.deepcopy(self.data))['analysis']
        candidates = result['tetrahedron_nonconstant_observables']
        self.assertEqual(len(candidates), 12)
        self.assertEqual(result['joint_nonconstant_observables'], [])
        for candidate in candidates:
            self.assertIsNotNone(candidate['first_octahedron_violation'])
            self.assertNotEqual(*candidate['violation_values'])
        reverse = result['reverse_pair_classes']
        self.assertEqual([reverse[j] for j in reverse], list(range(28)))
        self.assertEqual(result['joint_weight_nullity'], 11)
        self.assertEqual(sum(i != j for i, j in enumerate(reverse))//2, 10)
        for weight in result['joint_weights']:
            self.assertEqual(len({weight[i]+weight[j] for i, j in enumerate(reverse)}), 1)


if __name__ == '__main__':
    unittest.main()
