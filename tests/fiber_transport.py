import itertools
import json
import random
import subprocess
import sys
import unittest
from pathlib import Path

import numpy as np
from scipy.linalg import eigh

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'tools'))
from d4_loop_observer import ConnectionForest
from face_energy_obstruction import FiniteGroup
from fiber_mode_probe import lifted_graph_laplacian
from fiber_transport import FiberTransport, triangle_projective_audit
from transport_braid import pure_word, spectral_word_witness


class FiberTransportTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data = json.loads((ROOT/'data/fiber-transport-order.json').read_text())
        cls.experiments = [FiberTransport(cls.data['side'], r['fiber_vertices'], r['fiber_edges']) for r in cls.data['fibers']]

    @staticmethod
    def transformed(e, links, frames):
        g = e.geometry.group
        return [g.mul[frames[v]][g.mul[x][g.inv[frames[u]]]] for (u, v), x in zip(e.geometry.edges, links)]

    def test_generic_driver_derives_groups_and_replays_multiple_fibers(self):
        rng = random.Random(613902)
        examples = [(1, [], 1), (3, [(0, 1), (1, 2)], 2),
                    (3, [(0, 1), (1, 2), (2, 0)], 6),
                    (4, [(i, (i+1) % 4) for i in range(4)], 8),
                    (5, [(i, (i+1) % 5) for i in range(5)], 10),
                    (4, list(itertools.combinations(range(4), 2)), 24)]
        for size, edges, order in examples:
            e = FiberTransport(3, size, edges)
            g = e.geometry.group
            self.assertEqual(g.n, order)
            links = [rng.randrange(order) for _ in e.geometry.edges]
            frames = [rng.randrange(order) for _ in range(e.geometry.size)]
            frames[0] = order-1
            moved = self.transformed(e, links, frames)
            witness = e.frame_witness(links, moved)
            self.assertIsNotNone(witness)
            self.assertEqual(self.transformed(e, links, witness), moved)
            schedule = [rng.randrange(len(e.factor.pairs)) for _ in range(48)]
            a, b = e.compiled(links, schedule), e.compiled(moved, schedule)
            self.assertEqual(self.transformed(e, a['final_links'], frames), b['final_links'])
            self.assertEqual(a['final_signature'], b['final_signature'])
            self.assertEqual(e.compiled([g.identity]*len(links), schedule)['events'], [])

    def test_driver_rejects_invalid_protocol_before_producing_a_result(self):
        e = FiberTransport(3, 3, [(0, 1), (1, 2), (2, 0)])
        base = [3, 3, 3, 0, 0, 1, 0, 2, 1, 2]+[0]*len(e.geometry.edges)
        cases = [base[:-1], base+['unexpected']]
        for index, value in ((0, 25), (1, 7), (3, 100001), (5, 0), (9, 3), (10, 6)):
            changed = base[:]; changed[index] = value; cases.append(changed)
        changed = base[:]; changed[3] = 1; changed.append(len(e.factor.pairs)); cases.append(changed)
        cases.append([3, 6, 0, 0])
        for values in cases:
            result = subprocess.run(['build/wgphysics_fiber_transport'], input=' '.join(map(str, values))+'\n',
                                    capture_output=True, text=True, timeout=10)
            self.assertNotEqual(result.returncode, 0)
            self.assertEqual(result.stdout, '')
        with self.assertRaisesRegex(ValueError, 'schedule'):
            e.compiled([0]*len(e.geometry.edges), [-1])

    def test_complete_forest_is_independent_of_group_element_enumeration(self):
        e = self.experiments[0]
        g = e.geometry.group
        rng = random.Random(893051)
        links = [rng.randrange(g.n) for _ in e.geometry.edges]
        representative = e.forest.representative(e.forest.signature(links))
        self.assertTrue(e.factor.oracle.fan.gauge_equivalent(links, representative))
        reverse = FiniteGroup(list(reversed(g.elements)))
        relabeled = [reverse.elements.index(g.elements[x]) for x in links]
        forest = ConnectionForest(range(e.geometry.size), e.geometry.edges, reverse)
        self.assertEqual(forest.legacy_signature(relabeled), e.forest.legacy_signature(links))
        disconnected = ConnectionForest([0, 1, 2, 9], [(0, 1), (0, 2), (1, 2)], g)
        signature = disconnected.signature([1, 2, 5])
        self.assertEqual(disconnected.signature(disconnected.representative(signature)), signature)
        self.assertEqual([len(x) for x in signature], [1, 0])

    def test_saved_paths_are_identical_across_fibers_and_commuting_control(self):
        records = self.data['fibers']
        for e, r in zip(self.experiments, records):
            prep = r['preparation']
            self.assertEqual(e.compiled(prep['initial_links'], [x[2] for x in prep['events']]), prep['checked'])
            self.assertEqual(prep['paths'], records[0]['preparation']['paths'])
            self.assertEqual([c['schedule'] for c in r['circuits']], [c['schedule'] for c in records[0]['circuits']])
            for circuit in r['circuits']:
                self.assertEqual(len(circuit['schedule']), 32)
                for face, winding in circuit['windings'].items():
                    self.assertEqual(e.winding(circuit['walk'], int(face)), winding)

    def test_order_dependence_and_relative_loop_witness_on_full_connections(self):
        for index, (e, record) in enumerate(zip(self.experiments, self.data['fibers'])):
            finals = []
            for order in record['orders']:
                result = e.compiled(record['preparation']['links'], order['schedule'])
                self.assertEqual(result, order['checked'])
                finals.append(result['final_links'])
                self.assertEqual(e.geometry.sectors(result['final_links']), e.geometry.sectors(record['preparation']['links']))
            self.assertEqual(e.frame_witness(*finals) is not None, index != 0)
            self.assertEqual(e.factor.oracle.fan.gauge_equivalent(*finals), index != 0)
        witness = next(w for w in self.data['fibers'][0]['relative_loop_words'] if w['positions'] == [0, 2])
        self.assertEqual(witness['classes_by_order'], [3, 0])

    def test_exhausted_controlled_components_realize_a4_z2_and_trivial_actions(self):
        for e, record, size, order in zip(self.experiments, self.data['fibers'], (4, 2, 1), (12, 2, 1)):
            orbit = record['controlled_circuit_orbit']
            states = orbit['raw_representatives']
            self.assertTrue(orbit['exhausted'])
            self.assertEqual(len(states), size)
            self.assertEqual(len({e.forest.signature(s) for s in states}), size)
            self.assertEqual(orbit['permutation_group_order'], order)
            for source, row in enumerate(orbit['action_targets_by_state']):
                for action, target in enumerate(row):
                    checked = e.compiled(states[source], record['circuits'][action]['schedule'])
                    saved = orbit['edge_witnesses'][source][action]
                    self.assertEqual(checked, saved['checked'])
                    self.assertEqual(self.transformed(e, checked['final_links'], saved['frames_to_target_representative']), states[target])
            permutations = [tuple(p) for p in orbit['permutation_group']]
            compose = lambda a, b: tuple(a[b[i]] for i in range(size))
            self.assertTrue(all(compose(a, b) in permutations for a in permutations for b in permutations))
        triangle = self.data['fibers'][0]['controlled_circuit_orbit']
        self.assertTrue(triangle['equals_alternating_group_on_four_orbits'])
        self.assertEqual(triangle['generator_orders'], [3, 3])
        self.assertEqual(triangle['commutator_order'], 2)
        with self.assertRaisesRegex(ValueError, 'closure not proved'):
            self.experiments[0].orbit_census(self.data['fibers'][0]['preparation']['links'], self.data['fibers'][0]['circuits'], maximum=1)

    def test_deformed_programs_preserve_every_measured_orbit_action(self):
        for e, r in zip(self.experiments, self.data['fibers']):
            orbit = r['controlled_circuit_orbit']
            self.assertEqual(len(orbit['deformed_action_checks']), 2*len(orbit['raw_representatives']))
            for d in orbit['deformed_action_checks']:
                spec = orbit['chosen_winding_detours'][d['action']]
                original = orbit['edge_witnesses'][d['source_state']][d['action']]['checked']['final_links']
                self.assertEqual({i for i, x in enumerate(d['frames_from_reference_action']) if x}, {spec['primal_vertex']})
                self.assertEqual(self.transformed(e, original, d['frames_from_reference_action']), d['checked']['final_links'])
                self.assertEqual(e.forest.signature(d['checked']['final_links']),
                                 e.forest.signature(orbit['raw_representatives'][d['target_state']]))

    def test_inverse_pure_words_are_exact_inverses(self):
        for e in self.experiments[:2]:
            g = e.geometry.group
            for values in itertools.product(range(g.n), repeat=3):
                for i, j in itertools.combinations(range(3), 2):
                    moved = pure_word(g, values, i, j)
                    self.assertEqual(pure_word(g, moved, i, j, inverse=True), values)

    def test_fixed_based_loop_coordinates_intertwine_every_raw_circuit(self):
        for e, record in zip(self.experiments, self.data['fibers']):
            audit = e.loop_action_audit(record['preparation']['positions'], record['controlled_circuit_orbit'])
            self.assertEqual(json.loads(json.dumps(audit)), record['loop_action_audit'])
            self.assertEqual(audit['position_order'], (0, 1, 3, 2))
            self.assertEqual(audit['faces'], [52, 62, 208, 172])
            self.assertTrue(audit['neutral'])
            self.assertTrue(audit['injective_on_controlled_component'])
            for check in audit['checks']:
                self.assertEqual(check['actual_tuple'], check['predicted_tuple'])

    def test_neutral_triangle_tuples_are_a_projective_line_plus_commuting_point(self):
        e = self.experiments[0]
        audit = triangle_projective_audit(e.geometry.group)
        self.assertEqual(len(audit['neutral_raw_tuples']), 27)
        self.assertEqual(audit['neutral_gauge_orbits'], 5)
        self.assertEqual(audit['matrix_group_order'], 24)
        self.assertEqual(audit['projective_action_order'], 12)
        self.assertEqual(len(audit['projective_kernel']), 2)
        self.assertTrue(audit['equals_A4'])
        classes = {}
        for row in audit['neutral_raw_tuples']:
            classes.setdefault(row['projective_point'], []).append(row['tuple'])
        self.assertEqual(sorted(len(t) for t in classes.values()), [3, 6, 6, 6, 6])
        saved = self.data['fibers'][0]['loop_action_audit']['projective_points_by_state']
        self.assertEqual(set(map(tuple, saved)), set(audit['projective_points']))
        self.assertEqual(self.data['fibers'][2]['loop_action_audit']['projective_points_by_state'], [[0, 0]])
        # Relabeling group indices must not change the affine or projective result.
        reverse = triangle_projective_audit(FiniteGroup(list(reversed(e.geometry.group.elements))))
        for key in ('inverse_pure_matrices', 'projective_points', 'projective_generator_permutations', 'projective_kernel'):
            self.assertEqual(audit[key], reverse[key])
        with self.assertRaisesRegex(ValueError, 'triangle'):
            triangle_projective_audit(self.experiments[1].geometry.group)

    def test_full_bundle_construction_and_exact_order_dependent_spectrum(self):
        for size in (3, 4, 5):
            e = FiberTransport(3, size, [(i, (i+1) % size) for i in range(size)])
            matrix = lifted_graph_laplacian(range(9), e.geometry.edges, [0]*len(e.geometry.edges), e.geometry.group, e.geometry.fiber_edges)
            base = [6-2*(np.cos(x)+np.cos(y)+np.cos(x+y)) for x in 2*np.pi*np.arange(3)/3 for y in 2*np.pi*np.arange(3)/3]
            fiber = [2-2*np.cos(2*np.pi*k/size) for k in range(size)]
            np.testing.assert_allclose(eigh(np.asarray(matrix, dtype=float), eigvals_only=True), sorted(a+b for a in base for b in fiber), atol=1e-12)
        e, record = self.experiments[0], self.data['fibers'][0]
        states = [{'links': r['checked']['final_links']} for r in record['orders']]
        exact = spectral_word_witness(e.geometry, states, full_bundle=True)
        self.assertEqual(exact['first_distinct_trace_power'], 12)
        self.assertEqual(exact['trace_at_first_distinct_power_by_lap'], ['547642051586352', '547642051586712'])
        operators = [np.asarray(lifted_graph_laplacian(range(e.geometry.size), e.geometry.edges, s['links'], e.geometry.group, e.geometry.fiber_edges), dtype=np.int64) for s in states]
        count = len(operators[0])
        self.assertLess(count*16**12, 2**63)
        powers = [np.eye(count, dtype=np.int64) for _ in states]
        for k in range(1, 13):
            powers = [p@a for p, a in zip(powers, operators)]
            traces = [int(np.trace(p)) for p in powers]
            self.assertEqual(traces[0] == traces[1], k < 12)
        self.assertEqual(traces[1]-traces[0], 360)


if __name__ == '__main__':
    unittest.main()
