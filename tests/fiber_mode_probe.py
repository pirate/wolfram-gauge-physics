import copy
import random
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'tools'))
from d4_loop_observer import CentralExtension, conversion_witness
from face_energy_obstruction import derive_fiber_group
from fiber_mode_probe import (analyze, connection_laplacian, identity, moments, multiply,
                              representation, transpose, verify_graph_restriction)
from run_shared_edge import mesh_geometry


class FiberModeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.g = derive_fiber_group()
        cls.modes = representation(cls.g)

    def test_fiber_adjacency_derives_the_real_two_component_action(self):
        modes, e = self.modes, CentralExtension(self.g)
        self.assertEqual(multiply(modes['adjacency'], modes['basis_columns']), [[0, 0]]*4)
        r, s, z = [modes['matrices'][x] for x in (e.r, e.s, e.z)]
        self.assertEqual(z, [[-1, 0], [0, -1]])
        self.assertEqual(multiply(r, r), identity(2))
        self.assertEqual(multiply(s, s), identity(2))
        self.assertEqual(multiply(r, s), [[-x for x in row] for row in multiply(s, r)])
        with self.assertRaisesRegex(ValueError, 'dimensions'):
            multiply([[1, 2]], [[1]])

    def test_covariance_and_exact_restriction_of_the_actual_lifted_graph(self):
        rng = random.Random(5804321)
        vertices = list(range(9))
        edges, _, _ = mesh_geometry(3)
        for _ in range(8):
            links = [rng.randrange(8) for _ in edges]
            frames = [rng.randrange(8) for _ in vertices]
            changed = [self.g.mul[frames[v]][self.g.mul[x][self.g.inv[frames[u]]]]
                       for (u, v), x in zip(edges, links)]
            operator, full = verify_graph_restriction(vertices, edges, links, self.g, self.modes)
            moved, _ = verify_graph_restriction(vertices, edges, changed, self.g, self.modes)
            gauge = [[0]*18 for _ in range(18)]
            for v in vertices:
                for a in range(2):
                    for b in range(2):
                        gauge[2*v+a][2*v+b] = self.modes['matrices'][frames[v]][a][b]
            self.assertEqual(moved, multiply(multiply(gauge, operator), transpose(gauge)))
            self.assertEqual(operator, transpose(operator))
            self.assertEqual(moments(operator, 4), moments(moved, 4))
            self.assertTrue(all(sum(row) == 0 for row in full))
        corrupt = copy.deepcopy(self.modes)
        corrupt['matrices'][0] = [[-1, 0], [0, -1]]
        with self.assertRaisesRegex(ValueError, 'restriction'):
            verify_graph_restriction(vertices, edges, [0]*len(edges), self.g, corrupt)

    def test_quadratic_form_is_graph_incidence_not_a_fitted_action(self):
        vertices = list(range(9))
        edges, _, _ = mesh_geometry(3)
        links = [(i*5+i//3) % 8 for i in range(len(edges))]
        operator = connection_laplacian(vertices, edges, links, self.modes['matrices'])
        rng = random.Random(591037)
        for _ in range(12):
            field = [rng.randrange(-3, 4) for _ in range(18)]
            quadratic = sum(field[i]*operator[i][j]*field[j] for i in range(18) for j in range(18))
            differences = 0
            for (u, v), x in zip(edges, links):
                for a in range(2):
                    transported = sum(self.modes['matrices'][x][a][b]*field[2*u+b] for b in range(2))
                    differences += (field[2*v+a]-transported)**2
            self.assertEqual(quadratic, differences)
            self.assertGreaterEqual(quadratic, 0)
        flat = connection_laplacian(vertices, edges, [0]*len(edges), self.modes['matrices'])
        for component in range(2):
            constant = [[int(i % 2 == component)] for i in range(18)]
            self.assertEqual(multiply(flat, constant), [[0]]*18)

    def test_same_face_classes_have_different_mode_and_bundle_graph_spectra(self):
        witness = conversion_witness()
        result = analyze(witness)
        self.assertEqual(result['first_distinct_mode_trace_power'], 4)
        a, b = result['branches']
        self.assertEqual(a['mode_trace_powers_1_through_8'][:3], b['mode_trace_powers_1_through_8'][:3])
        self.assertEqual((a['mode_trace_powers_1_through_8'][3], b['mode_trace_powers_1_through_8'][3]), (43444, 43284))
        self.assertEqual(a['lifted_graph_trace_powers_1_through_8'][3]-b['lifted_graph_trace_powers_1_through_8'][3], 160)


if __name__ == '__main__':
    unittest.main()
