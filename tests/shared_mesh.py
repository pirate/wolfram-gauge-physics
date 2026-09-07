import copy
import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'tools'))
from run_mesh_experiments import analyze_result


class SharedMeshTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data = json.loads((ROOT/'data/d4-shared-mesh.json').read_text())
        census = json.loads((ROOT/'data/d4-involution-braid-search.json').read_text())
        cls.table = census['solutions'][cls.data['solution_id']]['table']

    def test_exhaustive_shared_face_algebra_and_charge_loss(self):
        result = analyze_result(copy.deepcopy(self.data), self.table)
        self.assertEqual(result['local_census']['independent_algebra_checks'], 4096)
        self.assertEqual(result['local_census']['nonclosed_incoming_sector_tuples'], 108)
        self.assertEqual(result['local_census']['additive_charge_basis'], [])
        self.assertEqual(result['stationary_reference']['noncommuting_pair']['stationary_active_fraction'], .875)
        self.assertEqual(result['stationary_reference']['reflection']['stationary_active_fraction'], .75)

    def test_missing_state_or_wrong_local_output_is_rejected(self):
        result = copy.deepcopy(self.data)
        result['local_census']['transitions'].pop()
        with self.assertRaisesRegex(ValueError, 'not exhaustive'):
            analyze_result(result, self.table)
        result = copy.deepcopy(self.data)
        result['local_census']['transitions'][0][4] = 1
        with self.assertRaisesRegex(ValueError, 'independent based-loop algebra'):
            analyze_result(result, self.table)

    def test_fixed_boundary_witness_with_independent_link_algebra(self):
        data = self.data
        group = [tuple(g) for g in data['group_permutations']]
        n = len(group)
        mul = [[group.index(tuple(a[b[i]] for i in range(4))) for b in group] for a in group]
        inv = [next(b for b in range(n) if mul[a][b] == 0) for a in range(n)]
        def conj(g, a): return mul[g][mul[a][inv[g]]]
        sectors = [min(conj(g, a) for g in range(n)) for a in range(n)]
        local = data['local_census']
        witness = local['fixed_boundary_witness']
        edges = [tuple(edge) for edge in witness['edges']]
        def edge(values, u, v):
            value = values[edges.index(tuple(sorted((u, v))))]
            return value if u < v else inv[value]
        def transport(values, path):
            result = 0
            for u, v in zip(path, path[1:]): result = mul[edge(values, u, v)][result]
            return result
        boundary = local['outer_boundary']
        first, second = witness['states']
        for u, v in zip(boundary, boundary[1:]):
            expected = edge(first['initial_links'], u, v)
            for state in witness['states']:
                self.assertEqual(edge(state['initial_links'], u, v), expected)
                self.assertEqual(edge(state['final_links'], u, v), expected)
        before, after = [], []
        for state in witness['states']:
            values = state['initial_links'][:]
            before.append([sectors[transport(values, face)] for face in local['faces']])
            a = transport(values, local['patch']['first'])
            b = transport(values, local['patch']['second'])
            t = transport(values, local['patch']['connector'])
            c, d = divmod(self.table[a*n+conj(t, b)], n)
            d = conj(inv[t], d)
            for path, old, target in ((local['patch']['first'], a, c), (local['patch']['second'], b, d)):
                u, v = path[-2:]
                value = mul[mul[target][inv[old]]][edge(values, u, v)]
                values[edges.index(tuple(sorted((u, v))))] = value if u < v else inv[value]
            self.assertEqual(values, state['final_links'])
            after.append([sectors[transport(values, face)] for face in local['faces']])
        self.assertEqual(before[0], before[1])
        self.assertNotEqual(after[0], after[1])
        # Exhaust all root frames; the boundary spanning tree determines the rest.
        equivalent = False
        for root in range(n):
            frames = {boundary[0]: root}
            for u, v in zip(boundary[:-2], boundary[1:-1]):
                frames[v] = mul[edge(second['initial_links'], u, v)][mul[frames[u]][inv[edge(first['initial_links'], u, v)]]]
            if all(mul[frames[v]][mul[edge(first['initial_links'], u, v)][inv[frames[u]]]]
                   == edge(second['initial_links'], u, v) for u, v in edges):
                equivalent = True
        self.assertFalse(equivalent)


if __name__ == '__main__':
    unittest.main()
