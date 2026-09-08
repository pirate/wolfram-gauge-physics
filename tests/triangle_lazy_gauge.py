import builtins
import itertools
import json
import random
import sys
import unittest
from pathlib import Path
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'tools'))
from d4_loop_observer import ConnectionForest
from triangle_charge_current import CurrentProbe
from triangle_lazy_gauge import GaugeConnection, LazyFrames, audit
from triangle_reference import AxialReference


def nested_write_forest(experiment):
    """A second valid coordinate tree, including nested edges inside triangles."""
    edges, g = experiment.geometry.edges, experiment.geometry.group
    forest = ConnectionForest(range(experiment.geometry.size), edges, g)
    forced = sorted(experiment.factor.oracle.fan.writes[0])
    center = next(iter(set(edges[forced[0]]) & set(edges[forced[1]])))
    root = next(v for v in edges[forced[0]] if v != center)
    representatives = list(range(experiment.geometry.size))
    def find(v):
        while representatives[v] != v:
            v = representatives[v]
        return v
    tree = set()
    for edge in forced+sorted(set(range(len(edges)))-set(forced)):
        a, b = [find(v) for v in edges[edge]]
        if a != b:
            representatives[b] = a
            tree.add(edge)
    adjacency = {v: [] for v in range(experiment.geometry.size)}
    for edge, (u, v) in enumerate(edges):
        if edge in tree:
            adjacency[u].append((v, edge, False))
            adjacency[v].append((u, edge, True))
    queue, seen, parents = [root], {root}, {}
    for u in queue:
        for v, edge, reverse in adjacency[u]:
            if v not in seen:
                seen.add(v); queue.append(v)
                parents[v] = (u, edge, reverse)
    forest.components = [{'root': root, 'queue': queue, 'parents': parents,
                          'chords': sorted(set(range(len(edges)))-tree)}]
    return forest


class TriangleLazyGaugeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.bank = json.loads((ROOT/'data/triangle-feedback.json').read_text())['bank']
        cls.p = CurrentProbe(3, cls.bank)
        cls.g = cls.p.e.geometry.group

    def test_all_small_noncommuting_interval_pairs_match_explicit_vertex_frames(self):
        intervals = [(i, j) for i in range(4) for j in range(i+1, 5)]
        for first, second in itertools.product(intervals, repeat=2):
            for a, b in itertools.product(range(6), repeat=2):
                tree, expected = LazyFrames(4, self.g), [0]*4
                for (lo, hi), h in ((first, a), (second, b)):
                    tree.left_multiply(lo, hi, h)
                    expected[lo:hi] = [self.g.mul[h][v] for v in expected[lo:hi]]
                    self.assertEqual([tree.point(i) for i in range(4)], expected)
        self.assertNotEqual(self.g.mul[1][2], self.g.mul[2][1])

    def test_interleaved_lazy_pushes_and_queries_preserve_chronological_order(self):
        rng = random.Random(994702)
        for size in (1, 7, 32, 257):
            frames, expected = LazyFrames(size, self.g), [0]*size
            for _ in range(400):
                lo = rng.randrange(size); hi = rng.randrange(lo+1, size+1); value = rng.randrange(6)
                frames.left_multiply(lo, hi, value)
                expected[lo:hi] = [self.g.mul[value][v] for v in expected[lo:hi]]
                i = rng.randrange(size)
                self.assertEqual(frames.point(i), expected[i])
            self.assertEqual([frames.point(i) for i in range(size)], expected)
            self.assertEqual(frames.all_points(), expected)
        for bad in ((-1, 1, 0), (0, 0, 0), (0, 258, 0), (0, 1, 6), (False, 1, 0)):
            with self.assertRaises(ValueError):
                frames.left_multiply(*bad)
        for index in (-1, 257, 0.5, True):
            with self.assertRaises(ValueError):
                frames.point(index)

    def test_all_cpp_whole_mesh_loop_histories_and_reverse_schedules_reproduce(self):
        data = json.loads((ROOT/'data/triangle-lazy-gauge.json').read_text())
        self.assertEqual(audit(), data)
        self.assertEqual([r['independent_loops'] for r in data['records']], [19, 73, 289, 1153])
        for row in data['records']:
            work, height = row['work_maxima'], (row['vertices']-1).bit_length()
            self.assertLessEqual(work['point_queries'], 5)
            self.assertLessEqual(work['range_updates'], 2)
            self.assertLessEqual(work['point_nodes'], 5*(height+1))
            self.assertLessEqual(work['range_nodes'], 2*(4*height+1))

    def test_subtree_intervals_and_nested_tree_writes_match_an_alternative_forest(self):
        e, rng = self.p.e, random.Random(41802)
        forest = nested_write_forest(e)
        raw = [rng.randrange(6) for _ in e.geometry.edges]
        state = GaugeConnection(e, raw, forest)
        for vertex in forest.components[0]['queue']:
            descendants = {vertex}
            for child in forest.components[0]['queue']:
                if child in forest.components[0]['parents'] and forest.components[0]['parents'][child][0] in descendants:
                    descendants.add(child)
            self.assertEqual({v for v in state.entry if state.entry[vertex] <= state.entry[v] < state.exit[vertex]}, descendants)
        schedule = [rng.randrange(13*len(e.factor.pairs)) for _ in range(2000)]
        tree_pairs = 0
        initial = raw[:]
        for op in schedule:
            rule, support = divmod(op, len(e.factor.pairs))
            e.factor.oracle.update(raw, rule, support, e.tables[rule])
            step = state.step(op)
            tree_pairs += step['tree_writes'] == 2
            self.assertEqual(state.loops(), forest.based_loops(raw)[0])
        self.assertGreater(tree_pairs, 0)
        for op in reversed(schedule):
            state.step(op)
        self.assertEqual(state.loops(), forest.based_loops(initial)[0])

    def test_ancestor_first_restoration_has_a_detectable_nonabelian_failure(self):
        e, rng, found = self.p.e, random.Random(190480), False
        forest = nested_write_forest(e)
        for _ in range(32):
            initial = [rng.randrange(6) for _ in e.geometry.edges]
            for operator in range(13*len(e.factor.pairs)):
                rule, support = divmod(operator, len(e.factor.pairs))
                raw = initial[:]
                e.factor.oracle.update(raw, rule, support, e.tables[rule])
                state = GaugeConnection(e, initial, forest)
                with patch('triangle_lazy_gauge.sorted', side_effect=lambda items, reverse=False: builtins.sorted(items, reverse=False)):
                    step = state.step(operator)
                if step['tree_writes'] == 2 and state.loops() != forest.based_loops(raw)[0]:
                    found = True
                    break
            if found:
                break
        self.assertTrue(found, 'negative control must expose a noncommuting nested-tree error')

    def test_local_frame_transformations_leave_only_the_expected_root_conjugation(self):
        p, rng = CurrentProbe(6, self.bank), random.Random(481170)
        e, g = p.e, p.e.geometry.group
        initial = [rng.randrange(6) for _ in e.geometry.edges]
        frames = [rng.randrange(6) for _ in range(e.geometry.size)]
        frames[e.forest.components[0]['root']] = 3
        transformed = [g.mul[frames[v]][g.mul[x][g.inv[frames[u]]]] for (u, v), x in zip(e.geometry.edges, initial)]
        a, b = GaugeConnection(e, initial), GaugeConnection(e, transformed)
        for _ in range(500):
            op = rng.randrange(13*len(e.factor.pairs))
            self.assertEqual(a.step(op)['changed'], b.step(op)['changed'])
            self.assertEqual(b.loops(), tuple(g.conj[3][v] for v in a.loops()))

    def test_flat_connections_with_different_global_handles_remain_distinguishable(self):
        reference = AxialReference(3, self.bank, 0)
        trivial = reference.reconstruct(0, 0, [0]*18)
        twisted = reference.reconstruct(1, 0, [0]*18)
        self.assertEqual(self.p.charge(trivial), self.p.charge(twisted))
        self.assertEqual(self.p.charge(trivial), [0]*18)
        a, b = GaugeConnection(self.p.e, trivial), GaugeConnection(self.p.e, twisted)
        first, second = a.loops(), b.loops()
        self.assertNotEqual(first, second)
        for op in range(13*len(self.p.e.factor.pairs)):
            self.assertFalse(a.step(op)['changed'])
            self.assertFalse(b.step(op)['changed'])
        self.assertEqual((a.loops(), b.loops()), (first, second))

    def test_update_path_does_not_hide_a_whole_mesh_observation(self):
        e, rng = self.p.e, random.Random(785790)
        initial = [rng.randrange(6) for _ in e.geometry.edges]
        state = GaugeConnection(e, initial)
        fail = AssertionError('whole mesh observation called from primitive update')
        with patch.object(e.forest, 'based_loops', side_effect=fail), \
             patch.object(state, 'loops', side_effect=fail), \
             patch.object(state, 'materialize', side_effect=fail), \
             patch.object(state.frames, 'all_points', side_effect=fail):
            for _ in range(500):
                state.step(rng.randrange(13*len(e.factor.pairs)))
        for bad in (-1, 13*len(e.factor.pairs), 1.0, True):
            with self.assertRaises(ValueError):
                state.step(bad)


if __name__ == '__main__':
    unittest.main()
