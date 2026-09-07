import itertools
import json
import random
import sys
import unittest
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'tools'))
from d4_loop_observer import CentralExtension, LoopForest, compiled_observe, conversion_witness, local_census
from face_energy_obstruction import FiniteGroup, derive_fiber_group
from run_three_face import LinkOracle


class LoopObserverTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.g = derive_fiber_group()
        cls.extension = CentralExtension(cls.g)

    def test_extension_is_derived_nonsplit_and_label_independent(self):
        e = self.extension
        self.assertEqual((e.r, e.s, e.z), (1, 2, 5))
        quotient = list(itertools.product(range(2), repeat=2))
        add = lambda x, y: tuple(a ^ b for a, b in zip(x, y))
        cocycle = lambda x, y: x[1] & y[0]
        for x, y, z in itertools.product(quotient, repeat=3):
            self.assertEqual(cocycle(x, y) ^ cocycle(add(x, y), z), cocycle(y, z) ^ cocycle(x, add(y, z)))
        for values in itertools.product(range(2), repeat=3):
            section = dict(zip(quotient, (0, *values)))
            self.assertTrue(any(cocycle(x, y) != section[x] ^ section[y] ^ section[add(x, y)]
                                for x, y in itertools.product(quotient, repeat=2)))
        relabeled = CentralExtension(FiniteGroup(self.g.elements[::-1]))
        self.assertEqual(len(relabeled.bits), 8)
        cyclic = FiniteGroup([tuple((i+k) % 8 for i in range(8)) for k in range(8)])
        with self.assertRaisesRegex(ValueError, 'center'):
            CentralExtension(cyclic)

    def test_complete_observer_and_exact_hidden_fiber_cardinalities(self):
        rows = local_census(self.extension)
        self.assertEqual([r['gauge_orbits'] for r in rows], [1, 5, 28, 176, 1216])
        for r in rows:
            m = r['based_loops']
            self.assertEqual(r['gauge_orbits'], (8**m+3*4**m)//4)
        fibers = {}
        for values in itertools.product(range(8), repeat=4):
            p = self.extension.probe(values)
            fibers.setdefault(p['classes'], set()).add(p['relation_bits'])
        for classes, bits in fibers.items():
            self.assertEqual(len(bits), 2**self.extension.probe(classes)['hidden_bits_beyond_individual_loop_classes'])

    def test_three_loop_words_are_necessary_not_just_pairwise_classes(self):
        left, right = (1, 2, 3), (1, 2, 6)
        classes = lambda values: tuple(self.g.sectors[x] for x in values)
        self.assertEqual(classes(left), classes(right))
        pairs = lambda values: tuple(self.g.sectors[self.g.mul[values[i]][values[j]]]
                                     for i, j in itertools.combinations(range(3), 2))
        self.assertEqual(pairs(left), pairs(right))
        a, b = self.extension.probe(left), self.extension.probe(right)
        self.assertNotEqual(a['relation_bits'], b['relation_bits'])
        self.assertEqual(len(a['relation_words'][0]), 3)
        with self.assertRaisesRegex(ValueError, 'relation bits'):
            self.extension.reconstruct(a['classes'], [])

    def test_graph_orbits_match_cpp_and_legacy_under_local_frames(self):
        rng = random.Random(473029)
        graphs = [(range(5), [(0, 1), (0, 2), (0, 3), (0, 4), (1, 2), (2, 3), (3, 4)]),
                  ([2, 4, 5, 8, 9, 12, 99], [(2, 4), (2, 5), (4, 5), (8, 9), (8, 12), (9, 12)]),
                  ([], [])]
        for vertices, edges in graphs:
            forest = LoopForest(vertices, edges, self.g)
            states = []
            for _ in range(20):
                raw = [rng.randrange(8) for _ in edges]
                frames = {v: rng.randrange(8) for v in vertices}
                changed = [self.g.mul[frames[v]][self.g.mul[x][self.g.inv[frames[u]]]]
                           for x, (u, v) in zip(raw, edges)]
                self.assertEqual(forest.signature(raw), forest.signature(changed))
                self.assertEqual(forest.signature(raw), forest.signature(forest.representative(forest.signature(raw))))
                states.extend((raw, changed))
            compiled_observe(forest, states)
        with self.assertRaisesRegex(ValueError, 'duplicate'):
            LoopForest([0, 0], [], self.g)
        with self.assertRaisesRegex(ValueError, 'canonical'):
            LoopForest([0, 1], [(1, 0)], self.g)

    def test_cpp_observer_exhausts_four_loop_tuples_against_legacy(self):
        edges = [(0, 1), (0, 2), (0, 3), (0, 4), (1, 2), (1, 3), (2, 3), (3, 4)]
        forest = LoopForest(range(5), edges, self.g)
        batch, signatures = [], set()
        for values in itertools.product(range(8), repeat=4):
            batch.append(forest.representative([values]))
            if len(batch) == 128:
                for observed in compiled_observe(forest, batch):
                    signatures.add(tuple(observed['normalized'][0]))
                batch = []
        self.assertFalse(batch)
        self.assertEqual(len(signatures), 1216)

    def test_fundamental_and_relation_walks_are_actual_edge_transports(self):
        oracle = LinkOracle(4, self.g)
        forest = LoopForest(range(16), oracle.edges, self.g)
        links = [(i*3+i//5) % 8 for i in range(len(oracle.edges))]
        ids = {edge: i for i, edge in enumerate(oracle.edges)}
        transport = lambda walk: oracle.transport(links, [(ids[tuple(sorted((a, b)))], a > b)
                                                         for a, b in zip(walk, walk[1:])])
        based = forest.based_loops(links)[0]
        walks = [forest.fundamental_walk(0, i) for i in range(len(based))]
        self.assertEqual(tuple(transport(w) for w in walks), based)
        p = forest.extension.probe(based)
        for word, bit in zip(p['relation_words'], p['relation_bits']):
            # Path traversal composes on the left, opposite to algebraic word order.
            walk = [0]
            for i, direction in reversed(word):
                piece = walks[i] if direction == 1 else walks[i][::-1]
                walk.extend(piece[1:])
            self.assertEqual(transport(walk), self.extension.z if bit else self.g.identity)
        compiled_observe(forest, [links])

    def test_raw_conversion_hides_a_real_bit_and_uniform_averaging_erases_its_feedback(self):
        witness = conversion_witness()
        a, b = witness['branches']
        self.assertEqual(a['face_classes'], b['face_classes'])
        self.assertNotEqual(a['local_loop_probe']['relation_bits'], b['local_loop_probe']['relation_bits'])
        self.assertTrue(witness['not_gauge_equivalent'])
        self.assertTrue(witness['uniform_next_class_distribution_identical'])
        self.assertNotEqual(witness['specified_rule_responses'][0]['face_classes'], witness['specified_rule_responses'][1]['face_classes'])
        self.assertTrue(all(b['compiled_raw_inverse'] for b in witness['branches']))


if __name__ == '__main__':
    unittest.main()
