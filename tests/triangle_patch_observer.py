import itertools
import json
import random
import sys
import unittest
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'tools'))
from triangle_charge_current import CurrentProbe
from triangle_patch_observer import PatchUnion, audit


class TrianglePatchObserverTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.bank = json.loads((ROOT/'data/triangle-feedback.json').read_text())['bank']
        cls.data = json.loads((ROOT/'data/triangle-patch-observer.json').read_text())
        cls.p = CurrentProbe(6, cls.bank)
        cls.u = PatchUnion(cls.p.e, (0, 2))

    def test_complete_census_and_transition_evidence_reproduce(self):
        self.assertEqual(audit(), self.data)
        self.assertEqual([r['gauge_orbits'] for r in self.data['word_censuses']], [49, 251, 1393])
        self.assertEqual(self.data['raw_union_connections'], 6**11)

    def test_every_raw_tuple_has_the_correct_short_word_and_anchored_orbit(self):
        words = self.u.words
        for rank in range(1, 6):
            signatures, anchored = {}, {}
            for raw in itertools.product(range(6), repeat=rank):
                canonical = words.canonical(raw)
                signature, anchor = words.signature(raw), words.anchored(raw)
                self.assertEqual(signatures.setdefault(signature, canonical), canonical)
                self.assertEqual(anchored.setdefault(anchor, canonical), canonical)
                self.assertEqual(words.canonical(words.reconstruct(anchor)), canonical)
            self.assertEqual(len(signatures), len(words.representatives(rank)))
            self.assertEqual(len(anchored), len(signatures))

    def test_each_extra_triple_word_is_necessary_with_individual_classes_fixed(self):
        w = self.u.words
        for row in self.data['triple_product_necessity']:
            left, right = row['tuples']
            a, b = w.signature(left), w.signature(right)
            omitted = row['omitted_word_index']
            self.assertNotEqual(w.canonical(left), w.canonical(right))
            self.assertEqual(a[:omitted]+a[omitted+1:], b[:omitted]+b[omitted+1:])
            self.assertNotEqual(a[omitted], b[omitted])

    def test_anchored_coordinates_are_linear_size_and_reject_invalid_labels(self):
        w, rng = self.u.words, random.Random(759120)
        raw = tuple(rng.randrange(6) for _ in range(100000))
        signature = w.anchored(raw)
        self.assertEqual(sum(map(len, signature)), 2*len(raw))
        self.assertEqual(w.canonical(w.reconstruct(signature)), w.canonical(raw))
        for row in w.g.conj:
            self.assertEqual(w.anchored(tuple(row[v] for v in raw)), signature)
        for bad in (((0,), (1,)), ((2,), (1,)), ((1,), (2,)), ((1,), ()), ((4,), (0,))):
            with self.assertRaises(ValueError):
                w.reconstruct(bad)

    def test_double_cosets_reconstruct_every_union_orbit_without_extra_frame_choices(self):
        data, u = self.data, self.u
        self.assertEqual(len(data['gluing_specs']), 845)
        self.assertEqual(Counter(len(states) for _, states in data['separate_observer_fibers']),
                         {1: 432, 2: 305, 3: 99, 6: 9})
        self.assertEqual(len({tuple(k) for k in data['sewn_coordinates']}), 1393)
        for spec in data['gluing_specs']:
            self.assertEqual(sorted(x for cell in spec['double_cosets'] for x in cell), spec['common_centralizer'])
        for i, values in enumerate(u.reps):
            self.assertEqual(u.observe(u.lift(values)), i)
            self.assertEqual(list(u.sewn(u.lift(values))), data['sewn_coordinates'][i])
        # Orbit-stabilizer weights reconstruct all raw assignments to 11 links,
        # rather than silently treating gauge orbits as equally sized.
        stabilizers = [sum(tuple(row[v] for v in state) == state for row in u.g.conj) for state in u.reps]
        self.assertEqual(sum(6**len(u.vertices)//s for s in stabilizers), 6**len(u.edges))
        for table in data['union_transition_tables']:
            self.assertEqual([stabilizers[j] for j in table], stabilizers)

    def test_every_one_face_overlap_uses_actual_connector_transport(self):
        e, rng, checked = self.p.e, random.Random(891330), 0
        candidates = [(0, b) for b, reads in enumerate(e.factor.oracle.fan.reads)
                      if e.factor.oracle.fan.writes[0] & reads]
        candidates.append((7, 0))  # Check the reverse direction of a nontrivial connector.
        for pair in candidates:
            if len(set(e.fan_faces[pair[0]]) & set(e.fan_faces[pair[1]])) != 1:
                continue
            u = PatchUnion(e, pair)
            self.assertEqual(len({u.sewn(u.lift(v)) for v in u.reps}), 1393)
            checked += 1
            for _ in range(3):
                raw = [rng.randrange(6) for _ in e.geometry.edges]
                frames = [rng.randrange(6) for _ in range(e.geometry.size)]
                moved = [u.g.mul[frames[v]][u.g.mul[x][u.g.inv[frames[a]]]]
                         for (a, v), x in zip(e.geometry.edges, raw)]
                self.assertEqual(u.observe(raw), u.observe(moved))
                self.assertEqual(u.sewn(raw), u.sewn(moved))
        self.assertEqual(checked, 13)

    def test_omitting_the_connector_is_detectably_not_gauge_invariant(self):
        u = PatchUnion(self.p.e, (0, 7))
        self.assertTrue(u.connector)
        found = False
        for values in u.reps:
            raw = u.lift(values)
            if u.triples(raw)[0][u.shared_positions[0]] != 0:
                continue
            correct = u.sewn(raw)
            for h in range(1, 6):
                frames = [0]*self.p.e.geometry.size
                frames[u.patch_roots[1]] = h
                moved = [u.g.mul[frames[b]][u.g.mul[x][u.g.inv[frames[a]]]]
                         for (a, b), x in zip(self.p.e.geometry.edges, raw)]
                self.assertEqual(u.observe(raw), u.observe(moved))
                self.assertEqual(correct, u.sewn(moved))
                connector, u.connector = u.connector, []
                try:
                    found = u.sewn(raw) != u.sewn(moved)
                finally:
                    u.connector = connector
                if found:
                    break
            if found:
                break
        self.assertTrue(found, 'negative control must change under a pure local frame transformation')

    def test_all_61292_quotient_transitions_match_cpp_raw_operator_inverse_scans(self):
        u, e, data = self.u, self.p.e, self.data
        self.assertEqual(len(data['channels']), 44)
        self.assertEqual(sum(rule == 0 for _, rule in data['channels']), 8)
        encoded = [rule*len(e.factor.pairs)+patch for patch, rule in data['channels']]
        scan = [op for op in encoded for _ in range(2)]
        checked = 0
        for start in range(0, len(u.reps), 16):
            values = u.reps[start:start+16]
            inputs = [u.lift(v) for v in values]
            runs = e.compiled_bank(inputs, scan, len(scan), capture_links=True)['runs']
            for local, initial in enumerate(inputs):
                run = runs[2*local+1]
                events = {event[0]: links for event, links in zip(run['events'], run['raw_event_links'])}
                for channel in range(len(encoded)):
                    target = events.get(2*channel+1, initial)
                    observed = u.observe(target)
                    self.assertEqual(observed, data['union_transition_tables'][channel][start+local])
                    self.assertEqual(list(u.sewn(target)), data['sewn_coordinates'][observed])
                    self.assertEqual(events.get(2*channel+2, initial), initial)
                    checked += 1
        self.assertEqual(checked, 1393*44)

    def test_gauge_transformed_raw_inputs_follow_the_same_kernel(self):
        u, e, rng = self.u, self.p.e, random.Random(24617)
        inputs, ids = [], [0, 7, 500, 931, 937, 1100, 1392]
        for state in ids:
            raw = u.lift(u.reps[state])
            frames = [rng.randrange(6) for _ in range(e.geometry.size)]
            inputs.append([u.g.mul[frames[v]][u.g.mul[x][u.g.inv[frames[a]]]]
                           for (a, v), x in zip(e.geometry.edges, raw)])
        encoded = [rule*len(e.factor.pairs)+patch for patch, rule in self.data['channels']]
        scan = [op for op in encoded for _ in range(2)]
        runs = e.compiled_bank(inputs, scan, len(scan), capture_links=True)['runs']
        for i, state in enumerate(ids):
            run = runs[2*i+1]
            events = {event[0]: links for event, links in zip(run['events'], run['raw_event_links'])}
            for channel in range(len(encoded)):
                self.assertEqual(u.observe(events.get(2*channel+1, inputs[i])),
                                 self.data['union_transition_tables'][channel][state])

    def test_three_step_mean_difference_uses_unquotiented_raw_path_counts(self):
        data, e = self.data, self.p.e
        expected, witness = data['mean_response_witness'], data['gluing_witness']
        faces, means = expected['face_ids'], []
        for initial in witness['raw_links']:
            paths, rows = Counter({tuple(initial): 1}), []
            for tick in range(4):
                self.assertEqual(sum(paths.values()), 44**tick)
                rows.append([sum(weight*self.p.charge(raw)[f] for raw, weight in paths.items()) for f in faces])
                if tick == 3:
                    break
                after = Counter()
                for raw, count in paths.items():
                    for patch, rule in data['channels']:
                        table = e.tables[rule]
                        code = e.factor.oracle.code(raw, rule, patch)
                        moved = list(raw)
                        if code != table[code]:
                            e.factor.oracle.update(moved, rule, patch, table)
                        after[tuple(moved)] += count
                paths = after
            means.append(rows)
        for row, a, b in zip(expected['differences'], *means):
            self.assertEqual(row['difference_numerator'], [x-y for x, y in zip(a, b)])
        self.assertEqual(expected['first_mean_difference_attempt'], 3)
        self.assertEqual(expected['differences'][-1]['difference_numerator'], [24, -24, 24, -48, 24])


if __name__ == '__main__':
    unittest.main()
