import itertools
import json
import random
import sys
import unittest
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'tools'))
from bank_equilibrium import compiled_run
from transport_braid import TransportExperiment, algebra_census, hurwitz, pure_word, spectral_word_witness, triangle_candidate


class TransportBraidTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data = json.loads((ROOT/'data/d4-transport-braid.json').read_text())
        cls.experiment = TransportExperiment(cls.data['side'])

    def transform(self, links, frames):
        e = self.experiment
        g = e.geometry.group
        return [g.mul[frames[v]][g.mul[x][g.inv[frames[u]]]] for (u, v), x in zip(e.geometry.edges, links)]

    def test_complete_pure_braid_orbit_census(self):
        self.assertEqual(json.loads(json.dumps(algebra_census())), self.data['algebra'])
        self.assertEqual(self.data['algebra']['census'][-1]['orbit_size_counts_over_raw_input_tuples'],
                         {'1': 1696, '2': 1440, '4': 960})

    def test_next_fiber_probe_derives_noncommuting_gauge_orbit_actions(self):
        probe = triangle_candidate()
        self.assertEqual(json.loads(json.dumps(probe)), self.data['next_fiber_algebra_probe'])
        self.assertEqual(len(probe['automorphisms_derived_from_adjacency']), 6)
        self.assertEqual(probe['paired_neutral_seed_and_generator_pair_cases'], 540)
        self.assertGreater(probe['noncommuting_on_gauge_orbits'], 0)
        w = probe['witness']
        self.assertNotEqual(w['canonical_first_then_second'], w['canonical_second_then_first'])

    def test_actual_hurwitz_words_preserve_products_and_pure_actions_commute(self):
        g = self.experiment.geometry.group
        for values in itertools.product(range(8), repeat=3):
            for i in (0, 1):
                changed = hurwitz(g, values, i)
                self.assertEqual(hurwitz(g, changed, i, inverse=True), values)
                product = lambda t: g.mul[t[0]][g.mul[t[1]][t[2]]]
                self.assertEqual(product(changed), product(values))
            pairs = list(itertools.combinations(range(3), 2))
            for a in pairs:
                self.assertEqual(pure_word(g, pure_word(g, values, *a), *a), values)
                for b in pairs:
                    self.assertEqual(pure_word(g, pure_word(g, values, *a), *b),
                                     pure_word(g, pure_word(g, values, *b), *a))

    def test_saved_histories_replay_on_raw_links_and_reverse_exactly(self):
        e = self.experiment
        for record in self.data['experiments']:
            links = record['initial_links'][:]
            checkpoints = {record['preparation_events']+12*i: s['links'] for i, s in enumerate(record['states'])}
            for number, (tick, rule, patch, code, target) in enumerate(record['events'], 1):
                self.assertEqual(tick, number)
                self.assertEqual(rule, 0)
                self.assertEqual(e.factor.oracle.update(links, rule, patch, e.factor.tables[rule]), (code, target))
                if number in checkpoints:
                    self.assertEqual(links, checkpoints[number])
            for _, rule, patch, code, target in reversed(record['events']):
                self.assertEqual(e.factor.oracle.update(links, rule, patch, e.factor.tables[rule]), (target, code))
            self.assertEqual(links, record['initial_links'])
            self.assertTrue(record['same_face_classes'])
            self.assertTrue(all(s['face_classes'] == record['states'][0]['face_classes'] for s in record['states']))

    def test_local_frame_covariance_of_the_entire_compiled_transport_history(self):
        e = self.experiment
        rng = random.Random(611704)
        for record in self.data['experiments']:
            frames = [rng.randrange(8) for _ in range(e.geometry.size)]
            initial = self.transform(record['initial_links'], frames)
            schedule = [event[2] for event in record['events']]
            checked = compiled_run(e.factor, initial, schedule, len(schedule))[0]
            self.assertEqual(checked['final_links'], self.transform(record['states'][-1]['links'], frames))
            self.assertTrue(checked['exact_link_inverse'])
            self.assertTrue(checked['independent_link_replay'])
            for s in record['states']:
                self.assertEqual(e.forest.signature(self.transform(s['links'], frames)), e.forest.signature(s['links']))

    def test_closed_windings_controls_and_explicit_frame_return_witnesses(self):
        e = self.experiment
        for number, record in enumerate(self.data['experiments']):
            initial = record['states'][0]['links']
            expected = [True, False, True] if number == 0 else [True]*3
            self.assertEqual(record['gauge_return_by_lap'], expected)
            for s, frames, equal in zip(record['states'], record['frames_from_lap_zero'], expected):
                self.assertEqual(e.frame_witness(initial, s['links']) is not None, equal)
                if equal:
                    self.assertEqual(self.transform(initial, frames), s['links'])
                else:
                    self.assertIsNone(frames)
            for face, winding in record['interior_windings'].items():
                self.assertEqual(e.winding(record['cycle'], int(face)), winding)
            if number == 2:
                self.assertTrue(all(e.winding(record['cycle'], int(face)) == 0 for face in record['interior_windings']))
        main = self.data['experiments'][0]
        self.assertFalse(main['raw_return_by_lap'][2])  # frame equivalence, not a raw echo
        self.assertEqual([r['central_holonomy_by_lap'] for r in main['same_class_loop_relations']], [[0, 5, 0], [5, 0, 5]])
        self.assertEqual(json.loads(json.dumps(e.loop_relations(main['states']))), main['same_class_loop_relations'])

    def test_all_saved_empty_star_deformations_are_single_vertex_frame_changes(self):
        e = self.experiment
        self.assertEqual([len(r['empty_star_deformations']) for r in self.data['experiments']], [12, 12, 24])
        for record in self.data['experiments']:
            reference = record['states'][1]['links']
            for detour in record['empty_star_deformations']:
                frames = detour['frames_from_reference_lap']
                self.assertEqual({v for v, x in enumerate(frames) if x}, {detour['primal_vertex']})
                self.assertEqual(self.transform(reference, frames), detour['final_links'])
                links, events = record['states'][0]['links'][:], []
                for u, v in zip(detour['walk'], detour['walk'][1:]):
                    e.step(links, u, v, events)
                self.assertEqual(events, detour['events'])
                self.assertEqual(links, detour['final_links'])
                for face, winding in record['interior_windings'].items():
                    self.assertEqual(e.winding(detour['walk'], int(face)), winding)

    def test_exact_spectral_memory_has_an_independent_integer_walk_moment(self):
        e = self.experiment
        for record in self.data['experiments']:
            self.assertEqual(spectral_word_witness(e.geometry, record['states']), record['exact_spectral_witness'])
            self.assertEqual(spectral_word_witness(e.geometry, record['states'], full_bundle=True),
                             record['exact_full_bundle_spectral_witness'])
        main = self.data['experiments'][0]
        self.assertEqual(main['exact_spectral_witness']['first_distinct_trace_power'], 14)
        # Signed row absolute sum is 12; all intermediate matrix powers and
        # traces through power 14 fit int64. This is independent of FLINT/Newton.
        self.assertLess(2*e.geometry.size*12**14, 2**63)
        powers = [np.eye(2*e.geometry.size, dtype=np.int64) for _ in range(2)]
        operators = [e.geometry.operator(s['links']).toarray().astype(np.int64) for s in main['states'][:2]]
        for k in range(1, 15):
            powers = [p@a for p, a in zip(powers, operators)]
            traces = [int(np.trace(p)) for p in powers]
            self.assertEqual(traces[0] == traces[1], k < 14)
        self.assertEqual(traces[0]-traces[1], 25200)
        full = main['exact_full_bundle_spectral_witness']
        self.assertEqual(full['first_distinct_trace_power'], 14)
        self.assertEqual(int(full['trace_at_first_distinct_power_by_lap'][0])-
                         int(full['trace_at_first_distinct_power_by_lap'][1]), 25200)


if __name__ == '__main__':
    unittest.main()
