import json
import random
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'tools'))
from fiber_transport import FiberTransport
from transport_braid import spectral_word_witness
from unrouted_transport import commuting_return_audit, detect_return, worldline_audit


class UnroutedTransportTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data = json.loads((ROOT/'data/unrouted-transport.json').read_text())

    def test_detector_reproduces_all_seeds_and_bounded_negatives(self):
        for r in self.data['runs']:
            e = FiberTransport(r['side'], 3, [(0, 1), (1, 2), (2, 0)])
            observed = detect_return(e, r['initial_links'], r['initial_positions'], r['seed'], r['attempt_budget'])
            for key, value in observed.items():
                if key == 'witness' and value is not None:
                    for k, v in value.items():
                        self.assertEqual(v, r[key][k])
                else:
                    self.assertEqual(value, r[key])
            self.assertEqual(e.compiled(r['initial_links'], r['schedule']), r['full_history_check'])
            if r['witness'] is None:
                self.assertEqual(r['attempts_used'], r['attempt_budget'])

    def test_passive_worldlines_match_the_changing_update_clock(self):
        for r in self.data['runs']:
            e = FiberTransport(r['side'], 3, [(0, 1), (1, 2), (2, 0)])
            whole = worldline_audit(e, r['initial_positions'], r['schedule'])
            self.assertEqual([[x['tick'], x['patch']] for x in whole['events']], [x[:2] for x in r['full_history_check']['events']])
            if r['witness'] is None:
                continue
            w = r['witness']; segment = r['schedule'][w['start_tick']:w['end_tick']]
            audit = worldline_audit(e, w['positions'], segment)
            self.assertEqual(audit, w['worldlines'])
            self.assertTrue(audit['all_tags_return'])
            self.assertEqual(audit['displacements'], [[x*3*r['side'] for x in d] for d in audit['torus_winding_by_tag']])
            for c in w['comparisons']:
                self.assertEqual([[x['tick'], x['patch']] for x in audit['events']], [x[:2] for x in c['segment_check']['events']])

    def test_controls_and_full_spectral_witnesses_are_independently_replayed(self):
        returns = []
        for r in self.data['runs']:
            if r['witness'] is None:
                continue
            e = FiberTransport(r['side'], 3, [(0, 1), (1, 2), (2, 0)])
            w = r['witness']; returns.append(w['comparisons'][1]['gauge_equivalent'])
            for c in w['comparisons']:
                before, after = c['prefix_check']['final_links'], c['segment_check']['final_links']
                self.assertEqual(e.geometry.sectors(before), e.geometry.sectors(after))
                self.assertEqual(e.factor.oracle.fan.gauge_equivalent(before, after), c['gauge_equivalent'])
                self.assertEqual(e.compiled(before, r['schedule'][w['start_tick']:w['end_tick']]), c['segment_check'])
                exact = spectral_word_witness(e.geometry, [{'links': before}, {'links': after}], full_bundle=True)
                self.assertEqual(exact, c['exact_full_bundle_spectrum'])
        self.assertIn(True, returns)
        self.assertIn(False, returns)

    def test_detection_is_covariant_under_local_frames(self):
        r = self.data['runs'][0]
        e = FiberTransport(r['side'], 3, [(0, 1), (1, 2), (2, 0)])
        g = e.geometry.group; rng = random.Random(5840)
        frames = [rng.randrange(g.n) for _ in range(e.geometry.size)]
        transform = lambda links: [g.mul[frames[v]][g.mul[x][g.inv[frames[u]]]] for (u, v), x in zip(e.geometry.edges, links)]
        changed = detect_return(e, transform(r['initial_links']), r['initial_positions'], r['seed'], r['attempt_budget'])
        self.assertEqual(changed['schedule'], r['schedule'])
        for key in ('start_tick', 'end_tick', 'positions'):
            self.assertEqual(changed['witness'][key], r['witness'][key])
        for key in ('before_links', 'after_links'):
            self.assertEqual(changed['witness'][key], transform(r['witness'][key]))

    def test_absorbing_vacuum_and_explicit_detector_bounds(self):
        e = FiberTransport(3, 3, [(0, 1), (1, 2), (2, 0)])
        flat = [0]*len(e.geometry.edges)
        result = detect_return(e, flat, [], 0, 1000)
        self.assertEqual(result['attempts_used'], 1000)
        self.assertEqual(result['changing_updates'], 0)
        self.assertIsNone(result['witness'])
        self.assertEqual(detect_return(e, flat, [], 0, 0)['schedule'], [])
        with self.assertRaisesRegex(ValueError, '100000'):
            detect_return(e, flat, [], 0, 100001)
        with self.assertRaisesRegex(ValueError, 'tags'):
            detect_return(e, flat, [0], 0, 10)

    def test_commuting_return_is_explained_by_torus_intersection_parity(self):
        for r in self.data['runs']:
            if r['witness'] is None:
                continue
            e = FiberTransport(r['side'], 3, [(0, 1), (1, 2), (2, 0)])
            w = r['witness']; control = w['comparisons'][1]
            audit = commuting_return_audit(e, control['prefix_check']['final_links'], control['segment_check']['final_links'], r['seed_generators'][0], w['worldlines'])
            self.assertEqual(audit, w['commuting_topology_audit'])
            self.assertEqual(audit['measured_period_bits'], audit['predicted_period_bits'])
            self.assertEqual(audit['measured_period_bits'] == [0, 0], control['gauge_equivalent'])


if __name__ == '__main__':
    unittest.main()
