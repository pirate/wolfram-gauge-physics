import copy
import json
import random
import sys
import unittest
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'tools'))
from triangle_spreading import SpreadingProbe, audit, digest, distances, run_case


class TriangleSpreadingTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.bank = json.loads((ROOT/'data/triangle-feedback.json').read_text())['bank']
        cls.data = json.loads((ROOT/'data/triangle-spreading.json').read_text())
        cls.p = SpreadingProbe(6, cls.bank)

    def test_published_small_mesh_case_reproduces_full_compiled_and_gauge_controls(self):
        self.assertEqual(run_case(6, 0, self.bank), self.data['cases'][0])

    def test_three_size_protocol_and_exact_waiting_law(self):
        cases = self.data['cases']
        self.assertEqual([(c['side'], c['trial']) for c in cases], [(s, t) for s in (6, 12, 24) for t in (0, 1)])
        self.assertEqual(sum(r['evolution']['event_count'] for c in cases for r in c['runs']), 319207)
        for case in cases:
            side, m = case['side'], case['attempted_operators']
            self.assertEqual(m, 78*side*side)
            self.assertEqual(case['attempts'], 12000*(side//6)**2)
            self.assertEqual(Fraction(*case['first_change_probability_per_attempt']), Fraction(96, m))
            self.assertEqual(Fraction(*case['geometric_mean_first_change_attempt']), Fraction(m, 96))
            self.assertTrue(case['frozen_and_reaction_disabled_controls_fixed'])
            self.assertEqual(len(case['runs']), 6)
            rng = random.Random(case['schedule_seed'])
            schedule = [rng.randrange(m) for _ in range(case['attempts'])]
            self.assertEqual(digest(schedule), case['schedule_sha256'])

    def test_saved_spatial_observations_reconstruct_from_faces_and_first_passage_times(self):
        for case in self.data['cases']:
            p = SpreadingProbe(case['side'], self.bank)
            for row in case['runs']:
                seed, run = row['seed'], row['evolution']
                d = distances(p.dual_neighbors, p.incident[seed['edge']])
                self.assertEqual(d, run['face_distances_from_seed'])
                self.assertEqual(distances(p.dependency_neighbors, [seed['edge']]), run['link_dependency_distances_from_seed'])
                self.assertEqual(p.measure(seed['links'])['changing_operators'], 96)
                first = run['first_charge_departure_tick']
                self.assertEqual(len(first), p.faces)
                previous_visited = 0
                for snapshot in run['snapshots']:
                    tick, q = snapshot['tick'], snapshot['charge_field']
                    visited = [f for f, t in enumerate(first) if t is not None and t <= tick]
                    self.assertGreaterEqual(len(visited), previous_visited)
                    previous_visited = len(visited)
                    self.assertEqual(sum(q), p.faces)
                    self.assertEqual(snapshot['activity']['populations'], [q.count(i) for i in range(3)])
                    self.assertFalse(snapshot['activity']['frozen'])
                    self.assertEqual(snapshot['ever_contrast_radius'], max((d[f] for f in visited), default=None))
                    self.assertEqual(snapshot['contrast_radius'], max((d[f] for f, h in enumerate(q) if h != 1), default=None))
                    for radius, size, n0, n1, n2, ever in snapshot['shells_radius_size_n0_n1_n2_ever']:
                        fs = [f for f, r in enumerate(d) if r == radius]
                        self.assertEqual(size, len(fs))
                        self.assertEqual([n0, n1, n2], [sum(q[f] == h for f in fs) for h in range(3)])
                        self.assertEqual(ever, sum(first[f] is not None and first[f] <= tick for f in fs))
                    seam = run['first_changing_support_touching_periodic_seam']
                    self.assertEqual(snapshot['before_any_changing_support_touches_seam'], seam is None or tick < seam)
                self.assertEqual(p.charge(run['final_links']), run['snapshots'][-1]['charge_field'])
                self.assertTrue(all(t is not None for t in first))
                self.assertEqual(p.measure(run['final_links']), run['snapshots'][-1]['activity'])
                self.assertTrue(run['each_changing_event_has_seed_ancestry'])
                self.assertTrue(run['exact_link_inverse'])

    def test_short_regression_accepts_a_window_with_no_changing_events(self):
        data = audit((6,), 1, 20)
        self.assertEqual(len(data['cases'][0]['runs']), 6)
        for row in data['cases'][0]['runs']:
            if row['evolution']['event_count'] == 0:
                self.assertIsNone(row['evolution']['first_event'])
                self.assertTrue(all(t is None for t in row['evolution']['first_charge_departure_tick']))
                self.assertFalse(row['evolution']['snapshots'][-1]['activity']['frozen'])

    def test_analyzer_rejects_unseeded_events_bad_ticks_and_bad_targets(self):
        p = self.p
        e = p.e
        background, seeds = p.initial_conditions()
        initial, edge = seeds[0]['links'], seeds[0]['edge']
        enabled = []
        for rule, table in enumerate(e.tables):
            for patch in range(len(e.factor.pairs)):
                code = e.factor.oracle.code(initial, rule, patch)
                if table[code] != code:
                    enabled.append((rule, patch, code, table[code]))
        self.assertEqual(len(enabled), 96)
        rule, patch, code, target = enabled[0]
        op = rule*len(e.factor.pairs)+patch
        schedule = [op]*20
        full = e.compiled_bank([initial], schedule, 1, capture_links=True)['runs'][1]
        p.analyze(background, initial, edge, schedule, full, 1)
        wrong = copy.deepcopy(full)
        wrong['events'][0][4] = code
        with self.assertRaisesRegex(ValueError, 'independent primitive replay'):
            p.analyze(background, initial, edge, schedule, wrong, 1)
        wrong = copy.deepcopy(full)
        wrong['events'][1][0] = wrong['events'][0][0]
        with self.assertRaisesRegex(ValueError, 'event tick'):
            p.analyze(background, initial, edge, schedule, wrong, 1)
        # Forge an event on a support that cannot read the initial perturbation.
        distant = next(k for k, read in enumerate(p.reads[1]) if edge not in read)
        wrong = copy.deepcopy(full)
        wrong['events'] = [[1, 1, distant, 0, 0]]
        with self.assertRaisesRegex(ValueError, 'damaged read dependency'):
            p.analyze(background, initial, edge, [len(e.factor.pairs)+distant]*20, wrong, 1)

    def test_dependency_metric_and_periodic_seams_are_from_actual_supports(self):
        p = self.p
        _, seeds = p.initial_conditions()
        edge = seeds[0]['edge']
        self.assertNotIn(edge, p.seam_edges)
        self.assertTrue(p.seam_edges)
        d = distances(p.dependency_neighbors, [edge])
        for reads, writes in zip(p.reads, p.writes):
            for read, write in zip(reads, writes):
                self.assertTrue(all(d[b] <= d[a]+1 for a in read for b in write))
        with self.assertRaisesRegex(ValueError, 'connected'):
            distances([set(), set()], [0])

    def test_seed_has_zero_mean_charge_drift_but_exact_positive_contrast_production(self):
        p, e = self.p, self.p.e
        _, seeds = p.initial_conditions()
        for seed in seeds:
            drift, contrast = [0]*p.faces, 0
            for rule, table in enumerate(e.tables):
                for patch in range(len(e.factor.pairs)):
                    code = e.factor.oracle.code(seed['links'], rule, patch)
                    if table[code] == code:
                        continue
                    self.assertGreater(rule, 0)
                    target = table[code]
                    values = [target//36, target//6 % 6, target % 6]
                    q = [e.charges[h] for h in values]
                    self.assertEqual(sorted(q), [0, 1, 2])
                    for f, value in zip(e.fan_faces[patch][::-1], q):
                        drift[f] += value-1
                        contrast += (value-1)**2
            self.assertEqual(drift, [0]*p.faces)
            self.assertEqual(contrast, 192)

    def test_bounds_and_initial_comparison_reject_invalid_inputs(self):
        for sides, trials in (((5,), 1), ((6, 6), 1), ((6,), 0)):
            with self.assertRaises(ValueError):
                audit(sides, trials)
        with self.assertRaises(ValueError):
            run_case(6, 0, self.bank, 21)
        base, seeds = self.p.initial_conditions()
        with self.assertRaisesRegex(ValueError, 'exactly the seed edge'):
            self.p.analyze(base, base, seeds[0]['edge'], [], {}, 1)


if __name__ == '__main__':
    unittest.main()
