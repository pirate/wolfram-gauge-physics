import itertools
import json
import random
import sys
import unittest
from collections import Counter
from fractions import Fraction
from pathlib import Path
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'tools'))
from triangle_event_clock import ActiveConnection, block_outcome, geometric_wait, summarize, trial
from triangle_spreading import SpreadingProbe


class IntegerStream:
    def __init__(self, values):
        self.values = iter(values)

    def randrange(self, stop):
        value = next(self.values)
        if not 0 <= value < stop:
            raise AssertionError('scripted integer outside uniform range')
        return value


class TriangleEventClockTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.bank = json.loads((ROOT/'data/triangle-feedback.json').read_text())['bank']
        cls.data = json.loads((ROOT/'data/triangle-first-passage.json').read_text())

    def test_exhaustive_integer_block_masses_are_exactly_geometric(self):
        for total in range(1, 7):
            for active in range(1, total+1):
                for length in range(1, 5):
                    counts = Counter(block_outcome(total, active, length, n) for n in range(total**length))
                    self.assertEqual(counts[None], (total-active)**length)
                    for t in range(1, length+1):
                        self.assertEqual(counts[t], active*(total-active)**(t-1)*total**(length-t))
                    self.assertEqual(sum(counts.values()), total**length)

    def test_multiple_blocks_and_right_censoring_preserve_exact_attempt_clock(self):
        counts = Counter(geometric_wait(IntegerStream([a, b]), 5, 2, 3)
                         for a in range(25) for b in range(5))
        self.assertEqual(counts, {1: 50, 2: 30, 3: 18, None: 27})
        self.assertIsNone(geometric_wait(IntegerStream([]), 5, 0, 10000))
        self.assertIsNone(geometric_wait(IntegerStream([]), 5, 5, 0))
        self.assertEqual(geometric_wait(IntegerStream([4]), 5, 5, 100), 1)
        for args in ((0, 0, 2), (5, 6, 1), (5, 2, -1), (5, True, 1)):
            with self.assertRaises(ValueError):
                geometric_wait(random.Random(0), *args)

    def test_residence_weighting_repairs_jump_chain_degree_bias(self):
        # Two involutions on a three-state path have degrees 1,2,1.
        transition = [[Fraction(1, 2), Fraction(1, 2), 0],
                      [Fraction(1, 2), 0, Fraction(1, 2)],
                      [0, Fraction(1, 2), Fraction(1, 2)]]
        self.assertTrue(all(sum(transition[i][j] for i in range(3)) == 1 for j in range(3)))
        jumps, mean_waits = [Fraction(1, 4), Fraction(1, 2), Fraction(1, 4)], [2, 1, 2]
        residence = [a*b for a, b in zip(jumps, mean_waits)]
        self.assertEqual([x/sum(residence) for x in residence], [Fraction(1, 3)]*3)

    def test_local_cache_matches_all_intermediate_cpp_links_and_full_oracles(self):
        p = SpreadingProbe(6, self.bank)
        e, g = p.e, p.e.geometry.group
        _, seeds = p.initial_conditions()
        rng = random.Random(9117)
        inputs = [seeds[0]['links'], seeds[3]['links']]
        inputs += [[rng.randrange(g.n) for _ in e.geometry.edges] for _ in range(2)]
        schedule = [rng.randrange(13*len(e.factor.pairs)) for _ in range(600)]
        compiled = e.compiled_bank(inputs, schedule, 600, capture_links=True)
        for condition, initial in enumerate(inputs):
            state = ActiveConnection(p, initial)
            expected = compiled['runs'][2*condition+1]
            seen, operators = [], []
            for tick, op in enumerate(schedule, 1):
                rule, patch_id = divmod(op, state.supports)
                code = e.factor.oracle.code(state.links, rule, patch_id)
                self.assertEqual(state.position[op] != -1, code != e.tables[rule][code])
                if state.position[op] == -1:
                    continue
                event, _, refreshed = state.execute(op)
                self.assertLessEqual(refreshed, 33)
                seen.append([tick, *event]); operators.append(op)
                self.assertEqual(state.links, expected['raw_event_links'][len(seen)-1])
                state.full_audit()
            self.assertEqual(seen, expected['events'])
            self.assertEqual(state.links, expected['final_links'])
            for op in reversed(operators):
                state.execute(op)
                state.full_audit()
            self.assertEqual(state.links, initial)

    def test_cache_audit_detects_corruption_and_fixed_states_remain_fixed(self):
        p = SpreadingProbe(6, self.bank)
        state = ActiveConnection(p, [0]*len(p.e.geometry.edges))
        self.assertEqual(state.active, [])
        with self.assertRaisesRegex(ValueError, 'not active'):
            state.execute(0)
        state.codes[0][0] = 1
        with self.assertRaisesRegex(ValueError, 'cached based tuple'):
            state.full_audit()

    def test_cached_sampling_commutes_with_nonconstant_local_gauge_frames(self):
        p = SpreadingProbe(6, self.bank)
        _, seeds = p.initial_conditions()
        g, rng = p.e.geometry.group, random.Random(8866)
        frames = [rng.randrange(g.n) for _ in range(p.e.geometry.size)]
        def transform(links):
            return [g.mul[frames[v]][g.mul[h][g.inv[frames[u]]]] for (u, v), h in zip(p.e.geometry.edges, links)]
        a, b = ActiveConnection(p, seeds[0]['links']), ActiveConnection(p, transform(seeds[0]['links']))
        for _ in range(100):
            self.assertEqual(a.active, b.active)
            op = a.active[rng.randrange(len(a.active))]
            a.execute(op); b.execute(op)
            self.assertEqual(transform(a.links), b.links)
            self.assertEqual(a.q, b.q)
        a.full_audit(); b.full_audit()

    def test_compiled_first_trials_reproduce_checked_in_evidence(self):
        for side in (12, 24):
            p = SpreadingProbe(side, self.bank)
            smaller = next((r for r in self.data['records'] if r['side'] == side//2 and r['trial'] == 0), None)
            observed = trial(p, 0, replay=True, smaller=smaller)
            expected = next(r for r in self.data['records'] if r['side'] == side and r['trial'] == 0)
            self.assertEqual(json.loads(json.dumps(observed)), expected)

    def test_clocks_do_not_change_the_seeded_jump_path_before_boundary_contact(self):
        p = SpreadingProbe(12, self.bank)
        normal = trial(p, 0)
        with patch('triangle_event_clock.geometric_wait', side_effect=lambda rng, m, k, h: 1 if h else None):
            fast = trial(p, 0)
        self.assertEqual(normal['spatial_events_sha256'], fast['spatial_events_sha256'])
        self.assertEqual(normal['final_links'], fast['final_links'])
        self.assertNotEqual(normal['events_sha256'], fast['events_sha256'])
        for (_, a), (_, b) in zip(normal['passages'], fast['passages']):
            if a is not None:
                self.assertEqual(a['event_number'], b['event_number'])
                self.assertEqual(a['conditional_mean_normalized_time'], b['conditional_mean_normalized_time'])

    def test_paired_sizes_have_exact_translated_raw_prefixes_and_holding_means(self):
        records = self.data['records']
        self.assertEqual([(r['side'], r['trial']) for r in records], [(s, i) for s in (12, 24) for i in range(32)])
        for large in records[32:]:
            small = records[large['trial']]
            check = large['smaller_mesh_prefix_check']
            self.assertTrue(check['same_spatial_event_digest'])
            self.assertTrue(check['same_translated_raw_difference'])
            self.assertEqual(check['events'], small['event_count'])
            for (r, a), (_, b) in zip(small['passages'], large['passages']):
                if a is not None:
                    self.assertIsNotNone(b)
                    self.assertEqual(a['event_number'], b['event_number'])
                    self.assertEqual(a['conditional_mean_normalized_time'], b['conditional_mean_normalized_time'])
        self.assertTrue(all(r['stopping'] == 'active_support_reaches_seam' for r in records))

    def test_informative_censoring_does_not_silently_drop_unresolved_trials(self):
        self.assertEqual(summarize(self.data['records'], (2, 4, 6, 8, 10)), self.data['summary'])
        records = [{'side': 12, 'attempted_operators': 100,
                    'passages': [[2, {'tick': 10, 'conditional_mean_normalized_time': [1, 10]}]]},
                   {'side': 12, 'attempted_operators': 100, 'passages': [[2, None]]}]
        row = summarize(records, (2,))[0]
        self.assertEqual((row['observed'], row['censored']), (1, 1))
        self.assertNotIn('mean_normalized_time', row)


if __name__ == '__main__':
    unittest.main()
