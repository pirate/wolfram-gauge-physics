import copy
import itertools
import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'tools'))
from face_energy_obstruction import derive_fiber_group
from triple_channels import analyze, positive_certificate
from triple_rule_search import make_table
from run_channel_banks import MixedOracle, experiment, verify_run


class ChannelTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.census = json.loads((ROOT/'data/d4-triple-rule-search.json').read_text())
        cls.channels = analyze(cls.census)

    def test_exact_positivity_and_channel_coverage(self):
        c = self.channels
        self.assertEqual(len(c['rules']), 945)
        self.assertEqual(sum(r['positive_certificate']['exists'] for r in c['rules']), 755)
        self.assertEqual(c['eligible_rules'], 144)
        self.assertEqual([len(f['rule_ids']) for f in c['families']], [48]*3)
        for r in c['rules']:
            cert = r['positive_certificate']
            row = cert.get('equation', cert.get('nonnegative_obstruction'))
            feasible = any(sum(a*b for a, b in zip(row, q)) == 0 for q in itertools.product(range(1, 4), repeat=4))
            self.assertEqual(cert['exists'], feasible)
        with self.assertRaisesRegex(ValueError, 'independent equation'):
            positive_certificate([(1, -1, 0), (0, 1, -1)], 3)
        with self.assertRaisesRegex(ValueError, 'dimensions'):
            positive_certificate([(1, -1)], 3)
        corrupt = copy.deepcopy(self.census)
        corrupt['minimal_rules'].pop()
        with self.assertRaisesRegex(ValueError, 'independent minimal closures'):
            analyze(corrupt)

    def test_charge_cones_and_sublattice_release(self):
        for f in self.channels['families']:
            self.assertEqual(len(f['charge_basis']), 3)
            self.assertEqual(len(f['triple_only_staggered_basis']), 4)
            self.assertEqual(len(f['shared_edge_transport_staggered_basis']), 3)
            for row in f['shared_edge_transport_staggered_basis']:
                self.assertEqual(row[:4], row[4:])
            self.assertEqual(f['reaction_species']['grand_canonical_occupation_ratio'], [4, 2])
        self.assertEqual(len(self.channels['family_incompatibilities']), 3)
        self.assertEqual(sorted(x['forced_zero_class_weights'] for x in self.channels['family_incompatibilities']), [[1], [2], [3]])

    def test_uniform_bank_factor_and_entropic_balance(self):
        g = derive_fiber_group()
        for f in self.channels['families']:
            self.assertIsNone(f['uniform_bank_class_rate_nonclosure'])
            rows = {tuple(r['input']): dict((tuple(out), count) for out, count in r['outputs']) for r in f['uniform_bank_class_rates']}
            self.assertEqual(len(rows), 125)
            size = lambda state: 2**sum(s in (1, 2, 3) for s in state)
            for state, outcomes in rows.items():
                self.assertEqual(sum(outcomes.values()), 48)
                for after, count in outcomes.items():
                    self.assertEqual(size(state)*count, size(after)*rows[after].get(state, 0))
            h = f['reaction_species']['heavy_class']
            a, b = f['reaction_species']['light_classes']
            grand_weights = {0: 1, h: 12, a: 4, b: 6, 5: 5}
            product_weight = lambda state: grand_weights[state[0]]*grand_weights[state[1]]*grand_weights[state[2]]
            for state, outcomes in rows.items():
                for after, count in outcomes.items():
                    self.assertEqual(product_weight(state)*count, product_weight(after)*rows[after].get(state, 0))
            tables = [make_table(self.channels['rules'][i]['transpositions']) for i in f['rule_ids']]
            for code in range(512):
                before = (code//64, (code//8) % 8, code % 8)
                observed = tuple(g.sectors[a] for a in before)
                actual = {}
                for table in tables:
                    out = table[code]
                    after = tuple(g.sectors[a] for a in (out//64, (out//8) % 8, out % 8))
                    actual[after] = actual.get(after, 0)+1
                self.assertEqual(actual, rows[observed])

    def test_raw_link_runs_controls_and_corruption(self):
        result = experiment(self.channels, 3, 2000, 100, 1, 6291871)
        self.assertEqual(len(result['runs']), 42)
        self.assertTrue(all(r['independent_link_replay'] for r in result['runs']))
        self.assertTrue(any(r['conversion_events'] for r in result['runs']))
        self.assertTrue(all(r['conversion_events'] == 0 for r in result['runs'] if r['condition'] in ('flat', 'r', 's', 'rotation')))
        bad = copy.deepcopy(next(r for r in result['runs'] if r['events']))
        bad['events'][0][3] ^= 1
        f = self.channels['families'][bad['family']]
        tables = [self.channels['shared_edge_vacancy_transport']]+[make_table(self.channels['rules'][i]['transpositions']) for i in f['rule_ids']]
        initial = result['initial_links'][result['conditions'].index(bad['condition'])]
        with self.assertRaisesRegex(ValueError, 'event log'):
            verify_run(MixedOracle(3, derive_fiber_group()), initial, result['schedules'][0], tables, f, bad, 100)


if __name__ == '__main__':
    unittest.main()
