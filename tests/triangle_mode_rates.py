import copy
import json
import sys
import unittest
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'tools'))
from triangle_mode_rates import audit, candidates, saturation_audit


class TriangleModeRateTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data = json.loads((ROOT/'data/triangle-mode-rates.json').read_text())
        cls.source = json.loads((ROOT/'data/triangle-feedback.json').read_text())
        cls.modes = json.loads((ROOT/'data/triangle-modes.json').read_text())

    def test_every_endpoint_and_full_compiled_kernel_reproduces(self):
        self.assertEqual(json.loads(json.dumps(audit())), self.data)

    def test_selection_covers_all_reactions_without_a_spectral_filter(self):
        selected = candidates(self.source, self.modes)
        self.assertEqual(len(selected), 28)
        self.assertEqual(selected, [{k: v for k, v in row.items() if k != 'kernel'}
                                    for row in self.data['candidates']])

    def test_all_noops_and_both_arities_are_in_probability_denominator(self):
        for row in self.data['candidates']:
            k = row['kernel']
            self.assertEqual(k['attempted_operators'], 2808)
            self.assertEqual(sum(n for _, n in k['rank_counts']), 2808)
            for arity, outcomes in enumerate(k['outcomes_transport_reaction']):
                self.assertEqual(sum(n for _, n in outcomes), 216*(12 if arity else 1))
                self.assertEqual(sum(n for state, n in outcomes if state[0] < k['observable'][0]),
                                 k['loss_operators_transport_reaction'][arity])
                self.assertEqual(sum(n for state, n in outcomes if state[0] > k['observable'][0]),
                                 k['gain_operators_transport_reaction'][arity])
            self.assertTrue(k['all_unique_targets_compiled_checked'])
            self.assertTrue(k['full_operator_inverse_scan_checked'])

    def test_nonclosure_is_inside_one_reachable_component_and_not_a_gauge_artifact(self):
        i, j = self.data['macrostate_nonclosure_witness']
        self.assertEqual([i, j], [3, 7])
        a, b = [self.data['candidates'][x]['kernel'] for x in (i, j)]
        self.assertEqual(a['observable'], [2, 2, 2, 4, 0])
        self.assertEqual(a['observable'], b['observable'])
        self.assertEqual(a['rank_counts'], [[1, 16], [2, 2792]])
        self.assertEqual(b['rank_counts'], [[1, 28], [2, 2780]])
        self.assertEqual([x['candidate_index'] for x in self.data['gauge_checks']], [i, j])
        self.assertTrue(all(x['exact_same_kernel'] for x in self.data['gauge_checks']))

    def test_saturation_bound_and_charge_two_inert_reaction_census(self):
        observed = saturation_audit(self.source['bank'])
        self.assertEqual(json.loads(json.dumps(observed)), self.data['saturation_obstruction'])
        self.assertEqual(len(observed['vacancy_charge_transfers']), 10)
        for row in self.data['candidates']:
            k = row['kernel']
            if k['observable'][0] == k['observable'][1] == k['observable'][2]:
                checks = k['saturated_vacancy_checks_operator_charge_count_bound']
                self.assertEqual(len(checks), k['changing_operators_transport_reaction'][0])
                self.assertEqual(k['loss_operators_transport_reaction'][0], len(checks))
                for _, charge, count, bound in checks:
                    self.assertGreater(charge, 0)
                    self.assertEqual(bound, k['observable'][0]-charge)
                    self.assertLessEqual(count, bound)

    def test_low_charge_obstruction_rejects_an_inserted_reaction(self):
        bad = copy.deepcopy(self.source['bank'])
        bad['tables'][0][0] = 1
        with self.assertRaisesRegex(ValueError, 'charge at most two'):
            saturation_audit(bad)

    def test_single_reflection_first_loss_clock_includes_noops(self):
        k = self.data['candidates'][1]['kernel']
        self.assertEqual(k['observable'], [1, 1, 1, 2, 0])
        self.assertEqual(k['rank_counts'], [[0, 8], [1, 2800]])
        self.assertEqual(k['changing_operators_transport_reaction'], [8, 0])
        self.assertEqual(Fraction(k['attempted_operators'], 8), 351)

    def test_one_step_stable_control_has_exact_two_step_escape(self):
        control = self.data['candidates'][2]['kernel']
        self.assertEqual(control['rank_counts'], [[1, 2808]])
        escape = self.data['one_step_stability_is_not_protection']
        self.assertEqual(escape['shortest_escape_length'], 2)
        self.assertEqual(escape['first_loss_by_two_attempts'], [1, 164268])
        self.assertEqual(escape['first_loss_by_two_changing_updates'], [13, 80])
        rows = escape['changed_first_step_targets']
        self.assertEqual(len(rows), 8)
        first = sum(len(row['first_operators']) for row in rows)
        self.assertEqual(first, 16)
        counted = sum(len(row['first_operators'])*sum(row['second_step_kernel']['loss_operators_transport_reaction']) for row in rows)
        self.assertEqual(counted, escape['losing_ordered_operator_pairs'])
        self.assertEqual(counted, 48)
        jump_probability = sum((Fraction(len(row['first_operators']), first)*Fraction(
            sum(row['second_step_kernel']['loss_operators_transport_reaction']),
            sum(row['second_step_kernel']['changing_operators_transport_reaction'])) for row in rows), Fraction(0))
        self.assertEqual(jump_probability, Fraction(13, 80))
        witness = escape['escape_witness']
        self.assertEqual(witness['final_observable'], [0, 0, 4, 4, 0])
        self.assertTrue(witness['compiled_exact_inverse'])


if __name__ == '__main__':
    unittest.main()
