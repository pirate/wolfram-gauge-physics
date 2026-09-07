import sys
import unittest
import itertools
import json
import copy
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'tools'))
from charge_quotient import induced_table, classify_exclusion, classify_block_exchange
from run_diffusion_experiment import fits, moments, validate_recording
from screen_braid_rules import nullspace
from exact_tag_distribution import tag_distribution, binomial_pmf


class EffectiveDynamicsTests(unittest.TestCase):
    def test_recorded_runs_have_consistent_independent_accumulators(self):
        root = Path(__file__).resolve().parents[1]
        for name in ('d4-diffusion.json', 'd4-diffusion-replication.json'):
            result = json.loads((root/'data'/name).read_text())
            validate_recording(result)
        corrupt = copy.deepcopy(result)
        corrupt['ensembles'][0]['block_raw_moment_sums'][0][0][0] += 1
        with self.assertRaisesRegex(ValueError, 'moments disagree'):
            validate_recording(corrupt)

    def test_hidden_state_can_prevent_a_closed_factor(self):
        # The two representatives of symbol zero do not predict the same result.
        labels = [0, 0, 1]
        rule = [2, 1, 0, 3, 4, 5, 6, 7, 8]
        self.assertIsNone(induced_table(rule, labels))

    def test_exact_factor_and_serialized_exclusion_classification(self):
        swap = [3*b+a for a in range(3) for b in range(3)]
        quotient = induced_table(swap, [0, 0, 1])
        self.assertEqual(quotient, [(0, 0), (1, 0), (0, 1), (1, 1)])
        self.assertEqual(classify_exclusion([list(pair) for pair in quotient], [0, 1]), 0)

    def test_block_exchange_classification_exhausts_five_symbol_swap_gates(self):
        # Every symmetric swap/no-swap choice on ten unordered distinct pairs.
        # The braid relation holds iff the non-exchange relation is a partition.
        alphabet = list(range(5))
        pairs = list(itertools.combinations(alphabet, 2))
        accepted = 0
        for choices in itertools.product((False, True), repeat=len(pairs)):
            swaps = {pair for pair, choice in zip(pairs, choices) if choice}
            table = [(b, a) if tuple(sorted((a, b))) in swaps else (a, b)
                     for a in alphabet for b in alphabet]
            braid = True
            for triple in itertools.product(alphabet, repeat=3):
                left, right = list(triple), list(triple)
                for state, order in ((left, (0, 1, 0)), (right, (1, 0, 1))):
                    for i in order: state[i:i+2] = table[state[i]*5+state[i+1]]
                if left != right:
                    braid = False
                    break
            blocks = classify_block_exchange(table, alphabet)
            self.assertEqual(blocks is not None, braid)
            accepted += braid
        self.assertEqual(accepted, 52)  # Bell number B_5.

    def test_exact_prediction_records_finite_time_bias(self):
        times = [0, 16, 32, 64]
        predicted = [tag_distribution(t, .25) for t in times]
        d, exponent = fits(times, predicted)
        self.assertAlmostEqual(d, 1.5, places=2)
        self.assertGreater(exponent, 1)  # Do not confuse a finite-time fit with alpha(infinity).

    def test_conservation_nullspace_rank_and_substitution(self):
        self.assertEqual(nullspace([[1, -1, 0], [0, 1, -1], [2, -2, 0]], 3), [[1, 1, 1]])
        self.assertEqual(nullspace([[1, 0], [0, 1]], 2), [])

    def test_moments_distinguish_drift_from_spreading(self):
        # One block with samples -1,+1; another with samples 1,3.
        result = moments([[[0, 2, 0, 2]], [[4, 10, 28, 82]]], 2)[0]
        self.assertEqual(result['mean'], 1)
        self.assertEqual(result['msd'], 3)
        self.assertEqual(result['variance'], 2)
        self.assertEqual(result['kurtosis'], 2)

    def test_diffusion_units_and_ballistic_control(self):
        times = [0, 4, 8, 16]
        d, alpha = fits(times, [{'variance': t} for t in times])
        self.assertAlmostEqual(d, .5)
        self.assertAlmostEqual(alpha, 1)
        _, alpha = fits(times, [{'variance': t*t} for t in times])
        self.assertAlmostEqual(alpha, 2)

    def test_exact_tag_distribution_against_exhaustive_initial_states(self):
        for rho in (.25, .5, .75):
            histogram = [0.]*5
            offsets = [i for i in range(-4, 5) if i]
            for bits in itertools.product((0, 1), repeat=8):
                weight = rho**sum(bits)*(1-rho)**(8-sum(bits))/2
                for phase in (0, 1):
                    state = [0]*24
                    state[12] = 1
                    for offset, bit in zip(offsets, bits): state[12+offset] = 2*bit
                    for tick in range(2):
                        for i in range((tick+phase) % 2, 23, 2):
                            if not state[i] or not state[i+1]: state[i], state[i+1] = state[i+1], state[i]
                    histogram[state.index(1)-12+2] += weight
            for measured, expected in zip(histogram, tag_distribution(2, rho)['pmf']):
                self.assertAlmostEqual(measured, expected, places=12)

    def test_exact_empty_and_full_controls(self):
        self.assertAlmostEqual(sum(binomial_pmf(1024, .75)), 1)
        self.assertAlmostEqual(tag_distribution(8, 0)['variance'], 64)
        self.assertAlmostEqual(tag_distribution(8, 1)['variance'], 0)


if __name__ == '__main__':
    unittest.main()
