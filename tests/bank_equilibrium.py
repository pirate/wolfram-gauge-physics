import copy
import itertools
import json
import sys
import unittest
from collections import Counter
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'tools'))
from bank_equilibrium import (BankFactor, canonical_counts, canonical_rates, compiled_run,
                              decode, reachable, relaxation, residence, spatial_reference, spatial_residence)
from run_channel_banks import seeds


class EquilibriumTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.channels = json.loads((ROOT/'data/d4-triple-channels.json').read_text())
        cls.factor = BankFactor(3, cls.channels)
        cls.initial = seeds(cls.factor.oracle, 3)[4]
        cls.partial = reachable(cls.factor, cls.initial, 1000)

    def test_canonical_counts_against_independent_small_enumeration(self):
        f = self.factor
        sectors = {}
        for state in itertools.product((0, 1, 2, 3, 5), repeat=5):
            counts = sectors.setdefault(f.signature(state), Counter())
            counts[state.count(f.h)] += 1
        for charges, counts in sectors.items():
            rows = canonical_counts(5, charges)
            self.assertEqual({r['heavy']: r['configurations'] for r in rows}, counts)
            self.assertEqual(sum(Fraction(*r['probability']) for r in rows), 1)
            for row in rows:
                h = row['heavy']
                self.assertEqual(row['weight_per_configuration'], 2**(charges[0]+charges[1]-h))
        rows = canonical_counts(18, (4, 2, 0))
        self.assertEqual([r['configurations'] for r in rows], [278460, 171360, 18360])
        self.assertEqual([r['probability'] for r in rows], [[182, 241], [56, 241], [3, 241]])
        self.assertEqual(canonical_rates(f, rows)['conversions_per_attempt'], [440, 86037])
        # Independent enumeration of whole five-face configurations, observing
        # one fixed fan. This checks the conditional hypergeometric rate formula
        # without drawing iid faces or assuming population-level Markov closure.
        weighted_rate, normalization = Fraction(0), 0
        for state in itertools.product((0, 1, 2, 3, 5), repeat=5):
            if f.signature(state) != (4, 2, 0):
                continue
            weight = f.weight(state)
            normalization += weight
            changing = sum(n for after, n in f.rates[state[:3]].items() if after != state[:3])
            weighted_rate += Fraction(weight*changing, 49)
        self.assertEqual(Fraction(*canonical_rates(f, canonical_counts(5, (4, 2, 0)))['conversions_per_attempt']),
                         weighted_rate/normalization)
        with self.assertRaisesRegex(ValueError, 'negative'):
            canonical_counts(18, (-1, 2, 0))

    def test_family_class_isomorphism_and_artifact_corruption(self):
        normalized = []
        for family in range(3):
            f = BankFactor(3, self.channels, family)
            rename = {0: 0, f.a: 1, f.b: 2, f.h: 3, 5: 4}
            normalized.append({tuple(rename[x] for x in before):
                               {tuple(rename[x] for x in after): n for after, n in outputs.items()}
                               for before, outputs in f.rates.items()})
            self.assertEqual(min(o[x] for x, o in f.rates.items()), 36)
        self.assertEqual(normalized[0], normalized[1])
        self.assertEqual(normalized[1], normalized[2])
        bad = copy.deepcopy(self.channels)
        bad['families'][0]['uniform_bank_class_rates'][0]['outputs'][0][1] -= 1
        with self.assertRaisesRegex(ValueError, 'saved class rates'):
            BankFactor(3, bad)
        bad = copy.deepcopy(self.channels)
        bad['families'][0]['rule_ids'][0] = bad['families'][0]['rule_ids'][1]
        with self.assertRaisesRegex(ValueError, 'distinct'):
            BankFactor(3, bad)

    def test_budget_is_inconclusive_and_cyclic_seed_is_a_real_restriction(self):
        f = self.factor
        self.assertFalse(self.partial['exhausted'])
        self.assertFalse(self.partial['all_formal_states_reachable'])
        self.assertEqual(self.partial['reachable_states'], 1000)
        with self.assertRaisesRegex(ValueError, 'exhaustive component'):
            relaxation(f, self.partial, 100, 20, 20, 1, 7)
        single = seeds(f.oracle, 3)[1]
        result = reachable(f, single, 1000)
        self.assertTrue(result['exhausted'])
        self.assertEqual(result['reachable_counts_by_heavy'], [(2, 153)])
        self.assertFalse(result['all_formal_states_reachable'])
        self.assertEqual(result['component_stationary_probabilities'], [[2, [1, 1]]])
        self.assertGreater(sum(r['configurations'] for r in result['formal_canonical_counts']), 153)
        flat = reachable(f, [0]*len(single), 1)
        self.assertTrue(flat['all_formal_states_reachable'])
        self.assertEqual(flat['reachable_states'], 1)

    def test_paths_lift_to_both_link_oracles_and_arbitrary_local_frames(self):
        f, g = self.factor, self.factor.group
        frames = [(v*3+1) % 8 for v in range(9)]
        def gauge(links):
            return [g.mul[frames[v]][g.mul[x][g.inv[frames[u]]]]
                    for x, (u, v) in zip(links, f.oracle.fan.edges)]
        for witness in self.partial['witnesses']:
            schedule, final = f.lift_path(gauge(self.initial), witness['path'])
            self.assertEqual(schedule, witness['schedule'])
            self.assertEqual(final, gauge(witness['final_links']))
            if schedule:
                run = compiled_run(f, gauge(self.initial), schedule, 1)[1]
                self.assertEqual(run['final_links'], final)
                self.assertTrue(run['independent_link_replay'])
        self.assertEqual([r['heavy'] for r in self.partial['witnesses']], [0, 1, 2])

    def test_exact_residence_clock_includes_null_attempts_and_burn_boundary(self):
        f = self.factor
        witness = next(w for w in self.partial['witnesses'] if w['heavy'] == 0)
        run = compiled_run(f, self.initial, witness['schedule'], 1)[1]
        events = copy.deepcopy(run['events'])
        self.assertEqual(len(events), 2)
        events[0][0], events[1][0] = 10, 20
        measured = residence(f, self.initial, events, 30, 5, 5)
        self.assertEqual(measured['occupation_attempts'], [11, 10, 4])
        self.assertEqual(measured['post_burn_conversions'], 2)
        self.assertEqual(residence(f, self.initial, events, 30, 10, 5)['occupation_attempts'], [11, 9, 0])
        self.assertEqual(residence(f, self.initial, events, 30, 10, 5)['post_burn_conversions'], 1)
        self.assertEqual(residence(f, self.initial, [], 30, 5, 5)['occupation_attempts'], [0, 0, 25])
        with self.assertRaisesRegex(ValueError, 'clock'):
            residence(f, self.initial, events[::-1], 30, 5, 5)
        with self.assertRaisesRegex(ValueError, 'clock'):
            residence(f, self.initial, events, 30, 5, 6)

    def test_spatial_exchangeability_and_event_observer(self):
        f = self.factor
        ref = spatial_reference(f, canonical_counts(18, (4, 2, 0)))
        self.assertEqual(sum(n for _, n in ref['ordered_face_pairs_by_distance']), 18*17)
        self.assertEqual(sum(Fraction(*p) for _, p in ref['expected_a_b_pairs_by_distance']),
                         Fraction(*ref['expected_total_a_b_pairs']))
        per_pair = Fraction(*ref['expected_a_b_per_ordered_face_pair'])
        for (d, n), (other, p) in zip(ref['ordered_face_pairs_by_distance'], ref['expected_a_b_pairs_by_distance']):
            self.assertEqual(d, other)
            self.assertEqual(Fraction(*p)/n, per_pair)
        witness = next(w for w in self.partial['witnesses'] if w['heavy'] == 0)
        run = compiled_run(f, self.initial, witness['schedule'], 1)[1]
        actual = spatial_residence(f, self.initial, run, 2, 0, ref['distances'])
        # Direct raw-link reconstruction at every clock tick, independent of the
        # cached face update used by the pair-distance measurement.
        expected = [0]*len(actual['a_b_pair_attempts_by_distance'])
        links = self.initial[:]
        for encoded in witness['schedule']:
            rule, patch = divmod(encoded, len(f.fans))
            f.oracle.update(links, rule, patch, f.tables[rule])
            classes = f.project(links)
            for a, x in enumerate(classes):
                for b, y in enumerate(classes):
                    if (x, y) == (f.a, f.b):
                        expected[ref['distances'][a][b]] += 1
        self.assertEqual(actual['a_b_pair_attempts_by_distance'], expected)
        self.assertEqual(sum(expected), 3+8)  # N_h=1 then N_h=0.
        bad = copy.deepcopy(run)
        bad['events'][0][3] = 0
        with self.assertRaisesRegex(ValueError, 'raw event input'):
            spatial_residence(f, self.initial, bad, 2, 0, ref['distances'])

    def test_raw_component_weights_are_unchanged_by_positive_generator_reweighting(self):
        f = self.factor
        labels = [tuple(f.group.sectors[x] for x in decode(code, 3)) for code in range(512)]
        unseen = set(range(512))
        coarse_biased = {}
        nonclosure = False
        while unseen:
            start = min(unseen)
            component, queue = {start}, [start]
            for code in queue:
                for table in f.tables[1:]:
                    target = table[code]
                    if target not in component:
                        component.add(target)
                        queue.append(target)
            unseen -= component
            observed = Counter(labels[code] for code in component)
            ratios = {Fraction(count, f.weight(state)) for state, count in observed.items()}
            self.assertEqual(len(ratios), 1)
            incoming = Counter()
            for code in component:
                outcomes = Counter()
                for weight, table in enumerate(f.tables[1:], 1):
                    incoming[table[code]] += weight
                    outcomes[labels[table[code]]] += weight
                if labels[code] in coarse_biased and coarse_biased[labels[code]] != outcomes:
                    nonclosure = True
                coarse_biased[labels[code]] = outcomes
            self.assertEqual(set(incoming.values()), {sum(range(1, 49))})
        # Reweighting changes the observer's memory/kinetics, despite leaving
        # uniform raw-component stationarity invariant.
        self.assertTrue(nonclosure)


if __name__ == '__main__':
    unittest.main()
