import json
import random
import sys
import unittest
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'tools'))
from triangle_reference import (ActivityProbe, AxialReference, audit, constant_connection, field_weight,
                                polynomial_power, reference, representation_audit, sector_audit, subgroup)
from triangle_charge_current import CurrentProbe


class TriangleReferenceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.bank = json.loads((ROOT/'data/triangle-feedback.json').read_text())['bank']
        cls.data = json.loads((ROOT/'data/triangle-reference.json').read_text())

    def test_exact_reference_sampler_and_compiled_controls_reproduce(self):
        self.assertEqual(json.loads(json.dumps(audit())), self.data)

    def test_characters_are_derived_and_complete_and_count_all_handle_pairs(self):
        p = CurrentProbe(3, self.bank)
        observed = representation_audit(p.e.geometry.group, p.e.charges)
        self.assertEqual(observed, self.data['representations'])
        self.assertEqual(observed['dimensions'], [1, 1, 2])
        self.assertEqual(observed['charge_polynomials'], [[1, 3, 2], [1, -3, 2], [1, 0, -1]])
        self.assertEqual(observed['direct_commutator_counts'], [18, 0, 0, 9, 9, 0])

    def test_character_polynomials_match_population_counts_and_topological_parity(self):
        for faces in range(1, 9):
            powers = [polynomial_power(p, faces) for p in self.data['representations']['charge_polynomials']]
            counts = [6*sum(p[q] for p in powers) for q in range(2*faces+1)]
            self.assertEqual(sum(counts), 6**(faces+1))
            self.assertEqual(counts[0], 18)
            for charge, count in enumerate(counts):
                if charge % 2:
                    self.assertEqual(count, 0)
                    with self.assertRaisesRegex(ValueError, 'empty'):
                        reference(faces, charge)
                else:
                    self.assertEqual(reference(faces, charge)['canonical_connections'], count)

    def test_subgroup_and_reflection_invariants_and_direct_four_face_counts(self):
        observed = sector_audit(self.bank)
        self.assertEqual(json.loads(json.dumps(observed)), self.data['sector_checks'])
        self.assertEqual(observed['local_rule_inputs_checked'], 2628)
        self.assertEqual(sorted(map(len, observed['all_subgroups'])), [1, 2, 2, 2, 3, 6])
        self.assertEqual(field_weight(14, 4, 0, 'all'), 972)
        self.assertEqual(field_weight(14, 4, 0, 'nonabelian_reflections'), 960)
        self.assertEqual(field_weight(16, 0, 2, 'nonabelian_reflections'), 0)

    def test_axial_map_is_a_bijection_on_unselected_arbitrary_raw_connections(self):
        for side in (3, 4, 6, 12):
            s = AxialReference(side, self.bank, 4)
            rng = random.Random(234+side)
            for _ in range(12):
                raw = [rng.randrange(6) for _ in s.e.geometry.edges]
                canonical, a, b, h, frames = s.canonicalize(raw)
                self.assertEqual(s.transform(canonical, frames), raw)
                self.assertEqual(s.reconstruct(a, b, h), canonical)
                self.assertEqual([s.e.charges[x] for x in h], s.probe.charge(raw))
                self.assertEqual(frames[0], 0)

    def test_wrap_faces_use_the_correct_basepoint_and_reject_invalid_data(self):
        s = AxialReference(3, self.bank, 4)
        record = s.sample(random.Random(180404))
        canonical = s.reconstruct(*record['handles'], record['cell_based_holonomies'])
        self.assertNotEqual(s.e.geometry.holonomies(canonical), record['cell_based_holonomies'])
        self.assertEqual(s.probe.charge(canonical), [s.e.charges[x] for x in record['cell_based_holonomies']])
        wrong = record['cell_based_holonomies'][:]
        wrong[0] = (wrong[0]+1) % 6
        with self.assertRaisesRegex(ValueError, 'commutator'):
            s.reconstruct(*record['handles'], wrong)
        with self.assertRaisesRegex(ValueError, 'handle'):
            s.reconstruct(-1, 0, record['cell_based_holonomies'])
        with self.assertRaisesRegex(ValueError, 'no torus'):
            AxialReference(3, self.bank, 3)
        with self.assertRaisesRegex(ValueError, 'empty'):
            AxialReference(3, self.bank, 0).sample(random.Random(0), 'nonabelian_reflections')

    def test_conditioned_sector_has_exact_population_and_spatial_predictions(self):
        for faces, expected in ((18, Fraction(9, 59)), (72, Fraction(9, 239))):
            r = reference(faces, 4, 'nonabelian_reflections')
            self.assertEqual(Fraction(*r['mean_rotation_count']), expected)
            self.assertEqual([x['populations'][2] for x in r['population_counts']], [0, 1])
        r = reference(72, 72, 'nonabelian_reflections')
        variance = Fraction(*r['single_face_charge_variance'])
        covariance = Fraction(*r['distinct_face_charge_covariance'])
        self.assertEqual(variance+71*covariance, 0)
        self.assertGreater(variance, 0)

    def test_nonabelian_frozen_states_survive_all_compiled_operators(self):
        frozen = self.data['frozen']
        self.assertEqual(len(frozen['constant_direction_homogeneous_charge_census']), 108)
        self.assertEqual(frozen['frozen_count'], 18)
        self.assertEqual(frozen['nonabelian_frozen_count'], 6)
        for witness in frozen['nonabelian_frozen_witnesses']:
            p = CurrentProbe(witness['side'], self.bank)
            self.assertEqual(p.charge(witness['links']), [1]*p.faces)
            self.assertEqual(len(subgroup(p.e.geometry.group, witness['global_based_loops'])), 6)
            self.assertTrue(witness['all_compiled_operators_fixed'])
            self.assertTrue(all(len(set(t)) == 3 for t in witness['root_star_fan_triples']))
        p = CurrentProbe(6, self.bank)
        raw = constant_connection(p, (1, 2, 5))
        g, rng = p.e.geometry.group, random.Random(99913)
        frames = [rng.randrange(6) for _ in range(p.e.geometry.size)]
        changed_frames = [g.mul[frames[v]][g.mul[x][g.inv[frames[u]]]] for (u, v), x in zip(p.e.geometry.edges, raw)]
        schedule = list(range(13*len(p.e.factor.pairs)))
        checked = p.e.compiled_bank([changed_frames], schedule, len(schedule))
        self.assertTrue(all(not row['events'] and row['final_links'] == changed_frames for row in checked['runs']))

    def test_exact_frozen_star_criterion_and_compiled_activity_degree(self):
        result = self.data['activity']
        self.assertEqual(result['local_pair_inputs_checked'], 36)
        self.assertEqual(result['local_triple_inputs_checked'], 216)
        self.assertEqual(result['reflection_star_inputs_checked'], 729)
        histogram = dict(result['reflection_star_active_fan_histogram'])
        self.assertEqual(sum(histogram.values()), 729)
        self.assertEqual(histogram[0], 3+(2**6+2))
        self.assertTrue(all(k % 2 == 0 for k in histogram))
        p = ActivityProbe(3, self.bank)
        for row in result['compiled_initial_states']:
            self.assertEqual(p.measure(row['links']), row['activity'])
            m = row['activity']
            if 0 < m['charge'] < p.faces:
                self.assertFalse(m['frozen'])
                self.assertGreater(m['vacancy_operators'], 0)
            self.assertEqual(m['frozen'], m['changing_operators'] == 0)

    def test_star_detector_is_gauge_invariant_at_multiple_sizes(self):
        for side in (3, 4, 6, 12):
            p = ActivityProbe(side, self.bank)
            rng, g = random.Random(7651+side), p.e.geometry.group
            for raw in (constant_connection(p, (1, 2, 5)),
                        [rng.randrange(g.n) for _ in p.e.geometry.edges]):
                frames = [rng.randrange(g.n) for _ in range(p.e.geometry.size)]
                transformed = [g.mul[frames[v]][g.mul[h][g.inv[frames[u]]]] for (u, v), h in zip(p.e.geometry.edges, raw)]
                self.assertEqual(p.measure(raw), p.measure(transformed))
                # Independently scan raw tuple targets, with no use of star decomposition.
                changing = sum(table[p.e.factor.oracle.code(raw, rule, patch)] != p.e.factor.oracle.code(raw, rule, patch)
                               for rule, table in enumerate(p.e.tables) for patch in range(len(p.e.factor.pairs)))
                self.assertEqual(p.measure(raw)['changing_operators'], changing)


if __name__ == '__main__':
    unittest.main()
