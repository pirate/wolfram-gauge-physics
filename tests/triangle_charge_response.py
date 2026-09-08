import itertools
import json
import random
import sys
import unittest
from fractions import Fraction
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'tools'))
from triangle_charge_response import (ChargeMarginal, ResponseGeometry, audit, coefficients,
                                      product_polynomial, sparse_formula, sparse_identity)
from triangle_reference import AxialReference
from triangle_response_memory import (audit as memory_audit, charge_targets, raw_targets, response)
from triangle_response_memory import activity_mask, joint_targets


class TriangleChargeResponseTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.bank = json.loads((ROOT/'data/triangle-feedback.json').read_text())['bank']
        cls.data = json.loads((ROOT/'data/triangle-charge-response.json').read_text())
        cls.geometry = ResponseGeometry(3, cls.bank)

    def test_exact_projection_and_exhaustive_census_reproduce(self):
        self.assertEqual(audit(), self.data)
        census = self.data['exhaustive_18_face_charge4_census']
        self.assertEqual(census['ordered_charge_fields'], 5508)
        self.assertEqual(census['canonical_connection_weight'], 3466368)

    def test_ordered_marginals_normalize_and_marginalize_with_raw_weights(self):
        for charge in (2, 4, 18):
            m = ChargeMarginal(18, charge, 'nonabelian_reflections')
            for size in range(5):
                total = 0
                for assignment in itertools.product(range(3), repeat=size):
                    counts = tuple(assignment.count(i) for i in range(3))
                    p = m.probability(counts)
                    self.assertGreaterEqual(p, 0)
                    total += p
                    extended = [m.probability(tuple(n+int(i == j) for j, n in enumerate(counts))) for i in range(3)]
                    self.assertEqual(sum(extended), p)
                self.assertEqual(total, 1)
        # Uniform charge fields would give the wrong mixture of populations.
        m = ChargeMarginal(18, 4, 'nonabelian_reflections')
        self.assertEqual(18*m.probability((0, 0, 1)), Fraction(9, 59))
        self.assertNotEqual(Fraction(2448, 5508), Fraction(9, 59))

    def test_disjoint_drift_terms_have_exact_zero_population_polynomial(self):
        for a, b in ((2, 2), (2, 3), (3, 2), (3, 3)):
            first = (tuple(range(a)), 2 if a == 2 else 4)
            second = (tuple(range(a, a+b)), 2 if b == 2 else 4)
            self.assertEqual(product_polynomial(first, second), ())
        self.assertTrue(product_polynomial(((0, 1), 2), ((0, 1), 2)))

    def test_actual_stencils_match_both_fourier_branches(self):
        for side in (3, 6, 12):
            g = self.geometry if side == 3 else ResponseGeometry(side, self.bank)
            self.assertTrue(np.array_equal(g.l2, 6*g.l1-g.l1@g.l1))
            eigenvalues = []
            for nx, ny in itertools.product(range(side), repeat=2):
                x, y = 2*np.pi*nx/side, 2*np.pi*ny/side
                magnitude = abs(1+np.exp(1j*x)+np.exp(-1j*y))
                eigenvalues.extend((3-magnitude, 3+magnitude))
            np.testing.assert_allclose(np.linalg.eigvalsh(g.l1), sorted(eigenvalues), atol=1e-12)
            self.assertTrue(all(len(terms) == 12 for terms in g.terms))

    def test_two_defect_mean_is_closed_for_every_small_mesh_position(self):
        g = self.geometry
        for occupied in itertools.combinations(range(g.faces), 2):
            q = np.zeros(g.faces, dtype=np.int64)
            q[list(occupied)] = 1
            self.assertEqual(g.drift(q), list(-2*g.l1@q))
        for r in self.data['records']:
            if r['charge'] == 2:
                self.assertEqual(r['nearest_rate'], [2, 1])
                self.assertEqual(r['next_nearest_rate'], [0, 1])
                self.assertEqual(r['closure_residual_covariance_orientation_rows'], [[], []])
                self.assertTrue(r['exact_linear_mean_closure'])

    def test_class_drift_matches_every_raw_operator_and_cpp_inverse(self):
        g = self.geometry
        for charge in (2, 4, 18):
            sampler = AxialReference(3, self.bank, charge)
            links = sampler.sample(random.Random(9082+charge), 'nonabelian_reflections')['links']
            expected = g.drift(g.p.charge(links))
            self.assertEqual(expected, g.p.analytic(links)['charge_drift_sum'])
            self.assertEqual(expected, g.p.enumerate(links, compiled=True)['charge_drift_sum'])

    def test_projection_orthogonality_and_two_step_excess_are_exact(self):
        g, c = self.geometry, self.data['exhaustive_18_face_charge4_census']
        convert = lambda key, denom: np.array([[Fraction(x, c[denom]) for x in row] for row in c[key]], dtype=object)
        s = convert('covariance_integer_numerator', 'covariance_denominator')
        k = convert('negative_drift_charge_cross_integer_numerator', 'cross_denominator')
        h = convert('drift_second_integer_numerator', 'drift_second_denominator')
        _, _, chi, k1, k2 = coefficients(ChargeMarginal(18, 4, 'nonabelian_reflections'))
        gamma = k1*g.l1+k2*g.l2
        self.assertTrue(np.array_equal(gamma@s, k))
        b = h-chi*gamma@gamma
        # This is E[(D+Gamma X)(D+Gamma X)^T], using independent census moments.
        gram = h-k@gamma.T-gamma@k.T+gamma@s@gamma.T
        self.assertTrue(np.array_equal(b, gram))
        m = 39*g.faces
        a = np.eye(g.faces, dtype=object)-gamma/m
        actual = s-2*k/m+h/(m*m)
        self.assertTrue(np.array_equal(actual-a@a@s, b/(m*m)))
        self.assertTrue(all(sum(row) == 0 for row in b))
        self.assertGreater(b[0, 0], 0)
        # Numerical PSD is a diagnostic; the exact weighted Gram identity is
        # the certificate, not an eigenvalue tolerance.
        self.assertGreater(np.linalg.eigvalsh(np.array(b, dtype=float)).min(), -1e-10)

    def test_sparse_rational_identities_and_strict_positivity(self):
        certificate = sparse_identity()
        self.assertEqual(certificate['cleared_identity_coefficients'], [0])
        self.assertTrue(all(x > 0 for x in certificate['defect_numerator_in_F_minus_18']))
        for faces in range(18, 101):
            f = sparse_formula(faces)
            a, b, chi, _, _ = coefficients(ChargeMarginal(faces, 4, 'nonabelian_reflections'))
            self.assertEqual((a, b, chi), (f['a'], f['b'], f['chi']))
            self.assertEqual(f['residual'], f['h']-(12*a*a+132*a*b+378*b*b)/chi)
            self.assertGreater(f['residual'], 0)
        for faces in (0, 17, 18.0, True):
            with self.assertRaises(ValueError):
                sparse_formula(faces)

    def test_all_216_local_charge_transition_multisets(self):
        e = self.geometry.e
        for values in itertools.product(range(6), repeat=3):
            code = 36*values[0]+6*values[1]+values[2]
            q = tuple(e.charges[x] for x in values)
            observed = {}
            for table in e.tables[1:]:
                after = table[code]
                target = tuple(e.charges[x] for x in (after//36, after//6 % 6, after % 6))
                if target != q:
                    observed[target] = observed.get(target, 0)+1
            expected = {}
            if sorted(q) == [0, 1, 2]:
                expected[(1, 1, 1)] = 4
            elif q == (1, 1, 1) and ((values[0] == values[1]) != (values[1] == values[2])):
                expected = {target: 2 for target in itertools.permutations((0, 1, 2))}
            self.assertEqual(observed, expected)

    def test_activity_augmented_generator_matches_actual_raw_charge_targets(self):
        rng, g = random.Random(175930), self.geometry
        for _ in range(3):
            links = [rng.randrange(6) for _ in g.e.geometry.edges]
            self.assertEqual(charge_targets(g, links)[0], raw_targets(g, links, compiled=True))
        # Independently apply L to the existing analytic raw-link drift.
        initial = g.p.analytic(links)['charge_drift_sum']
        second = [0]*g.faces
        for rule, table in enumerate(g.e.tables):
            for patch in range(len(g.e.factor.pairs)):
                code = g.e.factor.oracle.code(links, rule, patch)
                if code == table[code]:
                    continue
                moved = links[:]
                g.e.factor.oracle.update(moved, rule, patch, table)
                after = g.p.analytic(moved)['charge_drift_sum']
                second = [x+a-b for x, a, b in zip(second, after, initial)]
        self.assertEqual(second, response(g, links)['drift_of_drift_sum'])

    def test_two_step_hidden_gauge_response_and_frame_checks_reproduce(self):
        data = json.loads((ROOT/'data/triangle-response-memory.json').read_text())
        self.assertEqual(memory_audit(), data)
        a, b = data['records']
        self.assertEqual(a['charge_field'], b['charge_field'])
        self.assertEqual(a['charge_drift_sum'], b['charge_drift_sum'])
        self.assertNotEqual(a['drift_of_drift_sum'], b['drift_of_drift_sum'])
        self.assertEqual(sum(v for _, v in data['drift_of_drift_difference']), 0)
        self.assertEqual([x-y for x, y in zip(a['two_step_mean_charge_numerator'], b['two_step_mean_charge_numerator'])],
                         [x-y for x, y in zip(a['drift_of_drift_sum'], b['drift_of_drift_sum'])])

    def test_activity_itself_is_not_closed_and_third_mean_response_is_raw_verified(self):
        witness = json.loads((ROOT/'data/triangle-response-memory.json').read_text())['activity_closure_witness']
        g, e = self.geometry, self.geometry.e
        self.assertEqual(witness['sample_indices'], [181, 213])
        self.assertEqual(witness['same_active_reflection_fans'], [])
        self.assertEqual(witness['joint_transition_total_variation'], [2, 351])
        def raw_second(links):
            q = g.p.charge(links)
            initial = g.drift(q)
            result = [0]*g.faces
            for target, count in raw_targets(g, links).items():
                after = g.drift(target)
                result = [x+count*(a-b) for x, a, b in zip(result, after, initial)]
            return result
        third, initial_responses = [], []
        for index, links in enumerate(witness['links']):
            self.assertEqual(g.p.charge(links), witness['same_charge_field'])
            self.assertEqual(activity_mask(g, links), ())
            baseline = raw_second(links)
            initial_responses.append(response(g, links))
            self.assertEqual(baseline, initial_responses[-1]['drift_of_drift_sum'])
            result, cache = [0]*g.faces, {}
            # Independent two-level raw surgery; no activity-based generator is
            # used to compute this third-response numerator.
            for rule, table in enumerate(e.tables):
                for patch in range(len(e.factor.pairs)):
                    code = e.factor.oracle.code(links, rule, patch)
                    if code == table[code]:
                        continue
                    moved = links[:]
                    e.factor.oracle.update(moved, rule, patch, table)
                    key = tuple(moved)
                    if key not in cache:
                        cache[key] = raw_second(moved)
                    result = [x+a-b for x, a, b in zip(result, cache[key], baseline)]
            third.append(result)
            rng, group = random.Random(878+index), e.geometry.group
            frames = [rng.randrange(group.n) for _ in range(e.geometry.size)]
            transformed = [group.mul[frames[v]][group.mul[x][group.inv[frames[u]]]]
                           for (u, v), x in zip(e.geometry.edges, links)]
            self.assertEqual(joint_targets(g, links, compiled=True), joint_targets(g, transformed, compiled=True))
        self.assertEqual(initial_responses[0], initial_responses[1])
        self.assertEqual([a-b for a, b in zip(*third)], witness['three_step_mean_difference_numerator'])
        self.assertEqual(witness['three_step_mean_difference_denominator'], 702**3)
        self.assertTrue(any(witness['three_step_mean_difference_numerator']))


if __name__ == '__main__':
    unittest.main()
