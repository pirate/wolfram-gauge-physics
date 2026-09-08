import copy
import json
import sys
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'tools'))
from fiber_transport import TransportGeometry
from fiber_mode_probe import lifted_graph_laplacian
from triangle_modes import (PAIRS, exact_band_count, exact_forms, finite_study, operator,
                            perturbation, reaction_audit, restriction, verify_compact)


class TriangleModeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.certificate = json.loads((ROOT/'data/triangle-compact-modes.json').read_text())
        cls.data = json.loads((ROOT/'data/triangle-modes.json').read_text())

    def test_exact_infinite_certificates_do_not_call_an_eigensolver(self):
        with patch('triangle_modes.eigh', side_effect=AssertionError('not a verifier')), patch('triangle_modes.eigsh', side_effect=AssertionError('not a verifier')):
            result = verify_compact(self.certificate)
        self.assertEqual(result, self.certificate['verification'])
        self.assertEqual([r['exact_above_band_count'] for r in result], [0, 1, 1, 2])
        for r in result[1:]:
            self.assertEqual(r['strict_lower_edge'], '12 + 1/20')

    def test_certificate_rejects_tampering_and_overflow(self):
        for field in ('integer_rayleigh_matrix', 'integer_mass_matrix', 'perturbation_matrix'):
            record = copy.deepcopy(self.certificate)
            record['cases'][-1][field][0][0] += 1
            with self.assertRaises(ValueError):
                verify_compact(record)
        record = copy.deepcopy(self.certificate)
        record['cases'][-1]['integer_trial_columns'][0][0] = 10**100
        with self.assertRaisesRegex(ValueError, 'overflow'):
            verify_compact(record)
        record = copy.deepcopy(self.certificate); record['cases'].pop()
        with self.assertRaisesRegex(ValueError, 'four'):
            verify_compact(record)
        record = copy.deepcopy(self.certificate); record['verification'][-1]['exact_above_band_count'] = 3
        with self.assertRaisesRegex(ValueError, 'summary'):
            verify_compact(record)

    def test_full_graph_operator_matches_independent_lift_including_wrapped_links(self):
        for side in (3, 6, 16):
            geometry = TransportGeometry(side, 3, [(0, 1), (1, 2), (2, 0)])
            coords = [(x, y) for y in range(side) for x in range(side)]
            for pair in (*PAIRS, (3, 4)):
                links = geometry.paired_seed(*pair)[0]
                independent = lifted_graph_laplacian(range(geometry.size), geometry.edges, links, geometry.group, geometry.fiber_edges)
                np.testing.assert_array_equal(operator(coords, pair, period=side).toarray(), independent)

    def test_principal_restriction_keeps_infinite_boundary_degree(self):
        coords, laplacian = restriction(6, (0, 0))
        self.assertEqual(laplacian.shape, (507, 507))
        np.testing.assert_array_equal(laplacian.diagonal(), np.full(507, 8))
        self.assertTrue(np.any(np.asarray(laplacian.sum(axis=1)).ravel() > 0))
        self.assertEqual(int(laplacian[3*coords.index((0, 0))].sum()), 0)
        with self.assertRaisesRegex(ValueError, 'radius'):
            restriction(1, (1, 2))

    def test_perturbation_is_supported_on_three_endpoints_and_has_exact_positive_index(self):
        coords, flat = restriction(6, (0, 0))
        ids = [3*coords.index((x, 0))+a for x in range(3) for a in range(3)]
        for pair, expected in zip(PAIRS, (0, 1, 1, 2)):
            _, laplacian = restriction(6, pair)
            difference = laplacian-flat
            self.assertLessEqual(set(difference.nonzero()[0]) | set(difference.nonzero()[1]), set(ids))
            matrix, inertia = perturbation(pair)
            np.testing.assert_array_equal(difference[ids][:, ids].toarray(), matrix)
            self.assertEqual(inertia['positive'], expected)

    def test_compact_integer_forms_embed_unchanged_in_large_tori(self):
        local, _ = restriction(6, (1, 2))
        case = self.certificate['cases'][-1]
        for side in (16, 24, 48):
            coords = [(x, y) for y in range(side) for x in range(side)]
            laplacian = operator(coords, case['seed_pair'], period=side)
            q = np.zeros((3*side*side, 2), dtype=np.int64)
            for i, (x, y) in enumerate(local):
                v = ((y+side//2) % side)*side+(x+side//2) % side
                q[3*v:3*v+3] = case['integer_trial_columns'][3*i:3*i+3]
            gram, mass = exact_forms(laplacian, q.tolist(), 2)
            self.assertEqual(gram, case['integer_rayleigh_matrix'])
            self.assertEqual(mass, case['integer_mass_matrix'])

    def test_numerical_profiles_are_separate_from_exact_mode_counts(self):
        observed = finite_study(self.certificate, [16])
        for row, saved in zip(observed, self.data['static'][:4]):
            self.assertEqual(row['certified_above_band_count'], saved['certified_above_band_count'])
            for mode, reference in zip(row['numerical_modes'], saved['numerical_modes']):
                self.assertAlmostEqual(mode['eigenvalue_estimate'], reference['eigenvalue_estimate'], places=8)
                self.assertAlmostEqual(mode['base_participation_volume'], reference['base_participation_volume'], places=5)
        for row in self.data['static']:
            for mode in row['numerical_modes']:
                if 'base_projector_density' in mode:
                    density = mode['base_projector_density']
                    self.assertAlmostEqual(sum(density), 1, places=10)
                    self.assertAlmostEqual(1/sum(p*p for p in density), mode['base_participation_volume'], places=8)
        with self.assertRaisesRegex(ValueError, 'room'):
            finite_study(self.certificate, [12])

    def test_exact_dynamic_counts_reproduce_all_changing_events_and_reaction_endpoints(self):
        observed = reaction_audit()
        self.assertEqual(json.loads(json.dumps(observed)), self.data['dynamics'])
        self.assertEqual(len(observed['reactions']), 12)

    def test_orientation_capacity_forces_each_observed_reaction_loss(self):
        dynamics = self.data['dynamics']
        losses = []
        for r in dynamics['reactions']:
            counts = [x['positive'] for x in r['exact_inertias']]
            for n, charges in zip(counts, r['orientation_charges_before_after']):
                self.assertLessEqual(n, min(charges))
                self.assertEqual(sum(charges), 4)
            if counts[1] < counts[0]:
                losses.append(r['event'][0])
                self.assertTrue(r['decrease_forced_by_capacity'])
        self.assertEqual(losses, [545, 1975, 12579, 18296, 12937])
        for _, _, rank, up, down in dynamics['seed_zero_full_changing_event_counts']:
            self.assertEqual(up+down, 4)
            self.assertLessEqual(rank, min(up, down))

    def test_complete_attempt_interval_accounting_does_not_claim_particle_lifetimes(self):
        dynamics = self.data['dynamics']
        intervals = dynamics['seed_zero_constant_rank_intervals']
        self.assertEqual(intervals[0]['start'], 0)
        self.assertEqual(intervals[-1]['end'], 100000)
        self.assertTrue(all(a['end'] == b['start'] for a, b in zip(intervals, intervals[1:])))
        self.assertEqual(sum(dynamics['seed_zero_attempt_interval_residence'].values()), 100000)
        self.assertEqual(dynamics['seed_zero_attempt_interval_residence']['0'], 26940)
        self.assertTrue(any(rank == 0 and up > 0 and down > 0 for _, _, rank, up, down in dynamics['seed_zero_full_changing_event_counts']))

    def test_exact_polynomial_counter_requires_real_symmetric_integer_input(self):
        with self.assertRaisesRegex(ValueError, 'symmetric'):
            exact_band_count([[0, 1], [0, 0]])
        with self.assertRaisesRegex(ValueError, 'integer'):
            exact_band_count([[0.1]])
        self.assertEqual(exact_band_count([[11, 0, 0], [0, 12, 0], [0, 0, 13]]), {'positive': 1, 'negative': 1, 'zero': 1})


if __name__ == '__main__':
    unittest.main()
