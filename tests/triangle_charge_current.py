import json
import random
import subprocess
import sys
import unittest
from fractions import Fraction
from pathlib import Path
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'tools'))
from triangle_charge_current import CurrentProbe, audit, local_census


class TriangleChargeCurrentTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.source = json.loads((ROOT/'data/triangle-feedback.json').read_text())
        cls.data = json.loads((ROOT/'data/triangle-charge-current.json').read_text())

    def test_all_216_local_moment_identities_reproduce(self):
        rows = local_census(self.source['bank'])
        self.assertEqual(json.loads(json.dumps(rows)), self.data['local_census'])
        self.assertEqual(len(rows), 216)
        active = [r for r in rows if r['charges'] == [1, 1, 1] and r['charge_raw_second_sum'][0][0]]
        self.assertEqual(len(active), 12)
        for r in active:
            self.assertEqual(r['charge_drift_sum'], [0, 0, 0])
            self.assertEqual(r['charge_raw_second_sum'], [[8, -4, -4], [-4, 8, -4], [-4, -4, 8]])
            self.assertEqual(r['current_raw_second_sum'], [[8, 4], [4, 8]])

    def test_whole_current_and_activation_audit_reproduces(self):
        self.assertEqual(json.loads(json.dumps(audit())), self.data)

    def test_arbitrary_raw_connections_match_all_primitive_moments(self):
        for side in (3, 6):
            probe = CurrentProbe(side, self.source['bank'])
            rng = random.Random(71288+side)
            links = [rng.randrange(6) for _ in probe.e.geometry.edges]
            self.assertEqual(probe.analytic(links), probe.enumerate(links))
            self.assertEqual(probe.projected(links, [1]*probe.faces), (0, 0, 13*6*side*side))

    def test_current_paths_cross_exactly_the_written_primal_links(self):
        for side in (3, 6, 12):
            probe = CurrentProbe(side, self.source['bank'])
            self.assertEqual(len(set(probe.dual_primal_edges)), 3*side*side)
            for arity, paths in enumerate(probe.paths):
                for p, path in enumerate(paths):
                    expected = probe.e.factor.oracle.fan.writes[p] if arity else {probe.e.factor.oracle.pairs[p][0][0][0]}
                    self.assertEqual({probe.dual_primal_edges[edge] for edge, _ in path}, expected)

    def test_same_full_charge_field_and_drift_can_have_different_noise(self):
        readout = self.data['readout']
        a, b = readout['same_boundary_witness_states']
        self.assertEqual([a, b], [1, 2])
        left, right = [readout['records'][i]['moments'] for i in (a, b)]
        self.assertEqual(left['charge_drift_sum'], right['charge_drift_sum'])
        self.assertEqual(left['current_drift_sum'], right['current_drift_sum'])
        self.assertNotEqual(left['charge_raw_second_sum'], right['charge_raw_second_sum'])
        faces = self.source['readout']['encounter_faces']
        expected = [[a, b, 8 if a == b else -4] for a in faces for b in faces]
        self.assertEqual(readout['raw_second_sum_difference'], expected)
        self.assertTrue(all(r['exact_same_moments'] for r in readout['gauge_checks']))

    def test_interval_compensators_equal_expanded_attempt_by_attempt_accounting(self):
        for run in self.data['evolution']['runs']:
            denominator = run['denominator']
            drift_sum = variance_sum = residual_squares = raw_squares = raw_expected = 0
            previous = 0
            x = run['initial_region_charge']
            for start, end, before, dx, drift, second in run['constant_state_intervals_start_end_X_deltaX_driftSum_secondSum']:
                self.assertEqual((start, before), (previous, x))
                for tick in range(start+1, end+1):
                    step = dx if tick == end else 0
                    drift_sum += drift
                    variance_sum += denominator*second-drift**2
                    residual_squares += (denominator*step-drift)**2
                    raw_squares += step**2
                    raw_expected += second
                x += dx
                previous = end
            self.assertEqual(previous, run['attempts'])
            self.assertEqual(x, run['final_region_charge'])
            self.assertEqual(drift_sum, run['drift_compensator_numerator'])
            self.assertEqual(variance_sum, run['predictable_variance_numerator_over_denominator_squared'])
            self.assertEqual(residual_squares, run['realized_centered_square_numerator_over_denominator_squared'])
            self.assertEqual(raw_squares, run['realized_raw_squared_changes'])
            self.assertEqual(raw_expected, run['raw_square_compensator_numerator'])
            self.assertEqual(denominator*(x-run['initial_region_charge'])-drift_sum, run['martingale_residual_numerator'])

    def test_homogeneous_activation_and_all_fixed_controls(self):
        activation = self.data['homogeneous_activation']
        self.assertEqual([(r['side'], r['seed']) for r in activation['runs']], [(6, 0), (6, 1), (6, 2), (6, 3), (12, 0)])
        clocks = set()
        for run in activation['runs']:
            probe = CurrentProbe(run['side'], self.source['bank'])
            for initial in (run['initial_frozen_links'], run['initial_active_links']):
                self.assertEqual(probe.charge(initial), [1]*probe.faces)
            clocks.add(Fraction(run['attempts'], 13*6*run['side']**2))
            for control in run['checked']['runs'][:3]:
                self.assertFalse(control['events'])
            self.assertTrue(run['checked']['runs'][3]['events'][0][1])
            for snapshot in run['snapshots']:
                q = probe.charge(snapshot['links'])
                self.assertEqual(q, snapshot['charge_field'])
                n = snapshot['populations_0_1_2']
                self.assertEqual(n, [q.count(i) for i in range(3)])
                self.assertEqual(n[0], n[2])
                self.assertEqual(snapshot['charge_contrast_norm_squared'], n[0]+n[2])
                divergence = [0]*probe.faces
                for edge, flux in snapshot['integrated_dual_currents']:
                    a, b = probe.dual[edge]
                    divergence[a] -= flux; divergence[b] += flux
                self.assertEqual(q, [1+x for x in divergence])
        self.assertEqual(len(clocks), 1)

    def test_uniform_reflection_ensemble_activation_is_derived_not_fitted(self):
        a = self.data['homogeneous_activation']
        self.assertEqual(a['local_reflection_assignments'], 2187)
        self.assertEqual(len(a['based_triple_multiplicities']), 27)
        self.assertEqual({n for _, n in a['based_triple_multiplicities']}, {81})
        self.assertEqual(a['active_assignments'], 972)
        self.assertEqual(Fraction(a['active_assignments'], 2187)*Fraction(12, 13), Fraction(16, 39))
        self.assertEqual(a['initial_ensemble_changing_probability_per_attempt'], [16, 39])
        self.assertEqual(a['initial_ensemble_expected_charge_contrast_increment'], [32, 39])

    def test_raw_capture_preserves_default_protocol_and_rejects_bad_intermediate_links(self):
        e = CurrentProbe(6, self.source['bank']).e
        initial = self.source['evolution']['inputs'][2]
        schedule = self.source['evolution']['runs'][0]['schedule'][:1000]
        plain = e.compiled_bank([initial], schedule, 1000)
        captured = e.compiled_bank([initial], schedule, 1000, capture_links=True)
        for run in captured['runs']:
            self.assertEqual(len(run.pop('raw_event_links')), len(run['events']))
        self.assertEqual(captured, plain)
        actual_run = subprocess.run
        def tamper(*args, **kwargs):
            result = actual_run(*args, **kwargs)
            output = json.loads(result.stdout)
            event = next(run['events'][0] for run in output['runs'] if run['events'])
            event[-1][0] = (event[-1][0]+1) % 6
            result.stdout = json.dumps(output)
            return result
        with patch('triangle_feedback.subprocess.run', side_effect=tamper):
            with self.assertRaisesRegex(ValueError, 'intermediate raw connection'):
                e.compiled_bank([initial], schedule, 1000, capture_links=True)


if __name__ == '__main__':
    unittest.main()
