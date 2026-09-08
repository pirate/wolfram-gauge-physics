import itertools
import json
import random
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'tools'))
from fiber_transport import FiberTransport, TransportGeometry, VacancyFactor
from triangle_feedback import FeedbackExperiment, derive_bank
from triple_rule_search import algebra, assess
from triangle_spectral_capacity import capacity_certificate
from fiber_mode_probe import lifted_graph_laplacian
from spectral_transition import integer_inertia


class TriangleFeedbackTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data = json.loads((ROOT/'data/triangle-feedback.json').read_text())
        cls.bank = cls.data['bank']

    def test_independent_complete_minimal_census_and_derived_charge(self):
        g = FiberTransport(3, 3, [(0, 1), (1, 2), (2, 0)]).geometry.group
        observed = derive_bank(g)
        self.assertEqual(json.loads(json.dumps(observed)), self.bank)
        self.assertEqual(observed['minimal_rule_count'], 144)
        self.assertEqual(len(observed['tables']), 12)
        self.assertEqual(observed['charge_basis'], [[1, 2]])
        self.assertEqual(observed['element_charges'], observed['triangle_edge_transposition_lengths'])
        for table in observed['tables']:
            audit = assess(g, table)
            self.assertEqual(audit['moved_states'], 24)
            self.assertIsNotNone(audit['same_classes_and_boundary_feedback'])

    def test_uniform_bank_has_exact_rate_feedback_not_a_selected_rule_effect(self):
        g = FiberTransport(3, 3, [(0, 1), (1, 2), (2, 0)]).geometry.group
        states, products, actions, _ = algebra(g)
        witness = self.bank['uniform_bank_class_boundary_rate_witness']
        a, b = witness['input_codes']
        self.assertEqual(products[a], products[b])
        self.assertEqual(tuple(g.sectors[x] for x in states[a]), tuple(g.sectors[x] for x in states[b]))
        self.assertFalse(any(row[a] == b for row in actions))
        rates = [self.bank['uniform_bank_rate_rows'][i]['outgoing_class_counts'] for i in (a, b)]
        self.assertNotEqual(rates[0], rates[1])
        self.assertEqual(sum(n for _, n in rates[0]), 12)
        self.assertEqual([n for _, n in rates[0]], [2]*6)
        self.assertEqual(rates[1], [[[1, 1, 1], 12]])

    def test_all_single_reflection_subgroups_are_nonreacting_controls(self):
        e = FeedbackExperiment(3, self.bank); g = e.geometry.group
        for reflection in [x for x, q in enumerate(e.charges) if q == 1]:
            for values in itertools.product((g.identity, reflection), repeat=3):
                code = (values[0]*g.n+values[1])*g.n+values[2]
                self.assertTrue(all(table[code] == code for table in e.tables[1:]))

    def test_derived_charge_bounds_full_bundle_spectral_capacity(self):
        e = FeedbackExperiment(3, self.bank)
        certificate = capacity_certificate(e.geometry.group, e.charges)
        self.assertEqual(json.loads(json.dumps(certificate)), self.bank['full_bundle_spectral_capacity'])
        self.assertEqual(len(certificate['face_connection_census']), 216)
        self.assertTrue(certificate['assembly_check']['exact_identity'])
        self.assertTrue(certificate['assembly_check']['both_single_orientation_identities'])
        face_charges = certificate['assembly_check']['face_negative_inertias']
        self.assertEqual(certificate['assembly_check']['charges_by_face_orientation'],
                         [sum(face_charges[c::2]) for c in range(2)])
        self.assertEqual(certificate['flat_full_bundle_upper_edge'], 12)
        for row in certificate['face_connection_census']:
            self.assertEqual(row['inertia']['negative'], e.charges[row['holonomy']])

    def test_full_graph_band_counts_have_independent_exact_polynomial_checks(self):
        from flint import fmpz_mat
        e = FeedbackExperiment(3, self.bank)
        for pair in ((0, 0), (1, 0), (1, 1), (1, 2)):
            links = e.geometry.paired_seed(*pair)[0]
            matrix = lifted_graph_laplacian(range(e.geometry.size), e.geometry.edges, links, e.geometry.group, e.geometry.fiber_edges)
            for i in range(len(matrix)):
                matrix[i][i] -= 12
            inertia = integer_inertia(matrix)
            coefficients = list(fmpz_mat(matrix).charpoly())
            signs = [1 if x > 0 else -1 for x in reversed(coefficients) if x]
            positive = sum(a != b for a, b in zip(signs, signs[1:]))
            self.assertEqual(inertia['positive'], positive)
            charges = [e.charges[x] for x in e.geometry.holonomies(links)]
            self.assertLessEqual(positive, min(sum(charges[c::2]) for c in range(2)))

    def test_braid_memory_readout_and_full_compiled_generator_rates(self):
        e = FeedbackExperiment(12, self.bank)
        r = self.data['readout']
        states = [x['checked']['final_links'] for x in r['preparations']]
        self.assertEqual(len({tuple(e.geometry.sectors(s)) for s in states}), 1)
        self.assertEqual(len({e.forest.signature(s) for s in states}), 4)
        self.assertEqual([len(x['schedule']) for x in r['preparations']], [25]*4)
        self.assertEqual([x['class_changing_triple_operators'] for x in r['global_class_rates']], [12, 12, 0, 0])
        a, b = r['global_rate_witness']['states']
        self.assertEqual(r['encounter_products'][a], r['encounter_products'][b])
        scan = r['global_rate_scan_schedule']
        compiled = e.compiled_bank(states, scan, len(scan))
        self.assertEqual(compiled, r['checked_global_rate_scan'])
        for i, s in enumerate(states):
            rates, summary = e.global_class_rates(s)
            self.assertEqual(summary, r['global_class_rates'][i])
            self.assertEqual(rates, e.rates_from_compiled_scan(s, scan, compiled['runs'][2*i+1]))
            witness = r['global_rate_witness']
            if i in witness['states']:
                self.assertEqual(rates[tuple(witness['output_face_classes'])], witness['counts'][witness['states'].index(i)])

    def test_raw_readout_covariance_and_exact_inverse(self):
        e = FeedbackExperiment(12, self.bank); g = e.geometry.group
        r = self.data['readout']; initial = r['preparations'][0]['checked']['final_links']
        rng = random.Random(49901); frames = [rng.randrange(g.n) for _ in range(e.geometry.size)]
        transform = lambda links: [g.mul[frames[v]][g.mul[x][g.inv[frames[u]]]] for (u, v), x in zip(e.geometry.edges, links)]
        # One forward reaction is enough to test nontrivial gauge covariance;
        # inverse pairs are separately covered by the full readout scan.
        schedule = r['readout_schedule'][:1]
        records = e.compiled_bank([initial, transform(initial)], schedule, 1)['runs']
        self.assertTrue(records[1]['events'])
        self.assertEqual(records[1]['histograms'], records[3]['histograms'])
        self.assertEqual(transform(records[1]['final_links']), records[3]['final_links'])

    def test_unrouted_evolution_replays_all_trials_and_controls(self):
        data = self.data['evolution']; e = FeedbackExperiment(data['side'], self.bank)
        for r in data['runs']:
            checked = e.compiled_bank(data['inputs'], r['schedule'], data['attempts']//100)
            self.assertEqual(checked, r['checked'])
            for run in checked['runs']:
                for histogram in run['histograms']:
                    self.assertEqual(sum(q*n for q, n in zip(e.charges, histogram)), run['conserved_charge'])
                if run['condition'] < 2 or run['mode'] == 'transport':
                    self.assertEqual(run['conversions'], [])
        # These are explicitly bounded negatives, not omitted trials.
        self.assertEqual([r['checked']['runs'][-1]['conversions'] == [] for r in data['runs']], [False, True, True, False])

    def test_generic_cycle_runner_replays_orders_six_through_twelve(self):
        rng = random.Random(223)
        for size in (3, 4, 5, 6):
            geometry = TransportGeometry(3, size, [(i, (i+1) % size) for i in range(size)])
            factor = VacancyFactor(geometry); g = geometry.group
            triples = list(range(g.n**3)); links = [rng.randrange(g.n) for _ in geometry.edges]
            count = len(factor.pairs); schedule = [rng.randrange(2*count) for _ in range(80)]
            protocol = [3, 1, 1, len(schedule), len(schedule)]+factor.tables[0]+triples+links+schedule
            process = subprocess.run(['build/wgphysics_mixed_bank_experiments', '--cycle', str(size)],
                                     input=' '.join(map(str, protocol))+'\n', text=True, capture_output=True, check=True)
            result = json.loads(process.stdout)
            expected = links[:]
            for encoded in schedule:
                rule, patch = divmod(encoded, count)
                factor.oracle.update(expected, rule, patch, triples if rule else factor.tables[0])
            self.assertEqual(result['runs'][0]['final_links'], expected)
            self.assertEqual(result['runs'][1]['final_links'], expected)
            self.assertTrue(all(r['exact_link_inverse'] for r in result['runs']))

    def test_invalid_cycle_and_protocol_fail_before_output(self):
        for runner in ('wgphysics_triple_census', 'wgphysics_mixed_bank_experiments'):
            for args in (['--cycle', '2'], ['--cycle', '7'], ['--cycle', '3oops'], ['--cycle'], ['--unknown', '3']):
                p = subprocess.run(['build/'+runner]+args, input='', text=True, capture_output=True)
                self.assertNotEqual(p.returncode, 0)
                self.assertEqual(p.stdout, '')
        # Census bound is distinct from the runner's supported cycle range.
        p = subprocess.run(['build/wgphysics_triple_census', '--cycle', '5'], text=True, capture_output=True)
        self.assertNotEqual(p.returncode, 0)
        self.assertIn('order <=8', p.stderr)
        e = FeedbackExperiment(3, self.bank)
        with self.assertRaisesRegex(ValueError, 'schedule'):
            e.compiled_bank([[0]*len(e.geometry.edges)], [-1], 1)


if __name__ == '__main__':
    unittest.main()
