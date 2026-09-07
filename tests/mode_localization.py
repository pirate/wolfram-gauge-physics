import copy
import base64
import json
import itertools
import random
import sys
import unittest
import zlib
from pathlib import Path

import numpy as np
from scipy.linalg import eigh

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'tools'))
from fiber_mode_probe import connection_laplacian
from mode_localization import (ModeGeometry, band_report, density, full_profile, positive_rank_bound,
                               sparse_certified_profile, verify_positive_witness)
from run_shared_edge import mesh_geometry
from spectral_transition import (above_band_inertia, curvature_rank_bound,
                                 integer_inertia, triangle_matrix)


class LocalizationTests(unittest.TestCase):
    def test_compact_certificate_has_no_torus_or_induced_boundary_assumption(self):
        from compact_mode_certificate import restriction, verify
        record = json.loads((ROOT/'data/d4-compact-mode-certificate.json').read_text())
        self.assertEqual(verify(record)['above_band_eigenvalue_count'], 2)
        coordinates, operator, bound = restriction()
        g = ModeGeometry(24)
        indices = [2*((y+12)*24+x+12)+a for x, y in coordinates for a in range(2)]
        np.testing.assert_array_equal(operator.toarray(), g.operator(g.seed())[indices][:, indices].toarray())
        self.assertEqual(bound, 2)
        self.assertTrue(np.all(operator.diagonal() == 6))
        for key, value in (('integer_rayleigh_matrix', [['0', '0'], ['0', '0']]),
                           ('strict_gap_lower_bound', {'numerator': 1, 'denominator': 1})):
            bad = copy.deepcopy(record); bad[key] = value
            with self.assertRaises(ValueError):
                verify(bad)

    def test_geometry_patch_order_against_brute_force_incidence(self):
        for side in (3, 4, 6):
            edges, faces, patches = mesh_geometry(side)
            expected = []
            for edge in edges:
                adjacent = [i for i, f in enumerate(faces) if edge in
                            [tuple(sorted(e)) for e in zip(f, f[1:])]]
                self.assertEqual(len(adjacent), 2)
                for first in adjacent:
                    second = next(i for i in adjacent if i != first)
                    u, _ = next(e for e in zip(faces[first], faces[first][1:]) if tuple(sorted(e)) == edge)
                    loops = []
                    for i in (first, second):
                        f = faces[i][:-1]; k = f.index(u)
                        loops.append(f[k:]+f[:k]+(u,))
                    expected.append((*loops, first, second))
            self.assertEqual(patches, expected)

    def test_sparse_operator_matches_independent_integer_construction(self):
        rng = random.Random(103851)
        for side in (3, 4, 6):
            g = ModeGeometry(side)
            links = [rng.randrange(8) for _ in g.edges]
            expected = connection_laplacian(list(range(g.size)), g.edges, links, g.matrices.tolist())
            np.testing.assert_array_equal(g.operator(links).toarray(), expected)
            self.assertEqual(sum(x != 0 for x in g.seed()), 2)

    def test_flat_fourier_spectrum_and_exact_parallel_sections(self):
        for side in (3, 4, 6):
            g = ModeGeometry(side)
            flat = [0]*len(g.edges)
            ks = 2*np.pi*np.arange(side)/side
            expected = sorted(6-2*(np.cos(x)+np.cos(y)+np.cos(x+y)) for x in ks for y in ks for _ in range(2))
            np.testing.assert_allclose(eigh(g.operator(flat).toarray(), eigvals_only=True), expected, atol=1e-12)
            report = full_profile(g, flat)
            self.assertEqual(report['above_flat_band']['rank'], 0)
            self.assertEqual(report['parallel_section_dimension'], 2)
            self.assertAlmostEqual(report['ground']['ipr'], 1/g.size)
            self.assertEqual(full_profile(g, g.seed(True))['parallel_section_dimension'], 1)
            self.assertEqual(full_profile(g, g.seed())['parallel_section_dimension'], 0)

    def test_projector_density_ignores_local_frames_and_band_basis(self):
        g = ModeGeometry(6)
        rng = random.Random(41026)
        links = g.seed()
        frames = [rng.randrange(8) for _ in range(g.size)]
        moved = [g.group.mul[frames[v]][g.group.mul[x][g.group.inv[frames[u]]]]
                 for (u, v), x in zip(g.edges, links)]
        gauge = np.zeros((2*g.size, 2*g.size))
        for v, frame in enumerate(frames):
            gauge[2*v:2*v+2, 2*v:2*v+2] = g.matrices[frame]
        np.testing.assert_array_equal(g.operator(moved).toarray(), gauge@g.operator(links).toarray()@gauge.T)
        for source, changed in ((links, moved), ([0]*len(links),
                [g.group.mul[frames[v]][g.group.inv[frames[u]]] for u, v in g.edges])):
            a, b = full_profile(g, source), full_profile(g, changed)
            for band in ('ground', 'above_flat_band'):
                if a[band]['rank']:
                    np.testing.assert_allclose(a[band]['density'], b[band]['density'], atol=1e-10)
        values, vectors = eigh(g.operator(links).toarray())
        q = vectors[:, values > 9+1e-8]
        rotation = np.array([[3., -4.], [4., 3.]])/5
        np.testing.assert_allclose(density(q), density(gauge@q@rotation), atol=1e-14)
        background = [g.group.mul[frames[v]][g.group.inv[frames[u]]] for u, v in g.edges]
        self.assertEqual(positive_rank_bound(g, links), 2)
        self.assertEqual(positive_rank_bound(g, moved, background), 2)
        self.assertEqual(positive_rank_bound(g, g.seed(True)), 1)
        with self.assertRaisesRegex(ValueError, 'gauge-trivial'):
            positive_rank_bound(g, links, links)

    def test_exact_subspace_certificate_reconstruction_and_corruption(self):
        g = ModeGeometry(12)
        links = g.seed()
        sparse = sparse_certified_profile(g, links)
        dense = full_profile(g, links)
        np.testing.assert_allclose(sparse['above_flat_band']['eigenvalues'], dense['above_flat_band']['eigenvalues'], atol=1e-11)
        np.testing.assert_allclose(sparse['above_flat_band']['density'], dense['above_flat_band']['density'], atol=1e-10)
        certificate = sparse['exact_count_certificate']
        self.assertEqual(verify_positive_witness(g.operator(links), certificate), 2)
        for key, value in (('witness_sha256', '0'*64), ('directions', 1),
                           ('scaled_integer_rayleigh_matrix', [['0', '0'], ['0', '0']])):
            bad = copy.deepcopy(certificate); bad[key] = value
            with self.assertRaises(ValueError):
                verify_positive_witness(g.operator(links), bad)
        with self.assertRaises(ValueError):
            verify_positive_witness(g.operator([0]*len(links)), certificate)

    def test_exact_inertia_handles_singular_and_zero_diagonal_pivots(self):
        rng = np.random.default_rng(51303)
        inputs = [[[0, 1], [1, 0]], [[0, 0], [0, 0]], [[1, 1], [1, 1]],
                  [[0, 2, 1], [2, 0, 3], [1, 3, 0]]]
        for size in range(1, 12):
            for _ in range(8):
                a = rng.integers(-4, 5, size=(size, size))
                inputs.append((a+a.T).tolist())
                # Exact singular PSD and indefinite congruences exercise nullity.
                q = rng.integers(-2, 3, size=(size, max(1, size-2)))
                inputs.append((q@q.T).tolist())
        for a in inputs:
            result = integer_inertia(a)
            values = eigh(np.asarray(a, dtype=float), eigvals_only=True)
            self.assertEqual(result['positive'], int(np.sum(values > 1e-8)))
            self.assertEqual(result['negative'], int(np.sum(values < -1e-8)))
            self.assertEqual(result['zero'], int(np.sum(np.abs(values) <= 1e-8)))
        with self.assertRaises(ValueError):
            integer_inertia([[1, .5], [.5, 1]])

    def test_triangle_holonomy_inertia_and_global_curvature_bound(self):
        g = ModeGeometry(3)
        face = g.faces[0]
        expected = {0: 0, 1: 1, 2: 1, 3: 2, 5: 2}
        for values in itertools.product(range(8), repeat=3):
            links = [0]*len(g.edges)
            for (u, v), value in zip(zip(face, face[1:]), values):
                links[g.ids[tuple(sorted((u, v)))]] = value if u < v else g.group.inv[value]
            sector = g.sectors(links)[0]
            inertia = integer_inertia(triangle_matrix(g, links, face).tolist())
            self.assertEqual(inertia['negative'], expected[sector])
        rng = random.Random(962301)
        for links in ([0]*len(g.edges), g.seed(), [rng.randrange(8) for _ in g.edges]):
            assembled = np.zeros((2*g.size, 2*g.size), dtype=np.int64)
            for f in g.faces:
                indices = [2*v+a for v in f[:3] for a in range(2)]
                assembled[np.ix_(indices, indices)] += triangle_matrix(g, links, f)
            np.testing.assert_array_equal(assembled, 2*(9*np.eye(2*g.size)-g.operator(links).toarray()))
            self.assertLessEqual(above_band_inertia(g, links)['positive'],
                                 curvature_rank_bound(g, links)['above_band_count_upper_bound'])

    def test_saved_large_witnesses_recompute_without_eigensolver(self):
        data = json.loads((ROOT/'data/d4-mode-localization.json').read_text())
        for record in data['static']:
            if 'exact_count_certificate' not in record:
                continue
            g = ModeGeometry(record['side'])
            links = record['initial_links']
            count = verify_positive_witness(g.operator(links), record['exact_count_certificate'])
            self.assertEqual(count, positive_rank_bound(g, links))
            self.assertEqual(count, record['above_flat_band']['rank'])
            cert = record['exact_count_certificate']
            vectors = np.frombuffer(zlib.decompress(base64.b64decode(cert['witness_float64_le_zlib_base64'])), dtype='<f8').reshape(cert['witness_shape'])
            actual = band_report(g.operator(links), np.asarray(record['above_flat_band']['eigenvalues']), vectors, g.core_distances(links))
            np.testing.assert_allclose(actual['density'], record['above_flat_band']['density'], atol=1e-12)
            self.assertAlmostEqual(actual['ipr'], record['above_flat_band']['ipr'], places=12)
            audit = json.loads((ROOT/'data/d4-localization-density-audit.json').read_text())
            measured = next(r for r in audit['resolved_modes'] if r['side'] == g.side)
            for j, mode in enumerate(measured['modes']):
                p = density(vectors[:, j:j+1])
                self.assertAlmostEqual(1/float(p@p), mode['effective_vertices'], places=9)

    def test_saved_primitive_disappearance_has_exact_counts(self):
        from bank_equilibrium import BankFactor
        record = json.loads((ROOT/'data/d4-spectral-transition.json').read_text())
        g = ModeGeometry(record['side'])
        before, after = record['endpoints']
        factor = BankFactor(g.side, json.loads((ROOT/'data/d4-triple-channels.json').read_text()))
        links = before['links'][:]
        _, rule, patch, code, target = record['event']
        self.assertEqual(factor.oracle.update(links, rule, patch, factor.tables[rule]), (code, target))
        self.assertEqual(links, after['links'])
        for endpoint in record['endpoints']:
            self.assertEqual(above_band_inertia(g, endpoint['links']), endpoint['exact_inertia_of_L_minus_9I'])
        self.assertGreater(before['exact_inertia_of_L_minus_9I']['positive'], 0)
        self.assertEqual(after['exact_inertia_of_L_minus_9I']['positive'], 0)


if __name__ == '__main__':
    unittest.main()
