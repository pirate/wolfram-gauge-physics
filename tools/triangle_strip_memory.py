#!/usr/bin/env python3
"""Same-boundary-sector memory comparison on nested triangular strips.

These are one-cell-wide held-boundary strips, not two-dimensional bulk scaling.
The generator clock gives every contained primitive the same proposal rate.
"""
import itertools
import json
from collections import Counter
from pathlib import Path

import numpy as np
from scipy.sparse import coo_matrix, diags
from scipy.sparse.linalg import LinearOperator, cg, eigsh
from scipy.sparse.csgraph import connected_components

from triangle_elastic_scattering import ElasticExperiment, inverse
from triangle_reference import subgroup


class Strip:
    def __init__(self, length, shape='strip'):
        assert length in (3, 5, 7)
        assert shape in ('strip', 'compact') and (shape == 'strip' or length == 7)
        self.shape = shape
        bank = json.loads(Path('data/triangle-feedback.json').read_text())['bank']
        elastic = json.loads(Path('data/triangle-elastic-scattering.json').read_text())['pair_census']['elastic_tables']
        self.x = x = ElasticExperiment(8, bank, elastic)
        self.g, self.q = x.e.geometry.group, x.e.charges
        self.p = next(i for i, q in enumerate(self.q) if q == 1)
        cells = [2*(2*8+col) for col in range(1, 5)]
        self.faces = [face for a in cells for face in (a+1, a)][:length]
        if shape == 'compact':
            center = 3*8+3
            star = [f for f, vertices in enumerate(x.e.geometry.faces) if center in vertices]
            assert len(star) == 6
            star_edges = {e for f in star for e, _ in x.e.factor.oracle.fan.face_paths[f]}
            extra = next(f for f, path in enumerate(x.e.factor.oracle.fan.face_paths)
                         if f not in star and len({e for e, _ in path} & star_edges) == 1)
            self.faces = star+[extra]
        paths = [x.e.factor.oracle.fan.face_paths[f] for f in self.faces]
        incidence = Counter(e for path in paths for e, _ in path)
        outer = {e for e, c in incidence.items() if c == 1}
        internal = sorted(set(incidence)-outer)
        vertices = {v for e in incidence for v in x.e.geometry.edges[e]}
        boundary_vertices = {v for e in outer for v in x.e.geometry.edges[e]}
        interior_vertices = vertices-boundary_vertices
        assert len(interior_vertices) == int(shape == 'compact')
        self.interior_vertices = sorted(interior_vertices)
        self.fixed_gauge_edges = []
        self.interior_incident = []
        for v in self.interior_vertices:
            incident = [e for e in internal if v in x.e.geometry.edges[e]]
            fixed = next(e for e in incident if any(w in boundary_vertices for w in x.e.geometry.edges[e]))
            self.fixed_gauge_edges.append(fixed)
            self.interior_incident.append(incident)
        internal = [e for e in internal if e not in self.fixed_gauge_edges]
        assert len(internal) == length-1
        directed = [(a, b) for f in self.faces for a, b in zip(x.e.geometry.faces[f], x.e.geometry.faces[f][1:])
                    if x.e.geometry.ids[tuple(sorted((a, b)))] in outer]
        successor = dict(directed)
        root = min(successor)
        boundary = [root]
        for _ in outer:
            boundary.append(successor[boundary[-1]])
        assert boundary[-1] == root and len(set(boundary[:-1])) == len(outer)
        closing = x.e.geometry.ids[tuple(sorted(boundary[-2:]))]
        base = [self.g.identity]*len(x.e.geometry.edges)
        base[closing] = self.p
        self.keys, self.raw, self.charges = [], [], []
        for values in itertools.product(range(self.g.n), repeat=length-1):
            raw = base[:]
            for e, v in zip(internal, values):
                raw[e] = v
            q = tuple(self.q[x.e.factor.oracle.fan.transport(raw, path)] for path in paths)
            if sum(q) == length and len(subgroup(self.g, [self.p, *values])) == self.g.n:
                self.keys.append(values); self.raw.append(raw); self.charges.append(q)
        ids = {v: i for i, v in enumerate(self.keys)}
        partner = [ids[tuple(self.g.conj[self.p][a] for a in v)] for v in self.keys]
        assert all(j != i and partner[j] == i for i, j in enumerate(partner))
        representatives = [i for i, j in enumerate(partner) if i < j]
        self.representatives, self.paths = representatives, paths
        base_ids = {i: k for k, i in enumerate(representatives)}
        project = [base_ids[min(i, j)] for i, j in enumerate(partner)]
        oracle, edges = x.e.factor.oracle, set(incidence)
        pairs = [p for p, paths in enumerate(oracle.pairs) if {e for path in paths for e, _ in path} <= edges]
        fans = [p for p, reads in enumerate(oracle.fan.reads) if reads <= edges]
        self.channels = [(p, r) for p in pairs for r in range(3)]
        self.channels += [(p, r) for p in fans for r in range(3, len(x.tables))]
        assert all(not (oracle.fan.writes[p] & outer) for p in fans)
        assert all(oracle.pairs[p][0][0][0] not in outer for p in pairs)
        tables = []
        for p, r in self.channels:
            row = []
            for i in representatives:
                raw = self.raw[i][:]
                x.apply(raw, r*x.supports+p)
                self.normalize(raw)
                row.append(project[ids[tuple(raw[e] for e in internal)]])
            inverse(row)
            tables.append(row)
        assert all(inverse(t) in tables for t in tables)
        self.tables = np.array(tables)
        self.framed_states = len(self.keys)
        self.charges = np.array([self.charges[i] for i in representatives])
        n, m = len(representatives), len(tables)
        counts = coo_matrix((np.ones(n*m, dtype=np.int64),
                             (np.tile(np.arange(n), m), np.array(tables).reshape(-1))), shape=(n, n)).tocsr()
        assert (counts-counts.T).nnz == 0 and connected_components(counts, directed=False)[0] == 1
        self.laplacian = diags(np.full(n, m), dtype=float)-counts
        patterns = sorted(set(map(tuple, self.charges)))
        labels = {q: i for i, q in enumerate(patterns)}
        self.labels = np.array([labels[tuple(q)] for q in self.charges])
        self.sizes = np.bincount(self.labels)
        self.length, self.n, self.m = length, n, m

    def normalize(self, raw):
        # Only interior vertex frames change. Every boundary link/frame is held.
        for vertex, fixed, incident in zip(self.interior_vertices, self.fixed_gauge_edges,
                                           self.interior_incident):
            h = raw[fixed] if self.x.e.geometry.edges[fixed][0] == vertex else self.g.inv[raw[fixed]]
            for edge in incident:
                a, b = self.x.e.geometry.edges[edge]
                raw[edge] = (self.g.mul[raw[edge]][self.g.inv[h]] if a == vertex
                             else self.g.mul[h][raw[edge]])
            assert raw[fixed] == self.g.identity

    def project(self, v):
        return (np.bincount(self.labels, weights=v, minlength=len(self.sizes))/self.sizes)[self.labels]

    def spectrum(self):
        n, m, lap = self.n, self.m, self.laplacian
        k = min(8, n-2)
        values, vectors = eigsh(lap, k=k, which='SM', tol=1e-11, v0=np.linspace(1, 2, n))
        order = np.argsort(values); values, vectors = values[order], vectors[:, order]
        assert abs(values[0]) < 1e-8
        hidden = lambda v: v-self.project(v)
        def action(v):
            return hidden(lap@hidden(v))+2*m*self.project(v)
        op = LinearOperator((n, n), matvec=action, dtype=float)
        hv, hw = eigsh(op, k=min(k, n-len(self.sizes)), which='SM', tol=1e-11, v0=np.linspace(.3, 1.7, n)**2)
        order = np.argsort(hv); hv, hw = hv[order], hw[:, order]
        residues = [float(np.linalg.norm(self.project(lap@hw[:, i]))**2) for i in range(len(hv))]
        force = lap@(self.charges.astype(float)**2)
        force -= np.column_stack([self.project(force[:, i]) for i in range(self.length)])
        responses = hw.T@force
        force_weights = np.sum(responses**2, axis=1)/np.sum(force**2)
        # The full memory area uses all hidden modes, not just the slow tail.
        # In attempts-per-face units, (1/F) sum_j D^j = 45 (QLQ)^-1.
        solutions = []
        solve_errors = []
        for column in force.T:
            solution, info = cg(op, column, rtol=1e-11, atol=0)
            assert info == 0
            solve_errors.append(float(np.linalg.norm(action(solution)-column)/np.linalg.norm(column)))
            solutions.append(solution)
        inverse_force = np.column_stack(solutions)
        memory_area = float(np.sum(force*inverse_force))
        force_norm = float(np.sum(force**2))
        fields = self.charges.astype(float)**2
        direct_energy = float(np.sum(fields*(lap@fields)))
        slow_mode_areas = 45*force_weights/hv
        visible_ids = [i for i, r in enumerate(residues) if r > 1e-10]
        pair_profile = np.sum(responses[visible_ids[:2]]**2, axis=0)
        if pair_profile.sum() > 1e-20:
            pair_profile /= pair_profile.sum()
        rotation_counts = np.sum(self.charges == 2, axis=1)
        sector_weights = [{str(k): float(np.sum(hw[rotation_counts == k, i]**2))
                           for k in sorted(set(rotation_counts))} for i in visible_ids[:2]]
        errors = [float(np.linalg.norm(lap@vectors[:, i]-values[i]*vectors[:, i])) for i in range(k)]
        errors += [float(np.linalg.norm(action(hw[:, i])-hv[i]*hw[:, i])) for i in range(len(hv))]
        assert max(errors) < 1e-7
        visible = next((i for i, r in enumerate(residues) if r > 1e-10), None)
        assert visible is not None
        lifetime = lambda rate: float(-1/(self.length*np.log1p(-rate/(45*self.length))))
        return {'shape': self.shape, 'interior_vertices': self.interior_vertices,
                'faces': self.faces, 'face_count': self.length, 'framed_states': self.framed_states,
                'unframed_states': n, 'charge_patterns': len(self.sizes), 'contained_operator_slots': m,
                'full_generator_gap': float(values[1]),
                'full_gap_lifetime_attempts_per_face_padded_clock': lifetime(values[1]),
                'hidden_generator_rates': hv.tolist(), 'hidden_charge_coupling_squared': residues,
                'hidden_variance_force_spectral_fractions': force_weights.tolist(),
                'normalized_variance_memory_area_attempts_per_face': 45*memory_area/force_norm,
                'resolved_mode_memory_areas_attempts_per_face': slow_mode_areas.tolist(),
                'variance_force_norm_squared_per_state_per_face': force_norm/(n*self.length),
                'zero_frequency_variance_feedback_fraction_of_direct_energy': memory_area/direct_energy,
                'maximum_memory_solve_relative_residual': max(solve_errors),
                'two_slowest_visible_mode_variance_response_spatial_weights': pair_profile.tolist(),
                'two_slowest_visible_mode_rotation_population_weights': sector_weights,
                'slowest_resolved_charge_memory_rate': float(hv[visible]),
                'memory_lifetime_attempts_per_face_padded_clock': lifetime(hv[visible]),
                'maximum_eigen_residual': max(errors),
                'scope': 'Held-boundary disk, Q=face count, reflection boundary, full S3 image; not a bulk scaling exponent.'}


if __name__ == '__main__':
    for size in (3, 5, 7):
        print(json.dumps(Strip(size).spectrum()), flush=True)
