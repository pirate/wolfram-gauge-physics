#!/usr/bin/env python3
"""Release the former probe into the exact six-face dynamics.

Only the outer six-edge ring is held. The center-vertex frame is gauge fixed;
no internal loop, boundary between the old regions, or probe is held fixed.
"""
import itertools
import json
from collections import Counter

import numpy as np
from scipy.sparse import coo_matrix

from triangle_framed_overlap import FramedOverlap
from triangle_elastic_scattering import inverse
from triangle_reference import subgroup


class CoupledStar:
    def __init__(self):
        f = FramedOverlap()
        f.full_twist()
        self.f, self.x, self.g = f, f.x, f.g
        self.faces = f.order+[f.spectator_face]
        self.paths = [f.oracle.face_paths[i] for i in self.faces]
        self.read_edges = set().union(*({e for e, _ in p} for p in self.paths))
        incidence = Counter(e for p in self.paths for e, _ in p)
        self.outer = {e for e, count in incidence.items() if count == 1}
        self.internal = sorted(self.read_edges-self.outer)
        assert len(self.outer) == len(self.internal) == 6
        assert all(self.x.e.geometry.edges[e][0] == f.u.root for e in self.internal)
        self.fixed, self.variables = self.internal[0], self.internal[1:]
        i = f.ids[(1, 1, 1, 2, 2)]
        before = f.raw[i][:]
        e = next(e for e, _ in f.spectator_path if e not in f.u.edges)
        before[e] = None
        f.solve_edge(before, f.spectator_path, 2)
        after = before[:]
        for k in f.twist_channels:
            p, r = f.channels[k]
            self.x.apply(after, r*self.x.supports+p)
        self.normalize(before)
        self.normalize(after)
        assert all(before[e] == after[e] for e in self.outer)
        assert self.charge(before) == self.charge(after) == (1,)*6
        self.seed_raw = [before, after]
        self.keys, self.raw, self.charges = [], [], []
        for values in itertools.product(range(self.g.n), repeat=len(self.variables)):
            raw = before[:]
            for e, v in zip(self.variables, values):
                raw[e] = v
            word = tuple(f.oracle.transport(raw, p) for p in self.paths)
            if sum(f.q[v] for v in word) == 6 and len(subgroup(self.g, word)) == self.g.n:
                self.keys.append(values)
                self.raw.append(raw)
                self.charges.append(tuple(f.q[v] for v in word))
        self.ids = {key: i for i, key in enumerate(self.keys)}
        self.seeds = [self.ids[self.key(raw)] for raw in self.seed_raw]
        oracle = self.x.e.factor.oracle
        pairs = [p for p, paths in enumerate(oracle.pairs)
                 if {e for path in paths for e, _ in path} <= self.read_edges]
        fans = [p for p, reads in enumerate(oracle.fan.reads) if reads <= self.read_edges]
        self.channels = [(p, r) for p in pairs for r in range(3)]
        self.channels += [(p, r) for p in fans for r in range(3, len(self.x.tables))]
        self.tables = []
        for p, r in self.channels:
            table = []
            for raw in self.raw:
                moved = raw[:]
                self.x.apply(moved, r*self.x.supports+p)
                assert all(moved[e] == raw[e] for e in self.outer)
                self.normalize(moved)
                table.append(self.ids[self.key(moved)])
            inverse(table)
            self.tables.append(table)
        assert all(inverse(t) in self.tables for t in self.tables)
        n, m = len(self.keys), len(self.tables)
        self.counts = coo_matrix((np.ones(n*m, dtype=np.int64),
                                 (np.tile(np.arange(n), m), np.array(self.tables).reshape(-1))), shape=(n, n)).tocsr()
        assert (self.counts-self.counts.T).nnz == 0
        assert np.all(np.asarray(self.counts.sum(axis=1)) == m)

    def normalize(self, raw):
        # U_(0,v) -> U_(0,v) frame_0^{-1}; the boundary vertex frames stay fixed.
        h = self.g.inv[raw[self.fixed]]
        for e in self.internal:
            raw[e] = self.g.mul[raw[e]][h]
        assert raw[self.fixed] == self.g.identity

    def key(self, raw):
        return tuple(raw[e] for e in self.variables)

    def charge(self, raw):
        return tuple(self.f.q[self.f.oracle.transport(raw, p)] for p in self.paths)

    def exact_short_response(self, steps=12):
        # Integer backward powers include every ordinary scheduler word and
        # every no-op. This is not a selected successful interaction path.
        q = np.array(self.charges, dtype=object)
        features = np.concatenate((q, q*q), axis=1)
        rows = []
        for tick in range(steps+1):
            delta = features[self.seeds[1]]-features[self.seeds[0]]
            rows.append({'attempts': tick, 'denominator': len(self.tables)**tick,
                         'charge_mean_difference_numerator': [int(v) for v in delta[:6]],
                         'charge_square_difference_numerator': [int(v) for v in delta[6:]]})
            features = sum((features[t] for t in self.tables), np.zeros_like(features))
        return rows

    def variance_current_identity(self):
        n, clock = len(self.faces), len(self.tables)
        positions = {f: i for i, f in enumerate(self.faces)}
        laplacian = np.zeros((n, n), dtype=np.int64)
        for patch, rule in self.channels:
            if rule != 0:
                continue
            a, b = [positions[f] for f in self.x.p.supports[0][patch]]
            laplacian[a, a] += 1
            laplacian[b, b] += 1
            laplacian[a, b] -= 1
            laplacian[b, a] -= 1
        fans = sorted({p for p, r in self.channels if r >= 3})
        incidence = np.array([[int(f in self.x.e.fan_faces[p]) for p in fans] for f in self.faces])
        q = np.array(self.charges, dtype=np.int64)
        reflections = np.flatnonzero(np.all(q == 1, axis=1))
        first_mean = self.counts@q
        first_square = self.counts@(q*q)
        second_mean = self.counts@first_mean
        d = q-1
        twice_drift = (d*d-d)@laplacian.T
        for patch, rule in self.channels:
            if rule != 0:
                continue
            a, b = [positions[f] for f in self.x.p.supports[0][patch]]
            cubic = d[:, b]*d[:, a]**2-d[:, a]*d[:, b]**2
            twice_drift[:, a] += cubic
            twice_drift[:, b] -= cubic
        for patch in fans:
            fs = [positions[f] for f in self.x.e.fan_faces[patch]]
            mixed = np.all(np.sort(q[:, fs], axis=1) == (0, 1, 2), axis=1)
            twice_drift[:, fs] -= 8*d[:, fs]*mixed[:, None]
        assert np.array_equal(twice_drift, 2*(first_mean-clock*q))
        activity = []
        for i in reflections:
            row = []
            for p in fans:
                code = self.x.e.factor.oracle.code(self.raw[i], 1, p)
                a, b, c = code//36, code//6 % 6, code % 6
                row.append(int((a == b) != (b == c)))
            activity.append(row)
        activity = np.array(activity, dtype=np.int64)
        variance_numerator = first_square[reflections]-clock
        assert np.all(first_mean[reflections] == clock)
        assert np.array_equal(variance_numerator, 8*activity@incidence.T)
        assert np.array_equal(2*(second_mean[reflections]-clock**2), variance_numerator@laplacian.T)
        pair_census = []
        for a, b in itertools.product(range(3), repeat=2):
            d, e = a-1, b-1
            direct = b*int(a == 0)-a*int(b == 0)
            assert 2*direct == d*d-e*e-d+e+e*d*d-d*e*e
            pair_census.append([a, b, direct])
        return {'all_reflection_states_checked': len(reflections),
                'general_first_moment_identity_states_checked': len(q),
                'vacancy_weighted_laplacian': laplacian.tolist(),
                'contained_fan_ids': fans, 'face_fan_incidence': incidence.tolist(),
                'seed_activity_masks': [activity[list(reflections).index(i)].tolist() for i in self.seeds],
                'exact_identity': '2 M^2 (E[q(2)]-1) = L_v (M E[q(1)^2]-M) = 8 L_v B activity',
                'charge_pair_polynomial_census': pair_census}

    def response(self, times=(0, 1, 2, 4, 8, 16, 32, 64, 128, 256, 512)):
        q = np.array(self.charges, dtype=float)
        probability = np.zeros((len(self.keys), 2))
        probability[self.seeds[0], 0] = probability[self.seeds[1], 1] = 1
        result = []
        for tick in range(max(times)+1):
            if tick in times:
                result.append({'attempts': tick, 'total_variation': float(np.abs(probability[:, 1]-probability[:, 0]).sum()/2),
                               'mean_charge': (probability.T@q).tolist(),
                               'mean_charge_squared': (probability.T@(q*q)).tolist(),
                               'expected_rotation_count': (probability.T@(q == 2).sum(axis=1)).tolist()})
            probability = self.counts@probability/len(self.tables)
        return result


def calculation():
    s = CoupledStar()
    from scipy.sparse.csgraph import connected_components
    from fractions import Fraction
    populations = Counter(q.count(2) for q in s.charges)
    stationary_rotation = Fraction(sum(k*v for k, v in populations.items()), len(s.keys))
    return {'gauge_fixed_boundary_framed_states': len(s.keys), 'faces': s.faces,
            'operator_slots': len(s.tables), 'connected_components': connected_components(s.counts, directed=False)[0],
            'seed_state_ids': s.seeds, 'exact_response': s.exact_short_response(),
            'variance_current_identity': s.variance_current_identity(), 'response': s.response(),
            'stationary_rotation_count_populations': dict(populations),
            'stationary_expected_rotation_count': str(stationary_rotation)}


if __name__ == '__main__':
    print(json.dumps(calculation(), indent=2))
