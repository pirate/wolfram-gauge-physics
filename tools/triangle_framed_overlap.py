#!/usr/bin/env python3
"""Boundary-framed lift and classical alignment transport of the overlap.

Keep actual outer links and an actual exterior rotation loop fixed. The
residual centralizer sign is then observable, rather than quotiented away.
"""
import itertools
import json
from collections import Counter, deque
from functools import reduce

from flint import fmpq_mat

from triangle_overlap_algebra import build, components
from triangle_reference import subgroup
from triangle_elastic_scattering import inverse


class FramedOverlap:
    def __init__(self):
        self.x, self.u, self.channels, self.quotient_tables, faces, self.boundary = build()
        self.g = self.u.g
        self.q = self.x.e.charges
        self.oracle = self.x.e.factor.oracle.fan
        self.order = [70, 71, 10, 1, 0]
        assert set(self.order) == set(faces)
        self.paths = [self.oracle.face_paths[f] for f in self.order]
        self.p = next(i for i, q in enumerate(self.q) if q == 1)
        self.z = next(i for i, q in enumerate(self.q) if q == 2)
        self.internal = set(self.u.edges)-{edge for edge, _ in self.boundary}
        assert len(self.internal) == 4
        self.product = lambda values: reduce(lambda a, b: self.g.mul[a][b], values, self.g.identity)
        self.words = [w for w in itertools.product(range(self.g.n), repeat=5)
                      if sum(self.q[a] for a in w) == 5 and self.product(w) == self.p
                      and len(subgroup(self.g, w)) == self.g.n]
        assert len(self.words) == 560
        self.ids = {w: i for i, w in enumerate(self.words)}
        self.partner = [self.ids[tuple(self.g.conj[self.p][v] for v in w)] for w in self.words]
        assert all(j != i and self.partner[j] == i for i, j in enumerate(self.partner))
        # Actual exterior triangular loop at the same root, sharing only held
        # outer edges. Its one other edge can be set independently.
        choices = [f for f, path in enumerate(self.oracle.face_paths)
                   if self.x.e.geometry.faces[f][0] == self.u.root and
                   len({edge for edge, _ in path}-set(self.u.edges)) == 1 and
                   not self.internal.intersection(edge for edge, _ in path)]
        assert len(choices) == 1
        self.spectator_face = choices[0]
        self.spectator_path = self.oracle.face_paths[self.spectator_face]
        self.raw = [self.lift(w) for w in self.words]
        self.tables = []
        for p, r in self.channels:
            row = []
            for raw in self.raw:
                moved = raw[:]
                self.x.apply(moved, r*self.x.supports+p)
                assert all(moved[e] == raw[e] for e in range(len(raw)) if e not in self.internal)
                row.append(self.ids[self.read(moved)])
            inverse(row)
            assert all(row[self.partner[i]] == self.partner[row[i]] for i in range(len(row)))
            self.tables.append(row)
        assert len(components(self.tables)) == 1
        self.representatives = [i for i, j in enumerate(self.partner) if i < j]
        self.base_ids = {i: k for k, i in enumerate(self.representatives)}
        self.base = [self.base_ids[min(i, j)] for i, j in enumerate(self.partner)]
        self.orientation = [1 if i < j else -1 for i, j in enumerate(self.partner)]
        self.quotient = [self.u.observe(raw) for raw in self.raw]
        assert len(set(self.quotient)) == 280
        for t, qt in zip(self.tables, self.quotient_tables):
            assert all(self.quotient[t[i]] == qt[self.quotient[i]] for i in range(len(t)))

    def read(self, raw):
        return tuple(self.oracle.transport(raw, path) for path in self.paths)

    def solve_edge(self, raw, path, target):
        unknown = [i for i, (edge, _) in enumerate(path) if raw[edge] is None]
        assert len(unknown) == 1
        k = unknown[0]
        a = self.oracle.transport(raw, path[k+1:])
        b = self.oracle.transport(raw, path[:k])
        value = self.g.mul[self.g.inv[a]][self.g.mul[target][self.g.inv[b]]]
        edge, reverse = path[k]
        raw[edge] = self.g.inv[value] if reverse else value

    def lift(self, word):
        raw = [self.g.identity]*len(self.x.e.geometry.edges)
        edge, reverse = self.boundary[-1]
        raw[edge] = self.g.inv[self.p] if reverse else self.p
        for e in self.internal:
            raw[e] = None
        while any(raw[e] is None for e in self.internal):
            old = sum(raw[e] is None for e in self.internal)
            for path, target in zip(self.paths, word):
                if sum(raw[e] is None for e, _ in path) == 1:
                    self.solve_edge(raw, path, target)
            assert sum(raw[e] is None for e in self.internal) < old
        assert self.read(raw) == word
        spectator_edge = next(e for e, _ in self.spectator_path if e not in self.u.edges)
        raw[spectator_edge] = None
        self.solve_edge(raw, self.spectator_path, self.z)
        assert self.oracle.transport(raw, self.boundary) == self.p
        assert self.oracle.transport(raw, self.spectator_path) == self.z
        return raw

    def anchored_odd(self, word):
        # Gauge-invariant because z, p, and every face loop share an actual root.
        return tuple((self.q[self.product((self.z, self.p, a))]
                      -self.q[self.product((self.g.inv[self.z], self.p, a))])//2 for a in word)

    def full_twist(self):
        reflection = [i for i, w in enumerate(self.words) if all(self.q[v] == 1 for v in w)]
        assert len(reflection) == 80
        selected = []
        for position in range(4):
            expected = []
            for i in reflection:
                w = list(self.words[i])
                a, b = w[position:position+2]
                w[position:position+2] = [self.g.conj[a][b], a]
                expected.append(self.ids[tuple(w)])
            matching = [k for k, ((p, r), t) in enumerate(zip(self.channels, self.tables))
                        if r in (1, 2) and [t[i] for i in reflection] == expected]
            assert matching
            selected.append(matching[0])
        circuit = selected*5
        self.twist_channels = circuit
        for i in reflection:
            raw = self.raw[i][:]
            for k in circuit:
                p, r = self.channels[k]
                self.x.apply(raw, r*self.x.supports+p)
            assert self.read(raw) == self.words[self.partner[i]]
            assert self.anchored_odd(self.read(raw)) == tuple(-a for a in self.anchored_odd(self.words[i]))
            assert any(self.anchored_odd(self.words[i]))
            assert self.oracle.transport(raw, self.spectator_path) == self.z
            assert self.u.observe(raw) == self.quotient[i]
        example = reflection[0]
        return {'reflection_states': 80, 'channel_sequence': [self.channels[k] for k in circuit],
                'length': len(circuit), 'input_word': self.words[example],
                'output_word': self.words[self.partner[example]],
                'external_rotation_face': self.spectator_face,
                'external_rotation_holonomy': self.z, 'held_boundary_holonomy': self.p,
                'anchored_before': self.anchored_odd(self.words[example]),
                'anchored_after': self.anchored_odd(self.words[self.partner[example]]),
                'all_80_raw_circuits_return_projectively_and_invert_external_alignment': True}


    def crossing_readout(self):
        # Separate controlled readout experiment: replace the exterior rotation
        # by a reflection not commuting with p. The interior rules are unchanged.
        probe = next(a for a, q in enumerate(self.q) if q == 1 and a != self.p)
        outside = next(e for e, _ in self.spectator_path if e not in self.u.edges)
        edges = set(self.u.edges) | {outside}
        crossing = [p for p, reads in enumerate(self.oracle.reads)
                    if reads <= edges and not reads <= set(self.u.edges)]
        for i, word in enumerate(self.words):
            if not all(self.q[v] == 1 for v in word):
                continue
            before = self.raw[i][:]
            before[outside] = None
            self.solve_edge(before, self.spectator_path, probe)
            after = before[:]
            for k in self.twist_channels:
                p, r = self.channels[k]
                self.x.apply(after, r*self.x.supports+p)
            assert self.read(after) == self.words[self.partner[i]]
            assert self.x.p.charge(before) == self.x.p.charge(after)
            for p in crossing:
                for r in range(3, len(self.x.tables)):
                    a, b = before[:], after[:]
                    old_a, target_a = self.x.apply(a, r*self.x.supports+p)
                    old_b, target_b = self.x.apply(b, r*self.x.supports+p)
                    qa, qb = self.x.p.charge(a), self.x.p.charge(b)
                    if qa != qb:
                        fs = self.x.p.supports[1][p]
                        distributions = []
                        for state in (before, after):
                            counts = Counter()
                            for rule in range(3, len(self.x.tables)):
                                moved = state[:]
                                self.x.apply(moved, rule*self.x.supports+p)
                                q = self.x.p.charge(moved)
                                counts[tuple(q[f] for f in fs)] += 1
                            distributions.append(counts)
                        if distributions[0] == distributions[1]:
                            continue
                        squared = []
                        for distribution in distributions:
                            counts = Counter()
                            for q, count in distribution.items():
                                counts[sum(a*a for a in q)] += count
                            squared.append(dict(counts))
                        return {'reflection_probe_face': self.spectator_face, 'probe_holonomy': probe,
                                'input_word': word, 'twisted_word': self.words[self.partner[i]],
                                'crossing_channel': [p, r], 'affected_faces_table_order': fs,
                                'primitive_input_codes': [old_a, old_b], 'primitive_output_codes': [target_a, target_b],
                                'charge_outputs': [[q[f] for f in fs] for q in (qa, qb)],
                                'uniform_12_rule_charge_histograms': [[[list(q), count] for q, count in sorted(d.items())] for d in distributions],
                                'uniform_12_rule_squared_charge_histograms': squared,
                                'identical_all_face_charges_before_readout': True,
                                'identical_internal_projective_state_before_readout': True}
        raise ValueError('no boundary-crossing charge readout found')

    def return_alignment(self):
        n, clock = len(self.representatives), len(self.channels)
        plus, minus = [[0]*n for _ in range(n)], [[0]*n for _ in range(n)]
        for i, raw_id in enumerate(self.representatives):
            for t in self.tables:
                j = t[raw_id]
                plus[i][self.base[j]] += 1
                minus[i][self.base[j]] += self.orientation[j]
        matrices = [fmpq_mat(a) for a in (plus, minus)]
        assert all(a == a.transpose() for a in matrices)
        reflection = [i for i, j in enumerate(self.representatives) if all(self.q[a] == 1 for a in self.words[j])]
        mixed = [i for i in range(n) if i not in reflection]
        sub = lambda a, rows, cols: fmpq_mat([[a[i, j] for j in cols] for i in rows])
        identity = fmpq_mat([[clock*int(i == j) for j in mixed] for i in mixed])
        traces, excursions = [], []
        for a in matrices:
            rr, re, er, ee = sub(a, reflection, reflection), sub(a, reflection, mixed), sub(a, mixed, reflection), sub(a, mixed, mixed)
            excursion = re*(identity-ee).inv()*er/clock
            excursions.append(excursion)
            traces.append(rr/clock+excursion)
        even, odd = traces
        for i in range(40):
            for j in range(40):
                assert even[i, j]+odd[i, j] >= 0 and even[i, j]-odd[i, j] >= 0
                assert excursions[0][i, j]+excursions[1][i, j] >= 0 and excursions[0][i, j]-excursions[1][i, j] >= 0
        def diagonal(a):
            return sum(a[i, i] for i in range(40))/40
        summary = {}
        for name, (a, b) in (('first_positive_return', traces), ('mixed_excursion_only', excursions)):
            preserve, invert = (diagonal(a)+diagonal(b))/2, (diagonal(a)-diagonal(b))/2
            summary[name] = {'same_projective_point_same_alignment_probability': float(preserve),
                             'same_projective_point_reversed_alignment_probability': float(invert),
                             'conditional_reversal_among_same_projective_point_returns': float(invert/(preserve+invert)),
                             'same_alignment_probability_exact': str(preserve),
                             'reversed_alignment_probability_exact': str(invert)}
        import numpy as np
        eigenvalues = [np.linalg.eigvalsh(np.array(a, dtype=float)/clock) for a in (plus, minus)]
        summary['numerical_nonconstant_spectral_edges'] = {'even': float(eigenvalues[0][-2]), 'odd': float(eigenvalues[1][-1])}
        summary['odd_mixed_excursion_matrix_nonzero'] = any(excursions[1][i, j] != 0 for i in range(40) for j in range(40))
        return summary


def calculation():
    f = FramedOverlap()
    return {'framed_states': len(f.words), 'unframed_states': len(f.representatives),
            'all_framed_transitions_project_to_previous_unframed_kernel': True,
            'full_twist': f.full_twist(), 'crossing_readout': f.crossing_readout(),
            'boundary_alignment_return': f.return_alignment()}


if __name__ == '__main__':
    print(json.dumps(calculation(), indent=2))
