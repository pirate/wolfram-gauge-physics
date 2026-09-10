#!/usr/bin/env python3
"""Joint exact affine evolution of two inherited holonomy constraints.

No particles are labeled and no independent-carrier closure is assumed.
Initial all-reflection charge fields are matched across contact controls.
"""
import itertools
import json
import math
from pathlib import Path

import numpy as np
from flint import fmpz_mat

from fiber_constraint_dynamics import canonical, inverse, local_rule, loop_row, path_word, product
from fiber_refinement_limit import geometry, kind, primitive, read
from fiber_reaction_bursts import write_fan


def all_reflection_signs(oracle, fe, seed):
    rng = np.random.default_rng(seed)
    bits = rng.integers(2, size=len(oracle.edges), dtype=np.int64)
    incident = [[] for _ in oracle.edges]
    for f, edges in enumerate(fe):
        for e in edges:
            incident[e].append(f)
    adjacent = [[] for _ in fe]
    for e, (a, b) in enumerate(incident):
        adjacent[a].append((b, e)); adjacent[b].append((a, e))
    parent, order = {0: (0, -1)}, [0]
    for a in order:
        for b, e in adjacent[a]:
            if b not in parent:
                parent[b] = a, e; order.append(b)
    mismatch = np.bitwise_xor.reduce(bits[fe], axis=1)^1
    for f in reversed(order[1:]):
        if mismatch[f]:
            a, e = parent[f]
            bits[e] ^= 1; mismatch[a] ^= 1; mismatch[f] ^= 1
    assert not np.any(mismatch)
    return 1-2*bits, adjacent


def parameterize(rows):
    rank, E = len(rows), len(rows[0])
    if rank == 1:
        pivots = (next(j for j, x in enumerate(rows[0]) if abs(x) == 1),)
        inverse_pivot = [[rows[0][pivots[0]]]]
    elif rank == 2:
        pivots = next((i, j) for i, j in itertools.combinations(range(E), 2)
                      if abs(rows[0][i]*rows[1][j]-rows[0][j]*rows[1][i]) == 1)
        i, j = pivots
        det = rows[0][i]*rows[1][j]-rows[0][j]*rows[1][i]
        inverse_pivot = [[rows[1][j]//det, -rows[0][j]//det],
                         [-rows[1][i]//det, rows[0][i]//det]]
    else:
        support = [e for e in range(E) if any(row[e] for row in rows)]
        pivots = next(columns for columns in itertools.combinations(support, rank)
                      if abs(fmpz_mat([[row[e] for e in columns] for row in rows]).det()) == 1)
        pivot = fmpz_mat([[row[e] for e in pivots] for row in rows]).inv()
        inverse_pivot = [[int(x) for x in row] for row in pivot.tolist()]
    free = [i for i in range(E) if i not in pivots]
    forms = [[0]*(len(free)+1) for _ in range(E)]
    for k, e in enumerate(free):
        forms[e][k] = 1
        for i, p in enumerate(pivots):
            forms[p][k] = -sum(inverse_pivot[i][j]*rows[j][e] for j in range(rank))
    return list(map(tuple, forms))


class AffineProcess:
    def __init__(self, condition, side, seed):
        self.oracle, self.args = geometry(side)
        self.pe, self.pr, self.fe, self.fr, self.supports = self.args
        self.E, self.F = len(self.oracle.edges), len(self.fe)
        signs, adjacent = all_reflection_signs(self.oracle, self.fe, seed)
        distance = {int(f): 0 for f in self.supports[0][:2]}
        queue = list(distance)
        for f in queue:
            for g, _ in adjacent[f]:
                if g not in distance:
                    distance[g] = distance[f]+1; queue.append(g)
        far = max(range(len(self.pe)), key=lambda p: min(distance[int(f)] for f in self.supports[p][:2]))
        first_pair = set(map(int, self.supports[0][:2]))
        near = min((p for p in range(len(self.pe))
                    if first_pair.isdisjoint(map(int, self.supports[p][:2]))),
                   key=lambda p: min(distance[int(f)] for f in self.supports[p][:2]))
        selections = {'single_left': [(0, 0)], 'single_right': [(0, 1)],
                      'contact': [(0, 0), (0, 1)], 'separated': [(0, 0), (far, 0)],
                      'near_disjoint': [(0, 0), (near, 0)]}
        self.selection = selections[condition]
        self.rank = len(self.selection)
        rows = []
        for p, i in self.selection:
            a = loop_row(signs, self.pe[p, 2-i], self.pr[p, 2-i])[1]
            b = loop_row(signs, self.pe[p, 1-i], self.pr[p, 1-i])[1]
            rows.append(tuple(x-y for x, y in zip(a, b)))
        forms = parameterize(rows)
        self.raw = [(int(s), v) for s, v in zip(signs, forms)]
        self.initial_rows = rows
        self.initial_signs = signs.tolist()
        self.initial_pair_faces = [self.supports[p][i:i+2].tolist() for p, i in self.selection]
        self.events = []
        assert self.charges() == [1]*self.F
        self.Q = self.F

    def word(self, p):
        return tuple(path_word(self.raw, self.pe[p, 2-i], self.pr[p, 2-i]) for i in range(3))

    def charges(self):
        holonomies = [path_word(self.raw, e, r) for e, r in zip(self.fe, self.fr)]
        return [1 if h[0] == -1 else 0 if not any(h[1]) else 2 for h in holonomies]

    def features(self):
        if self.rank not in (1, 2):
            raise ValueError('These coefficient-norm diagnostics are defined only for ranks one and two')
        signs = [h[0] for h in self.raw]
        q = self.charges()
        assert sum(q) == self.Q and q.count(0) == q.count(2) and q.count(0) <= self.rank
        matrix = fmpz_mat([list(h[1][:-1]) for h in self.raw])
        kernel, nullity = matrix.transpose().nullspace()
        assert nullity == self.rank
        normals = [[int(kernel[e, j]) for e in range(self.E)] for j in range(nullity)]
        offsets = [sum(x*h[1][-1] for x, h in zip(row, self.raw)) for row in normals]
        # Primitive projective Pluecker coordinates are independent of a
        # basis chosen for the two-dimensional constraint row space.
        if self.rank == 2:
            extended = [row+[-b] for row, b in zip(normals, offsets)]
            pluecker = [extended[0][i]*extended[1][j]-extended[0][j]*extended[1][i]
                        for i, j in itertools.combinations(range(self.E+1), 2)]
            divisor = math.gcd(*pluecker)
            complexity = sum(abs(x)//divisor for x in pluecker)
        else:
            divisor = math.gcd(*(normals[0]+[offsets[0]]))
            complexity = sum(abs(x)//divisor for x in normals[0]+[offsets[0]])
        active_rows, visible_rows = set(), set()
        for f in range(self.F):
            if q[f] == 0:
                visible_rows.add(canonical(loop_row(signs, self.fe[f], self.fr[f])[1])[0])
        total = parallel = 0
        for p, fs in enumerate(self.supports):
            w = self.word(p)
            rows = [loop_row(signs, self.pe[p, 2-i], self.pr[p, 2-i])[1] for i in range(3)]
            for i in (0, 1):
                pair = product(w[i], w[i+1])
                if pair[0] == 1 and not any(pair[1]):
                    visible_rows.add(canonical(tuple(a+w[i][0]*b for a, b in zip(rows[i], rows[i+1])))[0])
            qq = tuple(q[f] for f in fs)
            if qq == (1, 1, 1):
                ab = w[0][1] == w[1][1]
                bc = w[1][1] == w[2][1]
                parallel += int(ab and bc)
                if ab != bc:
                    total += 12
                    i = 0 if ab else 1
                    active_rows.add(canonical(tuple(a-b for a, b in zip(rows[i], rows[i+1])))[0])
            elif sorted(qq) == [0, 1, 2]:
                total += 4
                active_rows.add(canonical(rows[qq.index(0)])[0])
        active_rank = fmpz_mat(list(active_rows)).rank() if active_rows else 0
        visible_rank = fmpz_mat(list(visible_rows)).rank() if visible_rows else 0
        zero_rank = self.rank-int(any(offsets))
        assert active_rank <= visible_rank <= zero_rank
        for row in normals:
            divergence = [0]*(max(max(e) for e in self.oracle.edges)+1)
            for x, s, (u, v) in zip(row, signs, self.oracle.edges):
                divergence[u] -= s*x; divergence[v] += x
            assert not any(divergence)
        return [total, q.count(0), active_rank, visible_rank, zero_rank, parallel,
                sum(any(row[e] for row in normals) for e in range(self.E)), complexity,
                int(total == 0), len(self.events)]

    def step(self, p, channel, t, shadow, parameters, modulus):
        w = self.word(p)
        target = local_rule(w, channel, lambda v: not any(v))
        agreed = True
        if shadow is not None:
            actual = tuple(read(shadow, self.pe[p, 2-i], self.pr[p, 2-i], modulus) for i in range(3))
            predicted = tuple((sum(c*int(x) for c, x in zip(h[1][:-1], parameters))+h[1][-1]) % modulus
                              +(modulus if h[0] == -1 else 0) for h in target)
            actual_target = primitive(*actual, channel, modulus)
            agreed = predicted == actual_target
            if actual_target != actual:
                write_fan(shadow, p, *actual_target, modulus, self.pe, self.pr)
        if target == w:
            return agreed
        aa, bb, cc = target
        prefix = path_word(self.raw, self.pe[p, 0, :2], self.pr[p, 0, :2])
        rim = self.raw[self.pe[p, 1, 1]]
        if self.pr[p, 1, 1]:
            rim = inverse(rim)
        s2 = product(prefix, inverse(cc))
        s3 = product(product(rim, prefix), inverse(product(bb, cc)))
        self.raw[self.pe[p, 0, 2]] = s2 if self.pr[p, 0, 2] else inverse(s2)
        self.raw[self.pe[p, 1, 2]] = s3 if self.pr[p, 1, 2] else inverse(s3)
        if 3 <= channel < 15:
            self.events.append((t, p, int(all(h[0] == -1 for h in w))))
        return agreed


def experiment(condition, alpha=1, side=4, trials=32, seed=914003001):
    times = (0., .01, .05, .2, .5, 1.)
    observations, histories, mismatches, proposals = [], [], [], []
    initial_geometry = None
    for trial in range(trials):
        process = AffineProcess(condition, side, seed+trial)
        if initial_geometry is None:
            initial_geometry = {'pairs': process.initial_pair_faces, 'signs': process.initial_signs,
                                'constraint_rows': process.initial_rows}
        modulus = 1_000_000_007
        rng = np.random.default_rng(seed+10_000+trial)
        parameters = rng.integers(modulus, size=process.E-process.rank, dtype=np.int64)
        shadow = np.array([(sum(c*int(x) for c, x in zip(h[1][:-1], parameters))+h[1][-1]) % modulus
                           +(modulus if h[0] == -1 else 0) for h in process.raw], dtype=np.int64)
        clock = np.random.default_rng(seed+20_000+trial)
        rate = len(process.pe)*(15+2*alpha)
        t, tick, attempts, mismatch, row = 0., 0, 0, None, []
        while tick < len(times):
            t += clock.exponential(1/rate)
            while tick < len(times) and times[tick] <= t:
                row.append(process.features())
                if shadow is not None:
                    q = [kind(read(shadow, process.fe[f], process.fr[f], modulus), modulus) for f in range(process.F)]
                    if q != process.charges():
                        mismatch, shadow = times[tick], None
                tick += 1
            if tick == len(times):
                break
            p = int(clock.integers(len(process.pe)))
            slot = int(clock.integers(15+2*alpha))
            channel = slot if slot < 15 else 15+(slot-15)//alpha
            attempts += 1
            if not process.step(p, channel, t, shadow, parameters, modulus) and shadow is not None:
                mismatch, shadow = t, None
        observations.append(row); histories.append(process.events); mismatches.append(mismatch); proposals.append(attempts)
    x = np.array(observations, dtype=float)
    return {'condition': condition, 'alpha': alpha, 'side': side, 'trials': trials, 'seed': seed,
            'initial_charge': 2*side**2, 'initial_preparation': 'Uniform raw sign field conditional on all-reflection faces; generic joint zero-constraint plane',
            'times': times, 'initial_geometry': initial_geometry,
            'columns': ['reaction_rate', 'identity_count', 'active_constraint_rank', 'visible_flat_rank',
                        'zero_offset_rank', 'parallel_RRR_fans', 'joint_support_edges', 'primitive_extended_Pluecker_L1',
                        'reaction_inactive', 'reaction_count'],
            'means': x.mean(axis=0).tolist(), 'independent_trajectory_SE': (x.std(axis=0, ddof=1)/np.sqrt(trials)).tolist(),
            'observations_by_trajectory': observations, 'reaction_histories': histories,
            'finite_shadow_modulus': modulus, 'first_shadow_mismatches': mismatches, 'proposal_counts': proposals}


if __name__ == '__main__':
    result = {'scope': 'Controlled joint-constraint encounters; not a stationary two-event Palm ensemble', 'runs': []}
    for condition, alpha in (('single_left', 1), ('single_right', 1), ('contact', 1), ('separated', 1), ('contact', 0), ('separated', 0)):
        row = experiment(condition, alpha)
        result['runs'].append(row)
        Path('data/fiber-two-constraints.json').write_text(json.dumps(result, indent=2)+'\n')
        print(json.dumps({'condition': condition, 'alpha': alpha, 'initial': row['means'][0],
                          'final': row['means'][-1], 'shadow_mismatches': sum(x is not None for x in row['first_shadow_mismatches'])}), flush=True)
