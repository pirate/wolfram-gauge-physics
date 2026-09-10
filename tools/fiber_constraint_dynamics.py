#!/usr/bin/env python3
"""Exact one-constraint, fine-fiber Palm dynamics with integer arithmetic.

State: raw link signs, primitive covector ell, and integer offset b,
meaning ell dot a = b. Every update is derived from the original local
group words. A finite-modulus shadow checks the same event schedule.
"""
import json
import math
from pathlib import Path

import numpy as np

from fiber_refinement_limit import LAYOUTS, geometry, kind, primitive, read
from fiber_reaction_bursts import write_fan


def canonical(row, offset=0):
    divisor = math.gcd(*row)
    assert divisor and offset % divisor == 0
    if next(x for x in row if x) < 0:
        divisor = -divisor
    return tuple(x//divisor for x in row), offset//divisor


def product(a, b):
    return a[0]*b[0], tuple(x+a[0]*y for x, y in zip(a[1], b[1]))


def inverse(a):
    return a[0], tuple(-a[0]*x for x in a[1])


def path_word(elements, edges, reverse):
    h = (1, (0,)*len(elements[0][1]))
    for e, rev in zip(edges, reverse):
        h = product(inverse(elements[e]) if rev else elements[e], h)
    return h


def loop_row(signs, edges, reverse):
    row = [0]*len(signs)
    parity = 1
    for e, rev in zip(edges, reverse):
        s = int(signs[e])
        row = [s*x for x in row]
        row[e] += -s if rev else 1
        parity *= s
    return parity, tuple(row)


def local_rule(word, channel, zero):
    a, b, c = word
    ident = (1, (0,)*len(a[1]))
    equal = lambda x, y: x[0] == y[0] and zero(tuple(u-v for u, v in zip(x[1], y[1])))
    charge = lambda x: 1 if x[0] == -1 else 0 if zero(x[1]) else 2
    ka, kb, kc = map(charge, word)
    if channel == 0:
        return (b, a, c) if ka == 0 or kb == 0 else word
    if channel == 1:
        return product(product(a, b), inverse(a)), a, c
    if channel == 2:
        return b, product(product(inverse(b), a), b), c
    if channel >= 15:
        if (ka, kb) not in ((1, 2), (2, 1)):
            return word
        z = a if ka == 2 else b
        v = list(z[1])
        v[-1] += 1 if channel == 15 else -1
        target = (1, tuple(v))
        if zero(target[1]):
            return word
        P = product(a, b)
        return (product(P, inverse(target)), target, c) if ka == 1 else (target, product(inverse(target), P), c)
    index = channel-3
    layout, orientation = tuple(LAYOUTS[index//2]), index % 2 == 0

    def reverse(w):
        return tuple(inverse(x) for x in w[::-1])

    def forward(r, s):
        z = product(r, s) if orientation else product(s, r)
        zi, ri = layout.index(2), layout.index(1)
        h = product(inverse(z), s) if zi < ri else product(s, inverse(z))
        entries = (ident, h, z)
        return tuple(entries[k] for k in layout)

    def backward(w):
        s = product(product(w[0], w[1]), w[2])
        z = w[layout.index(2)]
        r = product(z, s) if orientation else product(s, z)
        return r, r, s

    if ka == kb == kc == 1:
        ab, bc = equal(a, b), equal(b, c)
        if ab and not bc:
            return forward(a, c)
        if bc and not ab:
            return reverse(forward(c, a))
    elif (ka, kb, kc) == layout:
        return backward(word)
    elif (ka, kb, kc) == layout[::-1]:
        return reverse(backward(reverse(word)))
    return word


class ConstraintProcess:
    def __init__(self, side, seed):
        self.oracle, self.args = geometry(side)
        self.pe, self.pr, self.fe, self.fr, self.supports = self.args
        self.E, self.F = len(self.oracle.edges), len(self.oracle.faces)
        rng = np.random.default_rng(seed)
        self.signs = rng.choice((-1, 1), size=self.E).astype(int)
        if rng.integers(2):
            self.initial_direction = 'backward'
            kinds = (1, 1, 1)
            pair = int(rng.integers(2))
        else:
            self.initial_direction = 'forward'
            kinds = tuple(LAYOUTS[rng.integers(6)])
            pair = None
        for t in range(3):
            e0, e1, e2 = self.pe[0, t]
            self.signs[e1] = self.signs[e0]*self.signs[e2]*(-1 if kinds[2-t] == 1 else 1)
        rows = [loop_row(self.signs, self.pe[0, 2-i], self.pr[0, 2-i])[1] for i in range(3)]
        normal = tuple(a-b for a, b in zip(rows[pair], rows[pair+1])) if pair is not None else rows[kinds.index(0)]
        self.ell, self.b = canonical(normal)
        self.initial_signs, self.initial_ell = self.signs.copy(), self.ell
        self.initial_q = self.charges()
        self.Q = sum(self.initial_q)
        self.events = []
        distance = {}
        for coefficient, edge in zip(self.ell, self.oracle.edges):
            if coefficient:
                for vertex in edge:
                    distance[vertex] = 0
        adjacent = {}
        for u, v in self.oracle.edges:
            adjacent.setdefault(u, []).append(v)
            adjacent.setdefault(v, []).append(u)
        queue = list(distance)
        for u in queue:
            for v in adjacent[u]:
                if v not in distance:
                    distance[v] = distance[u]+1
                    queue.append(v)
        self.edge_radii = [max(distance[u], distance[v]) for u, v in self.oracle.edges]
        self.neighborhoods = []
        for p in range(len(self.pe)):
            edges = sorted(self.oracle.reads[p])
            ids = {e: i for i, e in enumerate(edges)}
            local_pe = np.array([[ids[e] for e in face] for face in self.pe[p]])
            basis = [tuple(int(i == j) for j in range(len(edges)+1)) for i in range(len(edges))]
            self.neighborhoods.append((edges, ids, local_pe, basis))

    def charges(self):
        out = []
        for e, r in zip(self.fe, self.fr):
            s, row = loop_row(self.signs, e, r)
            out.append(1 if s == -1 else 0 if self.b == 0 and canonical(row)[0] == self.ell else 2)
        return out

    def features(self):
        q = self.charges()
        assert sum(q) == self.Q and q.count(0) <= 1
        divergence = [0]*(max(max(e) for e in self.oracle.edges)+1)
        for coefficient, s, (u, v) in zip(self.ell, self.signs, self.oracle.edges):
            divergence[u] -= int(s)*coefficient
            divergence[v] += coefficient
        assert not any(divergence), 'constraint must remain a twisted closed circulation'
        total = 0
        pairs = set()
        for p, fs in enumerate(self.supports):
            qq = tuple(q[f] for f in fs)
            if sorted(qq) == [0, 1, 2]:
                total += 4
            elif qq == (1, 1, 1) and self.b == 0:
                rows = [loop_row(self.signs, self.pe[p, 2-i], self.pr[p, 2-i])[1] for i in range(3)]
                active = []
                for i in (0, 1):
                    pair = canonical(tuple(a-b for a, b in zip(rows[i], rows[i+1])))[0]
                    active.append(pair == self.ell)
                    if active[-1]:
                        pairs.add(tuple(sorted((int(fs[i]), int(fs[i+1])))))
                assert not all(active)
                if any(active):
                    total += 12
        assert len(pairs) <= 1
        assert not (q.count(0) and pairs)
        return [total, q.count(0), len(pairs), sum(x != 0 for x in self.ell),
                sum(abs(x) for x in self.ell), max(abs(x) for x in self.ell), abs(self.b),
                int(self.b != 0), int(total == 0), len(self.events),
                max(r for r, x in zip(self.edge_radii, self.ell) if x)]

    def finite_start(self, n, seed):
        rng = np.random.default_rng(seed)
        pivot = next(i for i, x in enumerate(self.ell) if x)
        while True:
            raw = rng.integers(n, size=self.E, dtype=np.int64)
            other = sum(a*int(x) for i, (a, x) in enumerate(zip(self.ell, raw)) if i != pivot)
            raw[pivot] = ((self.b-other)*pow(self.ell[pivot], -1, n)) % n
            raw += n*(self.signs == -1)
            w = tuple(read(raw, self.pe[0, 2-i], self.pr[0, 2-i], n) for i in range(3))
            if self.initial_direction == 'backward':
                valid = (w[0] == w[1]) != (w[1] == w[2])
            else:
                valid = sorted(kind(h, n) for h in w) == [0, 1, 2]
            if valid:
                return raw

    def step(self, p, channel, t, shadow=None, modulus=None):
        edges, ids, pe, basis = self.neighborhoods[p]
        elements = [(int(self.signs[e]), v) for e, v in zip(edges, basis)]
        support_inside = all(not v or e in ids for e, v in enumerate(self.ell))
        pivot = next((i for i, e in enumerate(edges) if self.ell[e]), None)

        def zero(v):
            if not any(v[:-1]):
                return v[-1] == 0
            if not support_inside or pivot is None:
                return False
            den, num = self.ell[edges[pivot]], v[pivot]
            return v[-1]*den+num*self.b == 0 and all(x*den == num*self.ell[e] for x, e in zip(v[:-1], edges))

        w = tuple(path_word(elements, pe[2-i], self.pr[p, 2-i]) for i in range(3))
        target = local_rule(w, channel, zero)
        agreed = True
        if shadow is not None:
            actual_w = tuple(read(shadow, self.pe[p, 2-i], self.pr[p, 2-i], modulus) for i in range(3))
            actual_target = primitive(*actual_w, channel, modulus)
            predicted = tuple((sum(c*int(shadow[e] % modulus) for c, e in zip(h[1][:-1], edges))+h[1][-1]) % modulus
                              +(modulus if h[0] == -1 else 0) for h in target)
            agreed = predicted == actual_target
            if actual_target != actual_w:
                write_fan(shadow, p, *actual_target, modulus, self.pe, self.pr)
        changed = any(x[0] != y[0] or not zero(tuple(a-b for a, b in zip(x[1], y[1]))) for x, y in zip(w, target))
        if not changed:
            return agreed
        aa, bb, cc = target
        prefix = path_word(elements, pe[0, :2], self.pr[p, 0, :2])
        rim1 = elements[pe[1, 1]]
        if self.pr[p, 1, 1]:
            rim1 = inverse(rim1)
        s2 = product(prefix, inverse(cc))
        s3 = product(product(rim1, prefix), inverse(product(bb, cc)))
        written = [s2 if self.pr[p, 0, 2] else inverse(s2), s3 if self.pr[p, 1, 2] else inverse(s3)]
        W = [int(self.pe[p, i, 2]) for i in range(2)]
        new_signs = self.signs.copy()
        for e, h in zip(W, written):
            new_signs[e] = h[0]
        if channel in (1, 2) or channel >= 15:
            M = [[h[1][ids[e]] for e in W] for h in written]
            det = M[0][0]*M[1][1]-M[0][1]*M[1][0]
            assert abs(det) == 1
            a, b = (self.ell[e] for e in W)
            lnew = ((a*M[1][1]-b*M[1][0])//det, (-a*M[0][1]+b*M[0][0])//det)
            ell = list(self.ell)
            for e, value in zip(W, lnew):
                ell[e] = value
            for e in edges:
                if e not in W:
                    ell[e] -= sum(l*h[1][ids[e]] for l, h in zip(lnew, written))
            new_b = self.b+sum(l*h[1][-1] for l, h in zip(lnew, written))
            self.ell, self.b = canonical(ell, new_b)
        else:
            charges = [1 if h[0] == -1 else 0 if zero(h[1]) else 2 for h in target]
            rows = [loop_row(new_signs, self.pe[p, 2-i], self.pr[p, 2-i])[1] for i in range(3)]
            if 0 in charges:
                normal = rows[charges.index(0)]
            else:
                assert charges == [1, 1, 1]
                pair = next(i for i in (0, 1) if zero(tuple(a-b for a, b in zip(target[i][1], target[i+1][1]))))
                normal = tuple(a-b for a, b in zip(rows[pair], rows[pair+1]))
            self.ell, self.b = canonical(normal)
        self.signs = new_signs
        if 3 <= channel < 15:
            direction = int(all(h[0] == -1 for h in w))
            if self.events:
                assert direction != self.events[-1][2], 'one constraint cannot branch into two identities'
            self.events.append((t, p, direction))
        if shadow is not None and agreed:
            agreed = np.all(self.signs == np.where(shadow < modulus, 1, -1)) and (
                sum(a*int(x % modulus) for a, x in zip(self.ell, shadow))-self.b) % modulus == 0
        return bool(agreed)


def experiment(side=4, alpha=1, trials=64, seed=913003001):
    times = (0., .01, .05, .2, .5, 1., 2.)
    observations, mismatches, examples, histories, proposals = [], [], [], [], []
    for trial in range(trials):
        process = ConstraintProcess(side, seed+trial)
        modulus = 1_000_000_007
        shadow = process.finite_start(modulus, seed+10_000+trial)
        mismatch = None
        rng = np.random.default_rng(seed+20_000+trial)
        rate = len(process.pe)*(15+2*alpha)
        t, tick, row, attempts = 0., 0, [], 0
        while tick < len(times):
            t += rng.exponential(1/rate)
            while tick < len(times) and times[tick] <= t:
                row.append(process.features())
                if shadow is not None:
                    q = [kind(read(shadow, process.fe[f], process.fr[f], modulus), modulus) for f in range(process.F)]
                    if q != process.charges():
                        mismatch, shadow = times[tick], None
                tick += 1
            if tick == len(times):
                break
            p = int(rng.integers(len(process.pe)))
            attempts += 1
            slot = int(rng.integers(15+2*alpha))
            channel = slot if slot < 15 else 15+(slot-15)//alpha
            if not process.step(p, channel, t, shadow, modulus) and shadow is not None:
                mismatch, shadow = t, None
        observations.append(row)
        mismatches.append(mismatch)
        histories.append(process.events)
        proposals.append(attempts)
        if trial < 2:
            examples.append({'initial_direction': process.initial_direction,
                             'initial_signs': process.initial_signs.tolist(), 'initial_covector': process.initial_ell,
                             'final_signs': process.signs.tolist(), 'final_covector': process.ell,
                             'final_offset': process.b, 'conserved_charge': process.Q, 'reactions': process.events})
    x = np.array(observations, dtype=float)
    return {'side': side, 'alpha': alpha, 'trials': trials, 'seed': seed, 'times': times,
            'columns': ['reaction_rate', 'identity_count', 'reactive_pair_count', 'constraint_support_edges',
                        'constraint_L1', 'maximum_coefficient', 'absolute_offset', 'nonzero_offset', 'reaction_inactive',
                        'reaction_count', 'support_radius_in_primal_graph_hops'],
            'means': x.mean(axis=0).tolist(), 'independent_trajectory_SE': (x.std(axis=0, ddof=1)/np.sqrt(trials)).tolist(),
            'observations_by_trajectory': observations, 'finite_shadow_modulus': modulus,
            'first_shadow_mismatches': mismatches, 'reaction_histories': histories,
            'proposal_counts': proposals, 'examples': examples}


if __name__ == '__main__':
    result = {'scope': 'Exact integer one-constraint Palm limit; not a finite-n replacement', 'runs': []}
    for side, alpha, trials in ((4, 0, 64), (4, 1, 64), (4, 8, 64), (8, 1, 32)):
        row = experiment(side=side, alpha=alpha, trials=trials)
        result['runs'].append(row)
        Path('data/fiber-constraint-dynamics.json').write_text(json.dumps(result, indent=2)+'\n')
        print(json.dumps({'side': side, 'alpha': alpha, 'final_means': row['means'][-1],
                          'shadow_mismatches': sum(t is not None for t in row['first_shadow_mismatches'])}), flush=True)
