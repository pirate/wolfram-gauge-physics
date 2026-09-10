#!/usr/bin/env python3
"""Couple actual cycle-fiber links to the exclusion process they approach.

Affine coordinates are the exact derived cycle automorphisms, not a new
force law. Angular proposals use exact generator thinning of idle slots.
"""
import itertools
import json
import math
from pathlib import Path

import numpy as np
from numba import njit

from fiber_angular_transport import AngularFiber
from cycle_relational_dynamics import TransportRule, family
from run_three_face import LinkOracle

LAYOUTS = np.array(list(itertools.permutations((0, 1, 2))), dtype=np.int64)


@njit
def mul(a, b, n):
    sa, sb = (1 if a < n else -1), (1 if b < n else -1)
    return (a % n+sa*(b % n)) % n+(n if sa*sb == -1 else 0)


@njit
def inv(a, n):
    return (-a) % n if a < n else a


@njit
def kind(a, n):
    return 0 if a == 0 else 1 if a >= n else 2


@njit
def forward(r, s, layout, orientation, n):
    z = mul(r, s, n) if orientation else mul(s, r, n)
    zi, ri = 0, 0
    for i in range(3):
        if layout[i] == 2:
            zi = i
        elif layout[i] == 1:
            ri = i
    h = mul(inv(z, n), s, n) if zi < ri else mul(s, inv(z, n), n)
    out = np.zeros(3, dtype=np.int64)
    out[zi], out[ri] = z, h
    return out[0], out[1], out[2]


@njit
def primitive(a, b, c, channel, n):
    if channel == 0:
        return (b, a, c) if a == 0 or b == 0 else (a, b, c)
    if channel == 1:
        return mul(mul(a, b, n), inv(a, n), n), a, c
    if channel == 2:
        return b, mul(mul(inv(b, n), a, n), b, n), c
    if channel >= 15:
        z = a if kind(a, n) == 2 and kind(b, n) == 1 else b if kind(a, n) == 1 and kind(b, n) == 2 else 0
        if not z:
            return a, b, c
        target = (z+(1 if channel == 15 else -1)) % n
        if target == 0:
            return a, b, c
        product = mul(a, b, n)
        return (mul(product, inv(target, n), n), target, c) if a >= n else (target, mul(inv(target, n), product, n), c)
    index = channel-3
    layout, orientation = LAYOUTS[index//2], index % 2 == 0
    ka, kb, kc = kind(a, n), kind(b, n), kind(c, n)
    if ka == kb == kc == 1:
        if a == b and b != c:
            return forward(a, c, layout, orientation, n)
        if b == c and a != b:
            x, y, z = forward(c, a, layout, orientation, n)
            return inv(z, n), inv(y, n), inv(x, n)
    elif ka == layout[0] and kb == layout[1] and kc == layout[2]:
        s = mul(mul(a, b, n), c, n)
        z = a if ka == 2 else b if kb == 2 else c
        r = mul(z, s, n) if orientation else mul(s, z, n)
        return r, r, s
    elif ka == layout[2] and kb == layout[1] and kc == layout[0]:
        s = inv(mul(mul(a, b, n), c, n), n)
        z = inv(a if ka == 2 else b if kb == 2 else c, n)
        r = mul(z, s, n) if orientation else mul(s, z, n)
        return inv(s, n), inv(r, n), inv(r, n)
    return a, b, c


@njit
def read(raw, edges, reverse, n):
    h = 0
    for i in range(3):
        a = raw[edges[i]]
        h = mul(inv(a, n) if reverse[i] else a, h, n)
    return h


@njit
def evolve(raw, n, alpha, duration, seed, pe, pr, fe, fr, supports):
    np.random.seed(seed)
    q = np.array([kind(read(raw, fe[i], fr[i], n), n) for i in range(len(fe))])
    parity = np.array([int(read(raw, fe[i], fr[i], n) < n) for i in range(len(fe))])
    initial_identity = int(np.sum(q == 0))
    matched = initial_identity == 0
    reaction_events, vacancy_events = 0, 0
    channels = 15+2*alpha
    attempts = np.random.poisson(len(pe)*channels*duration)
    for _ in range(attempts):
        slot = np.random.randint(len(pe)*channels)
        p, channel = slot//channels, slot % channels
        if channel >= 15:
            channel = 15+(channel-15)//alpha
        if channel == 1 or channel == 2:
            i, j = supports[p, 0], supports[p, 1]
            parity[i], parity[j] = parity[j], parity[i]
        a = read(raw, pe[p, 2], pr[p, 2], n)
        b = read(raw, pe[p, 1], pr[p, 1], n)
        c = read(raw, pe[p, 0], pr[p, 0], n)
        aa, bb, cc = primitive(a, b, c, channel, n)
        if (a, b, c) != (aa, bb, cc):
            prefix = 0
            for k in range(2):
                v = raw[pe[p, 0, k]]
                prefix = mul(inv(v, n) if pr[p, 0, k] else v, prefix, n)
            v = raw[pe[p, 1, 1]]
            if pr[p, 1, 1]:
                v = inv(v, n)
            s2 = mul(prefix, inv(cc, n), n)
            s3 = mul(mul(v, prefix, n), inv(mul(bb, cc, n), n), n)
            raw[pe[p, 0, 2]] = s2 if pr[p, 0, 2] else inv(s2, n)
            raw[pe[p, 1, 2]] = s3 if pr[p, 1, 2] else inv(s3, n)
            q[supports[p, 0]], q[supports[p, 1]], q[supports[p, 2]] = kind(aa, n), kind(bb, n), kind(cc, n)
            reaction_events += int(3 <= channel < 15)
            vacancy_events += int(channel == 0)
        for face in supports[p]:
            if q[face] != 1+parity[face]:
                matched = False
    actual = np.array([kind(read(raw, fe[i], fr[i], n), n) for i in range(len(fe))])
    assert np.all(actual == q)
    assert matched or initial_identity or reaction_events
    return q, parity, matched, initial_identity, reaction_events, vacancy_events, raw


def geometry(side):
    oracle = LinkOracle(side, AngularFiber(3).g)
    p, f = np.array(oracle.patches), np.array(oracle.face_paths)
    ids = {frozenset(face[:3]): i for i, face in enumerate(oracle.faces)}
    supports = np.array([[ids[frozenset(face[:3])] for face in spec[::-1]] for spec in oracle.specs])
    return oracle, (p[..., 0], p[..., 1], f[..., 0], f[..., 1], supports)


def representation_agreement(n):
    f = AngularFiber(n)
    labels = [p[0]+(0 if (p[1]-p[0]) % n == 1 else n) for p in f.g.elements]
    rules = [TransportRule(f.g, i) for i in range(3)]+family(f.g)
    compared = 0
    for w in itertools.product(range(f.g.n), repeat=3):
        encoded = tuple(labels[h] for h in w)
        for channel, rule in enumerate(rules):
            assert primitive(*encoded, channel, n) == tuple(labels[h] for h in rule.apply(w))
            compared += 1
        # The two directional proposals reproduce the sum of all angular
        # channels, not their individual labels or their idle-attempt counts.
        old = sorted(tuple(labels[h] for h in rule.apply(w)) for rule in f.rules if rule.apply(w) != w)
        new = sorted(primitive(*encoded, channel, n) for channel in (15, 16)
                     if primitive(*encoded, channel, n) != encoded)
        assert old == new
    return {'cycle_vertices': n, 'old_channel_images_compared': compared,
            'angular_generator_tuple_count': (2*n)**3}


def stationary_coupling(n, alpha, trials=4096, side=4, duration=0.02, seed=831004001):
    oracle, args = geometry(side)
    rng = np.random.default_rng(seed)
    rows = []
    for i in range(trials):
        raw = rng.integers(2*n, size=len(oracle.edges), dtype=np.int64)
        q, parity, matched, identities, reactions, vacancies, _ = evolve(raw, n, alpha, duration, seed+i, *args)
        rows.append([int(not matched), identities, reactions, vacancies, np.abs(q-1-parity).mean()])
    x = np.array(rows)
    faces = len(oracle.faces)
    return {'n': n, 'side': side, 'faces': faces, 'angular_rate': alpha, 'trials': trials,
            'duration': duration, 'seed': seed,
            'columns': ['path_mismatch', 'initial_identities', 'reaction_events', 'vacancy_events', 'final_L1_error_per_face'],
            'means': x.mean(axis=0).tolist(), 'independent_trajectory_SE': (x.std(axis=0, ddof=1)/np.sqrt(trials)).tolist(),
            'exact_expected_reactions': 18*faces*duration*(n-1)/n**2,
            'path_mismatch_upper_bound': min(1., faces/(2*n)+18*faces*duration*(n-1)/n**2),
            'expected_L1_error_upper_bound': 1/n+(39/n-36/n**2)*duration}


if __name__ == '__main__':
    result = {'representation': [representation_agreement(n) for n in (3, 5)], 'stationary_couplings': []}
    path = Path('data/fiber-refinement-limit.json')
    for n in (9, 27, 81, 243, 729):
        for alpha in (0, 1, 8):
            row = stationary_coupling(n, alpha)
            result['stationary_couplings'].append(row)
            path.write_text(json.dumps(result, indent=2)+'\n')
            print(json.dumps(row), flush=True)
