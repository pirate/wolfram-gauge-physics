#!/usr/bin/env python3
"""Exact reaction-channel calculation, not a simulation or a fitted model.

Reconstruct the twelve minimal reaction maps from their free symmetry orbits,
then resolve the averaged channel at fixed ordered boundary product.
"""
import itertools
import json
from collections import Counter
from pathlib import Path

from fiber_transport import FiberTransport
from triangle_patch_observer import WordObserver


def calculation():
    bank = json.loads(Path('data/triangle-feedback.json').read_text())['bank']
    g = FiberTransport(3, 3, [(0, 1), (1, 2), (2, 0)]).geometry.group
    q, tables = bank['element_charges'], bank['tables']
    states = list(itertools.product(range(g.n), repeat=3))
    ids = {s: i for i, s in enumerate(states)}
    product = lambda s: g.mul[g.mul[s[0]][s[1]]][s[2]]
    reflection = lambda s: all(q[v] == 1 for v in s) and ((s[0] == s[1]) != (s[1] == s[2]))
    mixed = lambda s: sorted(q[v] for v in s) == [0, 1, 2]
    rs = [s for s in states if reflection(s)]
    es = [s for s in states if mixed(s)]

    def action(s, frame, reverse):
        values = tuple(g.conj[frame][v] for v in s)
        return tuple(g.inv[v] for v in values[::-1]) if reverse else values

    actions = list(itertools.product(range(g.n), (False, True)))
    source = rs[0]
    assert len({action(source, *a) for a in actions}) == len(rs) == 12
    generated = []
    for target in es:
        if product(source) != product(target):
            continue
        pairs = [(ids[action(source, *a)], ids[action(target, *a)]) for a in actions]
        assert len({y for _, y in pairs}) == 12
        table = list(range(len(states)))
        for x, y in pairs:
            assert product(states[x]) == product(states[y])
            table[x], table[y] = y, x
        assert all(table[table[i]] == i for i in range(len(states)))
        generated.append(table)
    assert len(generated) == 12 and sorted(generated) == sorted(tables)

    blocks = []
    for p in range(g.n):
        if q[p] != 1:
            continue
        r = [ids[s] for s in rs if product(s) == p]
        e = [ids[s] for s in es if product(s) == p]
        assert len(r) == 4 and len(e) == 12
        for x in r:
            assert Counter(t[x] for t in tables) == Counter(e)
        for x in e:
            assert Counter(t[x] for t in tables) == Counter(r) + Counter({x: 8})
        block = r + e
        counts = [[sum(t[x] == y for t in tables) for y in block] for x in block]
        # Integer eigenbasis for 12 P. No numerical eigenvalue tolerance.
        eigenvectors = [(12, [1]*16), (-4, [3]*4 + [-1]*12)]
        for offset, size, eigenvalue in ((0, 4, 0), (4, 12, 8)):
            for i in range(size-1):
                v = [0]*16
                v[offset+i], v[offset+size-1] = 1, -1
                eigenvectors.append((eigenvalue, v))
        for eigenvalue, v in eigenvectors:
            assert [sum(a*b for a, b in zip(row, v)) for row in counts] == [eigenvalue*a for a in v]
        # Conditioned E -> R -> E returns uniformly over ALL 12 E states.
        for x in e:
            two_jump = Counter(t2[t1[x]] for t1 in tables if t1[x] != x for t2 in tables)
            assert two_jump == Counter({y: 4 for y in e})
        blocks.append({'boundary_product': p, 'reflection_states': [states[i] for i in r],
                       'mixed_states': [states[i] for i in e],
                       'eigenvalues_of_12P': dict(Counter(str(a) for a, _ in eigenvectors))})

    observer = WordObserver(g, q)
    # Two active reflection orbits and six mixed orbits. Each rule is already
    # equivariant, so this quotient is an exact lumping, not an approximation.
    rr = sorted({observer.canonical(s) for s in rs})
    ee = sorted({observer.canonical(s) for s in es})
    assert len(rr) == 2 and len(ee) == 6
    for s in rr:
        assert Counter(observer.canonical(states[t[ids[s]]]) for t in tables) == Counter({y: 2 for y in ee})
    for s in ee:
        assert Counter(observer.canonical(states[t[ids[s]]]) for t in tables) == Counter({y: 2 for y in rr}) + Counter({s: 8})
    # Actual two-reaction channel with two unchanged common-root rotation
    # spectators. Conditional on restoring the mixed charge order, the
    # complete three-rotation trace retains exactly one of its three contrasts.
    rotations = [v for v in range(g.n) if q[v] == 2]
    spectator_cases = 0
    for s in es:
        code = ids[s]
        slot = next(i for i, v in enumerate(s) if q[v] == 2)
        targets = [states[t2[t1[code]]] for t1 in tables if t1[code] != code for t2 in tables]
        targets = [v for v in targets if [q[a] for a in v] == [q[a] for a in s]]
        assert len(targets) == 8
        for a, b in itertools.product(rotations, repeat=2):
            initial = observer.canonical((a, b, s[slot]))
            same = sum(observer.canonical((a, b, v[slot])) == initial for v in targets)
            assert same == 4
            assert 4*same-len(targets) == len(targets)  # mean kernel = 1; trace normalization = 3
            spectator_cases += 1
    return {'derived_minimal_maps': len(generated), 'active_raw_states': len(rs)+len(es),
            'fixed_raw_states': len(states)-len(rs)-len(es), 'blocks': blocks,
            'active_gauge_orbits': {'reflection': rr, 'mixed': ee},
            'changing_E_to_R': 'uniform over 4 states at fixed ordered product',
            'changing_R_to_E': 'uniform over 12 states at fixed ordered product',
            'independent_two_jump_exact_return_probability': '1/12',
            'same_charge_order_probability_after_two_jumps': '1/6',
            'exact_return_given_same_charge_order': '1/2',
            'two_unchanged_rotation_spectator_cases': spectator_cases,
            'normalized_rotation_trace_given_restored_charge_order': '1/3'}


if __name__ == '__main__':
    print(json.dumps(calculation(), indent=2))
