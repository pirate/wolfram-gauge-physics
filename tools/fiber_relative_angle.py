#!/usr/bin/env python3
"""When independent holonomies expose angular motion to charge reactions.

Enumerate fixed-boundary fans and use the unchanged local channel bank.
No reaction rate is fitted to angle, and no physical energy is introduced.
"""
import itertools
import json
from collections import Counter
from fractions import Fraction
from pathlib import Path

import numpy as np
from scipy.sparse.linalg import LinearOperator, cg

from fiber_angular_transport import AngularFiber
from fiber_charge_coupling import build, integrated_correlation, lift_response_witness


def states_in_sector(f, faces, total_charge, boundary):
    g = f.g
    pools = [[g.identity], f.reflections, f.rotations]
    charge = {h: q for q, pool in enumerate(pools) for h in pool}
    states = []
    for layout in itertools.product(range(3), repeat=faces):
        if sum(layout) != total_charge or 1 not in layout:
            continue
        # Solve one reflection from L*h*R=P; this avoids scanning a redundant label.
        pivot = layout.index(1)
        positions = [i for i in range(faces) if i != pivot]
        for entries in itertools.product(*(pools[layout[i]] for i in positions)):
            w = [g.identity]*faces
            for i, h in zip(positions, entries):
                w[i] = h
            left, right = g.identity, g.identity
            for h in w[:pivot]:
                left = g.mul[left][h]
            for h in w[pivot+1:]:
                right = g.mul[right][h]
            w[pivot] = g.mul[g.mul[g.inv[left]][boundary]][g.inv[right]]
            if charge[w[pivot]] == 1:
                states.append(tuple(w))
    assert len(states) == len(set(states))
    return states


def algebraic_witness(f, flat_boundary=False):
    g, n = f.g, f.n
    r = lambda k: g.elements.index(tuple((k-v) % n for v in range(n)))
    t = lambda k: g.elements.index(tuple((v+k) % n for v in range(n)))
    before = (r(1), r(0), r(0), t(1))
    after = (r(1), r(0), r(1), t(2))
    if flat_boundary:
        before, after = before+(r(0),), after+(r(0),)
    assert f.rules[0].apply_pair(before[2:4]) == after[2:4]
    witness = {'before_tuple': before, 'after_tuple': after, 'pair_position': 2,
               'angular_channel_index': 0,
               'before_labels': ['r1', 'r0', 'r0', 't1']+(['r0'] if flat_boundary else []),
               'after_labels': ['r1', 'r0', 'r1', 't2']+(['r0'] if flat_boundary else [])}
    witness['shared_link_realization'] = lift_response_witness(n, witness)
    return witness


def charge_drift(q):
    """Exact full-bank drift: it depends on charges, not on holonomy labels."""
    b = [0]*len(q)
    for i in range(len(q)-1):
        current = (2+(0 in q[i:i+2]))*(q[i+1]-q[i])
        b[i] += current
        b[i+1] -= current
    for i in range(len(q)-2):
        if sorted(q[i:i+3]) == [0, 1, 2]:
            for j in range(i, i+3):
                b[j] += 4*(1-q[j])
    return b


def first_reaction(f, states, charges, A, B, flat_boundary):
    # Killing, not disabling channels: every original channel runs until
    # the first transition out of the one-rotation inventory.
    selected = np.array([i for i, q in enumerate(charges) if q.count(2) == 1])
    ids = {states[i]: j for j, i in enumerate(selected)}
    n, g = f.n, f.g
    r = lambda k: g.elements.index(tuple((k-v) % n for v in range(n)))
    t = lambda k: g.elements.index(tuple((v+k) % n for v in range(n)))
    seeds = [(r(1), r(0), r(k-1), t(k))+((r(0),) if flat_boundary else ()) for k in range(1, n)]
    records = {}
    for name, L in (('angular_off', A), ('angular_on', A+B)):
        K = L[selected][:, selected]
        rhs = np.ones(len(selected))
        M = LinearOperator(K.shape, matvec=lambda x: x/K.diagonal())
        u, info = cg(K, rhs, M=M, rtol=1e-12, atol=0, maxiter=20000)
        assert info == 0
        residual = np.linalg.norm(K@u-rhs)/np.linalg.norm(rhs)
        records[name] = {'uniform_mean_time': float(u.mean()),
                         'initial_uniform_hazard': float((K@rhs).mean()),
                         'seed_mean_times_by_relative_exponent': [float(u[ids[w]]) for w in seeds],
                         'poisson_relative_residual': float(residual)}
    return {'reactant_states': len(selected), 'relative_exponents': list(range(1, n)), **records}


def activity_closure_witness(f, states, charges, arows, bedges, flat_boundary):
    if f.n < 5:
        return None
    g, n = f.g, f.n
    r = lambda k: g.elements.index(tuple((k-v) % n for v in range(n)))
    t = lambda k: g.elements.index(tuple((v+k) % n for v in range(n)))
    def key(i):
        w = states[i]
        activity = tuple(int(all(h in f.reflections for h in w[j:j+3]) and
                             ((w[j] == w[j+1]) != (w[j+1] == w[j+2]))) for j in range(len(w)-2))
        return tuple(charges[i]), activity
    ids = {w: i for i, w in enumerate(states)}
    seeds = [(r(1), r(0), r(k-1), t(k))+((r(0),) if flat_boundary else ()) for k in (2, 3)]
    source_ids = [ids[w] for w in seeds]
    source = key(source_ids[0])
    assert key(source_ids[1]) == source
    target = (source[0], (1,)+(0,)*(len(seeds[0])-3))
    rates = []
    for i in source_ids:
        legacy = sum(rate for j, rate in arows[i].items() if key(j) == target)
        angular = sum(1 for a, b, _, _ in bedges if a == i and key(b) == target)
        rates.append({'legacy_rate': legacy, 'angular_rate': angular})
    assert rates == [{'legacy_rate': 2, 'angular_rate': 1}, {'legacy_rate': 0, 'angular_rate': 0}]
    return {'relative_exponents': [2, 3], 'source_tuples': seeds, 'common_source_charge_and_activity': source,
            'target_charge_and_activity': target, 'transition_rates': rates}


def calculation(n, faces, total_charge, flat_boundary):
    f = AngularFiber(n)
    g = f.g
    boundary = g.identity if flat_boundary else f.reflections[0]
    states = states_in_sector(f, faces, total_charge, boundary)
    f, states, arows, bedges, A, B = build(n, faces, states)
    N = len(states)
    charge = lambda h: 0 if h == g.identity else 1 if h in f.reflections else 2
    charges = [[charge(h) for h in w] for w in states]
    records = []
    for name, initial in (
        ('endpoint_charge', [faces*q[0]-total_charge for q in charges]),
        ('charge_dipole', [sum((2*i-faces+1)*v for i, v in enumerate(q)) for q in charges]),
        ('rotation_count', [q.count(2) for q in charges])):
        values, energies, witness = initial[:], [], None
        for depth in range(9):
            numerator = sum((values[i]-values[j])**2 for i, j, _, _ in bedges)
            energies.append(str(Fraction(numerator, 2*N)))
            if numerator:
                i, j, pair, channel = next(e for e in bedges if values[e[0]] != values[e[1]])
                witness = {'before_tuple': states[i], 'after_tuple': states[j],
                           'pair_position': pair, 'angular_channel_index': channel,
                           'charge_field': charges[i], 'response_before': values[i], 'response_after': values[j]}
                assert charges[i] == charges[j]
                break
            values = [sum(w*(values[i]-values[j]) for j, w in row.items()) for i, row in enumerate(arows)]
        observable = np.array(initial, dtype=float)
        observable -= observable.mean()
        off, on = integrated_correlation(A, observable), integrated_correlation(A+B, observable)
        records.append({'observable': name, 'first_visible_depth': depth if witness else None,
                        'first_changed_autocorrelation_derivative': 2*depth+1 if witness else None,
                        'exact_angular_dirichlet_energies': energies, 'witness': witness,
                        'angular_off': off, 'angular_on': on,
                        'relative_integrated_correlation_reduction': 1-on['normalized_integrated_autocorrelation']/off['normalized_integrated_autocorrelation']})
    witness = algebraic_witness(f, flat_boundary)
    ids = {w: i for i, w in enumerate(states)}
    rotations = [q.count(2) for q in charges]
    witness['positive_generator_on_rotation_count'] = [
        sum(rate*(rotations[i]-rotations[j]) for j, rate in arows[i].items())
        for i in (ids[tuple(witness['before_tuple'])], ids[tuple(witness['after_tuple'])])]
    assert witness['positive_generator_on_rotation_count'] == [-12, 0]
    moments = []
    for i in (ids[tuple(witness['before_tuple'])], ids[tuple(witness['after_tuple'])]):
        drift = [sum(rate*(charges[j][a]-charges[i][a]) for j, rate in arows[i].items()) for a in range(faces)]
        covariance = [[sum(rate*(charges[j][a]-charges[i][a])*(charges[j][b]-charges[i][b])
                           for j, rate in arows[i].items()) for b in range(faces)] for a in range(faces)]
        assert drift == charge_drift(charges[i])
        second = [sum(rate*(charge_drift(charges[j])[a]-drift[a]) for j, rate in arows[i].items()) for a in range(faces)]
        moments.append({'charge_drift': drift, 'instantaneous_charge_covariance': covariance,
                        'second_generator_power_on_charge': second})
    witness['conditional_charge_moments'] = moments
    return {'cycle_vertices': n, 'faces': faces, 'total_charge': total_charge,
            'boundary': 'identity' if flat_boundary else 'reflection', 'framed_states': N,
            'inventory_counts': dict(Counter(''.join('ERZ'[v] for v in q) for q in charges)),
            'observables': records, 'algebraic_witness': witness,
            'first_reaction': first_reaction(f, states, charges, A, B, flat_boundary),
            'charge_plus_activity_nonclosure': activity_closure_witness(f, states, charges, arows, bedges, flat_boundary),
            'scope': 'Invariant sector containing at least one reflection. Each unchanged channel has its own rate-one clock.'}


if __name__ == '__main__':
    output = Path('data/fiber-relative-angle.json')
    report = {'scope': 'Relative-holonomy coupling to charge, without a physical energy or continuum force.', 'calculations': []}
    for args in ((3, 4, 5, False), (5, 4, 5, False), (7, 4, 5, False), (11, 4, 5, False),
                 (3, 5, 6, True), (5, 5, 6, True), (7, 5, 6, True)):
        row = calculation(*args)
        report['calculations'].append(row)
        output.write_text(json.dumps(report, indent=2)+'\n')
        print(json.dumps({k: row[k] for k in ('cycle_vertices', 'faces', 'total_charge', 'framed_states')} |
                         {'observables': [{k: o[k] for k in ('observable', 'first_visible_depth', 'relative_integrated_correlation_reduction')} for o in row['observables']]}), flush=True)
