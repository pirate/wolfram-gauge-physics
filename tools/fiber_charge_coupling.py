#!/usr/bin/env python3
"""Detect when adjacency-assisted fiber dynamics becomes visible to charge.

Exact integer response depths and common-measure autocorrelation integrals;
no fitted charge equation or altered physical Hamiltonian is assumed.
"""
import itertools
import json
from collections import Counter
from fractions import Fraction
from pathlib import Path

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.sparse.linalg import LinearOperator, cg

from fiber_angular_transport import AngularFiber
from cycle_relational_dynamics import TransportRule, family


def single_fan_charge_closure():
    """Exact six-layout restriction; the constant layout vector is excluded."""
    layouts = list(itertools.permutations((0, 1, 2)))
    K = [[Fraction(4 if i == j else 0) for j in range(6)] for i in range(6)]
    for i, layout in enumerate(layouts):
        for position in (0, 1):
            target = list(layout)
            target[position], target[position+1] = target[position+1], target[position]
            weight = 3 if 0 in layout[position:position+2] else 2
            K[i][i] += weight
            K[i][layouts.index(tuple(target))] -= weight
    results = {}
    for name, v in (
        ('endpoint_charge', [3*w[0]-3 for w in layouts]),
        ('charge_dipole', [2*(w[2]-w[0]) for w in layouts])):
        augmented = [row[:]+[Fraction(value)] for row, value in zip(K, v)]
        for i in range(6):
            pivot = augmented[i][i]
            augmented[i] = [x/pivot for x in augmented[i]]
            for j in range(6):
                if i != j:
                    factor = augmented[j][i]
                    augmented[j] = [x-factor*y for x, y in zip(augmented[j], augmented[i])]
        results[name] = str(sum(x*row[-1] for x, row in zip(v, augmented))/sum(x*x for x in v))
    return {'layouts_as_charges': layouts, 'restricted_positive_generator': [[int(x) for x in row] for row in K],
            'exact_normalized_integrated_autocorrelations': results,
            'scope': 'Zero on RRR, zero-sum functions of the six ERZ layouts. Valid for every odd cycle n >= 3.'}


def lift_response_witness(n, witness):
    """Realize the response-distinguishable pair on shared links, exterior fixed."""
    from run_three_face import LinkOracle
    f = AngularFiber(n)
    g, oracle = f.g, LinkOracle(4, f.g)
    before, after = tuple(witness['before_tuple']), tuple(witness['after_tuple'])
    faces, boundary = len(before), g.identity
    for h in before:
        boundary = g.mul[boundary][h]
    assert 3 <= faces <= 5
    root, successor = 0, {}
    for face in oracle.faces:
        if root in face[:3]:
            i = face[:3].index(root)
            successor[face[(i+1) % 3]] = face[(i+2) % 3]
    rim = [min(successor)]
    for _ in range(faces):
        rim.append(successor[rim[-1]])
    ids = {edge: i for i, edge in enumerate(oracle.edges)}
    def path(vertices):
        return [(ids[tuple(sorted((u, v)))], u > v) for u, v in zip(vertices, vertices[1:])]
    native = [path((root, rim[i], rim[i+1], root)) for i in range(faces)]
    ordered = native[::-1]
    raw, spoke = [g.identity]*len(oracle.edges), g.identity
    for i, h in enumerate(before[::-1]):
        rim_value = boundary if i == faces-1 else g.identity
        edge, reverse = path((rim[i], rim[i+1]))[0]
        raw[edge] = g.inv[rim_value] if reverse else rim_value
        spoke = g.mul[g.mul[rim_value][spoke]][g.inv[h]]
        edge, reverse = path((root, rim[i+1]))[0]
        raw[edge] = g.inv[spoke] if reverse else spoke
    assert spoke == g.identity
    assert tuple(oracle.transport(raw, p) for p in ordered) == before
    start = min(witness['pair_position'], faces-3)
    patch = next(i for i, p in enumerate(oracle.patches) if p[::-1] == ordered[start:start+3])
    code = lambda w: (w[0]*g.n+w[1])*g.n+w[2]
    left, right = code(before[start:start+3]), code(after[start:start+3])
    evolved = raw[:]
    oracle.update(evolved, patch, {left: right})
    assert tuple(oracle.transport(evolved, p) for p in ordered) == after
    changed = [i for i, (a, b) in enumerate(zip(raw, evolved)) if a != b]
    assert len(changed) == 1 and set(changed) <= oracle.writes[patch]
    charge = lambda h: 0 if h == g.identity else 1 if h in f.reflections else 2
    initial_charges = [charge(oracle.transport(raw, p)) for p in oracle.face_paths]
    final_charges = [charge(oracle.transport(evolved, p)) for p in oracle.face_paths]
    assert initial_charges == final_charges
    assert sum(initial_charges) == sum(charge(h) for h in before)+charge(boundary)
    restored = evolved[:]
    oracle.update(restored, patch, {right: left})
    assert restored == raw
    return {'base_mesh_side': 4, 'root': root, 'rim_vertices': rim, 'native_three_face_patch': patch,
            'edges': oracle.edges, 'initial_links': raw, 'final_links': evolved,
            'changed_edge_indices': changed, 'whole_mesh_charge': sum(initial_charges),
            'exterior_fixed': True, 'inverse_restores_links': True}


def states_at_charge_three(f, faces):
    g, p = f.g, f.reflections[0]
    states = []
    for positions in itertools.combinations(range(faces), 3):
        for a in f.reflections:
            for b in f.reflections:
                c = g.mul[g.inv[g.mul[a][b]]][p]
                if a == b == c == p:
                    continue
                w = [g.identity]*faces
                for position, h in zip(positions, (a, b, c)):
                    w[position] = h
                states.append(tuple(w))
    for i, j in itertools.permutations(range(faces), 2):
        for z in f.rotations:
            r = g.mul[p][g.inv[z]] if i < j else g.mul[g.inv[z]][p]
            w = [g.identity]*faces
            w[i], w[j] = r, z
            states.append(tuple(w))
    assert len(set(states)) == len(states)
    return states


def build(n, faces, states=None):
    f = AngularFiber(n)
    g = f.g
    states = states_at_charge_three(f, faces) if states is None else states
    ids, N = {w: i for i, w in enumerate(states)}, len(states)
    transports, reactions = [TransportRule(g, k) for k in range(3)], family(g)
    arows, brows, angular_edges = [], [], []
    for i, w in enumerate(states):
        a, b = Counter(), Counter()
        for position in range(faces-1):
            pair = w[position:position+2]
            for rule in transports:
                target_pair = rule.apply((*pair, g.identity))[:2]
                target = w[:position]+target_pair+w[position+2:]
                j = ids[target]
                if i != j:
                    a[j] += 1
            z = pair[0] if pair[0] in f.distance and pair[1] in f.reflections else pair[1] if pair[1] in f.distance and pair[0] in f.reflections else None
            if z is not None:
                d = f.distance[z]
                channels = {d-1} | ({d-2} if d > 1 else set())
                for channel in channels:
                    target_pair = f.rules[channel].apply_pair(pair)
                    target = w[:position]+target_pair+w[position+2:]
                    j = ids[target]
                    assert i != j
                    b[j] += 1
                    angular_edges.append((i, j, position, channel))
        for position in range(faces-2):
            triple = w[position:position+3]
            for rule in reactions:
                target = w[:position]+rule.apply(triple)+w[position+3:]
                j = ids[target]
                if i != j:
                    a[j] += 1
        arows.append(a); brows.append(b)
    def matrix(rates):
        rows, columns, values = [], [], []
        for i, row in enumerate(rates):
            rows.append(i); columns.append(i); values.append(sum(row.values()))
            for j, value in row.items():
                rows.append(i); columns.append(j); values.append(-value)
        L = coo_matrix((values, (rows, columns)), shape=(N, N), dtype=float).tocsr()
        assert (L-L.T).nnz == 0
        return L
    return f, states, arows, angular_edges, matrix(arows), matrix(brows)


def integrated_correlation(L, observable):
    count, labels = connected_components(L, directed=False)
    sizes = np.bincount(labels)
    projection = lambda v: (np.bincount(labels, weights=v)/sizes)[labels]
    assert np.max(np.abs(projection(observable))) < 1e-12
    operator = LinearOperator(L.shape, matvec=lambda v: L@v+projection(v))
    preconditioner = LinearOperator(L.shape, matvec=lambda v: v/(L.diagonal()+1/sizes[labels]))
    solution, info = cg(operator, observable, M=preconditioner, rtol=1e-12, atol=0, maxiter=20000)
    assert info == 0
    residual = np.linalg.norm(operator@solution-observable)/np.linalg.norm(observable)
    return {'normalized_integrated_autocorrelation': float(observable@solution/(observable@observable)),
            'poisson_relative_residual': float(residual), 'stationary_components': int(count)}


def calculation(n, faces):
    f, states, arows, bedges, A, B = build(n, faces)
    charge = lambda h: 0 if h == f.g.identity else 1 if h in f.reflections else 2
    N = len(states)
    records = []
    for name, initial in (
        ('endpoint_charge', [faces*charge(w[0])-3 for w in states]),
        ('charge_dipole', [sum((2*i-faces+1)*charge(h) for i, h in enumerate(w)) for w in states])):
        value, energies, witness = initial[:], [], None
        for depth in range(17):
            numerator = sum((value[i]-value[j])**2 for i, j, _, _ in bedges)
            energies.append(str(Fraction(numerator, 2*N)))
            if numerator:
                i, j, position, channel = next(e for e in bedges if value[e[0]] != value[e[1]])
                witness = {'before_tuple': states[i], 'after_tuple': states[j], 'pair_position': position,
                           'angular_channel_index': channel, 'response_before': value[i], 'response_after': value[j],
                           'same_charge_field': [charge(h) for h in states[i]]}
                assert [charge(h) for h in states[i]] == [charge(h) for h in states[j]]
                break
            value = [sum(weight*(value[i]-value[j]) for j, weight in row.items()) for i, row in enumerate(arows)]
        v = np.array(initial, dtype=float)
        off, on = integrated_correlation(A, v), integrated_correlation(A+B, v)
        assert on['normalized_integrated_autocorrelation'] <= off['normalized_integrated_autocorrelation']+1e-10
        if witness and n == 9 and faces == 5:
            witness['shared_link_realization'] = lift_response_witness(n, witness)
        records.append({'observable': name, 'first_nonzero_angular_response_depth': depth if witness else None,
                        'angular_dirichlet_energy_of_A_power_observable': energies,
                        'first_changed_autocorrelation_derivative_order': 2*depth+1 if witness else None,
                        'angular_off': off, 'angular_on': on, 'response_witness': witness})
    return {'cycle_vertices': n, 'faces_in_ordered_fan': faces, 'framed_states': N,
            'response_depth_search_maximum': 16,
            'legacy_channels': 3*(faces-1)+12*(faces-2), 'angular_channels': f.m*(faces-1),
            'observables': records,
            'scope': 'Fixed boundary reflection and charge three; angular clocks added without renormalizing legacy clocks.'}


if __name__ == '__main__':
    output = Path('data/fiber-charge-coupling.json')
    report = {'scope': 'Charge-response coupling to angular primitives, not a derived spacetime transport coefficient.',
              'single_fan_exact_closure': single_fan_charge_closure(), 'calculations': []}
    for n, faces in ((9, 3), (9, 5), (27, 5), (81, 5), (9, 4), (5, 5), (7, 5), (25, 4), (49, 4)):
        row = calculation(n, faces)
        report['calculations'].append(row)
        output.write_text(json.dumps(report, indent=2)+'\n')
        print(json.dumps({'cycle': n, 'faces': faces, 'states': row['framed_states'],
                          'observables': [{k: o[k] for k in ('observable', 'first_nonzero_angular_response_depth',
                                                           'first_changed_autocorrelation_derivative_order', 'angular_off', 'angular_on')}
                                          for o in row['observables']]}), flush=True)
