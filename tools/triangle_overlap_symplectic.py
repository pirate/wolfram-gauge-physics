#!/usr/bin/env python3
"""Derive the finite symplectic action inside the five-reflection sector.

The alternating form is solved from primitive matrices, not supplied as an
equation of motion. Its coordinates are internal fiber labels, not spacetime.
"""
import itertools
import json
from collections import deque
from functools import reduce

from triangle_overlap_algebra import build


def multiply(a, b):
    return tuple(sum(a[4*i+k]*b[4*k+j] for k in range(4)) % 3 for i in range(4) for j in range(4))


def apply(a, v):
    return tuple(sum(a[4*i+j]*v[j] for j in range(4)) % 3 for i in range(4))


def transpose(a):
    return tuple(a[4*j+i] for i in range(4) for j in range(4))


def projective(v):
    return min(tuple(v), tuple(-a % 3 for a in v))


def calculation():
    x, u, channels, tables, faces, path = build()
    g = u.g
    adjacency = {f: [h for h in faces if h != f and
                    len(set(x.e.geometry.faces[f][:-1]) & set(x.e.geometry.faces[h][:-1])) == 2]
                 for f in faces}
    start = min(f for f in faces if len(adjacency[f]) == 1)
    order = [start]
    while len(order) < len(faces):
        order.append(next(h for h in adjacency[order[-1]] if h not in order))
    assert all(x.e.geometry.faces[f][0] == u.root for f in faces)
    product = lambda values: reduce(lambda a, b: g.mul[a][b], values, g.identity)
    raw_faces = [x.e.geometry.holonomies(u.lift(v)) for v in u.reps]
    outer = [x.e.factor.oracle.fan.transport(u.lift(v), path) for v in u.reps]
    possibilities = [o for o in (order, order[::-1])
                     if all(product([row[f] for f in o]) == p for row, p in zip(raw_faces, outer))]
    assert len(possibilities) == 1
    order = possibilities[0]
    reflections = [i for i, q in enumerate(x.e.charges) if q == 1]
    labels = {r: g.elements[r][0] for r in reflections}
    assert all(g.elements[r] == tuple((labels[r]-y) % 3 for y in range(3)) for r in reflections)
    positions = [i for i, row in enumerate(raw_faces) if all(row[f] in reflections for f in order)]
    coordinates = {}
    for i in positions:
        a = [labels[raw_faces[i][f]] for f in order]
        p = labels[outer[i]]
        assert sum((-1)**j*a[j] for j in range(5)) % 3 == p
        coordinates[i] = projective([(v-p) % 3 for v in a[:4]])
    assert len(set(coordinates.values())) == len(positions) == 41

    # Restrict the exact Hurwitz move (a,b)->(2a-b,a) to the hyperplane
    # v0-v1+v2-v3+v4=0. First four entries are independent coordinates.
    matrices = []
    for j in range(4):
        columns = []
        for k in range(4):
            v = [int(i == k) for i in range(4)]
            v.append((-v[0]+v[1]-v[2]+v[3]) % 3)
            v[j:j+2] = [(2*v[j]-v[j+1]) % 3, v[j]]
            assert sum((-1)**i*v[i] for i in range(5)) % 3 == 0
            columns.append(v[:4])
        matrices.append(tuple(columns[k][i] for i in range(4) for k in range(4)))
    identity = tuple(int(i == j) for i in range(4) for j in range(4))
    inverse_matrices = [multiply(m, m) for m in matrices]
    assert all(multiply(a, b) == identity for a, b in zip(matrices, inverse_matrices))
    predicted = {tuple(projective(apply(m, coordinates[i])) for i in positions)
                 for m in [identity]+matrices+inverse_matrices}
    actual = {tuple(coordinates[t[i]] for i in positions)
              for (p, r), t in zip(channels, tables) if r < 3}
    assert actual == predicted and len(actual) == 9

    pairs = list(itertools.combinations(range(4), 2))
    forms = []
    for coefficients in itertools.product(range(3), repeat=6):
        form = [0]*16
        for (i, j), a in zip(pairs, coefficients):
            form[4*i+j], form[4*j+i] = a, -a % 3
        form = tuple(form)
        if all(multiply(transpose(m), multiply(form, m)) == form for m in matrices):
            forms.append(form)
    assert len(forms) == 3
    form = min(f for f in forms if any(f))
    # Nondegeneracy is an exact kernel calculation over the 81 vectors.
    vectors = list(itertools.product(range(3), repeat=4))
    assert [v for v in vectors if not any(apply(form, v))] == [(0, 0, 0, 0)]
    generated, queue = {identity}, [identity]
    for m in queue:
        for a in matrices:
            target = multiply(a, m)
            if target not in generated:
                generated.add(target)
                queue.append(target)
    expected = 3**4*(3**2-1)*(3**4-1)
    assert len(generated) == expected == 51840
    assert tuple(-a % 3 for a in identity) in generated
    quotient = {min(m, tuple(-a % 3 for a in m)) for m in generated}
    assert len(quotient) == 25920
    orbit, queue = {(1, 0, 0, 0)}, deque([(1, 0, 0, 0)])
    while queue:
        v = queue.popleft()
        for m in matrices:
            w = apply(m, v)
            if w not in orbit:
                orbit.add(w)
                queue.append(w)
    assert len(orbit) == 80 and len({projective(v) for v in orbit}) == 40
    return {'ordered_face_ids': order, 'reflection_affine_labels': labels,
            'reflection_gauge_states': len(positions), 'generators': [list(m) for m in matrices],
            'invariant_alternating_forms_dimension': 1, 'alternating_form': list(form),
            'matrix_group_order': len(generated), 'projective_group_order': len(quotient),
            'nonzero_vector_orbit': len(orbit), 'nonzero_projective_orbit': 40,
            'actual_shared_edge_primitive_actions_equal_projective_matrix_actions': True,
            'scope': 'Elastic all-reflection sector only; reactions leave this chart. No physical symplectic phase space or quantum amplitudes asserted.'}


if __name__ == '__main__':
    print(json.dumps(calculation(), indent=2))
