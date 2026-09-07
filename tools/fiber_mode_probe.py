#!/usr/bin/env python3
"""Graph-native mode diagnostic; this does not define physical state updates."""
import argparse
import json
from pathlib import Path

from face_energy_obstruction import derive_fiber_group
from run_shared_edge import mesh_geometry
from screen_braid_rules import nullspace


def transpose(matrix):
    return [list(row) for row in zip(*matrix)]


def multiply(left, right):
    if (not left or not right or len(left[0]) != len(right) or
            any(len(row) != len(left[0]) for row in left) or
            any(len(row) != len(right[0]) for row in right)):
        raise ValueError('incompatible matrix dimensions')
    columns = transpose(right)
    return [[sum(a*b for a, b in zip(row, column)) for column in columns] for row in left]


def identity(size):
    return [[int(i == j) for j in range(size)] for i in range(size)]


def representation(group):
    edges = ((0, 1), (1, 2), (2, 3), (3, 0))
    adjacency = [[0]*4 for _ in range(4)]
    for u, v in edges:
        adjacency[u][v] = adjacency[v][u] = 1
    basis = nullspace(adjacency, 4)
    q = transpose(basis)
    if len(basis) != 2 or multiply(basis, q) != [[2, 0], [0, 2]]:
        raise ValueError('unexpected square adjacency zero-mode basis')
    matrices = []
    for p in group.elements:
        perm = [[int(i == p[j]) for j in range(4)] for i in range(4)]
        projected = multiply(multiply(basis, perm), q)
        if any(x % 2 for row in projected for x in row):
            raise ValueError('mode action has nonintegral coefficients')
        r = [[x//2 for x in row] for row in projected]
        if multiply(perm, q) != multiply(q, r) or multiply(transpose(r), r) != identity(2):
            raise ValueError('fiber modes are not an invariant orthogonal subspace')
        matrices.append(r)
    for a in range(group.n):
        for b in range(group.n):
            if multiply(matrices[a], matrices[b]) != matrices[group.mul[a][b]]:
                raise ValueError('mode action is not the derived group representation')
    return {'fiber_edges': edges, 'adjacency': adjacency, 'basis_columns': q,
            'basis_gram': [[2, 0], [0, 2]], 'matrices': matrices}


def connection_laplacian(vertices, edges, links, matrices):
    if len(links) != len(edges) or any(not 0 <= x < len(matrices) for x in links):
        raise ValueError('invalid diagnostic link vector')
    ids = {v: i for i, v in enumerate(vertices)}
    if len(ids) != len(vertices) or len(set(edges)) != len(edges) or any(u >= v or u not in ids or v not in ids for u, v in edges):
        raise ValueError('invalid fixed labeled base graph')
    result = [[0]*(2*len(vertices)) for _ in range(2*len(vertices))]
    for (u, v), x in zip(edges, links):
        u, v = ids[u], ids[v]
        for a in range(2):
            result[2*u+a][2*u+a] += 1
            result[2*v+a][2*v+a] += 1
            for b in range(2):
                result[2*v+a][2*u+b] -= matrices[x][a][b]
                result[2*u+b][2*v+a] -= matrices[x][a][b]
    return result


def lifted_graph_laplacian(vertices, edges, links, group, fiber_edges):
    ids = {v: i for i, v in enumerate(vertices)}
    lifted = set()
    for v in range(len(vertices)):
        for a, b in fiber_edges:
            lifted.add(tuple(sorted((4*v+a, 4*v+b))))
    for (u, v), x in zip(edges, links):
        for a, b in enumerate(group.elements[x]):
            lifted.add(tuple(sorted((4*ids[u]+a, 4*ids[v]+b))))
    result = [[0]*(4*len(vertices)) for _ in range(4*len(vertices))]
    for a, b in lifted:
        result[a][a] += 1
        result[b][b] += 1
        result[a][b] -= 1
        result[b][a] -= 1
    return result


def verify_graph_restriction(vertices, edges, links, group, modes):
    operator = connection_laplacian(vertices, edges, links, modes['matrices'])
    full = lifted_graph_laplacian(vertices, edges, links, group, modes['fiber_edges'])
    embedding = [[0]*(2*len(vertices)) for _ in range(4*len(vertices))]
    for v in range(len(vertices)):
        for a in range(4):
            for b in range(2):
                embedding[4*v+a][2*v+b] = modes['basis_columns'][a][b]
    shifted = [row[:] for row in operator]
    for i in range(len(shifted)):
        shifted[i][i] += 2  # Vertical C4 Laplacian eigenvalue in this actual subspace.
    if multiply(full, embedding) != multiply(embedding, shifted):
        raise ValueError('mode operator is not the restriction of the lifted-graph Laplacian')
    return operator, full


def moments(matrix, powers=8):
    current, result = identity(len(matrix)), []
    for _ in range(powers):
        current = multiply(current, matrix)
        result.append(sum(current[i][i] for i in range(len(matrix))))
    return result


def analyze(witness):
    group = derive_fiber_group()
    modes = representation(group)
    side = witness['side']
    if not 3 <= side <= 4:
        raise ValueError('dense exact spectral witness is bounded to side 3 or 4')
    vertices = list(range(side*side))
    edges, _, _ = mesh_geometry(side)
    branches = []
    for branch in witness['branches']:
        operator, full = verify_graph_restriction(vertices, edges, branch['final_links'], group, modes)
        branches.append({'rule': branch['rule'], 'mode_laplacian': operator,
                         'mode_trace_powers_1_through_8': moments(operator),
                         'lifted_graph_trace_powers_1_through_8': moments(full),
                         'exact_lifted_graph_intertwining': True})
    contrast = [i+1 for i, (a, b) in enumerate(zip(branches[0]['mode_trace_powers_1_through_8'],
                                                 branches[1]['mode_trace_powers_1_through_8'])) if a != b]
    if not contrast:
        raise ValueError('no exact spectral contrast in the supplied witness')
    return {'schema': 1, 'source': 'data/d4-loop-observer.json conversion_witness', 'representation': modes,
            'branches': branches, 'first_distinct_mode_trace_power': min(contrast),
            'scope': 'Adjacency-derived spectral diagnostic on an unweighted lifted graph. No Hamiltonian, wave equation, matter field, Born rule, or coupling is inserted into simulation updates.'}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, default=Path('out/fiber-mode-probe.json'))
    args = parser.parse_args()
    from d4_loop_observer import conversion_witness
    witness = conversion_witness()
    saved = json.loads(Path('data/d4-loop-observer.json').read_text())['conversion_witness']
    if json.loads(json.dumps(witness)) != saved:
        raise ValueError('saved witness differs from independent raw-link replay')
    result = analyze(witness)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2)+'\n')
    print('First distinct exact mode moment:', result['first_distinct_mode_trace_power'])


if __name__ == '__main__':
    main()
