#!/usr/bin/env python3
"""Independent Python screen of the least-support strict, subgroup-growing rules.

No spatial template or continuum target enters selection. Fractions recover all
additive class-function charges with q(identity)=0 from exact conservation equations.
"""
import argparse
import itertools
import json
import math
from fractions import Fraction
from pathlib import Path


def nullspace(matrix, columns):
    rows = [[Fraction(x) for x in row] for row in matrix if any(row)]
    pivots = []
    for column in range(columns):
        pivot = next((i for i in range(len(pivots), len(rows)) if rows[i][column]), None)
        if pivot is None:
            continue
        index = len(pivots)
        rows[index], rows[pivot] = rows[pivot], rows[index]
        scale = rows[index][column]
        rows[index] = [value / scale for value in rows[index]]
        for i in range(len(rows)):
            if i != index and rows[i][column]:
                factor = rows[i][column]
                rows[i] = [a - factor * b for a, b in zip(rows[i], rows[index])]
        pivots.append(column)
    basis = []
    for free in sorted(set(range(columns)) - set(pivots)):
        vector = [Fraction(0)] * columns
        vector[free] = 1
        for i, pivot in enumerate(pivots):
            vector[pivot] = -rows[i][free]
        scale = math.lcm(*(v.denominator for v in vector))
        basis.append([int(v * scale) for v in vector])
    return basis


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', type=Path, default=Path('data/d4-involution-braid-search.json'))
    parser.add_argument('--output', type=Path, default=Path('out/braid-rule-screen.json'))
    args = parser.parse_args()
    data = json.loads(args.input.read_text())
    edges = {frozenset(edge) for edge in ((0, 1), (1, 2), (2, 3), (3, 0))}
    group = [p for p in itertools.permutations(range(4))
             if {frozenset((p[a], p[b])) for a, b in ((0, 1), (1, 2), (2, 3), (3, 0))} == edges]
    n = len(group)
    identity = group.index(tuple(range(4)))
    multiply = [[group.index(tuple(a[b[i]] for i in range(4))) for b in group] for a in group]
    inverse = [next(b for b in range(n) if multiply[a][b] == identity) for a in range(n)]
    sector = [min(multiply[g][multiply[a][inverse[g]]] for g in range(n)) for a in range(n)]
    labels = sorted(set(sector) - {identity})
    candidates = [s for s in data['solutions'] if s['strict_braid'] and s['subgroup_grown_pairs']]
    complexity = lambda s: sum(i != y for i, y in enumerate(s['table']))
    minimum = min(map(complexity, candidates))
    selected = [s for s in candidates if complexity(s) == minimum]
    records = []
    for candidate in selected:
        table = candidate['table']
        if len(table) != n*n:
            raise ValueError('table and independently derived fiber group disagree')
        def update(state, i):
            state[i], state[i+1] = divmod(table[state[i]*n+state[i+1]], n)
        # Recheck every strict braid triple independently of the C++ search.
        for triple in itertools.product(range(n), repeat=3):
            left, right = list(triple), list(triple)
            for i in (0, 1, 0): update(left, i)
            for i in (1, 0, 1): update(right, i)
            if left != right:
                raise AssertionError('candidate fails independent braid check')
        equations = []
        for a, b in itertools.product(range(n), repeat=2):
            c, d = divmod(table[a*n+b], n)
            if table[c*n+d] != a*n+b:
                raise AssertionError('independent involution check failed')
            for g in range(n):
                conjugate = lambda value: multiply[g][multiply[value][inverse[g]]]
                if table[conjugate(a)*n+conjugate(b)] != conjugate(c)*n+conjugate(d):
                    raise AssertionError('independent gauge check failed')
            if multiply[a][b] != multiply[c][d]:
                raise AssertionError('product conservation failed')
            equations.append([int(sector[a] == s)+int(sector[b] == s)
                              -int(sector[c] == s)-int(sector[d] == s) for s in labels])
        basis = nullspace(equations, len(labels))
        if any(sum(a*b for a, b in zip(row, q)) for row in equations for q in basis):
            raise AssertionError('charge nullspace failed substitution')
        charge = [[0 if s == identity else q[labels.index(s)] for s in sector] for q in basis]
        runs = []
        for seed in range(n):
            if seed == identity: continue
            for phase in (0, 1):
                state = [identity]*256
                state[128] = seed
                initial_charge = [sum(q[a] for a in state) for q in charge]
                samples = []
                for tick in range(65):
                    if tick:
                        for i in range((tick-1+phase) % 2, 255, 2): update(state, i)
                    if [sum(q[a] for a in state) for q in charge] != initial_charge:
                        raise AssertionError('measured charge drift')
                    if tick % 8 == 0:
                        occupied = [i-128 for i, a in enumerate(state) if a != identity]
                        samples.append({'layer': tick, 'active': len(occupied),
                                        'left': min(occupied, default=0), 'right': max(occupied, default=0)})
                runs.append({'seed_element': seed, 'phase': phase, 'samples': samples})
        collisions = []
        # Distinguishable conserved species on opposite sublattices, including
        # stationary controls for laws under which the selected seeds do not move.
        # Single-defect references use identical sites and the identical schedule.
        for first_seed, second_seed in itertools.permutations(labels[:3], 2):
            together = [identity]*512
            alone_first, alone_second = together[:], together[:]
            together[224] = alone_first[224] = first_seed
            together[255] = alone_second[255] = second_seed
            positions = []
            for tick in range(97):
                if tick:
                    for state in (together, alone_first, alone_second):
                        for i in range((tick-1) % 2, 511, 2): update(state, i)
                row = [tick]
                for state, target in ((together, first_seed), (together, second_seed),
                                      (alone_first, first_seed), (alone_second, second_seed)):
                    sites = [i for i, value in enumerate(state) if sector[value] == target]
                    if len(sites) != 1:
                        raise AssertionError('collision did not preserve a single defect of each species')
                    row.append(sites[0])
                positions.append(row)
            shifts = [(p[1]-p[3], p[2]-p[4]) for p in positions[-16:]]
            contact = next((p[0] for p in positions if abs(p[1]-p[2]) == 1), None)
            collisions.append({'seeds': [first_seed, second_seed],
                               'first_contact_layer': contact,
                               'incoming_reference_velocity': [(positions[8][i]-positions[0][i])/8 for i in (3, 4)],
                               'outgoing_velocity': [(positions[-1][i]-positions[-17][i])/16 for i in (1, 2)],
                               'constant_late_shift': list(shifts[0]) if len(set(shifts)) == 1 else None,
                               'trajectory_columns': ['layer', 'first', 'second', 'first_alone', 'second_alone'],
                               'trajectory': positions[::4]})
        records.append({'solution_id': candidate['id'], 'charge_basis': basis,
                        'runs': runs, 'collisions': collisions})
    witness = data['spectator_counterexample']
    def equivalent(left, right):
        return any(all(multiply[g][multiply[a][inverse[g]]] == b for a, b in zip(left, right))
                   for g in range(n))
    if not equivalent(witness['left'][:3], witness['right'][:3]) or equivalent(witness['left'], witness['right']):
        raise AssertionError('spectator counterexample failed independent gauge test')
    result = {'schema': 1, 'selection': 'all strict subgroup-growing solutions with minimum moved-pair count',
              'moved_pairs': minimum, 'group_permutations': group, 'nonidentity_sector_labels': labels,
              'charge_convention': 'rational additive class functions, q(identity)=0',
              'cells': 256, 'layers': 64, 'selected_rules': records}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(f'{len(records)} minimal strict rules; {len(records)*14} seed runs; '
          f'{len(records)*6} two-defect experiments including no-contact controls; charges and braid checks passed')


if __name__ == '__main__':
    main()
