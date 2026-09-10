#!/usr/bin/env python3
"""Finite-group certificate for spectator-blind two-link face-energy obstruction.

For the four-face patch in shared-face-dynamics.md, all conserved real/rational
one-face sums are class functions of G/N, where N is the normal subgroup generated
by X A^-1 over F(A,B)=(X,Y). Every such observable is fixed face by face. This is
not an obstruction to larger-support observables or spectator-aware update laws.
"""
import argparse
import itertools
import json
from pathlib import Path
from screen_braid_rules import nullspace


class FiniteGroup:
    def __init__(self, permutations):
        self.elements = [tuple(p) for p in permutations]
        self.n = len(self.elements)
        degree = len(self.elements[0])
        if len(set(self.elements)) != self.n or any(sorted(p) != list(range(degree)) for p in self.elements):
            raise ValueError('expected distinct permutations of a common finite set')
        self.identity = self.elements.index(tuple(range(degree)))
        self.mul = [[self.elements.index(tuple(a[b[i]] for i in range(degree)))
                     for b in self.elements] for a in self.elements]
        self.inv = [next(b for b in range(self.n) if self.mul[a][b] == self.identity) for a in range(self.n)]
        self.conj = [[self.mul[g][self.mul[a][self.inv[g]]] for a in range(self.n)] for g in range(self.n)]
        self.sectors = [min(row[a] for row in self.conj) for a in range(self.n)]

    def normal_closure(self, generators):
        subgroup = {self.identity} | {row[a] for a in generators for row in self.conj}
        while True:
            expanded = subgroup | {self.mul[a][b] for a in subgroup for b in subgroup}
            if expanded == subgroup:
                return sorted(subgroup)
            subgroup = expanded

    def validate(self, table, strict=False):
        n = self.n
        if len(table) != n*n or any(not isinstance(x, int) or not 0 <= x < n*n for x in table):
            raise ValueError('pair table has invalid entries')
        # The theorem needs product preservation, not bijectivity or the braid law.
        for p, output in enumerate(table):
            a, b = divmod(p, n)
            x, y = divmod(output, n)
            if self.mul[a][b] != self.mul[x][y]:
                raise ValueError('ordered product is not preserved')
            if any(table[row[a]*n+row[b]] != row[x]*n+row[y] for row in self.conj):
                raise ValueError('pair table is not conjugation equivariant')
        if strict:
            if table[self.identity*n+self.identity] != self.identity*n+self.identity:
                raise ValueError('flat vacuum is not fixed')
            for p, output in enumerate(table):
                if table[output] != p:
                    raise ValueError('table is not an involution')
                reverse = lambda value: self.inv[value % n]*n+self.inv[value//n]
                if table[reverse(p)] != reverse(output):
                    raise ValueError('orientation covariance failed')
            for triple in itertools.product(range(n), repeat=3):
                left, right = list(triple), list(triple)
                for values, order in ((left, (0, 1, 0)), (right, (1, 0, 1))):
                    for i in order:
                        values[i:i+2] = divmod(table[values[i]*n+values[i+1]], n)
                if left != right:
                    raise ValueError('strict braid relation failed')


def derive_fiber_group():
    edges = ((0, 1), (1, 2), (2, 3), (3, 0))
    canonical = {frozenset(e) for e in edges}
    return FiniteGroup([p for p in itertools.permutations(range(4))
                        if {frozenset((p[u], p[v])) for u, v in edges} == canonical])


def certificate(group, table):
    group.validate(table)
    n, mul, inv = group.n, group.mul, group.inv
    # Varying an independent spectator C forces q(L^-1 C)-q(C) to be
    # constant. Summing over finite G makes that constant zero over R or Q.
    # Class invariance extends this to the normal closure of L=X A^-1.
    # Without independent spectators (or for torsion charges), this proof fails.
    increments = sorted({mul[output//n][inv[p//n]] for p, output in enumerate(table)})
    normal = group.normal_closure(increments)
    cosets = [min(mul[h][a] for h in normal) for a in range(n)]
    quotient_classes = [min(cosets[row[a]] for row in group.conj) for a in range(n)]
    labels = sorted(set(quotient_classes)-{quotient_classes[group.identity]})
    basis = [[int(label == q) for label in quotient_classes] for q in labels]
    for p, output in enumerate(table):
        a, b = divmod(p, n)
        x, y = divmod(output, n)
        if cosets[a] != cosets[x] or cosets[b] != cosets[y]:
            raise ValueError('derived quotient transport is not pointwise fixed')
    return {'increment_elements': increments, 'normal_subgroup': normal,
            'quotient_element_labels': cosets, 'quotient_conjugacy_labels': quotient_classes,
            'normalized_face_charge_dimension': len(basis), 'element_charge_basis': basis,
            'quotient_pair_components_pointwise_fixed': True}


def exhaustive_face_equations(group, table):
    """Independent four-holonomy equations, without using normal closure.

    Use the unsimplified based-pair formula, not the theorem's single L formula.
    Retain all distinct equations and check the proposed quotient basis against them.
    """
    n, mul, inv, sectors = group.n, group.mul, group.inv, group.sectors
    labels = sorted(set(sectors)-{sectors[group.identity]})
    equations = set()
    for a, b, c, d in itertools.product(range(n), repeat=4):
        transported_b = mul[inv[a]][mul[b][a]]
        x, y = divmod(table[a*n+transported_b], n)
        delta = mul[x][inv[a]]
        next_b = mul[a][mul[y][inv[a]]]
        after = (mul[a][delta], next_b, mul[inv[delta]][c], mul[inv[next_b]][mul[b][d]])
        before = (a, b, c, d)
        equations.add(tuple(sum(sectors[v] == label for v in after)
                            - sum(sectors[v] == label for v in before) for label in labels))
    basis = nullspace(sorted(equations), len(labels))
    predicted = certificate(group, table)
    lifted = [[0 if sector == sectors[group.identity] else q[labels.index(sector)] for sector in sectors]
              for q in basis]
    if len(basis) != predicted['normalized_face_charge_dimension']:
        raise ValueError('four-face nullity disagrees with normal-subgroup theorem')
    if any(q[mul[h][a]] != q[a] for q in lifted for h in predicted['normal_subgroup'] for a in range(n)):
        raise ValueError('four-face charge does not descend to the derived quotient')
    if any(sum(row[i]*q[labels[i]] for i in range(len(labels)))
           for row in equations for q in predicted['element_charge_basis']):
        raise ValueError('quotient charge fails independent four-face conservation')
    return {'states_checked': n**4, 'distinct_equations': len(equations),
            'nonidentity_sector_labels': labels, 'charge_basis': basis}


def analyze(census, minimal):
    group = derive_fiber_group()
    if (minimal['group_order'] != group.n or minimal['rule_count'] != len(minimal['rules'])
            or len({r['id'] for r in minimal['rules']}) != len(minimal['rules'])):
        raise ValueError('minimal-closure census metadata disagrees')
    strict = [s for s in census['solutions'] if s['strict_braid']]
    if not census['exhausted'] or len(strict) != census['strict_braid_solutions']:
        raise ValueError('requires exhausted strict-braid census')
    if len({s['id'] for s in strict}) != len(strict) or len({tuple(s['table']) for s in strict}) != len(strict):
        raise ValueError('duplicate strict rule ids or tables')
    classes = {}
    for solution in strict:
        group.validate(solution['table'], strict=True)
        result = certificate(group, solution['table'])
        key = tuple(result['normal_subgroup'])
        if key not in classes:
            classes[key] = {**result, 'solution_ids': [], 'subgroup_growing_rules': 0}
            # Increment lists vary between rules sharing a normal closure.
            del classes[key]['increment_elements']
        classes[key]['solution_ids'].append(solution['id'])
        classes[key]['subgroup_growing_rules'] += bool(solution['subgroup_grown_pairs'])
    checks = []
    for rule in minimal['rules']:
        checks.append({'minimal_rule_id': rule['id'], **certificate(group, rule['table']),
                       'independent_four_face_equations': exhaustive_face_equations(group, rule['table'])})
    return {'schema': 1, 'group_permutations': group.elements, 'strict_rules_checked': len(strict),
            'scope': 'unrestricted four-face spectator patch and its embeddings; common rational class function q with q(identity)=0',
            'theorem': 'conserved one-face sums are exactly class functions of G/N, and are fixed face by face; N is the normal closure of X A^-1',
            'strict_rule_classes': list(classes.values()), 'minimal_rule_controls': checks,
            'control_states_checked': sum(c['independent_four_face_equations']['states_checked'] for c in checks)}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, default=Path('out/face-energy-obstruction.json'))
    args = parser.parse_args()
    result = analyze(json.loads(Path('data/d4-involution-braid-search.json').read_text()),
                     json.loads(Path('data/d4-pair-involutions.json').read_text()))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2)+'\n')
    print('Strict rules:', result['strict_rules_checked'], 'independent control states:', result['control_states_checked'])
    for row in result['strict_rule_classes']:
        print('N:', row['normal_subgroup'], 'face-charge dimension:', row['normalized_face_charge_dimension'],
              'rules:', len(row['solution_ids']))


if __name__ == '__main__':
    main()
