#!/usr/bin/env python3
"""Classify every strict braid solution's maximal additive class-charge factor.

The input contains all isolated-triple solutions; no minimum-support or subgroup
growth filter is applied here. Charges are solved afresh for each exact solution.
"""
import argparse
import itertools
import json
from collections import Counter
from pathlib import Path
from charge_quotient import induced_table, classify_block_exchange
from screen_braid_rules import nullspace


def analyze(data, minimal):
    edges = {frozenset(e) for e in ((0, 1), (1, 2), (2, 3), (3, 0))}
    group = [p for p in itertools.permutations(range(4))
             if {frozenset((p[a], p[b])) for a, b in ((0, 1), (1, 2), (2, 3), (3, 0))} == edges]
    n = len(group)
    identity = group.index(tuple(range(4)))
    multiply = [[group.index(tuple(a[b[i]] for i in range(4))) for b in group] for a in group]
    inverse = [next(b for b in range(n) if multiply[a][b] == identity) for a in range(n)]
    conjugate = [[multiply[g][multiply[a][inverse[g]]] for a in range(n)] for g in range(n)]
    sectors = [min(row[a] for row in conjugate) for a in range(n)]
    nonidentity = sorted(set(sectors)-{identity})
    minimal_factors = []
    if minimal['group_order'] != n or minimal['rule_count'] != len(minimal['rules']):
        raise ValueError('minimal-closure census metadata disagrees')
    for rule in minimal['rules']:
        quotient = induced_table(rule['table'], sectors)
        if quotient is None:
            raise ValueError(f'minimal closure lacks sector factor: {rule["id"]}')
        minimal_factors.append({'minimal_rule_id': rule['id'], 'sector_factor': quotient})
    classes, counts = {}, Counter()
    strict = [s for s in data['solutions'] if s['strict_braid']]
    if not data['exhausted'] or len(strict) != data['strict_braid_solutions']:
        raise ValueError('requires the complete strict-braid census')
    for solution in strict:
        table = solution['table']
        if len(table) != n*n or sorted(table) != list(range(n*n)) or table[0] != 0:
            raise ValueError('invalid flat-fixing pair bijection')
        equations = []
        for a, b in itertools.product(range(n), repeat=2):
            c, d = divmod(table[a*n+b], n)
            if table[c*n+d] != a*n+b or multiply[a][b] != multiply[c][d]:
                raise ValueError('involution or product constraint failed')
            if table[inverse[b]*n+inverse[a]] != inverse[d]*n+inverse[c]:
                raise ValueError('orientation reversal constraint failed')
            if any(table[g[a]*n+g[b]] != g[c]*n+g[d] for g in conjugate):
                raise ValueError('conjugation covariance failed')
            equations.append([int(sectors[a] == s)+int(sectors[b] == s)
                              -int(sectors[c] == s)-int(sectors[d] == s) for s in nonidentity])
        for triple in itertools.product(range(n), repeat=3):
            left, right = list(triple), list(triple)
            for state, order in ((left, (0, 1, 0)), (right, (1, 0, 1))):
                for i in order:
                    state[i:i+2] = divmod(table[state[i]*n+state[i+1]], n)
            if left != right:
                raise ValueError('strict braid assertion failed')
        basis = nullspace(equations, len(nonidentity))
        vectors = [tuple(0 if s == identity else q[nonidentity.index(s)] for q in basis) for s in sectors]
        alphabet = sorted(set(vectors))
        labels = [alphabet.index(v) for v in vectors]
        quotient = induced_table(table, labels)
        blocks = classify_block_exchange(quotient, list(range(len(alphabet))))
        sector_closed = induced_table(table, sectors) is not None
        if not sector_closed or blocks is None:
            raise ValueError(f'classification counterexample: solution {solution["id"]}')
        signature = (tuple(map(tuple, basis)), tuple(map(tuple, quotient)))
        if signature not in classes:
            classes[signature] = {'charge_basis': basis, 'charge_vectors': alphabet,
                                  'element_symbols': labels, 'quotient_table': quotient,
                                  'exchange_blocks': blocks, 'solution_ids': []}
        classes[signature]['solution_ids'].append(solution['id'])
        counts[len(basis), len(alphabet), tuple(sorted(map(len, blocks)))] += 1
    return {'schema': 1, 'scope': 'all exact braid solutions in the exhausted square-fiber involution census',
            'claim': 'every maximal additive class-charge projection is an autonomous block-exchange gate',
            'strict_rules_checked': len(strict), 'sector_closed_rules': len(strict),
            'minimal_closure_sector_factors': minimal_factors,
            'extension_to_all_admitted_involutions':
                'every admitted table is a disjoint union of the minimal closures; '
                'a disjoint union equals their composition, and autonomous factors compose',
            'group_permutations': group, 'element_sectors': sectors,
            'nonidentity_sector_labels': nonidentity,
            'charge_convention': 'all rational additive class functions with q(identity)=0',
            'counts': [{'charge_rank': rank, 'symbols': symbols, 'block_sizes': sizes, 'rules': count}
                       for (rank, symbols, sizes), count in sorted(counts.items())],
            'quotient_classes': list(classes.values())}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', type=Path, default=Path('data/d4-involution-braid-search.json'))
    parser.add_argument('--output', type=Path, default=Path('out/braid-quotient-classification.json'))
    parser.add_argument('--minimal', type=Path, default=Path('data/d4-pair-involutions.json'))
    args = parser.parse_args()
    result = analyze(json.loads(args.input.read_text()), json.loads(args.minimal.read_text()))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(f'{result["strict_rules_checked"]} strict rules independently verified; '
          f'{len(result["quotient_classes"])} labeled charge quotients; all block exchanges')


if __name__ == '__main__':
    main()
