#!/usr/bin/env python3
"""Exact autonomous-factor test for measured one-cell observables."""
import argparse
import json
from pathlib import Path


def induced_table(table, labels):
    n = len(labels)
    alphabet = sorted(set(labels))
    targets = {}
    for a in range(n):
        for b in range(n):
            c, d = divmod(table[a*n+b], n)
            key, value = (labels[a], labels[b]), (labels[c], labels[d])
            if key in targets and targets[key] != value:
                return None
            targets[key] = value
    return [targets[a, b] for a in alphabet for b in alphabet]


def classify_exclusion(table, alphabet):
    table = [tuple(pair) for pair in table]
    for distinguished in alphabet:
        expected = [(b, a) if (a == distinguished or b == distinguished) else (a, b)
                    for a in alphabet for b in alphabet]
        if table == expected:
            return distinguished
    return None


def classify_block_exchange(table, alphabet):
    """Return blocks iff F(a,b) swaps precisely between distinct blocks."""
    if table is None:
        return None
    lookup = dict(zip(((a, b) for a in alphabet for b in alphabet), map(tuple, table)))
    blocks = {a: tuple(b for b in alphabet if lookup[a, b] == (a, b)) for a in alphabet}
    if any(lookup[a, b] != ((a, b) if blocks[a] == blocks[b] else (b, a))
           for a in alphabet for b in alphabet):
        return None
    return sorted(set(blocks.values()))


def analyze(census, screen):
    group = [tuple(p) for p in screen['group_permutations']]
    n = len(group)
    multiply = [[group.index(tuple(a[b[i]] for i in range(len(a)))) for b in group] for a in group]
    identity = group.index(tuple(range(len(group[0]))))
    inverse = [next(b for b in range(n) if multiply[a][b] == identity) for a in range(n)]
    sectors = [min(multiply[g][multiply[a][inverse[g]]] for g in range(n)) for a in range(n)]
    results = []
    for record in screen['selected_rules']:
        table = census['solutions'][record['solution_id']]['table']
        vectors = [tuple(0 if sector == identity else q[screen['nonidentity_sector_labels'].index(sector)]
                         for q in record['charge_basis']) for sector in sectors]
        alphabet = sorted(set(vectors))
        labels = [alphabet.index(v) for v in vectors]
        quotient = induced_table(table, labels)
        full_sector = induced_table(table, sectors)
        results.append({'solution_id': record['solution_id'], 'charge_vectors': alphabet,
                        'element_symbols': labels, 'sector_factor_closed': full_sector is not None,
                        'charge_factor_closed': quotient is not None,
                        'quotient_table': quotient,
                        'exclusion_distinguished_symbol': classify_exclusion(quotient, range(len(alphabet)))
                        if quotient is not None else None})
    return {'schema': 1, 'claim': 'exact local factor maps; therefore valid for every sequence of supported cell updates',
            'rules': results}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, default=Path('out/charge-quotient.json'))
    args = parser.parse_args()
    result = analyze(json.loads(Path('data/d4-involution-braid-search.json').read_text()),
                     json.loads(Path('data/d4-braid-rule-screen.json').read_text()))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(f'{len(result["rules"])} exact quotient checks; '
          f'{sum(r["exclusion_distinguished_symbol"] is not None for r in result["rules"])} exclusion factors')


if __name__ == '__main__':
    main()
