#!/usr/bin/env python3
"""Independent minimal triple-closure search and derived-charge feedback audit."""
import argparse
import itertools
import json
import subprocess
from collections import Counter
from pathlib import Path
from face_energy_obstruction import derive_fiber_group
from screen_braid_rules import nullspace


def algebra(group):
    states = list(itertools.product(range(group.n), repeat=3))
    ids = {state: i for i, state in enumerate(states)}
    products = [group.mul[group.mul[a][b]][c] for a, b, c in states]
    action = [[ids[tuple(row[a] for a in state)] for state in states] for row in group.conj]
    reverse = [ids[tuple(group.inv[a] for a in state[::-1])] for state in states]
    return states, products, action, reverse


def minimal_closures(group):
    states, products, action, reverse = algebra(group)
    operations = action+[[reverse[i] for i in row] for row in action]
    flat = states.index((group.identity,)*3)
    found = set()
    for x in range(len(states)):
        for y in range(x+1, len(states)):
            if flat in (x, y) or products[x] != products[y]: continue
            pairs = tuple(sorted({tuple(sorted((op[x], op[y]))) for op in operations}))
            if len({z for pair in pairs for z in pair}) == 2*len(pairs):
                found.add(pairs)
    return sorted(found)


def make_table(pairs, size=512):
    table = list(range(size))
    for a, b in pairs:
        if not 0 <= a < b < size or table[a] != a or table[b] != b:
            raise ValueError('transpositions overlap or have invalid entries')
        table[a], table[b] = b, a
    return table


def assess(group, table):
    states, products, action, reverse = algebra(group)
    size = len(states)
    if len(table) != size or sorted(table) != list(range(size)):
        raise ValueError('triple map is not a permutation')
    if table[states.index((group.identity,)*3)] != states.index((group.identity,)*3):
        raise ValueError('flat vacuum moved')
    for x, y in enumerate(table):
        if table[y] != x or products[x] != products[y] or table[reverse[x]] != reverse[y]:
            raise ValueError('triple involution, product, or orientation constraint failed')
        if any(table[op[x]] != op[y] for op in action):
            raise ValueError('triple gauge covariance failed')
    labels = sorted(set(group.sectors)-{group.identity})
    sectors = [tuple(group.sectors[a] for a in s) for s in states]
    equations = sorted({tuple(sectors[x].count(s)-sectors[y].count(s) for s in labels)
                        for x, y in enumerate(table)})
    basis = nullspace(equations, len(labels))
    charges = [[0 if sector == group.identity else q[labels.index(sector)] for sector in group.sectors] for q in basis]
    vectors = [tuple(q[a] for q in charges) for a in range(group.n)]
    observed = [tuple(vectors[a] for a in s) for s in states]
    def witness(keys):
        seen = {}
        for x, key in enumerate(keys):
            output = observed[table[x]]
            if key in seen and seen[key][1] != output:
                left = seen[key][0]
                if any(op[left] == x for op in action):
                    raise ValueError('feedback witness consists of gauge copies')
                return {'inputs': [left, x], 'outputs': [table[left], table[x]],
                        'input_triples': [states[left], states[x]], 'output_triples': [states[table[left]], states[table[x]]],
                        'input_charge_vectors': [observed[left], observed[x]],
                        'output_charge_vectors': [observed[table[left]], output]}
            seen[key] = (x, output)
        return None
    return {'charge_basis': basis, 'element_charges': charges,
            'moved_states': sum(x != y for x, y in enumerate(table)),
            'transported_charge': any(observed[x] != observed[y] for x, y in enumerate(table)),
            'charge_nonclosure': witness(observed),
            'boundary_charge_nonclosure': witness(list(zip(observed, products))),
            'same_classes_and_boundary_feedback': witness(list(zip(sectors, products)))}


def analyze(raw):
    group = derive_fiber_group()
    closures = minimal_closures(group)
    recorded = [tuple(map(tuple, r['transpositions'])) for r in raw['minimal_rules']]
    if raw['group_order'] != group.n or recorded != closures or [r['id'] for r in raw['minimal_rules']] != list(range(len(closures))):
        raise ValueError('C++ breadth-first closures disagree with independent Python action orbits')
    records, counts = [], Counter()
    for i, pairs in enumerate(closures):
        result = assess(group, make_table(pairs))
        counts[len(result['charge_basis']), result['transported_charge'], result['charge_nonclosure'] is None] += 1
        records.append({'id': i, 'transpositions': pairs, **result})
    eligible = [r for r in records if r['same_classes_and_boundary_feedback'] is not None]
    selected = min(eligible, key=lambda r: (r['moved_states'], r['transpositions']))
    # A separate vacuum-transport control, composed only after candidate selection.
    n = group.n
    transport_pairs = [(a, a*n*n) for a in range(1, n)]
    transport = make_table(transport_pairs)
    combined_pairs = sorted(list(selected['transpositions'])+transport_pairs)
    combined = make_table(combined_pairs)
    return {'schema': 1, 'scope': 'all nonidentity minimal product-preserving gauge/orientation transposition closures on G^3; not all composite triple rules',
            'group_permutations': group.elements, 'nonidentity_sector_labels': sorted(set(group.sectors)-{group.identity}),
            'minimal_rules': records, 'minimal_rule_count': len(records),
            'counts': [{'charge_rank': rank, 'transported': transported, 'charge_observer_closed': closed, 'rules': count}
                       for (rank, transported, closed), count in sorted(counts.items())],
            'same_classes_and_boundary_feedback_rules': len(eligible),
            'selection': 'minimum moved-state count, then lexicographic transpositions, among same-class same-boundary outgoing-charge feedback witnesses',
            'selected_rule': {**selected, 'table': make_table(selected['transpositions'])},
            'transport_control': {'transpositions': transport_pairs, 'table': transport, **assess(group, transport)},
            'combined_rule': {'transpositions': combined_pairs, 'table': combined, **assess(group, combined)},
            'limitations': 'no triple braid or overlapping-update schedule-consistency theorem; derived charges are not identified with physical energy'}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, default=Path('out/triple-rule-search.json'))
    args = parser.parse_args()
    process = subprocess.run(['build/wgphysics_triple_census'], text=True, capture_output=True, check=True)
    result = analyze(json.loads(process.stdout))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2)+'\n')
    print('Minimal closures:', result['minimal_rule_count'], 'strong feedback witnesses:', result['same_classes_and_boundary_feedback_rules'])
    print('Selected pairs:', result['selected_rule']['transpositions'], 'charges:', result['selected_rule']['charge_basis'])


if __name__ == '__main__':
    main()
