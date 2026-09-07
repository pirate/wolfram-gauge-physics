#!/usr/bin/env python3
"""Derive positivity certificates and activation channels across minimal rules."""
import argparse
import itertools
import json
import math
from collections import Counter
from pathlib import Path

from face_energy_obstruction import derive_fiber_group
from triple_rule_search import algebra, assess, make_table, minimal_closures
from screen_braid_rules import nullspace


def subgroup(group, generators):
    found = {group.identity, *generators}
    while True:
        following = found | {group.inv[a] for a in found} | {group.mul[a][b] for a in found for b in found}
        if following == found:
            return tuple(sorted(found))
        found = following


def primitive_row(row):
    scale = math.gcd(*row)
    if not scale:
        return tuple(row)
    row = tuple(x//scale for x in row)
    return row if next(x for x in row if x) > 0 else tuple(-x for x in row)


def positive_certificate(equations, width):
    """Exact rank-at-most-one feasibility, not a bounded integer mass search.

    Minimal closures here have one population equation up to sign. A strictly
    positive solution exists iff that row has both signs (or vanishes).
    """
    if width <= 0 or any(len(row) != width for row in equations):
        raise ValueError('invalid conservation-equation dimensions')
    rows = sorted({primitive_row(row) for row in equations if any(row)})
    if len(rows) > 1:
        raise ValueError('positivity oracle requires at most one independent equation')
    if not rows:
        return {'exists': True, 'weights': [1]*width, 'equation': [0]*width}
    row = rows[0]
    positive, negative = sum(x for x in row if x > 0), -sum(x for x in row if x < 0)
    if not negative:
        return {'exists': False, 'weights': None, 'nonnegative_obstruction': row}
    weights = [negative if x > 0 else positive if x < 0 else 1 for x in row]
    scale = math.gcd(*weights)
    weights = [x//scale for x in weights]
    if any(x <= 0 for x in weights) or sum(a*b for a, b in zip(row, weights)):
        raise ValueError('invalid exact positivity witness')
    return {'exists': True, 'weights': weights, 'equation': row}


def analyze(census):
    g = derive_fiber_group()
    states, products, action, _ = algebra(g)
    closures = minimal_closures(g)
    if [tuple(map(tuple, r['transpositions'])) for r in census['minimal_rules']] != closures:
        raise ValueError('channel census does not cover the independent minimal closures')
    labels = sorted(set(g.sectors)-{g.identity})
    populations = [tuple(sum(g.sectors[a] == s for a in state) for s in labels) for state in states]
    supports = [sum(a != g.identity for a in state) for state in states]
    generated = [subgroup(g, state) for state in states]
    central = {a for a in range(g.n) if all(g.mul[a][b] == g.mul[b][a] for b in range(g.n))}
    stabilizers = [tuple(h for h, op in enumerate(action) if op[x] == x) for x in range(len(states))]
    records = []
    for i, pairs in enumerate(closures):
        table = make_table(pairs)
        result = assess(g, table)
        equations = [tuple(a-b for a, b in zip(populations[x], populations[y])) for x, y in pairs]
        positive = positive_certificate(equations, len(labels))
        channels = []
        for x, y in pairs:
            if stabilizers[x] != stabilizers[y]:
                raise ValueError('a rule loses its exact frame stabilizer')
            if populations[x] == populations[y]:
                continue
            for source, target in ((x, y), (y, x)):
                channels.append({'input': source, 'output': target,
                                 'input_triple': states[source], 'output_triple': states[target],
                                 'input_nonflat': supports[source], 'output_nonflat': supports[target],
                                 'input_population': populations[source], 'output_population': populations[target],
                                 'input_subgroup': generated[source], 'output_subgroup': generated[target],
                                 'input_noncentral': sum(a not in central for a in states[source]),
                                 'stabilizer': stabilizers[source]})
        two_body = [c for c in channels if c['input_nonflat'] == c['input_noncentral'] == 2
                    and len(c['input_subgroup']) == g.n]
        records.append({'id': i, 'moved_states': result['moved_states'], 'transpositions': pairs,
                        'charge_basis': result['charge_basis'], 'positive_certificate': positive,
                        'class_conversion_channels': channels,
                        'noncommuting_two_noncentral_inputs': two_body,
                        'same_classes_and_boundary_feedback': result['same_classes_and_boundary_feedback']})
    eligible = [r for r in records if r['positive_certificate']['exists'] and r['noncommuting_two_noncentral_inputs']]
    minimum = min(r['moved_states'] for r in eligible) if eligible else None
    selected = [r['id'] for r in eligible if r['moved_states'] == minimum]
    counts = Counter((r['positive_certificate']['exists'], bool(r['class_conversion_channels']),
                      bool(r['noncommuting_two_noncentral_inputs'])) for r in records)
    families = []
    for equation in sorted({tuple(r['positive_certificate']['equation']) for r in eligible}):
        members = [r for r in eligible if tuple(r['positive_certificate']['equation']) == equation]
        equations, staggered = [], []
        for r in members:
            for x, y in r['transpositions']:
                equations.append(tuple(a-b for a, b in zip(populations[x], populations[y])))
                for phase in (0, 1):
                    parts = (phase, 1-phase, phase)
                    staggered.append(tuple(sum((g.sectors[a] == s)-(g.sectors[b] == s)
                                               for a, b, p in zip(states[x], states[y], parts) if p == part)
                                           for part in (0, 1) for s in labels))
        # A distinct, lower-arity transport primitive: (1,g) <-> (g,1).
        pair_transport = list(range(g.n*g.n))
        for a in range(1, g.n):
            pair_transport[a], pair_transport[a*g.n] = a*g.n, a
        mixed = staggered[:]
        for x, y in enumerate(pair_transport):
            before, after = divmod(x, g.n), divmod(y, g.n)
            mixed.append(tuple((g.sectors[before[p]] == s)-(g.sectors[after[p]] == s)
                               for p in (0, 1) for s in labels))
        families.append({'id': len(families), 'equation': equation,
                         'rule_ids': [r['id'] for r in members],
                         'charge_basis': nullspace(equations, len(labels)),
                         'positive_certificate': positive_certificate(equations, len(labels)),
                         'triple_only_staggered_basis': nullspace(staggered, 2*len(labels)),
                         'shared_edge_transport_staggered_basis': nullspace(mixed, 2*len(labels))})
        tables = [make_table(r['transpositions']) for r in members]
        classes = [tuple(g.sectors[a] for a in state) for state in states]
        rates, rate_witness = {}, None
        for x, observed in enumerate(classes):
            counts_out = tuple(sorted(Counter(classes[table[x]] for table in tables).items()))
            if observed in rates and rates[observed][1] != counts_out:
                if rate_witness is None:
                    previous, old = rates[observed]
                    rate_witness = {'inputs': [previous, x], 'input_triples': [states[previous], states[x]],
                                    'class_tuple': observed, 'outgoing_counts': [old, counts_out]}
            else:
                rates.setdefault(observed, (x, counts_out))
        families[-1]['uniform_bank_class_rate_nonclosure'] = rate_witness
        families[-1]['uniform_bank_class_rates'] = None if rate_witness else [
            {'input': observed, 'outputs': rate[1]} for observed, rate in sorted(rates.items())]
        if rate_witness is None:
            sizes = Counter(g.sectors)
            degeneracy = lambda values: math.prod(sizes[a] for a in values)
            checked = 0
            for before, (_, outputs) in rates.items():
                for after, count in outputs:
                    back = dict(rates[after][1]).get(before, 0)
                    if degeneracy(before)*count != degeneracy(after)*back:
                        raise ValueError('coarse channel rates violate degeneracy-weighted detailed balance')
                    checked += 1
            families[-1]['coarse_detailed_balance'] = {'class_sizes': dict(sorted(sizes.items())),
                                                       'directed_entries_checked': checked,
                                                       'weight': 'product of conjugacy-class sizes'}
            nonzero = [j for j, a in enumerate(equation) if a]
            heavy = next(j for j in nonzero if sum((equation[j] > 0) == (equation[k] > 0) for k in nonzero) == 1)
            light = [j for j in nonzero if j != heavy]
            families[-1]['reaction_species'] = {'heavy_class': labels[heavy], 'light_classes': [labels[j] for j in light],
                                                 'grand_canonical_occupation_ratio': [math.prod(sizes[labels[j]] for j in light), sizes[labels[heavy]]*sizes[g.identity]],
                                                 'ratio_convention': 'rho_a rho_b / (rho_h rho_flat) for the stationary product measure; not products of correlated canonical marginals'}
    incompatibilities = []
    for a, b in itertools.combinations(families, 2):
        basis = nullspace([a['equation'], b['equation']], len(labels))
        forced_zero = [labels[j] for j in range(len(labels)) if all(q[j] == 0 for q in basis)]
        if not forced_zero:
            raise ValueError('different family cones do not have the claimed strict-positive obstruction')
        incompatibilities.append({'families': [a['id'], b['id']], 'joint_charge_basis': basis,
                                  'forced_zero_class_weights': forced_zero})
    return {'schema': 1, 'group_permutations': g.elements, 'class_labels': labels,
            'selection': 'all least-support minimal rules with a strictly positive face weight and a class-converting channel from two noncentral holonomies generating the full group',
            'scope': 'positive weight is a feasibility certificate, not a selected mass or energy; activation counts occupied faces, not particles',
            'counts': [{'positive': p, 'class_conversion': c, 'two_noncentral_full_group': t, 'rules': n}
                       for (p, c, t), n in sorted(counts.items())],
            'eligible_rules': len(eligible), 'minimum_moved_states': minimum, 'selected_rule_ids': selected,
            'families': families, 'family_incompatibilities': incompatibilities,
            'shared_edge_vacancy_transport': pair_transport,
            'rules': records}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, default=Path('out/triple-channels.json'))
    args = parser.parse_args()
    result = analyze(json.loads(Path('data/d4-triple-rule-search.json').read_text()))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2)+'\n')
    print('Counts:', result['counts'])
    print('Eligible:', result['eligible_rules'], 'minimum support:', result['minimum_moved_states'])
    for family in result['families']:
        print('Family:', family['id'], 'rules:', len(family['rule_ids']), 'charge equation:', family['equation'],
              'staggered dimensions:', len(family['triple_only_staggered_basis']),
              'with shared-edge transport:', len(family['shared_edge_transport_staggered_basis']),
              'uniform class-rate closed:', family['uniform_bank_class_rate_nonclosure'] is None)


if __name__ == '__main__':
    main()
