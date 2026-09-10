#!/usr/bin/env python3
"""Exact prime-power resonance arithmetic and sparse sharp-arity reactions.

Scientific witnesses only: the spectral diagnostic is not installed as an
energy law or a replacement rule bank in the evolution engine.
"""
import itertools
import json
from collections import Counter, defaultdict
from pathlib import Path

from cycle_charge_compatibility import cycle_group
from run_three_face import LinkOracle
from screen_braid_rules import nullspace


def arithmetic(p, exponent):
    assert p >= 3 and all(p % d for d in range(2, int(p**0.5)+1))
    assert exponent >= 2
    n, m = p**exponent, p**(exponent-1)
    g = cycle_group(n)
    t = lambda k: g.elements.index(tuple((j+k) % n for j in range(n)))
    rotations = {t(k): k for k in range(n)}

    def energy(word):
        coefficients = [0]*n
        for h in word:
            coefficients[0] += 2
            if h in rotations:
                k = rotations[h]
                coefficients[k] -= 1
                coefficients[-k % n] -= 1
        # Phi_(p^a)(z) = sum_(j=0)^(p-1) z^(j*m). Independent
        # remainders for each residue mod m; integer arithmetic throughout.
        return tuple(coefficients[k]-coefficients[(p-1)*m+k % m]
                     for k in range((p-1)*m))

    def product(word):
        value = g.identity
        for h in word:
            value = g.mul[value][h]
        return value

    special = sorted({min(j*m, n-j*m) for j in range(1, p)})
    bands = [sorted({min((r+j*m) % n, -(r+j*m) % n) for j in range(p)})
             for r in range(1, (m+1)//2)]
    assert len(special) == (p-1)//2 and all(len(b) == p for b in bands)
    assert sorted(special+[k for band in bands for k in band]) == list(range(1, (n+1)//2))
    reflection = next(h for h in range(g.n) if h not in rotations)
    columns = [energy((reflection,))]+[energy((t(k),)) for k in range(1, (n+1)//2)]
    kernel = nullspace(list(zip(*columns)), len(columns))
    assert len(kernel) == (m+1)//2
    basis = []
    for band in [None]+bands:
        row = [-p]+[0]*((n-1)//2)
        for k in special if band is None else band:
            row[k] = 2 if band is None else 1
        assert all(sum(a*b for a, b in zip(row, coordinate)) == 0 for coordinate in zip(*columns))
        basis.append(row)
    # The disjoint supports of special/band coordinates prove independence.
    return g, t, energy, product, {
        'prime': p, 'exponent': exponent, 'cycle_vertices': n,
        'automorphisms_from_adjacency': g.n,
        'special_rotation_classes': special, 'generic_rotation_bands': bands,
        'spectral_population_kernel_dimension': len(kernel),
        'integer_energy_kernel_basis_R_then_rotation_pm_k': basis,
        'ordered_product_parity': 'sum of all integer basis coefficients must be even',
    }


def close_exchange(g, source, target, energy, product):
    reverse = lambda w: tuple(g.inv[h] for h in w[::-1])
    pairs = set()
    for row in g.conj:
        a, b = tuple(row[h] for h in source), tuple(row[h] for h in target)
        pairs.add(tuple(sorted((a, b))))
        pairs.add(tuple(sorted((reverse(a), reverse(b)))))
    if len({w for pair in pairs for w in pair}) != 2*len(pairs):
        return None
    table = {a: b for x, y in pairs for a, b in ((x, y), (y, x))}
    for x, y in table.items():
        assert table[y] == x and product(x) == product(y) and energy(x) == energy(y)
        assert table[reverse(x)] == reverse(y)
        assert all(table[tuple(row[h] for h in x)] == tuple(row[h] for h in y) for row in g.conj)
    assert (g.identity,)*len(source) not in table
    return tuple(sorted(pairs))


def raw_triangle_witness(g, source, target, pairs, energy):
    oracle, patch = LinkOracle(4, g), 0
    paths = oracle.patches[patch][::-1]
    writes = sorted(oracle.writes[patch])
    outside_edge = min(oracle.reads[patch]-set(writes))
    word = lambda raw: tuple(oracle.transport(raw, path) for path in paths)
    initial = None
    for outer, a, b in itertools.product(range(g.n), repeat=3):
        raw = [g.identity]*len(oracle.edges)
        raw[outside_edge], raw[writes[0]], raw[writes[1]] = outer, a, b
        if word(raw) == source:
            initial = raw
            break
    assert initial is not None
    code = lambda w: (w[0]*g.n+w[1])*g.n+w[2]
    table = list(range(g.n**3))
    for a, b in pairs:
        table[code(a)], table[code(b)] = code(b), code(a)
    moved = initial[:]
    oracle.update(moved, patch, table)
    assert word(moved) == target
    assert all(initial[e] == moved[e] for e in range(len(initial)) if e not in writes)
    total_energy = lambda raw: energy(tuple(oracle.transport(raw, path) for path in oracle.face_paths))
    assert total_energy(initial) == total_energy(moved)
    replay = moved[:]
    oracle.update(replay, patch, table)
    assert replay == initial
    return {'base_mesh_side': 4, 'patch': patch, 'written_internal_edges': writes,
            'initial_links': initial, 'final_links': moved,
            'whole_mesh_spectral_energy_remainder': total_energy(initial),
            'exterior_links_fixed': True, 'exact_raw_inverse': True}


def c9_low_arity_census(g, energy, product):
    records = []
    for arity in (1, 2, 3):
        buckets = defaultdict(list)
        for w in itertools.product(range(g.n), repeat=arity):
            buckets[(product(w), energy(w))].append(w)
        closures = set()
        for words in buckets.values():
            for i, source in enumerate(words):
                for target in words[i+1:]:
                    if Counter(g.sectors[h] for h in source) == Counter(g.sectors[h] for h in target):
                        continue
                    pairs = close_exchange(g, source, target, energy, product)
                    if pairs is not None:
                        closures.add(pairs)
        energies = sorted({energy(pairs[0][0]) for pairs in closures})
        records.append({'arity': arity, 'all_tuples_examined': g.n**arity,
                        'reactive_minimal_closures': len(closures),
                        'exact_reactive_energy_remainders': energies})
    assert records[0]['reactive_minimal_closures'] == records[1]['reactive_minimal_closures'] == 0
    assert records[2]['reactive_minimal_closures'] > 0
    return records


def witness(p, exponent):
    g, t, energy, product, record = arithmetic(p, exponent)
    n, m = p**exponent, p**(exponent-1)
    s = ((p*p-1)//4) % p
    r = s*(m//p)
    source_exponents = [(r+j*m) % n for j in range(p)]
    target_exponents = [j*m for j in range(1, (p+1)//2) for _ in range(2)]+[0]
    source, target = tuple(map(t, source_exponents)), tuple(map(t, target_exponents))
    pairs = close_exchange(g, source, target, energy, product)
    assert pairs is not None and len(pairs) == 4
    assert energy(source) == (2*p,)+(0,)*((p-1)*m-1)
    assert product(source) == product(target) == t(s*m)
    record.update({'minimum_reactive_arity': p, 'minimum_reactive_spectral_sum': 2*p,
                   'source_rotation_exponents': source_exponents,
                   'target_rotation_exponents': target_exponents,
                   'source_tuple': source, 'target_tuple': target,
                   'transpositions': pairs, 'moved_tuples': 8})
    if n == 9:
        record['complete_low_arity_census'] = c9_low_arity_census(g, energy, product)
        record['raw_link_witness'] = raw_triangle_witness(g, source, target, pairs, energy)
    return record


if __name__ == '__main__':
    report = {'scope': 'Exact conditional spectral-conservation results, not derived physical energy or deployed reaction dynamics.',
              'prime_power_witnesses': [witness(p, a) for p, a in ((3, 2), (3, 3), (3, 4), (5, 2), (7, 2))]}
    output = Path('data/cycle-prime-power-reactions.json')
    output.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({'output': str(output), 'results': [
        {k: row[k] for k in ('cycle_vertices', 'minimum_reactive_arity', 'minimum_reactive_spectral_sum',
                            'spectral_population_kernel_dimension', 'source_rotation_exponents', 'target_rotation_exponents')}
        for row in report['prime_power_witnesses']]}), flush=True)
