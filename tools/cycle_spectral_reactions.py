#!/usr/bin/env python3
"""Sparse sharp-arity witnesses for a cycle-derived spectral observable.

The observable is a diagnostic, not a replacement energy in the running
model. No continuum gauge group or numerical energy tolerance is supplied.
"""
import json
from pathlib import Path

from cycle_charge_compatibility import cycle_group, pentagon_census
from triangle_reference import subgroup


def spectral_witness(p):
    assert p >= 3 and all(p % d for d in range(2, int(p**0.5)+1))
    g = cycle_group(p)
    t = lambda k: g.elements.index(tuple((i+k) % p for i in range(p)))
    r = lambda k: g.elements.index(tuple((k-i) % p for i in range(p)))
    rotations = {t(k): k for k in range(p)}
    def energy(word):
        # E(t_k)=2-zeta^k-zeta^-k, E(reflection)=2, E(identity)=0.
        # Exact reduction by Phi_p(z)=1+...+z^(p-1); no floating-point cosines.
        coefficients = [0]*p
        for h in word:
            if h == g.identity:
                continue
            coefficients[0] += 2
            if h in rotations:
                k = rotations[h]
                coefficients[k] -= 1
                coefficients[-k] -= 1
        return tuple(c-coefficients[-1] for c in coefficients[:-1])
    def product(word):
        value = g.identity
        for h in word:
            value = g.mul[value][h]
        return value
    reverse = lambda word: tuple(g.inv[h] for h in word[::-1])
    source = (r(0),)*(2*p-1)+(r(1),)*2
    target = (r(0),)+tuple(h for k in range(1, (p+1)//2) for h in (t(k), t(-k), t(k), t(-k)))+(g.identity,)*2
    assert len(source) == len(target) == 2*p+1
    pairs = set()
    for row in g.conj:
        a, b = tuple(row[h] for h in source), tuple(row[h] for h in target)
        pairs.add(tuple(sorted((a, b))))
        pairs.add(tuple(sorted((reverse(a), reverse(b)))))
    pairs = sorted(pairs)
    assert len(pairs) == 4*p and len({w for pair in pairs for w in pair}) == 8*p
    table = {a: b for x, y in pairs for a, b in ((x, y), (y, x))}
    assert (g.identity,)*(2*p+1) not in table
    for x, y in table.items():
        assert table[y] == x and product(x) == product(y) and energy(x) == energy(y)
        assert table[reverse(x)] == reverse(y)
        for row in g.conj:
            assert table[tuple(row[h] for h in x)] == tuple(row[h] for h in y)
    assert product(source) == product(target) == r(0)
    assert len(subgroup(g, source)) == len(subgroup(g, target)) == 2*p
    assert energy(source) == (4*p+2,)+(0,)*(p-2)
    return {'prime_cycle_vertices': p, 'graph_derived_automorphisms': g.n,
            'minimum_reactive_arity': 2*p+1, 'source_tuple': source, 'target_tuple': target,
            'class_population_change': {'reflections': -2*p, 'each_nontrivial_rotation_class': 4, 'identity': 2},
            'spectral_energy_of_each_tuple': 4*p+2, 'fixed_boundary_element': r(0),
            'conjugation_and_reversal_closed_transpositions': pairs, 'moved_tuples': 8*p,
            'scope': 'Exact algebraic multi-loop involution attaining the proved bound. Not an installed shared-mesh update or a selected physical Hamiltonian.'}


if __name__ == '__main__':
    census = pentagon_census()
    twice_cost = ((4, 0), (5, -1), (5, 1))
    residuals = []
    for family in census['charge_equation_families']:
        row = family['row']
        value = [sum(row[i]*twice_cost[i][j] for i in range(3)) for j in range(2)]
        assert any(value)
        residuals.append({'equation': row, 'twice_spectral_cost_change_rational_sqrt5': value})
    report = {'C5_spectral_costs_twice_in_Q_sqrt5': twice_cost,
              'C5_reactive_triple_family_cost_residuals': residuals,
              'sharp_arity_witnesses': [spectral_witness(p) for p in (3, 5, 7, 11, 17)]}
    output = Path('data/cycle-spectral-reactions.json')
    output.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({'output': str(output), 'witnesses': [
        {k: row[k] for k in ('prime_cycle_vertices', 'minimum_reactive_arity', 'moved_tuples', 'spectral_energy_of_each_tuple')}
        for row in report['sharp_arity_witnesses']]}), flush=True)
