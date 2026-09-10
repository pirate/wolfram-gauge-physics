#!/usr/bin/env python3
"""Exact word algebra, spectral Fourier coefficients, and lookup-rule witnesses.

No dynamics or energy filter is changed. Integer frequencies describe the
adjacency-derived two-dimensional fiber action, not quantum amplitudes.
"""
import itertools
import json
from collections import Counter
from pathlib import Path

from cycle_charge_compatibility import cycle_group
from screen_braid_rules import nullspace
from triangle_reference import subgroup


def reduce_word(word):
    result = []
    for letter in word:
        if result and result[-1] == -letter:
            result.pop()
        else:
            result.append(letter)
    return tuple(result)


def substitute(word, replacements):
    return reduce_word([letter for index in word for letter in
                        (replacements[index-1] if index > 0 else
                         tuple(-x for x in replacements[-index-1][::-1]))])


def affine_word(word, signs):
    """Composition in x -> sign*x + angle, with formal input angles."""
    parity, coefficients = 1, [0]*len(signs)
    for letter in word:
        i = abs(letter)-1
        coefficients[i] += parity*(1 if letter > 0 else -signs[i])
        parity *= signs[i]
    return parity, tuple(coefficients)


def spectral_fourier(words, signs):
    # Twice (sum cos(output rotations) - sum cos(input rotations)).
    frequencies = Counter()
    for word in words:
        parity, coefficients = affine_word(word, signs)
        if parity == 1:
            frequencies[coefficients] += 1
            frequencies[tuple(-x for x in coefficients)] += 1
    for i, sign in enumerate(signs):
        if sign == 1:
            basis = tuple(int(j == i) for j in range(len(signs)))
            frequencies[basis] -= 1
            frequencies[tuple(-x for x in basis)] -= 1
    return {key: value for key, value in sorted(frequencies.items()) if value}


def evaluate(g, words, values):
    def evaluate_word(word):
        value = g.identity
        for letter in word:
            h = values[abs(letter)-1]
            value = g.mul[value][h if letter > 0 else g.inv[h]]
        return value
    return tuple(evaluate_word(word) for word in words)


def finite_spectral_fourier(words, signs, n):
    frequencies = Counter()
    for vector, coefficient in spectral_fourier(words, signs).items():
        frequencies[tuple(k % n for k in vector)] += coefficient
    return {key: value for key, value in frequencies.items() if value}


def primitive_report(name, words, inverse_words):
    k = len(words)
    assert reduce_word(sum(words, ())) == tuple(range(1, k+1))
    assert tuple(substitute(w, words) for w in inverse_words) == tuple((i,) for i in range(1, k+1))
    assert tuple(substitute(w, inverse_words) for w in words) == tuple((i,) for i in range(1, k+1))
    patterns = []
    for signs in itertools.product((1, -1), repeat=k):
        frequencies = spectral_fourier(words, signs)
        patterns.append({'input_signs': signs,
                         'output_affine_forms': [affine_word(w, signs) for w in words],
                         'spectral_residual_frequencies': [[v, c] for v, c in frequencies.items()]})
    charges = []
    for n in (3, 4, 5, 9):
        g = cycle_group(n)
        labels = sorted(set(g.sectors)-{g.identity})
        equations = set()
        for values in itertools.product(range(g.n), repeat=k):
            target = evaluate(g, words, values)
            equations.add(tuple(sum(g.sectors[h] == s for h in values)-sum(g.sectors[h] == s for h in target)
                                for s in labels))
        basis = nullspace(sorted(equations), len(labels))
        finite_spectral = all(not finite_spectral_fourier(words, signs, n)
                              for signs in itertools.product((1, -1), repeat=k))
        charges.append({'cycle_vertices': n, 'nonidentity_class_representatives': labels,
                        'conserved_additive_class_basis': basis,
                        'finite_spectral_conservation': finite_spectral})
    return {'name': name, 'words': words, 'inverse_words': inverse_words,
            'all_spectral_residuals_zero': all(not row['spectral_residual_frequencies'] for row in patterns),
            'sign_pattern_fourier_calculation': patterns, 'finite_additive_charge_spaces': charges}


def triangle_quotient_witnesses():
    bank = json.loads(Path('data/triangle-feedback.json').read_text())['bank']
    g = cycle_group(3)
    assert [list(p) for p in g.elements] == bank['group_automorphisms']
    signs = [(-1)**sum(p[i] > p[j] for i in range(3) for j in range(i+1, 3)) for p in g.elements]
    states = list(itertools.product(range(g.n), repeat=3))
    reflections = [h for h in range(g.n) if signs[h] == -1]
    parallel = (reflections[0],)*3
    mixed = (reflections[0], reflections[0], reflections[1])
    a, b = states.index(parallel), states.index(mixed)
    assert all(table[a] == a and table[b] != b for table in bank['tables'])
    records = []
    for rule_id, table in zip(bank['selected_rule_ids'], bank['tables']):
        target_a, target_b = states[table[a]], states[table[b]]
        sigma = lambda w: tuple(signs[h] for h in w)
        assert sigma(parallel) == sigma(mixed) and sigma(target_a) != sigma(target_b)
        records.append({'rule_id': rule_id, 'same_quotient_input': sigma(parallel),
                        'parallel_input': parallel, 'mixed_input': mixed,
                        'parallel_output': target_a, 'mixed_output': target_b,
                        'parallel_quotient_output': sigma(target_a), 'mixed_quotient_output': sigma(target_b)})
    return {'quotient': 'sign: Aut(C3) -> C2', 'every_selected_rule_fails_quotient_descent': records,
            'uniform_bank_parallel_output_counts': [
                [list(k), v] for k, v in sorted(Counter(tuple(signs[h] for h in states[t[a]]) for t in bank['tables']).items())],
            'uniform_bank_mixed_output_counts': [
                [list(k), v] for k, v in sorted(Counter(tuple(signs[h] for h in states[t[b]]) for t in bank['tables']).items())]}


def prime_power_subgroup_witnesses():
    report = json.loads(Path('data/cycle-prime-power-reactions.json').read_text())
    records = []
    for row in report['prime_power_witnesses']:
        g = cycle_group(row['cycle_vertices'])
        source, target = tuple(row['source_tuple']), tuple(row['target_tuple'])
        before, after = subgroup(g, source), subgroup(g, target)
        assert len(before) == row['prime']**2 and len(after) == row['prime']
        assert set(after) < set(before)
        records.append({'cycle_vertices': row['cycle_vertices'], 'source_generated_subgroup_order': len(before),
                        'target_generated_subgroup_order': len(after),
                        'source_tuple': source, 'target_tuple': target})
    return records


def relational_gate(n):
    """Define a matching predicate first; solve conservation only afterward."""
    g = cycle_group(n)
    involutions = [h for h in range(g.n) if h != g.identity and g.inv[h] == h]
    def p(w):
        a, b, c = w
        return g.mul[a][b], g.inv[b], g.mul[b][c]
    def matching(w):
        a, b, c = w
        return all(h in involutions for h in w) and ((a == b) != (b == c))
    sources = {w for a, b in itertools.permutations(involutions, 2) for w in ((a, a, b), (b, a, a))}
    table = {}
    for source in sources:
        target = p(source)
        assert matching(source) and not matching(target) and p(target) == source
        table[source], table[target] = target, source
    assert len(table) == 4*n*(n-1)
    reverse = lambda w: tuple(g.inv[h] for h in w[::-1])
    labels = sorted(set(g.sectors)-{g.identity})
    equations = set()
    for source, target in table.items():
        assert table[target] == source and table[reverse(source)] == reverse(target)
        assert g.mul[g.mul[source[0]][source[1]]][source[2]] == g.mul[g.mul[target[0]][target[1]]][target[2]]
        for row in g.conj:
            assert table[tuple(row[h] for h in source)] == tuple(row[h] for h in target)
        equations.add(tuple(sum(g.sectors[h] == s for h in source)-sum(g.sectors[h] == s for h in target) for s in labels))
    basis = nullspace(sorted(equations), len(labels))
    expected = [1 if h in involutions else 2 for h in labels]
    assert basis == [expected]
    report = {'cycle_vertices': n, 'matching_source_count': len(sources), 'moved_tuples': len(table),
              'nonidentity_class_representatives': labels, 'derived_additive_charge_basis': basis}
    if n == 3:
        states = list(itertools.product(range(g.n), repeat=3))
        ids = {w: i for i, w in enumerate(states)}
        full = [ids[table.get(w, w)] for w in states]
        bank = json.loads(Path('data/triangle-feedback.json').read_text())['bank']
        matches = [rule for rule, entries in zip(bank['selected_rule_ids'], bank['tables']) if entries == full]
        assert len(matches) == 1
        report['exact_existing_triangle_rule_id'] = matches[0]
    return g, table, report


def relational_gate_family():
    instances = {n: relational_gate(n) for n in (3, 5, 9, 27)}
    embeddings = []
    for source_n, target_n in ((3, 9), (3, 27), (9, 27)):
        a, ta, _ = instances[source_n]
        b, tb, _ = instances[target_n]
        ratio = target_n//source_n
        mapping = []
        for perm in a.elements:
            shift = perm[0]
            direction = 1 if (perm[1]-perm[0]) % source_n == 1 else -1
            lifted = tuple((direction*j+ratio*shift) % target_n for j in range(target_n))
            mapping.append(b.elements.index(lifted))
        assert len(set(mapping)) == a.n
        assert all(mapping[a.mul[x][y]] == b.mul[mapping[x]][mapping[y]] for x in range(a.n) for y in range(a.n))
        for word in itertools.product(range(a.n), repeat=3):
            lifted = tuple(mapping[h] for h in word)
            assert tuple(mapping[h] for h in ta.get(word, word)) == tb.get(lifted, lifted)
        embeddings.append({'source_cycle': source_n, 'target_cycle': target_n,
                           'all_source_tuples': a.n**3, 'intertwines_exactly': True})
    return {'predicate': 'All three inputs are nonidentity involutions, with exactly one equal adjacent pair; also enable the inverse image of this set under P.',
            'word': '(A B, inverse(B), B C)', 'fibers': [instance[2] for instance in instances.values()],
            'injective_refinement_examples': embeddings}


if __name__ == '__main__':
    primitives = [primitive_report('Hurwitz', ((1, 2, -1), (1,)), ((2,), (-2, 1, 2))),
                  primitive_report('BoundaryShear', ((1, 2, 1), (-1,)), ((-2,), (2, 1, 2))),
                  primitive_report('PartialProductSwap', ((1, 2), (-2,), (2, 3)), ((1, 2), (-2,), (2, 3)))]
    assert [p['all_spectral_residuals_zero'] for p in primitives] == [True, False, False]
    assert all(not row['conserved_additive_class_basis'] for row in primitives[-1]['finite_additive_charge_spaces'])
    result = {'scope': 'Word-rule conservation and quotient obstructions, not a new energy law or deployed evolution.',
              'primitive_word_calculations': primitives, 'triangle_lookup_rule_quotient_witnesses': triangle_quotient_witnesses(),
              'prime_power_lookup_rule_subgroup_witnesses': prime_power_subgroup_witnesses(),
              'relationally_gated_word_family': relational_gate_family()}
    output = Path('data/fiber-word-conservation.json')
    output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({'output': str(output), 'primitive_summary': [
        {'name': p['name'], 'spectral_conserving': p['all_spectral_residuals_zero'],
         'additive_charge_dimensions': [len(c['conserved_additive_class_basis']) for c in p['finite_additive_charge_spaces']]}
        for p in primitives],
        'triangle_nonword_rules': len(result['triangle_lookup_rule_quotient_witnesses']['every_selected_rule_fails_quotient_descent']),
        'prime_power_subgroup_orders': [[r['cycle_vertices'], r['source_generated_subgroup_order'], r['target_generated_subgroup_order']]
                                       for r in result['prime_power_lookup_rule_subgroup_witnesses']],
        'relational_gate_fibers': result['relationally_gated_word_family']['fibers']}), flush=True)
