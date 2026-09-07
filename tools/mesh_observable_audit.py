#!/usr/bin/env python3
"""Exact coarse-memory audit and a finite neighboring-loop invariant search."""
import argparse
import itertools
import json
import subprocess
from collections import Counter, defaultdict
from fractions import Fraction
from pathlib import Path
from screen_braid_rules import nullspace


def rational(value):
    return [value.numerator, value.denominator]


def verify_link_records(result, table, permutations):
    """Independent path multiplication, gauge normalization, and feature counts."""
    group = [tuple(g) for g in permutations]
    n = len(group)
    mul = [[group.index(tuple(a[b[i]] for i in range(4))) for b in group] for a in group]
    inv = [next(b for b in range(n) if mul[a][b] == 0) for a in range(n)]
    def conj(g, a): return mul[g][mul[a][inv[g]]]
    canonical_pair = [min(conj(g, a)*n+conj(g, b) for g in range(n)) for a in range(n) for b in range(n)]
    pair_reps = sorted(set(canonical_pair))
    pair_classes = [pair_reps.index(p) for p in canonical_pair]
    if pair_classes != result['pair_class_by_raw_pair'] or pair_reps != result['pair_class_representatives']:
        raise ValueError('based-pair gauge quotient disagrees with independent conjugation')
    def compiled(edges, specifications):
        keys = {(tuple(p['first']), tuple(p['second']), tuple(p['connector'])) for p in specifications}
        if len(keys) != len(specifications) or len(keys) != 8*len(edges):
            raise ValueError('closed triangular mesh lacks the complete ordered patch family')
        for first, second, connector in keys:
            if (second, first, tuple(reversed(connector))) not in keys:
                raise ValueError('patch family is not closed under exchanging the ordered faces')
            if (len(first) != 4 or len(second) != 4 or first[0] != first[-1] or second[0] != second[-1]
                    or len(set(first)) != 3 or len(set(second)) != 3
                    or connector[0] != second[0] or connector[-1] != first[0]
                    or len(set(connector)) != len(connector)):
                raise ValueError('invalid triangular patch or connector')
            left = {tuple(sorted(e)) for e in zip(first, first[1:])}
            right = {tuple(sorted(e)) for e in zip(second, second[1:])}
            writes = {tuple(sorted(first[-2:])), tuple(sorted(second[-2:]))}
            if (len(left & right) != 1 or writes & (left & right)
                    or any(tuple(sorted(e)) not in (left | right)-writes for e in zip(connector, connector[1:]))):
                raise ValueError('closing-link/connector support is invalid')
        def path(vertices):
            return [(edges.index(tuple(sorted((u, v)))), u > v) for u, v in zip(vertices, vertices[1:])]
        return [(path(p['first']), path(p['second']), path(p['connector'])) for p in specifications]
    def edge(values, link): return inv[values[link[0]]] if link[1] else values[link[0]]
    def transport(values, path):
        value = 0
        for link in path: value = mul[edge(values, link)][value]
        return value
    def feature(values, patches):
        counts = [0]*len(pair_reps)
        for first, second, connector in patches:
            a, b, t = transport(values, first), transport(values, second), transport(values, connector)
            counts[pair_classes[a*n+conj(t, b)]] += 1
        return counts
    def update(values, patch):
        first, second, connector = patch
        a, b, t = transport(values, first), transport(values, second), transport(values, connector)
        c, d = divmod(table[a*n+conj(t, b)], n)
        result = values[:]
        for path, old, new in ((first, a, c), (second, b, conj(inv[t], d))):
            link = path[-1]
            value = mul[mul[new][inv[old]]][edge(values, link)]
            result[link[0]] = inv[value] if link[1] else value
        return result
    tetra = result['tetrahedron']
    edges = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
    patches = compiled(edges, tetra['patch_specs'])
    if tetra['patches'] != len(patches):
        raise ValueError('tetrahedron generator count disagrees with patch specifications')
    face_paths = [[(edges.index(tuple(sorted((u, v)))), u > v) for u, v in zip(face, face[1:])]
                  for face in ((0, 1, 2, 0), (0, 2, 3, 0), (0, 3, 1, 0), (1, 3, 2, 1))]
    canonical = lambda state: min(tuple(conj(g, a) for a in state) for g in range(n))
    expected = {canonical(r) for r in itertools.product(range(n), repeat=3)}
    lookup = {tuple(s['rooted_chords']): s['id'] for s in tetra['states']}
    if set(lookup) != expected or len(lookup) != len(tetra['states']):
        raise ValueError('tetrahedron gauge quotient is incomplete')
    for index, state in enumerate(tetra['states']):
        if state['id'] != index or len(state['successors']) != len(patches):
            raise ValueError('tetrahedron state indexing or generator count is invalid')
        values = [0, 0, 0]+state['rooted_chords']
        expected_sectors = [min(conj(g, transport(values, path)) for g in range(n)) for path in face_paths]
        stabilizer = [g for g in range(n) if all(conj(g, a) == a for a in state['rooted_chords'])]
        if state['face_sectors'] != expected_sectors or state['stabilizer'] != stabilizer:
            raise ValueError('tetrahedron face sectors or stabilizer disagree with link algebra')
        if feature(values, patches) != state['pair_features']:
            raise ValueError('tetrahedron feature vector disagrees with independent link algebra')
        for p, target in enumerate(state['successors']):
            after = update(values, patches[p])
            rooted = [mul[inv[after[v-1]]][mul[after[index]][after[u-1]]] for index, (u, v) in enumerate(edges) if u]
            if lookup[canonical(rooted)] != target:
                raise ValueError('tetrahedron successor disagrees with independent link algebra')
    octa = result['octahedron']
    patches = compiled([tuple(e) for e in octa['edges']], octa['patches'])
    if len(octa['checks']) != octa['samples']*octa['steps_per_sample']:
        raise ValueError('octahedron sample coverage is incomplete')
    previous_after = None
    for index, check in enumerate(octa['checks']):
        sample, step = divmod(index, octa['steps_per_sample'])
        if (check['sample'], check['step']) != (sample, step):
            raise ValueError('octahedron sample indexing is invalid')
        before = check['before_links']
        if step and before != previous_after:
            raise ValueError('octahedron trajectory continuity failed')
        after = update(before, patches[check['patch']])
        if feature(before, patches) != check['before_features'] or feature(after, patches) != check['after_features']:
            raise ValueError('octahedron feature vector disagrees with independent link algebra')
        previous_after = after
    center = [a for a in range(n) if all(mul[a][g] == mul[g][a] for g in range(n))]
    frozen_gate = all((table[a*n+b] == a*n+b) if ((a in center) == (b in center)) else
                      ((table[a*n+b]//n in center) == (b in center)
                       and (table[a*n+b]%n in center) == (a in center))
                      for a in range(n) for b in range(n))
    return {'tetrahedron_successors': len(tetra['states'])*len(tetra['patch_specs']),
            'octahedron_feature_transitions': len(octa['checks']),
            'center_elements': center, 'central_type_frozen_gate_certificate': frozen_gate}


def memory_audit(local):
    raw = {tuple(row[:4]): tuple(row[4:]) for row in local['transitions']}
    n = len(local['element_sectors'])
    if (set(raw) != set(itertools.product(range(n), repeat=4))
            or len(raw) != len(local['transitions'])
            or any(raw.get(after) != before for before, after in raw.items())):
        raise ValueError('the microscopic two-step return is not an exhaustive involution')
    sectors = local['element_sectors']
    counts, mass = defaultdict(Counter), Counter()
    for before, after in raw.items():
        a, b = tuple(sectors[x] for x in before), tuple(sectors[x] for x in after)
        counts[a][b] += 1
        mass[a] += 1
    if any(count != counts[b][a] for a, row in counts.items() for b, count in row.items()):
        raise ValueError('coarse counts violate exact detailed balance of this involution')
    returns = {a: sum(Fraction(count, mass[a])*Fraction(counts[b][a], mass[b]) for b, count in row.items())
               for a, row in counts.items()}
    average = sum(Fraction(mass[a], len(raw))*p for a, p in returns.items())
    # The (current observation, next observation) refinement is closed for this
    # one involution: the next augmented observation is exactly the swapped pair.
    augmented = {(tuple(sectors[x] for x in a), tuple(sectors[x] for x in b)) for a, b in raw.items()}
    if any((b, a) not in augmented for a, b in augmented):
        raise ValueError('two-observation memory refinement failed')
    return {'microscopic_states': len(raw), 'coarse_states': len(counts),
            'exact_two_step_return_probability': [1, 1],
            'naive_markov_minimum_two_step_return': rational(min(returns.values())),
            'naive_markov_stationary_average_two_step_return': rational(average),
            'coarse_states_with_false_two_step_spreading': sum(p != 1 for p in returns.values()),
            'one_generator_memory_refinement_states': len(augmented),
            'measure': 'uniform gauge-fixed rooted holonomies, not uniform physical gauge orbits',
            'rows': [{'input_sectors': a, 'rooted_multiplicity': mass[a],
                      'outcomes': [{'output_sectors': b, 'count': count} for b, count in sorted(row.items())],
                      'markov_two_step_return': rational(returns[a])} for a, row in sorted(counts.items())]}


def predictive_partition(states):
    labels = [tuple(s['face_sectors']) for s in states]
    def compress(values):
        ids = {}
        return [ids.setdefault(value, len(ids)) for value in values]
    labels = compress(labels)
    counts = [len(set(labels))]
    while True:
        refined = compress([(labels[i], tuple(labels[j] for j in state['successors'])) for i, state in enumerate(states)])
        if len(set(refined)) == len(set(labels)):
            return {'refinement_sizes': counts, 'labels': refined}
        labels = refined
        counts.append(len(set(labels)))


def independent_observables(weights, feature_groups):
    """Remove coefficient null directions and constants separately on each mesh."""
    pivots, selected = {}, []
    for weight in weights:
        values = [sum(a*b for a, b in zip(weight, row)) for features in feature_groups for row in features]
        centered = []
        cursor = 0
        for features in feature_groups:
            centered.extend(Fraction(v-values[cursor]) for v in values[cursor:cursor+len(features)])
            cursor += len(features)
        for pivot, basis in sorted(pivots.items()):
            scale = centered[pivot]
            if scale:
                centered = [a-scale*b for a, b in zip(centered, basis)]
        pivot = next((i for i, value in enumerate(centered) if value), None)
        if pivot is not None:
            scale = centered[pivot]
            pivots[pivot] = [v/scale for v in centered]
            selected.append({'weights': weight, 'values': values})
    return selected


def analyze_orbits(result):
    states = result['tetrahedron']['states']
    size = len(states)
    patches = result['tetrahedron']['patches']
    for p in range(patches):
        if sorted(s['successors'][p] for s in states) != list(range(size)):
            raise ValueError('generator is not a permutation of physical gauge states')
        for i, s in enumerate(states):
            j = s['successors'][p]
            if states[j]['successors'][p] != i or len(states[j]['stabilizer']) != len(s['stabilizer']):
                raise ValueError('involution or gauge-orbit stabilizer size changed')
    if sum(8//len(s['stabilizer']) for s in states) != result['tetrahedron']['rooted_states']:
        raise ValueError('physical orbit multiplicities fail to recover the rooted state count')
    central = result['independent_verification']['center_elements']
    frozen_counts = Counter()
    for s in states:
        predicted_frozen = len({x in central for x in s['face_sectors']}) == 1
        fixed = all(i == s['id'] for i in s['successors'])
        if predicted_frozen != fixed:
            raise ValueError('central-type frozen criterion fails on tetrahedron')
        if fixed: frozen_counts['all_central' if s['face_sectors'][0] in central else 'all_noncentral'] += 1
    remaining, components = set(range(size)), []
    while remaining:
        stack = [min(remaining)]
        remaining.remove(stack[0])
        for i in stack:
            for j in states[i]['successors']:
                if j in remaining: remaining.remove(j); stack.append(j)
        components.append(sorted(stack))
    weights = sorted({len(states[i]['stabilizer']) for component in components for i in component})
    cycle_counts = {}
    for name, order in (('index_order', list(range(patches))),
                        ('permuted_order', sorted(range(patches), key=lambda p: (17*p) % patches))):
        successor = []
        for i in range(size):
            for p in order: i = states[i]['successors'][p]
            successor.append(i)
        unseen, cycles = set(range(size)), Counter()
        while unseen:
            i, length = min(unseen), 0
            while i in unseen: unseen.remove(i); length += 1; i = successor[i]
            cycles[length] += 1
        cycle_counts[name] = {'order': order, 'cycle_counts': [[length, count] for length, count in sorted(cycles.items())]}
    features = [s['pair_features'] for s in states]
    dimension = len(result['pair_class_representatives'])
    equations = sorted({tuple(b-a for a, b in zip(s['pair_features'], states[j]['pair_features']))
                        for s in states for j in s['successors']})
    basis = nullspace(equations, dimension)
    observables = independent_observables(basis, [features])
    oct_checks = result['octahedron']['checks']
    oct_equations = {tuple(b-a for a, b in zip(row['before_features'], row['after_features'])) for row in oct_checks}
    joint_basis = nullspace(sorted(set(equations) | oct_equations), dimension)
    oct_features = [row[key] for row in oct_checks for key in ('before_features', 'after_features')]
    survivors = independent_observables(joint_basis, [features, oct_features])
    # Reversing an ordered patch exchanges its two commonly based holonomies
    # up to simultaneous conjugation. Antisymmetric coefficients sum to zero
    # on EVERY mesh with the complete ordered operator family, not just samples.
    pair_class = result['pair_class_by_raw_pair']
    reverse_class = [pair_class[(p % 8)*8+p//8] for p in result['pair_class_representatives']]
    trivial_certificate = all(len({Fraction(w[i]+w[reverse_class[i]], 2) for i in range(dimension)}) == 1
                              for w in joint_basis)
    if trivial_certificate and survivors:
        raise ValueError('globally constant/antisymmetric weights produced a nonconstant observable')
    for observable in observables:
        weight = observable['weights']
        observable['first_octahedron_violation'] = next((i for i, row in enumerate(oct_checks)
            if sum(w*(b-a) for w, a, b in zip(weight, row['before_features'], row['after_features']))), None)
        if observable['first_octahedron_violation'] is not None:
            row = oct_checks[observable['first_octahedron_violation']]
            observable['violation_values'] = [sum(a*b for a, b in zip(weight, row[key]))
                                              for key in ('before_features', 'after_features')]
    result['analysis'] = {
        'physical_states': size, 'generator_count': patches, 'stabilizer_sizes': weights,
        'frozen_gauge_states': dict(frozen_counts),
        'components': components, 'component_size_counts': sorted(Counter(map(len, components)).items()),
        'predictive_face_partition': predictive_partition(states), 'periodic_sweeps': cycle_counts,
        'pair_ansatz': 'sum over all ordered adjacent-face closing-link patches of one common function of the based-pair gauge orbit',
        'coefficient_dimension': dimension, 'tetrahedron_weight_nullity': len(basis),
        'tetrahedron_nonconstant_observables': observables,
        'joint_weight_nullity': len(joint_basis), 'joint_nonconstant_observables': survivors,
        'reverse_pair_classes': reverse_class,
        'global_triviality_certificate': {'constant_symmetric_part_for_every_joint_weight': trivial_certificate,
                                          'antisymmetric_dimension': sum(i != reverse_class[i] for i in range(dimension))//2},
        'tetrahedron_equations': equations, 'joint_weights': joint_basis,
        'scope': 'tetrahedron exhaustive; octahedron sampled falsification only; not a classification of all local energies'}
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--binary', type=Path, default=Path('build/wgphysics_mesh_orbits'))
    parser.add_argument('--output', type=Path, default=Path('out/mesh-observable-audit.json'))
    parser.add_argument('--reanalyze', type=Path)
    args = parser.parse_args()
    mesh = json.loads(Path('data/d4-shared-mesh.json').read_text())
    census = json.loads(Path('data/d4-involution-braid-search.json').read_text())
    table = census['solutions'][mesh['solution_id']]['table']
    if args.reanalyze:
        result = json.loads(args.reanalyze.read_text())
        if result['solution_id'] != mesh['solution_id']: raise ValueError('rule mismatch')
    else:
        completed = subprocess.run([str(args.binary)], input=' '.join(map(str, table))+'\n', text=True,
                                   stdout=subprocess.PIPE, check=True, timeout=1800)
        result = json.loads(completed.stdout)
    result['solution_id'] = mesh['solution_id']
    result['independent_verification'] = verify_link_records(result, table, mesh['group_permutations'])
    result['memory'] = memory_audit(mesh['local_census'])
    analyze_orbits(result)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True)+'\n')
    print('memory return:', result['memory']['naive_markov_stationary_average_two_step_return'])
    print('components:', result['analysis']['component_size_counts'])
    print('predictive refinement:', result['analysis']['predictive_face_partition']['refinement_sizes'])
    print('nonconstant pair invariants:', len(result['analysis']['tetrahedron_nonconstant_observables']),
          'after larger-mesh checks:', len(result['analysis']['joint_nonconstant_observables']))


if __name__ == '__main__':
    main()
