#!/usr/bin/env python3
"""Charge compatibility under graph-derived cycle-fiber refinement.

No rule is added to the evolution engine. This is a census of the existing
minimal-closure selection criterion and its common conserved-charge space.
"""
import itertools
import json
import math
from collections import Counter
from pathlib import Path

from face_energy_obstruction import FiniteGroup
from screen_braid_rules import nullspace
from triple_channels import positive_certificate, primitive_row
from triple_rule_search import algebra, assess, make_table, minimal_closures
from run_three_face import LinkOracle


def cycle_group(n):
    adjacency = [{(i-1) % n, (i+1) % n} for i in range(n)]
    found = []
    def extend(path):
        if len(path) == n:
            if path[0] in adjacency[path[-1]]:
                found.append(tuple(path))
            return
        for target in sorted(adjacency[path[-1]]-set(path)):
            extend(path+[target])
    for first in range(n):
        extend([first])
    assert all({p[j] for j in adjacency[i]} == adjacency[p[i]] for p in found for i in range(n))
    return FiniteGroup(sorted(found))


def explicit_obstruction(n):
    g = cycle_group(n)
    states, products, actions, reverse = algebra(g)
    operations = actions+[[reverse[i] for i in row] for row in actions]
    labels = sorted(set(g.sectors)-{g.identity})
    rotation = lambda k: g.elements.index(tuple((i+k) % n for i in range(n)))
    reflection = lambda k: g.elements.index(tuple((k-i) % n for i in range(n)))
    kinds = {g.sectors[reflection(0)]: 'reflection'}
    kinds.update({g.sectors[rotation(k)]: f'rotation_pm_{k}' for k in range(1, (n+1)//2)})
    records, selected_rows = [], []
    examples = [(f'reflection_conversion_{k}', (reflection(0), reflection(0), reflection(k)),
                 (reflection(0), rotation(-k), g.identity)) for k in range(1, (n+1)//2)]
    examples.append(('rotation_fusion_1_1_to_2', (rotation(1), rotation(1), g.identity),
                     (rotation(2), g.identity, g.identity)))
    for name, source, target in examples:
        x, y = states.index(source), states.index(target)
        assert products[x] == products[y]
        pairs = tuple(sorted({tuple(sorted((op[x], op[y]))) for op in operations}))
        assert len({z for p in pairs for z in p}) == 2*len(pairs)
        table = make_table(pairs, g.n**3)
        audit = assess(g, table)
        equations = [tuple(sum(g.sectors[h] == s for h in states[a])-sum(g.sectors[h] == s for h in states[b])
                           for s in labels) for a, b in pairs]
        positive = positive_certificate(equations, len(labels))
        if positive['exists']:
            selected_rows.extend(equations)
        records.append({'name': name, 'input_triple': source, 'output_triple': target,
                        'transpositions': pairs, 'moved_states': audit['moved_states'],
                        'positive_charge_certificate': positive})
    raw_witnesses = [raw_link_witness(g, r) for r in records] if n == 5 else []
    return {'cycle_vertices': n, 'automorphisms_from_adjacency': g.n,
            'class_coordinate_order': [kinds[s] for s in labels], 'explicit_rules': records,
            'raw_link_witnesses': raw_witnesses,
            'common_charge_basis_for_individually_positive_examples': nullspace(selected_rows, len(labels))}


def raw_link_witness(g, record):
    oracle, patch = LinkOracle(4, g), 0
    paths = oracle.patches[patch][::-1]
    writes = sorted(oracle.writes[patch])
    outside_edge = min(oracle.reads[patch]-set(writes))
    source, target = tuple(record['input_triple']), tuple(record['output_triple'])
    word = lambda raw: tuple(oracle.transport(raw, path) for path in paths)
    initial = None
    for outer, a, b in itertools.product(range(g.n), repeat=3):
        raw = [g.identity]*len(oracle.edges)
        raw[outside_edge], raw[writes[0]], raw[writes[1]] = outer, a, b
        if word(raw) == source:
            initial = raw
            break
    assert initial is not None
    table = make_table(record['transpositions'], g.n**3)
    moved = initial[:]
    oracle.update(moved, patch, table)
    assert word(moved) == target and all(initial[e] == moved[e] for e in range(len(initial)) if e not in writes)
    n = len(g.elements[0])
    reflection = g.elements.index(tuple(-i % n for i in range(n)))
    coarse = [0 if h == g.identity else 1 if g.sectors[h] == g.sectors[reflection] else 2 for h in range(g.n)]
    charge = lambda raw: sum(coarse[oracle.transport(raw, p)] for p in oracle.face_paths)
    reversed_raw = moved[:]
    oracle.update(reversed_raw, patch, table)
    assert reversed_raw == initial
    return {'rule': record['name'], 'base_mesh_side': 4, 'patch': patch, 'internal_written_edges': writes,
            'initial_links': initial, 'final_links': moved, 'based_input': source, 'based_output': target,
            'whole_mesh_coarse_charge_change': charge(moved)-charge(initial),
            'all_exterior_links_unchanged': True, 'raw_inverse_returns_exactly': True}


def pentagon_census():
    g = cycle_group(5)
    states, _, _, _ = algebra(g)
    labels = sorted(set(g.sectors)-{g.identity})
    populations = [tuple(sum(g.sectors[a] == s for a in state) for s in labels) for state in states]
    closures = minimal_closures(g)
    equations, families = [], Counter()
    for pairs in closures:
        rows = [tuple(a-b for a, b in zip(populations[x], populations[y])) for x, y in pairs]
        if not any(sum(row) for row in rows):
            continue
        certificate = positive_certificate(rows, len(labels))
        if certificate['exists']:
            row = primitive_row(certificate['equation'])
            families[row] += 1
            equations.append(row)
    rays = set()
    rows = sorted(families)
    for a, b in itertools.combinations(rows, 2):
        v = (a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0])
        if all(x > 0 for x in v) or all(x < 0 for x in v):
            sign = 1 if v[0] > 0 else -1
            divisor = math.gcd(*v)
            rays.add(tuple(sign*x//divisor for x in v))
    compatible = []
    for ray in sorted(rays):
        active = [r for r in rows if sum(a*b for a, b in zip(ray, r)) == 0]
        compatible.append({'primitive_positive_charge': ray, 'conserved_equation_families': active,
                           'compatible_minimal_rules': sum(families[r] for r in active)})
    uncovered = [r for r in rows if not any(sum(a*b for a, b in zip(ray, r)) == 0 for ray in rays)]
    return {'cycle_vertices': 5, 'derived_group_order': g.n,
            'class_representatives': labels,
            'class_representative_permutations': [g.elements[i] for i in labels],
            'minimal_closures': len(closures), 'individually_positive_occupied_changing_rules': sum(families.values()),
            'charge_equation_families': [{'row': row, 'rules': families[row]} for row in rows],
            'joint_charge_basis': nullspace(equations, len(labels)),
            'maximal_jointly_positive_rays': compatible,
            'positive_equation_planes_without_an_intersection_ray': uncovered,
            'scope': 'Complete minimal transposition-closure census at C5, not all composite updates or a physical selection of a charge ray.'}


if __name__ == '__main__':
    result = {'odd_cycle_certificates': [explicit_obstruction(n) for n in (3, 5, 7, 9)],
              'pentagon_census': pentagon_census()}
    output = Path('data/cycle-charge-compatibility.json')
    output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({'output': str(output), 'odd_cycle_joint_dimensions': [
        [r['cycle_vertices'], len(r['common_charge_basis_for_individually_positive_examples'])]
        for r in result['odd_cycle_certificates']], 'pentagon_census': result['pentagon_census']}), flush=True)
