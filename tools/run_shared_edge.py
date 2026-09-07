#!/usr/bin/env python3
"""Derive mesh charges for the one-shared-edge lift; optionally replay every event."""
import argparse
import itertools
import json
import subprocess
from pathlib import Path
from charge_quotient import induced_table, classify_block_exchange
from face_energy_obstruction import derive_fiber_group
from screen_braid_rules import nullspace


def mesh_geometry(side):
    vertex = lambda x, y: (y % side)*side+x % side
    edges, faces = set(), []
    for y in range(side):
        for x in range(side):
            a, b, c, d = vertex(x, y), vertex(x+1, y), vertex(x+1, y+1), vertex(x, y+1)
            for face in ((a, b, c), (a, c, d)):
                face = min(face[k:]+face[:k] for k in range(3))
                faces.append(face+(face[0],))
                edges.update(tuple(sorted(e)) for e in zip(faces[-1], faces[-1][1:]))
    edges = sorted(edges)
    incidence = [[f for f, face in enumerate(faces) if edge in
                  [tuple(sorted(e)) for e in zip(face, face[1:])]] for edge in edges]
    # Match the documented sorted-edge, face-insertion, two-root ordering.
    patches = []
    for edge, adjacent in zip(edges, incidence):
        if len(adjacent) != 2:
            raise ValueError('expected a closed two-manifold')
        for first in adjacent:
            second = next(f for f in adjacent if f != first)
            face = faces[first]
            u, v = next(e for e in zip(face, face[1:]) if tuple(sorted(e)) == edge)
            def rotate(f):
                loop = faces[f][:-1]
                k = loop.index(u)
                loop = loop[k:]+loop[:k]
                return loop+(u,)
            a, b = rotate(first), rotate(second)
            if b[-2] != v:
                raise ValueError('shared face orientations do not cancel')
            patches.append((a, b, first, second))
    return edges, faces, patches


def class_memory_model(group, table, charges):
    """Derive, then exhaustively verify, a two-state vacancy-memory factor."""
    labels = sorted(set(group.sectors))
    zero = [s for s in labels if all(q[s] == 0 for q in charges)]
    if len(zero) != 2 or group.identity not in zero:
        return None
    zero = [group.identity, next(s for s in zero if s != group.identity)]
    factor = induced_table(table, group.sectors)
    if factor is None:
        return None
    lookup = dict(zip(itertools.product(labels, repeat=2), map(tuple, factor)))
    flips = {}
    for a, b in itertools.product(labels, repeat=2):
        x, y = lookup[a, b]
        if (a in zero) == (b in zero):
            if (x, y) != (a, b): return None
        else:
            occupied = b if a in zero else a
            old_empty = a if a in zero else b
            new_empty = y if a in zero else x
            if (x if a in zero else y) != occupied or new_empty not in zero:
                return None
            bit = zero.index(old_empty) ^ zero.index(new_empty)
            if occupied in flips and flips[occupied] != bit:
                return None
            flips[occupied] = bit
    flip_species = sorted(s for s, bit in flips.items() if bit)
    for a, b in itertools.product(labels, repeat=2):
        x, y = lookup[a, b]
        for first_part in (0, 1):
            value = lambda u, v: ((u == zero[1])+(v == zero[1])
                                  +first_part*(u in flip_species)+(1-first_part)*(v in flip_species)) % 2
            if value(a, b) != value(x, y):
                raise ValueError('derived staggered parity fails a local transition')
    return {'sector_labels': labels, 'sector_table': factor, 'vacancy_sectors': zero,
            'vacancy_flip_by_species': [[s, bit] for s, bit in sorted(flips.items())],
            'flip_species': flip_species, 'local_staggered_parity_checks': 2*len(labels)**2}


def bipartition(size, face_pairs):
    adjacent = [set() for _ in range(size)]
    for a, b in face_pairs:
        adjacent[a].add(b); adjacent[b].add(a)
    colors = [None]*size
    for start in range(size):
        if colors[start] is not None: continue
        colors[start] = 0
        queue = [start]
        for a in queue:
            for b in sorted(adjacent[a]):
                if colors[b] is None: colors[b] = 1-colors[a]; queue.append(b)
                elif colors[b] == colors[a]: raise ValueError('dual graph is not bipartite')
    return colors


def replay(result, table, group, memory):
    edges, faces, specs = mesh_geometry(result['side'])
    n, mul, inv = group.n, group.mul, group.inv
    def compile_path(loop):
        return [(edges.index(tuple(sorted(e))), e[0] > e[1]) for e in zip(loop, loop[1:])]
    face_paths = [compile_path(f) for f in faces]
    patches = [(compile_path(a), compile_path(b), f, g) for a, b, f, g in specs]
    colors = bipartition(len(faces), [p[2:] for p in patches]) if memory is not None else None
    if len(edges) != result['edges'] or len(faces) != result['faces'] or len(patches) != result['patches']:
        raise ValueError('exported mesh dimensions disagree with independent construction')
    def transport(values, path):
        h = group.identity
        for edge, reverse in path:
            h = mul[inv[values[edge]] if reverse else values[edge]][h]
        return h
    count = 0
    for run in result['runs']:
        hops, flips, first_parity = 0, 0, None
        schedule = run['schedule']
        if (len(schedule) != run['schedule_layers'] or
                sorted(p for layer in schedule for p in layer) != list(range(len(patches)))):
            raise ValueError('schedule drops or duplicates operators')
        for layer in schedule:
            used = set()
            for p in layer:
                pair = set(patches[p][2:])
                if used & pair:
                    raise ValueError('shared-edge layer contains a face conflict')
                used |= pair
        values = run['initial_links'][:]
        if len(values) != len(edges) or any(not 0 <= v < n for v in values):
            raise ValueError('invalid initial link vector')
        for tick, row in enumerate(run['trajectory']):
            if tick:
                for p in schedule[(tick-1) % len(schedule)]:
                    first, second, _, _ = patches[p]
                    a, b = transport(values, first), transport(values, second)
                    if memory is not None:
                        sa, sb = group.sectors[a], group.sectors[b]
                        zero = memory['vacancy_sectors']
                        if (sa in zero) != (sb in zero):
                            hops += 1
                            flips += (sb if sa in zero else sa) in memory['flip_species']
                    x, y = divmod(table[a*n+b], n)
                    # Exterior-path oracle, unlike the compiled S A^-1 X kernel.
                    q = transport(values, first[1:])
                    new = mul[inv[q]][x]
                    edge, reverse = first[0]
                    values[edge] = inv[new] if reverse else new
                    if transport(values, first) != x or transport(values, second) != y:
                        raise ValueError('independent shared-edge targets failed')
                    count += 1
            histogram = [0]*n
            face_sectors = [group.sectors[transport(values, face)] for face in face_paths]
            for sector in face_sectors: histogram[sector] += 1
            if memory is not None:
                parity = (histogram[memory['vacancy_sectors'][1]]
                          +sum(color for color, sector in zip(colors, face_sectors) if sector in memory['flip_species'])) % 2
                if tick == 0: first_parity = parity
                if parity != first_parity:
                    raise ValueError('derived staggered parity drifted in link-level replay')
            if histogram != row[5:]:
                raise ValueError('independent layer histogram differs from export')
        if values != run['final_links']:
            raise ValueError('independent final links differ from export')
        run['independent_transport'] = {'charge_hops': hops, 'vacancy_flip_hops': flips,
                                        'conserved_staggered_parity': first_parity} if memory is not None else None
    if count != sum(run['updates'] for run in result['runs']):
        raise ValueError('update count disagrees with independent replay')
    return count


def analyze(result, table, full_replay=False):
    if result.get('lift') != 'shared-edge':
        raise ValueError('this analyzer requires the explicitly selected shared-edge lift')
    group = derive_fiber_group()
    group.validate(table, strict=True)
    if [list(g) for g in group.elements] != result['group_permutations']:
        raise ValueError('exported fiber group differs from square adjacency automorphisms')
    n, sectors = group.n, group.sectors
    labels = sorted(set(sectors)-{group.identity})
    rows = result['local_census']['transitions']
    if len(rows) != n**4 or {tuple(r[:4]) for r in rows} != set(itertools.product(range(n), repeat=4)):
        raise ValueError('local census is not exhaustive')
    equations = set()
    for a, b, c, d, aa, bb, cc, dd in rows:
        x, y = divmod(table[b*n+a], n)
        if (aa, bb, cc, dd) != (y, x, c, d):
            raise ValueError('shared-edge census changes a spectator or misses a face target')
        equations.add(tuple(sum(sectors[v] == label for v in (aa, bb, cc, dd))
                            - sum(sectors[v] == label for v in (a, b, c, d)) for label in labels))
    basis = nullspace(sorted(equations), len(labels))
    charges = [[0 if s == group.identity else q[labels.index(s)] for s in sectors] for q in basis]
    symbols = [tuple(q[v] for q in charges) for v in range(n)]
    alphabet = sorted(set(symbols))
    factor = induced_table(table, [alphabet.index(s) for s in symbols])
    memory = class_memory_model(group, table, charges)
    commutators = [group.mul[group.mul[a][b]][group.mul[group.inv[a]][group.inv[b]]]
                   for a, b in itertools.product(range(n), repeat=2)]
    derived = group.normal_closure(commutators)
    cosets = [min(group.mul[h][a] for h in derived) for a in range(n)]
    expected = {(trial, condition) for trial in range(result['trials'])
                for condition in ('flat', 'reflection', 'rotation', 'noncommuting_pair', 'random_links')}
    if len(result['runs']) != len(expected) or {(r['trial'], r['condition']) for r in result['runs']} != expected:
        raise ValueError('incomplete condition/trial coverage')
    for run in result['runs']:
        run.pop('independent_transport', None)
        if [r[0] for r in run['trajectory']] != list(range(result['layers']+1)):
            raise ValueError('missing trajectory layers')
        totals = []
        for row in run['trajectory']:
            hist = row[5:]
            if (len(hist) != n or any(v < 0 for v in hist) or sum(hist) != result['faces']
                    or any(hist[v] for v in range(n) if v not in sectors)
                    or row[1] != result['faces']-hist[group.identity] or row[2] > row[0]):
                raise ValueError('invalid face histogram or one-hop propagation bound')
            totals.append([sum(q[v]*hist[v] for v in range(n)) for q in charges])
            product = group.identity
            for a, multiplicity in enumerate(hist):
                for _ in range(multiplicity): product = group.mul[a][product]
            if cosets[product] != cosets[group.identity]:
                raise ValueError('closed mesh violates abelianized curvature neutrality')
        if any(t != totals[0] for t in totals):
            raise ValueError('derived global face charge drifted')
        if not run['exact_inverse_replay'] or not run['reverse_histogram_echo']:
            raise ValueError('microscopic inverse or intermediate echo failed')
        run['derived_charge_totals'] = totals[0]
    result['analysis'] = {'local_states_checked': len(rows), 'nonidentity_sector_labels': labels,
                          'additive_charge_basis': basis, 'element_charges': charges,
                          'charge_vectors': alphabet, 'charge_factor_table': factor,
                          'exchange_blocks': classify_block_exchange(factor, list(range(len(alphabet)))),
                          'derived_class_memory': memory,
                          'commutator_subgroup': derived, 'abelianization_cosets': cosets,
                          'global_charge_conservation_checked': True, 'abelianized_curvature_neutrality_checked': True,
                          'forward_updates': sum(r['updates'] for r in result['runs']),
                          'independent_replayed_updates': replay(result, table, group, memory) if full_replay else None,
                          'scope': 'fixed oriented manifold; derived conserved observables, not physical energy or electric charge'}
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--side', type=int, default=12)
    parser.add_argument('--layers', type=int, default=256)
    parser.add_argument('--trials', type=int, default=8)
    parser.add_argument('--seed', type=int, default=819031)
    parser.add_argument('--output', type=Path, default=Path('out/shared-edge.json'))
    parser.add_argument('--reanalyze', type=Path)
    parser.add_argument('--full-replay', action='store_true')
    args = parser.parse_args()
    source = json.loads(Path('data/d4-shared-mesh.json').read_text())
    census = json.loads(Path('data/d4-involution-braid-search.json').read_text())
    table = census['solutions'][source['solution_id']]['table']
    if args.reanalyze:
        result = json.loads(args.reanalyze.read_text())
        if result['solution_id'] != source['solution_id']:
            raise ValueError('rule mismatch')
    else:
        command = ['build/wgphysics_mesh_experiments', '--lift', 'shared-edge']
        for name in ('side', 'layers', 'trials', 'seed'): command += ['--'+name, str(getattr(args, name))]
        output = subprocess.run(command, input=' '.join(map(str, table))+'\n', text=True,
                                stdout=subprocess.PIPE, check=True, timeout=1800)
        result = json.loads(output.stdout)
        result['solution_id'] = source['solution_id']
    analyze(result, table, args.full_replay)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True)+'\n')
    print('Shared-edge charge rank:', len(result['analysis']['additive_charge_basis']),
          'updates:', result['analysis']['forward_updates'],
          'independently replayed:', result['analysis']['independent_replayed_updates'])
    for condition in ('reflection', 'rotation', 'noncommuting_pair'):
        runs = [r for r in result['runs'] if r['condition'] == condition]
        print(condition, 'charges:', runs[0]['derived_charge_totals'],
              'final nonflat range:', min(r['trajectory'][-1][1] for r in runs), max(r['trajectory'][-1][1] for r in runs))


if __name__ == '__main__':
    main()
