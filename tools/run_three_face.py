#!/usr/bin/env python3
"""Controlled three-face interactions, exact charges, and independent link replay."""
import argparse
import itertools
import json
import subprocess
from pathlib import Path
from face_energy_obstruction import derive_fiber_group
from run_shared_edge import mesh_geometry, bipartition
from triple_rule_search import assess
from screen_braid_rules import nullspace

MODES = ('collision', 'transport', 'combined')
CONDITIONS = ('flat', 'reflection_link', 'rotation_link', 'reflection_pair', 'opposite_reflection_pair', 'witness_left', 'witness_right', 'random_links')


def fan_geometry(side):
    edges, faces, _ = mesh_geometry(side)
    links = {}
    for face in faces:
        for i in range(3):
            u, a, b = face[i], face[(i+1) % 3], face[(i+2) % 3]
            if a in links.setdefault(u, {}): raise ValueError('duplicate oriented fan link')
            links[u][a] = b
    patches = []
    for u, next_vertex in sorted(links.items()):
        for start in sorted(next_vertex):
            rim = [start]
            while len(rim) < 4 and rim[-1] in next_vertex:
                rim.append(next_vertex[rim[-1]])
            if len(rim) == 4 and len(set(rim)) == 4:
                patches.append([(u, rim[i], rim[i+1], u) for i in range(3)])
    return edges, faces, patches


class LinkOracle:
    def __init__(self, side, group):
        self.group = group
        self.edges, self.faces, self.specs = fan_geometry(side)
        ids = {edge: i for i, edge in enumerate(self.edges)}
        def path(loop): return [(ids[tuple(sorted((u, v)))], u > v) for u, v in zip(loop, loop[1:])]
        self.face_paths = [path(face) for face in self.faces]
        self.patches = [[path(face) for face in patch] for patch in self.specs]
        self.reads = [{e for face in patch for e, _ in face} for patch in self.patches]
        self.writes = [{patch[0][-1][0], patch[1][-1][0]} for patch in self.patches]
        incident = [[] for _ in self.edges]
        for f, path in enumerate(self.face_paths):
            for edge, _ in path: incident[edge].append(f)
        if any(len(fs) != 2 for fs in incident): raise ValueError('expected a closed manifold')
        self.colors = bipartition(len(self.faces), incident)

    def transport(self, values, path):
        g = self.group
        result = g.identity
        for edge, reverse in path:
            result = g.mul[g.inv[values[edge]] if reverse else values[edge]][result]
        return result

    def sectors(self, values):
        return [self.group.sectors[self.transport(values, face)] for face in self.face_paths]

    def update(self, values, p, table):
        g, patch = self.group, self.patches[p]
        before = tuple(self.transport(values, face) for face in patch[::-1])
        index = (before[0]*g.n+before[1])*g.n+before[2]
        target = table[index]
        a, b, c = target//(g.n*g.n), (target//g.n) % g.n, target % g.n
        # Two unchanged boundary-prefix paths; no reuse of the new first spoke.
        exterior = [patch[0][0], patch[0][1]]
        s2 = g.mul[self.transport(values, exterior)][g.inv[c]]
        s3 = g.mul[self.transport(values, exterior+[patch[1][1]])][g.inv[g.mul[b][c]]]
        for closing, spoke in ((patch[0][-1], s2), (patch[1][-1], s3)):
            edge, reverse = closing
            values[edge] = spoke if reverse else g.inv[spoke]
        if tuple(self.transport(values, face) for face in patch[::-1]) != (a, b, c):
            raise ValueError('independent fan reconstruction missed a face target')
        return index, target

    def gauge_equivalent(self, left, right):
        g = self.group
        adjacent = {}
        for i, (u, v) in enumerate(self.edges):
            adjacent.setdefault(u, []).append((v, i, False))
            adjacent.setdefault(v, []).append((u, i, True))
        for root_frame in range(g.n):
            frames, queue = {0: root_frame}, [0]
            for u in queue:
                for v, i, reverse in adjacent[u]:
                    if v in frames: continue
                    a = g.inv[left[i]] if reverse else left[i]
                    b = g.inv[right[i]] if reverse else right[i]
                    frames[v] = g.mul[b][g.mul[frames[u]][g.inv[a]]]
                    queue.append(v)
            if all(g.mul[frames[v]][g.mul[left[i]][g.inv[frames[u]]]] == right[i]
                   for i, (u, v) in enumerate(self.edges)):
                return True
        return False


def analyze(result, census, full_replay=False):
    g = derive_fiber_group()
    tables = dict(zip(MODES, (census[key]['table'] for key in ('selected_rule', 'transport_control', 'combined_rule'))))
    derived = {mode: assess(g, table) for mode, table in tables.items()}
    charges = derived['combined']['element_charges']
    weights = [sum(q[a] for q in charges) for a in range(g.n)]
    labels = sorted(set(g.sectors)-{g.identity})
    staggered_equations = set()
    for code, target in enumerate(tables['combined']):
        before = (code//64, (code//8) % 8, code % 8)
        after = (target//64, (target//8) % 8, target % 8)
        for phase in (0, 1):
            parts = (phase, 1-phase, phase)
            staggered_equations.add(tuple(sum((g.sectors[a] == s)-(g.sectors[b] == s)
                                               for a, b, p in zip(before, after, parts) if p == part)
                                           for part in (0, 1) for s in labels))
    staggered_basis = nullspace(sorted(staggered_equations), 2*len(labels))
    if weights[g.identity] != 0 or any(weights[a] <= 0 for a in range(g.n) if a != g.identity):
        raise ValueError('selected derived basis does not supply the claimed positive weight')
    for table in tables.values():
        for code, target in enumerate(table):
            before = (code//64, (code//8) % 8, code % 8)
            after = (target//64, (target//8) % 8, target % 8)
            if any(sum(q[a] for a in before) != sum(q[a] for a in after) for q in charges):
                raise ValueError('a control fails to preserve the shared derived charges')
    oracle = LinkOracle(result['side'], g)
    if (result['edges'], result['faces'], result['patches']) != (len(oracle.edges), len(oracle.faces), len(oracle.patches)):
        raise ValueError('independent fan geometry disagrees with export')
    if len(result['schedules']) != result['trials']:
        raise ValueError('missing trial schedules')
    for schedule in result['schedules']:
        if not schedule or sorted(p for layer in schedule for p in layer) != list(range(len(oracle.patches))):
            raise ValueError('schedule drops or duplicates fan operators')
        for layer in schedule:
            for i, p in enumerate(layer):
                for q in layer[:i]:
                    if oracle.writes[p] & oracle.reads[q] or oracle.writes[q] & oracle.reads[p]:
                        raise ValueError('claimed parallel fan layer has a read/write conflict')
    expected = set(itertools.product(range(result['trials']), CONDITIONS, MODES))
    lookup = {(r['trial'], r['condition'], r['mode']): r for r in result['runs']}
    if set(lookup) != expected or len(lookup) != len(result['runs']):
        raise ValueError('incomplete or duplicate control runs')
    checked = 0
    for (trial, condition, mode), run in lookup.items():
        if any(run['initial_links'] != lookup[trial, condition, other]['initial_links'] for other in MODES):
            raise ValueError('controls do not share their initial raw connections')
        if [r[0] for r in run['trajectory']] != list(range(result['layers']+1)):
            raise ValueError('missing fan trajectory layers')
        values = run['initial_links'][:]
        if len(values) != len(oracle.edges) or any(not 0 <= v < g.n for v in values):
            raise ValueError('invalid initial link vector')
        initial_sectors = oracle.sectors(values)
        conserved = [sum(q[s] for s in initial_sectors) for q in charges]
        part_totals = lambda sectors: [[sum(q[s] for s, color in zip(sectors, oracle.colors) if color == part)
                                       for q in charges] for part in (0, 1)]
        conserved_parts = part_totals(initial_sectors)
        total_weight = sum(conserved)
        updates = 0
        for tick, row in enumerate(run['trajectory']):
            if len(row) != 13 or any(v < 0 for v in row) or row[4] > 2*tick:
                raise ValueError('invalid trajectory row or causal radius')
            hist = row[5:]
            if sum(hist) != result['faces'] or row[3] != result['faces']-hist[g.identity]:
                raise ValueError('invalid fan face histogram')
            if [sum(q[a]*hist[a] for a in range(g.n)) for q in charges] != conserved:
                raise ValueError('derived global fan charge drifted')
            if row[3] > total_weight:
                raise ValueError('positive conserved weight failed to bound nonflat support')
            if full_replay:
                collisions = transports = 0
                if tick:
                    order = [result['prefix_patch']] if tick == 1 else result['schedules'][trial][(tick-2) % len(result['schedules'][trial])]
                    for p in order:
                        before, after = oracle.update(values, p, tables[mode])
                        if before != after:
                            if tables['collision'][before] != before: collisions += 1
                            else: transports += 1
                        checked += 1; updates += 1
                observed = oracle.sectors(values)
                if part_totals(observed) != conserved_parts:
                    raise ValueError('bipartite sublattice charges drifted')
                if [observed.count(a) for a in range(g.n)] != hist or [collisions, transports] != row[1:3]:
                    raise ValueError('independent link replay disagrees with histogram or event classification')
        if full_replay and (values != run['final_links'] or updates != run['updates']):
            raise ValueError('independent final links or update count disagree')
        final_sectors = oracle.sectors(run['final_links'])
        if [final_sectors.count(a) for a in range(g.n)] != run['trajectory'][-1][5:]:
            raise ValueError('final raw links disagree with reported face classes')
        if not run['exact_inverse_replay'] or not run['reverse_histogram_echo']:
            raise ValueError('microscopic inverse or histogram echo failed')
        run['derived_charge_totals'] = conserved
        run['sublattice_charge_totals'] = conserved_parts
        run['positive_conserved_weight'] = total_weight
    left = lookup[0, 'witness_left', 'combined']['initial_links']
    right = lookup[0, 'witness_right', 'combined']['initial_links']
    p = result['prefix_patch']
    if oracle.sectors(left) != oracle.sectors(right):
        raise ValueError('witness initial global face classes differ')
    if any(a != b for i, (a, b) in enumerate(zip(left, right)) if i not in oracle.writes[p]):
        raise ValueError('witness connections differ outside the internal write set')
    if oracle.gauge_equivalent(left, right):
        raise ValueError('witness connections are only gauge copies')
    out_left, out_right = left[:], right[:]
    oracle.update(out_left, p, tables['combined']); oracle.update(out_right, p, tables['combined'])
    density = lambda values: [[q[s] for q in charges] for s in oracle.sectors(values)]
    if density(out_left) == density(out_right):
        raise ValueError('the link-level witness has no conserved-density feedback')
    comparisons = []
    for trial, condition in itertools.product(range(result['trials']), CONDITIONS):
        a, b = lookup[trial, condition, 'combined'], lookup[trial, condition, 'transport']
        sector_difference = sum(x != y for x, y in zip(oracle.sectors(a['final_links']), oracle.sectors(b['final_links'])))
        comparisons.append({'trial': trial, 'condition': condition,
                            'combined_collision_events': sum(r[1] for r in a['trajectory']),
                            'combined_transport_events': sum(r[2] for r in a['trajectory']),
                            'final_face_class_differences_from_transport_control': sector_difference})
    result['analysis'] = {'charge_basis': derived['combined']['charge_basis'], 'element_charges': charges,
                          'bipartite_charge_basis': staggered_basis,
                          'positive_weight_per_element': weights,
                          'independent_replayed_updates': checked if full_replay else None,
                          'forward_updates': sum(r['updates'] for r in result['runs']),
                          'controlled_comparisons': comparisons,
                          'fixed_exterior_witness': {'prefix_patch': p, 'initial_left': left, 'initial_right': right,
                                                     'after_left': out_left, 'after_right': out_right,
                                                     'same_global_initial_face_classes': True,
                                                     'same_all_exterior_links': True, 'gauge_equivalent': False},
                          'scope': 'finite classical gauge model on supplied geometry; positive conserved weight is not a calibrated physical energy'}
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--side', type=int, default=12)
    parser.add_argument('--layers', type=int, default=256)
    parser.add_argument('--trials', type=int, default=8)
    parser.add_argument('--seed', type=int, default=6290831)
    parser.add_argument('--full-replay', action='store_true')
    parser.add_argument('--reanalyze', type=Path)
    parser.add_argument('--output', type=Path, default=Path('out/three-face.json'))
    args = parser.parse_args()
    census = json.loads(Path('data/d4-triple-rule-search.json').read_text())
    if args.reanalyze:
        result = json.loads(args.reanalyze.read_text())
        if result['selected_rule_id'] != census['selected_rule']['id']: raise ValueError('selected rule differs')
    else:
        inputs = [str(x) for key in ('selected_rule', 'transport_control', 'combined_rule') for x in census[key]['table']]
        command = ['build/wgphysics_three_face_experiments']
        for name in ('side', 'layers', 'trials', 'seed'): command += ['--'+name, str(getattr(args, name))]
        process = subprocess.run(command, input=' '.join(inputs)+'\n', text=True, stdout=subprocess.PIPE, check=True, timeout=1800)
        result = json.loads(process.stdout)
        result['selected_rule_id'] = census['selected_rule']['id']
    analyze(result, census, args.full_replay)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True)+'\n')
    print('Three-face charge basis:', result['analysis']['charge_basis'], 'updates:', result['analysis']['forward_updates'])
    print('Independent replay:', result['analysis']['independent_replayed_updates'])
    for condition in CONDITIONS:
        rows = [r for r in result['analysis']['controlled_comparisons'] if r['condition'] == condition]
        print(condition, 'collision range:', min(r['combined_collision_events'] for r in rows), max(r['combined_collision_events'] for r in rows),
              'final control difference:', [r['final_face_class_differences_from_transport_control'] for r in rows])


if __name__ == '__main__':
    main()
