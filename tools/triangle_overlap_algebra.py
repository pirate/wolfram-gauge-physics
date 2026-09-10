#!/usr/bin/env python3
"""Finite overlap algebra on the actual five-face shared-edge disk.

All primitive read supports lie inside the disk. This is a boundary-held
experiment, not a closed quotient of the open whole-mesh dynamics.
"""
import json
from collections import Counter, defaultdict, deque
from pathlib import Path

from flint import fmpq_mat

from triangle_elastic_scattering import ElasticExperiment, inverse
from triangle_patch_observer import PatchUnion
from triangle_reference import subgroup


def build():
    bank = json.loads(Path('data/triangle-feedback.json').read_text())['bank']
    elastic = json.loads(Path('data/triangle-elastic-scattering.json').read_text())['pair_census']['elastic_tables']
    x = ElasticExperiment(6, bank, elastic)
    u = PatchUnion(x.e, (0, 2))
    oracle, edges = x.e.factor.oracle, set(u.edges)
    pairs = [p for p, paths in enumerate(oracle.pairs)
             if {edge for path in paths for edge, _ in path} <= edges]
    fans = [p for p, reads in enumerate(oracle.fan.reads) if reads <= edges]
    channels = [(p, r) for p in pairs for r in range(3)]
    channels += [(p, r) for p in fans for r in range(3, len(x.tables))]
    tables = []
    for p, r in channels:
        row = []
        for rep in u.reps:
            raw = u.lift(rep)
            x.apply(raw, r*x.supports+p)
            row.append(u.observe(raw))
        inverse(row)
        tables.append(row)
    assert all(inverse(t) in tables for t in tables)
    faces = sorted(set().union(*(x.e.fan_faces[p] for p in u.patches)))
    directed = [(a, b) for f in faces for a, b in zip(x.e.geometry.faces[f], x.e.geometry.faces[f][1:])]
    counts = Counter(tuple(sorted(edge)) for edge in directed)
    outer = [(a, b) for a, b in directed if counts[tuple(sorted((a, b)))] == 1]
    successor = dict(outer)
    assert len(successor) == len(outer)
    root = min(successor)
    boundary = [root]
    for _ in outer:
        boundary.append(successor[boundary[-1]])
    assert boundary[-1] == root and len(set(boundary[:-1])) == len(outer)
    path = [(x.e.geometry.ids[tuple(sorted((a, b)))], a > b) for a, b in zip(boundary, boundary[1:])]
    writes = set().union(*(oracle.fan.writes[p] for p in fans),
                         *({oracle.pairs[p][0][0][0]} for p in pairs))
    assert not writes.intersection(edge for edge, _ in path)
    return x, u, channels, tables, faces, path


def components(tables):
    seen, result = set(), []
    for start in range(len(tables[0])):
        if start in seen:
            continue
        block = [start]
        seen.add(start)
        for i in block:
            for row in tables:
                j = row[i]
                if j not in seen:
                    seen.add(j)
                    block.append(j)
        result.append(sorted(block))
    return result


def reflection_return_channel(u, channels, tables, sector, reflection):
    """Exact first-return elimination of the 240 mixed states, with its clock."""
    r, e = sorted(reflection), sorted(sector-reflection)
    counts = lambda left, right: fmpq_mat([[sum(t[i] == j for t in tables) for j in right] for i in left])
    clock = len(tables)
    rr, re, er, ee = counts(r, r), counts(r, e), counts(e, r), counts(e, e)
    identity = fmpq_mat([[int(i == j)*clock for j in e] for i in e])
    green = (identity-ee).inv()
    hitting = green*er
    trace = (rr+re*hitting)/clock
    ones_e, ones_r = fmpq_mat([[1] for _ in e]), fmpq_mat([[1] for _ in r])
    assert hitting*ones_r == ones_e and trace*ones_r == ones_r
    assert trace.transpose() == trace
    assert all(trace[i, j] >= 0 for i in range(len(r)) for j in range(len(r)))
    mean = ones_r+re*green*ones_e
    average = sum(mean[i, 0] for i in range(len(r)))/len(r)
    assert average == 7 and len(sector) == 7*len(r)
    exits = sum(re[i, j] for i in range(len(r)) for j in range(len(e)))
    rids = {state: i for i, state in enumerate(r)}
    witness = None
    for channel, t in zip(channels, tables):
        if channel[1] not in (1, 2):
            continue
        perm = [rids[t[i]] for i in r]
        for i in range(len(r)):
            for j in range(len(r)):
                if trace[i, j] != trace[perm[i], perm[j]]:
                    witness = {'channel': channel, 'reflection_state_ids': [r[i], r[j]],
                               'conjugated_state_ids': [r[perm[i]], r[perm[j]]],
                               'probabilities': [str(trace[i, j]), str(trace[perm[i], perm[j]])]}
                    break
            if witness:
                break
        if witness:
            break
    assert witness is not None
    return {'reflection_states': len(r), 'eliminated_mixed_states': len(e),
            'reflection_state_ids': r,
            'trace_kernel_rational': [[str(trace[i, j]) for j in range(len(r))] for i in range(len(r))],
            'uniform_reflection_mean_first_positive_return_attempts': str(average),
            'uniform_reflection_departure_probability': str(exits/(clock*len(r))),
            'departure_flux_weighted_mixed_hitting_time': str(clock*(len(sector)-len(r))/exits),
            'statewise_mean_return_attempts': [str(mean[i, 0]) for i in range(len(r))],
            'trace_not_invariant_under_elastic_group_witness': witness,
            'exact_nonnegative_symmetric_stochastic_return_kernel': True}


def calculation():
    x, u, channels, tables, faces, path = build()
    labels, charges = [], []
    for rep in u.reps:
        raw = u.lift(rep)
        q = x.p.charge(raw)
        charges.append(tuple(q[f] for f in faces))
        total = sum(charges[-1])
        boundary = x.e.factor.oracle.fan.transport(raw, path)
        labels.append((total, x.e.charges[boundary], len(subgroup(u.g, rep))))
    assert all(labels[i] == labels[t[i]] for t in tables for i in range(len(u.reps)))
    blocks = components(tables)
    by_label = defaultdict(list)
    for block in blocks:
        by_label[labels[block[0]]].append(block)
    assert all(len(v) == 1 for k, v in by_label.items() if k[0] <= 5)
    # A constructive route: convert to all reflections without ever increasing
    # the number of rotations, then use the elastic reflection orbit.
    sector = {i for i, label in enumerate(labels) if label == (5, 1, 6)}
    level = [q.count(2) for q in charges]
    reflection = {i for i in sector if level[i] == 0}
    distances = {i: 0 for i in reflection}
    parents, queue = {}, deque(reflection)
    inverse_slots = [tables.index(inverse(t)) for t in tables]
    while queue:
        i = queue.popleft()
        for k, table in enumerate(tables):
            j = table[i]
            if j in sector and j not in distances and level[j] >= level[i]:
                distances[j] = distances[i]+1
                parents[j] = (i, inverse_slots[k])
                queue.append(j)
    assert set(distances) == sector
    for state in sector-reflection:
        raw, i = u.lift(u.reps[state]), state
        while i not in reflection:
            target, channel = parents[i]
            p, r = channels[channel]
            x.apply(raw, r*x.supports+p)
            assert u.observe(raw) == target and level[target] <= level[i]
            i = target
        assert all(x.p.charge(raw)[f] == 1 for f in faces)
    elastic_orbits = components([t for (p, r), t in zip(channels, tables) if r < 3])
    assert reflection in [set(b) for b in elastic_orbits]
    return {'gauge_states': len(u.reps), 'faces': faces, 'outer_path': path,
            'channels': channels, 'components': len(blocks),
            'invariant_labels': 'total charge, outer boundary charge, complete based holonomy image order',
            'charge5_constructive_reduction': {'states': len(sector), 'all_reflection_states': len(reflection),
                'rotation_count_populations': dict(Counter(level[i] for i in sector)),
                'shortest_monotone_distances': dict(sorted(Counter(distances.values()).items())),
                'maximum_monotone_distance': max(distances.values()),
                'all_monotone_paths_replayed_on_actual_raw_links': True},
            'reflection_first_return': reflection_return_channel(u, channels, tables, sector, reflection),
            'sectors': [{'label': k, 'sizes': [len(b) for b in v],
                         'representatives': [u.reps[b[0]] for b in v]}
                        for k, v in sorted(by_label.items())]}


if __name__ == '__main__':
    print(json.dumps(calculation(), indent=2))
