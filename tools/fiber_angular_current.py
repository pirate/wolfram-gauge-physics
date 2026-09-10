#!/usr/bin/env python3
"""Exact angular readout of winding-current memory on a real periodic mesh."""
import itertools
import json
from collections import Counter, deque
from pathlib import Path

import numpy as np

from fiber_bulk_transport import Experiment, current_drift


def based_activity(e, raw, p):
    w = [e.oracle.transport(raw, path) for path in e.oracle.patches[p][::-1]]
    return int(all(h in e.f.reflections for h in w) and ((w[0] == w[1]) != (w[1] == w[2])))


def current_response(e, raw):
    q = np.array([e.charge[e.oracle.transport(raw, p)] for p in e.oracle.face_paths])
    b = current_drift(q, e.supports, e.steps)
    result = np.zeros(2, dtype=np.int64)
    for p, support in enumerate(e.supports):
        a, b0, c = [e.oracle.transport(raw, path) for path in e.oracle.patches[p][::-1]]
        code = (a*e.g.n+b0)*e.g.n+c
        targets = Counter(int(table[code]) for table in e.tables[:15])
        original = q[support].copy()
        for target, multiplicity in targets.items():
            q[support] = e.charge[[target//(e.g.n**2), target//e.g.n % e.g.n, target % e.g.n]]
            result += multiplicity*(b-current_drift(q, e.supports, e.steps))
        q[support] = original
    return result


def witness(n=5, side=12, seed=519786001):
    e = Experiment(n, side)
    rng = np.random.default_rng(seed)
    raw = rng.integers(e.g.n, size=len(e.oracle.edges)).tolist()
    q = np.array([e.charge[e.oracle.transport(raw, p)] for p in e.oracle.face_paths])
    b = current_drift(q, e.supports, e.steps)
    weights = np.zeros((len(e.supports), 2), dtype=np.int64)
    for p, support in enumerate(e.supports):
        if np.all(q[support] == 1):
            for layout in itertools.permutations((0, 1, 2)):
                changed_q = q.copy()
                changed_q[support] = layout
                weights[p] -= 2*(current_drift(changed_q, e.supports, e.steps)-b)
    activities = np.array([based_activity(e, raw, p) for p in range(len(e.supports))])
    for p in range(len(e.supports)):
        for channel, rule in enumerate(e.f.rules):
            after = raw[:]
            e.oracle.update(after, p, rule)
            changed_edges = {edge for edge in e.oracle.writes[p] if raw[edge] != after[edge]}
            if not changed_edges:
                continue
            affected = [g for g, reads in enumerate(e.oracle.reads) if reads & changed_edges]
            delta = sum((weights[g]*(based_activity(e, after, g)-activities[g]) for g in affected),
                        np.zeros(2, dtype=np.int64))
            if delta[0] == 0:
                continue
            assert len(changed_edges) == 1
            after_q = np.array([e.charge[e.oracle.transport(after, path)] for path in e.oracle.face_paths])
            assert np.all(after_q == q)
            u0, u1 = current_response(e, raw), current_response(e, after)
            assert np.all(u1-u0 == delta)
            # Exact read set for this difference, including the boundary frames.
            reads = set(e.oracle.reads[p])
            for g in affected:
                reads.update(e.oracle.reads[g])
                faces = set(e.supports[g])
                for support in e.supports:
                    if faces & set(support):
                        for face in support:
                            reads.update(edge for edge, _ in e.oracle.face_paths[face])
            adjacent = {}
            for edge in reads:
                a, c = e.oracle.edges[edge]
                va, vc = np.array((a % side, a//side)), np.array((c % side, c//side))
                shift = (vc-va+side//2) % side-side//2
                adjacent.setdefault(a, []).append((c, shift))
                adjacent.setdefault(c, []).append((a, -shift))
            root = min(adjacent)
            lifted, queue = {root: np.zeros(2, dtype=np.int64)}, deque([root])
            while queue:
                a = queue.popleft()
                for c, shift in adjacent[a]:
                    target = lifted[a]+shift
                    if c in lifted:
                        assert np.all(lifted[c] == target)
                    else:
                        lifted[c] = target
                        queue.append(c)
            assert len(lifted) == len(adjacent)
            coordinates = np.array(list(lifted.values()))
            extent = coordinates.max(axis=0)-coordinates.min(axis=0)
            assert np.all(extent < side)
            return {'cycle_vertices': n, 'side': side, 'faces': len(q), 'seed': seed,
                    'angular_patch': p, 'angular_channel_index': channel, 'changed_edges': sorted(changed_edges),
                    'initial_links': raw, 'final_links': after, 'edges': e.oracle.edges,
                    'charge_field_before_and_after': q.tolist(), 'winding_drift_times_three_before_and_after': b.tolist(),
                    'old_generator_on_winding_drift_times_three_before': u0.tolist(),
                    'old_generator_on_winding_drift_times_three_after': u1.tolist(),
                    'difference': delta.tolist(), 'read_edge_indices': sorted(reads),
                    'read_edge_count': len(reads), 'unwrapped_read_extent': extent.tolist(),
                    'read_assignment': [[edge, raw[edge]] for edge in sorted(reads)],
                    'bulk_lower_bound_on_A_b_B_A_b_per_face': {
                        'numerator': int(delta[0])**2, 'denominator_factor': 36,
                        'group_order': e.g.n, 'group_order_exponent': len(reads)},
                    'scope': 'Exact local witness for angular sensitivity of winding-current memory. The lower bound concerns a short-time Taylor coefficient under uniform raw connections, not a positive DC mobility difference.'}
    raise ValueError('The specified reference draw has no angular current witness')


def bulk_dirichlet_moment(n=5, side=12, samples=20000, seed=620891001):
    e = Experiment(n, side)
    rng = np.random.default_rng(seed)
    affected_by_edge = [[p for p, reads in enumerate(e.oracle.reads) if edge in reads]
                        for edge in range(len(e.oracle.edges))]
    statistics = []
    nonzero = 0
    for _ in range(samples):
        raw = rng.integers(e.g.n, size=len(e.oracle.edges)).tolist()
        p, channel = int(rng.integers(len(e.supports))), int(rng.integers(e.f.m))
        after = raw[:]
        e.oracle.update(after, p, e.f.rules[channel])
        changed = [edge for edge in e.oracle.writes[p] if raw[edge] != after[edge]]
        delta = np.zeros(2, dtype=np.int64)
        if changed:
            assert len(changed) == 1
            q = np.array([e.charge[e.oracle.transport(raw, path)] for path in e.oracle.face_paths])
            b = current_drift(q, e.supports, e.steps)
            for g in affected_by_edge[changed[0]]:
                if not np.all(q[e.supports[g]] == 1):
                    continue
                difference = based_activity(e, after, g)-based_activity(e, raw, g)
                if difference:
                    for layout in itertools.permutations((0, 1, 2)):
                        changed_q = q.copy()
                        changed_q[e.supports[g]] = layout
                        delta -= 2*difference*(current_drift(changed_q, e.supports, e.steps)-b)
        nonzero += int(np.any(delta))
        # 3F rooted patches, m angular channels, Dirichlet factor 1/2,
        # and b3=3b: (3m/2)/9 = m/6 per face.
        statistics.append(e.f.m*np.array([delta[0]**2, delta[0]*delta[1], delta[1]**2])/6.)
    statistics = np.array(statistics)
    return {'cycle_vertices': n, 'side': side, 'independent_local_samples': samples, 'seed': seed,
            'nonzero_current_response_samples': nonzero,
            'A_b_B_A_b_per_face_xx_xy_yy': statistics.mean(axis=0).tolist(),
            'independent_sample_SE': (statistics.std(axis=0, ddof=1)/np.sqrt(samples)).tolist(),
            'scope': 'Direct stationary local Dirichlet moment; coefficient of short-time winding response, not DC mobility.'}


if __name__ == '__main__':
    result = witness()
    Path('data/fiber-angular-current.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k: result[k] for k in ('faces', 'difference', 'read_edge_count', 'unwrapped_read_extent',
                                           'winding_drift_times_three_before_and_after',
                                           'old_generator_on_winding_drift_times_three_before',
                                           'old_generator_on_winding_drift_times_three_after')}))
