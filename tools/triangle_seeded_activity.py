#!/usr/bin/env python3
"""Autonomous spreading from a localized nonabelian seed in an inactive background.

The background is a gauge transform of a constant reflection connection.
Geometry is supplied; dual-face distances are diagnostics, not emergent space.
"""
import json
import random
from collections import deque
from pathlib import Path

import numpy as np

from triangle_elastic_scattering import ElasticExperiment
from triangle_reference import reference
from triangle_released_traps import ReleasedTraps


def growth(r, side=8, trials=32, times=(0, 1, 2, 4, 8, 16, 32, 64, 128)):
    # These two fixed preparations are the strongest and weakest initial
    # reaction-growth classes of the complete nine-state alignment census.
    prepared, _ = r.embed('parallel_reflections', 0)
    chosen = (0, 3)
    bank = json.loads(Path('data/triangle-feedback.json').read_text())['bank']
    x = r.x if side == 8 else ElasticExperiment(side, bank, r.x.tables[1:3])
    count = len(x.e.geometry.faces)
    offset = side//2-3
    vertex = lambda v: (v//8+offset)*side+(v % 8+offset)
    initials = []
    for index in chosen:
        raw = [r.s.p]*len(x.e.geometry.edges)
        for edge, value in enumerate(prepared[index]):
            if value == r.s.p:
                continue
            a, b = r.x.e.geometry.edges[edge]
            assert abs(a % 8-b % 8) <= 1 and abs(a//8-b//8) <= 1
            raw[x.e.geometry.ids[tuple(sorted((vertex(a), vertex(b))))]] = value
        q = x.p.charge(raw)
        assert sum(q) == count and q.count(0) == q.count(2) == 1
        initials.append(raw)
    assert x.p.charge(initials[0]) == x.p.charge(initials[1])
    region = [2*((f//2//8+offset)*side+(f//2 % 8+offset))+f % 2 for f in r.s.faces]
    neighbors = [set() for _ in range(count)]
    for a, b in x.p.supports[0]:
        neighbors[a].add(b); neighbors[b].add(a)
    distance = np.full(count, -1, dtype=int)
    distance[region] = 0
    queue = deque(region)
    while queue:
        a = queue.popleft()
        for b in neighbors[a]:
            if distance[b] < 0:
                distance[b] = distance[a]+1
                queue.append(b)
    samples = []
    charges = x.e.charges
    for trial in range(trials):
        rng = random.Random(103725000+trial)
        schedule = [rng.randrange(x.operator_count) for _ in range(count*times[-1])]
        rows = x.compiled_output(initials, schedule, count, raw_events=False)['runs']
        pair = []
        for condition, initial in enumerate(initials):
            row = rows[2*condition+1]
            assert row['condition'] == condition and row['mode'] == 'combined'
            q = np.array(x.p.charge(initial))
            event_index, values = 0, []
            for time in times:
                while event_index < len(row['events']) and row['events'][event_index][0] <= time*count:
                    _, rule, patch, code, target = row['events'][event_index]
                    faces = x.p.supports[int(rule >= 3)][patch]
                    before = [charges[(code//r.s.g.n**k) % r.s.g.n] for k in reversed(range(len(faces)))]
                    after = [charges[(target//r.s.g.n**k) % r.s.g.n] for k in reversed(range(len(faces)))]
                    assert q[list(faces)].tolist() == before
                    q[list(faces)] = after
                    event_index += 1
                active = q != 1
                rotations = int(np.sum(q == 2))
                assert rotations == np.sum(q == 0)
                histogram = row['histograms'][time]
                assert rotations == sum(c for a, c in enumerate(histogram) if charges[a] == 2)
                values.append([rotations, float(np.mean(distance[active]**2)) if active.any() else 0,
                               int(np.max(distance[active])) if active.any() else 0])
            pair.append(values)
        samples.append(pair)
    samples = np.array(samples)
    mean = samples.mean(axis=0)
    se = samples.std(axis=0, ddof=1)/np.sqrt(trials)
    equilibrium = reference(count, count, 'nonabelian_reflections')['mean_rotation_count']
    return {'side': side, 'faces': count, 'trials': trials,
            'preparation_indices': chosen, 'times_attempts_per_face': times,
            'metrics': ['rotation_count', 'nonreflection_mean_squared_dual_distance_from_seed_region',
                        'nonreflection_maximum_dual_distance_from_seed_region'],
            'mean': mean.tolist(), 'SE': se.tolist(),
            'paired_rotation_count_difference': (samples[:, 0, :, 0]-samples[:, 1, :, 0]).mean(axis=0).tolist(),
            'paired_rotation_count_difference_SE': ((samples[:, 0, :, 0]-samples[:, 1, :, 0]).std(axis=0, ddof=1)/np.sqrt(trials)).tolist(),
            'uniform_sector_reference_mean_rotation_count': equilibrium[0]/equilibrium[1],
            'scope': 'Prepared local gauge quench in a supplied periodic triangular mesh; not an equilibrium seed ensemble or a front-speed exponent.'}


if __name__ == '__main__':
    r = ReleasedTraps()
    for side in (8, 16):
        print(json.dumps(growth(r, side)), flush=True)
