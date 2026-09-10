#!/usr/bin/env python3
"""Nonstationary raw-link seeds and charge-orthogonal Wilson relaxation.

Initial preparation changes nine explicit links; only existing primitive
rules govern all later evolution. Whole trajectories are independent units.
"""
import argparse
import itertools
import json
import random
import time
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor
from fractions import Fraction
from pathlib import Path

import numpy as np

from triangle_reference import reference, subgroup
from triangle_wilson_dynamics import LoopObserver


def seeds(observer):
    g, geometry = observer.g, observer.x.e.geometry
    reflections = [h for h in range(g.n) if observer.x.e.charges[h] == 1]
    result = []
    for name, spacing, parallel in (('cluster_S3', 1, False), ('spread_S3', 2, False), ('cluster_C2', 1, True)):
        raw, prepared = [g.identity]*len(geometry.edges), []
        for j, i in itertools.product(range(3), repeat=2):
            x, y = spacing*i, spacing*j
            u, v = y*observer.side+x, y*observer.side+(x+1) % observer.side
            edge = geometry.ids[tuple(sorted((u, v)))]
            value = reflections[0 if parallel else (i+2*j) % 3]
            raw[edge] = value
            prepared.append([edge, value])
        q = observer.x.p.charge(raw)
        image = subgroup(g, observer.x.e.forest.based_loops(raw)[0])
        assert sum(q) == 18 and q.count(1) == 18 and q.count(2) == 0
        assert len(image) == (2 if parallel else 6)
        result.append({'name': name, 'spacing_cells': spacing, 'holonomy_group_order': len(image),
                       'initial_links': raw, 'prepared_edge_values': prepared, 'initial_charge_field': q})
    assert result[0]['initial_charge_field'] == result[2]['initial_charge_field']
    return result


def conditional_count(theory, faces=4):
    """Independent finite sum grouped by prescribed inside/outside charges."""
    g = theory.g
    handles = defaultdict(list)
    for a, b in itertools.product(range(g.n), repeat=2):
        handles[g.mul[g.inv[b]][g.mul[g.inv[a]][g.mul[b][a]]]].append((a, b))
    closure = {mask: len(subgroup(g, [h for h in range(g.n) if mask >> h & 1]))
               for mask in range(1 << g.n)}
    counts, traces = Counter(), Counter()
    for word in itertools.product(range(g.n), repeat=faces):
        q = [theory.q[h] for h in word]
        if 1 not in q:
            continue
        product = [g.identity]
        for h in word:
            product.append(g.mul[product[-1]][h])
        mask = sum(1 << h for h in set(word))
        for a, b in handles[product[-1]]:
            if closure[mask | (1 << a) | (1 << b)] != g.n:
                continue
            for area in range(faces+1):
                key = (tuple(q[:area].count(i) for i in range(3)), tuple(q[area:].count(i) for i in range(3)))
                counts[key] += 1
                traces[key] += theory.characters[2][product[area]]
    for key, count in counts.items():
        assert theory.conditional_standard(*key) == Fraction(traces[key], 2*count)
    return {'faces': faces, 'distinct_population_conditions': len(counts),
            'all_conditional_characters_equal_direct_count': True}


def projected_loops(observer, q, sector):
    result = []
    for region in observer.regions:
        values = []
        for field in q:
            n = np.bincount(field, minlength=3)
            inside_q = field[region]
            if sector == 2:
                values.append(np.mean((1+(-1)**inside_q.sum(axis=1))/2))
                continue
            counts = np.stack([(inside_q == k).sum(axis=1) for k in range(3)], axis=1)
            values.append(np.mean([float(observer.theory.conditional_standard(tuple(int(x) for x in c),
                                      tuple(int(x) for x in n-c))) for c in counts]))
        result.append(values)
    return np.array(result).T


def trajectory(observer, prepared, schedule_seed, times):
    rng = random.Random(schedule_seed)
    attempts = times[-1]*observer.faces
    schedule = [rng.randrange(observer.x.operator_count) for _ in range(attempts)]
    row = observer.x.compiled_output([prepared['initial_links']], schedule, attempts, raw_events=True)['runs'][1]
    frames, cursor, current = [prepared['initial_links']], 0, prepared['initial_links']
    for t in times[1:]:
        while cursor < len(row['events']) and row['events'][cursor][0] <= t*observer.faces:
            current = row['events'][cursor][5]
            cursor += 1
        frames.append(current)
    assert current == row['final_links']
    q = np.array([observer.x.p.charge(raw) for raw in frames])
    assert np.all(q.sum(axis=1) == 18)
    measured = observer.observe(frames)
    projection = projected_loops(observer, q, prepared['holonomy_group_order'])
    residual = measured[:, :, 1]-projection
    if prepared['holonomy_group_order'] == 2:
        assert np.all(q != 2) and np.max(np.abs(residual)) < 1e-14
    return {'wilson': measured.tolist(), 'charge_projection': projection.tolist(), 'gauge_residual': residual.tolist(),
            'rotation_population': (q == 2).sum(axis=1).tolist(), 'charge_fields': q.tolist(),
            'changing_events': len(row['events'])}


def measure(observer, prepared, trials, seed, times, workers=4):
    started = time.monotonic()
    records = []
    run = lambda i: trajectory(observer, prepared, seed+i, times)
    with ThreadPoolExecutor(max_workers=workers) as pool:
        for i, row in enumerate(pool.map(run, range(trials)), 1):
            records.append(row)
            if i % 8 == 0:
                print(json.dumps({'seed': prepared['name'], 'progress': i, 'trials': trials,
                                  'elapsed_seconds': time.monotonic()-started}), flush=True)
    summary = {}
    for key in ('wilson', 'charge_projection', 'gauge_residual', 'rotation_population', 'charge_fields'):
        values = np.array([r[key] for r in records], dtype=float)
        summary[key] = {'mean': values.mean(axis=0).tolist(),
                        'trajectory_standard_error': (values.std(axis=0, ddof=1)/np.sqrt(trials)).tolist()}
    q = np.array([r['charge_fields'] for r in records], dtype=float)
    mean_q = q.mean(axis=0)
    # Unbiased ensemble mean-field inhomogeneity; finite-sample estimates can
    # be negative near zero. This is not the variance of a single configuration.
    inhomogeneity = ((mean_q-18/observer.faces)**2-q.var(axis=0, ddof=1)/trials).mean(axis=1)
    return {'preparation': prepared, 'trials': trials, 'schedule_seed_base': seed, 'times_attempts_per_face': times,
            'observables': summary, 'charge_mean_field_inhomogeneity_unbiased': inhomogeneity.tolist(),
            'combined_forward_attempts': trials*times[-1]*observer.faces,
            'combined_changing_events': sum(r['changing_events'] for r in records), 'trajectory_records': records,
            'scope': 'Nonstationary finite-seed ensemble. Agreement of selected observables is not whole-sector ergodicity or molecular binding.'}


def c2_means(faces, charge, area):
    import math
    sign = sum((-1)**k*math.comb(charge, k)*math.comb(faces-charge, area-k)
               for k in range(max(0, area-(faces-charge)), min(area, charge)+1))
    sign = Fraction(sign, math.comb(faces, area))
    return [float(sign), float((1+sign)/2)]


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--trials', type=int, default=128)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    observer = LoopObserver(6)
    times = [0, 1, 4, 16, 64, 256, 1024, 4096]
    report = {'conditional_count': conditional_count(observer.theory), 'shapes': observer.shapes,
              'S3_stationary_loop_means': [[float(observer.theory.mean(72, 18, s['area_faces'], t)) for t in (1, 2)]
                                           for s in observer.shapes],
              'C2_stationary_loop_means': [c2_means(72, 18, s['area_faces']) for s in observer.shapes],
              'S3_stationary_mean_rotations': reference(72, 18, 'nonabelian_reflections')['mean_rotation_count'],
              'experiments': []}
    for initial, seed in zip(seeds(observer), (329994001, 341105001, 352216001)):
        report['experiments'].append(measure(observer, initial, args.trials, seed, times))
        Path(args.output).write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({'output': args.output, 'complete': True}), flush=True)
