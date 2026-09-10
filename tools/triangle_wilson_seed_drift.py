#!/usr/bin/env python3
"""Exact first-step primitive contributions to the nonstationary loop probes."""
import json
from fractions import Fraction
from pathlib import Path

import numpy as np

from triangle_wilson_quench import LoopObserver, seeds


def observables(observer, raw, sector, shape=4):
    q = np.array(observer.x.p.charge(raw))
    n = np.bincount(q, minlength=3)
    loops = observer.holonomies(raw, shape)[0]
    actual = Fraction(sum(observer.theory.characters[2][int(h)] for h in loops), 2*len(loops))
    projected = []
    for field in q[observer.regions[shape]]:
        if sector == 2:
            projected.append(Fraction(1+(-1)**int(field.sum()), 2))
        else:
            inside = tuple(int((field == k).sum()) for k in range(3))
            outside = tuple(int(n[k])-inside[k] for k in range(3))
            projected.append(observer.theory.conditional_standard(inside, outside))
    projected = sum(projected)/len(projected)
    return (actual, projected, actual-projected, Fraction(int(n[2])))


def calculation():
    observer = LoopObserver(6)
    result = []
    labels = ('standard_loop_area8', 'charge_projection_area8', 'gauge_residual_area8', 'rotation_population')
    for prepared in seeds(observer):
        raw, sector = prepared['initial_links'], prepared['holonomy_group_order']
        before = observables(observer, raw, sector)
        sums = {name: [Fraction(0)]*4 for name in ('vacancy', 'elastic', 'reaction')}
        changing = dict.fromkeys(sums, 0)
        for op in range(observer.x.operator_count):
            moved = raw[:]
            a, b = observer.x.apply(moved, op)
            if a == b:
                continue
            rule = op//observer.x.supports
            family = 'vacancy' if rule == 0 else 'elastic' if rule < 3 else 'reaction'
            after = observables(observer, moved, sector)
            changing[family] += 1
            for i, (new, old) in enumerate(zip(after, before)):
                sums[family][i] += new-old
        # A derivative at the initial state means F times the exact one-attempt
        # expectation increment, not a fit to a finite elapsed-time interval.
        scale = Fraction(observer.faces, observer.x.operator_count)
        record = {'seed': prepared['name'], 'initial': {k: str(v) for k, v in zip(labels, before)},
                  'changing_slots': changing, 'total_slots': observer.x.operator_count,
                  'scaled_first_step_by_family': {name: {k: str(scale*v) for k, v in zip(labels, values)}
                                                  for name, values in sums.items()},
                  'scaled_first_step_total': {k: str(scale*sum(v[i] for v in sums.values())) for i, k in enumerate(labels)}}
        if sector == 2:
            lap = np.zeros((observer.faces, observer.faces), dtype=int)
            for u, v in observer.x.p.supports[0]:
                lap[u, u] += 1; lap[v, v] += 1
                lap[u, v] -= 1; lap[v, u] -= 1
            assert np.all(np.diag(lap) == 6) and set(lap[lap < 0]) == {-2}
            eigenvalues, vectors = np.linalg.eigh(lap)
            assert eigenvalues.min() > -1e-10
            eigenvalues = np.maximum(eigenvalues, 0)
            q0 = np.array(prepared['initial_charge_field'], dtype=float)
            times = [0, 1, 4, 16, 64, 256, 1024, 4096]
            fields = [vectors@((1-eigenvalues/observer.x.operator_count)**(t*observer.faces)*(vectors.T@q0))
                      for t in times]
            record['C2_exact_mean_heat_equation'] = {'times_attempts_per_face': times, 'mean_fields': [f.tolist() for f in fields],
                  'mean_field_inhomogeneity': [float(np.mean((f-18/observer.faces)**2)) for f in fields],
                  'nearest_neighbor_coefficient_per_attempt': str(Fraction(2, observer.x.operator_count))}
        result.append(record)
    return {'observables': labels, 'shape': observer.shapes[4], 'initial_drifts': result}


if __name__ == '__main__':
    result = calculation()
    output = Path('data/triangle-wilson-seed-drift.json')
    output.write_text(json.dumps(result, indent=2)+'\n')
    for row in result['initial_drifts']:
        print(json.dumps({k: v for k, v in row.items() if k != 'C2_exact_mean_heat_equation'}), flush=True)
