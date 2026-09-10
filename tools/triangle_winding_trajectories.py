#!/usr/bin/env python3
"""Winding fluctuations under the existing microscopic raw-link dynamics.

Exact stationary starts, independent uniform proposal schedules, and currents
reconstructed from the existing engine's actual events. No diffusion law is
supplied. Windows within a trajectory are clustered for uncertainty estimates.
"""
import json
import itertools
import random
import time
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

from triangle_elastic_scattering import ElasticExperiment
from triangle_reference import AxialReference
from triangle_charge_response import ChargeMarginal


class Winding:
    def __init__(self, side):
        self.bank = json.loads(Path('data/triangle-feedback.json').read_text())['bank']
        elastic = json.loads(Path('data/triangle-elastic-scattering.json').read_text())['pair_census']['elastic_tables']
        self.x = ElasticExperiment(side, self.bank, elastic)
        self.faces, self.side = self.x.geometry.faces, side
        self.coordinates = np.array([[3*(f//2 % side)+(2 if f % 2 == 0 else 1),
                                      3*(f//2//side)+(1 if f % 2 == 0 else 2)]
                                     for f in range(self.faces)])
        period = 3*side
        displacements = np.array([((self.coordinates[b]-self.coordinates[a]+period//2) % period)-period//2
                                  for a, b in self.x.p.dual])
        self.steps = [np.array([[displacements[e]*sign for e, sign in path]
                                for path in paths]) for paths in self.x.p.paths]
        self.charges = np.array(self.x.e.charges)

    def one_step_moments(self):
        fourth_counts = Counter()
        for arity in (0, 1):
            for steps in self.steps[arity]:
                states = ((0, 1), (1, 0), (0, 2), (2, 0)) if not arity else itertools.permutations((0, 1, 2))
                for q in states:
                    j = (int(steps[0, 0])*(q[0]-q[1]) if not arity else
                         int(steps[0, 0])*(q[0]-1)+int(steps[1, 0])*(1-q[2]))
                    # Four death slots, with their reverse births equal in
                    # stationary fourth moment; no conditional gate average.
                    fourth_counts[tuple(q.count(k) for k in range(3))] += (1 if not arity else 8)*j**4
        marginal = ChargeMarginal(self.faces, self.faces, 'nonabelian_reflections')
        fourth = marginal.expectation(fourth_counts.items())/(81*45*self.faces)
        second = (12*marginal.probability((1, 1, 0))+48*marginal.probability((1, 0, 1))
                  +480*marginal.probability((1, 1, 1)))/(9*45)
        return {'faces': self.faces, 'second_moment': float(second), 'fourth_moment': float(fourth),
                'j3_fourth_coefficients_per_face': [[list(k), v//self.faces] for k, v in fourth_counts.items()],
                'iid_excess_kurtosis_at_one_attempt_per_face': float((fourth/second**2-3)/self.faces),
                'scope': 'Exact one-step moments and an independent-increment null prediction, not the actual multi-step process.'}

    def trajectory(self, initial, schedule_seed, duration, windows, check_charge_events=False):
        rng = random.Random(schedule_seed)
        n = self.faces*duration
        schedule = [rng.randrange(self.x.operator_count) for _ in range(n)]
        output = self.x.compiled_output([initial], schedule, n, raw_events=False)
        row = output['runs'][1]
        assert row['mode'] == 'combined' and row['condition'] == 0
        events = np.array(row['events'], dtype=np.int64).reshape(-1, 5)
        currents = np.zeros((len(events), 2), dtype=np.int64)
        q = np.array(self.x.p.charge(initial))
        for arity in (0, 1):
            chosen = np.flatnonzero((events[:, 1] >= 3) == bool(arity))
            ticks, rules, patches, before, after = events[chosen].T
            first = 36 if arity else 6
            j0 = self.charges[before//first]-self.charges[after//first]
            currents[chosen] = self.steps[arity][patches, 0]*j0[:, None]
            if arity:
                j1 = self.charges[after % 6]-self.charges[before % 6]
                currents[chosen] += self.steps[arity][patches, 1]*j1[:, None]
        if check_charge_events:
            for _, rule, patch, before, after in events:
                faces = self.x.p.supports[int(rule >= 3)][patch]
                old = [self.charges[(before//6**k) % 6] for k in reversed(range(len(faces)))]
                new = [self.charges[(after//6**k) % 6] for k in reversed(range(len(faces)))]
                assert q[list(faces)].tolist() == old
                q[list(faces)] = new
            assert q.tolist() == self.x.p.charge(row['final_links'])
        initial_q = np.array(self.x.p.charge(initial))
        final_q = np.array(self.x.p.charge(row['final_links']))
        assert initial_q.sum() == final_q.sum() == self.faces
        # A winding cochain and bounded polarization agree modulo torus periods.
        assert np.all((currents.sum(axis=0)-(final_q-initial_q)@self.coordinates) % (3*self.side) == 0)
        unit = np.zeros((duration, 2), dtype=np.int64)
        np.add.at(unit, (events[:, 0]-1)//self.faces, currents)
        stats = []
        for width in windows:
            assert duration % width == 0
            block = unit.reshape(duration//width, width, 2).sum(axis=1)/3.0
            a, b = block.T
            divisor = 2*self.faces*width
            stats.append([float(np.mean(a*a)/divisor), float(np.mean(b*b)/divisor),
                          float(np.mean(a*b)/divisor), float(np.mean(a*a)), float(np.mean(a**4))])
        return {'window_statistics': stats, 'total_current': (currents.sum(axis=0)/3.0).tolist(),
                'changing_events': len(events)}


def measure(side=6, trials=128, duration=4096, seed=184651001,
            windows=(1, 4, 16, 64, 256, 1024), workers=4):
    started = time.monotonic()
    observer = Winding(side)
    reference = AxialReference(side, observer.bank, observer.faces)
    initials = [reference.sample(random.Random(seed+2*i), 'nonabelian_reflections')['links']
                for i in range(trials)]
    run = lambda i: observer.trajectory(initials[i], seed+2*i+1, duration, windows, i == 0)
    results = []
    with ThreadPoolExecutor(max_workers=workers) as pool:
        for result in pool.map(run, range(trials)):
            results.append(result)
            if len(results) % 8 == 0 or len(results) == trials:
                print(json.dumps({'completed_independent_trajectories': len(results),
                                  'elapsed_seconds': time.monotonic()-started}), flush=True)
    samples = np.array([row['window_statistics'] for row in results])
    mean = samples.mean(axis=0)
    se = samples.std(axis=0, ddof=1)/np.sqrt(trials) if trials > 1 else np.full_like(mean, np.nan)
    rows = []
    for i, width in enumerate(windows):
        q2, q4 = mean[i, 3:5]
        influence = ((samples[:, i, 4]-q4)/q2**2
                     -2*q4*(samples[:, i, 3]-q2)/q2**3)
        rows.append({'attempts_per_face': width, 'windows_per_trajectory': duration//width,
                     'sigma_xx': mean[i, 0], 'sigma_xx_cluster_SE': se[i, 0],
                     'sigma_yy': mean[i, 1], 'sigma_yy_cluster_SE': se[i, 1],
                     'sigma_xy': mean[i, 2], 'sigma_xy_cluster_SE': se[i, 2],
                     'x_excess_kurtosis': q4/q2**2-3,
                     'kurtosis_delta_method_cluster_SE': float(influence.std(ddof=1)/np.sqrt(trials)) if trials > 1 else None})
    total = np.array([row['total_current'] for row in results])
    return {'side': side, 'faces': observer.faces, 'independent_trajectories': trials,
            'duration_attempts_per_face': duration, 'total_forward_attempts': trials*duration*observer.faces,
            'initial_seed_base': seed, 'seed_rule': 'initial seed+2*i, schedule seed+2*i+1',
            'windows': rows, 'mean_total_current': total.mean(axis=0).tolist(),
            'mean_total_current_SE': (total.std(axis=0, ddof=1)/np.sqrt(trials)).tolist() if trials > 1 else None,
            'total_changing_events': sum(row['changing_events'] for row in results),
            'elapsed_seconds': time.monotonic()-started,
            'scope': 'Finite-time stationary winding fluctuations under the unchanged microscopic bank. Independent trajectory clusters, empirical errors; not proof of convergence, a hydrodynamic limit, or molecular/quantum behavior.'}


if __name__ == '__main__':
    print(json.dumps(measure()), flush=True)
