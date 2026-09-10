#!/usr/bin/env python3
"""Stationary current memory on actual periodic shared-link geometries.

Uniform raw-link starts, unchanged local rules, Poisson channel clocks.
The compiled loop only accelerates group multiplication and link rewriting.
"""
import itertools
import json
import time
from pathlib import Path

import numpy as np
from numba import njit

from cycle_relational_dynamics import TransportRule, family
from fiber_angular_transport import AngularFiber
from run_three_face import LinkOracle


@njit
def transport(raw, edges, reverse, mul, inv):
    h = 0
    for i in range(3):
        v = raw[edges[i]]
        if reverse[i]:
            v = inv[v]
        h = mul[v, h]
    return h


@njit
def current_drift(q, supports, steps):
    out = np.zeros(2, dtype=np.int64)
    for p in range(len(supports)):
        a, b, c = q[supports[p, 0]], q[supports[p, 1]], q[supports[p, 2]]
        j0 = (2+int(a == 0 or b == 0))*(a-b)
        j1 = 0
        if a+b+c == 3 and a != b and b != c and a != c:
            j0 += 4*(a-1)
            j1 += 4*(1-c)
        for d in range(2):
            out[d] += steps[p, 0, d]*j0+steps[p, 1, d]*j1
    return out


@njit
def trajectory(raw, seed, alpha, dt, intervals, mul, inv, charge, tables,
               patch_edges, patch_reverse, face_edges, face_reverse, supports, steps):
    np.random.seed(seed)
    order, patches = len(inv), len(supports)
    channels = 15+alpha*(len(tables)-15)
    q = np.empty(len(face_edges), dtype=np.int64)
    for i in range(len(q)):
        q[i] = charge[transport(raw, face_edges[i], face_reverse[i], mul, inv)]
    winding = np.zeros((intervals+1, 2), dtype=np.int64)
    drift = np.zeros((intervals+1, 2), dtype=np.int64)
    drift[0] = current_drift(q, supports, steps)
    events = np.zeros(3, dtype=np.int64)
    proposals = 0
    for tick in range(intervals):
        winding[tick+1] = winding[tick]
        attempts = np.random.poisson(patches*channels*dt)
        proposals += attempts
        for _ in range(attempts):
            slot = np.random.randint(patches*channels)
            p, channel = slot//channels, slot % channels
            if channel >= 15:
                channel = 15+(channel-15)//alpha
            # Native fan order is opposite to ordered boundary multiplication.
            a = transport(raw, patch_edges[p, 2], patch_reverse[p, 2], mul, inv)
            b = transport(raw, patch_edges[p, 1], patch_reverse[p, 1], mul, inv)
            c = transport(raw, patch_edges[p, 0], patch_reverse[p, 0], mul, inv)
            old = (a*order+b)*order+c
            new = tables[channel, old]
            if old == new:
                continue
            aa, bb, cc = new//(order*order), new//order % order, new % order
            # Reconstruct the two internal spokes using unchanged rim prefixes.
            v0 = raw[patch_edges[p, 0, 0]]
            v1 = raw[patch_edges[p, 0, 1]]
            v2 = raw[patch_edges[p, 1, 1]]
            if patch_reverse[p, 0, 0]:
                v0 = inv[v0]
            if patch_reverse[p, 0, 1]:
                v1 = inv[v1]
            if patch_reverse[p, 1, 1]:
                v2 = inv[v2]
            prefix = mul[v1, v0]
            s2 = mul[prefix, inv[cc]]
            s3 = mul[mul[v2, prefix], inv[mul[bb, cc]]]
            raw[patch_edges[p, 0, 2]] = s2 if patch_reverse[p, 0, 2] else inv[s2]
            raw[patch_edges[p, 1, 2]] = s3 if patch_reverse[p, 1, 2] else inv[s3]
            dq0, dq2 = charge[aa]-charge[a], charge[cc]-charge[c]
            assert dq0+charge[bb]-charge[b]+dq2 == 0
            for d in range(2):
                winding[tick+1, d] += -steps[p, 0, d]*dq0+steps[p, 1, d]*dq2
            q[supports[p, 0]] = charge[aa]
            q[supports[p, 1]] = charge[bb]
            q[supports[p, 2]] = charge[cc]
            events[0 if channel < 3 else 1 if channel < 15 else 2] += 1
        drift[tick+1] = current_drift(q, supports, steps)
    actual_q = np.empty(len(q), dtype=np.int64)
    for i in range(len(q)):
        actual_q[i] = charge[transport(raw, face_edges[i], face_reverse[i], mul, inv)]
    assert np.all(actual_q == q)
    return raw, q, winding, drift, events, proposals


class Experiment:
    def __init__(self, n, side):
        self.n, self.side = n, side
        f = AngularFiber(n)
        self.f, self.g = f, f.g
        g, oracle = self.g, LinkOracle(side, self.g)
        self.oracle = oracle
        assert g.identity == 0
        rules = [TransportRule(g, i) for i in range(3)]+family(g)+f.rules
        self.tables = np.array([[((a*g.n+b)*g.n+c) for a, b, c in
                                 (rule.apply(w) for w in itertools.product(range(g.n), repeat=3))]
                                for rule in rules], dtype=np.int64)
        self.charge = np.array([0 if h == g.identity else 1 if h in f.reflections else 2 for h in range(g.n)], dtype=np.int64)
        path_array = np.array(oracle.patches, dtype=np.int64)
        face_array = np.array(oracle.face_paths, dtype=np.int64)
        ids = {frozenset(face[:3]): i for i, face in enumerate(oracle.faces)}
        self.supports = np.array([[ids[frozenset(face[:3])] for face in spec[::-1]] for spec in oracle.specs])
        # Exactly one rooted pair per three-fan anchor; no duplicate first/last slots.
        assert len({tuple(row[:2]) for row in self.supports}) == len(self.supports)
        self.coordinates = np.array([[3*(i//2 % side)+(2 if i % 2 == 0 else 1),
                                      3*(i//2//side)+(1 if i % 2 == 0 else 2)] for i in range(len(oracle.faces))])
        period = 3*side
        self.steps = np.array([[(self.coordinates[b]-self.coordinates[a]+period//2) % period-period//2
                               for a, b in zip(support, support[1:])] for support in self.supports])
        self.arguments = (np.array(g.mul), np.array(g.inv), self.charge, self.tables,
                          path_array[..., 0], path_array[..., 1], face_array[..., 0], face_array[..., 1],
                          self.supports, self.steps)
        old = np.arange(g.n**3)
        before = self.charge[np.column_stack((old//(g.n*g.n), old//g.n % g.n, old % g.n))]
        self.shot_numerator = np.zeros((2, 2), dtype=np.int64)
        for table in self.tables:
            after = self.charge[np.column_stack((table//(g.n*g.n), table//g.n % g.n, table % g.n))]
            j0, j1 = before[:, 0]-after[:, 0], after[:, 2]-before[:, 2]
            for steps in self.steps:
                j = j0[:, None]*steps[0]+j1[:, None]*steps[1]
                self.shot_numerator += j.T@j
        self.bare_mobility = self.shot_numerator/(18*len(oracle.faces)*g.n**3)

    def run(self, alpha, trials=64, dt=0.001, duration=8., seed=318472001):
        intervals = round(duration/dt)
        cutoffs = np.array([0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0])
        lag_ids = np.rint(cutoffs/dt).astype(int)
        maxlag = max(lag_ids)
        samples, winding_samples, event_counts, quadrature_differences = [], [], [], []
        records = []
        started = time.monotonic()
        faces = len(self.oracle.faces)
        for trial in range(trials):
            rng = np.random.default_rng(seed+2*trial)
            raw = rng.integers(self.g.n, size=len(self.oracle.edges), dtype=np.int64)
            initial_q = np.array([self.charge[self.oracle.transport(raw, p)] for p in self.oracle.face_paths])
            final, q, winding, drift, events, proposals = trajectory(raw.copy(), seed+2*trial+1, alpha, dt, intervals, *self.arguments)
            assert q.sum() == initial_q.sum()
            assert np.all((winding[-1]-(q-initial_q)@self.coordinates) % (3*self.side) == 0)
            # Ensemble mean current drift is exactly zero; do not subtract a time mean.
            length = 1 << (2*len(drift)-1).bit_length()
            spectrum = np.fft.rfft(drift, n=length, axis=0)
            correlation = np.fft.irfft(abs(spectrum)**2, n=length, axis=0)[:maxlag+1]
            correlation /= (len(drift)-np.arange(maxlag+1))[:, None]
            correlation /= 9*faces
            trapezoid = np.cumsum((correlation[:-1]+correlation[1:])*dt/2, axis=0)[lag_ids-1]
            assert np.all(lag_ids % 2 == 0)
            integrals = np.array([(correlation[0]+correlation[lag]+4*correlation[1:lag:2].sum(axis=0)
                                   +2*correlation[2:lag:2].sum(axis=0))*dt/3 for lag in lag_ids])
            samples.append(integrals)
            quadrature_differences.append(trapezoid-integrals)
            blockstats = []
            for width in (0.1, 0.5, 1., 2.):
                size = round(width/dt)
                blocks = winding[size::size]-winding[:-size:size]
                blockstats.append(np.mean(blocks*blocks, axis=0)/(18*faces*width))
            winding_samples.append(blockstats)
            event_counts.append(events)
            if trial == 0:
                records.append({'initial_links': raw.tolist(), 'final_links': final.tolist(),
                                'initial_total_charge': int(initial_q.sum()), 'final_total_charge': int(q.sum()),
                                'integrated_winding_times_three': winding[-1].tolist(), 'proposals': proposals})
            if (trial+1) % max(16, trials//8) == 0:
                print(json.dumps({'side': self.side, 'alpha': alpha, 'trajectories': trial+1,
                                  'elapsed_seconds': time.monotonic()-started}), flush=True)
        samples, ws = np.array(samples), np.array(winding_samples)
        return {'cycle_vertices': self.n, 'side': self.side, 'faces': faces, 'angular_rate': alpha,
                'initial_and_schedule_seed_base': seed,
                'independent_trajectories': trials, 'duration_per_channel_time': duration, 'sampling_dt': dt,
                'cutoffs': cutoffs.tolist(), 'bare_mobility_tensor': self.bare_mobility.tolist(),
                'integrated_current_memory_per_face': samples.mean(axis=0).tolist(),
                'quadrature': 'Composite Simpson on sampled correlations',
                'trapezoid_minus_Simpson': np.mean(quadrature_differences, axis=0).tolist(),
                'trajectory_cluster_SE': (samples.std(axis=0, ddof=1)/np.sqrt(trials)).tolist(),
                'trajectory_memory_samples': samples.tolist(),
                'winding_windows': [0.1, 0.5, 1., 2.], 'finite_window_mobility': ws.mean(axis=0).tolist(),
                'finite_window_mobility_cluster_SE': (ws.std(axis=0, ddof=1)/np.sqrt(trials)).tolist(),
                'changed_pair_reaction_angular_events': np.sum(event_counts, axis=0).tolist(),
                'raw_trajectory_records': records,
                'scope': 'Uniform independent raw connections: stationary mixture over conserved charges and components. Finite-cutoff memory integrals are not exact asymptotic mobility; sampling quadrature and tail errors are not in cluster SE.'}


if __name__ == '__main__':
    path = Path('data/fiber-bulk-transport-confirmation.json')
    output = {'design': 'Fresh stationary starts; matched per-channel clocks; sides 4,8,12 and angular rates 0,1,8. Same initial connection within trajectory pairs across rates.',
              'calculations': []}
    for side in (4, 8, 12):
        experiment = Experiment(5, side)
        for alpha in (0, 1, 8):
            result = experiment.run(alpha, trials=512, dt=0.001, duration=32., seed=419683001)
            output['calculations'].append(result)
            path.write_text(json.dumps(output, indent=2)+'\n')
            print(json.dumps({k: result[k] for k in ('side', 'angular_rate', 'integrated_current_memory_per_face', 'trajectory_cluster_SE')}), flush=True)
