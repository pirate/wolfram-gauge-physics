#!/usr/bin/env python3
"""Released event-conditioned states and spontaneous full-link reaction logs.

Palm starts sample the state seen after a typical marked stationary
reaction. No links are frozen and no reaction or angular rate is changed.
Separate unconditioned trajectories observe spontaneous events directly.
"""
import json
from pathlib import Path

import numpy as np
from numba import njit

from fiber_refinement_limit import LAYOUTS, geometry, inv, kind, mul, primitive, read


@njit
def set_fan(raw, word, p, n, pe, pr):
    # Expose the three distinct rim links last. Each is a free bijective
    # coordinate for one based holonomy; all other links stay uniform.
    for t in range(3):
        first, last = raw[pe[p, t, 0]], raw[pe[p, t, 2]]
        if pr[p, t, 0]:
            first = inv(first, n)
        if pr[p, t, 2]:
            last = inv(last, n)
        rim = mul(mul(inv(last, n), word[2-t], n), inv(first, n), n)
        raw[pe[p, t, 1]] = inv(rim, n) if pr[p, t, 1] else rim
    for t in range(3):
        assert read(raw, pe[p, 2-t], pr[p, 2-t], n) == word[t]


@njit
def write_fan(raw, p, aa, bb, cc, n, pe, pr):
    first, rim0, rim1 = raw[pe[p, 0, 0]], raw[pe[p, 0, 1]], raw[pe[p, 1, 1]]
    if pr[p, 0, 0]:
        first = inv(first, n)
    if pr[p, 0, 1]:
        rim0 = inv(rim0, n)
    if pr[p, 1, 1]:
        rim1 = inv(rim1, n)
    prefix = mul(rim0, first, n)
    s2 = mul(prefix, inv(cc, n), n)
    s3 = mul(mul(rim1, prefix, n), inv(mul(bb, cc, n), n), n)
    raw[pe[p, 0, 2]] = s2 if pr[p, 0, 2] else inv(s2, n)
    raw[pe[p, 1, 2]] = s3 if pr[p, 1, 2] else inv(s3, n)


@njit
def activity(raw, p, n, pe, pr):
    a = read(raw, pe[p, 2], pr[p, 2], n)
    b = read(raw, pe[p, 1], pr[p, 1], n)
    c = read(raw, pe[p, 0], pr[p, 0], n)
    if a >= n and b >= n and c >= n:
        return 12 if (a == b) != (b == c) else 0
    qa, qb, qc = kind(a, n), kind(b, n), kind(c, n)
    return 4 if qa+qb+qc == 3 and qa != qb and qb != qc and qa != qc else 0


def palm_raw(rng, n, oracle, args, marked=0):
    raw = rng.integers(2*n, size=len(oracle.edges), dtype=np.int64)
    if rng.integers(2):
        r = int(rng.integers(n))
        s = (r+int(rng.integers(1, n))) % n
        word = (n+r, n+r, n+s) if rng.integers(2) else (n+s, n+r, n+r)
    else:
        values = (0, n+int(rng.integers(n)), int(rng.integers(1, n)))
        layout = LAYOUTS[rng.integers(6)]
        word = tuple(values[k] for k in layout)
    set_fan(raw, word, marked, n, args[0], args[1])
    return raw


def local_geometry(side):
    oracle, args = geometry(side)
    supports = args[-1]
    root = set(supports[0])
    selected = []
    for p, fs in enumerate(supports):
        if len(root & set(fs)) != 1:
            continue
        outside = [int(f) for f in fs if f not in root]
        A, B = ({e for e, _ in oracle.face_paths[f]} for f in outside)
        private = [sorted(A-oracle.reads[0]-B), sorted(B-oracle.reads[0]-A)]
        if all(private):
            selected.append((p, private))
    # A fixed geometrically observed target, not chosen by measured effect.
    target, private = next((p, e) for p, e in selected if p == 7)
    interfering = sum(bool(W & oracle.reads[target]) for W in oracle.writes)
    assert interfering == 19 and len(selected) == 12
    return oracle, args, target, {'marked_fan': 0, 'target_fan': target,
                               'marked_faces': supports[0].tolist(), 'target_faces': supports[target].tolist(),
                               'target_private_edges': private, 'interfering_roots': interfering,
                               'single_face_private_edge_targets': [p for p, _ in selected]}


@njit
def trajectory(raw, n, alpha, seed, times, target, pe, pr, fe, fr, supports):
    np.random.seed(seed)
    channels = 15+2*alpha
    rate = len(pe)*channels
    q = np.array([kind(read(raw, fe[f], fr[f], n), n) for f in range(len(fe))])
    initial_Q = q.sum()
    output = np.zeros((len(times), 11))
    counts = np.zeros(4, dtype=np.int64)
    tick, t = 0, 0.
    logs = []
    while tick < len(times):
        next_t = t+np.random.exponential(1/rate)
        while tick < len(times) and times[tick] <= next_t:
            total = 0
            for p in range(len(pe)):
                total += activity(raw, p, n, pe, pr)
            output[tick, :4] = (activity(raw, 0, n, pe, pr), activity(raw, target, n, pe, pr), total, np.sum(q == 0))
            output[tick, 4:8] = counts
            output[tick, 8:] = counts[1:] > 0
            tick += 1
        if tick == len(times):
            break
        t = next_t
        slot = np.random.randint(len(pe)*channels)
        p, channel = slot//channels, slot % channels
        if channel >= 15:
            channel = 15+(channel-15)//alpha
        a = read(raw, pe[p, 2], pr[p, 2], n)
        b = read(raw, pe[p, 1], pr[p, 1], n)
        c = read(raw, pe[p, 0], pr[p, 0], n)
        aa, bb, cc = primitive(a, b, c, channel, n)
        if (a, b, c) == (aa, bb, cc):
            continue
        write_fan(raw, p, aa, bb, cc, n, pe, pr)
        oldq = (kind(a, n), kind(b, n), kind(c, n))
        newq = (kind(aa, n), kind(bb, n), kind(cc, n))
        assert sum(oldq) == sum(newq)
        for i in range(3):
            q[supports[p, i]] = newq[i]
        if 3 <= channel < 15:
            exterior = False
            for i in range(3):
                outside = True
                for f in supports[0]:
                    if supports[p, i] == f:
                        outside = False
                if outside and newq[i] != oldq[i]:
                    exterior = True
            counts += np.array((1, int(p == 0), int(p == target), int(exterior)))
            logs.append((t, float(p), float(channel), float(int(oldq[0] == oldq[1] == oldq[2] == 1))))
    actual = np.array([kind(read(raw, fe[f], fr[f], n), n) for f in range(len(fe))])
    assert np.all(q == actual) and q.sum() == initial_Q
    return output, logs


def summarize(x):
    return {'means': x.mean(axis=0).tolist(),
            'independent_trajectory_SE': (x.std(axis=0, ddof=1)/np.sqrt(len(x))).tolist()}


def released_experiment(n, alpha, trials=2048, side=4, seed=911701001):
    oracle, args, target, geo = local_geometry(side)
    times = np.array((0., .002, .01, .05, .2, 1.))
    rng = np.random.default_rng(seed)
    rows = []
    for i in range(trials):
        raw = palm_raw(rng, n, oracle, args)
        observed, _ = trajectory(raw, n, alpha, seed+i, times, target, *args)
        rows.append(observed)
    return {'n': n, 'alpha': alpha, 'side': side, 'faces': len(oracle.faces), 'trials': trials,
            'seed': seed, 'times': times.tolist(), 'geometry': geo,
            'columns': ['marked_rate', 'neighbor_rate', 'total_rate', 'identity_faces',
                        'all_reactions', 'marked_reactions', 'neighbor_reactions', 'exterior_charge_reactions',
                        'any_marked_reaction', 'any_neighbor_reaction', 'any_exterior_charge_reaction'],
            'exact_stationary_rate_per_fan': 6*(n-1)/n**2,
            'exact_post_event_marked_rate': 8.,
            'exact_post_event_neighbor_rate': 1/3+16*(n-1)/(3*n**2),
            **summarize(np.array(rows))}


def window_counts(logs, duration, windows, supports):
    # Each trajectory, not each event, is the independent uncertainty unit.
    out = np.zeros((len(windows), 6), dtype=np.int64)
    for wi, delta in enumerate(windows):
        for i, event in enumerate(logs):
            t, p = event[0], int(event[1])
            if t+delta > duration:
                continue
            out[wi, 0] += 1
            parent = set(supports[p])
            any_event = same = single_overlap = False
            count = 0
            for after in logs[i+1:]:
                if after[0] > t+delta:
                    break
                any_event, count = True, count+1
                same |= int(after[1]) == p
                single_overlap |= len(parent & set(supports[int(after[1])])) == 1
            out[wi, 1:] += (any_event, same, single_overlap, count, int(event[3]))
    return out


def spontaneous_experiment(n, alpha=1, trials=32, side=4, seed=911905001):
    oracle, args, target, _ = local_geometry(side)
    duration = n/4
    windows = np.array((.01, .05, .2, 1.))
    rng = np.random.default_rng(seed)
    counts, all_windows, sample = [], [], None
    for i in range(trials):
        raw = rng.integers(2*n, size=len(oracle.edges), dtype=np.int64)
        initial = raw.copy() if i == 0 else None
        _, logs = trajectory(raw, n, alpha, seed+i, np.array((duration,)), target, *args)
        counts.append(len(logs))
        all_windows.append(window_counts(logs, duration, windows, args[-1]))
        if i == 0:
            sample = {'initial_raw_links': initial.tolist(), 'final_raw_links': raw.tolist(),
                      'reaction_log_columns': ['time', 'rooted_fan', 'channel', 'RRR_to_ERZ'],
                      'reactions': logs}
    x = np.array(all_windows, dtype=float)
    denominator = x[:, :, 0].mean(axis=0)
    ratios = x[:, :, 1:].mean(axis=0)/denominator[:, None]
    influence = (x[:, :, 1:]-x[:, :, :1]*ratios[None, ...])/denominator[None, :, None]
    lam = 18*len(oracle.faces)*(n-1)/n**2
    return {'n': n, 'alpha': alpha, 'side': side, 'faces': len(oracle.faces),
            'duration': duration, 'trials': trials, 'seed': seed, 'windows': windows.tolist(),
            'exact_stationary_total_rate': lam, 'expected_reactions_per_trajectory': lam*duration,
            'reaction_counts_by_trajectory': counts, 'window_counts_by_trajectory': x.astype(int).tolist(),
            'ratio_columns': ['any_reaction', 'same_marked_fan', 'single_face_overlap', 'mean_reaction_count', 'parent_forward_fraction'],
            'event_conditioned_ratios': ratios.tolist(),
            'independent_trajectory_ratio_SE': (influence.std(axis=0, ddof=1)/np.sqrt(trials)).tolist(),
            'independent_Poisson_any_event_prediction': (1-np.exp(-lam*windows)).tolist(),
            'sample_unconditioned_trajectory': sample}


if __name__ == '__main__':
    result = {'scope': 'Full raw-link process: all anchors released, original rates retained',
              'released_Palm_runs': [], 'spontaneous_runs': []}
    path = Path('data/fiber-reaction-bursts.json')
    for n, alpha in ((27, 1), (243, 1), (2187, 1), (1_000_000_007, 1), (1_000_000_007, 8)):
        row = released_experiment(n, alpha)
        result['released_Palm_runs'].append(row)
        path.write_text(json.dumps(result, indent=2)+'\n')
        print(json.dumps({'released_n': n, 'alpha': alpha, 'initial_rates': row['means'][0][:3],
                          'probabilities_at_0.2': row['means'][4][8:]}), flush=True)
    for n, alpha in ((243, 1), (729, 1), (2187, 1), (729, 8)):
        row = spontaneous_experiment(n, alpha)
        result['spontaneous_runs'].append(row)
        path.write_text(json.dumps(result, indent=2)+'\n')
        print(json.dumps({'spontaneous_n': n, 'alpha': alpha, 'reactions': sum(row['reaction_counts_by_trajectory']),
                          'window_0.05': row['event_conditioned_ratios'][1],
                          'SE': row['independent_trajectory_ratio_SE'][1],
                          'Poisson': row['independent_Poisson_any_event_prediction'][1]}), flush=True)
