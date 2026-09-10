#!/usr/bin/env python3
"""Release the exact activation planes; no reference or rule is frozen."""
import json
from pathlib import Path

import numpy as np
from flint import fmpz_mat

from fiber_constraint_dynamics import canonical
from fiber_constraint_intersection import angular_rate_drift
from fiber_refinement_limit import kind, read
from fiber_two_constraint_witness import rates
from fiber_two_constraints import AffineProcess, parameterize


TIMES = (0., .01, .05, .2, .5, 1., 2.)
ALPHAS = (0, 1, 8)
COLUMNS = ['reaction_rate', 'forward_rate', 'backward_rate', 'identity_count',
           'reaction_count', 'joint_support_edges', 'zero_offset_rank',
           'angular_rate_drift', 'rate_inside_initial_radius_two', 'rate_outside_initial_radius_two',
           'active_fans', 'flat_RR_pairs']


def prepare(source):
    state = AffineProcess('contact', source['side'], source['seed'])
    rows, signs = source['before_angular_constraints'], source['before_angular_signs']
    state.rank = 4
    state.initial_rows = rows
    state.initial_signs = signs
    state.raw = list(zip(signs, parameterize(rows)))
    return state


def distances(state, root):
    incident = [[] for _ in state.oracle.edges]
    for f, edges in enumerate(state.fe):
        for e in edges:
            incident[e].append(f)
    neighbors = [[] for _ in range(state.F)]
    for a, b in incident:
        neighbors[a].append(b); neighbors[b].append(a)
    distance, queue = {root: 0}, [root]
    for f in queue:
        for g in neighbors[f]:
            if g not in distance:
                distance[g] = distance[f]+1; queue.append(g)
    return distance


def observe(state, distance):
    intensity = rates(state)
    q = state.charges()
    k = q.count(0)
    assert k == q.count(2) and k <= 4 and sum(q) == state.F
    backwards = sum(rate for rate, fs in zip(intensity, state.supports) if any(q[int(f)] == 0 for f in fs))
    B = fmpz_mat([list(h[1][:-1]) for h in state.raw])
    kernel, rank = B.transpose().nullspace()
    assert rank == 4
    rows = [canonical([int(kernel[e, i]) for e in range(state.E)])[0] for i in range(rank)]
    offsets = [sum(x*h[1][-1] for x, h in zip(row, state.raw)) for row in rows]
    inside = sum(rate for rate, fs in zip(intensity, state.supports) if distance[int(fs[1])] <= 2)
    flat = set()
    for p, fs in enumerate(state.supports):
        if q[int(fs[0])] == q[int(fs[1])] == 1:
            word = state.word(p)
            if word[0][1] == word[1][1]:
                flat.add(tuple(sorted(map(int, fs[:2]))))
    return [sum(intensity), sum(intensity)-backwards, backwards, k, len(state.events),
            sum(any(row[e] for row in rows) for e in range(state.E)), rank-int(any(offsets)),
            angular_rate_drift(state, intensity)[0], inside, sum(intensity)-inside,
            sum(rate > 0 for rate in intensity), len(flat)]


def experiment(source, trials):
    observations = [[] for _ in ALPHAS]
    histories = [[] for _ in ALPHAS]
    mismatch_counts, proposal_counts = [0]*len(ALPHAS), []
    for trial in range(trials):
        states = [prepare(source) for _ in ALPHAS]
        distance = distances(states[0], source['flat_star_hub'])
        rng = np.random.default_rng(915101001+trial)
        modulus = 1_000_000_007
        parameters = rng.integers(modulus, size=states[0].E-4)
        raw = np.array([sum(c*int(z) for c, z in zip(h[1][:-1], parameters)) % modulus
                        +(modulus if h[0] == -1 else 0) for h in states[0].raw], dtype=np.int64)
        shadows = [raw.copy() for _ in ALPHAS]
        rows = [[] for _ in ALPHAS]
        clock = np.random.default_rng(915201001+trial)
        total_rate = len(states[0].pe)*31
        t, tick, proposals = 0., 0, 0
        while tick < len(TIMES):
            t += clock.exponential(1/total_rate)
            while tick < len(TIMES) and TIMES[tick] <= t:
                for i, state in enumerate(states):
                    rows[i].append(observe(state, distance))
                    if shadows[i] is not None:
                        q = [kind(read(shadows[i], e, r, modulus), modulus) for e, r in zip(state.fe, state.fr)]
                        if q != state.charges():
                            mismatch_counts[i] += 1; shadows[i] = None
                tick += 1
            if tick == len(TIMES):
                break
            p, slot = int(clock.integers(len(states[0].pe))), int(clock.integers(31))
            proposals += 1
            for i, (alpha, state) in enumerate(zip(ALPHAS, states)):
                if slot >= 15 and (slot-15) % 8 >= alpha:
                    continue
                channel = slot if slot < 15 else 15+(slot-15)//8
                if not state.step(p, channel, t, shadows[i], parameters, modulus) and shadows[i] is not None:
                    mismatch_counts[i] += 1; shadows[i] = None
        for i in range(len(ALPHAS)):
            observations[i].append(rows[i]); histories[i].append(states[i].events)
        proposal_counts.append(proposals)
        if (trial+1) % 8 == 0:
            print(json.dumps({'side': source['side'], 'birth_distance': source['birth_distance'],
                              'finished_trials': trial+1, 'trials': trials}), flush=True)
    output = []
    for i, alpha in enumerate(ALPHAS):
        x = np.array(observations[i], dtype=float)
        output.append({'side': source['side'], 'faces': states[0].F, 'edges': states[0].E,
                       'birth_distance': source['birth_distance'], 'alpha': alpha, 'trials': trials,
                       'initial_sign_seed': source['seed'], 'clock_seed': 915201001,
                       'times': TIMES, 'columns': COLUMNS, 'observations': observations[i],
                       'mean': x.mean(axis=0).tolist(), 'SE': (x.std(axis=0, ddof=1)/np.sqrt(trials)).tolist(),
                       'probability_angular_enhancement': (x[:, :, 7] > 0).mean(axis=0).tolist(),
                       'probability_reaction_active': (x[:, :, 0] > 0).mean(axis=0).tolist(),
                       'histories': histories[i], 'shadow_mismatches': mismatch_counts[i],
                       'maximal_clock_proposal_counts': proposal_counts,
                       'initial_radius_two_faces': sum(d <= 2 for d in distance.values())})
    return output


if __name__ == '__main__':
    sources = json.loads(Path('data/fiber-angular-activation.json').read_text())['runs']
    result = {'scope': 'Released fixed rank-four post-reaction planes with coupled original Poisson schedules', 'runs': []}
    for side, birth, trials in ((4, 1, 32), (4, 2, 32), (8, 2, 16)):
        source = next(row for row in sources if row['side'] == side and row['birth_distance'] == birth)
        result['runs'].extend(experiment(source, trials))
        Path('data/fiber-released-activation.json').write_text(json.dumps(result, indent=2)+'\n')
        for run in result['runs'][-3:]:
            print(json.dumps({k: run[k] for k in ('side', 'birth_distance', 'alpha', 'trials',
                                                  'probability_angular_enhancement', 'probability_reaction_active', 'shadow_mismatches')}), flush=True)
