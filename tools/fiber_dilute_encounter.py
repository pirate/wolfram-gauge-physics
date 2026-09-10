#!/usr/bin/env python3
"""Exact position generator during the one-constraint E/Z phase at Q=F.

Rates are counted from the original rooted channels. Translation orbits
reduce ordered face pairs, without altering their exclusion or killing.
"""
import json
import math
from collections import Counter
from fractions import Fraction
from pathlib import Path

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import eigsh, spsolve

from fiber_refinement_limit import geometry


def experiment(side):
    oracle, args = geometry(side)
    F = len(oracle.faces)
    incident = [[] for _ in oracle.edges]
    for f, path in enumerate(oracle.face_paths):
        for e, _ in path:
            incident[e].append(f)
    adjacent = [[] for _ in range(F)]
    for a, b in incident:
        adjacent[a].append(b); adjacent[b].append(a)
    fan_counts = Counter()
    for faces in args[-1]:
        for a in faces:
            for b in faces:
                if a != b:
                    fan_counts[int(a), int(b)] += 1

    def relative(a, b):
        ax, ay = a//2 % side, a//2//side
        bx, by = b//2 % side, b//2//side
        return a % 2, 2*(((by-ay) % side)*side+(bx-ax) % side)+b % 2

    states = [(a, b) for a in (0, 1) for b in range(F) if a != b]
    ids = {state: i for i, state in enumerate(states)}
    rows, cols, data, killing = [], [], [], []
    for i, (a, b) in enumerate(states):
        transitions = Counter()
        for other in adjacent[a]:
            if other != b:
                transitions[ids[relative(other, b)]] += 6
        for other in adjacent[b]:
            if other != a:
                transitions[ids[relative(a, other)]] += 4
        if b in adjacent[a]:
            transitions[ids[relative(b, a)]] += 6
        kappa = 4*fan_counts[a, b]
        distance_two = {c for v in adjacent[a] for c in adjacent[v]}-{a}
        assert kappa == (16 if b in adjacent[a] else 4 if b in distance_two else 0)
        for j, rate in transitions.items():
            rows.append(i); cols.append(j); data.append(-rate)
        rows.append(i); cols.append(i); data.append(sum(transitions.values())+kappa)
        killing.append(kappa)
    K = coo_matrix((data, (rows, cols)), shape=(len(states), len(states)), dtype=float).tocsr()
    assert (K-K.transpose()).nnz == 0
    kappa = np.array(killing, dtype=float)
    assert kappa.sum() == 144
    assert np.array_equal(K @ np.ones(len(states)), kappa)
    birth_counts = Counter()
    for (a, b), count in fan_counts.items():
        birth_counts[ids[relative(a, b)]] += count
    birth = np.array([birth_counts[i] for i in range(len(states))], dtype=float)
    birth /= birth.sum()
    assert np.array_equal(birth, kappa/144)
    mean_wait = spsolve(K, np.ones(len(states)))
    principal = eigsh(K, k=1, sigma=0, which='LM', return_eigenvectors=False)[0]
    fixed_parent_mean = np.mean([mean_wait[ids[relative(int(a), int(b))]]
                                 for a in args[-1][0] for b in args[-1][0] if a != b])
    return {'side': side, 'faces': F, 'full_ordered_pair_states': F*(F-1),
            'translation_orbit_states': len(states), 'E_R_exchange_rate': 6, 'Z_R_exchange_rate': 4,
            'E_Z_exchange_rate': 6, 'killing_rate_at_distance_one': 16, 'killing_rate_at_distance_two': 4,
            'sum_of_relative_state_killing_rates': int(kappa.sum()),
            'exact_birth_mean_wait': str(Fraction(F-1, 72)),
            'solved_birth_mean_wait': float(birth @ mean_wait),
            'fixed_parent_birth_mean_wait': float(fixed_parent_mean),
            'uniform_position_mean_wait': float(mean_wait.mean()),
            'uniform_mean_over_F_log_F': float(mean_wait.mean()/(F*math.log(F))),
            'principal_killing_rate': float(principal),
            'mean_equation_residual': float(np.max(np.abs(K @ mean_wait-1)))}


if __name__ == '__main__':
    result = {'scope': 'Exact stopped E/Z phase in the one-constraint, all-reflection-background Q=F limit', 'runs': []}
    for side in (4, 8, 16, 32, 64):
        row = experiment(side)
        result['runs'].append(row)
        Path('data/fiber-dilute-encounter.json').write_text(json.dumps(result, indent=2)+'\n')
        print(json.dumps(row), flush=True)
