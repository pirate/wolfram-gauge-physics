#!/usr/bin/env python3
"""Nonuniform density relaxation through full raw-link fiber updates.

The heat equation is a prediction, never an update rule. A simultaneous
Hurwitz-clock exclusion process is retained only as a coupled observable.
"""
import json
from pathlib import Path

import numpy as np

from fiber_refinement_limit import evolve, geometry, read


def dual_geometry(oracle, side):
    incident = [[] for _ in oracle.edges]
    for f, path in enumerate(oracle.face_paths):
        for e, _ in path:
            incident[e].append(f)
    adjacent = [[] for _ in oracle.faces]
    for e, (a, b) in enumerate(incident):
        adjacent[a].append((b, e))
        adjacent[b].append((a, e))
    parent = np.full(len(adjacent), -1, dtype=np.int64)
    parent_edge = parent.copy()
    parent[0] = 0
    order = [0]
    for a in order:
        for b, e in adjacent[a]:
            if parent[b] < 0:
                parent[b], parent_edge[b] = a, e
                order.append(b)
    assert len(order) == len(adjacent)
    k = 2*np.pi/side
    phi = 2+np.exp(1j*k)
    branch = np.conjugate(phi)/abs(phi)
    mode = np.array([np.real(np.exp(1j*k*x)*v)
                     for y in range(side) for x in range(side) for v in (1, branch)])
    eigenvalue = 3-abs(phi)
    lap_mode = np.array([3*mode[a]-sum(mode[b] for b, _ in neighbors)
                         for a, neighbors in enumerate(adjacent)])
    residual = float(np.max(np.abs(lap_mode-eigenvalue*mode)))
    assert residual < 1e-12
    assert abs(mode.mean()) < 1e-12
    # Two rooted pair descriptions, each with H and H inverse at rate one.
    pair_counts = {}
    _, args = geometry(side)
    for a, b, _ in args[-1]:
        pair = tuple(sorted((int(a), int(b))))
        pair_counts[pair] = pair_counts.get(pair, 0)+1
    assert len(pair_counts) == len(incident) and set(pair_counts.values()) == {2}
    return mode, eigenvalue, residual, (parent, parent_edge, order)


def conditioned_means(probability):
    factors = 1-2*probability
    product = np.prod(factors)
    return np.array([p*(1-np.prod(np.delete(factors, i)))/(1+product)
                     for i, p in enumerate(probability)])


def prepare(rng, probability, n, oracle, args, tree):
    # On this even-face closed mesh the reflection count, hence also the
    # even-holonomy count, must be even. Rejection gives the exact law.
    eta = rng.random(len(probability)) < probability
    while int(eta.sum()) % 2:
        eta = rng.random(len(probability)) < probability
    signs = rng.integers(2, size=len(oracle.edges), dtype=np.int64)
    fe = args[2]
    mismatch = np.bitwise_xor.reduce(signs[fe], axis=1) ^ (1-eta.astype(np.int64))
    parent, parent_edge, order = tree
    for f in reversed(order[1:]):
        if mismatch[f]:
            signs[parent_edge[f]] ^= 1
            mismatch[f] ^= 1
            mismatch[parent[f]] ^= 1
    assert not np.any(mismatch)
    # Non-tree signs stay independent fair bits; the remaining signs are
    # the unique solution. Shifts are uniform conditional on all signs.
    raw = rng.integers(n, size=len(signs), dtype=np.int64)+n*signs
    return raw, eta


def experiment(side, alpha, trials, n=1_000_000_007, tau=0.02, seed=911306001):
    oracle, args = geometry(side)
    mode, eigenvalue, residual, tree = dual_geometry(oracle, side)
    probability = .5+.2*mode
    initial_mean = float(mode @ (conditioned_means(probability)-.5)/len(mode))
    predicted = initial_mean*np.exp(-4*eigenvalue*side**2*tau)
    rng = np.random.default_rng(seed)
    rows, sample = [], None
    for i in range(trials):
        raw, eta = prepare(rng, probability, n, oracle, args, tree)
        initial_raw = raw.copy() if i == 0 else None
        q, shadow, matched, identities, reactions, vacancies, raw = evolve(
            raw, n, alpha, side**2*tau, seed+i, *args)
        rows.append([mode @ (eta-.5)/len(mode), mode @ (q-1.5)/len(mode),
                     mode @ (shadow-.5)/len(mode), int(not matched), identities,
                     reactions, vacancies, np.abs(q-1-shadow).mean()])
        if i == 0 and side == 8:
            sample = {'initial_raw_links': initial_raw.tolist(), 'final_raw_links': raw.tolist(),
                      'initial_even_holonomy': eta.astype(int).tolist(),
                      'final_charge': q.tolist(), 'final_coupled_exclusion': shadow.tolist()}
    x = np.array(rows)
    means, se = x.mean(axis=0), x.std(axis=0, ddof=1)/np.sqrt(trials)
    F, T = len(mode), side**2*tau
    result = {'side': side, 'faces': F, 'cycle_vertices': n, 'angular_rate': alpha,
              'trials': trials, 'seed': seed, 'diffusive_time': tau, 'microscopic_time': T,
              'mode_eigenvalue': eigenvalue, 'mode_eigenvector_residual': residual,
              'scaled_decay_rate': 4*eigenvalue*side**2,
              'continuum_decay_rate': 16*np.pi**2/3,
              'exact_initial_exclusion_mode_mean': initial_mean,
              'predicted_final_exclusion_mode_mean': float(predicted),
              'observed_raw_mode_minus_prediction_in_SE': float((means[1]-predicted)/se[1]),
              'path_mismatch_bound': min(1., F/n+72*F*T*(n-1)/n**2),
              'columns': ['initial_exclusion_mode', 'final_raw_charge_mode', 'final_exclusion_mode',
                          'path_mismatch', 'initial_identity_count', 'reaction_events',
                          'vacancy_events', 'final_L1_error_per_face'],
              'means': means.tolist(), 'independent_trajectory_SE': se.tolist()}
    if sample is not None:
        result['sample_raw_trajectory_endpoints'] = sample
    return result


if __name__ == '__main__':
    result = {'model': 'Full raw-link process; no continuum law in updates',
              'clock': 'Each old channel rate one; each angular channel rate alpha',
              'runs': []}
    for side, alpha, trials in ((4, 1, 2048), (8, 1, 1024), (16, 1, 512), (8, 8, 1024)):
        row = experiment(side, alpha, trials)
        result['runs'].append(row)
        Path('data/fiber-refinement-heat.json').write_text(json.dumps(result, indent=2)+'\n')
        print(json.dumps({k: v for k, v in row.items() if k != 'sample_raw_trajectory_endpoints'}), flush=True)
