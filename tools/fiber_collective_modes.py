#!/usr/bin/env python3
"""Full three-face reaction/transport/angular component and its collective modes.

Generators are assembled from the actual primitive maps. Sector rates are
conditional stationary flows, not assumed Markov closure or a field equation.
"""
import itertools
import json
import math
import random
import statistics
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.sparse.linalg import eigsh, spsolve

from fiber_angular_transport import AngularFiber, AngularLastPair, holonomy_image, seed_pair
from cycle_relational_dynamics import TransportRule, family


class LastTransport:
    def __init__(self, g, k):
        self.base = TransportRule(g, k)

    def apply(self, w):
        b, c, _ = self.base.apply((w[1], w[2], w[0]))
        return w[0], b, c


def component(f):
    g, n = f.g, f.n
    p = f.reflections[0]
    reflections, rotations = f.reflections, f.rotations
    states = [(a, b, g.mul[g.inv[g.mul[a][b]]][p]) for a in reflections for b in reflections]
    states.remove((p, p, p))
    for z in rotations:
        for layout in itertools.permutations(('E', 'R', 'Z')):
            r = g.mul[g.inv[z]][p] if layout.index('Z') < layout.index('R') else g.mul[p][g.inv[z]]
            values = {'E': g.identity, 'R': r, 'Z': z}
            states.append(tuple(values[k] for k in layout))
    assert len(set(states)) == len(states) == (n-1)*(n+7)
    ids = {w: i for i, w in enumerate(states)}
    rotation_labels = {z: k for k, z in enumerate(f.powers[f.steps[0]])}
    reflection_labels = {g.mul[z][p]: k for z, k in rotation_labels.items()}
    sectors = []
    for w in states:
        if all(h in reflection_labels for h in w):
            divisor = math.gcd(n, reflection_labels[w[0]], reflection_labels[w[1]])
        else:
            z = next(h for h in w if h in f.distance)
            divisor = math.gcd(n, rotation_labels[z])
        sectors.append(n//divisor)
    return p, states, ids, np.array(sectors, dtype=int)


def raw_first_exit(f, states, ids, sectors, q, expected, trials=64):
    g = f.g
    oracle, prepared, _ = seed_pair(f, f.powers[f.steps[0]][f.n//q])
    paths = oracle.patches[0][::-1]
    rules = ([TransportRule(g, k) for k in range(3)]+[LastTransport(g, k) for k in range(3)]+family(g)
             +f.rules+[AngularLastPair(rule) for rule in f.rules])
    assert len(rules) == f.n+17
    starts = np.flatnonzero(sectors == q).tolist()
    code = lambda w: (w[0]*g.n+w[1])*g.n+w[2]
    original = tuple(oracle.transport(prepared, path) for path in paths)
    times, example = [], None
    for trial in range(trials):
        rng = random.Random(252222000+1000*f.n+trial)
        word = states[rng.choice(starts)]
        raw = prepared[:]
        oracle.update(raw, 0, {code(original): code(word)})  # Uniform conditioned initial state, not evolution.
        assert len(holonomy_image(oracle, raw)) == 2*q
        initial = raw[:]
        for tick in range(1, 10000001):
            rule = rules[rng.randrange(len(rules))]
            word = tuple(oracle.transport(raw, path) for path in paths)
            target = rule.apply(word)
            if target == word:
                continue
            oracle.update(raw, 0, {code(word): code(target)})
            if sectors[ids[target]] != q:
                assert len(holonomy_image(oracle, raw)) == 2*f.n
                times.append(tick/len(rules))
                if example is None:
                    example = {'initial_links': initial, 'final_links': raw, 'exit_proposal': tick,
                               'before_exit_tuple': word, 'after_exit_tuple': target}
                break
        else:
            raise ValueError('First-exit observation reached its bound')
    return {'initial_rotation_subgroup_order': q, 'independent_trajectories': trials,
            'initial_measure': 'Uniform framed states in the selected fixed-boundary Q=3 subgroup sector.',
            'clock': 'Discrete proposals divided by n+17; expectation equals rate-one-per-channel backward equation.',
            'mean_first_exit_time': statistics.mean(times), 'standard_error': statistics.stdev(times)/trials**0.5,
            'backward_equation_prediction': expected, 'trajectory_exit_times': times, 'raw_example': example}


def calculation(n):
    f = AngularFiber(n)
    g, m = f.g, f.m
    p, states, ids, sectors = component(f)
    N = len(states)
    legacy = [TransportRule(g, k) for k in range(3)]+[LastTransport(g, k) for k in range(3)]+family(g)
    rows, columns, values = [], [], []
    old_rows, old_columns = [], []
    cross = Counter()
    angular_exit = np.zeros(N)
    activity = Counter()
    for i, w in enumerate(states):
        # Only four elastic channels act on a generic RRR tuple. Excluding
        # analytically inactive slots changes no generator or proposal clock.
        is_rrr = all(h in f.reflections for h in w)
        candidates = legacy if not is_rrr or ((w[0] == w[1]) != (w[1] == w[2])) else [legacy[k] for k in (1, 2, 4, 5)]
        for rule in candidates:
            j = ids[rule.apply(w)]
            if i != j:
                assert sectors[i] == sectors[j]
                rows.extend((i, i)); columns.extend((i, j)); values.extend((1, -1))
                old_rows.append(i); old_columns.append(j)
                activity['legacy_changing_slots'] += 1
        if is_rrr:
            continue
        for position in (0, 1):
            a, b = w[position:position+2]
            z = a if a in f.distance and b in f.reflections else b if b in f.distance and a in f.reflections else None
            if z is None:
                continue
            d = f.distance[z]
            gate_indices = {d-1}
            if d > 1:
                gate_indices.add(d-2)
            for index in gate_indices:
                rule = f.rules[index] if position == 0 else AngularLastPair(f.rules[index])
                j = ids[rule.apply(w)]
                assert i != j
                rows.extend((i, i)); columns.extend((i, j)); values.extend((1, -1))
                activity['angular_changing_slots'] += 1
                if sectors[i] != sectors[j]:
                    cross[int(sectors[i]), int(sectors[j])] += 1
                    angular_exit[i] += 1
    L = coo_matrix((values, (rows, columns)), shape=(N, N), dtype=float).tocsr()
    assert (L-L.T).nnz == 0 and np.max(np.abs(np.asarray(L.sum(axis=1)))) == 0
    count, _ = connected_components(L, directed=False)
    assert count == 1
    old = coo_matrix((np.ones(len(old_rows)), (old_rows, old_columns)), shape=(N, N)).tocsr()
    old_count, old_labels = connected_components(old, directed=False)
    assert old_count == len(set(sectors))
    assert all(len(set(sectors[old_labels == j])) == 1 for j in range(old_count))

    tau = np.array([ids[tuple(g.conj[p][h] for h in w)] for w in states])
    assert np.all(tau != np.arange(N)) and np.array_equal(tau[tau], np.arange(N))
    assert (L[tau][:, tau]-L).nnz == 0
    representatives = [i for i in range(N) if i < tau[i]]
    quotient_ids = {rep: j for j, rep in enumerate(representatives)}
    qids = [quotient_ids[min(i, tau[i])] for i in range(N)]
    B = coo_matrix((np.ones(N), (np.arange(N), qids)), shape=(N, N//2)).tocsr()
    Q = (B.T@L@B)*0.5
    if N//2 <= 20:
        eigenvalues, eigenvectors = np.linalg.eigh(Q.toarray())
        eigenvalues, eigenvectors = eigenvalues[:4], eigenvectors[:, :4]
    else:
        eigenvalues, eigenvectors = eigsh(Q, k=4, which='SM', tol=1e-10,
                                         v0=np.linspace(0.1, 1, N//2), maxiter=100000)
        order = np.argsort(eigenvalues)
        eigenvalues, eigenvectors = eigenvalues[order], eigenvectors[:, order]
    assert abs(eigenvalues[0]) < 1e-8 and eigenvalues[1] > 0
    residual = max(np.linalg.norm(Q@eigenvectors[:, j]-eigenvalues[j]*eigenvectors[:, j]) for j in range(4))
    sizes = Counter(map(int, sectors))
    qs = sorted(sizes)
    flows = []
    coarse = np.zeros((len(qs), len(qs)))
    for i, q in enumerate(qs):
        row = {'rotation_subgroup_order': q, 'framed_states': sizes[q], 'stationary_probability': sizes[q]/N,
               'mean_exit_rate': float(np.mean(angular_exit[sectors == q])),
               'within_sector_exit_rate_variance': float(np.var(angular_exit[sectors == q]))}
        phi = sum(math.gcd(j, q) == 1 for j in range(1, q+1))
        jordan = sum(math.gcd(q, a, b) == 1 for a in range(q) for b in range(q))
        assert sizes[q] == jordan+6*phi
        if q < n:
            assert cross[q, n] == 8*phi and all(target == n for source, target in cross if source == q)
            row['exact_exit_rate_numerator_denominator'] = [8*phi, sizes[q]]
            row['indicator_rayleigh_bound'] = (8*phi/sizes[q])/(1-sizes[q]/N)
            members = np.flatnonzero(sectors[representatives] == q)
            killed = Q[members][:, members]
            waiting = spsolve(killed, np.ones(len(members)))
            second = 2*spsolve(killed, waiting)
            row['mean_first_exit_from_uniform_sector'] = float(np.mean(waiting))
            row['first_exit_standard_deviation'] = math.sqrt(float(np.mean(second))-float(np.mean(waiting))**2)
            row['inverse_stationary_exit_rate'] = sizes[q]/(8*phi)
            row['first_exit_backward_equation_residual'] = float(np.max(np.abs(killed@waiting-1)))
            row['poisson_trial_rayleigh_bound'] = float(np.mean(waiting))/(float(np.mean(waiting**2))-(sizes[q]/N)*float(np.mean(waiting))**2)
        flows.append(row)
        for j, r in enumerate(qs):
            if i != j:
                coarse[i, j] = -cross[q, r]/math.sqrt(sizes[q]*sizes[r])
                coarse[i, i] += cross[q, r]/sizes[q]
    coarse_values = np.linalg.eigvalsh(coarse)
    v = eigenvectors[:, 1]
    qsectors = sectors[representatives]
    overlap = sum(float(np.sum(v[qsectors == q]))**2/np.sum(qsectors == q) for q in qs)
    projected = np.zeros_like(v)
    for q in qs:
        projected[qsectors == q] = np.mean(v[qsectors == q])
    residual_by_type = Counter()
    for j, state_id in enumerate(representatives):
        w = states[state_id]
        if all(h in f.reflections for h in w):
            kind = 'reaction_enabled_RRR' if ((w[0] == w[1]) != (w[1] == w[2])) else 'other_RRR'
        else:
            kind = 'rotation_present'
        residual_by_type[kind] += float((v[j]-projected[j])**2)
    lifted = v[np.array(qids)]/math.sqrt(2)
    phase = np.array([any(h in f.distance for h in w) for w in states])
    energy = Counter()
    entries = L.tocoo()
    for i, j, weight in zip(entries.row, entries.col, entries.data):
        if i < j:
            kind = 'cross_subgroup' if sectors[i] != sectors[j] else 'reaction_within_subgroup' if phase[i] != phase[j] else 'other_within_subgroup'
            energy[kind] += float(-weight*(lifted[i]-lifted[j])**2)
    assert abs(sum(energy.values())-eigenvalues[1]) < 1e-9
    assert activity['angular_changing_slots'] == 8*(n-2)
    largest_proper = max(q for q in qs if q < n)
    observation = raw_first_exit(f, states, ids, sectors, largest_proper,
                                next(row['mean_first_exit_from_uniform_sector'] for row in flows if row['rotation_subgroup_order'] == largest_proper))
    return {'cycle_vertices': n, 'framed_component_states': N, 'gauge_quotient_states': N//2,
            'old_bank_components': old_count, 'augmented_components': count,
            'operator_slots': n+17, 'clock': 'One rate-one clock per channel; divide rates by n+17 for uniform proposal time.',
            'stationary_rotation_present_fraction': 6/(n+7), 'stationary_angular_activity': 8*(n-2)/N,
            'sector_flows': flows, 'directed_cross_sector_slots': [[a, b, c] for (a, b), c in sorted(cross.items())],
            'lowest_gauge_even_decay_rates': eigenvalues.tolist(), 'eigenpair_residual_norm': float(residual),
            'sector_averaged_decay_rates_not_a_closed_generator': coarse_values.tolist(),
            'slowest_mode_squared_projection_on_subgroup_sectors': float(overlap),
            'slow_mode_omitted_squared_norm_by_configuration': dict(residual_by_type),
            'slow_mode_dirichlet_energy_by_transition_type': dict(energy),
            'projected_slow_mode_rayleigh_quotient': float(projected@(Q@projected)/(projected@projected)),
            'raw_first_exit_observation': observation,
            'scope': 'Complete fixed-boundary Q=3 three-face component with both adjacent supports, not an extended spatial mesh or a continuum field.'}


if __name__ == '__main__':
    output = Path('data/fiber-collective-modes.json')
    report = {'scope': 'Full local bank spectral and sector-flow calculation; no physical Hamiltonian or Markovian sector closure assumed.', 'fibers': []}
    for n in (9, 27, 81, 243):
        row = calculation(n)
        report['fibers'].append(row)
        output.write_text(json.dumps(report, indent=2)+'\n')
        print(json.dumps({'cycle': n, 'states': row['framed_component_states'],
                          'rates': row['lowest_gauge_even_decay_rates'],
                          'coarse_rates': row['sector_averaged_decay_rates_not_a_closed_generator'],
                          'sector_overlap': row['slowest_mode_squared_projection_on_subgroup_sectors']}), flush=True)
