#!/usr/bin/env python3
"""Exact charge projection and its unresolved-gauge memory operator.

Projection is conditional expectation in the actual stationary reference.
The derived memory is not replaced by a Markov update on charge fields.
"""
import json
import itertools
from collections import Counter

import numpy as np
from flint import fmpz_mat

from triangle_coupled_star import CoupledStar
from triangle_patch_observer import WordObserver


def calculation():
    s = CoupledStar()
    n, clock = len(s.keys), len(s.tables)
    patterns = sorted(set(s.charges))
    ids = {q: i for i, q in enumerate(patterns)}
    labels = np.array([ids[q] for q in s.charges])
    sizes = np.bincount(labels)
    membership = np.zeros((n, len(patterns)), dtype=np.int64)
    membership[np.arange(n), labels] = 1
    basis = membership/np.sqrt(sizes)
    assert np.allclose(basis.T@basis, np.eye(len(patterns)))
    outgoing = s.counts@membership
    fiber_sums = membership.T@outgoing
    residual_integer = outgoing*sizes[labels, None]-fiber_sums[labels]
    distinct = sorted(set(map(tuple, residual_integer.tolist())))
    rank = fmpz_mat([list(row) for row in distinct]).rank()
    assert not np.any(residual_integer@np.array(patterns, dtype=np.int64))

    # Factor the hidden coupling through actual reaction activities. Three
    # dependencies in the output map come from disjoint reaction diamonds.
    fans = sorted({p for p, rule in s.channels if rule >= 3})
    face_ids = {f: i for i, f in enumerate(s.faces)}
    supports = {p: [face_ids[f] for f in s.x.p.supports[1][p]] for p in fans}
    ports = [(j, p) for j, q in enumerate(patterns) for p in fans
             if all(q[f] == 1 for f in supports[p])]
    port_ids = {v: i for i, v in enumerate(ports)}
    raw_activity = np.zeros((n, len(ports)), dtype=np.int64)
    for i, raw in enumerate(s.raw):
        for p in fans:
            key = (int(labels[i]), p)
            if key not in port_ids:
                continue
            code = s.x.e.factor.oracle.code(raw, 1, p)
            a, b, c = code//36, code//6 % 6, code % 6
            raw_activity[i, port_ids[key]] = int((a == b) != (b == c))
    centered = raw_activity*sizes[labels, None]-(membership.T@raw_activity)[labels]
    def charge_targets(j, p):
        for values in itertools.permutations((0, 1, 2)):
            q = list(patterns[j])
            for f, a in zip(supports[p], values):
                q[f] = a
            yield ids[tuple(q)]
    output = np.zeros((len(ports), len(patterns)), dtype=np.int64)
    for k, (j, p) in enumerate(ports):
        output[k, j] -= 12
        for target in charge_targets(j, p):
            output[k, target] += 2
    assert np.array_equal(centered@output, residual_integer)
    activity_rank = fmpz_mat([list(row) for row in sorted(set(map(tuple, centered.tolist())))]).rank()
    output_rank = fmpz_mat(output.tolist()).rank()
    diamonds = []
    flat = ids[(1,)*6]
    for p, h in itertools.combinations(fans, 2):
        if set(supports[p]) & set(supports[h]):
            continue
        oracle = s.x.e.factor.oracle.fan
        assert not (oracle.writes[p] & oracle.reads[h] or oracle.writes[h] & oracle.reads[p])
        relation = np.zeros(len(ports), dtype=np.int64)
        relation[port_ids[flat, p]], relation[port_ids[flat, h]] = 6, -6
        for j in charge_targets(flat, p):
            relation[port_ids[j, h]] += 1
        for j in charge_targets(flat, h):
            relation[port_ids[j, p]] -= 1
        assert not np.any(relation@output)
        diamonds.append({'fan_pair': [p, h], 'relation_coefficients': relation.tolist()})
    assert fmpz_mat([d['relation_coefficients'] for d in diamonds]).rank() == len(ports)-output_rank
    assert activity_rank == len(ports) == 42 and output_rank == rank == 39
    dropped = [port_ids[flat, d['fan_pair'][0]] for d in diamonds]
    kept = [k for k in range(len(ports)) if k not in dropped]
    expansion6 = np.zeros((len(ports), len(kept)), dtype=np.int64)
    expansion6[kept, np.arange(len(kept))] = 6
    for k, diamond in zip(dropped, diamonds):
        relation = np.array(diamond['relation_coefficients'])
        assert relation[k] == 6 and all(relation[j] == 0 for j in dropped if j != k)
        expansion6[k] = -relation[kept]
    assert np.array_equal(expansion6@output[kept], 6*output)
    memory_input = (centered@expansion6)/(6*clock*sizes[labels, None])
    memory_output = output[kept]/np.sqrt(sizes)

    step = lambda v: s.counts@v/clock
    projected_step = basis.T@step(basis)
    coupling = step(basis)-basis@projected_step
    assert np.allclose(coupling, memory_input@memory_output, atol=1e-13)
    hidden_step = lambda v: step(v)-basis@(basis.T@step(v))
    covariance, memories = [np.eye(len(patterns))], []
    current = basis.copy()
    hidden = memory_input.copy()
    for k in range(9):
        current = step(current)
        covariance.append(basis.T@current)
        reduced_memory = memory_input.T@hidden
        memories.append(memory_output.T@reduced_memory@memory_output)
        hidden = hidden_step(hidden)
    assert np.allclose(memories[0], coupling.T@coupling, atol=1e-13)
    errors = []
    for t in range(1, 9):
        predicted = projected_step@covariance[t]
        for j in range(t):
            predicted += memories[j]@covariance[t-1-j]
        errors.append(float(np.max(np.abs(predicted-covariance[t+1]))))
    assert max(errors) < 1e-11

    prepared = np.zeros(n)
    prepared[s.seeds[1]], prepared[s.seeds[0]] = 1, -1
    assert np.all(basis.T@prepared == 0)
    hidden_preparation = prepared.copy()
    visible = [np.zeros(len(patterns))]
    preparation_errors = []
    for t in range(9):
        forcing = memory_output.T@(memory_input.T@hidden_preparation)
        predicted = projected_step@visible[-1]+forcing
        for j in range(t):
            predicted += memories[j]@visible[t-1-j]
        visible.append(predicted)
        prepared = step(prepared)
        preparation_errors.append(float(np.max(np.abs(predicted-basis.T@prepared))))
        hidden_preparation = hidden_step(hidden_preparation)
    assert max(preparation_errors) < 1e-11
    mean_difference_two = visible[2]@(np.array(patterns)*np.sqrt(sizes)[:, None])
    assert np.allclose(mean_difference_two, np.array([-16, 8, -8, 16, -8, 8])/clock**2, atol=1e-13)

    times = (1, 2, 4, 8, 16, 32, 64)
    sampled = {0: np.eye(len(patterns))}
    current = basis.copy()
    requested = set(times) | {2*t for t in times}
    for t in range(1, 2*max(times)+1):
        current = step(current)
        if t in requested:
            sampled[t] = basis.T@current
    gains = []
    for t in times:
        gram = sampled[2*t]-sampled[t]@sampled[t]
        gain = float(np.sqrt(max(0, np.linalg.eigvalsh((gram+gram.T)/2)[-1])))
        if t % 2 == 0:
            assert gain <= .5+1e-10
        gains.append({'attempts': t, 'gauge_to_complete_charge_L2_gain': gain})

    words = WordObserver(s.g, s.f.q)
    unframed = [words.canonical(tuple(s.f.oracle.transport(raw, path) for path in s.paths)) for raw in s.raw]
    multiplicities = sorted(set(Counter(unframed).values()))
    assert multiplicities == [3]
    return {'framed_states': n, 'complete_charge_patterns': len(patterns),
            'charge_fiber_sizes_and_counts': [[int(size), int(np.sum(sizes == size))] for size in sorted(set(sizes))],
            'exact_one_step_hidden_coupling_rank': rank,
            'activity_factorization': {'activity_ports': len(ports), 'centered_activity_rank': activity_rank,
                'charge_output_rank': output_rank, 'disjoint_reaction_diamonds': diamonds,
                'exact_reduced_memory_dimension': len(kept), 'eliminated_ports': [ports[k] for k in dropped]},
            'individual_charge_means_annihilate_one_step_hidden_coupling': True,
            'unframed_states': len(set(unframed)), 'framed_lifts_per_unframed_state': 3,
            'memory_recurrence_maximum_error': max(errors),
            'prepared_alignment_recurrence_maximum_error': max(preparation_errors),
            'prepared_alignment_two_attempt_mean_difference': mean_difference_two.tolist(),
            'first_memory_trace': float(np.trace(memories[0])),
            'memory_operator_norms': [float(np.linalg.norm(k, 2)) for k in memories],
            'gauge_to_charge_transfer': gains,
            'scope': 'Exact projection definition and integer coupling rank; floating-point memory and gain evaluation. No autonomous charge closure or physical instability inferred.'}


if __name__ == '__main__':
    print(json.dumps(calculation(), indent=2))
