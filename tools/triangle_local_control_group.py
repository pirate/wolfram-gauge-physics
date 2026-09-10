#!/usr/bin/env python3
"""Exact ordered-move algebra on a boundary-framed three-face disk.

All transitions are read from actual shared links. Group generation concerns
controllability by chosen words, not an autonomously selected physical law.
"""
import itertools
import json
import math
from collections import Counter
from pathlib import Path

from triangle_elastic_scattering import ElasticExperiment, inverse
from triangle_reference import subgroup


def compose(a, b):
    return tuple(a[b[i]] for i in range(len(b)))


def parity(p):
    return sum(p[i] > p[j] for i in range(len(p)) for j in range(i+1, len(p))) % 2


def binary_rank(words):
    basis = {}
    for value in words:
        while value:
            pivot = value.bit_length()-1
            if pivot not in basis:
                basis[pivot] = value
                break
            value ^= basis[pivot]
    return len(basis)


def calculation():
    bank = json.loads(Path('data/triangle-feedback.json').read_text())['bank']
    elastic = json.loads(Path('data/triangle-elastic-scattering.json').read_text())['pair_census']['elastic_tables']
    x = ElasticExperiment(6, bank, elastic)
    oracle, g, q = x.e.factor.oracle, x.e.geometry.group, x.e.charges
    fan = 0
    edges, internal = set(oracle.fan.reads[fan]), set(oracle.fan.writes[fan])
    boundary = edges-internal
    p = next(a for a in range(6) if q[a] == 1)
    word = lambda raw: tuple(oracle.fan.transport(raw, path) for path in oracle.fan.patches[fan][::-1])
    initial = [g.identity]*len(x.e.geometry.edges)
    initial[min(boundary)] = p
    # Keep an actual noncommuting exterior loop at the same root. Otherwise
    # the residual boundary conjugation could still be a global gauge change.
    root = oracle.fan.specs[fan][0][0]
    spectator = next(i for i, path in enumerate(oracle.fan.face_paths)
                     if x.e.geometry.faces[i][0] == root
                     and not internal.intersection(e for e, _ in path)
                     and any(e not in edges for e, _ in path))
    spectator_path = oracle.fan.face_paths[spectator]
    k = next(i for i, (e, _) in enumerate(spectator_path) if e not in edges)
    z = next(a for a in range(6) if q[a] == 2)
    a = oracle.fan.transport(initial, spectator_path[k+1:])
    b = oracle.fan.transport(initial, spectator_path[:k])
    value = g.mul[g.inv[a]][g.mul[z][g.inv[b]]]
    edge, reverse = spectator_path[k]
    initial[edge] = g.inv[value] if reverse else value
    assert oracle.fan.transport(initial, spectator_path) == z
    raw_states, words = [], []
    for values in itertools.product(range(6), repeat=2):
        raw = initial[:]
        for edge, value in zip(sorted(internal), values):
            raw[edge] = value
        w = word(raw)
        assert g.mul[g.mul[w[0]][w[1]]][w[2]] == p
        if sum(q[a] for a in w) == 3 and len(subgroup(g, w)) == 6:
            raw_states.append(raw)
            words.append(w)
    ids = {w: i for i, w in enumerate(words)}
    assert len(ids) == len(words) == 20
    tau = tuple(ids[tuple(g.conj[p][a] for a in w)] for w in words)
    assert all(tau[i] != i and tau[tau[i]] == i for i in range(20))
    odd = lambda w: tuple((q[g.mul[z][g.mul[p][a]]]-q[g.mul[g.inv[z]][g.mul[p][a]]])//2 for a in w)
    assert all(any(odd(w)) and odd(words[tau[i]]) == tuple(-v for v in odd(w)) for i, w in enumerate(words))
    pairs = [i for i, paths in enumerate(oracle.pairs) if {e for path in paths for e, _ in path} <= edges]
    fans = [i for i, reads in enumerate(oracle.fan.reads) if reads <= edges]
    channels = [(i, r) for i in pairs for r in range(3)]+[(i, r) for i in fans for r in range(3, 15)]
    tables = []
    for patch, rule in channels:
        table = []
        for raw in raw_states:
            moved = raw[:]
            x.apply(moved, rule*x.supports+patch)
            assert all(moved[e] == raw[e] for e in range(len(raw)) if e not in internal)
            assert oracle.fan.transport(moved, spectator_path) == z
            table.append(ids[word(moved)])
        inverse(table)
        tables.append(tuple(table))
    assert all(compose(t, tau) == compose(tau, t) for t in tables)
    reps = [i for i in range(20) if i < tau[i]]
    base = {i: k for k, rep in enumerate(reps) for i in (rep, tau[rep])}
    quotient = [tuple(base[t[i]] for i in reps) for t in tables]
    assert not any(parity(t) for t in tables+quotient)
    selected = []
    for position in (0, 1):
        target = []
        for w in words:
            moved = list(w)
            a, b = moved[position:position+2]
            if q[a] == q[b] == 1 and a != b:
                moved[position:position+2] = [g.conj[a][b], a]
            target.append(ids[tuple(moved)])
        selected.append(tables.index(tuple(target)))
    cycle = quotient[selected[0]]
    assert sum(i != j for i, j in enumerate(cycle)) == 3
    conjugates = {cycle: []}
    queue = [cycle]
    for old in queue:
        for k, t in enumerate(quotient):
            new = compose(t, compose(old, inverse(t)))
            if new not in conjugates:
                conjugates[new] = conjugates[old]+[k]
                queue.append(new)
    all_three_cycles = set()
    for a, b, c in itertools.combinations(range(10), 3):
        t = list(range(10)); t[a], t[b], t[c] = b, c, a
        all_three_cycles.update((tuple(t), tuple(inverse(t))))
    quotient_is_alternating = set(conjugates) == all_three_cycles
    twist = tuple(range(20))
    circuit = selected*3
    for k in circuit:
        twist = compose(tables[k], twist)
    reflections = [i for i, w in enumerate(words) if all(q[a] == 1 for a in w)]
    assert all(twist[i] == (tau[i] if i in reflections else i) for i in range(20))
    assert all(base[twist[i]] == base[i] for i in range(20))
    flips = sum(1 << k for k, i in enumerate(reps) if twist[i] == tau[i])
    sign_orbit, sign_queue = {flips}, [flips]
    for old in sign_queue:
        for t in quotient:
            new = sum(1 << t[i] for i in range(10) if old >> i & 1)
            if new not in sign_orbit:
                sign_orbit.add(new); sign_queue.append(new)
    rank = binary_rank(sign_orbit)
    star_words = []
    if quotient_is_alternating:
        for i in range(2, 10):
            target = list(range(10)); target[0], target[1], target[i] = 1, i, 0
            star_words.append({'cycle': [0, 1, i], 'conjugator_channel_indices': conjugates[tuple(target)]})
    # A short controlled conversion cycle, replayed on actual shared links.
    controlled_cycle = None
    if star_words:
        w = star_words[0]['conjugator_channel_indices']
        reverse = [tables.index(tuple(inverse(tables[k]))) for k in reversed(w)]
        local_circuit = reverse+[selected[0]]+w
        permutation = tuple(range(20))
        for k in local_circuit:
            permutation = compose(tables[k], permutation)
        for i, raw in enumerate(raw_states):
            moved = raw[:]
            for k in local_circuit:
                patch, rule = channels[k]
                x.apply(moved, rule*x.supports+patch)
            assert word(moved) == words[permutation[i]]
            assert all(moved[e] == raw[e] for e in range(len(raw)) if e not in internal)
        orbit = [reps[0], permutation[reps[0]], permutation[permutation[reps[0]]]]
        assert permutation[orbit[-1]] == orbit[0]
        controlled_cycle = {'quotient_cycle': [0, 1, 2], 'channel_indices': local_circuit,
                            'channels': [channels[k] for k in local_circuit],
                            'based_word_cycle': [words[i] for i in orbit],
                            'charge_cycle': [[q[a] for a in words[i]] for i in orbit],
                            'framed_states_moved': sum(i != j for i, j in enumerate(permutation)),
                            'all_20_inputs_replayed_on_actual_links': True}
    return {'fan': fan, 'boundary_holonomy': p, 'boundary_link_values': [[e, initial[e]] for e in sorted(boundary)],
            'exterior_reference_face': spectator, 'exterior_reference_holonomy': z,
            'internal_edges': sorted(internal), 'framed_states': words, 'frame_partner': tau,
            'representative_state_ids': reps, 'channels': channels, 'framed_tables': tables,
            'quotient_tables': quotient, 'all_raw_and_quotient_generators_even': True,
            'hurwitz_channel_indices': selected, 'three_cycle_conjugacy_orbit_size': len(conjugates),
            'quotient_equals_A10': quotient_is_alternating, 'star_three_cycle_conjugators': star_words,
            'central_twist_channel_indices': circuit, 'central_twist_flipped_pairs': flips.bit_count(),
            'sign_orbit_size': len(sign_orbit), 'sign_kernel_binary_rank': rank,
            'controlled_cycle': controlled_cycle,
            'framed_group_order': 2**rank*math.factorial(10)//2 if quotient_is_alternating else None,
            'classification': 'Even-weight binary flips semidirect A10' if quotient_is_alternating and rank == 9 else 'Not classified',
            'scope': 'Exact boundary-held ordered control group, not a spontaneous schedule, a Hamiltonian, a quantum amplitude space, or molecular motion.'}


if __name__ == '__main__':
    print(json.dumps(calculation()), flush=True)
