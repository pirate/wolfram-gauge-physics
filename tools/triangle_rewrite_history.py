#!/usr/bin/env python3
"""Actual closed rewrite walks versus auxiliary spin holonomy.

This analyzes both labeled changing events and their raw-state adjacency
skeleton. Removing idle attempts here is a geometric diagnostic, not a
change to the simulation clock or a declaration of physical equivalence.
"""
import json
import itertools
from collections import Counter

from triangle_local_control_group import calculation, compose
from triangle_projective_lift import projective_calculation


def reduce_walk(walk):
    stack = []
    for edge in walk:
        if stack and stack[-1] == -edge:
            stack.pop()
        else:
            stack.append(edge)
    return stack


def history_calculation():
    local, spin = calculation(), projective_calculation()
    tables, channels = local['framed_tables'], local['channels']
    count = len(local['framed_states'])
    inverses = [channels.index((p, 3-r if r in (1, 2) else r)) for p, r in channels]
    assert all(compose(tables[k], tables[inverses[k]]) == tuple(range(count)) for k in range(len(channels)))
    simple_edges = sorted({tuple(sorted((i, t[i]))) for t in tables for i in range(count) if t[i] != i})
    simple_ids = {e: j+1 for j, e in enumerate(simple_edges)}
    event_edges = sorted({min((i, k), (t[i], inverses[k]))
                          for k, t in enumerate(tables) for i in range(count) if t[i] != i})
    event_ids = {e: j+1 for j, e in enumerate(event_edges)}

    # A spanning tree gives an exact free-group word for a based walk in
    # either graph. No spanning-tree choice changes whether it is trivial.
    def chord_ids(labeled):
        edges = event_edges if labeled else simple_edges
        parent = list(range(count))
        def root(i):
            while parent[i] != i:
                i = parent[i]
            return i
        chords = set()
        for j, (u, v) in enumerate(edges, 1):
            if labeled:
                v = tables[v][u]
            a, b = root(u), root(v)
            if a == b:
                chords.add(j)
            else:
                parent[a] = b
        assert len({root(i) for i in range(count)}) == 1
        return chords
    chords = {'event': chord_ids(True), 'state': chord_ids(False)}

    start = spin['fixed_classical_state_id']
    def walk(word):
        state = start
        result = {'states': [state], 'event': [], 'state': [], 'idle_attempts': 0}
        for k in word:
            target = tables[k][state]
            if target == state:
                result['idle_attempts'] += 1
            else:
                arc, reverse = (state, k), (target, inverses[k])
                result['event'].append((1 if arc < reverse else -1)*event_ids[min(arc, reverse)])
                result['state'].append((1 if state < target else -1)*simple_ids[tuple(sorted((state, target)))])
            state = target
            result['states'].append(state)
        assert state == start
        for kind in ('event', 'state'):
            result[kind+'_reduced'] = reduce_walk(result[kind])
            result[kind+'_free_word'] = reduce_walk([e for e in result[kind] if abs(e) in chords[kind]])
            result[kind+'_net_edge_counts'] = {str(e): n for e, n in sorted(Counter(
                {abs(v): result[kind].count(abs(v))-result[kind].count(-abs(v)) for v in result[kind]}).items()) if n}
        return result

    a, b = spin['commuting_classical_words']
    reverse_word = lambda w: [inverses[k] for k in reversed(w)]
    commutator = a+b+reverse_word(a)+reverse_word(b)
    walks = {'A': walk(a), 'B': walk(b), 'commutator': walk(commutator)}
    assert all(not walks['commutator'][kind+'_net_edge_counts'] for kind in ('event', 'state'))
    reps, tau = local['representative_state_ids'], local['frame_partner']
    def odd_trace(table):
        return sum(1 if table[r] == r else -1 if table[r] == tau[r] else 0 for r in reps)
    identity = tuple(range(count))
    involutions = [k for k, t in enumerate(tables) if compose(t, t) == identity]
    idle_spin_witnesses = []
    for k, ell in itertools.combinations(involutions, 2):
        product = compose(tables[k], tables[ell])
        fixed = [i for i in range(count) if tables[k][i] == tables[ell][i] == i]
        if product != compose(tables[ell], tables[k]) or not fixed:
            continue
        # For commuting orthogonal involutions, this integer is the dimension
        # of the intersection of their negative eigenspaces. Both negative
        # dimensions here are even, so odd intersection gives spin commutator -1.
        numerator = count//2-odd_trace(tables[k])-odd_trace(tables[ell])+odd_trace(product)
        assert numerator % 4 == 0
        intersection = numerator//4
        if intersection % 2:
            idle_spin_witnesses.append({'channel_indices': [k, ell],
                'channels': [channels[k], channels[ell]], 'common_fixed_states': fixed,
                'odd_traces': [odd_trace(tables[k]), odd_trace(tables[ell]), odd_trace(product)],
                'negative_plane_intersection_dimension': intersection})
    idle = next(w for w in idle_spin_witnesses if start in w['common_fixed_states'])
    k, ell = idle['channel_indices']
    idle['based_word'] = local['framed_states'][start]
    idle['AB_state_history'] = [start, tables[k][start], tables[ell][tables[k][start]]]
    idle['BA_state_history'] = [start, tables[ell][start], tables[k][tables[ell][start]]]
    assert idle['AB_state_history'] == idle['BA_state_history'] == [start]*3
    # Each global involution fixes some states and has square -1 in its spin
    # lift. A connection on the actual state graph instead cancels retracing.
    witness_k = next(k for k, record in enumerate(spin['primitive_odd_actions'])
                     if record['spin_lift_square_sign_if_involution'] == -1)
    witness_x = next(i for i in range(count) if tables[witness_k][i] != i)
    witness_y = tables[witness_k][witness_x]
    assert tables[witness_k][witness_y] == witness_x
    return {'vertices': count,
            'changing_event_edges': len(event_edges), 'event_graph_free_rank': len(chords['event']),
            'raw_state_adjacency_edges': len(simple_edges), 'state_graph_free_rank': len(chords['state']),
            'edge_id_conventions': {'event': event_edges, 'state': simple_edges},
            'primitive_channel_indices': {'A': a, 'B': b, 'commutator': commutator},
            'walks': walks,
            'commuting_primitive_pairs_with_anticommuting_spin_lifts': idle_spin_witnesses,
            'two_idle_attempts_spin_witness': idle,
            'retracing_witness': {'channel_index': witness_k, 'channel': channels[witness_k],
                                  'state_path': [witness_x, witness_y, witness_x],
                                  'global_spin_lift_square': -1,
                                  'inverse_compatible_connection_holonomy': 1},
            'scalar_inverse_compatible_commutator_holonomy': 1,
            'auxiliary_spin_commutator': spin['projective_commutator_sign'],
            'scope': 'Exact walk graphs without attached two-cells. Not a spacetime topology, a physical history quotient, or derived amplitudes.'}


if __name__ == '__main__':
    print(json.dumps(history_calculation()), flush=True)
