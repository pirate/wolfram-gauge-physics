#!/usr/bin/env python3
"""Construct net angular activation from the old rules and four constraints."""
import json
from collections import Counter
from pathlib import Path

import numpy as np

from fiber_constraint_dynamics import loop_row
from fiber_constraint_intersection import angular_rate_drift, clone
from fiber_reaction_bursts import activity, write_fan
from fiber_refinement_limit import primitive, read
from fiber_three_constraints import diagnostics
from fiber_two_constraint_witness import rates
from fiber_two_constraints import AffineProcess, parameterize


def construct(side=4, seed=914407001, birth_distance=1):
    state = AffineProcess('contact', side, seed)
    incident = [[] for _ in state.oracle.edges]
    for f, edges in enumerate(state.fe):
        for e in edges:
            incident[e].append(f)
    neighbors = [[] for _ in range(state.F)]
    for a, b in incident:
        neighbors[a].append(b); neighbors[b].append(a)
    Eface, a, b = map(int, state.supports[0])
    Zface, leaf = (a, b) if birth_distance == 1 else (b, a)
    assert leaf in neighbors[Zface]
    hub = min(f for f in neighbors[leaf] if f not in (Eface, Zface))
    assert all(f not in (Eface, Zface) for f in neighbors[hub])
    signs = state.initial_signs.copy()
    dual_path = [Eface, Zface] if birth_distance == 1 else [Eface, leaf, Zface]
    for u, v in zip(dual_path, dual_path[1:]):
        shared = next(e for e, fs in enumerate(incident) if set(fs) == {u, v})
        signs[shared] *= -1
    rows = [loop_row(signs, state.fe[Eface], state.fr[Eface])[1]]
    for other in neighbors[hub]:
        p = next(p for p, fs in enumerate(state.supports) if set(map(int, fs[:2])) == {hub, other})
        a = loop_row(signs, state.pe[p, 2], state.pr[p, 2])[1]
        b = loop_row(signs, state.pe[p, 1], state.pr[p, 1])[1]
        rows.append(tuple(x-y for x, y in zip(a, b)))
    forms = parameterize(rows)
    state.raw = [(s, form) for s, form in zip(signs, forms)]
    state.rank = 4
    state.initial_rows = rows
    state.initial_signs = signs
    q = state.charges()
    assert q[Eface] == 0 and q[Zface] == 2 and q.count(0) == q.count(2) == 1
    drift, changes = angular_rate_drift(state)
    positive = next(change for change in changes if change['delta_reaction_rate'] > 0)
    after = clone(state)
    after.step(positive['fan'], positive['channel'], 0., None, (), 0)
    # Reverse an available original reaction, then apply that same
    # involution to obtain the constructed state from an all-R state.
    for channel in range(3, 15):
        initial = clone(state)
        initial.step(0, channel, 0., None, (), 0)
        if initial.charges() == [1]*state.F:
            break
    else:
        raise AssertionError('No original inverse reaction found')
    restored = clone(initial)
    restored.step(0, channel, 0., None, (), 0)
    assert restored.raw == state.raw
    initial_response = None
    if side == 4:
        response_counts = Counter()
        for p, rate in enumerate(rates(initial)):
            if not rate:
                continue
            for c in range(3, 15):
                middle = clone(initial)
                middle.step(p, c, 0., None, (), 0)
                response_counts[angular_rate_drift(middle)[0]] += 1
        initial_response = {'G0_A_lambda': sum(a*b for a, b in response_counts.items()),
                            'post_reaction_A_lambda_counts': sorted(response_counts.items())}
    # Actual finite raw links follow the same two primitive updates.
    n = 1_000_000_007
    parameters = np.random.default_rng(seed).integers(n, size=state.E-state.rank)
    def evaluate(s):
        return np.array([(sum(c*int(z) for c, z in zip(h[1][:-1], parameters))+h[1][-1]) % n
                         +(n if h[0] == -1 else 0) for h in s.raw], dtype=np.int64)
    raw = evaluate(initial)
    saved = [raw.tolist()]
    for p, c, target in ((0, channel, state), (positive['fan'], positive['channel'], after)):
        word = tuple(read(raw, state.pe[p, 2-i], state.pr[p, 2-i], n) for i in range(3))
        write_fan(raw, p, *primitive(*word, c, n), n, state.pe, state.pr)
        assert np.array_equal(raw, evaluate(target))
        assert rates(target) == [activity(raw, k, n, state.pe, state.pr) for k in range(len(state.pe))]
        saved.append(raw.tolist())
    def flat_graph(s):
        q = s.charges()
        flat = set()
        for p, faces in enumerate(s.supports):
            a, b = map(int, faces[:2])
            if q[a] == q[b] == 1:
                word = s.word(p)
                if word[0][1] == word[1][1]:
                    flat.add(tuple(sorted((a, b))))
        k = Counter(f for edge in flat for f in edge)
        d = {f: sum(q[g] == 1 for g in neighbors[f]) for f in range(s.F) if q[f] == 1}
        forward_rate = 12*sum(k[f]*(d[f]-k[f]) for f in d)
        backward_rate = 4*sum(sorted(q[int(f)] for f in faces) == [0, 1, 2] for faces in s.supports)
        assert forward_rate+backward_rate == sum(rates(s))
        return {'flat_edges': sorted(flat), 'forward_rate': forward_rate, 'backward_rate': backward_rate,
                'changed_face_R_degree': d[leaf], 'changed_face_flat_degree': k[leaf],
                'hub_R_degree': d[hub], 'hub_flat_degree': k[hub]}
    return {'side': side, 'faces': state.F, 'rank': state.rank, 'seed': seed, 'birth_distance': birth_distance,
            'identity_face': Eface, 'rotation_face': Zface, 'changed_reflection_face': leaf,
            'flat_star_hub': hub, 'flat_star_leaves': neighbors[hub],
            'before_angular_signs': signs, 'before_angular_constraints': rows,
            'all_reflection_initial': diagnostics(initial),
            'before_angular': diagnostics(state), 'after_angular': diagnostics(after),
            'A_lambda': drift, 'angular_changes': changes,
            'flat_graph_before': flat_graph(state), 'flat_graph_after': flat_graph(after),
            'all_first_reactions_response': initial_response,
            'reaction_then_angular_sequence': [[0, channel], [positive['fan'], positive['channel']]],
            'unchanged_charge_field': q, 'finite_modulus': n,
            'finite_raw_links_initial_before_after': saved}


if __name__ == '__main__':
    result = {'scope': 'Prepared four-constraint witness; no new transition, energy, or attraction', 'runs': []}
    for side in (4, 8, 16):
      for birth_distance in (1, 2):
        row = construct(side, birth_distance=birth_distance)
        result['runs'].append(row)
        Path('data/fiber-angular-activation.json').write_text(json.dumps(result, indent=2)+'\n')
        print(json.dumps({k: v for k, v in row.items() if k in ('side', 'faces', 'birth_distance', 'all_first_reactions_response', 'all_reflection_initial', 'before_angular',
                                                              'after_angular', 'A_lambda', 'angular_changes',
                                                              'reaction_then_angular_sequence')}), flush=True)
