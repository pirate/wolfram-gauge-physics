#!/usr/bin/env python3
"""Exact angular-response census after adding a third local relation."""
import json
from collections import Counter
from pathlib import Path

from flint import fmpz_mat

from fiber_constraint_dynamics import canonical, loop_row
from fiber_constraint_intersection import angular_rate_drift, clone, intersection_matrix
from fiber_two_constraint_witness import rates
from fiber_two_constraints import AffineProcess, parameterize


def prepare(base, p, i):
    state = clone(base)
    a = loop_row(base.initial_signs, base.pe[p, 2-i], base.pr[p, 2-i])[1]
    b = loop_row(base.initial_signs, base.pe[p, 1-i], base.pr[p, 1-i])[1]
    rows = base.initial_rows+[tuple(x-y for x, y in zip(a, b))]
    if fmpz_mat(rows).rank() != 3:
        return None
    forms = parameterize(rows)
    state.raw = [(s, form) for s, form in zip(base.initial_signs, forms)]
    state.rank = 3
    state.initial_rows = rows
    state.initial_pair_faces = base.initial_pair_faces+[base.supports[p][i:i+2].tolist()]
    return state


def diagnostics(state):
    B = fmpz_mat([list(h[1][:-1]) for h in state.raw])
    kernel, rank = B.transpose().nullspace()
    assert rank == state.rank
    rows = [list(canonical([int(kernel[e, i]) for e in range(state.E)])[0]) for i in range(rank)]
    signs = [h[0] for h in state.raw]
    N = fmpz_mat(rows)
    offsets = [sum(x*h[1][-1] for x, h in zip(row, state.raw)) for row in rows]
    q = state.charges()
    assert q.count(0) == q.count(2) and sum(q) == state.F
    return {'identity_count': q.count(0), 'reaction_rate': sum(rates(state)),
            'zero_offset_rank': rank-int(any(offsets)),
            'intersection_rank': (N*intersection_matrix(state.oracle, signs)*N.transpose()).rank(),
            'all_equal_RRR_fans': sum(all(h[0] == -1 for h in state.word(p))
                                     and len({h[1] for h in state.word(p)}) == 1
                                     for p in range(len(state.pe))),
            'joint_support': sum(any(row[e] for row in rows) for e in range(state.E))}


def census(side=4, seed=914407001):
    base = AffineProcess('contact', side, seed)
    seen, results = set(), []
    for p in range(len(base.pe)):
        for i in (0, 1):
            pair = tuple(sorted(map(int, base.supports[p][i:i+2])))
            if pair in seen:
                continue
            seen.add(pair)
            state = prepare(base, p, i)
            if state is None:
                continue
            initial = rates(state)
            angular_deltas, drifts = Counter(), Counter()
            coefficient, reaction_count, example = 0, 0, None
            for fan, intensity in enumerate(initial):
                if not intensity:
                    continue
                for channel in range(3, 15):
                    middle = clone(state)
                    middle.step(fan, channel, 0., None, (), 0)
                    reaction_count += 1
                    drift, changes = angular_rate_drift(middle)
                    coefficient += drift
                    drifts[drift] += 1
                    for change in changes:
                        delta = change['delta_reaction_rate']
                        angular_deltas[delta] += 1
                        if example is None and delta > 0:
                            after = clone(middle)
                            after.step(change['fan'], change['channel'], 0., None, (), 0)
                            example = {'first_reaction': [fan, channel], 'angular_step': change,
                                       'before': diagnostics(middle), 'after': diagnostics(after),
                                       'unchanged_charge_field': middle.charges(),
                                       'initial_signs': state.initial_signs, 'initial_rows': state.initial_rows}
                            assert middle.charges() == after.charges()
            assert reaction_count == sum(initial)
            row = {'side': side, 'seed': seed, 'third_pair_selection': [p, i],
                   'initial_pairs': state.initial_pair_faces, 'initial': diagnostics(state),
                   'G0_A_lambda': coefficient, 'angular_delta_counts': sorted(angular_deltas.items()),
                   'post_reaction_A_lambda_counts': sorted(drifts.items()), 'positive_step_example': example}
            results.append(row)
            Path('data/fiber-three-constraints.json').write_text(json.dumps({'scope': 'All distinct third adjacent-pair constraints added to a fixed overlapping pair',
                                                                           'runs': results}, indent=2)+'\n')
            print(json.dumps({k: v for k, v in row.items() if k not in ('positive_step_example', 'post_reaction_A_lambda_counts')}), flush=True)
    return results


if __name__ == '__main__':
    census()
