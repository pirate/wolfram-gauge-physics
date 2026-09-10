#!/usr/bin/env python3
"""Find and retain a legal released rank-two angular-feedback witness."""
import copy
import json
from pathlib import Path

import numpy as np

from fiber_two_constraints import AffineProcess
from fiber_reaction_bursts import activity


def clone(state):
    out = copy.copy(state)
    out.raw = state.raw.copy()
    out.events = state.events.copy()
    return out


def rates(state):
    out = []
    for p in range(len(state.pe)):
        w = state.word(p)
        q = [1 if h[0] == -1 else 0 if not any(h[1]) else 2 for h in w]
        if q == [1, 1, 1]:
            out.append(12 if (w[0][1] == w[1][1]) != (w[1][1] == w[2][1]) else 0)
        else:
            out.append(4 if sorted(q) == [0, 1, 2] else 0)
    return out


def search():
    seed = 914407001
    for condition in ('contact', 'separated'):
        start = AffineProcess(condition, 4, seed)
        original_rates = rates(start)
        for p in range(len(start.pe)):
            if original_rates[p] == 0:
                continue
            for channel in range(3, 15):
                middle = clone(start)
                middle.step(p, channel, 0., None, (), 0)
                before_rates = rates(middle)
                if not middle.events:
                    continue
                for q in range(len(start.pe)):
                    for angular in (15, 16):
                        after = clone(middle)
                        after.step(q, angular, 1., None, (), 0)
                        after_rates = rates(after)
                        if before_rates == after_rates:
                            continue
                        assert middle.charges() == after.charges()
                        modulus = 1_000_000_007
                        parameters = np.random.default_rng(seed).integers(modulus, size=start.E-start.rank)
                        def evaluate(state):
                            return np.array([(sum(c*int(x) for c, x in zip(h[1][:-1], parameters))+h[1][-1]) % modulus
                                             +(modulus if h[0] == -1 else 0) for h in state.raw], dtype=np.int64)
                        finite_before, finite_after = evaluate(middle), evaluate(after)
                        assert before_rates == [activity(finite_before, k, modulus, middle.pe, middle.pr) for k in range(len(middle.pe))]
                        assert after_rates == [activity(finite_after, k, modulus, after.pe, after.pr) for k in range(len(after.pe))]
                        changed = [i for i, (a, b) in enumerate(zip(before_rates, after_rates)) if a != b]
                        return {'condition': condition, 'seed': seed, 'side': 4, 'rank': 2,
                                'initial_signs': start.initial_signs, 'initial_constraints': start.initial_rows,
                                'initial_pair_faces': start.initial_pair_faces, 'conserved_charge': start.Q,
                                'sequence': [{'fan': p, 'channel': channel, 'kind': 'original reaction'},
                                             {'fan': q, 'channel': angular, 'kind': 'original angular generator direction'}],
                                'initial_rates': original_rates, 'before_angular_rates': before_rates,
                                'after_angular_rates': after_rates, 'changed_fans': changed,
                                'changed_fan_faces': [start.supports[i].tolist() for i in changed],
                                'unchanged_charge_field': middle.charges(), 'finite_modulus': modulus,
                                'finite_before_angular_links': finite_before.tolist(), 'finite_after_angular_links': finite_after.tolist(),
                                'before_features': middle.features(), 'after_features': after.features()}
    raise RuntimeError('No two-step witness found; do not infer an effect without evidence')


if __name__ == '__main__':
    result = search()
    Path('data/fiber-two-constraint-witness.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k: v for k, v in result.items() if k in ('condition', 'sequence', 'changed_fans', 'changed_fan_faces',
                                                             'before_features', 'after_features')}), flush=True)
