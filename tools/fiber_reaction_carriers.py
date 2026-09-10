#!/usr/bin/env python3
"""Event-conditioned reaction intensity and local flatness carriers."""
import json
from pathlib import Path

import numpy as np
from numba import njit

from fiber_reaction_bursts import activity, local_geometry, palm_raw, summarize, write_fan
from fiber_refinement_limit import kind, read, primitive, mul


@njit
def scan(raw, n, target, neighbors, pe, pr, fe, fr, supports):
    total = 0
    for p in range(len(pe)):
        total += activity(raw, p, n, pe, pr)
    q = np.array([kind(read(raw, fe[f], fr[f], n), n) for f in range(len(fe))])
    forward = np.any(q[supports[0]] == 0)
    created = 0
    if forward:
        for p in neighbors:
            shared = -1
            for f in supports[p]:
                for g in supports[0]:
                    if f == g:
                        shared = f
            if q[shared] == 0 and activity(raw, p, n, pe, pr) == 4:
                # Parent input was RRR; the exterior face classes do not
                # change in that parent rewrite. This ERZ gate was inactive.
                created += 1
    return np.array((activity(raw, 0, n, pe, pr), activity(raw, target, n, pe, pr),
                     total, int(forward), total*int(forward), total*int(not forward), created))


def experiment(side, trials, n=1_000_000_007, seed=912305001):
    oracle, args, target, geo = local_geometry(side)
    neighbors = np.array(geo['single_face_private_edge_targets'])
    root = args[-1][0]
    face_containers = [[int(p) for p, fs in enumerate(args[-1]) if f in fs] for f in root]
    pair_containers = [[int(p) for p, fs in enumerate(args[-1]) if a in fs and b in fs]
                       for a, b in zip(root, root[1:])]
    assert [len(x) for x in face_containers] == [9, 9, 9]
    assert [len(x) for x in pair_containers] == [4, 4]
    rng = np.random.default_rng(seed)
    x = np.array([scan(palm_raw(rng, n, oracle, args), n, target, neighbors, *args)
                  for _ in range(trials)])
    f = x[:, 3].mean()
    rates = np.array((x[:, 4].mean()/f, x[:, 5].mean()/(1-f)))
    influences = np.column_stack(((x[:, 4]-rates[0]*x[:, 3])/f,
                                  (x[:, 5]-rates[1]*(1-x[:, 3]))/(1-f)))
    return {'side': side, 'faces': len(oracle.faces), 'n': n, 'trials': trials, 'seed': seed,
            'geometry': geo, 'fans_per_root_face': face_containers, 'fans_per_root_adjacent_pair': pair_containers,
            'columns': ['marked_rate', 'neighbor_rate', 'total_rate', 'parent_forward',
                        'total_rate_times_forward', 'total_rate_times_backward', 'new_single_overlap_gates'],
            **summarize(x), 'conditional_total_rates_forward_backward': rates.tolist(),
            'conditional_rate_SE': (influences.std(axis=0, ddof=1)/np.sqrt(trials)).tolist(),
            'predicted_large_n_total_rate': 25., 'predicted_large_n_conditional_total_rates': [20., 30.],
            'exact_marked_rate': 8., 'exact_neighbor_rate': 1/3+16*(n-1)/(3*n**2),
            'exact_new_single_overlap_gates_mean': (n-1)/n}


def transport_witness(n=5, seed=912407001):
    oracle, args, target, geo = local_geometry(4)
    pe, pr, fe, fr, supports = args
    shared = next(iter(set(supports[0]) & set(supports[target])))
    rng = np.random.default_rng(seed)
    charges = lambda raw: np.array([kind(read(raw, fe[f], fr[f], n), n) for f in range(len(fe))])
    word = lambda raw, p: tuple(read(raw, pe[p, 2-i], pr[p, 2-i], n) for i in range(3))
    while True:
        middle = palm_raw(rng, n, oracle, args)
        qm = charges(middle)
        if qm[shared] == 0 and activity(middle, target, n, pe, pr) == 4:
            break
    wp = word(middle, 0)
    parent_channel = next(ch for ch in range(3, 15) if primitive(*wp, ch, n) != wp)
    before = middle.copy()
    write_fan(before, 0, *primitive(*wp, parent_channel, n), n, pe, pr)
    before_word = word(before, 0)
    assert all(h >= n for h in before_word)
    assert (before_word[0] == before_word[1]) != (before_word[1] == before_word[2])
    pair = 0 if before_word[0] == before_word[1] else 1
    assert mul(before_word[pair], before_word[pair+1], n) == 0
    replay = before.copy()
    write_fan(replay, 0, *primitive(*before_word, parent_channel, n), n, pe, pr)
    assert np.array_equal(replay, middle)
    assert activity(before, target, n, pe, pr) == 0
    wq = word(middle, target)
    target_channel = next(ch for ch in range(3, 15) if primitive(*wq, ch, n) != wq)
    after = middle.copy()
    write_fan(after, target, *primitive(*wq, target_channel, n), n, pe, pr)
    qb, qa = charges(before), charges(after)
    difference = qa-qb
    assert sorted(difference[difference != 0].tolist()) == [-1, 1]
    destination = int(np.where(difference == 1)[0][0])
    source = int(np.where(difference == -1)[0][0])
    assert source not in supports[0] and destination in supports[0]
    assert int(qb.sum()) == int(qm.sum()) == int(qa.sum())
    return {'n': n, 'seed': seed, 'geometry': geo, 'parent_channel': parent_channel,
            'target_channel': target_channel, 'flat_pair_faces_before': supports[0][pair:pair+2].tolist(),
            'shared_identity_face_between_events': int(shared), 'net_source_face': source,
            'net_destination_face': destination, 'conserved_charge': int(qb.sum()),
            'before_raw_links': before.tolist(), 'between_raw_links': middle.tolist(), 'after_raw_links': after.tolist(),
            'before_charge': qb.tolist(), 'between_charge': qm.tolist(), 'after_charge': qa.tolist(),
            'parent_word_before': list(before_word), 'parent_word_after': list(wp),
            'neighbor_word_before': list(word(before, target)), 'neighbor_word_between': list(wq),
            'neighbor_word_after': list(word(after, target))}


if __name__ == '__main__':
    result = {'scope': 'Exact Palm sampling of stationary marked reactions; no extra dynamics',
              'two_reaction_transport_witness': transport_witness(), 'runs': []}
    for side, trials in ((4, 65536), (8, 32768)):
        row = experiment(side, trials)
        result['runs'].append(row)
        Path('data/fiber-reaction-carriers.json').write_text(json.dumps(result, indent=2)+'\n')
        print(json.dumps({k: v for k, v in row.items()
                          if k in ('side', 'means', 'independent_trajectory_SE', 'conditional_total_rates_forward_backward', 'conditional_rate_SE')}), flush=True)
