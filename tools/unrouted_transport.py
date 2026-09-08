#!/usr/bin/env python3
"""Detect gauge memory under a state-independent scheduler, without routed circuits."""
import argparse
import itertools
import json
import random
from pathlib import Path

from fiber_transport import FiberTransport
from transport_braid import spectral_word_witness


def detect_return(e, initial, positions, seed, attempts):
    """First gauge-changing return to the most recent identical tagged placement.

    Tags are passive provenance labels; they never enter the scheduler or link rule.
    Retaining only the latest visit is an explicit detection policy, not an exhaustive
    search over all pairs of times. A miss is a bounded negative, never a no-go result.
    """
    if not 0 <= attempts <= 100000:
        raise ValueError('return detector supports at most 100000 attempts')
    g = e.geometry.group
    if len(set(positions)) != len(positions) or set(positions) != {i for i, x in enumerate(e.geometry.holonomies(initial)) if x != g.identity}:
        raise ValueError('tags must cover each initial nonflat face exactly once')
    links, placement = initial[:], list(positions)
    tag_at = {face: tag for tag, face in enumerate(positions)}
    rng = random.Random(seed)
    signature = e.forest.signature(links)
    seen = {tuple(placement): (0, links[:], signature)}
    schedule, changed, witness = [], 0, None
    for tick in range(1, attempts+1):
        patch = rng.randrange(len(e.factor.pairs))
        schedule.append(patch)
        code, target = e.factor.oracle.update(links, 0, patch, e.factor.tables[0])
        if code == target:
            continue
        u, v = e.factor.pairs[patch]
        if (u in tag_at) == (v in tag_at):
            raise ValueError('tag propagation disagrees with the vacancy update')
        if v in tag_at:
            u, v = v, u
        tag = tag_at.pop(u)
        tag_at[v] = tag
        placement[tag] = v
        changed += 1
        key = tuple(placement)
        signature = e.forest.signature(links)
        previous = seen.get(key)
        if previous is not None and signature != previous[2]:
            start, before, _ = previous
            witness = {'start_tick': start, 'end_tick': tick, 'positions': placement[:],
                       'before_links': before, 'after_links': links[:]}
            break
        seen[key] = (tick, links[:], signature)
    return {'seed': seed, 'attempt_budget': attempts, 'attempts_used': len(schedule),
            'changing_updates': changed, 'distinct_tagged_placements': len(seen),
            'schedule': schedule, 'final_links': links, 'witness': witness}


def worldline_audit(e, initial_positions, schedule):
    """Independent occupancy-only replay, with lifted torus displacements.

    Integer coordinates are three times the face barycenter. Minimum periodic
    displacements of adjacent dual faces are unambiguous for sides >= 3.
    Nonzero winding diagnoses a torus contribution; zero abelian winding does
    not prove a path is contractible in a punctured torus.
    """
    side = e.geometry.side
    period = 3*side
    def point(f):
        y, x = divmod(f//2, side)
        return (3*x+2-f % 2, 3*y+1+f % 2)
    placement = list(initial_positions)
    at = {face: tag for tag, face in enumerate(placement)}
    displacements = [[0, 0] for _ in placement]
    events = []
    for tick, patch in enumerate(schedule, 1):
        u, v = e.factor.pairs[patch]
        if (u in at) == (v in at):
            continue
        if v in at:
            u, v = v, u
        tag = at.pop(u); at[v] = tag; placement[tag] = v
        raw = [b-a for a, b in zip(point(u), point(v))]
        delta = [(x+period//2) % period-period//2 for x in raw]
        if any(abs(x) > 2 for x in delta):
            raise ValueError('worldline crosses nonadjacent dual faces')
        for j in range(2):
            displacements[tag][j] += delta[j]
        events.append({'tick': tick, 'patch': patch, 'tag': tag, 'source': u, 'target': v,
                       'displacement': delta, 'crosses_periodic_seam': raw != delta})
    closed = placement == list(initial_positions)
    if closed and any(x % period for d in displacements for x in d):
        raise ValueError('closed torus walk has nonintegral winding')
    return {'events': events, 'final_positions': placement, 'all_tags_return': closed,
            'displacements': displacements,
            'torus_winding_by_tag': [[x//period for x in d] for d in displacements] if closed else None,
            'seam_crossing_moves': sum(x['crosses_periodic_seam'] for x in events)}


def commuting_return_audit(e, before, after, reflection, worldlines):
    """Z_2 control: closed dual worldlines pair with primal torus periods.

    A vacancy move toggles its common edge by the sole nontrivial subgroup
    element. Once every tagged face returns, the difference is a flat Z_2
    cochain. It is a vertex-frame change iff its two torus periods vanish.
    """
    g, side = e.geometry.group, e.geometry.side
    if reflection == g.identity or g.mul[reflection][reflection] != g.identity:
        raise ValueError('commuting recurrence audit requires a nontrivial involution')
    if not worldlines['all_tags_return'] or any(x not in (g.identity, reflection) for x in before+after):
        raise ValueError('expected a closed return within one reflection subgroup')
    difference = [g.mul[y][g.inv[x]] for x, y in zip(before, after)]
    if any(x != g.identity for x in e.geometry.holonomies(difference)):
        raise ValueError('closed commuting return has a nonflat difference cochain')
    walks = [list(range(side))+[0], list(range(0, side*side, side))+[0]]
    periods = []
    for walk in walks:
        value = g.identity
        for u, v in zip(walk, walk[1:]):
            x = difference[e.geometry.ids[tuple(sorted((u, v)))]]
            value = g.mul[g.inv[x] if u > v else x][value]
        periods.append(int(value == reflection))
    winding = [sum(d[j] for d in worldlines['torus_winding_by_tag']) for j in range(2)]
    predicted = [winding[1] % 2, winding[0] % 2]
    if periods != predicted:
        raise ValueError('commuting loop periods disagree with passive-worldline intersection parity')
    equivalent = e.frame_witness(before, after) is not None
    if equivalent != (periods == [0, 0]):
        raise ValueError('commuting gauge return disagrees with the two-period criterion')
    return {'reflection': reflection, 'difference_links': difference, 'difference_flat': True,
            'primal_period_walks': walks, 'measured_period_bits': periods,
            'total_dual_winding': winding, 'predicted_period_bits': predicted,
            'gauge_equivalent_iff_both_periods_zero': True}


def experiment(side, seed, attempts):
    e = FiberTransport(side, 3, [(0, 1), (1, 2), (2, 0)])
    g = e.geometry.group
    a, b = next((a, b) for a, b in itertools.product(range(g.n), repeat=2)
                if a != g.identity and b != g.identity and g.mul[a][a] == g.identity and
                g.mul[b][b] == g.identity and g.mul[a][b] != g.mul[b][a])
    initial, positions = e.geometry.paired_seed(a, b)
    result = detect_return(e, initial, positions, seed, attempts)
    result.update(side=side, initial_links=initial, initial_positions=positions,
                  derived_automorphisms=g.elements, seed_generators=[a, b])
    # Even a bounded negative gets independent C++ replay of its entire history.
    result['full_history_check'] = e.compiled(initial, result['schedule'])
    if result['full_history_check']['final_links'] != result['final_links']:
        raise ValueError('compiled full history differs from the online detector')
    whole = worldline_audit(e, positions, result['schedule'])
    result['occupancy_only_changing_updates'] = len(whole['events'])
    if len(whole['events']) != result['changing_updates']:
        raise ValueError('unrouted occupancy-only clock disagrees with primitive evolution')
    witness = result['witness']
    if witness is None:
        return result
    begin, end = witness['start_tick'], witness['end_tick']
    prefix, segment = result['schedule'][:begin], result['schedule'][begin:end]
    witness['worldlines'] = worldline_audit(e, witness['positions'], segment)
    if not witness['worldlines']['all_tags_return']:
        raise ValueError('discovered witness is not a tagged return')
    comparisons = []
    for name, generators in (('noncommuting', (a, b)), ('commuting', (a, a))):
        start, control_positions = e.geometry.paired_seed(*generators)
        if control_positions != positions:
            raise ValueError('control seed has different tagged placement')
        before_check = e.compiled(start, prefix)
        before = before_check['final_links']
        after_check = e.compiled(before, segment)
        after = after_check['final_links']
        if e.geometry.sectors(before) != e.geometry.sectors(after):
            raise ValueError('tagged recurrence changed the face-class field')
        frames = e.frame_witness(before, after)
        if (frames is not None) != (before_check['final_signature'] == after_check['final_signature']):
            raise ValueError('unrouted gauge-equivalence checks disagree')
        spectrum = spectral_word_witness(e.geometry, [{'links': before}, {'links': after}], full_bundle=True)
        comparisons.append({'condition': name, 'prefix_check': before_check, 'segment_check': after_check,
                            'gauge_equivalent': frames is not None, 'frames_between_endpoints': frames,
                            'same_face_classes': True, 'exact_full_bundle_spectrum': spectrum})
    if comparisons[0]['gauge_equivalent'] or comparisons[0]['prefix_check']['final_links'] != witness['before_links'] or comparisons[0]['segment_check']['final_links'] != witness['after_links']:
        raise ValueError('discovered raw witness failed independent compiled replay')
    # Both controls have identical occupancy evolution, including which attempts act.
    if [x[:2] for x in comparisons[0]['segment_check']['events']] != [x[:2] for x in comparisons[1]['segment_check']['events']]:
        raise ValueError('same-occupancy control changed the update clock')
    witness['comparisons'] = comparisons
    witness['commuting_topology_audit'] = commuting_return_audit(
        e, comparisons[1]['prefix_check']['final_links'], comparisons[1]['segment_check']['final_links'], a, witness['worldlines'])
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--attempts', type=int, default=10000)
    parser.add_argument('--sides', type=int, nargs='+', default=[3, 4, 6])
    parser.add_argument('--seeds', type=int, nargs='+', default=[0, 1])
    parser.add_argument('--output', type=Path, default=Path('out/unrouted-transport.json'))
    args = parser.parse_args()
    result = {'scope': 'First gauge-changing most-recent tagged recurrence under uniform state-independent rooted-patch sampling. Supplied torus and classical vacancy dynamics, not autonomous binding or emergent geometry.', 'runs': []}
    for side, seed in itertools.product(args.sides, args.seeds):
        row = experiment(side, seed, args.attempts)
        result['runs'].append(row)
        w = row['witness']
        print('side', side, 'seed', seed, 'attempts', row['attempts_used'],
              'return', None if w is None else (w['start_tick'], w['end_tick']),
              'control returns', None if w is None else [c['gauge_equivalent'] for c in w['comparisons']], flush=True)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, separators=(',', ':'))+'\n')


if __name__ == '__main__':
    main()
