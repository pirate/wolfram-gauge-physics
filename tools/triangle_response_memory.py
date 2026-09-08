#!/usr/bin/env python3
"""Exact charge-observable transition counts, including relative-gauge activity.

This is an observer of the existing raw-link rules, not a replacement evolution
on charge fields: the activity mask is read from actual based holonomies.
"""
import argparse
import itertools
import json
import random
from collections import Counter
from fractions import Fraction
from pathlib import Path

from triangle_charge_response import ResponseGeometry


def activity_mask(geometry, links, q=None):
    q = geometry.p.charge(links) if q is None else q
    active = []
    e, n = geometry.e, geometry.e.geometry.group.n
    for patch, faces in enumerate(geometry.p.supports[1]):
        if all(q[f] == 1 for f in faces):
            code = e.factor.oracle.code(links, 1, patch)
            a, b, c = code//n**2, code//n % n, code % n
            if (a == b) != (b == c):
                active.append(patch)
    return tuple(active)


def joint_targets(geometry, links, compiled=False):
    """Full one-attempt (q,A) target law, including no-ops and multiplicities."""
    e, p = geometry.e, geometry.p
    q = tuple(p.charge(links))
    initial = (q, activity_mask(geometry, links, q))
    targets, cache = Counter(), {}
    changing = 0
    for rule, table in enumerate(e.tables):
        for patch in range(len(e.factor.pairs)):
            code = e.factor.oracle.code(links, rule, patch)
            if code == table[code]:
                continue
            moved = links[:]
            e.factor.oracle.update(moved, rule, patch, table)
            key = tuple(moved)
            if key not in cache:
                after = tuple(p.charge(moved))
                cache[key] = (after, activity_mask(geometry, moved, after))
            targets[cache[key]] += 1
            changing += 1
    targets[initial] += 39*geometry.faces-changing
    if compiled:
        scan = [op for op in range(39*geometry.faces) for _ in range(2)]
        run = e.compiled_bank([links], scan, len(scan), capture_links=True)['runs'][1]
        observed = Counter()
        for state in run['raw_event_links'][::2]:
            after = tuple(p.charge(state))
            observed[after, activity_mask(geometry, state, after)] += 1
        observed[initial] += 39*geometry.faces-len(run['raw_event_links'])//2
        if observed != targets:
            raise ValueError('C++ joint charge/activity targets disagree')
    return targets


def observed_charge_targets(geometry, q, active):
    """Charge generator given actual observed (q,A), not autonomous evolution."""
    p, e = geometry.p, geometry.e
    q, active = list(q), set(active)
    targets = Counter()
    for a, b in e.factor.pairs:
        if (q[a] == 0) != (q[b] == 0):
            moved = q[:]
            moved[a], moved[b] = q[b], q[a]
            targets[tuple(moved)] += 1
    for patch, faces in enumerate(p.supports[1]):
        values = [q[f] for f in faces]
        if sorted(values) == [0, 1, 2]:
            moved = q[:]
            for f in faces:
                moved[f] = 1
            targets[tuple(moved)] += 4
        elif values == [1, 1, 1]:
            if patch in active:
                for assignment in itertools.permutations((0, 1, 2)):
                    moved = q[:]
                    for f, v in zip(faces, assignment):
                        moved[f] = v
                    targets[tuple(moved)] += 2
    return targets


def charge_targets(geometry, links):
    """Nontrivial charge targets with original rooted-operator multiplicities."""
    q = geometry.p.charge(links)
    active = activity_mask(geometry, links, q)
    return observed_charge_targets(geometry, q, active), list(active)


def second_drift(geometry, q, active):
    d = geometry.drift(q)
    result = [0]*geometry.faces
    for target, count in observed_charge_targets(geometry, q, active).items():
        after = geometry.drift(target)
        result = [x+count*(a-b) for x, a, b in zip(result, after, d)]
    return result


def raw_targets(geometry, links, compiled=False):
    """Independent enumeration through raw connection surgery, optionally C++."""
    e, p = geometry.e, geometry.p
    initial = tuple(p.charge(links))
    targets = Counter()
    for rule, table in enumerate(e.tables):
        for patch in range(len(e.factor.pairs)):
            code = e.factor.oracle.code(links, rule, patch)
            if code == table[code]:
                continue
            moved = links[:]
            e.factor.oracle.update(moved, rule, patch, table)
            target = tuple(p.charge(moved))
            if target != initial:
                targets[target] += 1
    if compiled:
        # An involution followed by itself resets the initial state. This scans
        # every operator, rather than only a seeded trajectory's selected ones.
        scan = [op for op in range(13*len(e.factor.pairs)) for _ in range(2)]
        run = e.compiled_bank([links], scan, len(scan), capture_links=True)['runs'][1]
        observed = Counter(tuple(p.charge(s)) for s in run['raw_event_links'][::2])
        observed.pop(initial, None)
        if observed != targets:
            raise ValueError('C++ raw targets disagree with microscopic charge counts')
    return targets


def response(geometry, links, verify_raw=False):
    q = geometry.p.charge(links)
    d = geometry.drift(q)
    targets, active = charge_targets(geometry, links)
    if verify_raw and targets != raw_targets(geometry, links, compiled=True):
        raise ValueError('charge/activity generator disagrees with raw transitions')
    first, second = [0]*geometry.faces, [0]*geometry.faces
    for target, count in targets.items():
        after = geometry.drift(target)
        for f in range(geometry.faces):
            first[f] += count*(target[f]-q[f])
            second[f] += count*(after[f]-d[f])
    if first != d or sum(second):
        raise ValueError('response hierarchy fails first drift or conservation')
    m = 39*geometry.faces
    return {'charge_field': q, 'charge_drift_sum': d,
            'drift_of_drift_sum': second, 'active_reflection_fans': active,
            'attempted_operators': m,
            'two_step_mean_charge_numerator': [m*m*x+2*m*y+z for x, y, z in zip(q, d, second)],
            'two_step_mean_charge_denominator': m*m,
            'distinct_charge_targets': len(targets),
            'changing_charge_operators': sum(targets.values())}


def audit(compiled=True):
    source = json.loads(Path('data/triangle-charge-current.json').read_text())['readout']
    bank = json.loads(Path('data/triangle-feedback.json').read_text())['bank']
    geometry = ResponseGeometry(12, bank)
    ids = source['same_boundary_witness_states']
    links = [source['records'][i]['links'] for i in ids]
    records = [response(geometry, state, verify_raw=compiled) for state in links]
    a, b = records
    if a['charge_field'] != b['charge_field'] or a['charge_drift_sum'] != b['charge_drift_sum']:
        raise ValueError('witness does not hold charge and first drift fixed')
    difference = [[f, x-y] for f, (x, y) in enumerate(zip(a['drift_of_drift_sum'], b['drift_of_drift_sum'])) if x != y]
    if not difference:
        raise ValueError('prepared states do not distinguish two-step mean response')
    gauge_checks = []
    e, g = geometry.e, geometry.e.geometry.group
    for i, (state, expected) in enumerate(zip(links, records)):
        rng = random.Random(219440+i)
        frames = [rng.randrange(g.n) for _ in range(e.geometry.size)]
        moved = [g.mul[frames[v]][g.mul[x][g.inv[frames[u]]]] for (u, v), x in zip(e.geometry.edges, state)]
        if response(geometry, moved, verify_raw=compiled) != expected:
            raise ValueError('local frame transformation changes second response')
        gauge_checks.append({'frames': frames, 'same_response': True})
    return {'side': 12, 'source_states': ids, 'records': records,
            'drift_of_drift_difference': difference, 'gauge_checks': gauge_checks,
            'activity_closure_witness': activity_closure_audit(bank, compiled),
            'scope': 'Charge plus drift fails at two-attempt mean response. Charge plus actual based-fan activity suffices for one-step charge transition counts but is not a closed joint observer, witnessed by exact raw transitions on an 18-face mesh.'}


def activity_closure_audit(bank, compiled=True):
    """Bounded reproducible search; exact operator census certifies its witness."""
    from triangle_reference import AxialReference, fraction
    geometry = ResponseGeometry(3, bank)
    sampler = AxialReference(3, bank, 4)
    seed, limit = 426880, 1000
    rng, seen = random.Random(seed), {}
    for index in range(limit):
        links = sampler.sample(rng, 'nonabelian_reflections')['links']
        q = tuple(geometry.p.charge(links))
        key = (q, activity_mask(geometry, links, q))
        if key not in seen:
            seen[key] = (index, links)
            continue
        old_index, old = seen[key]
        left, right = joint_targets(geometry, old), joint_targets(geometry, links)
        if left == right:
            continue
        if compiled and (left != joint_targets(geometry, old, True) or right != joint_targets(geometry, links, True)):
            raise ValueError('compiled joint witness differs')
        difference = [(target, left[target]-right[target]) for target in sorted(set(left) | set(right)) if left[target] != right[target]]
        # Marginalizing A must recover the same next-charge law, exactly.
        marginal = Counter()
        for (charge, _), count in difference:
            marginal[charge] += count
        if any(marginal.values()) or sum(count for _, count in difference):
            raise ValueError('activity witness does not preserve the next-charge marginal')
        m = 39*geometry.faces
        third = [0]*geometry.faces
        for (charge, active), count in difference:
            e = second_drift(geometry, charge, active)
            third = [x+count*y for x, y in zip(third, e)]
        if sum(third):
            raise ValueError('third response violates conservation')
        # The original pair has the same q,A, hence the same D and L D.
        # Delta(P^3 q) = Delta(L^3 q)/M^3 = Delta(L^2 D)/M^3.
        return {'side': 3, 'charge': 4, 'search_seed': seed, 'search_limit': limit,
                'sample_indices': [old_index, index], 'links': [old, links],
                'same_charge_field': list(q), 'same_active_reflection_fans': list(key[1]),
                'joint_target_count_difference': [[list(charge), list(active), count] for (charge, active), count in difference],
                'joint_transition_total_variation': fraction(Fraction(sum(abs(count) for _, count in difference), 2*m)),
                'three_step_mean_difference_numerator': third,
                'three_step_mean_difference_denominator': m**3}
    raise ValueError('bounded activity-closure search found no witness')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, default=Path('out/triangle-response-memory.json'))
    args = parser.parse_args()
    data = audit()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(data, separators=(',', ':'))+'\n')


if __name__ == '__main__':
    main()
