#!/usr/bin/env python3
"""Exact discrete residence times and dependency-updated active primitive sets."""
import argparse
import hashlib
import json
import random
from collections import Counter
from fractions import Fraction
from pathlib import Path

from triangle_spreading import SpreadingProbe, digest, distances


def block_outcome(total, active, length, integer):
    """None means all `length` attempts were no-ops; otherwise first change 1..B.

    For uniform integer in [0,M**B), the masses are s**B for no change and
    K*s**(t-1)*M**(B-t) for first change at t, where s=M-K. No floats.
    """
    if not 1 <= active <= total or length < 1 or not 0 <= integer < total**length:
        raise ValueError('invalid exact geometric block')
    denominator, inactive = total**length, total-active
    survival = inactive**length
    if integer < survival:
        return None
    integer -= survival
    lo, hi = 1, length
    while lo < hi:
        middle = (lo+hi)//2
        cumulative = denominator-inactive**middle*total**(length-middle)
        if integer < cumulative:
            hi = middle
        else:
            lo = middle+1
    return lo


def geometric_wait(rng, total, active, horizon):
    """Exact Geom(K/M), right-censored at the given remaining attempt horizon."""
    if any(type(x) is not int for x in (total, active, horizon)) or total < 1 or not 0 <= active <= total or horizon < 0:
        raise ValueError('invalid discrete residence-time dimensions')
    if not active or not horizon:
        return None
    block = min(256, max(1, total//active))
    elapsed = 0
    while elapsed < horizon:
        length = min(block, horizon-elapsed)
        outcome = block_outcome(total, active, length, rng.randrange(total**length))
        if outcome is not None:
            return elapsed+outcome
        elapsed += length
    return None


class ActiveConnection:
    def __init__(self, probe, links):
        self.p, self.e, self.links = probe, probe.e, links[:]
        e = self.e
        self.supports = len(e.factor.pairs)
        self.total = len(e.tables)*self.supports
        self.codes = [[e.factor.oracle.code(links, arity, patch) for patch in range(self.supports)] for arity in (0, 1)]
        self.enabled = [[[0] if e.tables[0][code] != code else [] for code in range(36)],
                        [[r for r, table in enumerate(e.tables[1:], 1) if table[code] != code] for code in range(216)]]
        self.seam = [[bool(read & probe.seam_edges) for read in reads] for reads in probe.reads]
        self.active, self.position, self.seam_active = [], [-1]*self.total, 0
        for rule, table in enumerate(e.tables):
            arity = int(bool(rule))
            for patch, code in enumerate(self.codes[arity]):
                if table[code] != code:
                    self.add(rule*self.supports+patch)
        self.readers = [set() for _ in links]
        for arity, reads in enumerate(probe.reads):
            for patch, read in enumerate(reads):
                for edge in read:
                    self.readers[edge].add((arity, patch))
        self.q = probe.charge(links)
        self.charge = sum(self.q)
        self.full_audit()

    def add(self, op):
        if self.position[op] != -1:
            raise ValueError('active primitive inserted twice')
        self.position[op] = len(self.active)
        self.active.append(op)
        rule, patch = divmod(op, self.supports)
        self.seam_active += self.seam[int(bool(rule))][patch]

    def remove(self, op):
        i = self.position[op]
        if i == -1:
            raise ValueError('inactive primitive removed')
        last = self.active.pop()
        if i < len(self.active):
            self.active[i] = last
            self.position[last] = i
        self.position[op] = -1
        rule, patch = divmod(op, self.supports)
        self.seam_active -= self.seam[int(bool(rule))][patch]

    def full_audit(self):
        e = self.e
        expected = set()
        for arity in (0, 1):
            for patch in range(self.supports):
                code = e.factor.oracle.code(self.links, arity, patch)
                if code != self.codes[arity][patch]:
                    raise ValueError('cached based tuple differs from full raw reconstruction')
                expected.update(rule*self.supports+patch for rule in self.enabled[arity][code])
        if expected != set(self.active) or len(expected) != len(self.active):
            raise ValueError('active primitive set is incomplete or duplicated')
        if any(self.position[op] != i for i, op in enumerate(self.active)) or {op for op, i in enumerate(self.position) if i != -1} != expected:
            raise ValueError('active primitive reverse index is inconsistent')
        if self.q != self.p.charge(self.links) or sum(self.q) != self.charge:
            raise ValueError('cached charges differ from full raw holonomies')
        if len(self.active) != self.p.measure(self.links)['changing_operators']:
            raise ValueError('active primitive count differs from independent star formula')
        if self.seam_active != sum(self.seam[int(bool(op//self.supports))][op % self.supports] for op in self.active):
            raise ValueError('active boundary-support count is inconsistent')

    def execute(self, op):
        if type(op) is not int or not 0 <= op < self.total or self.position[op] == -1:
            raise ValueError('selected primitive is not active')
        e, p = self.e, self.p
        rule, patch = divmod(op, self.supports)
        arity = int(bool(rule))
        code, target = self.codes[arity][patch], e.tables[rule][self.codes[arity][patch]]
        writes = p.writes[arity][patch]
        old = {edge: self.links[edge] for edge in writes}
        if e.factor.oracle.update(self.links, rule, patch, e.tables[rule]) != (code, target):
            raise ValueError('cached choice differs from independent raw primitive')
        changed = {edge for edge in writes if old[edge] != self.links[edge]}
        if not changed:
            raise ValueError('selected active primitive did not change raw links')
        touched_faces = {f for edge in changed for f in p.incident[edge]}
        before_charge = sum(self.q[f] for f in touched_faces)
        for f in touched_faces:
            self.q[f] = e.charges[e.factor.oracle.fan.transport(self.links, e.factor.oracle.fan.face_paths[f])]
        if sum(self.q[f] for f in touched_faces) != before_charge:
            raise ValueError('local event changed total charge')
        affected = sorted({support for edge in changed for support in self.readers[edge]})
        for kind, index in affected:
            before = self.codes[kind][index]
            after = e.factor.oracle.code(self.links, kind, index)
            if before == after:
                continue
            for r in self.enabled[kind][before]:
                self.remove(r*self.supports+index)
            self.codes[kind][index] = after
            for r in self.enabled[kind][after]:
                self.add(r*self.supports+index)
        return [rule, patch, code, target], touched_faces, len(affected)


def relative_link(probe, edge):
    side = probe.e.geometry.side
    return [[v % side-side//2, v//side-side//2] for v in probe.e.geometry.edges[edge]]


def difference_digest(probe, links, background):
    return digest([[relative_link(probe, edge), value] for edge, value in enumerate(links) if value != background[edge]])


def trial(probe, trial_id, radii=(2, 4, 6, 8, 10), horizon=None, replay=False, smaller=None):
    side = probe.e.geometry.side
    horizon = horizon if horizon is not None else 12000*(side//6)**2
    background, seeds = probe.initial_conditions()
    seed = seeds[0]
    state = ActiveConnection(probe, seed['links'])
    distance = distances(probe.dual_neighbors, probe.incident[seed['edge']])
    jump_seed, clock_seed = 443000+trial_id, 844000+trial_id
    jump_rng, clock_rng = random.Random(jump_seed), random.Random(clock_seed)
    tick, residence, events = 0, Fraction(), []
    first = [None]*probe.faces
    passage = {r: None for r in radii}
    touched_counts = Counter()
    stopping = 'attempt_horizon'
    censored_at = horizon
    spatial_hash = hashlib.sha256()
    prefix = None
    for event_number in range(20000):
        if state.seam_active:
            stopping, censored_at = 'active_support_reaches_seam', tick
            break
        k = len(state.active)
        if not k:
            raise ValueError('active component entered a fixed state')
        wait = geometric_wait(clock_rng, state.total, k, horizon-tick)
        if wait is None:
            break
        next_tick = tick+wait
        op = state.active[jump_rng.randrange(k)]
        rule, patch = divmod(op, state.supports)
        if probe.reads[int(bool(rule))][patch] & probe.seam_edges:
            raise ValueError('seam-active count missed a boundary operator')
        tick = next_tick
        residence += Fraction(1, k)
        event, affected, count = state.execute(op)
        events.append([tick, k, *event])
        paths = probe.e.factor.oracle.fan.patches[patch] if rule else probe.e.factor.oracle.pairs[patch]
        description = [k, rule, [[[relative_link(probe, edge), reverse] for edge, reverse in path] for path in paths], *event[2:]]
        spatial_hash.update((json.dumps(description, separators=(',', ':'))+'\n').encode())
        if smaller is not None and len(events) == smaller['event_count']:
            prefix = {'events': len(events),
                      'same_spatial_event_digest': spatial_hash.hexdigest() == smaller['spatial_events_sha256'],
                      'same_translated_raw_difference': difference_digest(probe, state.links, background) == smaller['final_relative_difference_sha256']}
        touched_counts[count] += 1
        for f in affected:
            if state.q[f] != 1 and first[f] is None:
                first[f] = tick
                for r in radii:
                    if distance[f] >= r and passage[r] is None:
                        passage[r] = {'tick': tick, 'event_number': len(events),
                                      'conditional_mean_normalized_time': [residence.numerator, residence.denominator]}
        if len(events) % 256 == 0:
            state.full_audit()
    else:
        stopping, censored_at = 'event_budget', tick
    state.full_audit()
    verification = None
    if replay and events:
        schedule = [rule*state.supports+patch for _, _, rule, patch, _, _ in events]
        # No-op attempts are irrelevant to raw replay; residence times are
        # independently certified by the exact geometric law.
        result = probe.e.compiled_bank([seed['links']], schedule, 1)
        full = result['runs'][1]
        expected = [[i+1, *row[2:]] for i, row in enumerate(events)]
        if full['events'] != expected or full['final_links'] != state.links:
            raise ValueError('event-driven history differs from C++ raw-link replay')
        verification = {'all_changed_events_replayed_by_cpp': True, 'exact_inverse': full['exact_link_inverse']}
    return {'side': side, 'trial': trial_id, 'jump_seed': jump_seed, 'clock_seed': clock_seed,
            'attempt_horizon': horizon, 'attempted_operators': state.total,
            'stopping': stopping, 'censored_at_tick': censored_at,
            'passages': [[r, passage[r]] for r in radii], 'event_count': len(events),
            'events_sha256': digest(events), 'first_charge_departures': first,
            'spatial_events_sha256': spatial_hash.hexdigest(),
            'final_relative_difference_sha256': difference_digest(probe, state.links, background),
            'smaller_mesh_prefix_check': prefix,
            'read_supports_refreshed_histogram': sorted(touched_counts.items()),
            'final_links': state.links, 'final_activity': probe.measure(state.links),
            'compiled_event_replay': verification}


def summarize(records, radii):
    rows = []
    for side in sorted({r['side'] for r in records}):
        selected = [r for r in records if r['side'] == side]
        for radius in radii:
            values = [dict(r['passages'])[radius] for r in selected]
            observed = [v for v in values if v is not None]
            # No complete-case average when boundary stopping censors a trial.
            # Censoring is informative, so no Kaplan-Meier independence claim.
            row = {'side': side, 'radius': radius, 'trials': len(values), 'observed': len(observed),
                   'censored': len(values)-len(observed)}
            if len(observed) == len(values):
                xs = [Fraction(v['tick'], selected[0]['attempted_operators']) for v in observed]
                mean = sum(xs)/len(xs)
                variance = sum((x-mean)**2 for x in xs)/(len(xs)-1) if len(xs) > 1 else Fraction()
                rb = sum(Fraction(*v['conditional_mean_normalized_time']) for v in observed)/len(observed)
                row.update(mean_normalized_time=[mean.numerator, mean.denominator],
                           sample_variance_normalized_time=[variance.numerator, variance.denominator],
                           mean_conditional_residence=[rb.numerator, rb.denominator])
            rows.append(row)
    return rows


def audit(sides=(12, 24), trials=32):
    if not sides or any(s not in (6, 12, 24) for s in sides) or len(set(sides)) != len(sides) or not 1 <= trials <= 128:
        raise ValueError('invalid bounded first-passage audit')
    bank = json.loads(Path('data/triangle-feedback.json').read_text())['bank']
    radii, records = (2, 4, 6, 8, 10), []
    for side in sides:
        p = SpreadingProbe(side, bank)
        for i in range(trials):
            smaller = next((r for r in records if r['trial'] == i and r['side'] == side//2), None)
            row = trial(p, i, radii, replay=i < 2, smaller=smaller)
            records.append(row)
            print('first passage', side, i, row['event_count'], row['stopping'],
                  'observed', [r for r, value in row['passages'] if value is not None], flush=True)
    return {'records': records, 'summary': summarize(records, radii),
            'scope': 'Independent trial seeds within each size; matching seeds across sizes are paired. Exact original discrete uniform-attempt law, stopped when any active support reaches the periodic seam. Informative censoring; no fitted speed exponent or emergent spacetime claim.'}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--sides', type=int, nargs='+', default=[12, 24])
    parser.add_argument('--trials', type=int, default=32)
    parser.add_argument('--output', type=Path, default=Path('out/triangle-first-passage.json'))
    args = parser.parse_args()
    data = audit(tuple(args.sides), args.trials)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(data, separators=(',', ':'))+'\n')


if __name__ == '__main__':
    main()
