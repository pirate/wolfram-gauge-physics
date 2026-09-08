#!/usr/bin/env python3
"""Seeded activity spreading on actual connections, with exact damage ancestry.

The one-link perturbation is supplied as an initial condition, not generated from
the frozen background by an allowed move. No propagation equation is inserted.
"""
import argparse
import hashlib
import json
import random
from collections import deque
from fractions import Fraction
from pathlib import Path

from triangle_reference import ActivityProbe, constant_connection, subgroup


def digest(value):
    return hashlib.sha256(json.dumps(value, separators=(',', ':')).encode()).hexdigest()


def distances(adjacent, sources):
    result = [None]*len(adjacent)
    queue = deque()
    for source in sources:
        result[source] = 0
        queue.append(source)
    while queue:
        u = queue.popleft()
        for v in sorted(adjacent[u]):
            if result[v] is None:
                result[v] = result[u]+1
                queue.append(v)
    if any(d is None for d in result):
        raise ValueError('spreading diagnostic requires a connected support graph')
    return result


class SpreadingProbe(ActivityProbe):
    def __init__(self, side, bank):
        super().__init__(side, bank)
        e, oracle = self.e, self.e.factor.oracle
        self.reads = [[{edge for path in pair for edge, _ in path} for pair in oracle.pairs], oracle.fan.reads]
        self.writes = [[{pair[0][0][0]} for pair in oracle.pairs], oracle.fan.writes]
        self.incident = [[] for _ in e.geometry.edges]
        for f, path in enumerate(oracle.fan.face_paths):
            for edge, _ in path:
                self.incident[edge].append(f)
        self.dual_neighbors = [set() for _ in range(self.faces)]
        for a, b in self.incident:
            self.dual_neighbors[a].add(b); self.dual_neighbors[b].add(a)
        # Directed links-to-links graph of possible primitive dependencies. This
        # is NOT a Euclidean distance or an inferred physical light cone.
        self.dependency_neighbors = [set() for _ in e.geometry.edges]
        for reads, writes in zip(self.reads, self.writes):
            for read, write in zip(reads, writes):
                for edge in read:
                    self.dependency_neighbors[edge].update(write)
        self.seam_edges = {i for i, (u, v) in enumerate(e.geometry.edges)
                           if abs(u % side-v % side) > 1 or abs(u//side-v//side) > 1}

    def initial_conditions(self):
        e, side = self.e, self.e.geometry.side
        background = constant_connection(self, (1, 2, 5))
        u = (side//2)*side+side//2
        seeds = []
        reflections = [h for h, q in enumerate(e.charges) if q == 1]
        for direction, v in zip(('horizontal', 'vertical', 'diagonal'), (u+1, u+side, u+side+1)):
            edge = e.geometry.ids[tuple(sorted((u, v)))]
            for value in reflections:
                if value == background[edge]:
                    continue
                links = background[:]; links[edge] = value
                m = self.measure(links)
                if m['populations'] != [0, self.faces, 0] or m['changing_operators'] != 96:
                    raise ValueError('single-link seed fails its homogeneous charge/activity certificate')
                if len(subgroup(e.geometry.group, e.forest.based_loops(links)[0])) != 6:
                    raise ValueError('seed does not preserve full nonabelian holonomy')
                seeds.append({'direction': direction, 'edge': edge,
                              'old': background[edge], 'new': value, 'links': links})
        if not self.measure(background)['frozen']:
            raise ValueError('unperturbed reference is not globally frozen')
        return background, seeds

    def analyze(self, background, initial, seed_edge, schedule, run, stride):
        e, oracle = self.e, self.e.factor.oracle
        if len(background) != len(e.geometry.edges) or len(initial) != len(background):
            raise ValueError('invalid coupled initial connection sizes')
        if [i for i, (a, b) in enumerate(zip(background, initial)) if a != b] != [seed_edge]:
            raise ValueError('coupled initial states must differ at exactly the seed edge')
        if self.charge(background) != [1]*self.faces or self.charge(initial) != [1]*self.faces:
            raise ValueError('spreading comparison requires identical homogeneous charge fields')
        if not self.measure(background)['frozen']:
            raise ValueError('damage ancestry requires the specified background to be fixed')
        face_distance = distances(self.dual_neighbors, self.incident[seed_edge])
        dependency_distance = distances(self.dependency_neighbors, [seed_edge])
        links, q = initial[:], self.charge(initial)
        # A non-null depth belongs only to a currently different raw link.
        # Identically reframing both runs preserves every equality comparison.
        depth = [None]*len(links); depth[seed_edge] = 0
        first_departure = [None]*self.faces
        first_departure_depth = [None]*self.faces
        snapshots, max_depth, index, first_seam = [], 0, 0, None
        events = run['events']
        previous_tick = 0
        for event in events:
            tick, rule, patch, _, _ = event
            if not previous_tick < tick <= len(schedule) or schedule[tick-1] != rule*len(e.factor.pairs)+patch:
                raise ValueError('event tick does not match the supplied primitive schedule')
            previous_tick = tick
        for checkpoint in range(0, len(schedule)+1, stride):
            while index < len(events) and events[index][0] <= checkpoint:
                tick, rule, patch, before, target = events[index]
                arity = int(bool(rule))
                read, write = self.reads[arity][patch], self.writes[arity][patch]
                if first_seam is None and read & self.seam_edges:
                    first_seam = tick
                parents = [depth[edge] for edge in read if links[edge] != background[edge]]
                if not parents or any(parent is None for parent in parents):
                    raise ValueError('changing event lacks a currently damaged read dependency')
                # A support identical to the frozen run must have an identity
                # action. Also check its actual table entry independently.
                base_code = oracle.code(background, rule, patch)
                if e.tables[rule][base_code] != base_code:
                    raise ValueError('claimed fixed background has an enabled operator')
                if oracle.update(links, rule, patch, e.tables[rule]) != (before, target):
                    raise ValueError('event differs from independent primitive replay')
                child_depth = max(parents)+1
                for edge in write:
                    depth[edge] = child_depth if links[edge] != background[edge] else None
                    if depth[edge] is not None and dependency_distance[edge] > depth[edge]:
                        raise ValueError('damage crossed its exact primitive-dependency bound')
                max_depth = max(max_depth, child_depth)
                for f in (e.fan_faces[patch] if rule else e.factor.pairs[patch]):
                    h = oracle.fan.transport(links, oracle.fan.face_paths[f])
                    q[f] = e.charges[h]
                    if q[f] != 1 and first_departure[f] is None:
                        first_departure[f] = tick
                        first_departure_depth[f] = child_depth
                index += 1
            if q != self.charge(links) or sum(q) != self.faces:
                raise ValueError('incremental charge field differs from full raw holonomies')
            if any((d is None) != (h == b) for d, h, b in zip(depth, links, background)):
                raise ValueError('current damage support differs from ancestry bookkeeping')
            m = self.measure(links)
            if m['frozen']:
                raise ValueError('invertible dynamics entered a globally fixed state')
            shells = []
            for radius in range(max(face_distance)+1):
                fs = [f for f, d in enumerate(face_distance) if d == radius]
                shells.append([radius, len(fs), *[sum(q[f] == a for f in fs) for a in range(3)],
                               sum(first_departure[f] is not None for f in fs)])
            current = [face_distance[f] for f in range(self.faces) if q[f] != 1]
            visited = [face_distance[f] for f in range(self.faces) if first_departure[f] is not None]
            snapshots.append({'tick': checkpoint, 'charge_field': q[:], 'activity': m,
                              'damaged_links': sum(d is not None for d in depth),
                              'events': index, 'maximum_damage_lineage_depth': max_depth,
                              'before_any_changing_support_touches_seam': first_seam is None,
                              'contrast_radius': max(current, default=None),
                              'ever_contrast_radius': max(visited, default=None),
                              'shells_radius_size_n0_n1_n2_ever': shells})
        if index != len(events) or links != run['final_links']:
            raise ValueError('damage replay did not reproduce the compiled final connection')
        return {'face_distances_from_seed': face_distance,
                'link_dependency_distances_from_seed': dependency_distance,
                'first_charge_departure_tick': first_departure,
                'first_charge_departure_lineage_depth': first_departure_depth,
                'first_changing_support_touching_periodic_seam': first_seam,
                'snapshots': snapshots, 'events_sha256': digest(events),
                'first_event': events[0] if events else None, 'event_count': len(events),
                'final_links': links,
                'each_changing_event_has_seed_ancestry': True,
                'independent_raw_link_replay': run['independent_raw_link_replay'],
                'exact_link_inverse': run['exact_link_inverse']}


def run_case(side, trial, bank, attempts=None):
    p = SpreadingProbe(side, bank)
    e, g = p.e, p.e.geometry.group
    background, seeds = p.initial_conditions()
    attempts = attempts if attempts is not None else 12000*(side//6)**2
    if attempts < 20 or attempts % 20 or attempts > 1000000:
        raise ValueError('attempts must be a multiple of 20 in [20, 1000000]')
    stride, schedule_seed = attempts//20, 76100+trial
    rng = random.Random(schedule_seed)
    schedule = [rng.randrange(13*len(e.factor.pairs)) for _ in range(attempts)]
    inputs = [background]+[seed['links'] for seed in seeds]
    compiled = e.compiled_bank(inputs, schedule, stride)
    records = []
    for condition, seed in enumerate(seeds, 1):
        control, full = compiled['runs'][2*condition:2*condition+2]
        if control['events'] or control['final_links'] != seed['links']:
            raise ValueError('reaction-disabled homogeneous-charge control moved')
        analyzed = p.analyze(background, seed['links'], seed['edge'], schedule, full, stride)
        # Before the first change, every attempted operator still sees the
        # initial connection; directly predict that event from the schedule.
        predicted = None
        for tick, op in enumerate(schedule, 1):
            rule, patch = divmod(op, len(e.factor.pairs))
            code = e.factor.oracle.code(seed['links'], rule, patch)
            if e.tables[rule][code] != code:
                predicted = [tick, rule, patch, code, e.tables[rule][code]]
                break
        if predicted != analyzed['first_event']:
            raise ValueError('first-event waiting time disagrees with exact initial activity')
        records.append({'seed': seed, 'evolution': analyzed})
        last = analyzed['snapshots'][-1]
        print('spread', side, trial, condition, 'events', analyzed['event_count'],
              'populations', last['activity']['populations'], 'ever radius', last['ever_contrast_radius'],
              'diameter from seed', max(analyzed['face_distances_from_seed']), flush=True)
    if any(run['events'] or run['final_links'] != background for run in compiled['runs'][:2]):
        raise ValueError('unperturbed frozen control moved')
    gauge_check = None
    if side == 6 and trial == 0:
        rng = random.Random(99877)
        frames = [rng.randrange(g.n) for _ in range(e.geometry.size)]
        def transform(links):
            return [g.mul[frames[v]][g.mul[h][g.inv[frames[u]]]] for (u, v), h in zip(e.geometry.edges, links)]
        framed = e.compiled_bank([transform(background), transform(seeds[0]['links'])], schedule, stride)
        original = compiled['runs'][3]
        moved = framed['runs'][3]
        if [ev[:3] for ev in original['events']] != [ev[:3] for ev in moved['events']] or transform(original['final_links']) != moved['final_links']:
            raise ValueError('coupled seed evolution does not commute with local frame changes')
        result = p.analyze(transform(background), transform(seeds[0]['links']), seeds[0]['edge'], schedule, moved, stride)
        for key in ('snapshots', 'first_charge_departure_tick', 'first_charge_departure_lineage_depth',
                    'first_changing_support_touching_periodic_seam'):
            if result[key] != records[0]['evolution'][key]:
                raise ValueError('spreading diagnostic is not invariant under a common local reframing')
        gauge_check = {'frames': frames, 'same_event_locations': True,
                       'final_raw_connections_transform_exactly': True,
                       'identical_spreading_observations': True}
    probability = Fraction(96, 13*len(e.factor.pairs))
    return {'side': side, 'trial': trial, 'attempts': attempts, 'stride': stride,
            'schedule_seed': schedule_seed, 'schedule_sha256': digest(schedule),
            'schedule_rule': 'Independent uniform integer draws from all rooted operators, including no-ops; the six seeds share this schedule.',
            'initial_background': background, 'initial_charge': p.faces,
            'initial_changing_operators_each_seed': 96,
            'first_change_probability_per_attempt': [probability.numerator, probability.denominator],
            'geometric_mean_first_change_attempt': [probability.denominator, probability.numerator],
            'attempted_operators': 13*len(e.factor.pairs),
            'frozen_and_reaction_disabled_controls_fixed': True,
            'runs': records, 'common_local_frame_control': gauge_check}


def audit(sides=(6, 12, 24), trials=2, attempts=None):
    if not sides or any(s not in (6, 12, 24) for s in sides) or len(set(sides)) != len(sides) or not 1 <= trials <= 8:
        raise ValueError('bounded spreading audit supports distinct sides 6, 12, 24 and 1..8 trials')
    bank = json.loads(Path('data/triangle-feedback.json').read_text())['bank']
    return {'cases': [run_case(side, trial, bank, attempts) for side in sides for trial in range(trials)],
            'scope': 'Externally seeded classical activity on supplied geometry and random schedules. Charge is unchanged by initialization. No spontaneous nucleation, emergent geometry, quantum signal, ballistic speed, or particle persistence claim.'}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--sides', type=int, nargs='+', default=[6, 12, 24])
    parser.add_argument('--trials', type=int, default=2)
    parser.add_argument('--attempts', type=int)
    parser.add_argument('--output', type=Path, default=Path('out/triangle-spreading.json'))
    args = parser.parse_args()
    data = audit(tuple(args.sides), args.trials, args.attempts)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(data, separators=(',', ':'))+'\n')


if __name__ == '__main__':
    main()
