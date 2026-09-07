#!/usr/bin/env python3
"""Exact abelian-subgroup face factor and sparse encounter reachability.

This is a derived factor of the link updates, not replacement dynamics. An
existential schedule is not a scattering probability or a typical trajectory.
"""
import argparse
import itertools
import json
import random
import subprocess
from collections import Counter, deque
from pathlib import Path

from face_energy_obstruction import derive_fiber_group
from run_three_face import LinkOracle
from triple_rule_search import assess


def generated_subgroup(group, generators):
    found = {group.identity, *generators}
    while True:
        closure = found | {group.inv[a] for a in found} | {group.mul[a][b] for a in found for b in found}
        if closure == found:
            return tuple(sorted(found))
        found = closure


class AbelianFaceFactor:
    def __init__(self, side, group, table, subgroup):
        if group.n != 8 or group.identity != 0:
            raise ValueError('this finite-state encoding requires the derived order-eight group')
        self.side = side
        assess(group, table)
        self.oracle = LinkOracle(side, group)
        self.group, self.table, self.subgroup = group, table, tuple(subgroup)
        if generated_subgroup(group, subgroup) != self.subgroup:
            raise ValueError('not a closed subgroup')
        if any(group.mul[a][b] != group.mul[b][a] for a, b in itertools.product(subgroup, repeat=2)):
            raise ValueError('face-only factor requires an abelian link subgroup')
        for triple in itertools.product(subgroup, repeat=3):
            target = self.decode(table[self.encode(triple)])
            if any(a not in subgroup for a in target):
                raise ValueError('table does not preserve the link subgroup')
        # Cyclic rotation is harmless only because all subgroup transports commute.
        canonical = lambda face: min(face[i:]+face[:i] for i in range(len(face)))
        face_ids = {canonical(tuple(face[:-1])): i for i, face in enumerate(self.oracle.faces)}
        self.patches = [tuple(face_ids[canonical(tuple(face[:-1]))] for face in spec[::-1])
                        for spec in self.oracle.specs]
        self.incident_patches = [set() for _ in self.oracle.faces]
        for p, faces in enumerate(self.patches):
            for f in faces:
                self.incident_patches[f].add(p)

    @staticmethod
    def encode(triple):
        a, b, c = triple
        return (a*8+b)*8+c

    @staticmethod
    def decode(code):
        return code//64, (code//8) % 8, code % 8

    def project(self, links):
        if any(a not in self.subgroup for a in links):
            raise ValueError('link state is outside the proven subgroup factor')
        return bytes(self.oracle.transport(links, path) for path in self.oracle.face_paths)

    def step(self, state, patch):
        faces = self.patches[patch]
        code = self.encode(state[f] for f in faces)
        target = self.table[code]
        if target == code:
            return state
        result = bytearray(state)
        for f, a in zip(faces, self.decode(target)):
            result[f] = a
        return bytes(result)

    def candidate_patches(self, state):
        # Flat triples are fixed. This sparse frontier omits only proven idle events.
        return sorted(set().union(*(self.incident_patches[f] for f, a in enumerate(state) if a)))

    def signature(self, state):
        weights = {0: 0, 1: 1, 4: 1, 5: 2}
        totals, products = [0, 0], [0, 0]
        for a, part in zip(state, self.oracle.colors):
            totals[part] += weights[a]
            products[part] = self.group.mul[products[part]][a]
        return tuple(totals+products)

    def lift_faces(self, state):
        """Construct a link preimage by eliminating leaves of a dual spanning tree."""
        if len(state) != len(self.oracle.faces) or any(a not in self.subgroup for a in state):
            raise ValueError('invalid face-factor state')
        links = [0]*len(self.oracle.edges)
        incident = [[] for _ in links]
        for f, path in enumerate(self.oracle.face_paths):
            for e, reverse in path:
                incident[e].append((f, reverse))
        adjacent = [[] for _ in state]
        for e, ((f, reverse), (other, _)) in enumerate(incident):
            adjacent[f].append((other, e, reverse))
            adjacent[other].append((f, e, not reverse))
        parents, queue = {0: None}, [0]
        for f in queue:
            for other, e, reverse in adjacent[f]:
                if other not in parents:
                    parents[other] = (f, e, not reverse)
                    queue.append(other)
        for f in reversed(queue[1:]):
            _, e, reverse = parents[f]
            current = self.oracle.transport(links, self.oracle.face_paths[f])
            delta = self.group.mul[self.group.inv[current]][state[f]]
            links[e] = self.group.mul[links[e]][self.group.inv[delta] if reverse else delta]
        if self.project(links) != state:
            raise ValueError('face field has nonidentity global product and cannot be lifted')
        return links


def link_pair(factor, first=1, second=4):
    side = factor.side
    center = (side//2)*side+side//2
    other = (side//2)*side+(side//2+max(1, side//3)) % side
    end = (other//side)*side+(other+1) % side
    links = [factor.group.identity]*len(factor.oracle.edges)
    for u, v, value in ((center, center+1, first), (other, end, second)):
        edge = factor.oracle.edges.index(tuple(sorted((u, v))))
        links[edge] = value if u < v else factor.group.inv[value]
    return links


def search(factor, initial, collision, max_states=100000, exhaustive=False):
    """BFS all changing updates, retaining a shortest pre-collision witness.

    Inert events can be omitted: the boundary-fixed lift is the identity when
    all three target holonomies are unchanged. A truncated BFS never proves
    unreachability. In exhaustive mode, continue through collisions as well.
    """
    queue = deque([initial])
    parents = {initial: None}
    first = None
    active_edges = collision_edges = 0
    min_idle = len(factor.patches)
    invariant = factor.signature(initial)
    while queue:
        state = queue.popleft()
        active_here = 0
        for p in factor.candidate_patches(state):
            faces = factor.patches[p]
            code = factor.encode(state[f] for f in faces)
            target = factor.table[code]
            if code == target:
                continue
            active_edges += 1
            active_here += 1
            if collision[code] != code:
                collision_edges += 1
                if first is None:
                    path, ancestor = [p], state
                    while parents[ancestor] is not None:
                        ancestor, event = parents[ancestor]
                        path.append(event)
                    first = list(reversed(path))
                    if not exhaustive:
                        return {'complete': False, 'stop_reason': 'first collision found',
                                'states_discovered': len(parents), 'shortest_collision_path': first}
            following = factor.step(state, p)
            if factor.step(following, p) != state:
                raise ValueError('face-factor transition is not an involution')
            if following not in parents:
                if len(parents) >= max_states:
                    return {'complete': False, 'stop_reason': 'state budget',
                            'states_discovered': len(parents), 'shortest_collision_path': first}
                parents[following] = (state, p)
                if factor.signature(following) != invariant:
                    raise ValueError('reachable state left its charge/product sector')
                queue.append(following)
        min_idle = min(min_idle, len(factor.patches)-active_here)
    return {'complete': True, 'stop_reason': 'reachable component exhausted',
            'states_discovered': len(parents), 'shortest_collision_path': first,
            'directed_changing_events': active_edges, 'directed_collision_events': collision_edges,
            'minimum_idle_patches_per_state': min_idle,
            'central_face_histogram': dict(sorted(Counter(s.count(5) for s in parents).items()))}


def replay(factor, links, path, collision):
    original = links[:]
    state = factor.project(links)
    records = []
    for p in path:
        before, target = factor.oracle.update(links, p, factor.table)
        state = factor.step(state, p)
        if factor.project(links) != state:
            raise ValueError('face factor diverged from independent raw link replay')
        records.append({'patch': p, 'before': factor.decode(before), 'after': factor.decode(target),
                        'collision': collision[before] != before, 'face_holonomies': list(state)})
    final = links[:]
    for p in reversed(path):
        factor.oracle.update(links, p, factor.table)
    if links != original:
        raise ValueError('raw link encounter path did not invert exactly')
    return {'initial_links': original, 'final_links': final, 'steps': records,
            'independent_factor_replay': True, 'exact_link_inverse': True}


def compiled_replay(factor, initial, path, expected):
    inputs = [factor.side, *factor.table, *initial, len(path), *path]
    process = subprocess.run(['build/wgphysics_three_face_replay'], input=' '.join(map(str, inputs))+'\n',
                             text=True, capture_output=True, check=True, timeout=120)
    result = json.loads(process.stdout)
    if result['final_links'] != expected['final_links'] or not result['exact_link_inverse']:
        raise ValueError('compiled link replay disagrees with independent oracle')
    steps = [[r['patch'], factor.encode(r['before']), factor.encode(r['after'])] for r in expected['steps']]
    if result['steps'] != steps or bytes(result['final_face_holonomies']) != factor.project(expected['final_links']):
        raise ValueError('compiled face sequence disagrees with exact factor')
    return True


def random_schedule(factor, initial, collision, attempts, seed, control=None):
    """Uniform attempted patches, including idle events in occupation times."""
    rng = random.Random(seed)
    state = factor.project(initial)
    control_state = state
    control_path = []
    path, event_times, residence = [], [], Counter()
    for tick in range(attempts):
        p = rng.randrange(len(factor.patches))
        if control is not None:
            next_control = control.step(control_state, p)
            if next_control != control_state:
                control_path.append(p)
            control_state = next_control
        following = factor.step(state, p)
        if following != state:
            path.append(p)
            event_times.append(tick+1)
        state = following
        residence[state.count(5)] += 1
    raw = replay(factor, initial[:], path, collision)
    compiled_replay(factor, initial, path, raw)
    collision_times = [t for t, r in zip(event_times, raw['steps']) if r['collision']]
    result = {'seed': seed, 'attempted_updates': attempts, 'changing_events': len(path),
            'collision_events': len(collision_times), 'first_collision_attempt': next(iter(collision_times), None),
            'central_face_residence_attempts': dict(sorted(residence.items())),
            'changing_patch_path': path, 'changing_event_attempts': event_times,
            'collision_attempts': collision_times, 'final_links': raw['final_links'],
            'independent_raw_link_replay': True, 'compiled_raw_link_replay': True,
            'exact_link_inverse': True}
    if control is not None:
        control_raw = replay(control, initial[:], control_path, collision)
        compiled_replay(control, initial, control_path, control_raw)
        if not collision_times and control_raw['final_links'] != raw['final_links']:
            raise ValueError('no-collision run differs from its matched transport control')
        result['transport_control'] = {
            'same_attempted_patch_schedule': True, 'changing_events': len(control_path),
            'changing_patch_path': control_path, 'final_links': control_raw['final_links'],
            'final_face_class_difference': sum(factor.group.sectors[a] != factor.group.sectors[b]
                                               for a, b in zip(state, control_state)),
            'independent_raw_link_replay': True, 'compiled_raw_link_replay': True,
            'exact_link_inverse': True}
    return result


def verify_local_factor(factor):
    """Exhaust all subgroup assignments on a fan's seven distinct read links."""
    p = 0
    support = sorted(factor.oracle.reads[p])
    if len(support) != 7:
        raise ValueError('expected seven fan links')
    checked = 0
    for values in itertools.product(factor.subgroup, repeat=7):
        links = [factor.group.identity]*len(factor.oracle.edges)
        for e, a in zip(support, values):
            links[e] = a
        target = factor.step(factor.project(links), p)
        factor.oracle.update(links, p, factor.table)
        if factor.project(links) != target:
            raise ValueError('exhaustive local factor check failed')
        checked += 1
    return checked


def sector_count(factor):
    """Algebraic candidate states: each bipartite part has weight 2/product z."""
    g = factor.group
    # Establish the local conservation laws used by the combinatorial count.
    weight = {0: 0, 1: 1, 4: 1, 5: 2}
    for triple in itertools.product(factor.subgroup, repeat=3):
        target = factor.decode(factor.table[factor.encode(triple)])
        if target[1] != triple[1]:
            raise ValueError('middle-face conservation assumption failed')
        if (sum(weight[triple[i]] for i in (0, 2)) != sum(weight[target[i]] for i in (0, 2))
                or g.mul[triple[0]][triple[2]] != g.mul[target[0]][target[2]]):
            raise ValueError('sublattice weight or product is not conserved')
        if all(a in (0, 5) for a in triple) and any(a not in (0, 5) for a in target):
            raise ValueError('central-only sector is not closed')
    m = factor.side**2
    if Counter(factor.oracle.colors) != {0: m, 1: m}:
        raise ValueError('unequal dual sublattices')
    return {'faces_per_part': m, 'all_charge_and_product_compatible_states': m**4,
            'reflection_present_states': m**4-m*m,
            'central_only_states': m*m,
            'reflection_present_central_histogram': {0: m*m*(m-1)**2, 1: 2*m*m*(m-1)},
            'uniform_reflection_present_mean_central_faces': [2, m+1],
            'reflection_present_directed_collision_events': 48*m*(m-1),
            'uniform_collision_probability_per_attempt': [8, m*m*(m+1)],
            'uniform_expected_collisions_per_sweep': [48, m*(m+1)]}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--side', type=int, default=3)
    parser.add_argument('--max-states', type=int, default=100000)
    parser.add_argument('--exhaustive', action='store_true')
    parser.add_argument('--random-trials', type=int, default=0)
    parser.add_argument('--sweeps', type=int, default=100)
    parser.add_argument('--seed', type=int, default=5920941)
    parser.add_argument('--output', type=Path, default=Path('out/three-face-reachability.json'))
    args = parser.parse_args()
    if not 3 <= args.side <= 24 or args.max_states <= 0 or not 0 <= args.random_trials <= 128 or args.sweeps <= 0:
        parser.error('side must be 3..24 and state budget positive')
    census = json.loads(Path('data/d4-triple-rule-search.json').read_text())
    g = derive_fiber_group()
    factor = AbelianFaceFactor(args.side, g, census['combined_rule']['table'], generated_subgroup(g, [1, 4]))
    control = AbelianFaceFactor(args.side, g, census['transport_control']['table'], factor.subgroup)
    initial = link_pair(factor)
    if factor.signature(factor.project(initial)) != (2, 2, 5, 5):
        raise ValueError('two-link seed has the wrong sublattice charges or products')
    result = search(factor, factor.project(initial), census['selected_rule']['table'], args.max_states, args.exhaustive)
    if result['shortest_collision_path'] is not None:
        result['witness'] = replay(factor, initial, result['shortest_collision_path'], census['selected_rule']['table'])
        result['witness']['compiled_raw_link_replay'] = compiled_replay(factor, initial, result['shortest_collision_path'], result['witness'])
    result['sector_count'] = sector_count(factor)
    if result['complete'] and result['states_discovered'] != result['sector_count']['reflection_present_states']:
        raise ValueError('reachable component does not saturate the candidate reflection-present sector')
    if result['complete'] and result['directed_collision_events'] != result['sector_count']['reflection_present_directed_collision_events']:
        raise ValueError('reachable event count disagrees with combinatorial stationary rate')
    if args.exhaustive:
        central_links = [0]*len(factor.oracle.edges)
        central_links[0] = 5
        central = search(factor, factor.project(central_links), census['selected_rule']['table'], args.max_states, True)
        if central['complete'] and central['states_discovered'] != result['sector_count']['central_only_states']:
            raise ValueError('central-only component does not saturate its candidate sector')
        result['central_only_component'] = central
    result['random_schedule_runs'] = [random_schedule(factor, initial, census['selected_rule']['table'],
                                                     args.sweeps*len(factor.patches), args.seed+i*100003, control)
                                       for i in range(args.random_trials)]
    result.update(schema=1, side=args.side, subgroup=factor.subgroup, initial_links=initial,
                  selected_rule_id=census['selected_rule']['id'],
                  random_schedule='Python random.Random(seed).randrange(patch_count); each sweep is patch_count attempts',
                  scope='exact abelian subgroup face factor; arbitrary sequential schedules; no typicality or physical scattering claim')
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2)+'\n')
    print({k: v for k, v in result.items() if k not in ('witness', 'random_schedule_runs', 'initial_links')})
    if result['random_schedule_runs']:
        print('Random-schedule collision counts:', [r['collision_events'] for r in result['random_schedule_runs']])


if __name__ == '__main__':
    main()
