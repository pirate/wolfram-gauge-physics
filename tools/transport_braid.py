#!/usr/bin/env python3
"""Closed transport of actual gauge defects; no quantum braid amplitudes assumed."""
import argparse
import hashlib
import itertools
import json
from collections import Counter, deque
from pathlib import Path

from bank_equilibrium import BankFactor, compiled_run
from d4_loop_observer import CentralExtension, LoopForest, compiled_observe
from face_energy_obstruction import FiniteGroup, derive_fiber_group
from fiber_mode_probe import lifted_graph_laplacian, representation
from mode_localization import ModeGeometry, full_profile


def hurwitz(group, values, index, inverse=False):
    output = list(values)
    a, b = values[index:index+2]
    output[index:index+2] = ((b, group.mul[group.inv[b]][group.mul[a][b]]) if inverse else
                             (group.mul[a][group.mul[b][group.inv[a]]], a))
    return tuple(output)


def pure_word(group, values, i, j, inverse=False):
    """Execute adjacent Hurwitz moves, not the simplified central-bit formula."""
    if not 0 <= i < j < len(values):
        raise ValueError('invalid pure-braid pair')
    preparation = list(range(j-1, i, -1))
    for k in preparation:
        values = hurwitz(group, values, k)
    values = hurwitz(group, hurwitz(group, values, i, inverse=inverse), i, inverse=inverse)
    for k in reversed(preparation):
        values = hurwitz(group, values, k, inverse=True)
    return tuple(values)


def algebra_census(maximum=4):
    if not 0 <= maximum <= 5:
        raise ValueError('bounded pure-braid census supports zero through five entries')
    group = derive_fiber_group()
    extension = CentralExtension(group)
    rows = []
    for size in range(maximum+1):
        pairs = list(itertools.combinations(range(size), 2))
        distribution = Counter()
        for values in itertools.product(range(8), repeat=size):
            for i, j in pairs:
                a, b = values[i], values[j]
                commutator = group.mul[a][group.mul[b][group.mul[group.inv[a]][group.inv[b]]]]
                expected = list(values)
                expected[i] = group.mul[commutator][a]
                expected[j] = group.mul[commutator][b]
                if pure_word(group, values, i, j) != tuple(expected):
                    raise ValueError('actual Hurwitz word differs from central commutator law')
            bits = [extension.bits[x][:2] for x in values]
            noncentral = [q for q in bits if q != (0, 0)]
            rank = min(2, len(set(noncentral)))
            total = tuple(sum(q[i] for q in bits) % 2 for i in range(2))
            exponent = 0 if rank < 2 else len(noncentral)-(3 if total == (0, 0) else 2)
            expected_size = 2**exponent
            normalized = extension.normalized(values)[0]
            seen, queue = {normalized}, [normalized]
            for state in queue:
                for i, j in pairs:
                    target = extension.normalized(pure_word(group, state, i, j))[0]
                    if target not in seen:
                        seen.add(target)
                        queue.append(target)
            if len(seen) != expected_size:
                raise ValueError('pure-braid gauge-orbit count differs from the derived rank formula')
            distribution[len(seen)] += 1
        rows.append({'tuple_length': size, 'raw_tuples_checked': 8**size,
                     'orbit_size_counts_over_raw_input_tuples': dict(sorted(distribution.items()))})
    return {'census': rows, 'scope': 'Classical fixed ordered D4 holonomy tuples modulo simultaneous conjugation; not a count of all torus connections or quantum fusion states.',
            'pure_word_law': '(g_i,g_j) -> (kappa*g_i,kappa*g_j), kappa=[g_i,g_j] in {1,z}',
            'pure_action': 'Pair generators commute and square to identity; noncommuting flux group does not imply noncommuting pure-winding actions.',
            'orbit_formula': 'For quotient rank below two: 1. Otherwise 2^(n-3) if the total quotient is zero, and 2^(n-2) otherwise; n counts noncentral entries.'}


def triangle_candidate():
    """A graph-derived next-fiber criterion, deliberately separate from raw-link results."""
    edges = ((0, 1), (1, 2), (2, 0))
    canonical_edges = {frozenset(e) for e in edges}
    group = FiniteGroup([p for p in itertools.permutations(range(3))
                         if {frozenset((p[u], p[v])) for u, v in edges} == canonical_edges])
    canonical = lambda values: min(tuple(row[x] for x in values) for row in group.conj)
    def product(values):
        result = group.identity
        for x in values:
            result = group.mul[result][x]
        return result
    pairs = list(itertools.combinations(range(4), 2))
    checked, noncommuting, witness = 0, 0, None
    for a, b in itertools.product(range(group.n), repeat=2):
        initial = (a, group.inv[a], b, group.inv[b])
        for first, second in itertools.combinations(pairs, 2):
            ab = pure_word(group, pure_word(group, initial, *first), *second)
            ba = pure_word(group, pure_word(group, initial, *second), *first)
            checked += 1
            if product(ab) != group.identity or product(ba) != group.identity:
                raise ValueError('triangle-fiber pure words changed total flux')
            if canonical(ab) != canonical(ba):
                noncommuting += 1
                if witness is None:
                    witness = {'initial': initial, 'first_pair': first, 'second_pair': second,
                               'first_then_second': ab, 'second_then_first': ba,
                               'canonical_first_then_second': canonical(ab),
                               'canonical_second_then_first': canonical(ba)}
    return {'fiber_edges': edges, 'automorphisms_derived_from_adjacency': group.elements,
            'paired_neutral_seed_and_generator_pair_cases': checked,
            'noncommuting_on_gauge_orbits': noncommuting, 'witness': witness,
            'scope': 'Algebraic probe for a triangle fiber only. The current mixed-bank engine and the raw mesh experiments below still use the square fiber; no triangle-fiber mesh evolution is claimed.'}


def spectral_word_witness(geometry, states, full_bundle=False):
    """Exact characteristic-polynomial comparison and first different trace moment."""
    from flint import fmpz_mat
    def matrix(links):
        if full_bundle:
            return lifted_graph_laplacian(list(range(geometry.size)), geometry.edges, links,
                                          geometry.group, geometry.fiber_edges if hasattr(geometry, 'fiber_edges')
                                          else representation(geometry.group)['fiber_edges'])
        return geometry.operator(links).toarray().astype(int).tolist()
    coefficients = [list(reversed(list(fmpz_mat(matrix(s['links'])).charpoly()))) for s in states]
    first = next((i for i, (a, b) in enumerate(zip(coefficients[0], coefficients[1])) if a != b), None)
    traces = []
    if first is not None:
        for c in coefficients:
            moments = [len(c)-1]
            for k in range(1, first+1):
                moments.append(-k*int(c[k])-sum(int(c[i])*moments[k-i] for i in range(1, k)))
            traces.append(str(moments[first]))
    return {'method': 'Exact FLINT integer characteristic polynomials; Newton identities for traces',
            'same_characteristic_polynomial_as_lap_zero': [c == coefficients[0] for c in coefficients],
            'first_distinct_trace_power': first, 'trace_at_first_distinct_power_by_lap': traces,
            'characteristic_polynomial_sha256_by_lap': [hashlib.sha256(','.join(map(str, c)).encode()).hexdigest() for c in coefficients]}


class FixedMeshTransport:
    """Group-independent geometric paths and explicit frame witnesses on a connected mesh."""
    def compile_dual_topology(self):
        self.adjacency = [set() for _ in self.geometry.faces]
        self.patch = {}
        for p, (u, v) in enumerate(self.factor.pairs):
            self.adjacency[u].add(v)
            self.adjacency[v].add(u)
            self.patch.setdefault((u, v), p)
            self.patch.setdefault((v, u), p)

    def tree_transports(self, links):
        g = self.geometry.group
        spec = self.forest.components[0]
        transports = {spec['root']: g.identity}
        for v in spec['queue'][1:]:
            u, edge, reverse = spec['parents'][v]
            value = g.inv[links[edge]] if reverse else links[edge]
            transports[v] = g.mul[value][transports[u]]
        return transports

    def frame_witness(self, before, after):
        g = self.geometry.group
        a, b = self.tree_transports(before), self.tree_transports(after)
        for root_frame in range(g.n):
            frames = [g.mul[b[v]][g.mul[root_frame][g.inv[a[v]]]] for v in range(self.geometry.size)]
            transformed = [g.mul[frames[v]][g.mul[x][g.inv[frames[u]]]]
                           for (u, v), x in zip(self.geometry.edges, before)]
            if transformed == after:
                return frames
        return None

    def based_faces(self, links):
        g = self.geometry.group
        transports = self.tree_transports(links)
        values = []
        for face, path in zip(self.geometry.faces, self.geometry.paths):
            h = g.identity
            for edge, reverse in path:
                value = g.inv[links[edge]] if reverse else links[edge]
                h = g.mul[value][h]
            t = transports[face[0]]
            values.append(g.mul[g.inv[t]][g.mul[h][t]])
        return values

    def empty_hexagon_detours(self, links, walk):
        occupied = {i for i, c in enumerate(self.geometry.sectors(links)) if c != self.geometry.group.identity}
        candidates = []
        for k, (u, v) in enumerate(zip(walk, walk[1:])):
            shared = set(self.geometry.faces[u][:3]) & set(self.geometry.faces[v][:3])
            for center in sorted(shared):
                star = {i for i, f in enumerate(self.geometry.faces) if center in f[:3]}
                path = [u]
                previous = v
                while path[-1] != v:
                    choices = sorted((self.adjacency[path[-1]] & star)-{previous})
                    if len(choices) != 1:
                        raise ValueError('expected a degree-six primal-vertex star')
                    previous = path[-1]
                    path.append(choices[0])
                    if len(path) > 6:
                        raise ValueError('star detour failed to close')
                if len(path) != 6 or set(path[1:-1]) & (occupied-{walk[0]}):
                    continue
                altered = walk[:k]+path+walk[k+2:]
                if any(self.winding(altered, f) != self.winding(walk, f) for f in occupied-{walk[0]}):
                    raise ValueError('empty star detour changed a defect winding')
                candidates.append({'primal_vertex': center, 'replaced_step': k, 'walk': altered})
        return candidates

    def distances(self, origin):
        distances, queue = {origin: 0}, [origin]
        for u in queue:
            for v in sorted(self.adjacency[u]):
                if v not in distances:
                    distances[v] = distances[u]+1
                    queue.append(v)
        return distances

    def cycle(self, center):
        distances = self.distances(center)
        core = {v for v, d in distances.items() if d in (1, 2, 3)}
        while True:
            leaves = {v for v in core if len(self.adjacency[v] & core) < 2}
            if not leaves:
                break
            core -= leaves
        if len(core) != 12 or any(len(self.adjacency[v] & core) != 2 for v in core):
            raise ValueError('expected the contractible 12-face transport ring')
        walk = [min(core)]
        while len(walk) < len(core):
            choices = sorted((self.adjacency[walk[-1]] & core)-set(walk))
            if not choices:
                raise ValueError('transport ring is disconnected')
            walk.append(choices[0])
        if walk[0] not in self.adjacency[walk[-1]]:
            raise ValueError('transport ring does not close')
        return walk+[walk[0]]

    def point(self, face):
        side = self.geometry.side
        vertices = self.geometry.faces[face][:3]
        xs, ys = [v % side for v in vertices], [v//side for v in vertices]
        if max(xs)-min(xs) > 1 or max(ys)-min(ys) > 1:
            raise ValueError('planar winding audit refuses a face crossing the torus seam')
        # Three times the barycenter: exact integer polygon predicates.
        return sum(xs), sum(ys)

    def winding(self, walk, face):
        p = self.point(face)
        result = 0
        for u, v in zip(walk, walk[1:]):
            a, b = self.point(u), self.point(v)
            cross = (b[0]-a[0])*(p[1]-a[1])-(b[1]-a[1])*(p[0]-a[0])
            if cross == 0 and min(a[0], b[0]) <= p[0] <= max(a[0], b[0]) and min(a[1], b[1]) <= p[1] <= max(a[1], b[1]):
                raise ValueError('winding target lies on the transport polygon')
            if a[1] <= p[1] < b[1] and cross > 0:
                result += 1
            elif b[1] <= p[1] < a[1] and cross < 0:
                result -= 1
        return result

    def step(self, links, u, v, events):
        before = self.geometry.sectors(links)
        if before[u] == self.geometry.group.identity or before[v] != self.geometry.group.identity:
            raise ValueError('transport requires one nonflat source and an empty target')
        patch = self.patch[u, v]
        code, target = self.factor.oracle.update(links, 0, patch, self.factor.tables[0])
        expected = before[:]
        expected[u], expected[v] = expected[v], expected[u]
        if self.geometry.sectors(links) != expected:
            raise ValueError('raw transport failed the independent face-class swap')
        events.append([len(events)+1, 0, patch, code, target])

    def move(self, links, origin, target, events):
        occupied = {i for i, c in enumerate(self.geometry.sectors(links)) if c != self.geometry.group.identity and i != origin}
        parents, queue = {origin: None}, deque([origin])
        while queue and target not in parents:
            u = queue.popleft()
            for v in sorted(self.adjacency[u]-occupied):
                if v not in parents:
                    parents[v] = u
                    queue.append(v)
        if target not in parents:
            raise ValueError('no vacancy path to target')
        path = [target]
        while path[-1] != origin:
            path.append(parents[path[-1]])
        path.reverse()
        for u, v in zip(path, path[1:]):
            self.step(links, u, v, events)
        return path

class TransportExperiment(FixedMeshTransport):
    def __init__(self, side=12):
        if not 10 <= side <= 24:
            raise ValueError('transport audit supports sides 10 through 24')
        self.geometry = ModeGeometry(side)
        self.factor = BankFactor(side, json.loads(Path('data/d4-triple-channels.json').read_text()))
        self.compile_dual_topology()
        self.forest = LoopForest(range(self.geometry.size), self.geometry.edges, self.geometry.group)

    def loop_relations(self, states):
        g = self.geometry.group
        based = [self.based_faces(s['links']) for s in states]
        pairs = []
        occupied = [i for i, c in enumerate(states[0]['face_classes']) if c]
        for k, i in enumerate(occupied):
            for j in occupied[k+1:]:
                if states[0]['face_classes'][i] != states[0]['face_classes'][j]:
                    continue
                values = [g.mul[row[i]][g.inv[row[j]]] for row in based]
                if any(x not in (0, self.forest.extension.z) for x in values):
                    raise ValueError('same-class reflection relation is not central')
                pairs.append({'faces': [i, j], 'word': 'H_i H_j^-1 at the fixed spanning-tree root',
                              'central_holonomy_by_lap': values})
        return pairs

    def prepare(self, commuting=False, empty=False):
        links = self.geometry.seed(reflection_only=commuting)
        initial = links[:]
        sectors = self.geometry.sectors(links)
        center = next(i for i, c in enumerate(sectors) if c == (1 if commuting else 2))
        mover = next(i for i, c in enumerate(sectors) if c == 1 and i != center)
        walk = self.cycle(center)
        distances = self.distances(center)
        # Park spectators in the same planar chart but outside the ring. This
        # preparation uses the identical primitive, not direct connection edits.
        parking = []
        for f in range(len(sectors)):
            try:
                self.point(f)
            except ValueError:
                continue
            if distances[f] >= 7 and not sectors[f]:
                parking.append(f)
        events = []
        spectators = [i for i, c in enumerate(sectors) if c and i not in (center, mover)]
        if empty:
            spectators.append(center)
        for f in spectators:
            target = next(t for t in parking if not self.geometry.sectors(links)[t])
            self.move(links, f, target, events)
        self.move(links, mover, walk[0], events)
        occupied = {i for i, c in enumerate(self.geometry.sectors(links)) if c}
        if occupied & set(walk) != {walk[0]}:
            raise ValueError('prepared ring has an extra occupied face')
        winding = {i: self.winding(walk, i) for i in occupied if i != walk[0]}
        if {i for i, n in winding.items() if n} != (set() if empty else {center}):
            raise ValueError('transport ring encloses the wrong defects')
        return initial, links, events, walk, center, winding

    def run(self, commuting=False, empty=False):
        initial, links, events, walk, center, winding = self.prepare(commuting, empty)
        preparation = len(events)
        states = []
        for lap in range(3):
            states.append({'lap': lap, 'links': links[:], 'face_classes': self.geometry.sectors(links),
                           'gauge_signature': self.forest.signature(links),
                           'spectrum': full_profile(self.geometry, links)})
            if lap < 2:
                for u, v in zip(walk, walk[1:]):
                    self.step(links, u, v, events)
        schedule = [e[2] for e in events]  # rule zero's exact rooted placement ids
        checked = compiled_run(self.factor, initial, schedule, len(schedule))[0]
        if checked['events'] != events or checked['final_links'] != links:
            raise ValueError('closed transport differs from independent compiled raw-link replay')
        compiled_observe(self.forest, [s['links'] for s in states])
        occupied = [i for i, c in enumerate(states[0]['face_classes']) if c]
        based = [[self.based_faces(s['links'])[i] for i in occupied] for s in states]
        if not empty:
            pair = sorted((occupied.index(walk[0]), occupied.index(center)))
            if any(pure_word(self.geometry.group, a, *pair) != tuple(b) for a, b in zip(based, based[1:])):
                raise ValueError('based defect transport differs from the independently evaluated pure Hurwitz word')
        elif any(row != based[0] for row in based):
            raise ValueError('empty-interior transport changed the fixed based face tuple')
        deformations = []
        for spec in self.empty_hexagon_detours(states[0]['links'], walk):
            moved, local_events = states[0]['links'][:], []
            for u, v in zip(spec['walk'], spec['walk'][1:]):
                self.step(moved, u, v, local_events)
            witness = self.frame_witness(states[1]['links'], moved)
            if witness is None:
                raise ValueError('empty star path deformation changed the gauge orbit')
            if any(x for v, x in enumerate(witness) if v != spec['primal_vertex']):
                raise ValueError('empty star deformation requires a frame change outside its primal vertex')
            replay = compiled_run(self.factor, states[0]['links'], [e[2] for e in local_events], len(local_events))[0]
            if replay['final_links'] != moved or replay['events'] != local_events:
                raise ValueError('deformed transport differs from independent compiled replay')
            deformations.append(dict(spec, events=local_events, final_links=moved,
                                     frames_from_reference_lap=witness, exact_link_inverse=replay['exact_link_inverse']))
        return {'commuting_control': commuting, 'empty_interior_control': empty,
                'initial_links': initial, 'preparation_events': preparation, 'events': events,
                'cycle': walk, 'center_face': center, 'interior_windings': winding,
                'states': states, 'same_face_classes': all(s['face_classes'] == states[0]['face_classes'] for s in states),
                'gauge_return_by_lap': [s['gauge_signature'] == states[0]['gauge_signature'] for s in states],
                'raw_return_by_lap': [s['links'] == states[0]['links'] for s in states],
                'frames_from_lap_zero': [self.frame_witness(states[0]['links'], s['links']) for s in states],
                'same_class_loop_relations': self.loop_relations(states), 'empty_star_deformations': deformations,
                'based_defect_faces': occupied, 'based_defect_holonomies_by_lap': based,
                'actual_transport_matches_pure_hurwitz_word': True if not empty else None,
                'exact_spectral_witness': spectral_word_witness(self.geometry, states),
                'exact_full_bundle_spectral_witness': spectral_word_witness(self.geometry, states, full_bundle=True),
                'compiled_and_legacy_gauge_observers_agree': True,
                'exact_link_inverse': checked['exact_link_inverse'],
                'independent_link_replay': checked['independent_link_replay']}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--side', type=int, default=12)
    parser.add_argument('--output', type=Path, default=Path('out/transport-braid.json'))
    args = parser.parse_args()
    experiment = TransportExperiment(args.side)
    result = {'side': args.side, 'scope': 'Classical gauge-link holonomy memory under actual vacancy transport; not quantum statistics or braid amplitudes.',
              'algebra': algebra_census(), 'next_fiber_algebra_probe': triangle_candidate(), 'experiments': []}
    for commuting, empty in ((False, False), (True, False), (False, True)):
        record = experiment.run(commuting, empty)
        result['experiments'].append(record)
        print('Controls', commuting, empty, 'gauge returns', record['gauge_return_by_lap'],
              'band counts', [s['spectrum']['above_flat_band']['rank'] for s in record['states']], flush=True)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, separators=(',', ':'))+'\n')


if __name__ == '__main__':
    main()
