#!/usr/bin/env python3
"""Group-independent raw-link transport, with automorphisms derived from a supplied fiber."""
import argparse
import itertools
import json
import math
import subprocess
from pathlib import Path

from d4_loop_observer import ConnectionForest
from face_energy_obstruction import FiniteGroup
from run_channel_banks import MixedOracle
from run_shared_edge import mesh_geometry
from transport_braid import FixedMeshTransport, pure_word, spectral_word_witness


def triangle_projective_audit(group):
    """Exact coordinates for neutral four-reflection tuples, not a state update.

    Label the triangle vertices by F_3. Every derived reflection is y -> -y+x.
    Neutrality is a-b+c-d=0. Common affine conjugation removes a and identifies
    (u,v)=(b-a,c-a) with (-u,-v). Nonzero classes are P^1(F_3).
    """
    g = group
    if g.n != 6 or len(g.elements[0]) != 3:
        raise ValueError('projective audit requires the full triangle automorphism group')
    reflections = [x for x in range(g.n) if x != g.identity and g.mul[x][x] == g.identity]
    labels = {x: g.elements[x][0] for x in reflections}
    if len(reflections) != 3 or any(g.elements[x] != tuple((labels[x]-y) % 3 for y in range(3)) for x in reflections):
        raise ValueError('derived involutions are not the three affine reflections')
    canonical = lambda t: min(tuple(row[x] for x in t) for row in g.conj)
    projective = lambda p: min(tuple(p), tuple(-x % 3 for x in p))
    def coordinates(t):
        a, b, c, d = [labels[x] for x in t]
        if (a-b+c-d) % 3:
            raise ValueError('four-reflection tuple is not neutral')
        return ((b-a) % 3, (c-a) % 3)
    # Direct substitution in inverse Hurwitz squares gives these matrices.
    # Check them on EVERY neutral raw tuple, not just the four observed states.
    matrices = [((1, 0), (1, 1)), ((2, 2), (1, 0))]
    apply = lambda m, p: tuple(sum(a*b for a, b in zip(row, p)) % 3 for row in m)
    rows = []
    for t in itertools.product(reflections, repeat=4):
        product = g.identity
        for x in t:
            product = g.mul[product][x]
        if product != g.identity:
            continue
        p = coordinates(t)
        for m, j in zip(matrices, (1, 3)):
            if coordinates(pure_word(g, t, 0, j, inverse=True)) != apply(m, p):
                raise ValueError('projective matrix does not realize the Hurwitz action')
        for frame in g.conj:
            if projective(coordinates(tuple(frame[x] for x in t))) != projective(p):
                raise ValueError('projective point changes under a common gauge frame')
        rows.append({'tuple': t, 'canonical_tuple': canonical(t), 'coordinates': p, 'projective_point': projective(p)})
    correspondence = {(r['canonical_tuple'], r['projective_point']) for r in rows}
    if len(correspondence) != 5 or len({a for a, _ in correspondence}) != 5 or len({b for _, b in correspondence}) != 5:
        raise ValueError('projective coordinates fail to classify neutral tuples')
    product_matrix = lambda a, b: tuple(tuple(sum(a[i][k]*b[k][j] for k in range(2)) % 3 for j in range(2)) for i in range(2))
    closure, queue = {((1, 0), (0, 1))}, [((1, 0), (0, 1))]
    for m in queue:
        for a in matrices:
            target = product_matrix(a, m)
            if target not in closure:
                closure.add(target); queue.append(target)
    sl2 = {((a, b), (c, d)) for a, b, c, d in itertools.product(range(3), repeat=4) if (a*d-b*c) % 3 == 1}
    if closure != sl2:
        raise ValueError('derived matrices do not generate SL(2,F_3)')
    points = sorted({b for _, b in correspondence}-{(0, 0)})
    action = lambda m: tuple(points.index(projective(apply(m, p))) for p in points)
    permutations = {action(m) for m in closure}
    alternating = {p for p in itertools.permutations(range(4)) if sum(p[i] > p[j] for i in range(4) for j in range(i+1, 4)) % 2 == 0}
    kernel = sorted(m for m in closure if action(m) == tuple(range(4)))
    if permutations != alternating or kernel != [((1, 0), (0, 1)), ((2, 0), (0, 2))]:
        raise ValueError('projective action is not SL(2,F_3)/{I,-I} = A_4')
    return {'reflection_elements': reflections, 'affine_labels': [labels[x] for x in reflections],
            'neutral_raw_tuples': rows, 'neutral_gauge_orbits': len(correspondence),
            'projective_points': points, 'inverse_pure_matrices': matrices,
            'matrix_group_order': len(closure), 'projective_kernel': kernel,
            'projective_generator_permutations': [action(m) for m in matrices],
            'projective_action_order': len(permutations), 'equals_A4': True}


class TransportGeometry:
    def __init__(self, side, fiber_vertices, fiber_edges):
        if not 3 <= side <= 24 or not 1 <= fiber_vertices <= 6:
            raise ValueError('invalid bounded generic transport geometry')
        canonical = [tuple(sorted(e)) for e in fiber_edges]
        if len(set(canonical)) != len(canonical) or any(len(e) != 2 or not 0 <= e[0] < e[1] < fiber_vertices for e in canonical):
            raise ValueError('invalid simple fiber edges')
        self.fiber_vertices, self.fiber_edges = fiber_vertices, sorted(canonical)
        edge_set = set(canonical)
        self.group = FiniteGroup([p for p in itertools.permutations(range(fiber_vertices))
                                 if {tuple(sorted((p[u], p[v]))) for u, v in canonical} == edge_set])
        if self.group.n > 24:
            raise ValueError('bounded transport supports at most 24 derived automorphisms')
        self.side, self.size = side, side*side
        self.edges, self.faces, self.pair_specs = mesh_geometry(side)
        self.ids = {edge: i for i, edge in enumerate(self.edges)}
        self.paths = [[(self.ids[tuple(sorted((u, v)))], u > v) for u, v in zip(f, f[1:])] for f in self.faces]

    def holonomies(self, links):
        if len(links) != len(self.edges) or any(not isinstance(x, int) or not 0 <= x < self.group.n for x in links):
            raise ValueError('invalid generic link vector')
        result = []
        g = self.group
        for path in self.paths:
            h = g.identity
            for edge, reverse in path:
                value = g.inv[links[edge]] if reverse else links[edge]
                h = g.mul[value][h]
            result.append(h)
        return result

    def sectors(self, links):
        return [self.group.sectors[x] for x in self.holonomies(links)]

    def paired_seed(self, a, b):
        if any(not isinstance(x, int) or not 0 <= x < self.group.n for x in (a, b)):
            raise ValueError('invalid derived seed elements')
        y = x = self.side//2
        vertices = [y*self.side+(x+i) % self.side for i in range(3)]
        links = [self.group.identity]*len(self.edges)
        touched = []
        for (u, v), value in zip(zip(vertices, vertices[1:]), (a, b)):
            edge = tuple(sorted((u, v)))
            links[self.ids[edge]] = value if u < v else self.group.inv[value]
            touched.append(sorted(i for i, f in enumerate(self.faces) if edge in [tuple(sorted(t)) for t in zip(f, f[1:])]))
        return links, touched[0]+touched[1]


class VacancyFactor:
    def __init__(self, geometry):
        g = geometry.group
        self.oracle = MixedOracle(geometry.side, g)
        self.pairs = [tuple(spec[2:]) for spec in geometry.pair_specs]
        self.tables = [[b*g.n+a if (a == g.identity) != (b == g.identity) else a*g.n+b
                        for a, b in itertools.product(range(g.n), repeat=2)]]
        g.validate(self.tables[0], strict=True)


class FiberTransport(FixedMeshTransport):
    def __init__(self, side, fiber_vertices, fiber_edges):
        self.geometry = TransportGeometry(side, fiber_vertices, fiber_edges)
        self.factor = VacancyFactor(self.geometry)
        self.compile_dual_topology()
        self.forest = ConnectionForest(range(self.geometry.size), self.geometry.edges, self.geometry.group)

    def compiled(self, initial, schedule):
        g, geometry = self.geometry.group, self.geometry
        geometry.holonomies(initial)
        if len(schedule) > 100000 or any(not isinstance(p, int) or not 0 <= p < len(self.factor.pairs) for p in schedule):
            raise ValueError('invalid bounded generic schedule')
        inputs = [geometry.side, geometry.fiber_vertices, len(geometry.fiber_edges), len(schedule)]
        inputs += [v for edge in geometry.fiber_edges for v in edge]+list(initial)+list(schedule)
        process = subprocess.run(['build/wgphysics_fiber_transport'], input=' '.join(map(str, inputs))+'\n',
                                 text=True, capture_output=True, check=True, timeout=120)
        record = json.loads(process.stdout)
        if record['automorphisms'] != [list(p) for p in g.elements] or record['initial_signature'] != self.forest.legacy_signature(initial):
            raise ValueError('independently derived group or initial gauge quotient disagrees')
        links, expected = initial[:], []
        initial_faces = geometry.holonomies(links)
        for tick, patch in enumerate(schedule, 1):
            before, before_links = geometry.holonomies(links), links[:]
            code, target = self.factor.oracle.update(links, 0, patch, self.factor.tables[0])
            after = geometry.holonomies(links)
            affected = set(self.factor.pairs[patch])
            shared = self.factor.oracle.pairs[patch][0][0][0]
            if any(a != b for i, (a, b) in enumerate(zip(before, after)) if i not in affected):
                raise ValueError('generic vacancy transport changed a spectator holonomy')
            if any(a != b for i, (a, b) in enumerate(zip(before_links, links)) if i != shared):
                raise ValueError('generic vacancy transport changed a boundary link')
            expected_classes = [g.sectors[x] for x in before]
            u, v = self.factor.pairs[patch]
            if (expected_classes[u] == g.identity) != (expected_classes[v] == g.identity):
                expected_classes[u], expected_classes[v] = expected_classes[v], expected_classes[u]
            if [g.sectors[x] for x in after] != expected_classes:
                raise ValueError('generic transport class swap failed')
            if code != target:
                expected.append([tick, patch, code, target])
        if (record['events'] != expected or record['final_links'] != links or record['initial_faces'] != initial_faces or
                record['final_faces'] != geometry.holonomies(links) or record['final_signature'] != self.forest.legacy_signature(links)):
            raise ValueError('compiled generic transport differs from independent raw-link replay or gauge quotient')
        for patch in reversed(schedule):
            self.factor.oracle.update(links, 0, patch, self.factor.tables[0])
        if links != initial or not record['exact_link_inverse']:
            raise ValueError('generic independent inverse failed')
        record.update(independent_link_replay=True, spectator_holonomies_unchanged=True, boundary_links_unchanged=True)
        return record

    def prepare_four(self, generators):
        links, origins = self.geometry.paired_seed(*generators)
        initial = links[:]
        targets = [2*(y*self.geometry.side+x) for x, y in ((2, 2), (7, 2), (2, 7), (8, 8))]
        if self.geometry.side < 12:
            raise ValueError('four-defect planar preparation requires side at least twelve')
        events, paths = [], []
        for origin, target in zip(origins, targets):
            paths.append(self.move(links, origin, target, events))
        checked = self.compiled(initial, [event[2] for event in events])
        if checked['final_links'] != links:
            raise ValueError('generic preparation replay failed')
        return {'initial_links': initial, 'links': links, 'origins': origins, 'positions': targets,
                'paths': paths, 'events': events, 'checked': checked}

    def circuit(self, links, mover, center):
        # The lead path is retraced after the circuit. All scheduling choices
        # use occupancy and geometry only, never a desired holonomy response.
        moved, events = links[:], []
        cycle = self.cycle(center)
        lead = self.move(moved, mover, cycle[0], events)
        for u, v in zip(cycle, cycle[1:]):
            self.step(moved, u, v, events)
        reverse = list(reversed(lead))
        for u, v in zip(reverse, reverse[1:]):
            self.step(moved, u, v, events)
        walk = lead+cycle[1:]+reverse[1:]
        occupied = {i for i, c in enumerate(self.geometry.sectors(links)) if c and i != mover}
        winding = {i: self.winding(walk, i) for i in sorted(occupied)}
        if {i for i, n in winding.items() if n} != {center}:
            raise ValueError('generic circuit winds around the wrong stationary defects')
        if self.geometry.sectors(moved) != self.geometry.sectors(links):
            raise ValueError('generic circuit did not restore every face class')
        return {'walk': walk, 'windings': winding, 'lead_steps': len(lead)-1,
                'cycle_steps': len(cycle)-1, 'schedule': [event[2] for event in events]}

    def compare_orders(self, generators):
        preparation = self.prepare_four(generators)
        initial = preparation['links']
        positions = preparation['positions']
        circuits = [self.circuit(initial, positions[0], center) for center in positions[1:3]]
        orders = []
        for order in ((0, 1), (1, 0)):
            schedule = [patch for i in order for patch in circuits[i]['schedule']]
            checked = self.compiled(initial, schedule)
            final = checked['final_links']
            if self.geometry.sectors(final) != self.geometry.sectors(initial):
                raise ValueError('ordered circuits did not restore every face class')
            orders.append({'order': order, 'schedule': schedule, 'checked': checked,
                           'based_defect_holonomies': [self.based_faces(final)[f] for f in positions]})
        frames = self.frame_witness(orders[0]['checked']['final_links'], orders[1]['checked']['final_links'])
        equivalent = self.factor.oracle.fan.gauge_equivalent(orders[0]['checked']['final_links'], orders[1]['checked']['final_links'])
        if equivalent != (frames is not None) or equivalent != (orders[0]['checked']['final_signature'] == orders[1]['checked']['final_signature']):
            raise ValueError('three independent full-connection gauge-equivalence checks disagree')
        pair_words = []
        g = self.geometry.group
        for i, j in itertools.combinations(range(4), 2):
            values = [g.mul[row['based_defect_holonomies'][i]][g.inv[row['based_defect_holonomies'][j]]] for row in orders]
            pair_words.append({'positions': [i, j], 'faces': [positions[i], positions[j]],
                               'root_based_word': 'H_i H_j^-1', 'holonomies_by_order': values,
                               'classes_by_order': [g.sectors[x] for x in values]})
        spectra = spectral_word_witness(self.geometry, [{'links': row['checked']['final_links']} for row in orders], full_bundle=True)
        spectra = {k.replace('as_lap_zero', 'as_first_order').replace('by_lap', 'by_order'): v for k, v in spectra.items()}
        orbit = self.orbit_census(initial, circuits)
        basis = self.loop_action_audit(positions, orbit)
        return {'fiber_vertices': self.geometry.fiber_vertices, 'fiber_edges': self.geometry.fiber_edges,
                'group_automorphisms': self.geometry.group.elements, 'generators': generators,
                'preparation': preparation, 'circuits': circuits, 'orders': orders,
                'orders_gauge_equivalent': equivalent, 'frames_between_orders': frames,
                'relative_loop_words': pair_words, 'exact_full_bundle_spectral_comparison': spectra,
                'controlled_circuit_orbit': orbit,
                'loop_action_audit': basis,
                'same_face_classes': True, 'initial_based_defect_holonomies': [self.based_faces(initial)[f] for f in positions]}

    def loop_action_audit(self, positions, orbit):
        """Check a fixed perimeter-ordered observer against the actual raw action.

        The preparation parks faces lower-left, lower-right, upper-left, upper-right.
        This is an explicitly chosen based-loop convention, not a theorem that any
        spatial sorting produces a bouquet basis for arbitrary connections.
        """
        g = self.geometry.group
        ordering = (0, 1, 3, 2)
        faces = [positions[i] for i in ordering]
        tuples = [tuple(self.based_faces(links)[f] for f in faces) for links in orbit['raw_representatives']]
        canonical = lambda t: min(tuple(row[x] for x in t) for row in g.conj)
        if len({canonical(t) for t in tuples}) != len(tuples):
            raise ValueError('four-loop observer is not injective on the controlled component')
        checks = []
        for source, values in enumerate(tuples):
            product = g.identity
            for x in values:
                product = g.mul[product][x]
            if product != g.identity:
                raise ValueError('perimeter-ordered tuple is not neutral')
            for action, j in enumerate((1, 3)):
                witness = orbit['edge_witnesses'][source][action]
                actual = tuple(self.based_faces(witness['checked']['final_links'])[f] for f in faces)
                predicted = pure_word(g, values, 0, j, inverse=True)
                if actual != predicted:
                    raise ValueError('raw circuit does not equal its predicted inverse pure Hurwitz word')
                target = orbit['action_targets_by_state'][source][action]
                if canonical(actual) != canonical(tuples[target]):
                    raise ValueError('based tuple action disagrees with the full-connection action')
                checks.append({'source_state': source, 'action': action, 'target_state': target,
                               'predicted_tuple': predicted, 'actual_tuple': actual})
        result = {'position_order': ordering, 'faces': faces, 'pure_pairs': [(0, 1), (0, 3)],
                  'inverse_squares': True, 'based_tuples': tuples, 'checks': checks,
                  'neutral': True, 'injective_on_controlled_component': True,
                  'scope': 'Fixed spanning-tree based loops in perimeter order; exact equality on every controlled action edge, not a general torus reconstruction from four faces.'}
        if self.geometry.fiber_vertices == 3 and g.n == 6:
            algebra = triangle_projective_audit(g)
            point_by_tuple = {tuple(r['tuple']): r['projective_point'] for r in algebra['neutral_raw_tuples']}
            result['projective_points_by_state'] = [point_by_tuple[t] for t in tuples]
            result['triangle_projective_algebra'] = algebra
        return result

    def orbit_census(self, initial, circuits, maximum=64):
        """Exhaust the complete gauge orbits reachable by these two fixed programs.

        Representatives are actual reachable raw states. Canonical signatures
        deduplicate them; no canonicalized links are injected into the dynamics.
        """
        states = [initial[:]]
        lookup = {self.forest.signature(initial): 0}
        actions, witnesses = [], []
        for links in states:
            targets, records = [], []
            for circuit in circuits:
                checked = self.compiled(links, circuit['schedule'])
                final = checked['final_links']
                signature = self.forest.signature(final)
                if signature not in lookup:
                    if len(states) >= maximum:
                        raise ValueError('controlled orbit budget exhausted; closure not proved')
                    lookup[signature] = len(states)
                    states.append(final)
                target = lookup[signature]
                frames = self.frame_witness(final, states[target])
                if frames is None:
                    raise ValueError('orbit edge lacks an explicit frame-equivalence witness')
                targets.append(target)
                records.append({'checked': checked, 'frames_to_target_representative': frames})
            actions.append(targets)
            witnesses.append(records)
        count = len(states)
        permutations = [tuple(row[j] for row in actions) for j in range(len(circuits))]
        identity = tuple(range(count))
        if any(sorted(p) != list(identity) for p in permutations):
            raise ValueError('a closed reversible circuit failed to induce a permutation')
        compose = lambda a, b: tuple(a[b[i]] for i in range(count))
        closure, queue = {identity}, [identity]
        for p in queue:
            for generator in permutations:
                target = compose(generator, p)
                if target not in closure:
                    if len(closure) >= 10000:
                        raise ValueError('controlled permutation group budget exhausted')
                    closure.add(target)
                    queue.append(target)
        def order(p):
            visited, lengths = set(), []
            for start in identity:
                if start in visited:
                    continue
                v, length = start, 0
                while v not in visited:
                    visited.add(v); length += 1; v = p[v]
                lengths.append(length)
            return math.lcm(*lengths)
        inverse = lambda p: tuple(p.index(i) for i in identity)
        a, b = permutations
        commutator = compose(a, compose(b, compose(inverse(a), inverse(b))))
        alternating_four = count == 4 and closure == {
            p for p in itertools.permutations(range(4))
            if sum(p[i] > p[j] for i in range(4) for j in range(i+1, 4)) % 2 == 0}
        chosen = []
        for circuit in circuits:
            candidates = [d for d in self.empty_hexagon_detours(initial, circuit['walk'])
                          if circuit['lead_steps'] <= d['replaced_step'] < circuit['lead_steps']+circuit['cycle_steps']]
            if not candidates:
                raise ValueError('no empty-star deformation of the winding portion')
            chosen.append(dict(candidates[0], admissible_winding_detours=len(candidates)))
        deformed = []
        for state, links in enumerate(states):
            for action, spec in enumerate(chosen):
                moved, events = links[:], []
                for u, v in zip(spec['walk'], spec['walk'][1:]):
                    self.step(moved, u, v, events)
                checked = self.compiled(links, [event[2] for event in events])
                if checked['final_links'] != moved:
                    raise ValueError('deformed generic circuit replay disagrees')
                reference = witnesses[state][action]['checked']['final_links']
                frames = self.frame_witness(reference, moved)
                if frames is None or any(x for v, x in enumerate(frames) if v != spec['primal_vertex']):
                    raise ValueError('deformed generic circuit is not a single-vertex frame change')
                if self.forest.signature(moved) != self.forest.signature(states[actions[state][action]]):
                    raise ValueError('deformed generic circuit changes its orbit action')
                deformed.append({'source_state': state, 'action': action, 'target_state': actions[state][action],
                                 'checked': checked, 'frames_from_reference_action': frames})
        return {'exhausted': True, 'raw_representatives': states, 'action_targets_by_state': actions,
                'edge_witnesses': witnesses, 'generator_permutations': permutations,
                'generator_orders': [order(p) for p in permutations],
                'permutation_group_order': len(closure), 'permutation_group': sorted(closure),
                'commutator_permutation': commutator, 'commutator_order': order(commutator),
                'equals_alternating_group_on_four_orbits': alternating_four,
                'chosen_winding_detours': chosen, 'deformed_action_checks': deformed,
                'scope': 'Complete gauge-orbit component under these two controlled transport programs, not all microscopic updates or a quantum state space.'}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--side', type=int, default=12)
    parser.add_argument('--output', type=Path, default=Path('out/fiber-transport.json'))
    args = parser.parse_args()
    result = {'side': args.side, 'scope': 'Controlled classical raw-link transport with groups derived from fiber graphs. No quantum amplitudes, force law, or stable matter are asserted.', 'fibers': []}
    for size, commuting in ((3, False), (4, False), (3, True)):
        e = FiberTransport(args.side, size, [(i, (i+1) % size) for i in range(size)])
        g = e.geometry.group
        generators = next((a, b) for a, b in itertools.product(range(g.n), repeat=2)
                          if a != g.identity and b != g.identity and g.mul[a][a] == g.identity and
                          g.mul[b][b] == g.identity and g.mul[a][b] != g.mul[b][a])
        if commuting:
            generators = (generators[0], generators[0])
        record = e.compare_orders(generators)
        record['seed_condition'] = 'commuting' if commuting else 'noncommuting'
        result['fibers'].append(record)
        print('Fiber', size, 'derived group', g.n, 'orders gauge equivalent', record['orders_gauge_equivalent'],
              'based tuples', [r['based_defect_holonomies'] for r in record['orders']], flush=True)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, separators=(',', ':'))+'\n')


if __name__ == '__main__':
    main()
