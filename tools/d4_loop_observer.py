#!/usr/bin/env python3
"""Derived central-extension coordinates and complete based-loop observables.

These are classical gauge-invariant data, not quantum amplitudes or a Hamiltonian.
The spanning forest is a coordinate convention on a fixed labeled graph.
"""
import argparse
import itertools
import json
import subprocess
from collections import Counter
from pathlib import Path

from face_energy_obstruction import derive_fiber_group


class CentralExtension:
    def __init__(self, group):
        self.g = group
        e = group.identity
        central = [a for a in range(group.n) if all(group.mul[a][b] == group.mul[b][a] for b in range(group.n))]
        if group.n != 8 or len(central) != 2:
            raise ValueError('expected an order-eight group with a two-element center')
        self.z = next(a for a in central if a != e)
        pairs = [(r, s) for r in range(8) for s in range(8) if r not in central and s not in central
                 and group.mul[r][r] == group.mul[s][s] == e and group.mul[r][s] != group.mul[s][r]]
        if not pairs:
            raise ValueError('no noncommuting involutions for the square-fiber normal form')
        self.r, self.s = min(pairs)
        self.elements = {}
        for a, b, c in itertools.product(range(2), repeat=3):
            self.elements[a, b, c] = group.mul[self.z if c else e][
                group.mul[self.r if a else e][self.s if b else e]]
        self.bits = {value: key for key, value in self.elements.items()}
        if len(self.bits) != 8:
            raise ValueError('normal form does not cover the derived group')
        for x, y in itertools.product(range(8), repeat=2):
            a, b, c = self.bits[x]
            d, f, h = self.bits[y]
            if self.elements[a ^ d, b ^ f, c ^ h ^ (b & d)] != group.mul[x][y]:
                raise ValueError('central-extension multiplication law failed')
        for frame, value in itertools.product(range(8), repeat=2):
            u, v, _ = self.bits[frame]
            a, b, c = self.bits[value]
            if self.bits[group.conj[frame][value]] != (a, b, c ^ (u & b) ^ (v & a)):
                raise ValueError('central-extension frame action failed')

    def normalized(self, values):
        """Fix at most two pivot central bits, in O(number of based loops)."""
        bits = [self.bits[x] for x in values]
        pivots = []
        for i, (a, b, c) in enumerate(bits):
            if (a or b) and (not pivots or (a, b) != bits[pivots[0]][:2]):
                pivots.append(i)
                if len(pivots) == 2:
                    break
        u, v = next((u, v) for u, v in itertools.product(range(2), repeat=2)
                    if all(bits[i][2] ^ (u & bits[i][1]) ^ (v & bits[i][0]) == 0 for i in pivots))
        canonical = tuple(self.elements[a, b, c ^ (u & b) ^ (v & a)] for a, b, c in bits)
        return canonical, tuple(pivots)

    def word(self, values, factors):
        result = self.g.identity
        for i, direction in factors:
            if direction not in (-1, 1):
                raise ValueError('word exponent must be +1 or -1')
            value = values[i] if direction == 1 else self.g.inv[values[i]]
            result = self.g.mul[result][value]
        return result

    def probe(self, values):
        classes = tuple(self.g.sectors[x] for x in values)
        representatives, words = {}, []
        for i, x in enumerate(values):
            kind = self.bits[x][:2]
            if kind == (0, 0):
                continue
            if kind in representatives:
                words.append(((i, 1), (representatives[kind], -1)))
            else:
                representatives[kind] = i
        if len(representatives) == 3:
            i, j, k = representatives.values()
            words.append(((i, 1), (j, 1), (k, -1)))
        central_values = tuple(self.word(values, word) for word in words)
        if any(x not in (self.g.identity, self.z) for x in central_values):
            raise ValueError('relation word is not central')
        normalized, pivots = self.normalized(values)
        rank = len(pivots)
        noncentral = sum(self.bits[x][:2] != (0, 0) for x in values)
        if len(words) != noncentral-rank:
            raise ValueError('relation basis has the wrong number of bits')
        return {'classes': classes, 'relation_words': words,
                'relation_bits': tuple(int(x == self.z) for x in central_values),
                'normalized': normalized, 'quotient_rank': rank,
                'hidden_bits_beyond_individual_loop_classes': noncentral-rank,
                'coordinate_bits': 3*len(values)-rank}

    def reconstruct(self, classes, relation_bits):
        # The relation-word basis depends only on gauge-invariant class labels.
        if any(c not in set(self.g.sectors) for c in classes):
            raise ValueError('invalid conjugacy class')
        raw = [self.elements[self.bits[c][0], self.bits[c][1], 0] if self.bits[c][:2] != (0, 0) else c for c in classes]
        spec = self.probe(raw)
        if len(relation_bits) != len(spec['relation_words']) or any(bit not in (0, 1) for bit in relation_bits):
            raise ValueError('invalid relation bits')
        # The cross-type relation determines the third type's representative;
        # solve it before propagating same-type differences to other loops.
        equations = list(zip(spec['relation_words'], relation_bits))
        for word, bit in sorted(equations, key=lambda row: -len(row[0])):
            target = self.z if bit else self.g.identity
            index = word[2][0] if len(word) == 3 else word[0][0]
            a, b, _ = self.bits[raw[index]]
            candidates = []
            for c in range(2):
                raw[index] = self.elements[a, b, c]
                if self.word(raw, word) == target:
                    candidates.append(raw[index])
            if len(candidates) != 1:
                raise ValueError('relation word did not determine one central bit')
            raw[index] = candidates[0]
        if self.probe(raw)['relation_bits'] != tuple(relation_bits):
            raise ValueError('reconstructed relation data disagree')
        return self.normalized(raw)[0]


class ConnectionForest:
    def __init__(self, vertices, edges, group):
        self.vertices = tuple(sorted(vertices))
        self.edges = tuple(tuple(e) for e in edges)
        if len(set(self.vertices)) != len(self.vertices) or len(set(self.edges)) != len(self.edges):
            raise ValueError('duplicate graph entries')
        vertex_set = set(self.vertices)
        if any(u >= v or u not in vertex_set or v not in vertex_set for u, v in self.edges):
            raise ValueError('expected canonical edges with known distinct vertices')
        self.g = group
        adjacency = {v: [] for v in self.vertices}
        for i, (u, v) in enumerate(self.edges):
            adjacency[u].append((v, i, False))
            adjacency[v].append((u, i, True))
        seen, self.components = set(), []
        for root in self.vertices:
            if root in seen:
                continue
            seen.add(root)
            queue, parents, tree = [root], {}, set()
            for u in queue:
                for v, edge, reverse in sorted(adjacency[u]):
                    if v not in seen:
                        seen.add(v)
                        queue.append(v)
                        parents[v] = (u, edge, reverse)
                        tree.add(edge)
            members = set(queue)
            chords = [i for i, (u, v) in enumerate(self.edges) if u in members and i not in tree]
            self.components.append({'root': root, 'queue': queue, 'parents': parents, 'chords': chords})

    def based_loops(self, links):
        if len(links) != len(self.edges) or any(not isinstance(x, int) or not 0 <= x < self.g.n for x in links):
            raise ValueError('invalid raw link vector')
        result = []
        for component in self.components:
            transports = {component['root']: self.g.identity}
            for v in component['queue'][1:]:
                u, edge, reverse = component['parents'][v]
                value = self.g.inv[links[edge]] if reverse else links[edge]
                transports[v] = self.g.mul[value][transports[u]]
            loops = []
            for edge in component['chords']:
                u, v = self.edges[edge]
                loops.append(self.g.mul[self.g.inv[transports[v]]][self.g.mul[links[edge]][transports[u]]])
            result.append(tuple(loops))
        return result

    def signature(self, links):
        return tuple(min(tuple(row[x] for x in values) for row in self.g.conj)
                     for values in self.based_loops(links))

    def legacy_signature(self, links):
        result = []
        for values in self.based_loops(links):
            flattened = min(tuple(v for x in values for v in self.g.elements[row[x]]) for row in self.g.conj)
            result.extend((len(values), *flattened))
        return result

    def representative(self, signature):
        if len(signature) != len(self.components):
            raise ValueError('signature has wrong number of components')
        links = [self.g.identity]*len(self.edges)
        for component, values in zip(self.components, signature):
            if len(values) != len(component['chords']):
                raise ValueError('signature has wrong cycle rank')
            for edge, x in zip(component['chords'], values):
                if not isinstance(x, int) or not 0 <= x < self.g.n:
                    raise ValueError('invalid loop group value')
                links[edge] = x
        return links

    def fundamental_walk(self, component, chord):
        spec = self.components[component]
        def path(v):
            values = [v]
            while v != spec['root']:
                v = spec['parents'][v][0]
                values.append(v)
            return values[::-1]
        u, v = self.edges[spec['chords'][chord]]
        return path(u)+[v]+path(v)[-2::-1]


class LoopForest(ConnectionForest):
    """Square-fiber specialization retaining its fast central-coordinate normal form."""
    def __init__(self, vertices, edges, group):
        super().__init__(vertices, edges, group)
        self.extension = CentralExtension(group)

    def signature(self, links):
        return tuple(self.extension.normalized(values)[0] for values in self.based_loops(links))


def local_census(extension, maximum_loops=4):
    rows = []
    for size in range(maximum_loops+1):
        observed_to_orbit, orbit_to_observed = {}, {}
        hidden = Counter()
        for values in itertools.product(range(8), repeat=size):
            probe = extension.probe(values)
            observed = (probe['classes'], probe['relation_bits'])
            orbit = min(tuple(row[x] for x in values) for row in extension.g.conj)
            if observed in observed_to_orbit and observed_to_orbit[observed] != orbit:
                raise ValueError('complete observer merges distinct gauge orbits')
            if orbit in orbit_to_observed and orbit_to_observed[orbit] != observed:
                raise ValueError('observer splits a gauge orbit')
            observed_to_orbit[observed] = orbit
            orbit_to_observed[orbit] = observed
            if extension.reconstruct(probe['classes'], probe['relation_bits']) != probe['normalized']:
                raise ValueError('relation-word reconstruction failed')
        for classes, bits in observed_to_orbit:
            hidden[len(bits)] += 1
        rows.append({'based_loops': size, 'raw_tuples': 8**size, 'gauge_orbits': len(observed_to_orbit),
                     'orbit_counts_by_hidden_relation_bits': sorted(hidden.items()),
                     'complete_against_all_simultaneous_conjugations': True,
                     'relation_reconstruction': True})
    return rows


def compiled_observe(forest, states):
    inputs = [len(forest.vertices), len(forest.edges), len(states), *forest.vertices]
    inputs += [x for edge in forest.edges for x in edge]+[x for state in states for x in state]
    process = subprocess.run(['build/wgphysics_loop_observer'], input=' '.join(map(str, inputs))+'\n',
                             text=True, capture_output=True, check=True, timeout=120)
    results = json.loads(process.stdout)['states']
    if len(results) != len(states):
        raise ValueError('compiled observer dropped states')
    for links, result in zip(states, results):
        based = forest.based_loops(links)
        signature = forest.signature(links)
        expected_ranks = [len(forest.extension.normalized(x)[1]) for x in based]
        if (result['based_loops'] != [list(x) for x in based] or
                result['normalized'] != [list(x) for x in signature] or result['ranks'] != expected_ranks or
                result['representative_links'] != forest.representative(signature) or not result['legacy_quotient_agrees']):
            raise ValueError('C++ and Python complete observers disagree')
    return results


def conversion_witness():
    from bank_equilibrium import BankFactor, compiled_run, decode
    from run_channel_banks import seeds
    f = BankFactor(3, json.loads(Path('data/d4-triple-channels.json').read_text()))
    initial = seeds(f.oracle, 3)[4]
    forest = LoopForest(range(9), f.oracle.fan.edges, f.group)
    for patch, spec in enumerate(f.oracle.fan.specs):
        code = f.oracle.code(initial, 1, patch)
        by_classes = {}
        for rule in range(1, 49):
            target = f.tables[rule][code]
            if target == code:
                continue
            probe = forest.extension.probe(decode(target, 3))
            classes = probe['classes']
            if classes in by_classes and by_classes[classes][1]['relation_bits'] != probe['relation_bits']:
                previous, old_probe = by_classes[classes]
                branches = []
                for selected, observed in ((previous, old_probe), (rule, probe)):
                    links = initial[:]
                    f.oracle.update(links, selected, patch, f.tables[selected])
                    run = compiled_run(f, initial, [selected*len(f.fans)+patch], 1)[1]
                    if run['final_links'] != links:
                        raise ValueError('microscopic conversion witness disagrees with C++')
                    if any(a != b for i, (a, b) in enumerate(zip(initial, links)) if i not in f.oracle.fan.writes[patch]):
                        raise ValueError('conversion witness altered an exterior link')
                    branches.append({'rule': selected, 'census_rule_id': f.family['rule_ids'][selected-1],
                                     'schedule': [selected*len(f.fans)+patch], 'final_links': links,
                                     'output_triple': decode(f.tables[selected][code], 3),
                                     'face_classes': list(f.project(links)), 'local_loop_probe': observed,
                                     'compiled_raw_inverse': run['exact_link_inverse'],
                                     'independent_link_replay': run['independent_link_replay']})
                left, right = [b['final_links'] for b in branches]
                if f.project(left) != f.project(right) or f.oracle.fan.gauge_equivalent(left, right):
                    raise ValueError('witness is not a hidden gauge-invariant distinction')
                u = spec[0][0]
                exterior = [u, spec[0][1], spec[0][2], spec[1][2], spec[2][2], u]
                ids = {edge: i for i, edge in enumerate(forest.edges)}
                path = [(ids[tuple(sorted((a, b)))], a > b) for a, b in zip(exterior, exterior[1:])]
                boundary = [f.oracle.fan.transport(x, path) for x in (initial, left, right)]
                if len(set(boundary)) != 1:
                    raise ValueError('witness failed exact boundary preservation')
                observed = compiled_observe(forest, [initial, left, right])
                if observed[1]['normalized'] == observed[2]['normalized']:
                    raise ValueError('complete observer missed the hidden branch distinction')
                responses, averaged = [], []
                for branch in branches:
                    links = branch['final_links'][:]
                    response = compiled_run(f, links, [previous*len(f.fans)+patch], 1)[1]
                    f.oracle.update(links, previous, patch, f.tables[previous])
                    if links != response['final_links']:
                        raise ValueError('specified-rule response disagrees with independent links')
                    responses.append({'final_links': links, 'face_classes': list(f.project(links)),
                                      'compiled_raw_inverse': response['exact_link_inverse']})
                    branch_code = f.oracle.code(branch['final_links'], 1, patch)
                    averaged.append(Counter(tuple(f.group.sectors[x] for x in decode(t[branch_code], 3))
                                            for t in f.tables[1:]))
                if responses[0]['face_classes'] == responses[1]['face_classes'] or averaged[0] != averaged[1]:
                    raise ValueError('specified-rule feedback or uniform-bank averaging claim failed')
                return {'side': 3, 'initial_links': initial, 'fan_patch': patch, 'fan_faces': spec,
                        'input_triple': decode(code, 3), 'exterior_walk': exterior,
                        'exterior_holonomy': boundary[0], 'branches': branches,
                        'same_all_face_classes': True, 'same_exterior_links': True,
                        'not_gauge_equivalent': True, 'compiled_complete_observer': observed,
                        'specified_followup_rule': previous, 'specified_rule_responses': responses,
                        'uniform_next_class_counts': sorted(averaged[0].items()),
                        'uniform_next_class_distribution_identical': True}
            by_classes.setdefault(classes, (rule, probe))
    raise ValueError('no noncentral relative-loop conversion witness found')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, default=Path('out/d4-loop-observer.json'))
    args = parser.parse_args()
    extension = CentralExtension(derive_fiber_group())
    result = {'schema': 1, 'generators': {'r': extension.r, 's': extension.s, 'z': extension.z},
              'coordinates': [extension.bits[x] for x in range(8)],
              'multiplication': '(a,b,c)(d,f,h)=(a xor d,b xor f,c xor h xor (b and d))',
              'frame_action': '(a,b,c) -> (a,b,c xor (u and b) xor (v and a))',
              'census': local_census(extension),
              'conversion_witness': conversion_witness(),
              'scope': 'complete classical gauge observer on a fixed labeled graph; not graph-isomorphism canonicalization, quantum amplitudes, or a physical time/energy law'}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2)+'\n')
    print('Exact loop-orbit census:', [(r['based_loops'], r['gauge_orbits']) for r in result['census']])


if __name__ == '__main__':
    main()
