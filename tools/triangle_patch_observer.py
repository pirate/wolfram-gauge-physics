#!/usr/bin/env python3
"""Exact short-word gauge observers and overlapping-patch quotient kernels."""
import argparse
import itertools
import json
from collections import Counter, defaultdict
from functools import lru_cache
from pathlib import Path

from d4_loop_observer import ConnectionForest
from triangle_charge_current import CurrentProbe


class WordObserver:
    def __init__(self, group, charges):
        self.g, self.charges = group, charges
        self.conjugates = group.conj
        self.reflections = [x for x, q in enumerate(charges) if q == 1]
        self.rotations = [x for x, q in enumerate(charges) if q == 2]
        if group.n != 6 or len(self.reflections) != 3 or len(self.rotations) != 2:
            raise ValueError('short-word coordinates require the derived triangle group')

    def canonical(self, values):
        return min(tuple(row[x] for x in values) for row in self.conjugates)

    @lru_cache(maxsize=None)
    def representatives(self, rank):
        return tuple(sorted({self.canonical(values) for values in itertools.product(range(self.g.n), repeat=rank)}))

    def signature(self, values):
        """Charges of all increasing-index ordered subwords of lengths 1..3."""
        result = []
        for size in range(1, min(3, len(values))+1):
            for indices in itertools.combinations(range(len(values)), size):
                word = self.g.identity
                for i in indices:
                    word = self.g.mul[word][values[i]]
                result.append(self.charges[word])
        return tuple(result)

    def anchored(self, values):
        """Linear-size invariant coordinates, using only words of length <=3."""
        classes = tuple(self.charges[v] for v in values)
        rotations = [v for v in values if self.charges[v] == 2]
        reflections = [v for v in values if self.charges[v] == 1]
        labels, seen = [], {}
        z = rotations[0] if rotations else None
        r = reflections[0] if reflections else None
        for v, q in zip(values, classes):
            if q == 0:
                labels.append(0)
            elif q == 2:
                labels.append(int(self.charges[self.g.mul[z][v]] == 0))
            elif z is None:
                # Equality partition of reflections: common conjugation acts as
                # all permutations of the three labels. At most two comparisons
                # with earlier representatives recover each partition label.
                labels.append(seen.setdefault(v, len(seen)))
            else:
                pair = self.g.mul[r][v]
                labels.append(0 if self.charges[pair] == 0 else
                              1 if self.charges[self.g.mul[pair][z]] == 0 else 2)
        return classes, tuple(labels)

    def reconstruct(self, signature):
        classes, labels = signature
        if len(classes) != len(labels) or any(q not in (0, 1, 2) for q in classes):
            raise ValueError('invalid anchored coordinate dimensions or charges')
        z, r = min(self.rotations), min(self.reflections)
        has_rotation = 2 in classes
        reflected = [r]
        for label in (1, 2):
            reflected.append(next(v for v in self.reflections if v != r and
                                  (self.charges[self.g.mul[self.g.mul[r][v]][z]] == 0) == (label == 1)))
        result = []
        for q, label in zip(classes, labels):
            if q == 0 and label == 0:
                result.append(self.g.identity)
            elif q == 2 and label in (0, 1):
                result.append(z if label == 0 else self.g.inv[z])
            elif q == 1 and label in (0, 1, 2):
                result.append(reflected[label] if has_rotation else self.reflections[label])
            else:
                raise ValueError('invalid anchored coordinate label')
        result = tuple(result)
        if self.anchored(result) != (tuple(classes), tuple(labels)):
            raise ValueError('coordinate anchors are inconsistent with first occurrences')
        return result

    def census(self, rank):
        reps = self.representatives(rank)
        signatures = [self.signature(v) for v in reps]
        if len(set(signatures)) != len(reps):
            raise ValueError('short-word signature misses a simultaneous-conjugacy orbit')
        anchored = [self.anchored(v) for v in reps]
        if len(set(anchored)) != len(reps) or any(self.canonical(self.reconstruct(s)) != v for s, v in zip(anchored, reps)):
            raise ValueError('anchored short words do not reconstruct the exact orbit')
        # Burnside uses centralizers derived from the actual multiplication table.
        centralizers = [sum(self.g.mul[h][x] == self.g.mul[x][h] for x in range(self.g.n))
                        for h in range(self.g.n)]
        burnside = sum(c**rank for c in centralizers)//self.g.n
        if len(reps) != burnside:
            raise ValueError('direct orbit enumeration disagrees with Burnside count')
        return {'rank': rank, 'raw_tuples': self.g.n**rank, 'gauge_orbits': len(reps),
                'word_count': len(signatures[0]), 'distinct_signatures': len(set(signatures)),
                'centralizer_sizes': centralizers, 'anchored_signatures': len(set(anchored))}

    def triple_necessity(self):
        result = []
        reps = self.representatives(3)
        for omitted in range(3, 7):
            seen = {}
            for values in reps:
                signature = self.signature(values)
                reduced = signature[:omitted]+signature[omitted+1:]
                if reduced in seen:
                    other = seen[reduced]
                    result.append({'omitted_word_index': omitted, 'tuples': [list(other), list(values)],
                                   'full_signatures': [list(self.signature(other)), list(signature)]})
                    break
                seen[reduced] = values
            else:
                raise ValueError('claimed necessary product word is redundant')
        return result


class PatchUnion:
    """Tree gauge fixing on the union's actual read graph; no remote frame guess."""
    def __init__(self, experiment, patches):
        patches = tuple(patches)
        self.e, self.patches = experiment, patches
        self.g = experiment.geometry.group
        self.words = WordObserver(self.g, experiment.charges)
        self.edges = sorted(set().union(*(experiment.factor.oracle.fan.reads[p] for p in patches)))
        self.vertices = sorted({v for edge in self.edges for v in experiment.geometry.edges[edge]})
        self.forest = ConnectionForest(self.vertices, [experiment.geometry.edges[i] for i in self.edges], self.g)
        if len(self.forest.components) != 1:
            raise ValueError('overlap read graph is disconnected')
        component = self.forest.components[0]
        self.root = component['root']
        self.tree = [(u, v, self.edges[edge], reverse) for v in component['queue'][1:]
                     for u, edge, reverse in (component['parents'][v],)]
        self.chords = [self.edges[edge] for edge in component['chords']]
        self.rank = len(self.chords)
        self.reps = self.words.representatives(self.rank)
        self.ids = {state: i for i, state in enumerate(self.reps)}
        self.triple_ids = {state: i for i, state in enumerate(self.words.representatives(3))}
        faces = [experiment.fan_faces[p][::-1] for p in patches]
        shared = set(faces[0]) & set(faces[1]) if len(faces) == 2 else set()
        self.shared_positions = tuple(fs.index(next(iter(shared))) for fs in faces) if len(shared) == 1 else None
        self.patch_roots = [experiment.factor.oracle.fan.specs[p][0][0] for p in patches]
        self.connector = []
        if self.shared_positions is not None and self.patch_roots[0] != self.patch_roots[1]:
            a, b = self.patch_roots
            edge = experiment.geometry.ids[tuple(sorted((a, b)))]
            intersection = set.intersection(*(experiment.factor.oracle.fan.reads[p] for p in patches))
            if edge not in intersection:
                raise ValueError('shared-face connector is not in both actual read supports')
            self.connector = [(edge, a > b)]

    def based_loops(self, links):
        return self.forest.based_loops([links[edge] for edge in self.edges])[0]

    def observe(self, links):
        return self.ids[self.words.canonical(self.based_loops(links))]

    def lift(self, values):
        if len(values) != self.rank:
            raise ValueError('wrong number of independent union loops')
        links = [self.g.identity]*len(self.e.geometry.edges)
        for edge, value in zip(self.chords, values):
            links[edge] = value
        return links

    def separate(self, links):
        return tuple(self.triple_ids[self.words.canonical(v)] for v in self.triples(links))

    def triples(self, links):
        n = self.g.n
        codes = [self.e.factor.oracle.code(links, 1, p) for p in self.patches]
        values = [(code//n**2, code//n % n, code % n) for code in codes]
        if self.connector:
            k = self.e.factor.oracle.fan.transport(links, self.connector)
            values[1] = tuple(self.g.mul[self.g.inv[k]][self.g.mul[v][k]] for v in values[1])
        if self.shared_positions is not None:
            i, j = self.shared_positions
            if values[0][i] != values[1][j]:
                raise ValueError('connector transport fails to identify the actual shared face')
        return tuple(values)

    @lru_cache(maxsize=None)
    def gluing_spec(self, key):
        """Relative alignment is C(A) \\ C(h) / C(B), not a free frame choice."""
        if self.shared_positions is None:
            raise ValueError('this gluing coordinate requires exactly one shared face')
        a, b = [self.words.representatives(3)[i] for i in key]
        i, j = self.shared_positions
        matches = [row for row in self.g.conj if row[b[j]] == a[i]]
        if not matches:
            raise ValueError('patch orbits have incompatible shared-face holonomy')
        aligned = tuple(matches[0][v] for v in b)
        common = [h for h, row in enumerate(self.g.conj) if row[a[i]] == a[i]]
        left = [h for h, row in enumerate(self.g.conj) if tuple(row[v] for v in a) == a]
        right = [h for h, row in enumerate(self.g.conj) if tuple(row[v] for v in aligned) == aligned]
        cosets = sorted({tuple(sorted({self.g.mul[x][self.g.mul[h][y]] for x in left for y in right})) for h in common})
        if sorted(x for coset in cosets for x in coset) != common:
            raise ValueError('double cosets do not partition the common centralizer')
        return {'aligned_second_tuple': aligned, 'common_centralizer': common,
                'first_stabilizer': left, 'second_stabilizer': right, 'double_cosets': cosets}

    def sewn(self, links):
        first, second = self.triples(links)
        canonical, frame = min((tuple(row[v] for v in first), h) for h, row in enumerate(self.g.conj))
        key = self.triple_ids[canonical], self.triple_ids[self.words.canonical(second)]
        spec = self.gluing_spec(key)
        actual_second = tuple(self.g.conj[frame][v] for v in second)
        matches = [h for h in spec['common_centralizer']
                   if tuple(self.g.conj[h][v] for v in spec['aligned_second_tuple']) == actual_second]
        if not matches:
            raise ValueError('no compatible relative alignment of actual patch frames')
        cells = {i for i, coset in enumerate(spec['double_cosets']) if any(h in coset for h in matches)}
        if len(cells) != 1:
            raise ValueError('relative alignment does not determine one double coset')
        return (*key, next(iter(cells)))

    def kernel(self):
        """All existing operators whose complete read support lies in the union."""
        edges, oracle = set(self.edges), self.e.factor.oracle
        pairs = [p for p, paths in enumerate(oracle.pairs)
                 if {edge for path in paths for edge, _ in path} <= edges]
        fans = [p for p, reads in enumerate(oracle.fan.reads) if reads <= edges]
        channels = [(patch, 0) for patch in pairs]
        channels += [(patch, rule) for patch in fans for rule in range(1, len(self.e.tables))]
        tables = [[] for _ in channels]
        separate = []
        for values in self.reps:
            links = self.lift(values)
            separate.append(self.separate(links))
            for table, (patch, rule) in zip(tables, channels):
                moved = links[:]
                self.e.factor.oracle.update(moved, rule, patch, self.e.tables[rule])
                table.append(self.observe(moved))
        for table in tables:
            if sorted(table) != list(range(len(self.reps))) or any(table[table[i]] != i for i in range(len(self.reps))):
                raise ValueError('union gauge kernel is not a permutation involution')
        return channels, tables, separate


def audit():
    bank = json.loads(Path('data/triangle-feedback.json').read_text())['bank']
    e = CurrentProbe(6, bank).e
    union = PatchUnion(e, (0, 2))
    if union.rank != 5:
        raise ValueError('selected overlapping patches no longer have five independent loops')
    words = union.words
    channels, tables, separate = union.kernel()
    actor_tables = [table for (patch, rule), table in zip(channels, tables) if patch == union.patches[0] and rule]
    actor_channels = [channel for channel in channels if channel[0] == union.patches[0] and channel[1]]
    fibers = defaultdict(list)
    for state, key in enumerate(separate):
        fibers[key].append(state)
    sewn = [union.sewn(union.lift(v)) for v in union.reps]
    if len(set(sewn)) != len(union.reps):
        raise ValueError('sewn double-coset coordinates do not separate union states')
    gluing = []
    for key in itertools.product(range(len(words.representatives(3))), repeat=2):
        a, b = [words.representatives(3)[i] for i in key]
        i, j = union.shared_positions
        if e.charges[a[i]] != e.charges[b[j]]:
            continue
        spec = union.gluing_spec(key)
        if len(fibers[key]) != len(spec['double_cosets']):
            raise ValueError('independent double-coset count disagrees with union fiber')
        gluing.append({'patch_orbits': list(key),
                       'aligned_second_tuple': list(spec['aligned_second_tuple']),
                       'common_centralizer': spec['common_centralizer'],
                       'first_stabilizer': spec['first_stabilizer'],
                       'second_stabilizer': spec['second_stabilizer'],
                       'double_cosets': [list(c) for c in spec['double_cosets']]})
    witnesses = []
    for key, states in sorted(fibers.items()):
        for i, j in itertools.combinations(states, 2):
            # Select the actor's bank, conditioned on patch 0, retaining every
            # rule multiplicity. This is not the full-mesh averaged generator.
            left = Counter(separate[table[i]] for table in actor_tables)
            right = Counter(separate[table[j]] for table in actor_tables)
            if left == right:
                continue
            differing = [(k, left[k]-right[k]) for k in sorted(set(left) | set(right)) if left[k] != right[k]]
            rule = next(k for k, table in enumerate(actor_tables) if separate[table[i]] != separate[table[j]])
            witnesses.append({'input_union_states': [i, j], 'same_separate_patch_orbits': list(key),
                              'first_distinguishing_channel': list(actor_channels[rule]),
                              'separate_targets': [list(separate[actor_tables[rule][s]]) for s in (i, j)],
                              'actor_bank_target_count_difference': [[list(k), value] for k, value in differing],
                              'raw_links': [union.lift(union.reps[s]) for s in (i, j)]})
            break
        if witnesses:
            break
    if not witnesses:
        raise ValueError('no dynamic gluing obstruction found in exhaustive union census')
    forecast = mean_witness(union, tables, witnesses[0]['input_union_states'])
    return {'side': 6, 'patches': [0, 2], 'read_edges': union.edges, 'vertices': union.vertices,
            'patch_face_ids': [e.fan_faces[p] for p in union.patches],
            'face_vertices': [[f, list(e.geometry.faces[f][:3])] for f in forecast['face_ids']],
            'tree': [list(edge) for edge in union.tree], 'chords': union.chords,
            'patch_roots': union.patch_roots, 'connector': union.connector,
            'raw_union_connections': union.g.n**len(union.edges),
            'word_censuses': [words.census(r) for r in (3, 4, 5)],
            'triple_product_necessity': words.triple_necessity(),
            'union_representatives': [list(v) for v in union.reps],
            'channels': [list(v) for v in channels], 'union_transition_tables': tables,
            'separate_patch_orbits': [list(v) for v in separate],
            'sewn_coordinates': [list(v) for v in sewn], 'gluing_specs': gluing,
            'separate_observer_fibers': [[list(key), states] for key, states in sorted(fibers.items())],
            'gluing_witness': witnesses[0], 'mean_response_witness': forecast,
            'previous_activity_witness_resolution': previous_witness(bank, words),
            'scope': 'Exact gauge-complete local word observer and closed union quotient under all boundary-contained primitive operators. Separately complete patch orbits need not determine overlap evolution. No whole-mesh closure, emergent geometry, or quantum dynamics claimed.'}


def mean_witness(union, tables, states):
    """Exact powers of the finite operator; integers retain all no-op paths."""
    import numpy as np
    faces = sorted(set().union(*(union.e.fan_faces[p] for p in union.patches)))
    charge = lambda links: [union.e.charges[x] for x in union.e.geometry.holonomies(links)]
    values = np.array([[charge(union.lift(v))[f] for f in faces] for v in union.reps], dtype=object)
    totals = values.sum(axis=1)
    if any(not np.array_equal(totals, totals[t]) for t in tables):
        raise ValueError('contained bank fails to conserve union charge')
    rows = []
    for tick in range(9):
        difference = [int(x) for x in values[states[0]]-values[states[1]]]
        rows.append({'attempts': tick, 'difference_numerator': difference, 'denominator': len(tables)**tick})
        if any(difference):
            return {'face_ids': faces, 'contained_operator_count': len(tables), 'differences': rows,
                    'first_mean_difference_attempt': tick}
        values = sum((values[t] for t in tables), np.zeros_like(values))
    raise ValueError('gluing witness has no mean-charge distinction within the audit horizon')


def previous_witness(bank, words):
    """The new observables must distinguish the earlier (q,A) counterexample."""
    source = json.loads(Path('data/triangle-response-memory.json').read_text())['activity_closure_witness']
    e = CurrentProbe(source['side'], bank).e
    rows = []
    for patch in range(len(e.factor.pairs)):
        signatures = []
        for links in source['links']:
            code = e.factor.oracle.code(links, 1, patch)
            signatures.append(words.signature((code//36, code//6 % 6, code % 6)))
        if signatures[0] != signatures[1]:
            rows.append({'patch': patch, 'signatures': [list(s) for s in signatures]})
    if not rows:
        raise ValueError('patch words fail to resolve the prior charge/activity witness')
    return {'side': source['side'], 'source_sample_indices': source['sample_indices'], 'differing_patches': rows}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, default=Path('out/triangle-patch-observer.json'))
    args = parser.parse_args()
    result = audit()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, separators=(',', ':'))+'\n')


if __name__ == '__main__':
    main()
