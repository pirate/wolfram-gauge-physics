#!/usr/bin/env python3
"""Exact whole-mesh tree-gauge evolution with lazy nonabelian subtree frames.

The original primitive tables and attempted-update scheduler are unchanged.
Only their connection coordinates change. The remaining common-root conjugation
is retained internally; gauge-invariant signatures are separate observations.
"""
import argparse
import hashlib
import json
import random
from collections import Counter
from pathlib import Path

from triangle_charge_current import CurrentProbe


class LazyFrames:
    """Range LEFT multiplication and point lookup; multiplication need not commute."""
    def __init__(self, size, group):
        if type(size) is not int or size < 1:
            raise ValueError('frame tree needs a positive integer size')
        self.size, self.g = size, group
        self.tags = [group.identity]*(4*size)
        self.reset_counts()

    def reset_counts(self):
        self.counts = Counter()

    def left_multiply(self, begin, end, value):
        if any(type(x) is not int for x in (begin, end, value)) or not 0 <= begin < end <= self.size or not 0 <= value < self.g.n:
            raise ValueError('invalid frame interval or group element')
        if value == self.g.identity:
            return
        self.counts['range_updates'] += 1
        def visit(node, lo, hi):
            if end <= lo or hi <= begin:
                return
            self.counts['range_nodes'] += 1
            if begin <= lo and hi <= end:
                self.tags[node] = self.g.mul[value][self.tags[node]]
                return
            pending = self.tags[node]
            if pending != self.g.identity:
                for child in (2*node, 2*node+1):
                    self.tags[child] = self.g.mul[pending][self.tags[child]]
                self.tags[node] = self.g.identity
            mid = (lo+hi)//2
            visit(2*node, lo, mid)
            visit(2*node+1, mid, hi)
        visit(1, 0, self.size)

    def point(self, index):
        if type(index) is not int or not 0 <= index < self.size:
            raise ValueError('frame point is outside the tree')
        self.counts['point_queries'] += 1
        node, lo, hi, value = 1, 0, self.size, self.g.identity
        while True:
            self.counts['point_nodes'] += 1
            # An unpushed ancestor tag is newer than all descendant tags.
            value = self.g.mul[value][self.tags[node]]
            if hi-lo == 1:
                return value
            mid = (lo+hi)//2
            if index < mid:
                node, hi = 2*node, mid
            else:
                node, lo = 2*node+1, mid

    def all_points(self):
        """Explicit O(V) snapshot, never part of a primitive's update path."""
        result = [self.g.identity]*self.size
        stack = [(1, 0, self.size, self.g.identity)]
        while stack:
            node, lo, hi, outer = stack.pop()
            value = self.g.mul[outer][self.tags[node]]
            if hi-lo == 1:
                result[lo] = value
            else:
                mid = (lo+hi)//2
                stack.extend(((2*node, lo, mid, value), (2*node+1, mid, hi, value)))
        return result


class ConnectionView:
    """One primitive's bounded read cache and simultaneous write buffer."""
    def __init__(self, state):
        self.state, self.pending, self.frames = state, {}, {}

    def __len__(self):
        return len(self.state.e.geometry.edges)

    def frame(self, vertex):
        if vertex not in self.frames:
            self.frames[vertex] = self.state.frames.point(self.state.entry[vertex])
        return self.frames[vertex]

    def __getitem__(self, edge):
        s = self.state
        if edge in self.pending:
            return self.pending[edge]
        if edge in s.tree_edges:
            return s.g.identity
        u, v = s.e.geometry.edges[edge]
        return s.g.mul[self.frame(v)][s.g.mul[s.chords[edge]][s.g.inv[self.frame(u)]]]

    def __setitem__(self, edge, value):
        self.pending[edge] = value


class GaugeConnection:
    def __init__(self, experiment, links, forest=None):
        self.e, self.g = experiment, experiment.geometry.group
        forest = experiment.forest if forest is None else forest
        if forest.edges != tuple(experiment.geometry.edges) or forest.g.elements != self.g.elements:
            raise ValueError('coordinate forest differs from the primitive connection graph or group')
        self.forest = forest
        if len(forest.components) != 1:
            raise ValueError('this evolution adapter requires a connected fixed mesh')
        component = forest.components[0]
        self.root = component['root']
        self.chord_order = component['chords']
        self.chords = dict(zip(self.chord_order, forest.based_loops(links)[0]))
        self.tree_edges, children = {}, {v: [] for v in component['queue']}
        self.depth = {self.root: 0}
        for v in component['queue'][1:]:
            u, edge, reverse = component['parents'][v]
            self.tree_edges[edge] = (v, reverse)
            children[u].append(v)
            self.depth[v] = self.depth[u]+1
        self.entry, self.exit = {}, {}
        stack = [(self.root, False)]
        while stack:
            vertex, exiting = stack.pop()
            if exiting:
                self.exit[vertex] = len(self.entry)
            else:
                self.entry[vertex] = len(self.entry)
                stack.append((vertex, True))
                stack.extend((v, False) for v in reversed(children[vertex]))
        self.frames = LazyFrames(len(self.entry), self.g)

    def step(self, encoded):
        e, g = self.e, self.g
        supports = len(e.factor.pairs)
        if type(encoded) is not int or not 0 <= encoded < len(e.tables)*supports:
            raise ValueError('invalid rooted operator')
        rule, patch = divmod(encoded, supports)
        view = ConnectionView(self)
        self.frames.reset_counts()
        code = e.factor.oracle.code(view, rule, patch)
        target = e.tables[rule][code]
        if code != target:
            e.factor.oracle.update(view, rule, patch, e.tables[rule])
            # First write chords in the current frame. The subsequent tree
            # restoration must transform these NEW values, not their old ones.
            for edge, value in view.pending.items():
                if edge in self.chords:
                    u, v = e.geometry.edges[edge]
                    self.chords[edge] = g.mul[g.inv[view.frame(v)]][g.mul[value][view.frame(u)]]
            tree_writes = [(self.depth[self.tree_edges[edge][0]], edge, value)
                           for edge, value in view.pending.items() if edge in self.tree_edges]
            # Deepest first: restoring a descendant never changes an ancestor
            # write; a later ancestor conjugation preserves restored identities.
            for _, edge, value in sorted(tree_writes, reverse=True):
                child, reverse = self.tree_edges[edge]
                correction = value if reverse else g.inv[value]
                self.frames.left_multiply(self.entry[child], self.exit[child], correction)
        counts = dict(self.frames.counts)
        if counts.get('point_queries', 0) > (5 if rule else 4) or counts.get('range_updates', 0) > (2 if rule else 1):
            raise ValueError('primitive gauge work escaped its bounded local support')
        return {'changed': code != target, 'rule': rule, 'patch': patch,
                'tree_writes': sum(edge in self.tree_edges for edge in view.pending),
                'chord_writes': sum(edge in self.chords for edge in view.pending), 'work': counts}

    def loops(self):
        frames = self.frames.all_points()
        result = []
        for edge in self.chord_order:
            u, v = self.e.geometry.edges[edge]
            result.append(self.g.mul[frames[self.entry[v]]][
                self.g.mul[self.chords[edge]][self.g.inv[frames[self.entry[u]]]]])
        return tuple(result)

    def materialize(self):
        result = [self.g.identity]*len(self.e.geometry.edges)
        for edge, value in zip(self.chord_order, self.loops()):
            result[edge] = value
        return result


def trial(side, bank, attempts=1000):
    p = CurrentProbe(side, bank)
    e, rng = p.e, random.Random(554020+side)
    initial = [rng.randrange(e.geometry.group.n) for _ in e.geometry.edges]
    schedule = [rng.randrange(len(e.tables)*len(e.factor.pairs)) for _ in range(attempts)]
    run = e.compiled_bank([initial], schedule, attempts, capture_links=True)['runs'][1]
    events = {event[0]: raw for event, raw in zip(run['events'], run['raw_event_links'])}
    state, raw, digest = GaugeConnection(e, initial), initial, hashlib.sha256()
    totals, maxima, tree_patterns = Counter(), Counter(), Counter()
    for tick, operator in enumerate(schedule, 1):
        step = state.step(operator)
        if step['changed'] != (tick in events):
            raise ValueError('tree-gauge update changes the actual activity schedule')
        raw = events.get(tick, raw)
        # This compares every rooted loop, not just conjugacy classes or charge.
        expected = e.forest.based_loops(raw)[0]
        if state.loops() != expected:
            raise ValueError('incremental tree gauge differs from complete C++ connection loops')
        digest.update(bytes(expected))
        totals.update(step['work'])
        for key, value in step['work'].items():
            maxima[key] = max(maxima[key], value)
        if step['changed']:
            tree_patterns[step['tree_writes'], step['chord_writes']] += 1
    final = state.loops()
    for operator in reversed(schedule):
        state.step(operator)
    if state.loops() != e.forest.based_loops(initial)[0]:
        raise ValueError('reverse schedule fails in the complete tree-gauge state')
    return {'side': side, 'vertices': e.geometry.size, 'links': len(e.geometry.edges),
            'independent_loops': len(final), 'attempts': attempts, 'events': len(events),
            'loop_history_sha256': digest.hexdigest(), 'final_loops': list(final),
            'work_totals': dict(totals), 'work_maxima': dict(maxima),
            'write_patterns': [[a, b, n] for (a, b), n in sorted(tree_patterns.items())]}


def audit():
    bank = json.loads(Path('data/triangle-feedback.json').read_text())['bank']
    records = [trial(side, bank) for side in (3, 6, 12, 24)]
    return {'records': records,
            'scope': 'Exact original whole-mesh primitive evolution in root-fixed tree gauge. Subtree frames are coordinate changes, not physical nonlocal updates. Full loop snapshots remain linear-size observations; no speedup over raw-link evolution or emergent quantum physics is claimed.'}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, default=Path('out/triangle-lazy-gauge.json'))
    args = parser.parse_args()
    data = audit()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(data, separators=(',', ':'))+'\n')


if __name__ == '__main__':
    main()
