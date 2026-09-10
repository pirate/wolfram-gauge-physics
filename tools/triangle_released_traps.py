#!/usr/bin/env python3
"""Release the compact trap into the actual periodic primitive dynamics.

All nine preparations have the same complete initial charge field within
each shared exterior. Backgrounds are explicit preparations, not equilibrium
claims. First exit means leaving the local one-vacancy/one-rotation quiet set,
not annihilation or loss of all gauge information.
"""
import json
import random
from collections import Counter

import numpy as np

from triangle_strip_memory import Strip


class ReleasedTraps:
    def __init__(self):
        self.s = s = Strip(7, shape='compact')
        self.x, self.oracle = s.x, s.x.e.factor.oracle
        self.mesh_faces = len(s.x.e.geometry.faces)
        self.region_edges = {e for path in s.paths for e, _ in path}
        self.pairs = sorted({p for p, r in s.channels if r == 1})
        self.fans = sorted({p for p, r in s.channels if r == 3})
        self.reactive = {code for table in s.x.tables[3:] for code, target in enumerate(table) if code != target}
        rules = np.array([r for _, r in s.channels])
        quiet = (np.sum(s.charges == 2, axis=1) == 1) & np.all(
            s.tables[rules != 0] == np.arange(s.n)[None, :], axis=0)
        target = np.all(s.charges == (0, 1, 1, 1, 1, 2, 1), axis=1)
        selected = np.flatnonzero(quiet & target)
        continuation = np.sum((s.tables[:, selected] != selected[None, :]) & quiet[s.tables[:, selected]], axis=0)
        assert sorted(continuation) == [0]*8+[2]
        selected = selected[np.argsort(-continuation, kind='stable')]
        self.preparations = [s.raw[s.representatives[i]][:] for i in selected]
        assert all(self.is_quiet(raw) for raw in self.preparations)
        self.touching = []
        self.writes = []
        for encoded in range(self.x.operator_count):
            rule, patch = divmod(encoded, self.x.supports)
            writes = ({self.oracle.pairs[patch][0][0][0]} if rule < 3
                      else self.oracle.fan.writes[patch])
            touches = bool(writes & self.region_edges)
            self.writes.append(touches)
            if touches:
                self.touching.append(encoded)

    def is_quiet(self, raw):
        q = [self.s.q[self.oracle.fan.transport(raw, path)] for path in self.s.paths]
        if q.count(0) != 1 or q.count(2) != 1:
            return False
        for patch in self.pairs:
            code = self.oracle.code(raw, 0, patch)
            if self.x.tables[1][code] != code:
                return False
        return all(self.oracle.code(raw, 1, patch) not in self.reactive for patch in self.fans)

    def embed(self, background, seed):
        rng = random.Random(seed)
        if background == 'vacuum':
            exterior = [self.s.g.identity]*len(self.preparations[0])
        elif background == 'parallel_reflections':
            exterior = [self.s.p]*len(self.preparations[0])
            # Match the patch's boundary-tree frame by an actual gauge map.
            # Pasting identity boundary links into constant-reflection links
            # would instead create a ring of vacuum faces outside the patch.
            incidence = Counter(e for path in self.s.paths for e, _ in path)
            outer = {e for e, count in incidence.items() if count == 1}
            tree = [e for e in outer if self.preparations[0][e] == self.s.g.identity]
            assert len(tree) == len(outer)-1
            vertices = {v for e in outer for v in self.x.e.geometry.edges[e]}
            frames = {min(vertices): self.s.g.identity}
            while len(frames) < len(vertices):
                old_size = len(frames)
                for edge in tree:
                    a, b = self.x.e.geometry.edges[edge]
                    if a in frames and b not in frames:
                        frames[b] = self.s.g.mul[frames[a]][self.s.p]
                    elif b in frames and a not in frames:
                        frames[a] = self.s.g.mul[frames[b]][self.s.p]
                assert len(frames) > old_size
            for edge, (a, b) in enumerate(self.x.e.geometry.edges):
                fa, fb = frames.get(a, self.s.g.identity), frames.get(b, self.s.g.identity)
                exterior[edge] = self.s.g.mul[fb][self.s.g.mul[exterior[edge]][self.s.g.inv[fa]]]
            assert all(exterior[e] == self.preparations[0][e] for e in outer)
            assert self.x.p.charge(exterior) == [1]*self.mesh_faces
        elif background == 'haar':
            exterior = [rng.randrange(self.s.g.n) for _ in self.preparations[0]]
        else:
            raise ValueError(background)
        inputs = []
        for prepared in self.preparations:
            raw = exterior[:]
            for edge in self.region_edges:
                raw[edge] = prepared[edge]
            assert self.is_quiet(raw)
            inputs.append(raw)
        charges = [self.x.p.charge(raw) for raw in inputs]
        assert all(q == charges[0] for q in charges)
        if background == 'parallel_reflections':
            assert sum(charges[0]) == self.mesh_faces
            assert all(q == 1 for f, q in enumerate(charges[0]) if f not in self.s.faces)
        return inputs, sum(charges[0])

    def initial_escape(self, background, seed):
        inputs, charge = self.embed(background, seed)
        rows = []
        contained = {r*self.x.supports+p for p, r in self.s.channels}
        for raw in inputs:
            escapes = Counter()
            for encoded in self.touching:
                moved = raw[:]
                old, new = self.x.apply(moved, encoded)
                if old != new and not self.is_quiet(moved):
                    rule = encoded//self.x.supports
                    family = 'vacancy' if rule == 0 else 'elastic' if rule < 3 else 'reaction'
                    escapes[('contained_' if encoded in contained else 'crossing_')+family] += 1
            births = deaths = 0
            for patch in range(self.x.supports):
                code = self.oracle.code(raw, 1, patch)
                before = sum(self.s.q[(code//self.s.g.n**k) % self.s.g.n] == 2 for k in range(3))
                for table in self.x.tables[3:]:
                    target = table[code]
                    after = sum(self.s.q[(target//self.s.g.n**k) % self.s.g.n] == 2 for k in range(3))
                    births += max(0, after-before)
                    deaths += max(0, before-after)
            rows.append({'quiet_exit_slots': dict(escapes),
                         'global_rotation_birth_slots': births,
                         'global_rotation_death_slots': deaths,
                         'initial_rotation_growth_per_attempt_per_face': (births-deaths)/45})
        return {'background': background, 'total_initial_charge': charge,
                'full_operator_slots': self.x.operator_count, 'escape_slot_counts': rows}

    def experiment(self, background, trials=512, cap=32):
        durations, censored, totals = [], [], []
        for trial in range(trials):
            inputs, total = self.embed(background, 91523000+trial)
            rng = random.Random(92624000+trial)
            ticks = np.full(9, cap*self.mesh_faces, dtype=int)
            alive = list(range(9))
            for tick in range(1, cap*self.mesh_faces+1):
                encoded = rng.randrange(self.x.operator_count)
                remaining = []
                for i in alive:
                    old, new = self.x.apply(inputs[i], encoded)
                    if self.writes[encoded] and old != new and not self.is_quiet(inputs[i]):
                        ticks[i] = tick
                    else:
                        remaining.append(i)
                alive = remaining
                if not alive:
                    break
            # All operators, including those wholly outside the region, were
            # evolved until each trajectory's first exit. No boundary is held.
            durations.append(ticks/self.mesh_faces)
            censored.append(len(alive))
            totals.append(total)
        times = np.array(durations)
        slow, others = times[:, 0], times[:, 1:].mean(axis=1)
        delta = slow-others
        se = lambda x: float(x.std(ddof=1)/np.sqrt(len(x)))
        return {'background': background, 'trials': trials, 'cap_attempts_per_face': cap,
                'censored_trajectories': sum(censored),
                'initial_total_charge_range': [min(totals), max(totals)],
                'former_two_state_trap_mean_exit': float(slow.mean()), 'former_two_state_trap_SE': se(slow),
                'former_singletons_mean_exit': float(others.mean()), 'former_singletons_SE': se(others),
                'paired_mean_difference': float(delta.mean()), 'paired_difference_SE': se(delta),
                'meaning': 'Restricted first-exit mean if censored; local quiet-set exit, not gauge-memory erasure.'}


if __name__ == '__main__':
    r = ReleasedTraps()
    for background in ('vacuum', 'parallel_reflections', 'haar'):
        print(json.dumps(r.initial_escape(background, 91523000)), flush=True)
        print(json.dumps(r.experiment(background)), flush=True)
