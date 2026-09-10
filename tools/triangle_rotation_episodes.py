#!/usr/bin/env python3
"""Event-derived rotation episodes and exact transport through native face frames.

Episode labels are annotations. Ending an episode at a reversible reaction does
not assert that its information is erased or that physical particles exist.
"""
import itertools
import json
import random
from pathlib import Path

import numpy as np

from triangle_complete_memory import CompleteObserver
from triangle_elastic_scattering import ElasticExperiment
from triangle_encounter_memory import EncounterObserver
from triangle_reference import AxialReference
from triangle_relational_memory import replay_samples


class RotationEpisodes:
    def __init__(self, experiment, raw):
        self.x, self.g = experiment, experiment.e.geometry.group
        self.oracle = experiment.e.factor.oracle.fan
        self.tags = [-1]*experiment.geometry.faces
        self.records = []
        self.hops = self.births = self.deaths = self.return_hops = self.return_flips = 0
        self.charge = experiment.p.charge(raw)
        self.prefixes = []
        for a, b, fa, fb in experiment.e.geometry.pair_specs:
            root = a[0]
            if b[0] != root:
                raise ValueError('pair loops are not based in the same frame')
            self.prefixes.append([self.oracle.face_paths[f][:experiment.e.geometry.faces[f][:-1].index(root)]
                                  for f in (fa, fb)])
        for face, q in enumerate(self.charge):
            if q == 2:
                self.create(raw, face, 0)
        self.initial_count = len(self.records)
        self.births = 0

    def holonomy(self, raw, face):
        return self.oracle.transport(raw, self.oracle.face_paths[face])

    def create(self, raw, face, tick):
        if self.tags[face] != -1:
            raise ValueError('rotation birth would overwrite an existing episode')
        self.tags[face] = len(self.records)
        self.records.append({'birth_attempt': tick, 'birth_face': face,
                             'birth_holonomy': self.holonomy(raw, face),
                             'face': face, 'transport': self.g.identity,
                             'hops': 0, 'death_attempt': None})
        self.births += 1

    def event(self, before, raw, tick, op, code, target):
        x, g = self.x, self.g
        rule, patch = divmod(op, x.supports)
        arity = 2 if rule < 3 else 3
        decode = lambda value: [(value//g.n**i) % g.n for i in reversed(range(arity))]
        old, new = decode(code), decode(target)
        old_q, new_q = [x.e.charges[v] for v in old], [x.e.charges[v] for v in new]
        faces = x.p.supports[int(rule >= 3)][patch]
        if [self.charge[f] for f in faces] != old_q:
            raise ValueError('episode charges disagree with the actual primitive input')
        if rule < 3 and 2 in old_q:
            if rule != 0 or sorted(old_q) != [0, 2] or new_q != old_q[::-1]:
                raise ValueError('unexpected rotation-containing pair law')
            source, dest = old_q.index(2), new_q.index(2)
            f, h = faces[source], faces[dest]
            tag = self.tags[f]
            if tag < 0 or self.tags[h] != -1:
                raise ValueError('rotation transport has invalid episode occupancy')
            # A: old native face basepoint -> common event root, before update.
            # B: new native face basepoint -> common event root, after update.
            # Native-frame transport is B^{-1} A, in this noncommutative order.
            a = self.oracle.transport(before, self.prefixes[patch][source])
            b = self.oracle.transport(raw, self.prefixes[patch][dest])
            transport = g.mul[g.inv[b]][a]
            old_h, new_h = self.holonomy(before, f), self.holonomy(raw, h)
            if new_h != g.mul[transport][g.mul[old_h][g.inv[transport]]]:
                raise ValueError('native holonomy did not follow the derived event connector')
            record = self.records[tag]
            record['transport'] = g.mul[transport][record['transport']]
            record['hops'] += 1; record['face'] = h
            self.tags[f], self.tags[h] = -1, tag
            self.hops += 1
            if h == record['birth_face']:
                self.return_hops += 1
                flipped = new_h != record['birth_holonomy']
                if flipped != (x.e.charges[record['transport']] % 2 == 1):
                    raise ValueError('closed-return inversion disagrees with transport parity')
                self.return_flips += int(flipped)
        elif rule >= 3 and 2 in old_q:
            if sorted(old_q) != [0, 1, 2] or new_q != [1, 1, 1]:
                raise ValueError('unexpected rotation-destroying reaction')
            f = faces[old_q.index(2)]
            tag = self.tags[f]
            if tag < 0:
                raise ValueError('reaction consumed no rotation episode')
            self.records[tag]['death_attempt'] = tick
            self.records[tag]['face'] = -1
            self.tags[f] = -1; self.deaths += 1
        elif rule >= 3 and 2 in new_q:
            if old_q != [1, 1, 1] or sorted(new_q) != [0, 1, 2]:
                raise ValueError('unexpected rotation-creating reaction')
            self.create(raw, faces[new_q.index(2)], tick)
        for face, q in zip(faces, new_q):
            self.charge[face] = q

    def verify(self, raw):
        g = self.g
        if self.charge != self.x.p.charge(raw):
            raise ValueError('incremental episode charges differ from native face holonomies')
        live = [tag for tag in self.tags if tag >= 0]
        if len(set(live)) != len(live) or len(live) != self.initial_count+self.births-self.deaths:
            raise ValueError('episode accounting lost injectivity or birth/death balance')
        for face, (q, tag) in enumerate(zip(self.charge, self.tags)):
            if (q == 2) != (tag >= 0):
                raise ValueError('episode occupancy differs from rotation charge')
            if tag >= 0:
                r = self.records[tag]
                expected = g.mul[r['transport']][g.mul[r['birth_holonomy']][g.inv[r['transport']]]]
                if r['face'] != face or r['death_attempt'] is not None or self.holonomy(raw, face) != expected:
                    raise ValueError('cumulative native-frame holonomy transport failed')
        return {'initial_alive': sum(tag < self.initial_count for tag in live),
                'born_alive': sum(tag >= self.initial_count for tag in live),
                'births': self.births, 'deaths': self.deaths, 'rotation_hops': self.hops,
                'return_hops': self.return_hops, 'return_inversions': self.return_flips}


class EpisodeObserver(EncounterObserver):
    def __init__(self, experiment, complete):
        super().__init__(experiment, complete)
        self.episodes = None
        self.faces = np.array(experiment.e.fan_faces)

    def event(self, raw, tick, op, code, target):
        self.episodes.event(self.shadow, raw, tick, op, code, target)
        super().event(raw, tick, op, code, target)

    def measure(self, raw):
        sample = super().measure(raw)
        if self.episodes is None:
            self.episodes = RotationEpisodes(self.x, raw)
        census = self.episodes.verify(raw)
        return {**sample, 'episode_tags': np.array(self.episodes.tags)[self.faces], 'episode_census': census}

    def returns(self, first, later):
        # Classification concerns the original patch's three uninterrupted
        # rotation episodes. A reaction crossing does not imply information loss.
        initial = first['episode_tags']
        current = later['episode_tags']
        both = np.all(initial >= 0, axis=1) & np.all(current >= 0, axis=1)
        untouched = later['history_labels'] == 0
        exact = np.all(initial == current, axis=1)
        cohort = np.all(np.sort(initial, axis=1) == np.sort(current, axis=1), axis=1)
        origin = np.all((current >= 0) & (current < self.episodes.initial_count), axis=1)
        labels = np.where(untouched, 0, np.where(exact, 1, np.where(cohort, 2, np.where(origin, 3, 4))))
        same = first['states'] == later['states']
        counts = np.bincount(labels[both], minlength=5)
        returns = np.bincount(labels[both & same], minlength=5)
        contributions = (4*returns-counts)/(3*self.x.supports*float(self.rotation_probability))
        ordinary, _ = self.correlate(first, later)
        if not np.isclose(sum(contributions), sum(ordinary['local_contributions'][1]), atol=1e-11):
            raise ValueError('episode return bins do not sum to complete rotation memory')
        # The three orthogonal 222 contrasts are pairwise relative orientations.
        # Native fan order is opposite to the canonical tuple's table order.
        pairs = np.array([(0, 1), (0, 2), (1, 2)])
        reps = np.array(self.complete.reference.reps)
        a, b = reps[first['states']][:, ::-1], reps[later['states']][:, ::-1]
        phi = lambda v: 2*(v[:, pairs[:, 0]] == v[:, pairs[:, 1]]).astype(int)-1
        products = phi(a)*phi(b)
        if not np.array_equal(products[both].sum(axis=1), 4*same[both].astype(int)-1):
            raise ValueError('pairwise orientation contrasts disagree with the complete orbit kernel')
        retained = np.all(initial[:, pairs] == current[:, pairs], axis=2)
        selected = both & (labels == 4)
        pair_counts = [int(np.sum(mask[selected])) for mask in (retained, ~retained)]
        pair_sums = [int(np.sum((products*mask)[selected])) for mask in (retained, ~retained)]
        pair_contributions = np.array(pair_sums)/(3*self.x.supports*float(self.rotation_probability))
        if not np.isclose(sum(pair_contributions), contributions[4], atol=1e-11):
            raise ValueError('pair decomposition does not recover reaction-born-patch memory')
        return {'both_rotation_endpoint_counts': counts.tolist(), 'same_orbit_endpoint_counts': returns.tolist(),
                'local_return_contributions': contributions.tolist(), 'episode_census': later['episode_census'],
                'reaction_born_patch_pair_counts': pair_counts,
                'reaction_born_patch_pair_signed_sums': pair_sums,
                'reaction_born_patch_pair_contributions': pair_contributions.tolist()}


def pilot(side=3, attempts=10000):
    bank = json.loads(Path('data/triangle-feedback.json').read_text())['bank']
    elastic = json.loads(Path('data/triangle-elastic-scattering.json').read_text())['pair_census']['elastic_tables']
    x = ElasticExperiment(side, bank, elastic)
    c = CompleteObserver(x)
    observer = EpisodeObserver(x, c)
    raw = AxialReference(side, bank, x.geometry.faces).sample(random.Random(772991), 'nonabelian_reflections')['links']
    rng = random.Random(881992)
    samples, check = replay_samples(x, raw, [rng.randrange(x.operator_count) for _ in range(attempts)],
                                    [0, attempts], observer, observer.event)
    return {'returns': observer.returns(samples[0], samples[-1]), **check}


def ensemble(side=12, trials=64, times=(0, 1, 2, 4, 8, 16),
             initial_seed=58112000, schedule_seed=68212000):
    """Independent exact-reference trajectories, with uncertainty across runs."""
    bank = json.loads(Path('data/triangle-feedback.json').read_text())['bank']
    elastic = json.loads(Path('data/triangle-elastic-scattering.json').read_text())['pair_census']['elastic_tables']
    x = ElasticExperiment(side, bank, elastic)
    complete = CompleteObserver(x)
    reference = AxialReference(side, bank, x.geometry.faces)
    rows = []
    for index in range(trials):
        observer = EpisodeObserver(x, complete)
        raw = reference.sample(random.Random(initial_seed+index), 'nonabelian_reflections')['links']
        rng = random.Random(schedule_seed+index)
        ticks = [t*x.geometry.faces for t in times]
        samples, check = replay_samples(x, raw, [rng.randrange(x.operator_count) for _ in range(ticks[-1])],
                                        ticks, observer, observer.event)
        rows.append({'index': index, 'observations': [observer.returns(samples[0], s) for s in samples], **check})
    def estimate(key):
        a = np.array([[s[key] for s in r['observations']] for r in rows])
        return {'mean': a.mean(axis=0).tolist(), 'standard_error': (a.std(axis=0, ddof=1)/trials**.5).tolist()}
    return {'side': side, 'trials': trials, 'times': times, 'initial_seed': initial_seed,
            'schedule_seed': schedule_seed, 'episode_bins': ['untouched', 'same_episodes_same_faces',
            'same_cohort_permuted', 'other_initial_episodes', 'at_least_one_reaction_born_episode'],
            'reaction_born_pair_bins': ['same_two_original_episodes_same_faces', 'at_least_one_different_episode'],
            'estimates': {key: estimate(key) for key in ('local_return_contributions', 'reaction_born_patch_pair_contributions')},
            'replicates': rows}


if __name__ == '__main__':
    print(json.dumps(pilot(), indent=2))
