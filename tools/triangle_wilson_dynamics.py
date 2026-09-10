#!/usr/bin/env python3
"""Wilson characters on actual boundary links under unchanged primitive updates.

Exact stationary starts, independent schedules, fixed attempted-time sampling,
and trajectory-clustered errors. No fit selects the area-law prediction.
"""
import argparse
import json
import random
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

from triangle_elastic_scattering import ElasticExperiment
from triangle_reference import AxialReference
from triangle_wilson_reference import WilsonReference


class LoopObserver:
    def __init__(self, side):
        self.side = side
        self.bank = json.loads(Path('data/triangle-feedback.json').read_text())['bank']
        elastic = json.loads(Path('data/triangle-elastic-scattering.json').read_text())['pair_census']['elastic_tables']
        self.x = ElasticExperiment(side, self.bank, elastic)
        self.g, self.faces = self.x.e.geometry.group, 2*side*side
        self.theory = WilsonReference(self.g, self.x.e.charges)
        self.mul, self.inv = np.array(self.g.mul), np.array(self.g.inv)
        self.characters = np.array(self.theory.characters, dtype=float)/np.array(self.theory.dimensions)[:, None]
        self.shapes, self.paths, self.regions = [], [], []
        shapes = [('triangle', [(0, 0, 0)])]
        shapes += [(f'{w}x{h}_cells', [(x, y, t) for y in range(h) for x in range(w) for t in range(2)])
                   for w, h in ((1, 1), (2, 1), (3, 1), (4, 1), (2, 2), (3, 2))]
        geometry = self.x.e.geometry
        for name, offsets in shapes:
            paths, regions = [], []
            for y in range(side):
                for x in range(side):
                    region = [2*(((y+dy) % side)*side+(x+dx) % side)+t for dx, dy, t in offsets]
                    boundary = set()
                    for f in region:
                        face = geometry.faces[f]
                        for u, v in zip(face, face[1:]):
                            if (v, u) in boundary:
                                boundary.remove((v, u))
                            else:
                                boundary.add((u, v))
                    successor = dict(boundary)
                    assert len(successor) == len(boundary)
                    start = min(successor)
                    u, path = start, []
                    for _ in boundary:
                        v = successor[u]
                        path.append((geometry.ids[tuple(sorted((u, v)))], u > v))
                        u = v
                    assert u == start and len({e for e, _ in path}) == len(boundary)
                    paths.append(path); regions.append(region)
            self.shapes.append({'name': name, 'area_faces': len(offsets), 'perimeter_edges': len(paths[0]),
                                'translated_loops_per_frame': side*side})
            self.paths.append(np.array(paths, dtype=int))
            self.regions.append(np.array(regions, dtype=int))

    def holonomies(self, links, shape):
        states = np.asarray(links, dtype=int)
        if states.ndim == 1:
            states = states[None, :]
        paths = self.paths[shape]
        product = np.full((len(states), len(paths)), self.g.identity, dtype=int)
        for step in range(paths.shape[1]):
            values = states[:, paths[:, step, 0]]
            values = np.where(paths[:, step, 1], self.inv[values], values)
            product = self.mul[values, product]
        return product

    def observe(self, states):
        return np.stack([self.characters[1:, self.holonomies(states, shape)].mean(axis=-1).T
                         for shape in range(len(self.shapes))], axis=1)

    def hidden_loop_witness(self, reference, seed):
        oracle = self.x.e.factor.oracle
        for draw in range(32):
            raw = reference.sample(random.Random(seed+draw), 'nonabelian_reflections')['links']
            before_q = self.x.p.charge(raw)
            before = self.holonomies(raw, 1)[0]
            for translation, path in enumerate(self.paths[1]):
                reads = {int(edge) for edge, _ in path}
                for patch in range(self.x.supports):
                    if oracle.pairs[patch][0][0][0] not in reads:
                        continue
                    moved = raw[:]
                    code, target = self.x.apply(moved, self.x.supports+patch)
                    if code == target:
                        continue
                    after = int(self.holonomies(moved, 1)[0, translation])
                    if self.characters[2, before[translation]] != self.characters[2, after]:
                        assert self.x.p.charge(moved) == before_q
                        return {'seed': seed+draw, 'shape': self.shapes[1], 'translation_index': translation,
                                'path': path.tolist(), 'patch': patch, 'rule': 1, 'raw_tuple_before_after': [code, target],
                                'boundary_holonomy_before_after': [int(before[translation]), after],
                                'standard_character_before_after': [float(self.characters[2, before[translation]]),
                                                                    float(self.characters[2, after])],
                                'all_face_charges_unchanged': True, 'face_charges': before_q,
                                'initial_links': raw, 'final_links': moved}
        raise RuntimeError('no hidden loop witness found in the stated bounded search')

    def trajectory(self, initial, seed, duration, interval):
        rng = random.Random(seed)
        attempts = duration*self.faces
        schedule = [rng.randrange(self.x.operator_count) for _ in range(attempts)]
        row = self.x.compiled_output([initial], schedule, attempts, raw_events=True)['runs'][1]
        assert row['mode'] == 'combined'
        frames, cursor, current = [initial], 0, initial
        events = row['events']
        for tick in range(interval*self.faces, attempts+1, interval*self.faces):
            while cursor < len(events) and events[cursor][0] <= tick:
                current = events[cursor][5]
                cursor += 1
            frames.append(current)
        assert current == row['final_links']
        values = self.observe(frames)
        # Each trajectory contributes one unit of weight. Translations and
        # consecutive frames are never treated as independent replicates.
        return {'mean': values[1:].mean(axis=0).tolist(), 'initial': values[0].tolist(),
                'final': values[-1].tolist(), 'changing_events': len(events)}


def measure(side, charge, trials, duration, interval, seed, workers=4):
    started = time.monotonic()
    observer = LoopObserver(side)
    reference = AxialReference(side, observer.bank, charge)
    initials = [reference.sample(random.Random(seed+2*i), 'nonabelian_reflections')['links'] for i in range(trials)]
    # The sign character is exactly the parity of the enclosed scalar charge.
    # The standard character instead retains nonabelian alignment information.
    for raw in initials[:4]:
        q = np.array(observer.x.p.charge(raw))
        for shape, region in enumerate(observer.regions):
            assert np.array_equal(observer.characters[1, observer.holonomies(raw, shape)[0]], (-1)**q[region].sum(axis=1))
    run = lambda i: observer.trajectory(initials[i], seed+2*i+1, duration, interval)
    records = []
    with ThreadPoolExecutor(max_workers=workers) as pool:
        for i, record in enumerate(pool.map(run, range(trials)), 1):
            records.append(record)
            if i % 16 == 0:
                print(json.dumps({'progress': i, 'trials': trials, 'side': side, 'charge': charge,
                                  'elapsed_seconds': time.monotonic()-started}), flush=True)
    means = np.array([r['mean'] for r in records])
    final = np.array([r['final'] for r in records])
    initial = np.array([r['initial'] for r in records])
    se = lambda a: a.std(axis=0, ddof=1)/np.sqrt(len(a))
    rows = []
    for i, shape in enumerate(observer.shapes):
        exact = [observer.theory.mean(observer.faces, charge, shape['area_faces'], t) for t in (1, 2)]
        rows.append({**shape, 'exact_mean': [float(v) for v in exact], 'exact_mean_rational': [str(v) for v in exact],
                     'mean': means[:, i].mean(axis=0).tolist(), 'trajectory_standard_error': se(means[:, i]).tolist(),
                     'initial_mean': initial[:, i].mean(axis=0).tolist(),
                     'final_mean': final[:, i].mean(axis=0).tolist(),
                     'final_minus_initial_mean': (final[:, i]-initial[:, i]).mean(axis=0).tolist(),
                     'final_minus_initial_standard_error': se(final[:, i]-initial[:, i]).tolist()})
    difference = means[:, 4]-means[:, 5]
    equal_perimeter_difference = means[:, 4]-means[:, 6]
    return {'side': side, 'faces': observer.faces, 'charge': charge, 'trials': trials, 'duration_attempts_per_face': duration,
            'observation_interval_attempts_per_face': interval, 'seed': seed, 'representation_order': ['sign', 'standard_normalized'],
            'combined_forward_attempts': trials*duration*observer.faces,
            'combined_changing_events': sum(r['changing_events'] for r in records),
            'bulk_prediction': observer.theory.bulk(charge/observer.faces), 'loops': rows,
            'same_area_different_perimeter': {'shapes': [4, 5], 'area_faces': 8, 'perimeters': [10, 8],
                                             'paired_mean_difference': difference.mean(axis=0).tolist(),
                                             'paired_standard_error': se(difference).tolist(), 'exact_difference': [0, 0]},
            'same_perimeter_different_area': {'shapes': [4, 6], 'perimeter_edges': 10, 'areas': [8, 12],
                                             'paired_mean_difference': equal_perimeter_difference.mean(axis=0).tolist(),
                                             'paired_standard_error': se(equal_perimeter_difference).tolist(),
                                             'exact_difference': [rows[4]['exact_mean'][i]-rows[6]['exact_mean'][i] for i in range(2)]},
            'trajectory_records': records, 'elapsed_seconds': time.monotonic()-started,
            'scope': 'Stationary-start actual compiled primitive evolution. No equilibration claim, fitted force, physical string tension, or quantum dynamics.'}


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--trials', type=int, default=128)
    parser.add_argument('--duration', type=int, default=128)
    parser.add_argument('--confirmation', action='store_true')
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    observer = LoopObserver(6)
    ref = AxialReference(6, observer.bank, 72)
    report = {'exact_small_count': observer.theory.direct_small_torus_count(),
              'hidden_loop_witness': observer.hidden_loop_witness(ref, 252217001), 'experiments': []}
    cases = ((6, 18, 296661001), (8, 32, 307772001), (6, 36, 318883001)) if args.confirmation else (
             (6, 18, 263328001), (8, 32, 274439001), (6, 36, 285550001))
    report['confirmation'] = args.confirmation
    for side, charge, seed in cases:
        report['experiments'].append(measure(side, charge, args.trials, args.duration, 4, seed))
        Path(args.output).write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({'output': args.output, 'complete': True}), flush=True)
