#!/usr/bin/env python3
"""Complete three-loop relational memory with an exact conditional reference.

All functions of the 49 local gauge states are retained, after projection off
the entire charge field. Whitening is an observation metric, never a new rule.
"""
import argparse
import itertools
import json
import math
import random
from collections import Counter, defaultdict
from fractions import Fraction
from functools import cached_property, lru_cache
from pathlib import Path

import numpy as np

from triangle_elastic_scattering import ElasticExperiment
from triangle_patch_observer import WordObserver
from triangle_reference import AxialReference, field_weight, fraction, reference
from triangle_relational_memory import replay_samples
from triangle_spatial_memory import Moments, Translations


def rotation_gate(experiment):
    """Exhaust the primitive obstruction to changing a rotation without a vacancy."""
    q, g = experiment.e.charges, experiment.e.geometry.group
    rows = []
    for rule, table in enumerate(experiment.tables):
        arity = 2 if rule < 3 else 3
        tuples = list(itertools.product(range(g.n), repeat=arity))
        enabled, rotation_enabled, patterns = 0, 0, Counter()
        for code, values in enumerate(tuples):
            after = tuples[table[code]]
            enabled += values != after
            if values == after or not any(q[h] == 2 for h in values):
                continue
            rotation_enabled += 1
            pattern = tuple(q[h] for h in values)
            patterns[pattern] += 1
            if 0 not in pattern:
                raise ValueError('an active rotation-containing primitive lacks a vacancy')
            for a, b in zip(values, after):
                if q[a] == 2 and (q[b] == 2 or a == b):
                    raise ValueError('an active primitive retains its input rotation face')
        rows.append({'rule_slot': rule, 'arity': arity, 'raw_inputs': len(tuples),
                     'enabled_inputs': enabled, 'enabled_rotation_inputs': rotation_enabled,
                     'enabled_rotation_charge_patterns': [{'pattern': p, 'inputs': n} for p, n in sorted(patterns.items())]})
    return {'primitive_rows': rows,
            'statement': 'Every changing primitive containing a rotation face also contains a vacancy, and each input rotation face leaves its charge-two sector. In a vacancy-free rotation region, the first change must therefore enter through a support meeting its boundary.',
            'scope': 'Exact rule-bank kinetic constraint. Boundary updates can erode a rotation region; this is not binding or a protected propagating object.'}


class CompleteReference:
    def __init__(self, group, charges, faces, charge):
        self.g, self.q, self.faces = group, charges, faces
        self.words = WordObserver(group, charges)
        self.reps = self.words.representatives(3)
        self.ids = {rep: i for i, rep in enumerate(self.reps)}
        self.tuples = list(itertools.product(range(group.n), repeat=3))
        self.codes = np.array([self.ids[self.words.canonical(v)] for v in self.tuples])
        self.orbit_sizes = Counter(self.codes)
        self.patterns = [tuple(charges[h] for h in rep) for rep in self.reps]
        self.buckets = defaultdict(list)
        for i, pattern in enumerate(self.patterns):
            self.buckets[pattern].append(i)
        self.counts = reference(faces, charge, 'nonabelian_reflections')
        # Each charge bucket has one exact sum-to-zero constraint. Differences
        # from its last category give a deterministic, nonredundant basis.
        columns = []
        for pattern, ids in sorted(self.buckets.items()):
            for i in ids[:-1]:
                col = [0]*len(self.reps)
                col[i], col[ids[-1]] = 1, -1
                columns.append(col)
        self.contrasts = np.array(columns, dtype=int).T
        self.covariance = self.exact_covariance()
        h = self.contrasts.astype(object)
        self.gram = h.T@self.covariance@h
        self.dimension = len(columns)
        if (len(self.reps), len(self.buckets), self.dimension) != (49, 27, 22):
            raise ValueError('unexpected three-loop gauge/charge dimensions')

    @cached_property
    def whitener(self):
        # Tiny counting presentations may not populate all 22 contrasts. They
        # still have valid exact marginals, but cannot use this full-rank metric.
        # Do not cure a singular reference by adding an arbitrary regularizer.
        if any(self.gram[i, i] == 0 for i in range(self.dimension)):
            raise ValueError('reference does not populate all relational contrasts')
        lower = np.linalg.cholesky(np.array(self.gram, dtype=float))
        return self.contrasts@np.linalg.solve(lower.T, np.eye(self.dimension))

    def completion(self, values, populations):
        """Fixed based local tuple, remaining fixed charges, full-S3 handles.

        Character convolution counts the residual handle relation. Remove C2
        completions separately, including three subgroups for an all-flat patch.
        """
        counts = Counter(self.q[h] for h in values)
        remaining = [n-counts[q] for q, n in enumerate(populations)]
        if min(remaining) < 0:
            return 0
        _, r1, r2 = remaining
        p = self.g.identity
        for h in values:
            p = self.g.mul[p][h]
        character = sum(i == v for i, v in enumerate(self.g.elements[p]))-1
        sign = (-1)**self.q[p]
        weight = 6*3**r1*2**r2*(1+(-1)**r1*sign)+3*int(r1 == 0)*(-1)**r2*character
        if populations[2] == 0:
            weight -= 4*sum(all(h in (self.g.identity, r) for h in values)
                            for r in self.words.reflections)
        if weight < 0:
            raise ValueError('negative conditional completion count')
        return weight

    @lru_cache(maxsize=None)
    def distribution(self, populations, pattern):
        if sum(populations) != self.faces or min(populations) < 0:
            raise ValueError('invalid population dimensions')
        z = field_weight(*populations, 'nonabelian_reflections')
        if not z or any(pattern.count(q) > n for q, n in enumerate(populations)):
            raise ValueError('empty conditioning event')
        ids = self.buckets[pattern]
        weights = [self.orbit_sizes[i]*self.completion(self.reps[i], populations) for i in ids]
        if sum(weights) != z:
            raise ValueError('complete local marginal does not sum to the fixed-field count')
        return tuple(Fraction(w, z) for w in weights)

    def exact_covariance(self):
        result = np.full((len(self.reps), len(self.reps)), Fraction(), dtype=object)
        falling_faces = math.prod(self.faces-i for i in range(3))
        for row in self.counts['population_counts']:
            pop = tuple(row['populations'])
            population_weight = Fraction(row['canonical_connections'], self.counts['canonical_connections'])
            for pattern, ids in self.buckets.items():
                if any(pattern.count(q) > n for q, n in enumerate(pop)):
                    continue
                assignments = math.prod(math.prod(n-i for i in range(pattern.count(q))) for q, n in enumerate(pop))
                weight = population_weight*Fraction(assignments, falling_faces)
                p = self.distribution(pop, pattern)
                for i, a in enumerate(ids):
                    for j, b in enumerate(ids):
                        result[a, b] += weight*(p[i]*int(i == j)-p[i]*p[j])
        return result

    @lru_cache(maxsize=256)
    def lookup(self, populations):
        vectors = np.zeros((len(self.reps), self.dimension))
        valid = np.zeros(len(self.reps), dtype=bool)
        for pattern, ids in self.buckets.items():
            if any(pattern.count(q) > n for q, n in enumerate(populations)):
                continue
            p = np.array(self.distribution(populations, pattern), dtype=float)
            mean = p@self.whitener[ids]
            for i in ids:
                vectors[i] = self.whitener[i]-mean
                valid[i] = True
        return vectors, valid

    def report(self):
        return {'local_gauge_states': len(self.reps), 'local_charge_patterns': len(self.buckets),
                'conditional_relational_dimension': self.dimension,
                'representatives': self.reps, 'orbit_sizes': [int(self.orbit_sizes[i]) for i in range(len(self.reps))],
                'contrast_columns': self.contrasts.T.tolist(),
                'exact_contrast_covariance': [[fraction(v) for v in row] for row in self.gram],
                'scope': 'Complete fixed functions of a single three-loop disk gauge state, centered on the entire charge field. Not a complete observer of unions of patches or of history.'}


class CompleteObserver:
    def __init__(self, experiment):
        self.x = experiment
        self.reference = CompleteReference(experiment.e.geometry.group, experiment.e.charges,
                                           experiment.geometry.faces, experiment.geometry.faces)
        self.translations = Translations(experiment)
        self.layout = self.translations.layouts[2]
        _, eigenvalues, _ = self.translations.charge_blocks()
        level = min(v for v in eigenvalues[:, :, 0].flat if v > 1e-10)
        self.shell = [list(map(int, a)) for a in np.argwhere(abs(eigenvalues[:, :, 0]-level) < 1e-10)]

    def measure(self, raw):
        q = self.x.p.charge(raw)
        if sum(q) != self.x.geometry.faces:
            raise ValueError('complete-memory observer requires Q=F')
        populations = tuple(q.count(i) for i in range(3))
        oracle = self.x.e.factor.oracle
        codes = [oracle.code(raw, 1, p) for p in range(self.x.supports)]
        states = self.reference.codes[codes]
        lookup, valid = self.reference.lookup(populations)
        if not all(valid[states]):
            raise ValueError('raw patch observed an impossible charge event')
        # [orientation, relational contrast, y, x]; no orientation averaging.
        fields = lookup[states[self.layout]].transpose(0, 3, 1, 2)
        return {'populations': list(populations), 'states': states,
                'fields': fields, 'fourier': np.fft.fft2(fields)}

    def correlate(self, first, later):
        a, b = first['fourier'], later['fourier']
        v = self.x.e.geometry.size
        # Whitening uses SINGLE-PATCH covariance. Mode time-zero values can
        # differ from one due to overlap and spatial correlations.
        by_contrast = np.sum(a.conj()*b, axis=0)/(len(self.layout)*v)
        trace = by_contrast.mean(axis=0)
        shell_values = np.mean([by_contrast[:, y, x].real for y, x in self.shell], axis=0)
        y, x = self.shell[0]
        return {'lowest_shell_by_contrast': shell_values.tolist(),
                'lowest_shell_trace': float(shell_values.mean()),
                'oriented_trace_imaginary': float(trace[y, x].imag),
                'local_trace': float(np.fft.ifft2(trace)[0, 0].real)}, np.fft.ifft2(trace).real


def run_pair(observer, control, sampler, index, times):
    x, s = observer.x, observer.x.e.geometry.side
    seeds = (38100000+s*1000+index, 48200000+s*1000+index)
    raw = sampler.sample(random.Random(seeds[0]), 'nonabelian_reflections')['links']
    ticks = [t*x.geometry.faces for t in times]
    rng = random.Random(seeds[1])
    schedule = [rng.randrange(x.operator_count) for _ in range(ticks[-1])]
    runs, maps = {}, {}
    for name, experiment in (('baseline_padded', control), ('elastic', x)):
        samples, check = replay_samples(experiment, raw, schedule, ticks, observer)
        measured = [observer.correlate(samples[0], sample) for sample in samples]
        runs[name] = {'observations': [a for a, _ in measured],
                      'populations': [a['populations'] for a in samples], **check}
        maps[name] = np.array([b for _, b in measured])
    return {'index': index, 'initial_seed': seeds[0], 'schedule_seed': seeds[1], 'runs': runs}, maps


def audit(side, trials, bank, elastic, times=(0, 1, 2, 4, 8, 16, 32, 64)):
    if type(trials) is not int or trials < 2:
        raise ValueError('at least two independent replicates required')
    x = ElasticExperiment(side, bank, elastic)
    control = ElasticExperiment(side, bank, [list(range(36))]*2)
    observer = CompleteObserver(x)
    sampler = AxialReference(side, bank, x.geometry.faces)
    spatial = {name: Moments() for name in ('baseline_padded', 'elastic', 'paired_elastic_minus_baseline')}
    rows = []
    for i in range(trials):
        row, maps = run_pair(observer, control, sampler, i, times)
        rows.append(row)
        maps['paired_elastic_minus_baseline'] = maps['elastic']-maps['baseline_padded']
        for name, values in maps.items():
            spatial[name].add(values)
    estimates = {}
    for key in rows[0]['runs']['elastic']['observations'][0]:
        arrays = {name: np.array([[a[key] for a in row['runs'][name]['observations']] for row in rows])
                  for name in ('baseline_padded', 'elastic')}
        arrays['paired_elastic_minus_baseline'] = arrays['elastic']-arrays['baseline_padded']
        estimates[key] = {name: {'mean': a.mean(axis=0).tolist(),
                                 'standard_error': (a.std(axis=0, ddof=1)/math.sqrt(trials)).tolist()}
                          for name, a in arrays.items()}
    return {'side': side, 'faces': x.geometry.faces, 'trials': trials,
            'times_in_attempts_per_face': times, 'reference': observer.reference.report(),
            'rotation_gate': rotation_gate(x),
            'shell_wavevectors_yx': observer.shell, 'replicates': rows, 'estimates': estimates,
            'spatial_estimates': {name: m.result() for name, m in spatial.items()},
            'scope': 'Stationary Q=F full-S3 reflection-present reference. Supplied mesh and uniform inverse-paired scheduler. Whitening and Fourier phases are diagnostics, not a Hamiltonian or quantum amplitudes.'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--sides', nargs='+', type=int, default=[6, 12])
    parser.add_argument('--trials', type=int, default=64)
    parser.add_argument('--output', type=Path, default=Path('out/triangle-complete-memory.json'))
    args = parser.parse_args()
    bank = json.loads(Path('data/triangle-feedback.json').read_text())['bank']
    elastic = json.loads(Path('data/triangle-elastic-scattering.json').read_text())['pair_census']['elastic_tables']
    rows = []
    for side in args.sides:
        row = audit(side, args.trials, bank, elastic)
        rows.append(row)
        print('side', side, 'elastic complete shell', row['estimates']['lowest_shell_trace']['elastic'], flush=True)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps({'whole_mesh': rows}, indent=2)+'\n')
    print('Wrote', args.output)


if __name__ == '__main__':
    main()
