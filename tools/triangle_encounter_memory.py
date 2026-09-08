#!/usr/bin/env python3
"""Split gauge memory by actual encounters without changing any update rule.

Unconditional contributions, not survival-conditioned covariances: they sum
back to the original complete detector at every displacement and sample time.
"""
import argparse
import json
import math
import random
from fractions import Fraction
from pathlib import Path

import numpy as np

from triangle_complete_memory import CompleteObserver
from triangle_elastic_scattering import ElasticExperiment
from triangle_reference import AxialReference
from triangle_relational_memory import replay_samples
from triangle_spatial_memory import Moments


BINS = ('untouched_raw_patch', 'raw_touched_orbit_unchanged', 'orbit_reconfigured')


class EncounterObserver:
    def __init__(self, experiment, complete):
        self.x, self.complete = experiment, complete
        self.shadow = None
        self.readers = [set() for _ in experiment.e.geometry.edges]
        for patch, edges in enumerate(experiment.e.factor.oracle.fan.reads):
            for edge in edges:
                self.readers[edge].add(patch)
        self.first_raw = np.zeros(experiment.supports, dtype=np.int64)
        self.first_orbit = np.zeros(experiment.supports, dtype=np.int64)
        self.orbit_changes = np.zeros(experiment.supports, dtype=np.int64)
        self.groups = [list(range(complete.reference.dimension)), []]
        for j, column in enumerate(complete.reference.contrasts.T):
            state = next(i for i, v in enumerate(column) if v)
            if complete.reference.patterns[state] == (2, 2, 2):
                self.groups[1].append(j)
        if len(self.groups[1]) != 3:
            raise ValueError('expected three independent rotation-triple contrasts')
        self.rotation_states = complete.reference.buckets[(2, 2, 2)]
        covariance = complete.reference.covariance
        self.rotation_probability = Fraction(4, 3)*sum(covariance[i, i] for i in self.rotation_states)
        for i in self.rotation_states:
            for j in self.rotation_states:
                if covariance[i, j] != self.rotation_probability*(Fraction(int(i == j), 4)-Fraction(1, 16)):
                    raise ValueError('rotation conditional states are not the exact uniform four-state reference')

    def event(self, raw, tick, op, code, target):
        if self.shadow is None or tick <= 0 or code == target:
            raise ValueError('encounter callback requires an initialized changing event')
        x = self.x
        rule, patch = divmod(op, x.supports)
        writes = (x.e.factor.oracle.fan.writes[patch] if rule >= 3 else
                  {x.e.factor.oracle.pairs[patch][0][0][0]})
        changed = {edge for edge in writes if self.shadow[edge] != raw[edge]}
        if not changed:
            raise ValueError('changing event wrote no raw edge')
        affected = set().union(*(self.readers[edge] for edge in changed))
        for p in affected:
            old = int(self.current[p])
            new_code = x.e.factor.oracle.code(raw, 1, p)
            new = int(self.complete.reference.codes[new_code])
            if not self.first_raw[p]:
                self.first_raw[p] = tick
            if self.complete.reference.patterns[old] == (2, 2, 2):
                if self.complete.reference.patterns[new] == (2, 2, 2):
                    raise ValueError('a rotation-only fan was raw-touched without leaving charge 222')
            if old != new:
                if not self.first_orbit[p]:
                    self.first_orbit[p] = tick
                self.orbit_changes[p] += 1
            self.current[p] = new
        for edge in changed:
            self.shadow[edge] = raw[edge]

    def measure(self, raw):
        sample = self.complete.measure(raw)
        if self.shadow is None:
            self.shadow = raw[:]
            self.current = sample['states'].copy()
            self.initial_states = self.current.copy()
            self.initial_rotation = np.array([self.complete.reference.patterns[i] == (2, 2, 2)
                                               for i in self.current])
        if self.shadow != raw or not np.array_equal(self.current, sample['states']):
            raise ValueError('incremental visitation tracker disagrees with the full raw observation')
        if np.any((self.first_orbit > 0) & ((self.first_raw == 0) | (self.first_orbit < self.first_raw))):
            raise ValueError('physical orbit changed before any actual raw write')
        if not np.array_equal(self.first_raw[self.initial_rotation], self.first_orbit[self.initial_rotation]):
            raise ValueError('initial rotation fan survived its first raw encounter in charge 222')
        labels = np.where(self.first_orbit > 0, 2, np.where(self.first_raw > 0, 1, 0))
        fields = sample['fields']
        masks = np.array([(labels == i)[self.complete.layout] for i in range(3)])
        masked_fourier = [np.fft.fft2(fields*mask[:, None]) for mask in masks]
        return {**sample, 'masked_fourier': masked_fourier, 'history_labels': labels,
                'history_bin_counts': np.bincount(labels, minlength=3).tolist(),
                'initial_rotation_bin_counts': np.bincount(labels[self.initial_rotation], minlength=3).tolist()}

    def correlate(self, first, later):
        x, base = self.x, self.complete
        local, shell, zero, imaginary, maps = [], [], [], [], []
        for ids in self.groups:
            a = first['fourier'][:, ids]
            normalization = x.supports*len(ids)
            spectra = np.array([np.sum(a.conj()*b[:, ids], axis=(0, 1))/normalization
                                for b in later['masked_fourier']])
            complete_spectrum = np.sum(a.conj()*later['fourier'][:, ids], axis=(0, 1))/normalization
            np.testing.assert_allclose(spectra.sum(axis=0), complete_spectrum, atol=1e-11)
            spatial = np.fft.ifft2(spectra).real
            maps.append(spatial)
            local.append(spatial[:, 0, 0].tolist())
            shell.append(np.mean([spectra[:, y, x].real for y, x in base.shell], axis=0).tolist())
            zero.append(spectra[:, 0, 0].real.tolist())
            y, kx = base.shell[0]
            imaginary.append(spectra[:, y, kx].imag.tolist())
        n = x.supports
        # Independent integer return-count estimator. Uniformity of the four
        # rotation orbits gives k(i,j)=(4 delta_ij-1)/Pr(charge 222) exactly.
        rotation = np.isin(first['states'], self.rotation_states) & np.isin(later['states'], self.rotation_states)
        same = rotation & (first['states'] == later['states'])
        both_counts = np.bincount(later['history_labels'][rotation], minlength=3)
        same_counts = np.bincount(later['history_labels'][same], minlength=3)
        direct = (4*same_counts-both_counts)/(3*n*float(self.rotation_probability))
        np.testing.assert_allclose(local[1], direct, atol=1e-11)
        return {'local_contributions': local, 'lowest_shell_contributions': shell,
                'zero_wavevector_contributions': zero, 'oriented_imaginary_contributions': imaginary,
                'rotation_both_endpoint_counts': both_counts.tolist(),
                'rotation_same_endpoint_counts': same_counts.tolist(),
                'history_bin_fractions': [v/n for v in later['history_bin_counts']],
                'initial_rotation_bin_fractions_of_all_patches': [v/n for v in later['initial_rotation_bin_counts']]}, np.array(maps)

    def report(self):
        return {'initial_states': self.initial_states.tolist(), 'final_states': self.current.tolist(),
                'first_raw_change_attempt': self.first_raw.tolist(),
                'first_orbit_change_attempt': self.first_orbit.tolist(),
                'orbit_change_counts': self.orbit_changes.tolist()}


def run_pair(complete, control, sampler, index, times):
    x, s = complete.x, complete.x.e.geometry.side
    seeds = (58100000+s*1000+index, 68200000+s*1000+index)
    raw = sampler.sample(random.Random(seeds[0]), 'nonabelian_reflections')['links']
    ticks = [t*x.geometry.faces for t in times]
    rng = random.Random(seeds[1])
    schedule = [rng.randrange(x.operator_count) for _ in range(ticks[-1])]
    runs, maps = {}, {}
    for name, experiment in (('baseline_padded', control), ('elastic', x)):
        observer = EncounterObserver(experiment, complete)
        samples, check = replay_samples(experiment, raw, schedule, ticks, observer, observer.event)
        measured = [observer.correlate(samples[0], sample) for sample in samples]
        runs[name] = {'observations': [a for a, _ in measured],
                      'populations': [a['populations'] for a in samples],
                      'visit_history': observer.report(), **check}
        maps[name] = np.array([b for _, b in measured])
    return {'index': index, 'initial_seed': seeds[0], 'schedule_seed': seeds[1], 'runs': runs}, maps


def audit(side, trials, bank, elastic, times=(0, 1, 2, 4, 8, 16, 32, 64)):
    if type(trials) is not int or trials < 2:
        raise ValueError('at least two independent trajectory pairs are required')
    x = ElasticExperiment(side, bank, elastic)
    control = ElasticExperiment(side, bank, [list(range(36))]*2)
    complete = CompleteObserver(x)
    sampler = AxialReference(side, bank, x.geometry.faces)
    spatial = {name: Moments() for name in ('baseline_padded', 'elastic', 'paired_elastic_minus_baseline')}
    rows = []
    for index in range(trials):
        row, maps = run_pair(complete, control, sampler, index, times)
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
    # Total uncertainties must include covariance among bins: sum WITHIN each
    # replicate before computing an error bar, never add independent-bin errors.
    totals = {}
    for key in ('local_contributions', 'lowest_shell_contributions', 'zero_wavevector_contributions'):
        arrays = {name: np.array([[np.sum(a[key], axis=-1) for a in row['runs'][name]['observations']] for row in rows])
                  for name in ('baseline_padded', 'elastic')}
        arrays['paired_elastic_minus_baseline'] = arrays['elastic']-arrays['baseline_padded']
        totals[key] = {name: {'mean': a.mean(axis=0).tolist(),
                             'standard_error': (a.std(axis=0, ddof=1)/math.sqrt(trials)).tolist()}
                      for name, a in arrays.items()}
    probability = EncounterObserver(x, complete).rotation_probability
    return {'side': side, 'faces': x.geometry.faces, 'trials': trials,
            'times_in_attempts_per_face': times, 'history_bins': BINS,
            'feature_groups': ['complete_22', 'rotation_222_3'],
            'rotation_222_probability': [probability.numerator, probability.denominator],
            'reference': complete.reference.report(), 'shell_wavevectors_yx': complete.shell,
            'replicates': rows, 'estimates': estimates, 'total_estimates': totals,
            'spatial_estimates': {name: m.result() for name, m in spatial.items()},
            'scope': 'Fixed-attempt stationary measurements, Q=F full-S3 reflection-present reference. History bins refer to actual changed links and complete local gauge orbits; no survival-conditioned normalization or changed dynamics.'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--sides', nargs='+', type=int, default=[6, 12])
    parser.add_argument('--trials', type=int, default=64)
    parser.add_argument('--output', type=Path, default=Path('out/triangle-encounter-memory.json'))
    args = parser.parse_args()
    bank = json.loads(Path('data/triangle-feedback.json').read_text())['bank']
    elastic = json.loads(Path('data/triangle-elastic-scattering.json').read_text())['pair_census']['elastic_tables']
    rows = []
    for side in args.sides:
        row = audit(side, args.trials, bank, elastic)
        rows.append(row)
        print('side', side, 'elastic shell', row['estimates']['lowest_shell_contributions']['elastic'], flush=True)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps({'whole_mesh': rows}, indent=2)+'\n')
    print('Wrote', args.output)


if __name__ == '__main__':
    main()
