#!/usr/bin/env python3
"""Exact relational-memory calibration and stationary whole-mesh correlations.

No new update rules: this observes the separately identified elastic bank.
Independent initialized trajectories, not sites or time samples, are replicates.
"""
import argparse
import hashlib
import itertools
import json
import math
import random
from collections import Counter
from fractions import Fraction
from pathlib import Path

import numpy as np
from flint import fmpq, fmpq_mat, fmpz_mat

from triangle_charge_response import ChargeMarginal
from triangle_elastic_scattering import ElasticExperiment
from triangle_patch_observer import PatchUnion
from triangle_reference import AxialReference, fraction, reference


def refine_partition(counts, labels):
    """Coarsest strongly lumpable refinement of an observed partition."""
    levels = [labels[:]]
    while True:
        keys = [(labels[i], tuple(sum(row[j] for j in range(len(labels)) if labels[j] == block)
                                 for block in sorted(set(labels)))) for i, row in enumerate(counts)]
        ids = {key: i for i, key in enumerate(sorted(set(keys)))}
        refined = [ids[key] for key in keys]
        if len(set(refined)) == len(set(labels)):
            return levels
        labels = refined; levels.append(labels)


def spectral_measure(counts, values):
    """Exact autocovariance sum w_lambda lambda^t for a symmetric finite chain."""
    size, clock = len(counts), sum(counts[0])
    if any(sum(row) != clock for row in counts) or any(counts[i][j] != counts[j][i] for i in range(size) for j in range(size)):
        raise ValueError('spectral calibration requires symmetric transition counts')
    matrix = fmpz_mat(counts)
    factors = matrix.charpoly().factor()[1]
    if any(poly.degree() != 1 for poly, _ in factors):
        raise ValueError('this exact rational-pole calibration encountered a nonrational pole')
    poles = sorted((Fraction(-int(poly[0]), int(poly[1])*clock), multiplicity) for poly, multiplicity in factors)
    transition = fmpq_mat(matrix)/clock
    identity = fmpq_mat([[int(i == j) for j in range(size)] for i in range(size)])
    observable = fmpq_mat([[fmpq(str(v))] for v in values])
    measures = []
    for pole, multiplicity in poles:
        projector = identity
        for other, _ in poles:
            if pole != other:
                projector = projector*(transition-identity*fmpq(str(other)))/fmpq(str(pole-other))
        if projector*projector != projector:
            raise ValueError('exact spectral projector is not idempotent')
        weight = Fraction(str((observable.transpose()*projector*observable)[0, 0]))/size
        if weight < 0:
            raise ValueError('symmetric-chain spectral weight is negative')
        measures.append({'pole': fraction(pole), 'multiplicity': multiplicity, 'covariance_weight': fraction(weight)})
    if sum(Fraction(*row['covariance_weight']) for row in measures) != sum(v*v for v in values)/size:
        raise ValueError('spectral measure lost equal-time variance')
    return measures


def patch_calibration(bank, elastic, source):
    experiment = ElasticExperiment(3, bank, elastic)
    union = PatchUnion(experiment.e, (0,))
    kernel = source['patch_activation']
    states = kernel['component_mergers'][0]['combined_component']
    tables = kernel['transition_tables']
    padded = tables[:16]+[list(range(49)) for _ in range(8)]
    count = lambda rows: [[sum(t[i] == j for t in rows) for j in states] for i in states]
    baseline, candidate = count(padded), count(tables)
    charges = [tuple(experiment.p.charge(union.lift(union.reps[i]))[f] for f in experiment.e.fan_faces[0]) for i in states]
    q = [Fraction(values[0]-1) for values in charges]
    active = []
    for state in states:
        links = union.lift(union.reps[state])
        code = experiment.e.factor.oracle.code(links, 1, 0)
        a, b, c = code//36, code//6 % 6, code % 6
        active.append(int(all(experiment.e.charges[x] == 1 for x in (a, b, c)) and ((a == b) != (b == c))))
    # On this ten-state component, two of the four q=111 states are active.
    hidden = [Fraction(2*a-int(qs == (1, 1, 1))) for a, qs in zip(active, charges)]
    ids = {values: i for i, values in enumerate(sorted(set(charges)))}
    refinement = refine_partition(candidate, [ids[q] for q in charges])
    records = {}
    for name, counts in (('baseline_padded', baseline), ('elastic', candidate)):
        records[name] = {label: spectral_measure(counts, values) for label, values in (('charge', q), ('relational_residual', hidden))}
    return {'states': states, 'clock': 24, 'transition_counts': {'baseline_padded': baseline, 'elastic': candidate},
            'charge_values': [list(q) for q in charges], 'relational_residual': [int(x) for x in hidden],
            'charge_partition_refinement': refinement, 'spectral_measures': records,
            'scope': 'Exact isolated three-face component, uniformly weighted; not a whole-mesh relaxation law.'}


def conditional_marginals(g, charges, faces, total_charge):
    """Exact disk-patch marginals inside the full-S3, reflection-present torus reference.

    Fixing a local based tuple of product p leaves R faces and two handles.
    Sum the remaining product convolution against the handle commutator count.
    Local reflection tuples exclude proper subgroups only when all are the same
    reflection: subtract four C2 handle assignments times the binomial face count.
    """
    ref = reference(faces, total_charge, 'nonabelian_reflections')
    marginal = ChargeMarginal(faces, total_charge, 'nonabelian_reflections')
    dp = [[[int(p == g.identity) for p in range(g.n)]]]
    for length in range(1, faces):
        row = [[0]*g.n for _ in range(min(total_charge, 2*length)+1)]
        for q, products in enumerate(dp[-1]):
            for p, count in enumerate(products):
                for h, cost in enumerate(charges):
                    if q+cost < len(row):
                        row[q+cost][g.mul[p][h]] += count
        dp.append(row)
    commutators = Counter(g.mul[g.inv[b]][g.mul[g.inv[a]][g.mul[b][a]]] for a, b in itertools.product(range(g.n), repeat=2))
    reflection = next(h for h, q in enumerate(charges) if q == 1)
    rotation = next(h for h, q in enumerate(charges) if q == 2)

    def completion(k, p):
        remaining = total_charge-k
        if not 0 <= remaining < len(dp[faces-k]):
            return 0
        return sum(count*commutators[g.mul[p][v]] for v, count in enumerate(dp[faces-k][remaining]))

    records = {}
    for k, name in ((2, 'distinct_reflection_pair'), (3, 'active_reflection_fan')):
        remainder = total_charge-k
        proper = 4*math.comb(faces-k, remainder) if total_charge % 2 == 0 and 0 <= remainder <= faces-k else 0
        if k == 2:
            w0, wz = completion(2, g.identity), completion(2, rotation)
            all_weight, active_weight = 3*(w0-proper)+6*wz, 6*wz
        else:
            wr = completion(3, reflection)
            all_weight, active_weight = 27*wr-3*proper, 12*wr
        probability = Fraction(all_weight, ref['canonical_connections'])
        if probability != marginal.probability((0, k, 0)):
            raise ValueError('group-product disk marginal differs from independent charge-field counts')
        if not all_weight:
            raise ValueError('requested reflection conditioning event is empty')
        conditional = Fraction(active_weight, all_weight)
        variance = probability*conditional*(1-conditional)
        records[name] = {'conditioning_probability': fraction(probability),
                         'conditional_mean': fraction(conditional), 'residual_variance': fraction(variance),
                         'proper_subgroup_subtraction_per_equal_tuple': proper}
    return {'faces': faces, 'total_charge': total_charge, 'reference_connections': ref['canonical_connections'],
            'charge_variance': ref['single_face_charge_variance'], 'observables': records}


class RelationalObserver:
    def __init__(self, experiment, marginal):
        self.x = experiment
        self.alpha = Fraction(*marginal['observables']['active_reflection_fan']['conditional_mean'])
        self.beta = Fraction(*marginal['observables']['distinct_reflection_pair']['conditional_mean'])
        self.mean_charge = Fraction(marginal['total_charge'], marginal['faces'])
        self.variances = [Fraction(*marginal['charge_variance'])]
        self.variances += [Fraction(*marginal['observables'][key]['residual_variance'])
                          for key in ('distinct_reflection_pair', 'active_reflection_fan')]
        g, q = experiment.e.geometry.group, experiment.e.charges
        self.pair_values = []
        for a, b in itertools.product(range(g.n), repeat=2):
            rr = q[a] == q[b] == 1
            self.pair_values.append((int(rr and a != b), int(rr)))
        self.fan_values = []
        for a, b, c in itertools.product(range(g.n), repeat=3):
            rrr = q[a] == q[b] == q[c] == 1
            self.fan_values.append((int(rrr and ((a == b) != (b == c))), int(rrr)))

    def measure(self, links):
        x, oracle = self.x, self.x.e.factor.oracle
        return [np.array(x.p.charge(links), dtype=np.int64),
                np.array([self.pair_values[oracle.code(links, 0, p)] for p in range(x.supports)], dtype=np.int64),
                np.array([self.fan_values[oracle.code(links, 1, p)] for p in range(x.supports)], dtype=np.int64)]

    def correlation(self, first, later):
        """Integer sufficient statistics, exact centering; floats only at export."""
        moments = [[int(first[0]@later[0]), int(first[0].sum()), int(later[0].sum())]]
        values = [(Fraction(moments[0][0])-self.mean_charge*(moments[0][1]+moments[0][2]))/len(first[0])+self.mean_charge**2]
        for a, b, conditional in zip(first[1:], later[1:], (self.beta, self.alpha)):
            row = [int(a[:, 0]@b[:, 0]), int(a[:, 0]@b[:, 1]), int(a[:, 1]@b[:, 0]), int(a[:, 1]@b[:, 1])]
            moments.append(row)
            values.append((row[0]-conditional*(row[1]+row[2])+conditional**2*row[3])/len(a))
        return [float(v/variance) for v, variance in zip(values, self.variances)], moments

    def spatial_sums(self, sample):
        result = [Fraction(int(sample[0].sum()))-len(sample[0])*self.mean_charge]
        for values, conditional in zip(sample[1:], (self.beta, self.alpha)):
            result.append(int(values[:, 0].sum())-conditional*int(values[:, 1].sum()))
        return result

    def spatial_means(self, sample):
        return [float(x/len(values)) for x, values in zip(self.spatial_sums(sample), sample)]

    def integrated_correlation(self, first, later):
        """Sum of correlations over ALL support pairs, divided by N Var(h).

        Invariant under spatial permutations of either configuration. This is
        not normalized to one at t=0 and does not identify a bound object.
        """
        return [float(a*b/len(values)/variance) for a, b, values, variance in
                zip(self.spatial_sums(first), self.spatial_sums(later), first, self.variances)]


def replay_samples(experiment, initial, schedule, ticks, observer):
    """All local event targets + final raw links checked against C++; observe fixed attempts.

    C++ performs a full inverse echo. No event-time sampling or trajectory-derived
    centering enters the correlation estimator. Expensive full spectator/frame
    tests are already exhaustive for the primitive laws, not repeated per attempt.
    """
    if not ticks or ticks[0] != 0 or sorted(set(ticks)) != ticks or ticks[-1] != len(schedule):
        raise ValueError('sample ticks must be increasing, start at zero, and end at the schedule length')
    run = experiment.compiled_output([initial], schedule, len(schedule), raw_events=False)['runs'][1]
    if run['condition'] != 0 or run['mode'] != 'combined':
        raise ValueError('sampled replay received the wrong compiled condition')
    raw, event_index, samples = initial[:], 0, [observer.measure(initial)]
    tick_set = set(ticks[1:])
    for tick, op in enumerate(schedule, 1):
        code, target = experiment.apply(raw, op)
        if code != target:
            rule, patch = divmod(op, experiment.supports)
            if event_index >= len(run['events']) or run['events'][event_index] != [tick, rule, patch, code, target]:
                raise ValueError('sampled replay differs from C++ event targets')
            event_index += 1
        if tick in tick_set:
            if sum(experiment.p.charge(raw)) != sum(experiment.p.charge(initial)):
                raise ValueError('sampled replay changed total charge')
            samples.append(observer.measure(raw))
    if event_index != len(run['events']) or raw != run['final_links'] or not run['exact_link_inverse'] or not run['reverse_histogram_echo']:
        raise ValueError('sampled replay failed C++ final-state or inverse verification')
    return samples, {'events': event_index, 'final_links_sha256': hashlib.sha256(json.dumps(raw, separators=(',', ':')).encode()).hexdigest(),
                     'cpp_all_event_targets_final_state_and_inverse_verified': True}


def whole_mesh(side, bank, elastic, trials, times):
    x = ElasticExperiment(side, bank, elastic)
    control = ElasticExperiment(side, bank, [list(range(36)), list(range(36))])
    faces = x.geometry.faces
    marginal = conditional_marginals(x.e.geometry.group, x.e.charges, faces, faces)
    observer = RelationalObserver(x, marginal)
    sampler = AxialReference(side, bank, faces)
    ticks = [t*faces for t in times]
    rows = []
    for trial in range(trials):
        initial_seed, schedule_seed = 7701000+side*1000+trial, 9101000+side*1000+trial
        initial = sampler.sample(random.Random(initial_seed), 'nonabelian_reflections')['links']
        rng = random.Random(schedule_seed)
        schedule = [rng.randrange(x.operator_count) for _ in range(ticks[-1])]
        runs = {}
        for name, experiment in (('baseline_padded', control), ('elastic', x)):
            samples, check = replay_samples(experiment, initial, schedule, ticks, observer)
            measured = [observer.correlation(samples[0], snapshot) for snapshot in samples]
            runs[name] = {'normalized_autocovariances': [r[0] for r in measured],
                          'integer_moment_counts': [r[1] for r in measured],
                          'integrated_correlations': [observer.integrated_correlation(samples[0], sample) for sample in samples],
                          'spatial_indicator_totals': [[[int(v.sum())] if v.ndim == 1 else [int(v[:, 0].sum()), int(v[:, 1].sum())]
                                                        for v in sample] for sample in samples], **check}
        rows.append({'trial': trial, 'initial_seed': initial_seed, 'schedule_seed': schedule_seed,
                     'initial_spatial_means': observer.spatial_means(samples[0]), 'runs': runs})
    arrays = {name: np.array([row['runs'][name]['normalized_autocovariances'] for row in rows]) for name in ('baseline_padded', 'elastic')}
    arrays['paired_elastic_minus_baseline'] = arrays['elastic']-arrays['baseline_padded']
    estimates = {name: {'mean': a.mean(axis=0).tolist(), 'standard_error': (a.std(axis=0, ddof=1)/math.sqrt(trials)).tolist()}
                 for name, a in arrays.items()}
    integrated = {name: np.array([row['runs'][name]['integrated_correlations'] for row in rows]) for name in ('baseline_padded', 'elastic')}
    integrated['paired_elastic_minus_baseline'] = integrated['elastic']-integrated['baseline_padded']
    integrated_estimates = {name: {'mean': a.mean(axis=0).tolist(), 'standard_error': (a.std(axis=0, ddof=1)/math.sqrt(trials)).tolist()}
                            for name, a in integrated.items()}
    return {'side': side, 'faces': faces, 'charge': faces, 'independent_trials': trials,
            'times_in_attempts_per_face': list(times), 'attempt_ticks': ticks,
            'observable_order': ['charge', 'distinct_reflection_pair_residual', 'active_reflection_fan_residual'],
            'exact_marginals': marginal, 'replicates': rows, 'estimates': estimates, 'integrated_estimates': integrated_estimates,
            'scope': 'Exact stationary reference initializations; sites are spatial averages within a replicate, not independent samples. Error bars are pointwise replicate standard errors, not simultaneous intervals or a mixing/binding test.'}


def audit(trials=32, sides=(3, 6, 12), times=(0, 1, 2, 4, 8, 16, 32, 64, 128)):
    if type(trials) is not int or trials < 2:
        raise ValueError('at least two independent trials are required')
    bank = json.loads(Path('data/triangle-feedback.json').read_text())['bank']
    source = json.loads(Path('data/triangle-elastic-scattering.json').read_text())
    elastic = source['pair_census']['elastic_tables']
    return {'scope': 'Relational persistence diagnostics of the existing elastic candidate, not a new force or a quantum model.',
            'patch_calibration': patch_calibration(bank, elastic, source),
            'whole_mesh': [whole_mesh(side, bank, elastic, trials, times) for side in sides]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--trials', type=int, default=32)
    parser.add_argument('--sides', type=int, nargs='+', default=[3, 6, 12])
    parser.add_argument('--output', type=Path, default=Path('out/triangle-relational-memory.json'))
    args = parser.parse_args()
    result = audit(args.trials, args.sides)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2)+'\n')
    for row in result['whole_mesh']:
        estimates = row['estimates']['elastic']
        print('side', row['side'], 'late correlations', estimates['mean'][-1], 'SE', estimates['standard_error'][-1])
    print('Wrote', args.output)


if __name__ == '__main__':
    main()
