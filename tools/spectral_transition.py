#!/usr/bin/env python3
"""Find and exactly audit a primitive loss of above-band graph modes, not particles."""
import argparse
import hashlib
import json
import random
from collections import Counter
from pathlib import Path

import numpy as np

from bank_equilibrium import BankFactor, compiled_run
from mode_localization import ModeGeometry, full_profile


def integer_inertia(matrix):
    """Exact symmetric congruence elimination with fraction-free Schur complements.

    The working block is previous_pivot times the actual Schur complement.
    Hence the next congruence pivot has sign pivot/previous_pivot. Symmetric
    permutations and integer shears handle zero diagonal pivots without any
    perturbation or floating-point decision. The terminal zero block is nullity.
    """
    a = [[int(x) for x in row] for row in matrix]
    n = len(a)
    if any(len(row) != n for row in a) or any(
            a[i][j] != matrix[i][j] or a[i][j] != a[j][i]
            for i in range(n) for j in range(n)):
        raise ValueError('inertia requires an integer symmetric square matrix')
    previous, counts, signs = 1, {'positive': 0, 'negative': 0, 'zero': 0}, []
    for k in range(n):
        index = next((i for i in range(k, n) if a[i][i]), None)
        if index is None:
            pair = next(((i, j) for i in range(k, n) for j in range(i+1, n) if a[i][j]), None)
            if pair is None:
                counts['zero'] = n-k
                break
            i, j = pair
            # e_i -> e_i + e_j is unimodular and creates diagonal 2*a_ij.
            for t in range(k, n):
                a[i][t] += a[j][t]
            for t in range(k, n):
                a[t][i] += a[t][j]
            index = i
        if index != k:
            a[k], a[index] = a[index], a[k]
            for row in a:
                row[k], row[index] = row[index], row[k]
        pivot = a[k][k]
        sign = 1 if (pivot > 0) == (previous > 0) else -1
        signs.append(sign)
        counts['positive' if sign > 0 else 'negative'] += 1
        for i in range(k+1, n):
            for j in range(i, n):
                numerator = pivot*a[i][j]-a[i][k]*a[k][j]
                value, remainder = divmod(numerator, previous)
                if remainder:
                    raise ArithmeticError('fraction-free congruence division was not exact')
                a[i][j] = a[j][i] = value
        previous = pivot
    return dict(counts, congruence_pivot_signs=signs)


def above_band_inertia(geometry, links):
    matrix = geometry.operator(links).toarray().astype(np.int64)
    matrix -= 9*np.eye(len(matrix), dtype=np.int64)
    return integer_inertia(matrix.tolist())


def triangle_matrix(geometry, links, face):
    """I plus the transported adjacency of a triangle, in its three local frames."""
    vertices = face[:3]
    result = np.eye(6, dtype=np.int64)
    for i, u in enumerate(vertices):
        j = (i+1) % 3
        v = vertices[j]
        x = links[geometry.ids[tuple(sorted((u, v)))]]
        if u > v:
            x = geometry.group.inv[x]
        result[2*j:2*j+2, 2*i:2*i+2] = geometry.matrices[x]
        result[2*i:2*i+2, 2*j:2*j+2] = geometry.matrices[x].T
    return result


def curvature_rank_bound(geometry, links):
    """n_+(L-9I) <= sum_f n_-(S_f), from 9I-L = (1/2) sum_f S_f.

    Counts follow from the holonomy-twisted three-cycle, not a fit to global
    spectra. They are checked by exact inertia on every actual face below.
    """
    predicted = {0: 0, 1: 1, 2: 1, 3: 2, 5: 2}
    counts = []
    for face, sector in zip(geometry.faces, geometry.sectors(links)):
        actual = integer_inertia(triangle_matrix(geometry, links, face).tolist())['negative']
        if actual != predicted[sector]:
            raise ValueError('face inertia differs from derived holonomy-class count')
        counts.append(actual)
    return {'above_band_count_upper_bound': sum(counts), 'face_negative_inertias': counts,
            'class_negative_inertia': predicted}


def verify_characteristic_polynomials(record):
    """Independent exact audit: FLINT integer characteristic polynomials, no LDL.

    For a real-rooted polynomial, Descartes sign variations equal the positive
    root count exactly: positive and negative variations together cannot exceed
    its degree after removing zero roots, and each bounds its respective count.
    Symmetry guarantees real roots here.
    """
    from flint import fmpz_mat
    geometry = ModeGeometry(record['side'])
    results = []
    for endpoint in record['endpoints']:
        matrix = geometry.operator(endpoint['links']).toarray().astype(np.int64).tolist()
        for i in range(len(matrix)):
            matrix[i][i] -= 9
        coefficients = list(fmpz_mat(matrix).charpoly())  # ascending degree
        nullity = next(i for i, x in enumerate(coefficients) if x)
        signs = [1 if x > 0 else -1 for x in reversed(coefficients) if x]
        positive = sum(a != b for a, b in zip(signs, signs[1:]))
        result = {'positive': positive, 'negative': len(matrix)-positive-nullity, 'zero': nullity}
        if any(result[k] != endpoint['exact_inertia_of_L_minus_9I'][k] for k in result):
            raise ValueError('independent characteristic polynomial disagrees with exact inertia')
        results.append(result)
    return results


def find_transition(source):
    study = json.loads(source.read_text())
    experiment = study['evolution']
    # Selection is explicit: first stored run and first sampled disappearance
    # with an earlier positive sample. Bisection finds a crossing, not the first
    # crossing in time: presence need not be monotone between the samples.
    selected = next((run, i) for run in experiment['runs']
                    for i, snapshot in enumerate(run['snapshots'])
                    if i and snapshot['above_flat_band']['rank'] == 0 and
                    run['snapshots'][i-1]['above_flat_band']['rank'] > 0)
    run, sample = selected
    geometry = ModeGeometry(experiment['side'])
    factor = BankFactor(geometry.side, json.loads(Path('data/d4-triple-channels.json').read_text()))
    rng = random.Random(run['seed'])
    schedule = [rng.randrange(49*len(factor.fans)) for _ in range(experiment['attempts'])]
    checked = next(r for r in compiled_run(factor, run['initial_links'], schedule, experiment['stride']) if r['mode'] == run['mode'])
    if checked['events'] != run['events'] or checked['final_links'] != run['snapshots'][-1]['links']:
        raise ValueError('source event history differs from independent compiled replay')
    events = run['events']

    def state_at(count):
        links = run['initial_links'][:]
        for _, rule, patch, code, target in events[:count]:
            if factor.oracle.update(links, rule, patch, factor.tables[rule]) != (code, target):
                raise ValueError('raw transition replay mismatch')
        return links

    lo_tick, hi_tick = [run['snapshots'][i]['tick'] for i in (sample-1, sample)]
    lo, hi = [sum(e[0] <= tick for e in events) for tick in (lo_tick, hi_tick)]
    audit = []
    while hi-lo > 1:
        mid = (lo+hi)//2
        rank = full_profile(geometry, state_at(mid))['above_flat_band']['rank']
        audit.append({'events_applied': mid, 'above_band_rank': rank})
        if rank:
            lo = mid
        else:
            hi = mid
    before, after = state_at(lo), state_at(hi)
    endpoints = []
    for links in (before, after):
        profile = full_profile(geometry, links)
        exact = above_band_inertia(geometry, links)
        if exact['positive'] != profile['above_flat_band']['rank']:
            raise ValueError('numerical threshold and exact above-band count disagree')
        endpoints.append({'links': links, 'face_classes': geometry.sectors(links),
                          'exact_inertia_of_L_minus_9I': exact, 'profile': profile,
                          'curvature_rank_bound': curvature_rank_bound(geometry, links)})
    if not endpoints[0]['exact_inertia_of_L_minus_9I']['positive'] or endpoints[1]['exact_inertia_of_L_minus_9I']['positive']:
        raise ValueError('selected event is not an exact disappearance')
    histograms = [Counter(e['face_classes']) for e in endpoints]
    return {'schema': 1, 'source_sha256': hashlib.sha256(source.read_bytes()).hexdigest(),
            'side': geometry.side, 'trial': run['trial'], 'mode': run['mode'],
            'selection': 'First sampled positive-to-zero bracket in stored run order; bisection yields one crossing, not necessarily the first crossing.',
            'sample_bracket_ticks': [lo_tick, hi_tick], 'bisection': audit,
            'event_index_zero_based': lo, 'event': events[lo],
            'changed_edges': [i for i, (a, b) in enumerate(zip(before, after)) if a != b],
            'class_histogram_preserved': histograms[0] == histograms[1],
            'endpoints': endpoints, 'independent_full_schedule_replay': True,
            'exact_link_inverse': checked['exact_link_inverse'],
            'scope': 'Exact spectral-count change under one existing raw-link rule. Not wave evolution, particle annihilation, or an energy conservation claim.'}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--source', type=Path, default=Path('data/d4-mode-localization.json'))
    parser.add_argument('--output', type=Path, default=Path('out/spectral-transition.json'))
    parser.add_argument('--verify-with-flint', type=Path)
    args = parser.parse_args()
    if args.verify_with_flint:
        print(verify_characteristic_polynomials(json.loads(args.verify_with_flint.read_text())))
        return
    result = find_transition(args.source)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, separators=(',', ':'))+'\n')
    print('Event', result['event'], 'changed edges', result['changed_edges'],
          'exact above-band counts', [e['exact_inertia_of_L_minus_9I']['positive'] for e in result['endpoints']])


if __name__ == '__main__':
    main()
