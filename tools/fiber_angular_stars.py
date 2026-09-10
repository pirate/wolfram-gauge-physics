#!/usr/bin/env python3
"""Exact conditional angular orbits and their Fourier-reduced spectra.

The star coordinates are actual raw-link affine shifts, with link signs
and all inactive links held fixed by the existing angular primitives.
"""
import itertools
import json
import math
from collections import Counter
from fractions import Fraction
from pathlib import Path

import numpy as np
from scipy.linalg import eigvalsh_tridiagonal
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import eigsh

from fiber_angular_transport import AngularFiber
from run_three_face import LinkOracle
from fiber_charge_coupling import lift_response_witness


def affine_coordinates(g, n, raw):
    return [(1 if (g.elements[h][1]-g.elements[h][0]) % n == 1 else -1, g.elements[h][0]) for h in raw]


def describe_stars(f, oracle, raw):
    g, n = f.g, f.n
    affine = affine_coordinates(g, n, raw)
    q = [0 if (h := oracle.transport(raw, path)) == g.identity else 1 if h in f.reflections else 2
         for path in oracle.face_paths]
    incident = [[] for _ in oracle.edges]
    for face, path in enumerate(oracle.face_paths):
        for edge, _ in path:
            incident[edge].append(face)
    active = {e for e, faces in enumerate(incident) if sorted(q[i] for i in faces) == [1, 2]}
    stars = []
    for face, path in enumerate(oracle.face_paths):
        if q[face] != 2:
            continue
        coefficients, constant, suffix_sign = {}, 0, 1
        for edge, reverse in path[::-1]:
            sign, shift = affine[edge]
            coefficient = suffix_sign*(-sign if reverse else 1)
            if edge in active:
                coefficients[edge] = coefficient
            else:
                constant += coefficient*shift
            suffix_sign *= sign
        assert suffix_sign == 1
        z = (constant+sum(coefficient*affine[e][1] for e, coefficient in coefficients.items())) % n
        assert z != 0 and z == g.elements[oracle.transport(raw, path)][0]
        if coefficients:
            degree = len(coefficients)
            stars.append({'rotation_face': face, 'active_edges': list(coefficients),
                          'coefficients': list(coefficients.values()), 'constant': constant % n,
                          'rotation_shift': z, 'degree': degree, 'orbit_size': (n-1)*n**(degree-1)})
    assert sorted(e for star in stars for e in star['active_edges']) == sorted(active)
    return q, affine, active, stars


def raw_factorization(n=5, side=4, seed=721903001):
    f = AngularFiber(n)
    g, oracle = f.g, LinkOracle(side, f.g)
    rng = np.random.default_rng(seed)
    raw = rng.integers(g.n, size=len(oracle.edges)).tolist()
    q, affine, active, stars = describe_stars(f, oracle, raw)
    actual = Counter()
    for p in range(len(oracle.patches)):
        for rule in f.rules:
            target = raw[:]
            oracle.update(target, p, rule)
            if target != raw:
                changed = [e for e in oracle.writes[p] if target[e] != raw[e]]
                assert len(changed) == 1
                actual[(changed[0], target[changed[0]])] += 1
    expected = Counter()
    for star in stars:
        for edge, coefficient in zip(star['active_edges'], star['coefficients']):
            sign, a = affine[edge]
            for step in (-1, 1):
                if (star['rotation_shift']+coefficient*step) % n:
                    target = g.elements.index(tuple((sign*v+a+step) % n for v in range(n)))
                    expected[(edge, target)] = 2
    assert actual == expected
    # Uniform orbit averaging can be sampled exactly, without evolving B.
    resampled = raw[:]
    for star in stars:
        edges, coefficients = star['active_edges'], star['coefficients']
        shift_sum = star['constant']
        for edge, coefficient in zip(edges[:-1], coefficients[:-1]):
            a = int(rng.integers(n))
            sign = affine[edge][0]
            resampled[edge] = g.elements.index(tuple((sign*v+a) % n for v in range(n)))
            shift_sum += coefficient*a
        target_z = int(rng.integers(1, n))
        edge, coefficient = edges[-1], coefficients[-1]
        a = coefficient*(target_z-shift_sum) % n
        sign = affine[edge][0]
        resampled[edge] = g.elements.index(tuple((sign*v+a) % n for v in range(n)))
    new_q, new_affine, new_active, _ = describe_stars(f, oracle, resampled)
    assert new_q == q and new_active == active
    assert all(new_affine[e][0] == affine[e][0] for e in range(len(raw)))
    assert all(resampled[e] == raw[e] for e in range(len(raw)) if e not in active)
    # The factorization must remain true away from the initial point of its orbit.
    resampled_actual = Counter()
    for p in range(len(oracle.patches)):
        for rule in f.rules:
            target = resampled[:]
            oracle.update(target, p, rule)
            if target != resampled:
                changed = [e for e in oracle.writes[p] if target[e] != resampled[e]]
                assert len(changed) == 1
                resampled_actual[(changed[0], target[changed[0]])] += 1
    _, _, _, new_stars = describe_stars(f, oracle, resampled)
    resampled_expected = Counter()
    for star in new_stars:
        for edge, coefficient in zip(star['active_edges'], star['coefficients']):
            sign, a = new_affine[edge]
            for step in (-1, 1):
                if (star['rotation_shift']+coefficient*step) % n:
                    target = g.elements.index(tuple((sign*v+a+step) % n for v in range(n)))
                    resampled_expected[(edge, target)] = 2
    assert resampled_actual == resampled_expected
    return {'cycle_vertices': n, 'side': side, 'faces': len(q), 'seed': seed,
            'initial_links': raw, 'uniform_orbit_sample_links': resampled, 'edges': oracle.edges,
            'charge_field': q, 'stars': stars, 'star_degree_counts': dict(Counter(s['degree'] for s in stars)),
            'angular_component_size': str(math.prod(star['orbit_size'] for star in stars)),
            'nonidle_rooted_channel_slots_initial': sum(actual.values()),
            'nonidle_rooted_channel_slots_resampled': sum(resampled_actual.values()),
            'rate_per_allowed_coordinate_step': 2,
            'scope': 'A conditional raw-link angular orbit, not independent physical faces or a finite-rate replacement rule.'}


def spectral_gap(n, degree):
    length = n-1
    longitudinal = 4*degree*(1-math.cos(math.pi/length))
    if degree == 1:
        return {'n': n, 'degree': degree, 'states': length, 'gap': longitudinal,
                'n_squared_gap': n*n*longitudinal, 'slowest_sector': 'longitudinal'}
    best, best_k = -1., None
    for k in itertools.product(range(n), repeat=degree-1):
        if not any(k):
            continue
        rho = abs(1+sum(np.exp(2j*np.pi*j/n) for j in k))
        if rho > best:
            best, best_k = float(rho), k
    diagonal = np.full(length, 4.*degree)
    diagonal[[0, -1]] = 2.*degree
    transverse = float(eigvalsh_tridiagonal(diagonal, np.full(length-1, -2*best), select='i', select_range=(0, 0))[0])
    gap = min(longitudinal, transverse)
    lower = min(longitudinal, 2*(degree-1)/degree*(1-math.cos(2*math.pi/n)))
    assert gap >= lower-1e-12
    return {'n': n, 'degree': degree, 'states': length*n**(degree-1), 'gap': gap,
            'n_squared_gap': n*n*gap, 'longitudinal_gap': longitudinal, 'lowest_transverse_rate': transverse,
            'maximal_transverse_hopping': best, 'transverse_frequency': best_k,
            'rigorous_gap_lower_bound': lower,
            'slowest_sector': 'longitudinal' if longitudinal <= transverse else 'transverse'}


def full_star_spectrum(n, degree):
    states = [w for w in itertools.product(range(n), repeat=degree) if sum(w) % n]
    ids = {w: i for i, w in enumerate(states)}
    rows, cols, values = [], [], []
    for i, w in enumerate(states):
        count = 0
        for axis in range(degree):
            for step in (-1, 1):
                target = list(w)
                target[axis] = (target[axis]+step) % n
                if sum(target) % n:
                    rows.append(i); cols.append(ids[tuple(target)]); values.append(-2.)
                    count += 2
        rows.append(i); cols.append(i); values.append(float(count))
    L = coo_matrix((values, (rows, cols)), shape=(len(states), len(states))).tocsr()
    measured = np.linalg.eigvalsh(L.toarray()) if len(states) <= 3 else np.sort(eigsh(L, k=3, which='SM', tol=1e-11, return_eigenvectors=False))
    predicted = spectral_gap(n, degree)
    assert abs(measured[0]) < 1e-9 and abs(measured[1]-predicted['gap']) < 1e-9
    return {'n': n, 'degree': degree, 'three_lowest_full_graph_eigenvalues': measured.tolist(),
            'Fourier_predicted_gap': predicted['gap']}


def linear_solution_count(rows, right, variables, prime):
    """Exact affine-system count over F_p, used only for prime cycle sizes."""
    matrix = [[a % prime for a in row]+[b % prime] for row, b in zip(rows, right)]
    rank = 0
    for column in range(variables):
        pivot = next((i for i in range(rank, len(matrix)) if matrix[i][column]), None)
        if pivot is None:
            continue
        matrix[rank], matrix[pivot] = matrix[pivot], matrix[rank]
        inverse = pow(matrix[rank][column], -1, prime)
        matrix[rank] = [(a*inverse) % prime for a in matrix[rank]]
        for i in range(len(matrix)):
            if i != rank and matrix[i][column]:
                factor = matrix[i][column]
                matrix[i] = [(a-factor*b) % prime for a, b in zip(matrix[i], matrix[rank])]
        rank += 1
    if any(not any(row[:-1]) and row[-1] for row in matrix):
        return 0
    return prime**(variables-rank)


def conditional_activity(f, oracle, raw, patch):
    n = f.n
    assert all(n % divisor for divisor in range(2, math.isqrt(n)+1))
    q, affine, active, stars = describe_stars(f, oracle, raw)
    paths = oracle.patches[patch][::-1]
    if not all(oracle.transport(raw, path) in f.reflections for path in paths):
        return Fraction(0)
    directly_read = {edge for path in paths for edge, _ in path if edge in active}
    touched = [star for star in stars if directly_read & set(star['active_edges'])]
    variables = sorted(edge for star in touched for edge in star['active_edges'])
    ids = {edge: i for i, edge in enumerate(variables)}
    forms, constants = [], []
    for path in paths:
        row, constant, suffix = [0]*len(variables), 0, 1
        for edge, reverse in path[::-1]:
            sign, shift = affine[edge]
            coefficient = suffix*(-sign if reverse else 1)
            if edge in ids:
                row[ids[edge]] += coefficient
            else:
                constant += coefficient*shift
            suffix *= sign
        assert suffix == -1
        forms.append(row); constants.append(constant)
    equalities = [[a-b for a, b in zip(forms[i], forms[i+1])] for i in (0, 1)]
    equality_right = [constants[i+1]-constants[i] for i in (0, 1)]
    exclusions, exclusion_right = [], []
    for star in touched:
        row = [0]*len(variables)
        for edge, coefficient in zip(star['active_edges'], star['coefficients']):
            row[ids[edge]] = coefficient
        exclusions.append(row); exclusion_right.append(-star['constant'])
    denominator = math.prod(star['orbit_size'] for star in touched)
    def probability(selected):
        total = 0
        for mask in range(1 << len(touched)):
            zero_stars = [i for i in range(len(touched)) if mask & (1 << i)]
            rows = [equalities[i] for i in selected]+[exclusions[i] for i in zero_stars]
            right = [equality_right[i] for i in selected]+[exclusion_right[i] for i in zero_stars]
            total += (-1)**len(zero_stars)*linear_solution_count(rows, right, len(variables), n)
        return Fraction(total, denominator)
    return probability((0,))+probability((1,))-2*probability((0, 1))


def averaged_reaction_examples(n=5):
    f = AngularFiber(n)
    g, oracle = f.g, LinkOracle(4, f.g)
    r = lambda k: g.elements.index(tuple((k-v) % n for v in range(n)))
    t = lambda k: g.elements.index(tuple((v+k) % n for v in range(n)))
    examples = []
    common_charge = None
    for a, b, expected in ((0, 1, Fraction(0)), (1, 0, Fraction(1, n-1)),
                            (1, 1, Fraction(n-2, n-1)), (0, 0, Fraction(1))):
        w = (r(a), r(b), r(1-a+b), t(1))
        target = (r(a), r(b), r(2-a+b), t(2))
        realization = lift_response_witness(n, {'before_tuple': w, 'after_tuple': target, 'pair_position': 2})
        raw = realization['initial_links']
        q, affine, active, stars = describe_stars(f, oracle, raw)
        assert len(stars) == 1 and stars[0]['degree'] == 1
        common_charge = q if common_charge is None else common_charge
        assert q == common_charge
        root, rim = realization['root'], realization['rim_vertices']
        ordered_loops = [(root, rim[i], rim[i+1], root) for i in reversed(range(4))]
        patch = next(p for p, spec in enumerate(oracle.specs) if spec[::-1] == ordered_loops[:3])
        average = conditional_activity(f, oracle, raw, patch)
        assert average == expected
        total_activity = sum((conditional_activity(f, oracle, raw, p) for p in range(len(oracle.patches))), Fraction(0))
        reverse_encounters = sum(sorted(0 if (h := oracle.transport(raw, path)) == g.identity else 1 if h in f.reflections else 2
                                        for path in paths) == [0, 1, 2] for paths in oracle.patches)
        examples.append({'anchor_labels': [a, b], 'initial_links': raw, 'reaction_patch': patch,
                         'angular_orbit_size': stars[0]['orbit_size'], 'exact_averaged_activity': str(average),
                         'averaged_forward_reaction_rate': str(12*average),
                         'whole_mesh_averaged_activity_sum': str(total_activity),
                         'whole_mesh_averaged_rotation_count_drift': str(12*total_activity-4*reverse_encounters)})
    return {'cycle_vertices': n, 'common_full_charge_field': common_charge, 'examples': examples,
            'scope': 'Different actual angular components with the same whole charge field and boundary; averaging is conditional on the frozen raw data.'}


if __name__ == '__main__':
    report = {'raw_factorizations': [raw_factorization(n, side) for n, side in ((3, 4), (5, 4), (5, 8), (7, 4))],
              'small_full_star_spectra': [full_star_spectrum(n, d) for n in (3, 5, 7) for d in (1, 2, 3)],
              'spectral_gaps': [spectral_gap(n, d) for n in (3, 5, 9, 27, 81, 243) for d in (1, 2, 3)],
              'exact_averaged_reaction_examples': [averaged_reaction_examples(n) for n in (5, 7)]}
    Path('data/fiber-angular-stars.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({'components': [{k: row[k] for k in ('cycle_vertices', 'side', 'star_degree_counts', 'angular_component_size')}
                                     for row in report['raw_factorizations']], 'spectral_gaps': report['spectral_gaps']}, indent=2))
