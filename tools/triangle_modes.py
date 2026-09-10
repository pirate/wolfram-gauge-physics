#!/usr/bin/env python3
"""Full-bundle localized modes, exact compact certificates, and primitive crossings."""
import argparse
import json
from collections import Counter
from pathlib import Path

import numpy as np
from scipy.linalg import eigh
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import eigsh

from fiber_transport import TransportGeometry
from fiber_mode_probe import lifted_graph_laplacian
from spectral_transition import integer_inertia
from triangle_feedback import FeedbackExperiment

PAIRS = ((0, 0), (1, 0), (1, 1), (1, 2))
NAMES = ('flat', 'single_reflection', 'commuting_pair', 'noncommuting_pair')


def group():
    return TransportGeometry(3, 3, [(0, 1), (1, 2), (2, 0)]).group


def operator(coordinates, pair, period=None):
    """Full graph; a finite support retains infinite diagonal degree eight."""
    ids = {p: i for i, p in enumerate(coordinates)}
    if len(ids) != len(coordinates):
        raise ValueError('duplicate base coordinates')
    g = group()
    if len(pair) != 2 or any(type(x) is not int or not 0 <= x < g.n for x in pair):
        raise ValueError('invalid two-link seed')
    center = 0 if period is None else period//2
    defects = {((center, center), (center+1, center)): pair[0],
               ((center+1, center), (center+2, center)): pair[1]}
    if period is not None:
        if type(period) is not int or period < 3:
            raise ValueError('period must be an integer at least three')
        defects = {tuple(tuple(x % period for x in p) for p in edge): value for edge, value in defects.items()}
    count = 3*len(ids)
    rows, cols, values = list(range(count)), list(range(count)), [8]*count
    def edge(u, v):
        rows.extend((u, v)); cols.extend((v, u)); values.extend((-1, -1))
    for (x, y), i in ids.items():
        for a, b in ((0, 1), (1, 2), (2, 0)):
            edge(3*i+a, 3*i+b)
        for dx, dy in ((1, 0), (0, 1), (1, 1)):
            neighbor = (x+dx, y+dy)
            if period is not None:
                neighbor = tuple(a % period for a in neighbor)
            if neighbor not in ids:
                continue
            j = ids[neighbor]
            value = defects.get(((x, y), neighbor), g.identity)
            for a, b in enumerate(g.elements[value]):
                edge(3*i+a, 3*j+b)
    return coo_matrix((values, (rows, cols)), shape=(count, count), dtype=np.int64).tocsr()


def restriction(radius, pair):
    if not 2 <= radius <= 12:
        raise ValueError('compact support radius must be in 2..12')
    coordinates = [(x, y) for y in range(-radius, radius+1) for x in range(-radius, radius+1)]
    return coordinates, operator(coordinates, pair)


def perturbation(pair):
    """The exact finite-rank difference on the three defect endpoints."""
    g = group()
    matrix = np.zeros((9, 9), dtype=np.int64)
    for u, element in enumerate(pair):
        p = np.zeros((3, 3), dtype=np.int64)
        for a, b in enumerate(g.elements[element]):
            p[b, a] = 1
        difference = np.eye(3, dtype=np.int64)-p
        matrix[3*(u+1):3*(u+2), 3*u:3*(u+1)] += difference
        matrix[3*u:3*(u+1), 3*(u+1):3*(u+2)] += difference.T
    return matrix.tolist(), integer_inertia(matrix.tolist())


def positive(matrix):
    if not matrix:
        return True
    if len(matrix) == 1:
        return matrix[0][0] > 0
    if len(matrix) != 2:
        raise ValueError('compact verifier supports at most two trial columns')
    a, b = matrix
    return a[1] == b[0] and a[0] > 0 and a[0]*b[1]-a[1]*b[0] > 0


def exact_forms(laplacian, columns, rank):
    raw = np.asarray(columns, dtype=object)
    count = laplacian.shape[0]
    if raw.shape != (count, rank) or any(type(x) is not int for x in raw.flat):
        raise ValueError('invalid integer compact trial columns')
    maximum = max(map(abs, raw.flat), default=0)
    # Bound every absolute intermediate dot product before using signed int64.
    # The +12 also bounds the shifted operator and dominates the mass form.
    row_bound = int(abs(laplacian).sum(axis=1).max())+12
    if count*maximum*maximum*row_bound >= 2**63:
        raise ValueError('integer trial exceeds the verified no-overflow bound')
    q = np.asarray(raw, dtype=np.int64)
    gram = q.T@(laplacian@q-12*q)
    mass = q.T@q
    return gram.tolist(), mass.tolist()


def verify_compact(record):
    if record['fiber_automorphisms'] != [list(p) for p in group().elements]:
        raise ValueError('compact certificate uses a different derived fiber group')
    if [(case['condition'], tuple(case['seed_pair'])) for case in record['cases']] != list(zip(NAMES, PAIRS)):
        raise ValueError('compact certificate does not cover the four specified controls')
    results = []
    for case in record['cases']:
        coords, laplacian = restriction(record['support_radius'], case['seed_pair'])
        delta, inertia = perturbation(case['seed_pair'])
        # Since L_flat <= 12I, every above-band subspace is positive for Delta,
        # hence has dimension <= n_+(Delta). The compact trial form below gives
        # the matching lower bound by min-max after zero extension to infinity.
        rank = inertia['positive']
        if case['perturbation_matrix'] != delta or case['perturbation_inertia'] != inertia or case['certified_mode_count'] != rank:
            raise ValueError('compact perturbation upper bound differs')
        gram, mass = exact_forms(laplacian, case['integer_trial_columns'], rank)
        if gram != case['integer_rayleigh_matrix'] or mass != case['integer_mass_matrix']:
            raise ValueError('saved compact quadratic forms differ')
        # 20G-M > 0 makes every trial Rayleigh quotient > 12+1/20, not just
        # each column separately. Finite-rank Delta leaves the essential band
        # unchanged, so the certified modes are discrete infinite-graph modes.
        if not positive(gram) or not positive([[20*gram[i][j]-mass[i][j] for j in range(rank)] for i in range(rank)]):
            raise ValueError('compact gap certificate 12 + 1/20 failed')
        if any(sum(case['integer_trial_columns'][3*v+a][j] for a in range(3)) for v in range(len(coords)) for j in range(rank)):
            raise ValueError('trial columns do not have exact zero fiber sum')
        results.append({'condition': case['condition'], 'exact_above_band_count': rank,
                        'strict_lower_edge': '12 + 1/20' if rank else None})
    if 'verification' in record and record['verification'] != results:
        raise ValueError('saved compact verification summary disagrees')
    return results


def generate_compact():
    record = {'support_radius': 6, 'fiber_automorphisms': [list(p) for p in group().elements],
              'coordinate_order': 'y then x, -6..6, followed by fiber vertex 0..2',
              'scope': 'Exact counts above [0,12] for the full infinite supplied bundle graph with specified finite link defects; not physical energies or stable matter.', 'cases': []}
    for name, pair in zip(NAMES, PAIRS):
        coords, laplacian = restriction(record['support_radius'], pair)
        delta, inertia = perturbation(pair)
        rank = inertia['positive']
        q = np.zeros((laplacian.shape[0], rank), dtype=np.int64)
        if rank:
            _, vectors = eigh(laplacian.toarray().astype(float), subset_by_index=[len(q)-rank, len(q)-1])
            for j in range(rank):
                anchor = next(x for x in vectors[:, j] if abs(x) > max(abs(vectors[:, j]))/2)
                if anchor < 0:
                    vectors[:, j] *= -1
            q = np.rint(4096*vectors).astype(np.int64)
            q[2::3] = -q[0::3]-q[1::3]
        gram, mass = exact_forms(laplacian, q.tolist(), rank)
        record['cases'].append({'condition': name, 'seed_pair': list(pair),
                                'perturbation_matrix': delta, 'perturbation_inertia': inertia,
                                'certified_mode_count': rank, 'integer_trial_columns': q.tolist(),
                                'integer_rayleigh_matrix': gram, 'integer_mass_matrix': mass})
    record['verification'] = verify_compact(record)
    return record


def finite_study(certificate, sides):
    verify_compact(certificate)
    radius = certificate['support_radius']
    local, _ = restriction(radius, [1, 2])
    records = []
    for side in sides:
        if not 2*radius+2 <= side <= 192:
            raise ValueError('finite study requires room for the compact witness and side <=192')
        coordinates = [(x, y) for y in range(side) for x in range(side)]
        for case in certificate['cases']:
            laplacian = operator(coordinates, case['seed_pair'], period=side)
            rank = case['certified_mode_count']
            embedded = np.zeros((3*side*side, rank), dtype=np.int64)
            for i, (x, y) in enumerate(local):
                v = ((y+side//2) % side)*side+(x+side//2) % side
                embedded[3*v:3*v+3] = case['integer_trial_columns'][3*i:3*i+3]
            gram, mass = exact_forms(laplacian, embedded.tolist(), rank)
            if gram != case['integer_rayleigh_matrix'] or mass != case['integer_mass_matrix']:
                raise ValueError('torus witness does not reproduce the infinite principal restriction')
            modes = []
            if rank:
                rng = np.random.default_rng(40833)
                values, vectors = eigsh(laplacian.astype(float), k=rank, which='LA', tol=1e-11,
                                       v0=rng.normal(size=laplacian.shape[0]))
                for j, value in enumerate(values):
                    vector = vectors[:, j]
                    residual = float(np.linalg.norm(laplacian@vector-value*vector))
                    if value <= 12.05 or residual > 1e-8:
                        raise ValueError('numerical mode estimate fails its certificate or residual check')
                    density = np.sum(vector.reshape(-1, 3)**2, axis=1)
                    mode = {'eigenvalue_estimate': float(value), 'residual_norm': residual,
                            'base_participation_volume': float(1/np.sum(density**2))}
                    if side == 24:
                        mode['base_projector_density'] = density.tolist()
                    modes.append(mode)
            records.append({'side': side, 'condition': case['condition'], 'full_graph_vertices': 3*side*side,
                            'certified_above_band_count': rank, 'compact_embedding_verified': True,
                            'numerical_modes': modes})
            print('static', side, case['condition'], [m['eigenvalue_estimate'] for m in modes], flush=True)
    return records


def exact_band_count(laplacian):
    from flint import fmpz_mat
    matrix = laplacian.tolist() if isinstance(laplacian, np.ndarray) else [list(row) for row in laplacian]
    if any(len(row) != len(matrix) for row in matrix) or any(type(x) is not int or x != matrix[j][i] for i, row in enumerate(matrix) for j, x in enumerate(row)):
        raise ValueError('exact band count requires an integer symmetric matrix')
    for i in range(len(matrix)):
        matrix[i][i] -= 12
    coefficients = list(fmpz_mat(matrix).charpoly())
    # Symmetry makes all roots real. Descartes' bounds for p(x) and p(-x)
    # therefore saturate after zero roots are removed; sign variation is an
    # exact inertia count here, not a floating eigenvalue threshold.
    signs = [1 if x > 0 else -1 for x in reversed(coefficients) if x]
    positive_count = sum(a != b for a, b in zip(signs, signs[1:]))
    nullity = next(i for i, x in enumerate(coefficients) if x)
    return {'positive': positive_count, 'negative': len(matrix)-positive_count-nullity, 'zero': nullity}


def reaction_audit():
    source = json.loads(Path('data/triangle-feedback.json').read_text())
    e = FeedbackExperiment(6, source['bank'])
    initial = source['evolution']['inputs'][2]
    reactions, trajectory = [], []
    def matrix(links):
        return lifted_graph_laplacian(range(e.geometry.size), e.geometry.edges, links, e.geometry.group, e.geometry.fiber_edges)
    def oriented_charges(links):
        charges = [e.charges[x] for x in e.geometry.holonomies(links)]
        return [sum(charges[c::2]) for c in range(2)]
    for run in source['evolution']['runs']:
        links = initial[:]
        # Revalidate the full compiled prefix history, not only isolated reactions.
        checked = e.compiled_bank([initial], run['schedule'], 1000)['runs'][1]
        if checked != dict(run['checked']['runs'][-1], condition=0):
            raise ValueError('source evolution no longer reproduces')
        previous = exact_band_count(matrix(links))
        if run['seed'] == 0:
            trajectory.append([0, None, previous['positive'], *oriented_charges(links)])
        for tick, rule, patch, code, target in checked['events']:
            if rule:
                before = links[:]
                before_matrix = matrix(before)
                before_count = exact_band_count(before_matrix)
                before_charges = oriented_charges(before)
            actual = e.factor.oracle.update(links, rule, patch, e.tables[rule])
            if actual != (code, target):
                raise ValueError('primitive spectral audit diverges from source events')
            if rule or run['seed'] == 0:
                after_matrix = matrix(links)
                after_count = exact_band_count(after_matrix)
                after_charges = oriented_charges(links)
                if sum(after_charges) != checked['conserved_charge'] or after_count['positive'] > min(after_charges):
                    raise ValueError('observed mode count exceeds the exact orientation-resolved capacity')
                if run['seed'] == 0:
                    trajectory.append([tick, rule, after_count['positive'], *after_charges])
            if rule:
                inertia = []
                for endpoint, counts in ((before_matrix, before_count), (after_matrix, after_count)):
                    shifted = [[x-12*int(i == j) for j, x in enumerate(row)] for i, row in enumerate(endpoint)]
                    exact = integer_inertia(shifted)
                    if any(exact[k] != counts[k] for k in counts):
                        raise ValueError('independent integer inertia and characteristic polynomial disagree')
                    inertia.append(exact)
                replay = e.compiled_bank([before], [rule*len(e.factor.pairs)+patch], 1)['runs'][1]
                if replay['final_links'] != links:
                    raise ValueError('isolated C++ reaction disagrees with the recorded spectral endpoints')
                reactions.append({'seed': run['seed'], 'event': [tick, rule, patch, code, target],
                                  'before_links': before, 'after_links': links[:],
                                  'exact_inertias': inertia, 'independent_characteristic_polynomial_counts': [before_count, after_count],
                                  'orientation_charges_before_after': [before_charges, after_charges],
                                  'decrease_forced_by_capacity': after_count['positive'] < before_count['positive'] and min(after_charges) < before_count['positive'],
                                  'compiled_reaction_check': replay})
        if links != checked['final_links']:
            raise ValueError('spectral replay does not reach the compiled final connection')
    end = source['evolution']['attempts']
    if trajectory[-1][0] != end:
        trajectory.append([end, None, *trajectory[-1][2:]])
    residence = Counter()
    segments = []
    for a, b in zip(trajectory, trajectory[1:]):
        residence[a[2]] += b[0]-a[0]
        if segments and segments[-1]['rank'] == a[2]:
            segments[-1]['end'] = b[0]
        else:
            segments.append({'start': a[0], 'end': b[0], 'rank': a[2]})
    return {'source': 'data/triangle-feedback.json', 'reactions': reactions,
            'seed_zero_full_changing_event_counts': trajectory,
            'seed_zero_attempt_interval_residence': dict(sorted(residence.items())),
            'seed_zero_constant_rank_intervals': segments,
            'scope': 'Exact count at every changing event of seed zero and both endpoints of every reaction in all four runs. Counts are not tracked particle identities; time is supplied attempted-update ticks.'}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--generate-certificate', type=Path)
    parser.add_argument('--certificate', type=Path, default=Path('data/triangle-compact-modes.json'))
    parser.add_argument('--verify-only', action='store_true')
    parser.add_argument('--sides', type=int, nargs='+', default=[16, 24, 48, 96])
    parser.add_argument('--output', type=Path, default=Path('out/triangle-modes.json'))
    args = parser.parse_args()
    if args.generate_certificate:
        record = generate_compact()
        args.generate_certificate.parent.mkdir(parents=True, exist_ok=True)
        args.generate_certificate.write_text(json.dumps(record, separators=(',', ':'))+'\n')
        print(record['verification'])
        return
    record = json.loads(args.certificate.read_text())
    verified = verify_compact(record)
    if args.verify_only:
        print(verified)
        return
    result = {'compact_certificate': str(args.certificate), 'verification': verified,
              'static': finite_study(record, args.sides), 'dynamics': reaction_audit()}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, separators=(',', ':'))+'\n')


if __name__ == '__main__':
    main()
