#!/usr/bin/env python3
"""A finite integer witness for two above-band modes on the infinite supplied lattice."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.linalg import eigh
from scipy.sparse import coo_matrix

from face_energy_obstruction import derive_fiber_group
from fiber_mode_probe import representation
from mode_localization import exact_positive_witness, verify_positive_witness


def restriction(radius=6):
    """Principal restriction of the infinite degree-six operator, NOT induced degrees."""
    if not 2 <= radius <= 12:
        raise ValueError('invalid bounded compact support')
    coordinates = [(x, y) for y in range(-radius, radius+1) for x in range(-radius, radius+1)]
    ids = {p: i for i, p in enumerate(coordinates)}
    matrices = np.asarray(representation(derive_fiber_group())['matrices'], dtype=np.int64)
    result = 6*np.eye(2*len(ids), dtype=np.int64)
    defects = {((0, 0), (1, 0)): 1, ((1, 0), (2, 0)): 2}
    for (x, y), i in ids.items():
        for dx, dy in ((1, 0), (0, 1), (1, 1)):
            neighbor = x+dx, y+dy
            if neighbor not in ids:
                continue
            j = ids[neighbor]
            value = matrices[defects.get(((x, y), neighbor), 0)]
            result[2*j:2*j+2, 2*i:2*i+2] = -value
            result[2*i:2*i+2, 2*j:2*j+2] = -value.T
    exact_ranks = []
    for g in defects.values():
        a, b = (matrices[g]-matrices[0]).tolist()
        exact_ranks.append(2 if a[0]*b[1]-a[1]*b[0] else int(any(a+b)))
    return coordinates, coo_matrix(result, dtype=float).tocsr(), sum(exact_ranks)


def positive_two_by_two(matrix):
    a, b = matrix
    return a[1] == b[0] and a[0] > 0 and a[0]*b[1]-a[1]*b[0] > 0


def verify(record):
    coordinates, operator, upper = restriction(record['support_radius'])
    q = np.asarray(record['integer_trial_columns'], dtype=object)
    if q.shape != (2*len(coordinates), 2) or any(type(x) is not int for x in q.flat):
        raise ValueError('invalid integer trial columns')
    b = operator.toarray().astype(np.int64).astype(object)-9*np.eye(len(q), dtype=object)
    gram, mass = q.T@b@q, q.T@q
    if not positive_two_by_two(gram.tolist()):
        raise ValueError('compact above-band witness is not positive')
    denominator = record['strict_gap_lower_bound']['denominator']
    numerator = record['strict_gap_lower_bound']['numerator']
    if numerator != 1 or denominator != 250 or not positive_two_by_two((denominator*gram-numerator*mass).tolist()):
        raise ValueError('compact rational gap bound failed')
    for key, matrix in (('integer_rayleigh_matrix', gram), ('integer_mass_matrix', mass)):
        if record[key] != [[str(int(x)) for x in row] for row in matrix]:
            raise ValueError('compact saved matrix differs: '+key)
    if upper != 2 or record['above_band_rank_upper_bound'] != upper:
        raise ValueError('compact rank upper bound failed')
    # Integer columns are far below 2**53, so this conversion is exact.
    if np.max(np.abs(q)) >= 2**53:
        raise ValueError('trial integers exceed exact float64 representation')
    certificate = exact_positive_witness(operator, np.asarray(q, dtype=float))
    if certificate['witness_sha256'] != record['exact_positive_certificate']['witness_sha256']:
        raise ValueError('integer columns and dyadic witness differ')
    verify_positive_witness(operator, record['exact_positive_certificate'])
    return {'above_band_eigenvalue_count': upper, 'both_strictly_above': '9 + 1/250',
            'support_vertices': len(coordinates)}


def generate():
    coordinates, operator, upper = restriction()
    values, vectors = eigh(operator.toarray())
    # Numerical search proposes short integer trial columns; only exact
    # Rayleigh inequalities, not the numerical eigenvalues, prove existence.
    q = np.rint(4096*vectors[:, -2:]).astype(np.int64)
    integers = q.astype(object)
    gram = integers.T@(operator.toarray().astype(np.int64).astype(object)-9*np.eye(len(q), dtype=object))@integers
    mass = integers.T@integers
    result = {'schema': 1, 'support_radius': 6, 'support_vertices': len(coordinates),
              'coordinate_order': 'y outer then x, each from -radius through +radius; two components per vertex',
              'lattice': 'Z^2 with undirected steps (1,0), (0,1), (1,1); degree six retained at support boundary',
              'defects': [{'edge': [[0, 0], [1, 0]], 'group_index': 1},
                          {'edge': [[1, 0], [2, 0]], 'group_index': 2}],
              'integer_trial_columns': q.tolist(), 'above_band_rank_upper_bound': upper,
              'integer_rayleigh_matrix': [[str(int(x)) for x in row] for row in gram],
              'integer_mass_matrix': [[str(int(x)) for x in row] for row in mass],
              'strict_gap_lower_bound': {'numerator': 1, 'denominator': 250},
              'numerical_restricted_top_eigenvalues_not_used_as_proof': values[-2:].tolist(),
              'exact_positive_certificate': exact_positive_witness(operator, q.astype(float)),
              'scope': 'Two discrete square-summable spectral modes above the flat band on a supplied infinite triangular lattice with two specified fiber-link defects. Not emergent matter, geometry, or a physical evolution law.'}
    result['verification'] = verify(result)
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, default=Path('out/compact-mode-certificate.json'))
    parser.add_argument('--verify', type=Path)
    args = parser.parse_args()
    if args.verify:
        print(verify(json.loads(args.verify.read_text())))
        return
    result = generate()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, separators=(',', ':'))+'\n')
    print(result['verification'])


if __name__ == '__main__':
    main()
