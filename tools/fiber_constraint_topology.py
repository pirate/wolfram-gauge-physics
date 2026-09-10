#!/usr/bin/env python3
"""Exact rational cycle-space ranks for the sign-cover interpretation."""
import json
from pathlib import Path

from flint import fmpz_mat

from fiber_constraint_dynamics import canonical, loop_row
from fiber_refinement_limit import geometry


def analyze(side, signs, ell, offset, Q):
    oracle, args = geometry(side)
    V = max(max(e) for e in oracle.edges)+1
    E, F = len(oracle.edges), len(oracle.faces)
    incidence = [[0]*E for _ in range(V)]
    for e, ((u, v), s) in enumerate(zip(oracle.edges, signs)):
        incidence[u][e] = -s
        incidence[v][e] = 1
    rows = [loop_row(signs, edges, reverse) for edges, reverse in zip(args[2], args[3])]
    R = sum(s == -1 for s, _ in rows)
    boundaries = [list(row) for s, row in rows if s == 1]
    C, B = fmpz_mat(incidence), fmpz_mat(boundaries)
    assert C*B.transpose() == fmpz_mat(V, len(boundaries))
    assert C*fmpz_mat([[x] for x in ell]) == fmpz_mat(V, 1)
    rankC, rankB = C.rank(), B.rank()
    assert R > 0 and rankC == V and rankB == F-R
    homology_dimension = E-rankC-rankB
    carrier_class_nonzero = fmpz_mat(boundaries+[list(ell)]).rank() > rankB
    identities = sum(s == 1 and offset == 0 and canonical(row)[0] == tuple(ell) for s, row in rows)
    assert 2*F-R-2*identities == Q
    genus = 1+R//2
    assert genus+identities == 1+F-Q//2
    return {'vertices': V, 'edges': E, 'faces': F, 'reflection_faces': R, 'identity_faces': identities,
            'twisted_incidence_rank': rankC, 'even_face_boundary_rank': rankB,
            'anti_invariant_homology_dimension_over_Q': homology_dimension,
            'carrier_homology_class_nonzero': carrier_class_nonzero,
            'branched_cover_cell_counts': [2*V, 2*E, 2*F-R],
            'branched_cover_genus': genus, 'genus_plus_identity_count': genus+identities,
            'conserved_charge': Q}


if __name__ == '__main__':
    data = json.loads(Path('data/fiber-constraint-dynamics.json').read_text())
    out = []
    for run in data['runs']:
        for trial, sample in enumerate(run['examples']):
            for when in ('initial', 'final'):
                row = analyze(run['side'], sample[when+'_signs'], sample[when+'_covector'],
                              0 if when == 'initial' else sample['final_offset'], sample['conserved_charge'])
                out.append({'side': run['side'], 'alpha': run['alpha'], 'example': trial,
                            'when': when, 'initial_direction': sample['initial_direction'], **row})
    Path('data/fiber-constraint-topology.json').write_text(json.dumps({'states': out}, indent=2)+'\n')
    for row in out:
        print(json.dumps(row), flush=True)
