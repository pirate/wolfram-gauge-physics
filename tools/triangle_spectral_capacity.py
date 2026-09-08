#!/usr/bin/env python3
"""Exact link between derived reaction charge and full-graph spectral capacity."""
import itertools
import random

from fiber_mode_probe import lifted_graph_laplacian
from fiber_transport import TransportGeometry
from spectral_transition import integer_inertia


def face_matrix(group, transports):
    """I plus transported adjacency, in the full three-dimensional fiber space."""
    if len(transports) != 3 or len(group.elements[0]) != 3:
        raise ValueError('triangle capacity needs three links and a three-vertex fiber')
    matrix = [[int(i == j) for j in range(9)] for i in range(9)]
    for i, element in enumerate(transports):
        j = (i+1) % 3
        for a, b in enumerate(group.elements[element]):
            matrix[3*j+b][3*i+a] = matrix[3*i+a][3*j+b] = 1
    return matrix


def capacity_certificate(group, charges):
    if len(charges) != group.n:
        raise ValueError('invalid capacity charge vector')
    rows = []
    for links in itertools.product(range(group.n), repeat=3):
        holonomy = group.identity
        for x in links:
            holonomy = group.mul[x][holonomy]
        inertia = integer_inertia(face_matrix(group, links))
        if inertia['negative'] != charges[holonomy]:
            raise ValueError('derived charge does not equal full-fiber face negative inertia')
        rows.append({'links': links, 'holonomy': holonomy,
                     'inertia': {k: inertia[k] for k in ('positive', 'negative', 'zero')}})
    geometry = TransportGeometry(3, 3, [(0, 1), (1, 2), (2, 0)])
    if geometry.group.elements != group.elements:
        raise ValueError('capacity witness and supplied group enumerations disagree')
    rng = random.Random(59130)
    links = [rng.randrange(group.n) for _ in geometry.edges]
    count = 3*geometry.size
    laplacian = lifted_graph_laplacian(range(geometry.size), geometry.edges, links, group, geometry.fiber_edges)
    lhs = [[24*int(i == j)-2*laplacian[i][j] for j in range(count)] for i in range(count)]
    rhs = [[0]*count for _ in range(count)]
    oriented = [[[0]*count for _ in range(count)] for _ in range(2)]
    face_charges = []
    for f, (face, path) in enumerate(zip(geometry.faces, geometry.paths)):
        local = face_matrix(group, [group.inv[links[e]] if reverse else links[e] for e, reverse in path])
        ids = [3*v+a for v in face[:3] for a in range(3)]
        for i, u in enumerate(ids):
            for j, v in enumerate(ids):
                rhs[u][v] += local[i][j]
                oriented[f % 2][u][v] += local[i][j]
        face_charges.append(integer_inertia(local)['negative'])
    # 3I - L(C_3) = J_3 = 1 1^T is positive semidefinite.
    for v in range(geometry.size):
        for a, b in itertools.product(range(3), repeat=2):
            rhs[3*v+a][3*v+b] += 2
            for color in range(2):
                oriented[color][3*v+a][3*v+b] += 1
    if lhs != rhs:
        raise ValueError('full-bundle triangle decomposition failed')
    if any(2*color[i][j] != lhs[i][j] for color in oriented for i in range(count) for j in range(count)):
        raise ValueError('single-orientation triangle decomposition failed')
    return {'flat_full_bundle_upper_edge': 12,
            'face_connection_census': rows,
            'class_capacity_equals_derived_charge': charges,
            'identity': '2(12 I - L_bundle) = sum_f embedded(S_f) + 2 direct_sum_v J_3',
            'stronger_identity': '12 I - L_bundle = sum_{f of either one orientation} embedded(S_f) + direct_sum_v J_3',
            'assembly_check': {'side': 3, 'links': links, 'twice_12I_minus_L': lhs,
                               'face_negative_inertias': face_charges, 'exact_identity': True,
                               'both_single_orientation_identities': True,
                               'charges_by_face_orientation': [sum(face_charges[c::2]) for c in range(2)]},
            'conclusion': 'On the supplied triangular tori, n_+(L_bundle - 12 I) <= min(Q_up,Q_down) <= floor(Q/2). Only Q=Q_up+Q_down is conserved; neither the orientation-resolved capacity nor the actual mode count is generally conserved.'}
