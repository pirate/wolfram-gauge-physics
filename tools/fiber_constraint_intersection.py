#!/usr/bin/env python3
"""Intersection form of inherited cycles on the auxiliary sign cover.

The cyclic order comes from the actual oriented mesh, not a layout.
No motion or reaction channel is changed by this calculation.
"""
import copy
import json
from collections import Counter
from pathlib import Path

from flint import fmpz_mat

from fiber_constraint_dynamics import canonical, loop_row
from fiber_two_constraints import AffineProcess
from fiber_two_constraint_witness import rates


def intersection_matrix(oracle, signs):
    E = len(oracle.edges)
    edge_id = {edge: i for i, edge in enumerate(oracle.edges)}
    following = {}
    for face in oracle.faces:
        for i in range(3):
            v, a, b = face[i], face[(i+1) % 3], face[(i+2) % 3]
            following.setdefault(v, {})[a] = b
    J = [[0]*E for _ in range(E)]
    for v, successor in following.items():
        start = min(successor)
        order = [start]
        while successor[order[-1]] != start:
            order.append(successor[order[-1]])
        assert len(order) == len(successor)
        darts = []
        for w in order:
            e = edge_id[tuple(sorted((v, w)))]
            darts.append((e, signs[e] if v < w else -1))
        # Half the cyclic-order wedge on each sheet. The two sheets
        # contribute equally for anti-invariant cycles, cancelling 1/2.
        for i, (e, a) in enumerate(darts):
            for f, b in darts[i+1:]:
                J[e][f] += a*b
                J[f][e] -= a*b
    return fmpz_mat(J)


def normals(state):
    B = fmpz_mat([list(h[1][:-1]) for h in state.raw])
    kernel, rank = B.transpose().nullspace()
    return [list(canonical([int(kernel[e, i]) for e in range(state.E)])[0]) for i in range(rank)]


def analyze(state, ambient=False):
    signs = [h[0] for h in state.raw]
    J = intersection_matrix(state.oracle, signs)
    rows = normals(state)
    N = fmpz_mat(rows)
    restricted = N*J*N.transpose()
    out = {'identity_count': state.charges().count(0),
           'constraint_intersection_matrix': [[int(x) for x in row] for row in restricted.tolist()],
           'constraint_intersection_rank': restricted.rank(),
           'constraint_rows': rows,
           'features': state.features()}
    if ambient:
        V = max(max(e) for e in state.oracle.edges)+1
        C = [[0]*state.E for _ in range(V)]
        for e, ((u, v), s) in enumerate(zip(state.oracle.edges, signs)):
            C[u][e] = -s; C[v][e] = 1
        C = fmpz_mat(C)
        cycles, dimension = C.nullspace()
        Z = fmpz_mat([[int(cycles[e, k]) for k in range(dimension)] for e in range(state.E)])
        boundaries = [row for s, row in [loop_row(signs, edges, reverse)
                      for edges, reverse in zip(state.fe, state.fr)] if s == 1]
        assert C*N.transpose() == fmpz_mat(V, state.rank)
        if boundaries:
            assert fmpz_mat(boundaries)*J*Z == fmpz_mat(len(boundaries), dimension)
        induced = Z.transpose()*J*Z
        R = state.charges().count(1)
        assert induced.rank() == R
        out.update(cycle_space_dimension=dimension, even_face_boundaries=len(boundaries),
                   ambient_intersection_rank=induced.rank(), reflection_faces=R)
    return out


def clone(state):
    out = copy.copy(state)
    out.raw = state.raw.copy()
    out.events = state.events.copy()
    return out


def elastic_certificate():
    """Exact whole-cycle matrix identity for every sign choice in one read patch."""
    state = AffineProcess('contact', 4, 914407001)
    E = state.E
    V = max(max(e) for e in state.oracle.edges)+1
    read_edges = sorted(state.oracle.reads[0])
    cases = 0
    for mask in range(1 << len(read_edges)):
        signs = [1]*E
        for i, e in enumerate(read_edges):
            signs[e] = -1 if mask >> i & 1 else 1
        state.raw = [(s, tuple(int(e == j) for j in range(E))+(0,))
                     for e, s in enumerate(signs)]
        C = [[0]*E for _ in range(V)]
        for e, ((u, v), s) in enumerate(zip(state.oracle.edges, signs)):
            C[u][e] = -s; C[v][e] = 1
        kernel, dimension = fmpz_mat(C).nullspace()
        Z = fmpz_mat([[int(kernel[e, k]) for k in range(dimension)] for e in range(E)])
        J = intersection_matrix(state.oracle, signs)
        for channel in (1, 2):
            after = clone(state)
            after.step(0, channel, 0., None, (), 0)
            M = fmpz_mat([list(h[1][:-1]) for h in after.raw])
            assert abs(M.det()) == 1
            Jnew = intersection_matrix(after.oracle, [h[0] for h in after.raw])
            Minv = M.inv()
            defect = Z.transpose()*(Minv*Jnew*Minv.transpose()-J)*Z
            assert all(not x for row in defect.tolist() for x in row)
            cases += 1
    return {'patch': 0, 'read_edges': read_edges, 'sign_assignments': 1 << len(read_edges),
            'channels': [1, 2], 'whole_cycle_matrix_identities': cases,
            'outside_read_signs': 'positive', 'identity': 'Z^T (M^-1 J_new M^-T - J_old) Z = 0'}


def single_rate(state, p):
    w = state.word(p)
    q = [1 if h[0] == -1 else 0 if not any(h[1]) else 2 for h in w]
    if q == [1, 1, 1]:
        return 12 if (w[0][1] == w[1][1]) != (w[1][1] == w[2][1]) else 0
    return 4 if sorted(q) == [0, 1, 2] else 0


def angular_rate_drift(state, initial=None):
    initial = rates(state) if initial is None else initial
    charges = state.charges()
    drift, changes = 0, []
    for p in range(len(state.pe)):
        if tuple(charges[int(f)] for f in state.supports[p][:2]) not in ((1, 2), (2, 1)):
            continue
        for channel in (15, 16):
            after = clone(state)
            after.step(p, channel, 0., None, (), 0)
            if after.raw == state.raw:
                continue
            writes = {e for e in state.oracle.writes[p] if after.raw[e] != state.raw[e]}
            affected = [j for j, edges in enumerate(state.oracle.reads) if writes & edges]
            delta = sum(single_rate(after, j)-initial[j] for j in affected)
            drift += delta
            if delta:
                changes.append({'fan': p, 'channel': channel, 'delta_reaction_rate': delta})
    return drift, changes


def local_coefficients(condition, side, seed):
    """Local generator expansion, retaining all guards and their read supports."""
    state = AffineProcess(condition, side, seed)
    initial = rates(state)
    support = {e for e in range(state.E) if any(row[e] for row in state.initial_rows)}
    G_lambda, G0_A_lambda, reaction_channels, feedback_channels = 0, 0, 0, 0
    for p in range(len(state.pe)):
        # Outside the inherited support, elastic writes leave C and the
        # all-reflection sign field unchanged. Reactions are initially idle.
        if not (state.oracle.writes[p] & support):
            assert initial[p] == 0
            continue
        for channel in range(1, 15):
            if channel >= 3 and initial[p] == 0:
                continue
            after = clone(state)
            after.step(p, channel, 0., None, (), 0)
            if after.raw == state.raw:
                continue
            writes = {e for e in state.oracle.writes[p] if after.raw[e] != state.raw[e]}
            after_rates = initial.copy()
            for j, read_edges in enumerate(state.oracle.reads):
                if writes & read_edges:
                    after_rates[j] = single_rate(after, j)
            G_lambda += sum(after_rates)-sum(initial)
            if channel >= 3:
                reaction_channels += 1
                drift, changes = angular_rate_drift(after, after_rates)
                G0_A_lambda += drift
                feedback_channels += bool(changes)
    assert reaction_channels == sum(initial)
    return {'condition': condition, 'side': side, 'faces': state.F, 'seed': seed,
            'initial_pair_faces': state.initial_pair_faces,
            'lambda_0': sum(initial), 'G_lambda_0': G_lambda,
            'G0_A_lambda_0': G0_A_lambda, 'feedback_reaction_channels': feedback_channels}


def calculation():
    initial, transitions, coefficients = [], [], []
    for condition in ('contact', 'near_disjoint', 'separated'):
        state = AffineProcess(condition, 4, 914407001)
        before = analyze(state, True)
        initial.append({'condition': condition, 'initial_pair_faces': state.initial_pair_faces, **before})
        G_lambda, G0_A_lambda, feedback = 0, 0, []
        for p in range(len(state.pe)):
            for channel in range(17):
                after = clone(state)
                after.step(p, channel, 0., None, (), 0)
                if after.raw == state.raw:
                    continue
                row = analyze(after, True)
                G_lambda += row['features'][0]-before['features'][0]
                if 3 <= channel < 15:
                    drift, changes = angular_rate_drift(after)
                    G0_A_lambda += drift
                    if drift:
                        feedback.append({'reaction_fan': p, 'reaction_channel': channel,
                                         'A_lambda_after_reaction': drift, 'angular_changes': changes})
                transitions.append({'condition': condition, 'fan': p, 'channel': channel,
                                    'before_intersection_rank': before['constraint_intersection_rank'],
                                    **row})
        coefficients.append({'condition': condition, 'lambda_0': before['features'][0],
                             'G_lambda_0': G_lambda, 'G0_A_lambda_0': G0_A_lambda,
                             'angular_feedback_witnesses': feedback})
    return {'scope': 'Exact auxiliary surface intersection; not a physical spatial metric',
            'initial': initial, 'one_step_transitions': transitions,
            'short_time_reaction_count_coefficients': coefficients,
            'elastic_matrix_certificate': elastic_certificate()}


if __name__ == '__main__':
    import sys
    if '--volumes' in sys.argv:
        result = {'scope': 'Exact local reaction-count coefficients, not fitted trajectories', 'runs': []}
        for side in (4, 8, 16):
            for condition in ('contact', 'near_disjoint', 'separated'):
                row = local_coefficients(condition, side, 914407001)
                result['runs'].append(row)
                Path('data/fiber-constraint-proximity.json').write_text(json.dumps(result, indent=2)+'\n')
                print(json.dumps(row), flush=True)
        raise SystemExit(0)
    result = calculation()
    Path('data/fiber-constraint-intersection.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({'initial': [{k: v for k, v in row.items() if k != 'constraint_rows'} for row in result['initial']],
                      'transitions': len(result['one_step_transitions']),
                      'transition_counts': [list(key)+[count] for key, count in sorted(Counter(
                          (row['condition'], 'reaction' if 3 <= row['channel'] < 15 else 'elastic',
                           row['before_intersection_rank'], row['constraint_intersection_rank'])
                          for row in result['one_step_transitions']).items())]}), flush=True)
