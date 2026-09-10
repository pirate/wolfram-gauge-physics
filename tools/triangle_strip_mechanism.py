#!/usr/bin/env python3
"""Resolve primitive Dirichlet costs of the slow strip gauge-memory modes."""
import json

import numpy as np
from scipy.sparse.linalg import LinearOperator, eigsh
from scipy.sparse.csgraph import connected_components

from triangle_strip_memory import Strip


def kinetic_components(s):
    """Find all one-rotation states changed only by vacancy moves."""
    hidden = lambda v: v-s.project(v)
    op = LinearOperator((s.n, s.n), dtype=float,
                        matvec=lambda v: hidden(s.laplacian@hidden(v))+2*s.m*s.project(v))
    rates, modes = eigsh(op, k=2, which='SM', tol=1e-11,
                         v0=np.linspace(.3, 1.7, s.n)**2)
    one = np.sum(s.charges == 2, axis=1) == 1
    rules = np.array([r for _, r in s.channels])
    quiet = one & np.all(s.tables[rules != 0] == np.arange(s.n)[None, :], axis=0)
    ids = np.flatnonzero(quiet)
    killed = s.laplacian[ids][:, ids]
    number, labels = connected_components(killed, directed=False)
    results = []
    for component in range(number):
        chosen = ids[labels == component]
        block = s.laplacian[chosen][:, chosen].toarray()
        ordered = sorted(range(len(chosen)), key=lambda i: tuple(s.charges[chosen[i]]))
        results.append({'state_count': len(chosen),
                        'charge_patterns': s.charges[chosen[ordered]].tolist(),
                        'killed_generator': block[np.ix_(ordered, ordered)].tolist(),
                        'minimum_killed_rate': float(np.linalg.eigvalsh(block)[0]),
                        'mean_first_exit_from_each_state_attempts_per_face': (45*np.linalg.solve(block, np.ones(len(chosen))))[ordered].tolist(),
                        'slow_hidden_mode_masses': np.sum(modes[chosen]**2, axis=0).tolist()})
    results.sort(key=lambda row: row['minimum_killed_rate'])
    return {'shape': s.shape, 'faces': s.faces, 'hidden_rates': rates.tolist(),
            'quiet_states': len(ids), 'components': results}


def mechanism(length):
    s = Strip(length)
    hidden = lambda v: v-s.project(v)
    action = lambda v: hidden(s.laplacian@hidden(v))+2*s.m*s.project(v)
    op = LinearOperator((s.n, s.n), matvec=action, dtype=float)
    rates, modes = eigsh(op, k=2, which='SM', tol=1e-11,
                         v0=np.linspace(.3, 1.7, s.n)**2)
    order = np.argsort(rates)
    rates, modes = rates[order], modes[:, order]
    rotations = np.sum(s.charges == 2, axis=1)
    one = rotations == 1
    rpos = np.argmax(s.charges == 2, axis=1)
    vpos = np.argmax(s.charges == 0, axis=1)
    rules = np.array([r for _, r in s.channels])
    escape = np.sum(rotations[s.tables] != rotations[None, :], axis=0)
    holonomies = np.array([[s.x.e.factor.oracle.fan.transport(s.raw[i], path)
                           for path in s.paths] for i in s.representatives])
    parallel = np.array([len(set(h[q == 1])) == 1 for h, q in zip(holonomies, s.charges)])
    endpoint = (rpos == 0) | (rpos == length-1)
    quiet = one & endpoint & parallel & (np.abs(rpos-vpos) >= 3)
    quiet_ids = np.flatnonzero(quiet)
    killed = s.laplacian[quiet_ids][:, quiet_ids].toarray()
    npath = length-3
    assert len(quiet_ids) == 2*npath
    assert np.all(s.tables[rules != 0][:, quiet_ids] == quiet_ids[None, :])
    for side in (0, length-1):
        ordered = sorted((i for i in quiet_ids if rpos[i] == side),
                         key=lambda i: abs(rpos[i]-vpos[i]))
        block = s.laplacian[ordered][:, ordered].toarray()
        expected = np.diag([4]*(npath-1)+[2])-2*np.eye(npath, k=1)-2*np.eye(npath, k=-1)
        assert np.array_equal(block, expected)
        distance = np.arange(1, npath+1)
        passage = distance*(2*npath+1-distance)/4
        assert np.array_equal(block@passage, np.ones(npath))
    killed_rates, killed_modes = np.linalg.eigh(killed)
    trial = np.zeros(s.n)
    trial[quiet_ids] = killed_modes[:, 0]
    trial = hidden(trial)
    trial /= np.linalg.norm(trial)
    # An integer calculation of the conditional-centering penalty. Different
    # quiet configurations have different complete charge fields.
    d = 3**(length-3)
    assert np.all(s.sizes[s.labels[quiet_ids]] == d)
    assert len(set(s.labels[quiet_ids])) == len(quiet_ids)
    centered_integer = -(s.labels[:, None] == s.labels[quiet_ids][None, :]).astype(np.int64)
    centered_integer[quiet_ids, np.arange(len(quiet_ids))] += d
    integer_lap = s.laplacian.astype(np.int64)
    numerator = centered_integer.T@(integer_lap@centered_integer)
    face_ids = {f: i for i, f in enumerate(s.faces)}
    fan_supports = [[face_ids[f] for f in s.x.p.supports[1][p]]
                    for p in sorted({p for p, r in s.channels if r >= 3})]
    reflection_fans = np.array([sum(np.all(s.charges[i, support] == 1) for support in fan_supports)
                                for i in quiet_ids])
    expected_numerator = d*(d-1)*killed.astype(np.int64)+np.diag((16*d//3)*reflection_fans)
    assert np.array_equal(numerator, expected_numerator)
    centered_killed = numerator/(d*(d-1))
    centered_rates = np.linalg.eigvalsh(centered_killed)
    output = []
    for rate, v in zip(rates, modes.T):
        assert np.linalg.norm(action(v)-rate*v) < 1e-8
        gradient_squared = (v[None, :]-v[s.tables])**2/2
        costs = {name: float(gradient_squared[mask].sum())
                 for name, mask in [('vacancy', rules == 0),
                                    ('elastic', (rules == 1) | (rules == 2)),
                                    ('reaction', rules >= 3)]}
        assert abs(sum(costs.values())-rate) < 1e-8
        mass = np.zeros((length, length))
        np.add.at(mass, (rpos[one], vpos[one]), v[one]**2)
        exit_weights = {str(k): float(np.sum(v[one & (escape == k)]**2))
                        for k in sorted(set(escape[one]))}
        output.append({'rate': float(rate), 'dirichlet_costs': costs,
                       'parallel_reflection_one_rotation_mass': float(np.sum(v[one & parallel]**2)),
                       'endpoint_rotation_one_rotation_mass': float(np.sum(v[one & endpoint]**2)),
                       'quiet_endpoint_mass': float(np.sum(v[quiet]**2)),
                       'one_rotation_mass_by_rotation_row_vacancy_column': mass.tolist(),
                       'one_rotation_mass_by_population_change_slot_count': exit_weights,
                       'stationary_one_rotation_counts_by_population_change_slots': {
                           str(k): int(np.sum(one & (escape == k))) for k in sorted(set(escape[one]))}})
    return {'face_count': length, 'modes': output,
            'quiet_endpoint_state_count': len(quiet_ids),
            'quiet_killed_generator': killed.tolist(),
            'quiet_killed_minimum_rate': float(killed_rates[0]),
            'quiet_charge_fiber_size': d,
            'quiet_endpoint_mean_first_exit_attempts_per_face': 45*npath*(npath+1)/4,
            'quiet_all_reflection_fan_counts': reflection_fans.tolist(),
            'exact_centered_quiet_compression_formula': 'K + diag(16*g/(3*(3**(F-3)-1)))',
            'centered_quiet_variational_upper_bound': float(centered_rates[0]),
            'quiet_centered_trial_rayleigh_rate': float(trial@s.laplacian@trial),
            'quiet_centered_trial_slow_pair_overlap': float(np.sum((modes.T@trial)**2))}


if __name__ == '__main__':
    for length in (5, 7):
        print(json.dumps(mechanism(length)), flush=True)
