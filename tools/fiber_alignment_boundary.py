#!/usr/bin/env python3
"""Exact conditional alignment gate, not a replacement for full dynamics.

One degree-one angular star evolves; a marked RRR reaction kills at z=1.
The frozen anchors and omitted old transport are explicit restrictions.
"""
import json
from pathlib import Path

import numpy as np
from scipy.linalg import eigh_tridiagonal, solve_banded
from scipy.optimize import brentq


def experiment(n, c, kappa=12., gamma=1., gate_fraction=0.):
    N = n-1
    d = min(N, max(1, int(round(gate_fraction*N+.5))))
    alpha = c*N
    degree = np.full(N, 2.)
    degree[[0, -1]] = 1.
    diagonal = 2*alpha*degree
    diagonal[d-1] += kappa
    off = np.full(N-1, -2*alpha)
    band = np.zeros((3, N))
    band[0, 1:], band[1], band[2, :-1] = off, diagonal, off
    mean_from_state = solve_banded((1, 1), band, np.ones(N))
    j = np.arange(1, N+1, dtype=float)
    exact = N/kappa+np.where(j <= d, d*(d-1)-j*(j-1),
                             (j-d)*(2*N-j-d+1))/(4*alpha)
    relative_error = float(np.max(np.abs(mean_from_state-exact))/np.max(exact))
    eigenvalue = float(eigh_tridiagonal(diagonal, off, select='i', select_range=(0, 0),
                                      eigvals_only=True)[0])
    root = brentq(lambda k: k*(np.tan(k*gate_fraction)+np.tan(k*(1-gate_fraction)))-kappa/(2*c),
                  1e-10, np.pi/(2*max(gate_fraction, 1-gate_fraction))-1e-10)
    source = np.zeros(N)
    source[d-1] = kappa
    band[1] += gamma
    probability_before_loss = solve_banded((1, 1), band, source)
    # The rank-one resolvent formula uses the un-killed angular operator.
    band[1, d-1] -= kappa
    unit = source/kappa
    gdd = solve_banded((1, 1), band, unit)[d-1]
    probability_formula = kappa/(N*gamma*(1+kappa*gdd))
    return {'n': n, 'angular_rate': alpha, 'c': c, 'gate_rate': kappa,
            'gate_index': d, 'limiting_gate_fraction': gate_fraction,
            'scaled_uniform_mean_wait': float(mean_from_state.mean()/N),
            'exact_scaled_uniform_mean_wait': float(exact.mean()/N),
            'continuum_uniform_mean_wait': 1/kappa+(gate_fraction**3+(1-gate_fraction)**3)/(6*c),
            'mean_equation_relative_error': relative_error,
            'scaled_principal_decay_rate': N*eigenvalue,
            'continuum_principal_decay_rate': 2*c*root**2,
            'independent_anchor_loss_rate': gamma,
            'probability_reaction_before_anchor_loss': float(probability_before_loss.mean()),
            'rank_one_probability_formula': float(probability_formula),
            'probability_upper_bound': min(1., kappa/(N*gamma))}


if __name__ == '__main__':
    result = {'scope': 'Conditional killed degree-one angular star with frozen anchors; not full link dynamics',
              'runs': [experiment(n, c, gate_fraction=a)
                       for c, a in ((.1, 0.), (1., 0.), (10., 0.), (1., .25), (1., .5))
                       for n in (9, 33, 129, 513, 2049)]}
    Path('data/fiber-alignment-boundary.json').write_text(json.dumps(result, indent=2)+'\n')
    for r in result['runs']:
        print(json.dumps(r), flush=True)
