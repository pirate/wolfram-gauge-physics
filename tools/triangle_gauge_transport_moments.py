#!/usr/bin/env python3
"""Stationary estimates for a fixed, genuinely gauge-hidden trial corrector.

The training draws choose one coefficient. Independent evaluation draws
estimate its variational functional, not a fitted physical update law.
"""
import json
import random
from fractions import Fraction
from pathlib import Path

import numpy as np

from triangle_charge_response import ChargeMarginal, coefficients
from triangle_gauge_current_corrector import GaugeCurrent
from triangle_reference import AxialReference


def moments(side=6, training=128, evaluation=512, evaluation_seed=136059001):
    g = GaugeCurrent(side)
    bank = json.loads(Path('data/triangle-feedback.json').read_text())['bank']
    f, m = g.x.geometry.faces, g.x.operator_count
    reference = AxialReference(side, bank, f)
    marginal = ChargeMarginal(f, f, 'nonabelian_reflections')
    v = marginal.expectation(g.variance_counts.items())
    w = marginal.expectation(g.twice_energy_counts.items())/2
    a, b, _, _, _ = coefficients(marginal)
    bare = float((a+10*b)/135)
    base = bare-float(2*v*v/(810*f*w))
    c = float(v/w)
    batches = []
    for name, count, seed in [('training', training, 125948001), ('evaluation', evaluation, evaluation_seed)]:
        rows = []
        rng = random.Random(seed)
        for _ in range(count):
            raw = reference.sample(rng, 'nonabelian_reflections')['links']
            h, s = g.elastic_moments(raw)
            op = rng.randrange(m)
            moved = raw[:]
            old, new = g.x.apply(moved, op)
            other = g.hidden_corrector(moved) if old != new else h
            # S is summed over every elastic slot; T uses one independent
            # uniform full-bank slot. Both are unbiased stationary moments.
            t = Fraction(m, 2)*(other-h)**2
            rows.append((s/f, float(t)/f))
        samples = np.array(rows)
        batches.append(samples)
        print(json.dumps({'batch': name, 'draws': count, 'S_per_face': float(samples[:, 0].mean()),
                          'T_per_face': float(samples[:, 1].mean()),
                          'nonzero_T_draws': int(np.count_nonzero(samples[:, 1]))}), flush=True)
    train, observed = batches
    assert train[:, 1].mean() > 0
    d = -c*train[:, 0].mean()/train[:, 1].mean()
    changes = (4*c*d*observed[:, 0]+2*d*d*observed[:, 1])/810
    se = float(changes.std(ddof=1)/np.sqrt(evaluation))
    return {'side': side, 'faces': f, 'training_draws': training, 'evaluation_draws': evaluation,
            'training_seed': 125948001, 'evaluation_seed': evaluation_seed,
            'charge_corrector_coefficient': c, 'fixed_hidden_corrector_coefficient': d,
            'exact_one_function_upper_bound': base,
            'estimated_change_in_upper_bound': float(changes.mean()), 'change_standard_error': se,
            'estimated_two_function_trial_functional': base+float(changes.mean()),
            'relative_extra_reduction_of_bare': -float(changes.mean())/bare,
            'S_per_face_estimate': float(observed[:, 0].mean()),
            'T_per_face_estimate': float(observed[:, 1].mean()),
            'scope': 'Independent stationary evaluation of a fixed variational trial; numerical estimates are not rigorous confidence-certified bounds or a hydrodynamic diffusion coefficient.'}


if __name__ == '__main__':
    print(json.dumps(moments()), flush=True)
