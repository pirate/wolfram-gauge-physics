#!/usr/bin/env python3
"""Conserved-sector prediction, not a fitted rate or a dynamics change."""
import json
import math
from pathlib import Path

import numpy as np

from fiber_refinement_limit import geometry, kind, read


def falling(k, d):
    return math.prod(k-i for i in range(d))


def analyze(run):
    oracle, args = geometry(run['side'])
    n, F = run['n'], len(oracle.faces)
    rng = np.random.default_rng(run['seed'])
    charges, predictions = [], []
    for i in range(run['trials']):
        # This is the saved independent initial-state generator, not a
        # reconstructed history. Its draws are independent of event clocks.
        raw = rng.integers(2*n, size=len(oracle.edges), dtype=np.int64)
        if i == 0:
            assert raw.tolist() == run['sample_unconditioned_trajectory']['initial_raw_links']
        Q = sum(kind(read(raw, args[2][f], args[3][f], n), n) for f in range(F))
        assert F <= Q <= 2*F and Q % 2 == 0
        K = 2*F-Q
        charges.append(int(Q))
        predictions.append(144*F*run['duration']/n*falling(K, 3)/falling(F, 3))
    counts = np.array(run['reaction_counts_by_trajectory'])
    predicted = np.array(predictions)
    residual = counts-predicted
    weights = np.array([math.comb(F, k)/2**(F-1) if k % 2 == 0 else 0. for k in range(F+1)])
    m = np.array([144*F*run['duration']/n*falling(k, 3)/falling(F, 3) for k in range(F+1)])
    expected = weights @ m
    mixture_variance = weights @ (m-expected)**2
    return {'n': n, 'alpha': run['alpha'], 'faces': F, 'trials': len(counts),
            'conserved_charge_by_trajectory': charges,
            'leading_sector_mean_predictions': predictions,
            'observed_reaction_counts': counts.tolist(),
            'mean_observed_count': float(counts.mean()),
            'mean_predicted_count_for_sampled_sectors': float(predicted.mean()),
            'mean_residual': float(residual.mean()),
            'independent_trajectory_mean_residual_SE': float(residual.std(ddof=1)/np.sqrt(len(counts))),
            'raw_sample_count_variance_over_mean': float(counts.var(ddof=1)/counts.mean()),
            'sector_centered_sample_second_moment_over_mean': float(np.mean(residual**2)/counts.mean()),
            'limiting_sector_mixture_variance': float(mixture_variance),
            'limiting_sector_mixture_variance_over_mean': float(mixture_variance/expected),
            'scope': 'Leading n-to-infinity sector means; residual moment is not an exact finite-n conditional Fano factor'}


if __name__ == '__main__':
    source = json.loads(Path('data/fiber-reaction-bursts.json').read_text())
    result = {'runs': [analyze(run) for run in source['spontaneous_runs']]}
    Path('data/fiber-reaction-sector-rates.json').write_text(json.dumps(result, indent=2)+'\n')
    for r in result['runs']:
        print(json.dumps({k: v for k, v in r.items() if not isinstance(v, list)}), flush=True)
