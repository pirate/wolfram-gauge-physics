#!/usr/bin/env python3
"""Intrinsic charge-mode correlations reconstructed from autonomous raw events.

Laplacian eigenvectors are observables of the supplied mesh, not a dynamical
wave equation. Independent stationary trajectories supply every time series.
"""
import json
import random
import time
from concurrent.futures import ThreadPoolExecutor

import numpy as np

from triangle_charge_response import ChargeMarginal, coefficients
from triangle_reference import AxialReference
from triangle_winding_trajectories import Winding


def measure(side=6, trials=128, duration=4096, seed=218984001, workers=4):
    started = time.monotonic()
    observer = Winding(side)
    x, f = observer.x, observer.faces
    eigenvalues, eigenvectors = np.linalg.eigh(x.geometry.l1)
    bands = []
    for value in eigenvalues[eigenvalues > 1e-8]:
        if not bands or abs(value-bands[-1][0]) > 1e-8:
            if len(bands) == 2:
                break
            bands.append((float(value), np.flatnonzero(abs(eigenvalues-value) < 1e-8)))
    selected = np.concatenate([indices for _, indices in bands])
    vectors = eigenvectors[:, selected]
    slices = []
    start = 0
    for _, indices in bands:
        slices.append(slice(start, start+len(indices)))
        start += len(indices)
    marginal = ChargeMarginal(f, f, 'nonabelian_reflections')
    a, b, chi, _, _ = coefficients(marginal)
    chi = float(chi)
    reference = AxialReference(side, observer.bank, f)
    initials = [reference.sample(random.Random(seed+2*i), 'nonabelian_reflections')['links'] for i in range(trials)]
    lags = np.array([0, 1, 2, 4, 8, 16, 32, 64, 128, 256])
    paths = [np.array(supports) for supports in x.p.supports]

    def run(i):
        rng = random.Random(seed+2*i+1)
        schedule = [rng.randrange(x.operator_count) for _ in range(f*duration)]
        row = x.compiled_output([initials[i]], schedule, f*duration, raw_events=False)['runs'][1]
        events = np.array(row['events'], dtype=np.int64).reshape(-1, 5)
        increments = np.zeros((duration, len(selected)))
        for arity in (0, 1):
            chosen = events[(events[:, 1] >= 3) if arity else (events[:, 1] == 0)]
            ticks, _, patch, old, new = chosen.T
            powers = 6**np.arange(arity+1, -1, -1)
            delta = observer.charges[(new[:, None]//powers) % 6]-observer.charges[(old[:, None]//powers) % 6]
            changes = np.einsum('efm,ef->em', vectors[paths[arity][patch]], delta)
            np.add.at(increments, (ticks-1)//f, changes)
        initial_mode = (np.array(x.p.charge(initials[i]))-1)@vectors
        series = np.vstack([initial_mode, initial_mode+np.cumsum(increments, axis=0)])
        final_mode = (np.array(x.p.charge(row['final_links']))-1)@vectors
        assert np.max(abs(series[-1]-final_mode)) < 1e-9
        length = 1 << (2*len(series)-1).bit_length()
        transform = np.fft.rfft(series, n=length, axis=0)
        autocorrelation = np.fft.irfft(abs(transform)**2, n=length, axis=0)[lags]
        # No empirical time-mean subtraction: the exact ensemble mean is zero.
        return np.array([autocorrelation[:, band].mean(axis=1)/(len(series)-lags)/chi for band in slices])

    samples = []
    with ThreadPoolExecutor(max_workers=workers) as pool:
        for result in pool.map(run, range(trials)):
            samples.append(result)
            if len(samples) % 16 == 0:
                print(json.dumps({'side': side, 'independent_trajectories': len(samples),
                                  'elapsed_seconds': time.monotonic()-started}), flush=True)
    samples = np.array(samples)
    means = samples.mean(axis=0)
    errors = samples.std(axis=0, ddof=1)/np.sqrt(trials)
    rows = []
    for k, (value, indices) in enumerate(bands):
        initial_rate = float(((a+10*b)*value-b*value**2)/(45*chi))
        rows.append({'laplacian_eigenvalue': value, 'degeneracy': len(indices),
                     'normalized_correlation': means[k].tolist(), 'trajectory_cluster_SE': errors[k].tolist(),
                     'exact_initial_rate_per_attempt_per_face': initial_rate,
                     'memoryless_projection_curve': ((1-initial_rate/f)**(f*lags)).tolist()})
    endpoints = (8, 32) if side == 6 else (16, 64) if side == 8 else (32, 128)
    i, j = [lags.tolist().index(t) for t in endpoints]
    left, right = means[0, i], means[0, j]
    assert left > 0 and right > 0
    delta_t = endpoints[1]-endpoints[0]
    rate = -np.log(right/left)/delta_t
    influence = -(samples[:, 0, j]/right-samples[:, 0, i]/left)/delta_t
    rate_se = influence.std(ddof=1)/np.sqrt(trials)
    qk = (2*np.pi/side)**2
    return {'side': side, 'faces': f, 'trials': trials, 'duration_attempts_per_face': duration,
            'seed': seed, 'susceptibility': chi, 'lags_attempts_per_face': lags.tolist(), 'bands': rows,
            'lowest_band_fixed_interval': list(endpoints), 'lowest_band_effective_rate': float(rate),
            'lowest_band_effective_rate_SE': float(rate_se),
            'finite_wavevector_mobility_proxy': float(chi*rate/qk),
            'finite_wavevector_mobility_proxy_SE': float(chi*rate_se/qk),
            'laplacian_calibrated_mobility_proxy': float(chi*rate/(3*bands[0][0])),
            'laplacian_calibrated_mobility_proxy_SE': float(chi*rate_se/(3*bands[0][0])),
            'elapsed_seconds': time.monotonic()-started,
            'scope': 'Charge-mode correlation measurements, fixed-interval effective slopes and conditional diffusive mobility proxies; no hydrodynamic limit, fitted physical update, or quantum amplitudes.'}


if __name__ == '__main__':
    print(json.dumps(measure()), flush=True)
