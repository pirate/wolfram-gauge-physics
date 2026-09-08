#!/usr/bin/env python3
"""Exact patch memory and measured whole-mesh covariance, from saved evidence."""
import json
from fractions import Fraction
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np


data = json.loads(Path('data/triangle-relational-memory.json').read_text())
fig, axes = plt.subplots(1, 3, figsize=(15, 5.5), layout='constrained')
colors = ['#326aab', '#bb6327', '#27856a']
patch = data['patch_calibration']
steps = np.arange(41)
for mode, label, title, color, style in (
    ('baseline_padded', 'relational_residual', 'Original: frozen-state plateau', '#777777', '--'),
    ('elastic', 'relational_residual', 'Elastic: relational residual', colors[2], '-'),
    ('elastic', 'charge', 'Charge: identical in both banks', colors[0], '-')):
    rows = patch['spectral_measures'][mode][label]
    variance = sum(Fraction(*r['covariance_weight']) for r in rows)
    values = [float(sum(Fraction(*r['covariance_weight'])*Fraction(*r['pole'])**int(t) for r in rows)/variance) for t in steps]
    axes[0].plot(steps, values, style, color=color, label=title)
axes[0].set(title='Exact patch: a false persistence signal', xlabel='Attempts in the isolated 24-operator kernel', ylabel='Normalized autocovariance')
axes[0].legend(fontsize=8.5)

mesh = next(row for row in data['whole_mesh'] if row['side'] == 12)
times = np.array(mesh['times_in_attempts_per_face'])
labels = ['Charge', 'Pair alignment residual', 'Fan activity residual']
for column, label in enumerate(labels):
    mean = np.array(mesh['estimates']['elastic']['mean'])[:, column]
    error = np.array(mesh['estimates']['elastic']['standard_error'])[:, column]
    axes[1].plot(times, mean, '.-', label=label, color=colors[column])
    axes[1].fill_between(times, mean-error, mean+error, color=colors[column], alpha=.15)
axes[1].set(title='Whole mesh: local memory decays', ylabel='Single-support variance-normalized covariance')
axes[1].legend(fontsize=8.5)
for column in (1, 2):
    mean = np.array(mesh['integrated_estimates']['elastic']['mean'])[:, column]
    error = np.array(mesh['integrated_estimates']['elastic']['standard_error'])[:, column]
    axes[2].plot(times, mean, '.-', label=labels[column], color=colors[column])
    axes[2].fill_between(times, mean-error, mean+error, color=colors[column], alpha=.15)
axes[2].set(title='Spatial sums: late-time memory unresolved', ylabel=r'$\langle(\sum h_0)(\sum h_t)\rangle/[N\,\mathrm{Var}(h)]$')
axes[2].legend(fontsize=8.5)
for ax in axes[1:]:
    ax.set_xscale('symlog', linthresh=1)
    ax.set_xticks([0, 1, 4, 16, 64, 128], ['0', '1', '4', '16', '64', '128'])
    ax.set_xlabel('Whole-mesh attempts / number of faces')
for ax in axes:
    ax.axhline(0, color='#999999', linewidth=.7)
    ax.grid(alpha=.16)
fig.suptitle('Gauge memory: distinguish frozen sectors, local relaxation, and spatially summed information', fontsize=14)
fig.supxlabel('Right two panels: 288 faces, Q = F, exact stationary full-S3 reflection-present initializations; 32 independent trajectories.\n'
              'Shading is ±1 replicate standard error, not a simultaneous confidence band. Spatial-sum curves are not normalized to one at time zero. No binding claim.', fontsize=9)
output = Path('docs/images/triangle-relational-memory.png')
fig.savefig(output, dpi=160)
print(output)
