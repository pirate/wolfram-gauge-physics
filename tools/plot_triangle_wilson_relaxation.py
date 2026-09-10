#!/usr/bin/env python3
"""Nonstationary loop, gauge-residual, and reaction transients from saved runs."""
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np


data = json.loads(Path('data/triangle-wilson-quench.json').read_text())
fig, axes = plt.subplots(1, 3, figsize=(14, 4.6))
colors = ['#275b9b', '#cc7128', '#428568']
labels = ['Clustered S₃', 'Dispersed S₃', 'Clustered C₂']
for e, color, label in zip(data['experiments'], colors, labels):
    t = np.array(e['times_attempts_per_face'])
    o = e['observables']
    means = [np.array(o['wilson']['mean'])[:, 4, 1],
             np.array(o['gauge_residual']['mean'])[:, 4], np.array(o['rotation_population']['mean'])]
    errors = [np.array(o['wilson']['trajectory_standard_error'])[:, 4, 1],
              np.array(o['gauge_residual']['trajectory_standard_error'])[:, 4],
              np.array(o['rotation_population']['trajectory_standard_error'])]
    for ax, mean, error in zip(axes, means, errors):
        ax.errorbar(t, mean, yerr=2*error, color=color, marker='o', ms=4, lw=1.2, capsize=2,
                    label=label)

axes[0].axhline(data['S3_stationary_loop_means'][4][1], color='#666666', ls='--', lw=1,
                label='Sector reference means')
axes[0].axhline(data['C2_stationary_loop_means'][4][1], color='#666666', ls='--', lw=1)
axes[0].set(title='Area-eight standard Wilson loop', ylabel='Translation-averaged character')
axes[1].axhline(0, color='#777777', ls='--', lw=1)
axes[1].set(title='After removing the charge prediction', ylabel='Charge-orthogonal loop residual')
num, den = data['S3_stationary_mean_rotations']
axes[2].axhline(num/den, color='#666666', ls='--', lw=1)
axes[2].set(title='Reaction overshoot versus delayed encounters', ylabel='Number of rotation faces')
for ax in axes:
    ax.set_xscale('symlog', linthresh=1)
    ax.set_xlim(-0.05, 6000)
    ax.set_xticks([0, 1, 16, 256, 4096], ['0', '1', '16', '256', '4096'])
    ax.set_xlabel('Attempts per face')
    ax.spines[['top', 'right']].set_visible(False)
    ax.grid(alpha=0.15)
axes[0].legend(fontsize=8, loc='lower right', frameon=False)
fig.suptitle('Gauge statistics develop from elementary nonstationary link seeds', fontsize=15, y=0.99)
fig.text(0.5, 0.90, '72 supplied faces · fixed charge 18 · 128 independent schedules per seed · error bars ±2 SE · lines are guides, not fits',
         ha='center', fontsize=9.5, color='#555555')
fig.subplots_adjust(top=0.78, bottom=0.17, left=0.05, right=0.99, wspace=0.31)
output = Path('docs/images/triangle-wilson-relaxation.png')
fig.savefig(output, dpi=180)
print(output.resolve())
