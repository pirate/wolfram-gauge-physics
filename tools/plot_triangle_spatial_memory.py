#!/usr/bin/env python3
"""Plot saved spatial measurements without fitting or positivity clipping."""
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np


data = json.loads(Path('data/triangle-spatial-memory.json').read_text())
fig, axes = plt.subplots(1, 3, figsize=(15, 5.4), layout='constrained')
colors = ['#326aab', '#bb6327', '#27856a']
for row, color in zip(data['whole_mesh'], colors):
    t = np.array(row['times_in_attempts_per_face'])
    mean = np.array(row['shell_estimates']['elastic']['mean'])[:, 0]
    se = np.array(row['shell_estimates']['elastic']['standard_error'])[:, 0]
    axes[0].plot(t, mean, '.-', color=color, label=f"{row['faces']} faces: measured")
    axes[0].fill_between(t, mean-se, mean+se, color=color, alpha=.15)
    axes[0].plot(t, np.array(row['charge_jensen_lower_bounds'])[:, 0], '--', color=color,
                 label=f"{row['faces']} faces: exact lower bound")
axes[0].set(title='Slowest charge wavelength', ylabel='Static-variance-normalized covariance')
row = data['whole_mesh'][-1]
t = np.array(row['times_in_attempts_per_face'])
for index, label, color in zip((0, 3, 5), ('Charge', 'Pair alignment residual', 'Fan activity residual'), colors):
    mean = np.array(row['shell_estimates']['elastic']['mean'])[:, index]
    se = np.array(row['shell_estimates']['elastic']['standard_error'])[:, index]
    axes[1].plot(t, mean, '.-', color=color, label=label)
    axes[1].fill_between(t, mean-se, mean+se, color=color, alpha=.15)
axes[1].set(title='Same wavelength, different information', ylabel='Covariance (normalizations described below)')
for ax in axes[:2]:
    ax.axhline(0, color='#999999', lw=.7)
    ax.set_xscale('symlog', linthresh=1)
    ax.set_xticks([0, 1, 4, 16, 64, 128], ['0', '1', '4', '16', '64', '128'])
    ax.set_xlabel('Attempts / number of faces')
    ax.grid(alpha=.16)
    ax.legend(fontsize=8)
index = list(t).index(8)
values = np.fft.fftshift(np.array(row['spatial_correlations']['elastic']['mean'])[index, 0])
limit = max(abs(values.min()), abs(values.max()))
s = row['side']; extent = [-s//2-.5, s//2-.5, -s//2-.5, s//2-.5]
im = axes[2].imshow(values, origin='lower', extent=extent, cmap='RdBu_r', vmin=-limit, vmax=limit)
axes[2].set(title='Actual charge memory across displacement', xlabel='Periodic translation dx', ylabel='Periodic translation dy')
fig.colorbar(im, ax=axes[2], shrink=.7, label='Charge covariance / single-face variance')
fig.suptitle('Translation-resolved memory: slow charge relaxation, no resolved extra persistent gauge mode', fontsize=13)
fig.supxlabel('Elastic bank; Q = F; exact stationary full-S3 reference; 64 independent trajectories per mesh. Shading: ±1 replicate SE, not simultaneous bands.\n'
              'Middle/right: 288 faces. Right: 8 attempts per face, unsmoothed periodic displacement coordinates, not emergent space.\n'
              'Gauge curves use single-support variance, not mode variance: their time-zero value is not constrained to one. Dashed curves are bounds, not fits.', fontsize=8.5)
output = Path('docs/images/triangle-spatial-memory.png')
fig.savefig(output, dpi=160)
print(output)
