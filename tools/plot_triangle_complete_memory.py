#!/usr/bin/env python3
"""Plot the complete local gauge detector, including all charge-block traces."""
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np


data = json.loads(Path('data/triangle-complete-memory.json').read_text())
row = next(r for r in data['whole_mesh'] if r['side'] == 12)
t = np.array(row['times_in_attempts_per_face'])
fig, axes = plt.subplots(1, 3, figsize=(15, 5.5), layout='constrained')
for mode, color, label in (('baseline_padded', '#777777', 'Original, same clock'),
                           ('elastic', '#326aab', 'Elastic bank')):
    r = row['estimates']['lowest_shell_trace'][mode]
    mean, se = np.array(r['mean']), np.array(r['standard_error'])
    axes[0].plot(t, mean, '.-', color=color, label=label)
    axes[0].fill_between(t, mean-se, mean+se, color=color, alpha=.16)
axes[0].set(title='All 22 relational contrasts', ylabel='Mean whitened-contrast Fourier covariance')
q = [0, 1, 1, 2, 2, 1]
patterns = [tuple(q[h] for h in rep) for rep in row['reference']['representatives']]
groups = {}
for j, c in enumerate(row['reference']['contrast_columns']):
    p = patterns[next(i for i, value in enumerate(c) if value)]
    groups.setdefault(p, []).append(j)
values = np.array([[o['lowest_shell_by_contrast'] for o in r['runs']['elastic']['observations']]
                   for r in row['replicates']])
for p, color, label in (((1, 1, 1), '#27856a', '111: four reflection contrasts'),
                       ((2, 2, 2), '#bb6327', '222: three rotation contrasts')):
    a = values[:, :, groups[p]].mean(axis=2)
    mean, se = a.mean(axis=0), a.std(axis=0, ddof=1)/len(a)**.5
    axes[1].plot(t, mean, '.-', color=color, label=label)
    axes[1].fill_between(t, mean-se, mean+se, color=color, alpha=.16)
axes[1].set(title='Previously unobserved rotation relationships', ylabel='Within-block mean Fourier covariance')
all_means = np.array([values[:, :, ids].mean(axis=2).mean(axis=0) for ids in groups.values()])
limit = float(np.max(abs(all_means)))
im = axes[2].imshow(all_means, aspect='auto', cmap='RdBu_r', vmin=-limit, vmax=limit)
axes[2].set_xticks(range(len(t)), [str(v) for v in t])
axes[2].set_yticks(range(len(groups)), [''.join(map(str, p)) for p in groups])
axes[2].set(title='Every nontrivial charge block (no selection)', xlabel='Attempts / faces', ylabel='Ordered local charge pattern')
fig.colorbar(im, ax=axes[2], shrink=.8, label='Mean whitened-contrast Fourier covariance')
for ax in axes[:2]:
    ax.axhline(0, lw=.7, color='#999999')
    ax.set_xscale('symlog', linthresh=1)
    ax.set_xticks([0, 1, 4, 16, 64], ['0', '1', '4', '16', '64'])
    ax.set_xlabel('Attempts / number of faces')
    ax.legend(fontsize=8)
    ax.grid(alpha=.15)
fig.suptitle('Complete local gauge memory: additional transient structure, not a persistent bound object', fontsize=13)
fig.supxlabel('288 faces; Q = F; stationary full-S3 reference; 64 independent paired initializations. Lowest positive acoustic wavevector shell.\n'
              'Whitening uses exact single-patch covariance: Fourier time-zero values need not be one. Shading is ±1 replicate SE, not simultaneous bands.\n'
              'All six support orientations retained. Charge-block breakdown is exploratory; the complete trace was specified before simulation.', fontsize=8.5)
output = Path('docs/images/triangle-complete-memory.png')
fig.savefig(output, dpi=160)
print(output)
