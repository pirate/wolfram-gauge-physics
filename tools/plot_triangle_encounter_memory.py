#!/usr/bin/env python3
"""Measured encounter decomposition; no survival conditioning or curve fitting."""
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np


data = json.loads(Path('data/triangle-encounter-memory.json').read_text())
row = next(r for r in data['whole_mesh'] if r['side'] == 12)
t = np.array(row['times_in_attempts_per_face'])
fig, axes = plt.subplots(1, 3, figsize=(15, 5.4), layout='constrained')
colors = ['#326aab', '#bb6327', '#27856a']


def line(ax, key, group, category, label, color):
    r = row['estimates'][key]['elastic']
    mean = np.array(r['mean'])[:, group, category]
    se = np.array(r['standard_error'])[:, group, category]
    ax.plot(t, mean, '.-', color=color, label=label)
    ax.fill_between(t, mean-se, mean+se, color=color, alpha=.16)


for ax, key, title in ((axes[0], 'local_contributions', 'Rotation memory at the original patch'),
                       (axes[2], 'lowest_shell_contributions', 'Rotation memory at the longest shell')):
    line(ax, key, 1, 0, 'No actual raw touch', colors[0])
    line(ax, key, 1, 2, 'Local gauge state changed', colors[1])
    total = row['total_estimates'][key]['elastic']
    ax.plot(t, np.array(total['mean'])[:, 1], '--', color='#777777', label='Total, sum of history bins')
    ax.set(title=title, ylabel='Unconditional normalized contribution')
line(axes[1], 'local_contributions', 0, 2, 'All 22 contrasts: after reconfiguration', colors[2])
line(axes[1], 'local_contributions', 1, 2, 'Rotation contrasts: after reconfiguration', colors[1])
axes[1].set(title='Small local return excess remains', ylabel='Unconditional normalized contribution')
for ax in axes:
    ax.axhline(0, color='#999999', lw=.7)
    ax.set_xscale('symlog', linthresh=1)
    ax.set_xticks([0, 1, 4, 16, 64], ['0', '1', '4', '16', '64'])
    ax.set_xlabel('Attempts / number of faces')
    ax.grid(alpha=.15)
    ax.legend(fontsize=8)
fig.suptitle('Encounter-resolved gauge memory: distinguish untouched survival from post-interaction return', fontsize=13)
fig.supxlabel('Elastic bank; 288 faces; Q = F; exact stationary full-S3 reference; 64 independent trajectory pairs. Shading is ±1 replicate SE.\n'
              'Raw-touch / unchanged-orbit contribution is exactly zero for rotation triples. Returning to the initial state never resets visitation.\n'
              'Bins sum without survivor normalization. Post-encounter contributions may be negative; these are not quantum oscillations or evidence of binding.', fontsize=8.5)
output = Path('docs/images/triangle-encounter-memory.png')
fig.savefig(output, dpi=160)
print(output)
