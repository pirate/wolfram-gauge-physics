#!/usr/bin/env python3
"""Show actual based reflection stars from the checked raw-connection controls."""
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from triangle_reference import ActivityProbe, subgroup

data = json.loads(Path('data/triangle-reference.json').read_text())
bank = json.loads(Path('data/triangle-feedback.json').read_text())['bank']
probe = ActivityProbe(3, bank)
fig, axes = plt.subplots(1, 3, figsize=(13, 5.4))
colors = {1: '#386cb0', 2: '#d99026', 5: '#7757a5'}
labels = {1: 'A', 2: 'B', 5: 'C'}
angles = np.pi/2-np.arange(6)*np.pi/3
xy = np.array([np.cos(angles), np.sin(angles)]).T
for ax, index, title in zip(axes, (1, 2, 8),
        ('Equal reflections: frozen', 'Distinct neighbors: frozen', 'Equality boundary: active')):
    row = data['activity']['compiled_initial_states'][index]
    raw = row['links']
    values = [probe.e.factor.oracle.fan.transport(raw, path) for path in probe.stars[0]]
    if any(bank['element_charges'][h] != 1 for h in values):
        raise ValueError('illustrated star is not uniformly charge one')
    equal = [values[i] == values[(i+1) % 6] for i in range(6)]
    enabled = [(i+1) % 6 for i in range(6) if equal[i] != equal[(i+1) % 6]]
    for i in range(6):
        j = (i+1) % 6
        ax.plot(*np.array([xy[i], xy[j]]).T, color='#b4b8bd', lw=2, zorder=1)
        midpoint = 1.25*(xy[i]+xy[j])/2
        ax.text(*midpoint, '=' if equal[i] else '≠', ha='center', va='center', fontsize=15)
    ax.scatter(*xy.T, s=1000, c=[colors[v] for v in values], edgecolors='white', linewidths=2, zorder=3)
    if enabled:
        ax.scatter(*xy[enabled].T, s=1450, facecolors='none', edgecolors='#c53f38', linewidths=2.5, zorder=4)
    for point, value in zip(xy, values):
        ax.text(*point, labels[value], color='white', ha='center', va='center', fontsize=17, zorder=5)
    ax.text(0, 0, r'$q_f=1$'+'\non every face', ha='center', va='center', fontsize=12)
    group_order = len(subgroup(probe.e.geometry.group, probe.e.forest.based_loops(raw)[0]))
    ax.text(0, -1.46, f'{len(enabled)} active fans at this vertex\n'
            f"{row['activity']['changing_operators']} changing operators on the whole mesh\n"
            f'Global holonomy group order: {group_order}', ha='center', va='top', fontsize=10)
    ax.set(title=title, xlim=(-1.4, 1.4), ylim=(-2.1, 1.4), aspect='equal')
    ax.axis('off')
fig.suptitle('Identical charge fields can conceal different microscopic activity', fontsize=17, y=0.96)
fig.text(0.5, 0.11, 'Circles are the six actual face holonomies based at vertex 0; red rings mark active fan centers.\n'
         'A, B, C are the three reflections in a common frame. Relabeling that frame preserves equality and activity.',
         ha='center', fontsize=10)
fig.text(0.5, 0.025, 'Checked raw-link states on an 18-face supplied torus. Cyclic order is exact; drawing coordinates are schematic, not emergent space.',
         ha='center', fontsize=9, color='#50565e')
fig.subplots_adjust(top=0.82, bottom=0.19, wspace=0.2)
output = Path('docs/images/triangle-reference.png')
fig.savefig(output, dpi=160)
print(output)
