#!/usr/bin/env python3
"""Plot actual charge fields and exact noise differences, without a fitted PDE."""
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import PolyCollection
from matplotlib.colors import BoundaryNorm, ListedColormap

data = json.loads(Path('data/triangle-charge-current.json').read_text())
charges = json.loads(Path('data/triangle-feedback.json').read_text())['bank']['element_charges']
runs = data['homogeneous_activation']['runs']
large = runs[-1]
side = large['side']
fig, axes = plt.subplots(2, 2, figsize=(13, 11), layout='constrained')
colors = ListedColormap(['#315e91', '#e4e7e9', '#d56c3b'])
norm = BoundaryNorm([-0.5, 0.5, 1.5, 2.5], colors.N)
polygons = []
for y in range(side):
    for x in range(side):
        polygons.extend([[(x, y), (x+1, y), (x+1, y+1)],
                         [(x, y), (x+1, y+1), (x, y+1)]])
for ax, snapshot, title in zip(axes[0], [large['snapshots'][0], large['snapshots'][-1]],
        ['Both initial charge fields; frozen control at all times', 'Disordered reflection connection: evolved charge field']):
    field = np.asarray(snapshot['charge_field'])
    collection = PolyCollection(polygons, array=field, cmap=colors, norm=norm,
                                edgecolors='#687782', linewidths=0.25)
    ax.add_collection(collection)
    ax.set(xlim=(0, side), ylim=(0, side), aspect='equal', xlabel='Supplied mesh x index', ylabel='Supplied mesh y index')
    ax.set_title(title, fontsize=12)
    n = snapshot['populations_0_1_2']
    ax.text(0.02, 0.98, f'Charge populations: {n[0]} / {n[1]} / {n[2]}\nTotal charge = {2*side*side}',
            transform=ax.transAxes, va='top', fontsize=10,
            bbox={'facecolor': 'white', 'alpha': 0.9, 'edgecolor': 'none'})
fig.colorbar(collection, ax=axes[0].tolist(), ticks=[0, 1, 2], label='Derived scalar face charge q')

ax = axes[1, 0]
for run in runs:
    count = 2*run['side']**2
    denominator = 13*6*run['side']**2
    histograms = run['checked']['runs'][3]['histograms']
    stride = run['attempts']//(len(histograms)-1)
    recorded = {i*stride: sum((q-1)**2*n for q, n in zip(charges, histogram))/count
                for i, histogram in enumerate(histograms)}
    recorded.update({s['tick']: s['charge_contrast_norm_squared']/count for s in run['snapshots']})
    t = [tick/denominator for tick in sorted(recorded)]
    contrast = [recorded[tick] for tick in sorted(recorded)]
    ax.plot(t, contrast, 'o-', alpha=0.55 if run['side'] == 6 else 1,
            linewidth=1.3 if run['side'] == 6 else 2.5,
            label=f"{run['side']}×{run['side']}, seed {run['seed']}")
ax.axhline(0, color='#60666d', linestyle='--', label='Frozen / reaction-disabled controls')
ax.set(xlabel='Attempted updates per rooted operator (supplied clock)',
       ylabel=r'Mean charge contrast $\sum_f(q_f-1)^2/F$',
       title='Fluctuations activate despite zero initial mean current')
ax.legend(fontsize=8, loc='lower right')
ax.grid(alpha=0.2)

ax = axes[1, 1]
faces = sorted({a for a, _, _ in data['readout']['raw_second_sum_difference']})
matrix = np.zeros((3, 3), dtype=int)
for a, b, n in data['readout']['raw_second_sum_difference']:
    matrix[faces.index(a), faces.index(b)] = n
heatmap = ax.imshow(matrix, cmap='RdBu_r', vmin=-8, vmax=8)
for (i, j), value in np.ndenumerate(matrix):
    ax.text(j, i, str(value), ha='center', va='center', fontsize=20,
            color='white' if abs(value) == 8 else 'black')
ax.set(xticks=range(3), yticks=range(3), xticklabels=faces, yticklabels=faces,
       xlabel='Face id', ylabel='Face id', title='Stored braid memory changes noise, not mean drift')
ax.text(0.5, -0.16, 'Exact covariance difference × 11,232\nSame full charge field and encounter-boundary holonomy',
        transform=ax.transAxes, ha='center', va='top', fontsize=10)
fig.colorbar(heatmap, ax=ax, label='Integer covariance-difference numerator')
fig.suptitle('Conserved charge currents and gauge-dependent classical fluctuations', fontsize=18)
fig.supxlabel('Actual primitive histories; geometry and update clock are supplied. No quantum amplitudes, emergent matter, or phase transition is claimed.', fontsize=10)
output = Path('docs/images/triangle-charge-current.png')
fig.savefig(output, dpi=160)
print(output)
