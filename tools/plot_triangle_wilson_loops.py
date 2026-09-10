#!/usr/bin/env python3
"""Static scientific figure from the saved independent loop measurements."""
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from triangle_wilson_dynamics import LoopObserver


data = json.loads(Path('data/triangle-wilson-confirmation.json').read_text())
fig, (left, right) = plt.subplots(1, 2, figsize=(12, 4.6), gridspec_kw={'width_ratios': [1.1, 1]})
fig.suptitle('Spatial gauge flux depends on area, not just perimeter', fontsize=15, y=0.98)
fig.text(0.5, 0.905, 'Supplied triangular mesh · unchanged microscopic rules · stationary starts · no force fit',
         ha='center', fontsize=10, color='#555555')

case = data['experiments'][1]
theory = LoopObserver(8).theory
areas = np.arange(1, 13)
exact = [float(theory.mean(128, 32, int(a), 2)) for a in areas]
left.plot(areas, exact, color='#275b9b', label='Exact finite-volume prediction')
left.plot(areas, case['bulk_prediction']['character_area_factors'][2]**areas,
          '--', color='#999999', label='Infinite-volume area law')
for i, row in enumerate(case['loops']):
    offset = -0.12 if i == 4 else 0.12 if i == 5 else 0
    left.errorbar(row['area_faces']+offset, row['mean'][1], yerr=2*row['trajectory_standard_error'][1],
                  fmt='s' if i == 5 else 'o', color='#c66b22', ms=5, capsize=3,
                  label='Evolved links (±2 SE)' if i == 0 else None)
left.annotate('Same area 8\nperimeters 8 and 10', xy=(8, 0.10), xytext=(4.6, 0.065),
              arrowprops={'arrowstyle': '-', 'color': '#666666'}, fontsize=9)
left.set(yscale='log', xlabel='Enclosed area (faces)', ylabel='Mean normalized standard character',
         xlim=(0.6, 12.5), ylim=(0.02, 1.05), title='128 faces, total charge 32')
left.set_xticks([1, 2, 4, 6, 8, 12])
left.legend(fontsize=8.5, loc='upper right', frameon=False)

colors = ['#275b9b', '#c66b22', '#43866a']
for i, (case, color) in enumerate(zip(data['experiments'], colors)):
    offset = (i-1)*0.13
    for j, key in enumerate(('same_area_different_perimeter', 'same_perimeter_different_area')):
        comparison = case[key]
        right.errorbar(j+offset, comparison['paired_mean_difference'][1],
                       yerr=2*comparison['paired_standard_error'][1], fmt='o', color=color, capsize=3,
                       label=f"F={case['faces']}, Q={case['charge']}" if j == 0 else None)
        right.plot(j+offset, comparison['exact_difference'][1], marker='_', ms=14, color='black')
right.axhline(0, color='#aaaaaa', lw=0.8)
right.set_xticks([0, 1], ['Same area 8\nperimeters 10 − 8', 'Same perimeter 10\nareas 8 − 12'])
right.set(xlim=(-0.42, 1.42), ylim=(-0.009, 0.083), ylabel='Paired difference of mean characters',
          title='Matched geometric comparisons')
right.text(0.03, 0.96, 'Black ticks: exact predictions\nPoints: measurements ±2 SE',
           transform=right.transAxes, va='top', fontsize=9)
right.legend(loc='center left', fontsize=9, frameon=False)
for ax in (left, right):
    ax.spines[['top', 'right']].set_visible(False)
    ax.grid(axis='y', alpha=0.16)
fig.subplots_adjust(top=0.80, bottom=0.17, left=0.075, right=0.985, wspace=0.29)
output = Path('docs/images/triangle-wilson-loops.png')
fig.savefig(output, dpi=180)
print(output.resolve())
