#!/usr/bin/env python3
"""Actual overlap geometry and exact conditional mean response, not a mockup."""
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import PolyCollection
from matplotlib.colors import TwoSlopeNorm
from matplotlib.patches import Patch

data = json.loads(Path('data/triangle-patch-observer.json').read_text())
side = data['side']
def point(v):
    y, x = divmod(v, side)
    return ((x+side//2) % side-side//2, (y+side//2) % side-side//2)
polygons = [[point(v) for v in vertices] for _, vertices in data['face_vertices']]
faces = [f for f, _ in data['face_vertices']]
patches = [set(fs) for fs in data['patch_face_ids']]
palette = ['#569cc7', '#df9a52', '#9975b1']
colors = [palette[2 if f in patches[0] & patches[1] else 0 if f in patches[0] else 1] for f in faces]
fig, axes = plt.subplots(1, 2, figsize=(12, 6), layout='constrained')
axes[0].add_collection(PolyCollection(polygons, facecolors=colors, edgecolors='#333333', linewidths=1.2))
axes[0].set_title('Two complete patch descriptions still lose alignment')
axes[0].legend(handles=[Patch(color=c, label=s) for c, s in zip(palette, ('First fan only', 'Second fan only', 'Shared face'))],
               fontsize=9, loc='lower right')
centers = np.mean(np.array(polygons), axis=1)
for face, center in zip(faces, centers):
    axes[0].text(*center, f'face {face}', ha='center', va='center', fontsize=10)
values = data['mean_response_witness']['differences'][-1]['difference_numerator']
col = PolyCollection(polygons, array=np.array(values), cmap='RdBu_r',
                     norm=TwoSlopeNorm(vmin=-48, vcenter=0, vmax=48), edgecolors='#333333', linewidths=1.2)
axes[1].add_collection(col)
axes[1].set_title('Different gluing changes mean charge at attempt three')
for value, center in zip(values, centers):
    axes[1].text(*center, f'{value:+}', ha='center', va='center', fontsize=12,
                 color='white' if abs(value) >= 40 else 'black')
fig.colorbar(col, ax=axes[1], shrink=0.65, label=r'$44^3$ times the conditional mean difference')
for ax in axes:
    ax.scatter(*zip(*(point(v) for v in data['vertices'])), s=22, color='#222222', zorder=4)
    ax.set(xlim=(-1.3, 1.3), ylim=(-1.3, 1.3), aspect='equal',
           xlabel='Supplied, locally unwrapped mesh x', ylabel='Supplied, locally unwrapped mesh y')
fig.suptitle('Exact overlapping-patch gauge dynamics: keep the gluing information', fontsize=15)
fig.supxlabel('845 compatible pairs of patch orbits expand to 1,393 union states when relative alignment is retained.\n'
              'Right: exact controlled evolution with 44 boundary-contained operators; not the full-mesh scheduler or emergent geometry.', fontsize=10)
output = Path('docs/images/triangle-patch-observer.png')
fig.savefig(output, dpi=160)
print(output)
