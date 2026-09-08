#!/usr/bin/env python3
"""Plot exact response data; supplied mesh coordinates, no fitted dynamics."""
import json
from fractions import Fraction
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import PolyCollection
from matplotlib.colors import TwoSlopeNorm
from matplotlib.ticker import NullFormatter

data = json.loads(Path('data/triangle-charge-response.json').read_text())
memory = json.loads(Path('data/triangle-response-memory.json').read_text())
number = lambda pair: float(Fraction(*pair))
fig, axes = plt.subplots(1, 3, figsize=(16, 5.3), layout='constrained')
k = np.linspace(0, np.pi, 200)
lam = 3-np.abs(2+np.exp(1j*k))
for charge, color, label in ((2, '#287c8e', 'Two defects: exact mean closure'),
                              (72, '#bb6633', 'Reactive density Q/F = 1: projection only')):
    r = next(r for r in data['records'] if r['side'] == 6 and r['charge'] == charge)
    g = number(r['nearest_rate'])*lam+number(r['next_nearest_rate'])*(6*lam-lam**2)
    axes[0].plot(k, g, color=color, label=label, lw=2)
axes[0].set(title='Derived acoustic response symbol', xlabel=r'Cell wavevector $k_x$ ($k_y=0$)',
            ylabel=r'Projected decay $g_-(k)$ in scaled-operator units')
axes[0].set_xticks([0, np.pi/2, np.pi], ['0', r'$\pi/2$', r'$\pi$'])
axes[0].legend(fontsize=8, loc='upper left')
axes[0].grid(alpha=0.2)

rows = [r for r in data['records'] if r['charge'] == 4]
f = np.array([r['faces'] for r in rows])
residual = [number(r['closure_residual_squared_norm_per_face']) for r in rows]
limit = number(data['sparse_charge4_polynomial_certificate']['limit_F_cubed_times_defect'])
axes[1].loglog(f, residual, 'o-', color='#bb6633', lw=2, label='Exact four-charge residual')
axes[1].loglog(f, limit/f.astype(float)**3, '--', color='#555555', label=r'Derived leading term $(26568/5)F^{-3}$')
axes[1].set(title='Reactive correction becomes dilute', xlabel='Number of supplied faces F; fixed Q = 4',
            ylabel=r'Per-face unresolved squared drift $B_{ff}$')
axes[1].set_xticks(f, [str(x) for x in f])
axes[1].xaxis.set_minor_formatter(NullFormatter())
axes[1].legend(fontsize=8, loc='lower left')
axes[1].grid(alpha=0.2, which='both')

polygons = []
side = memory['side']
for y in range(side):
    for x in range(side):
        polygons.extend([[(x, y), (x+1, y), (x+1, y+1)],
                         [(x, y), (x+1, y+1), (x, y+1)]])
difference = np.zeros(2*side*side)
for face, value in memory['drift_of_drift_difference']:
    difference[face] = value
col = PolyCollection(polygons, array=difference, cmap='RdBu_r',
                     norm=TwoSlopeNorm(vmin=-48, vcenter=0, vmax=48),
                     edgecolors='#aaaaaa', linewidths=0.6)
axes[2].add_collection(col)
centers = np.mean(np.array(polygons), axis=1)
active = np.array([f for f, _ in memory['drift_of_drift_difference']])
for face in active:
    axes[2].text(*centers[face], f'{int(difference[face]):+}', ha='center', va='center', fontsize=9,
                 color='white' if abs(difference[face]) >= 40 else 'black')
low, high = centers[active].min(axis=0)-0.8, centers[active].max(axis=0)+0.8
axes[2].set(xlim=(low[0], high[0]), ylim=(low[1], high[1]), aspect='equal',
            title='Same charge and drift; different next response',
            xlabel='Supplied mesh x index', ylabel='Supplied mesh y index')
fig.colorbar(col, ax=axes[2], shrink=0.75, label=r'$M^2$ times the two-attempt mean difference')
fig.suptitle('From exact classical transport to a measurable hidden-gauge correction', fontsize=16)
fig.supxlabel('Exact finite-model calculations, not fitted trajectories. Geometry and random scheduling are supplied.\n'
              'The reactive symbol is an initial projection, not a closed diffusion law; the right panel uses actual prepared connections.', fontsize=10)
output = Path('docs/images/triangle-charge-response.png')
fig.savefig(output, dpi=160)
print(output)
