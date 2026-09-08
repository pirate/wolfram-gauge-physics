#!/usr/bin/env python3
"""Plot actual seeded charge histories and explicitly mark periodic seam contact."""
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import PolyCollection
from matplotlib.colors import BoundaryNorm, ListedColormap

data = json.loads(Path('data/triangle-spreading.json').read_text())
case = next(c for c in data['cases'] if c['side'] == 24 and c['trial'] == 0)
run = case['runs'][0]['evolution']
side, operators = case['side'], case['attempted_operators']
polygons = []
for y in range(side):
    for x in range(side):
        polygons.extend([[(x, y), (x+1, y), (x+1, y+1)],
                         [(x, y), (x+1, y+1), (x, y+1)]])
fig, axes = plt.subplots(2, 2, figsize=(13, 10), layout='constrained')
first = np.ma.array([t/operators if t is not None else 0 for t in run['first_charge_departure_tick']],
                    mask=[t is None for t in run['first_charge_departure_tick']])
col = PolyCollection(polygons, array=first, cmap='viridis', edgecolors='#657080', linewidths=0.12)
axes[0, 0].add_collection(col)
fig.colorbar(col, ax=axes[0, 0], label='First charge departure: attempts per rooted operator')
axes[0, 0].set_title('One-link seed: when each face first changes charge')
qmap = ListedColormap(['#315e91', '#eceeef', '#d56c3b'])
col = PolyCollection(polygons, array=np.array(run['snapshots'][-1]['charge_field']), cmap=qmap,
                     norm=BoundaryNorm([-0.5, 0.5, 1.5, 2.5], qmap.N), edgecolors='#657080', linewidths=0.12)
axes[0, 1].add_collection(col)
fig.colorbar(col, ax=axes[0, 1], ticks=[0, 1, 2], label='Derived face charge')
axes[0, 1].set_title('Final charge field; initially every face had q = 1')
for ax in axes[0]:
    ax.plot([side//2, side//2+1], [side//2, side//2], color='black', lw=3, label='Initially changed link')
    ax.set(xlim=(0, side), ylim=(0, side), aspect='equal', xlabel='Supplied mesh x index', ylabel='Supplied mesh y index')
    ax.legend(loc='upper left', fontsize=8)

colors = {6: '#4477aa', 12: '#cc8844', 24: '#228877'}
for c in data['cases']:
    n, m = 2*c['side']**2, c['attempted_operators']
    for i, row in enumerate(c['runs']):
        r = row['evolution']
        t = np.array([s['tick']/m for s in r['snapshots']])
        radius = np.array([s['ever_contrast_radius'] if s['ever_contrast_radius'] is not None else np.nan for s in r['snapshots']])
        seam = r['first_changing_support_touching_periodic_seam']
        before = np.array([seam is None or s['tick'] < seam for s in r['snapshots']])
        label = f'{n:,} faces' if c['trial'] == 0 and i == 0 else None
        color = colors[c['side']]
        axes[1, 0].plot(t, radius, ':', color=color, alpha=0.3, linewidth=1)
        axes[1, 0].plot(t[before], radius[before], '-', color=color, alpha=0.75, linewidth=1.5, label=label)
        fraction = [(s['activity']['populations'][0]+s['activity']['populations'][2])/n for s in r['snapshots']]
        axes[1, 1].plot(t, fraction, color=color, alpha=0.3, linewidth=1, label=label)
axes[1, 0].set(xlabel='Attempts per rooted operator (supplied clock)', ylabel='Largest visited dual-graph distance from seed',
               title='Spreading radius; dotted after first seam-touching event')
axes[1, 1].axhline(0, color='#555555', linestyle='--', label='Unperturbed / reaction-disabled controls')
axes[1, 1].set(xlabel='Attempts per rooted operator (supplied clock)', ylabel='Fraction of faces with charge different from one',
               title='Activity creates charge contrast without changing total charge')
for ax in axes[1]:
    ax.grid(alpha=0.2)
    ax.legend(fontsize=8)
fig.suptitle('A local gauge perturbation activates a frozen nonabelian background', fontsize=17)
fig.supxlabel('Actual primitive evolution, not an inserted wave law. Six seed variants share each schedule; only two independent schedules per mesh size.\n'
              'Periodic-box saturation is visible; these runs do not establish a speed law, a phase transition, or a localized particle.', fontsize=10)
output = Path('docs/images/triangle-spreading.png')
fig.savefig(output, dpi=160)
print(output)
