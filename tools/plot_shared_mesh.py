#!/usr/bin/env python3
"""Plot measured shared-face evolution and the exact fixed-boundary witness."""
import json
import math
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, Patch

data = json.loads(Path('data/d4-shared-mesh-replication.json').read_text())
colors = {'flat': '#999999', 'reflection': '#1768a4', 'rotation': '#4d8871',
          'noncommuting_pair': '#c35e27', 'random_links': '#8763a2'}
names = {'flat': 'Flat control', 'reflection': 'One reflection link', 'rotation': 'One rotation link',
         'noncommuting_pair': 'Two noncommuting seed links', 'random_links': 'Uniform random links'}
fig, axes = plt.subplots(1, 2, figsize=(11.6, 4.6), constrained_layout=True)
for condition, color in colors.items():
    runs = [r for r in data['runs'] if r['condition'] == condition]
    times = [row[0] for row in runs[0]['trajectory']]
    fractions = [[row[1]/data['faces'] for row in run['trajectory']] for run in runs]
    average = [sum(values)/len(values) for values in zip(*fractions)]
    axes[0].plot(times, average, color=color, label=names[condition], linewidth=1.8)
    axes[0].fill_between(times, [min(v) for v in zip(*fractions)], [max(v) for v in zip(*fractions)],
                         color=color, alpha=.1)
    if condition in ('reflection', 'noncommuting_pair'):
        entropy = [[-sum((count/data['faces'])*math.log(count/data['faces']) for count in row[5:] if count)
                    for row in run['trajectory']] for run in runs]
        mean = [sum(values)/len(values) for values in zip(*entropy)]
        axes[1].plot(times, mean, color=color, label=names[condition])
        if not all(r['reverse_histogram_echo'] for r in runs):
            raise ValueError('cannot display an unverified inverse echo')
        axes[1].plot([2*data['layers']-t for t in reversed(times)], list(reversed(mean)), '--', color=color)
        axes[1].axhline(data['stationary_reference'][condition]['stationary_sector_entropy_nats'],
                        color=color, linewidth=.6, alpha=.4)
axes[0].axhline(.75, color='#666666', linewidth=.8, linestyle=':')
axes[0].axhline(.875, color='#666666', linewidth=.8, linestyle=':')
axes[0].set(title='Local seeds spread on a shared-link mesh', xlabel='Conflict-free update layers', ylabel='Fraction of nonflat faces')
axes[0].legend(frameon=False, fontsize=8, loc='lower right')
axes[1].axvline(data['layers'], color='#777777', linewidth=.8, linestyle=':')
axes[1].text(data['layers'], .3, 'Reverse the microscopic schedule', rotation=90, va='bottom', ha='right', fontsize=8)
axes[1].set(title='Coarse spreading is exactly reversible', xlabel='Forward layers, then inverse layers', ylabel='Empirical face-sector entropy (nats)')
axes[1].legend(frameon=False, fontsize=8, loc='lower center')
for ax in axes: ax.grid(alpha=.18)
fig.suptitle('Same derived gauge table; now every incident face responds', fontsize=14)
fig.supxlabel('1,152 faces; eight fixed schedules. Shading: min–max across schedules, not a confidence interval.\n'
              'Reference lines come from invariant uniform-link measures. Local histogram agreement does not prove thermalization; geometry is supplied.', fontsize=8)
fig.savefig('docs/images/shared-mesh-relaxation.png', dpi=180)

local = data['local_census']
witness = local['fixed_boundary_witness']
positions = {0: (0, 0), 1: (1, 0), 2: (1, 1), 3: (0, 1), 4: (.5, -.8), 5: (.5, 1.8)}
palette = {0: '#e5e7eb', 1: '#74b3d5', 5: '#e7a66e'}
fig, axes = plt.subplots(1, 3, figsize=(10.6, 4.6), constrained_layout=True)
states = [witness['states'][0]['initial_face_sectors']]+[s['final_face_sectors'] for s in witness['states']]
titles = ['Identical incoming face classes', 'First connection after update', 'Second connection after update']
for ax, sectors, title in zip(axes, states, titles):
    for label, face, sector in zip('ABCD', local['faces'], sectors):
        points = [positions[v] for v in face[:-1]]
        ax.add_patch(Polygon(points, facecolor=palette[sector], edgecolor='#555555', linewidth=1))
        x, y = [sum(p[i] for p in points)/3 for i in range(2)]
        ax.text(x, y, label, ha='center', va='center', fontsize=13)
    for u, v in ((0, 1), (2, 3)):
        ax.plot([positions[u][0], positions[v][0]], [positions[u][1], positions[v][1]], color='#b82435', linewidth=3)
    for v, (x, y) in positions.items(): ax.text(x, y, str(v), ha='right', va='bottom', fontsize=8)
    ax.set(xlim=(-.15, 1.15), ylim=(-.95, 1.95), title=title)
    ax.set_aspect('equal'); ax.axis('off')
fig.legend(handles=[Patch(color=palette[0], label='Flat'), Patch(color=palette[1], label='Reflection sector'),
                    Patch(color=palette[5], label='Central half-turn')], loc='lower center', ncol=3, frameon=False,
           bbox_to_anchor=(.5, .045), fontsize=9)
fig.suptitle('Identical exterior links + identical face classes ≠ identical future face classes', fontsize=13)
fig.supxlabel('A and B are the selected pair; red links are written. C and D respond through shared incidence.\n'
              'The two initial connections are gauge-inequivalent. Positions show supplied incidence, not inferred physical geometry.', fontsize=8)
fig.savefig('docs/images/shared-face-witness.png', dpi=180)
print('docs/images/shared-mesh-relaxation.png\ndocs/images/shared-face-witness.png')
