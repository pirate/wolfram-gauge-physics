#!/usr/bin/env python3
"""Draw the actual transport polygon and its measured invariant outputs."""
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
import numpy as np

from transport_braid import TransportExperiment


def main():
    data = json.loads(Path('data/d4-transport-braid.json').read_text())
    record = data['experiments'][0]
    e = TransportExperiment(data['side'])
    side = e.geometry.side
    points = np.array([(v % side, v//side) for v in range(e.geometry.size)], dtype=float)
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.6), constrained_layout=True,
                             gridspec_kw={'width_ratios': [1.2, 1, 1]})
    ax = axes[0]
    for u, v in e.geometry.edges:
        a, b = points[u], points[v]
        if np.max(np.abs(a-b)) <= 1:
            ax.plot([a[0], b[0]], [a[1], b[1]], color='#dddddd', linewidth=.5, zorder=0)
    for f, c in enumerate(record['states'][0]['face_classes']):
        if c:
            ax.add_patch(Polygon(points[list(e.geometry.faces[f][:3])], facecolor='#3278b7' if c == 1 else '#da7926', alpha=.85))
    path = np.asarray([e.point(f) for f in record['cycle']], dtype=float)/3
    ax.plot(path[:, 0], path[:, 1], '.-', color='#7145a0', linewidth=2)
    for k in (0, 3, 6, 9):
        a, b = path[k:k+2]
        ax.annotate('', xy=b, xytext=a, arrowprops={'arrowstyle': '->', 'color': '#7145a0', 'lw': 1.8})
    ax.annotate('moving r', xy=path[0], xytext=(4.0, 2.5), arrowprops={'arrowstyle': '-', 'color': '#555'}, fontsize=9)
    center = np.asarray(e.point(record['center_face']))/3
    ax.annotate('enclosed s', xy=center, xytext=(8.0, 8.0), arrowprops={'arrowstyle': '-', 'color': '#555'}, fontsize=9)
    ax.annotate('r and s anchors', xy=(.5, .5), xytext=(1.5, 2.1), arrowprops={'arrowstyle': '-', 'color': '#555'}, fontsize=9)
    ax.set(xlim=(-.25, side-.75), ylim=(-.25, side-.75), aspect='equal', xlabel='Supplied mesh x index', ylabel='Supplied mesh y index',
           title='Twelve actual vacancy-transport updates')
    ax.text(.02, .98, 'Periodic seams omitted\nAll defect positions return', transform=ax.transAxes, va='top', fontsize=9)

    ax = axes[1]
    for relation, label, color in zip(record['same_class_loop_relations'], ('r-pair relation', 's-pair relation'), ('#3278b7', '#da7926')):
        ax.plot(range(3), [int(x != 0) for x in relation['central_holonomy_by_lap']], 'o-', color=color, label=label)
    ax.set(xticks=range(3), xlabel='Completed circuits', yticks=[0, 1], yticklabels=['identity', 'central half-turn z'],
           ylim=(-.2, 1.4), title='Gauge-invariant memory flips and returns')
    ax.legend(loc='upper center', fontsize=9)
    ax.grid(alpha=.2)

    ax = axes[2]
    values = np.array([s['spectrum']['above_flat_band']['eigenvalues'] for s in record['states']])
    for j in range(2):
        ax.plot(range(3), values[:, j], 'o-', label=f'Above-band eigenvalue {j+1}')
    ax.set(xticks=range(3), xlabel='Completed circuits', ylabel='Two-component connection eigenvalue',
           ylim=(9.025, 9.145), title='Same face classes; different exact spectrum')
    ax.legend(fontsize=8, loc='upper center')
    ax.grid(alpha=.2)
    ax.text(.04, .04, r'First exact trace difference: $\mathrm{tr}(L^{14})$'+'\nDifference after one circuit: −25,200',
            transform=ax.transAxes, fontsize=8)
    fig.suptitle('Classical holonomy memory from primitive link transport — not quantum braid amplitudes', fontsize=14)
    fig.savefig('docs/images/transport-braid.png', dpi=180)


if __name__ == '__main__':
    main()
