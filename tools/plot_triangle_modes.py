#!/usr/bin/env python3
"""Measured full-graph mode densities, size scaling, and exact primitive evolution."""
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np


def main():
    data = json.loads(Path('data/triangle-modes.json').read_text())
    rows = [r for r in data['static'] if r['condition'] == 'noncommuting_pair']
    sample = next(r for r in rows if r['side'] == 24)
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), layout='constrained')
    maximum = max(max(m['base_projector_density']) for m in sample['numerical_modes'])
    for ax, mode, label in zip(axes[0], sample['numerical_modes'], ('Lower', 'Upper')):
        density = np.array(mode['base_projector_density']).reshape(24, 24)
        plotted = ax.imshow(density, origin='lower', cmap='magma', norm=LogNorm(vmin=1e-6, vmax=maximum))
        ax.plot([12, 13, 14], [12, 12, 12], 'o-', c='#78e0ce', lw=2, markersize=4)
        ax.set(title=f'{label} full-graph mode: λ ≈ {mode["eigenvalue_estimate"]:.6f}',
               xlabel='Supplied mesh x index', ylabel='Supplied mesh y index')
        ax.text(.03, .97, f'Participation volume ≈ {mode["base_participation_volume"]:.2f}',
                transform=ax.transAxes, va='top', color='white', fontsize=9)
    fig.colorbar(plotted, ax=list(axes[0]), label='Normalized graph-mode density (log scale)', shrink=.9)

    ax = axes[1, 0]
    trace = np.array([[r[0], r[2], min(r[3:])] for r in data['dynamics']['seed_zero_full_changing_event_counts']])
    ax.step(trace[:, 0], trace[:, 2], where='post', color='#c08b32', lw=1.5, ls='--', label='Exact orientation capacity')
    ax.step(trace[:, 0], trace[:, 1], where='post', color='#3c78a1', lw=1.2, label='Exact above-band count')
    for r in data['dynamics']['reactions']:
        if r['seed'] == 0:
            a, b = [x['positive'] for x in r['exact_inertias']]
            ax.scatter(r['event'][0], b, marker='v' if b < a else '^' if b > a else 'o',
                       c='#b94037' if b < a else '#278158' if b > a else '#777', s=35, zorder=4)
    ax.set(xlim=(0, 24000), ylim=(-.15, 2.4), yticks=(0, 1, 2), xlabel='Attempted-update tick',
           ylabel='Mode count / upper bound', title='Reactions can force full-band mode loss')
    ax.legend(fontsize=8, loc='upper right')
    ax.grid(alpha=.15)

    ax = axes[1, 1]
    for j, label, color in ((0, 'Lower mode', '#8d6ab1'), (1, 'Upper mode', '#348c7e')):
        ax.plot([r['side']**2 for r in rows], [r['numerical_modes'][j]['base_participation_volume'] for r in rows],
                'o-', label=label, color=color)
    ax.set(xscale='log', xlabel='Base vertices in the finite torus', ylabel='Participation volume in base vertices',
           title='Densities remain concentrated as the graph grows')
    ax.legend(fontsize=9)
    ax.grid(alpha=.15)
    fig.suptitle('Isolated modes in the full triangle-fiber graph — static localization is not particle stability', fontsize=14)
    fig.supxlabel('The integer certificate proves exactly two infinite-graph modes above 12.05. Densities and decimal eigenvalues are numerical.\n'
                  'Unrouted seed zero: no above-band modes for 26,940 of 100,000 attempt intervals. No physical wave law is assumed.', fontsize=9)
    fig.savefig('docs/images/triangle-modes.png', dpi=180)


if __name__ == '__main__':
    main()
