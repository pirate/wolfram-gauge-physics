#!/usr/bin/env python3
"""Plot measured raw-link residence times; no synthesized class trajectories."""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from bank_equilibrium import BankFactor, residence


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', type=Path, default=Path('data/d4-bank-equilibrium.json'))
    parser.add_argument('--output', type=Path, default=Path('docs/images/bank-equilibrium.png'))
    args = parser.parse_args()
    data = json.loads(args.input.read_text())
    experiment = data['relaxation']
    factor = BankFactor(data['side'], json.loads(Path('data/d4-triple-channels.json').read_text()), data['family'])
    colors = ['#007f86', '#b85520', '#7252a3']
    fig, (left, right) = plt.subplots(1, 2, figsize=(11.8, 4.5), gridspec_kw={'width_ratios': [1.6, 1]})
    exact = sum(row['heavy']*row['probability'][0]/row['probability'][1] for row in data['formal_canonical_counts'])
    attempts = experiment['attempts_per_run']
    step = 1000
    for h, color in enumerate(colors):
        selected = [r for r in experiment['runs'] if r['initial_heavy'] == h]
        curves = []
        for run in selected:
            blocks = residence(factor, run['initial_links'], run['events'], attempts, 0, step)['block_occupation_attempts']
            total, curve = 0, []
            for i, row in enumerate(blocks, 1):
                total += sum(k*n for k, n in enumerate(row))
                curve.append(total/(i*step))
            curves.append(curve)
        mean = [sum(values)/len(values) for values in zip(*curves)]
        left.plot(range(step, attempts+1, step), mean, color=color, label=f'Start $N_h={h}$ (4 runs)')
    left.axhline(exact, color='#202c39', ls='--', lw=1.2, label=r'Exact $62/241$')
    left.axvline(experiment['burn_attempts'], color='#888888', ls=':', lw=1)
    left.set(xscale='log', xlabel='Attempted updates (idle slots included)',
             ylabel=r'Running time average of $N_h$', title='Raw-link trajectories approach the counting prediction')
    left.legend(fontsize=8, loc='best')
    for i, run in enumerate(experiment['runs']):
        right.scatter(i, 1000*run['conversion_rate_per_attempt'], color=colors[run['initial_heavy']], s=28)
    rate = data['canonical_reaction_rates']['conversions_per_attempt']
    right.axhline(1000*rate[0]/rate[1], color='#202c39', ls='--', lw=1.2, label=r'Exact $440/86037$ per attempt')
    right.set(xticks=[1.5, 5.5, 9.5], xticklabels=[r'$N_h=0$', r'$N_h=1$', r'$N_h=2$'],
              xlabel='Initial population; each dot is one independent run',
              ylabel='Conversions per 1,000 attempted updates', title='A second prediction: conversion frequency')
    right.legend(fontsize=8, loc='best')
    for ax in (left, right):
        ax.spines[['top', 'right']].set_visible(False)
        ax.grid(axis='y', alpha=.2)
    fig.suptitle('Finite nonabelian gauge-link statistics: 18 faces, fixed charges (4, 2, 0)', fontsize=13)
    fig.text(.5, .015, 'Right panel: 200,000-attempt burn, then 800,000 attempts per run. Supplied torus and random clock; no molecular or quantum claim.',
             ha='center', fontsize=8, color='#444444')
    fig.tight_layout(rect=(0, .055, 1, .94))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=180)


if __name__ == '__main__':
    main()
