#!/usr/bin/env python3
"""Scientific plots of recorded microscopic trajectories; no resimulation or fits hidden here."""
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

main = json.loads(Path('data/d4-diffusion.json').read_text())
replication = json.loads(Path('data/d4-diffusion-replication.json').read_text())
colors = {0.25: '#1768a4', 0.5: '#c35e27', 0.75: '#4d8871'}
fig, axes = plt.subplots(1, 3, figsize=(13.2, 4.7), constrained_layout=True)
for data, marker in ((main, 'o'), (replication, '^')):
    for ensemble in data['ensembles']:
        rho = ensemble['density']
        if rho not in colors:
            continue
        times = ensemble['times'][1:]
        variances = [o['variance'] for o in ensemble['observables'][1:]]
        axes[0].loglog(times, variances, marker, color=colors[rho], markersize=4,
                       label=fr'$\rho={rho}$' if data is main else None)
        if data is replication:
            axes[0].loglog(times, ensemble['finite_time_prediction']['variances'][1:],
                           '-', color=colors[rho], alpha=.6)
        fit = ensemble['late_fit']
        value, (low, high) = fit['diffusion'], fit['diffusion_ci95']
        axes[1].errorbar(rho, value, yerr=[[value-low], [high-value]], fmt=marker,
                         markerfacecolor='white' if data is main else colors[rho],
                         color=colors[rho], capsize=3)
densities = [.2+i*.006 for i in range(101)]
axes[1].plot(densities, [(1/rho-1)/2 for rho in densities], color='#444444', alpha=.7,
             label=r'$D=(1-\rho)/(2\rho)$')
axes[1].set_yscale('log')
axes[1].legend(frameon=False, fontsize=9)
ensemble = next(e for e in replication['ensembles'] if e['density'] == .5)
width = 4
hist = ensemble['displacement_histograms'][-1]
pmf = ensemble['exact_late_prediction']['pmf']
positions = [i-replication['layers']+(width-1)/2 for i in range(0, len(hist), width)]
axes[2].plot(positions, [sum(hist[i:i+width])/ensemble['trials'] for i in range(0, len(hist), width)],
             'o', markersize=3, color=colors[.5], label='Microscopic trajectories')
axes[2].plot(positions, [sum(pmf[i:i+width]) for i in range(0, len(pmf), width)],
             '-', color='#444444', label='Exact finite-time probability')
axes[2].set_xlim(-90, 90)
axes[2].legend(frameon=False, fontsize=8)
axes[0].set(title='Diffusive spreading', xlabel='Single matching layers', ylabel='Tag displacement variance')
axes[0].legend(frameon=False, fontsize=9)
axes[1].set(title='Density-dependent transport', xlabel=r'Background occupancy $\rho$', ylabel='Half the late variance slope')
axes[2].set(title=r'Full distribution: $\rho=0.5$, layer 512', xlabel='Tag displacement (cells)', ylabel='Probability per 4-cell bin')
for ax in axes:
    ax.grid(alpha=.18)
fig.suptitle('Deterministic gauge-table evolution → an exactly identified exclusion process', fontsize=14)
fig.supxlabel('Circles: 256-layer run. Triangles: independent 512-layer run. Left lines: exact finite-time variance.\n'
              'Error bars: pointwise 95% block bootstrap; 5/6 cover theory. '
              'Fixed 1D chain; random initial ensemble only. Empty/full controls: ballistic/blocked.', fontsize=8)
output = Path('docs/images/derived-diffusion.png')
fig.savefig(output, dpi=180)
print(output)
