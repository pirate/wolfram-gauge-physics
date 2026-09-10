#!/usr/bin/env python3
"""Scientific plot of actual winding-current measurements, with no fitted law."""
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np


data = json.loads(Path('data/triangle-winding-fluctuations.json').read_text())
fig, axes = plt.subplots(2, 1, figsize=(9, 7.5), sharex=True, layout='constrained')
for run, color, marker in zip(data['runs'], ('#2463a6', '#b34b32'), ('o', 's')):
    rows = np.array(run['windows'])
    x = rows[:, 0]
    label = f"{run['faces']} faces; {run['trajectories']} independent trajectories"
    axes[0].errorbar(x, rows[:, 1], yerr=rows[:, 2], fmt=marker+'-', color=color,
                     capsize=3, lw=1.1, ms=5, label=label)
    axes[0].axhline(run['bare_xx'], color=color, ls=':', alpha=.65, lw=1)
    axes[0].axhline(run['exact_one_function_upper_bound'], color=color, ls='--', alpha=.65, lw=1)
    axes[1].errorbar(x, rows[:, 7], yerr=rows[:, 8], fmt=marker+'-', color=color,
                     capsize=3, lw=1.1, ms=5)
    axes[1].plot(x, run['iid_kurtosis_at_one_attempt_per_face']/x, ':', color=color, alpha=.75)
for ax in axes:
    ax.set_xscale('log', base=4)
    ax.grid(True, alpha=.18)
    ax.spines[['top', 'right']].set_visible(False)
axes[0].set_title('Autonomous raw-link dynamics: winding-current fluctuations', loc='left', weight='bold')
axes[0].set_ylabel(r'$\mathbb{E}[J_N^2]/(2N)$')
axes[0].legend(loc='upper left', fontsize=9, frameon=False)
axes[0].text(.02, .04, 'Dotted: instantaneous estimate. Dashed: exact one-function mobility bound.',
             transform=axes[0].transAxes, fontsize=8)
axes[0].set_ylim(.0195, .028)
axes[1].axhline(0, color='#555555', lw=.8)
axes[1].set_ylabel('Excess kurtosis of winding current')
axes[1].set_xlabel(r'Window length $N/F$ (attempted updates per face)')
axes[1].set_title('Short-window correlations exceed an independent-increment null', loc='left', fontsize=11)
axes[1].text(.02, .04, 'Dotted: exact one-step marginal with independent increments. Zero: Gaussian kurtosis.',
             transform=axes[1].transAxes, fontsize=8)
axes[1].set_ylim(-.48, .4)
fig.suptitle('285 million forward proposals; error bars are one trajectory-cluster standard error.\n'
             'Supplied triangular geometry and finite fiber; no fitted diffusion law. Long-window convergence remains open.',
             fontsize=10)
output = Path('docs/images/triangle-winding-fluctuations.png')
fig.savefig(output, dpi=170)
print(output.resolve())
