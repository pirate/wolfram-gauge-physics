#!/usr/bin/env python3
"""Plot boundary-censored first-passage observations without fitting a speed law."""
import json
from fractions import Fraction
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

data = json.loads(Path('data/triangle-first-passage.json').read_text())
fig, axes = plt.subplots(1, 2, figsize=(13, 5.5), layout='constrained')
radii = [2, 4, 6, 8, 10]
for record in data['records']:
    if record['side'] != 24:
        continue
    ys = [float(Fraction(*v['conditional_mean_normalized_time'])) if v else np.nan for _, v in record['passages']]
    axes[0].plot(radii, ys, color='#228877', alpha=0.16, linewidth=0.8)
for side, color, marker in ((24, '#228877', 'o'), (12, '#cc8844', 's')):
    rows = [r for r in data['summary'] if r['side'] == side]
    observed = [r for r in rows if not r['censored']]
    axes[0].plot([r['radius'] for r in observed],
                 [float(Fraction(*r['mean_conditional_residence'])) for r in observed],
                 color=color, marker=marker, markersize=8 if side == 12 else 5,
                 markerfacecolor='none' if side == 12 else color,
                 linewidth=1.5, label=f'{2*side*side:,} faces: mean of 32 paths')
    offset = -0.18 if side == 12 else 0.18
    bars = axes[1].bar(np.arange(5)+offset, [r['observed'] for r in rows], width=0.34,
                      color=color, alpha=0.8, label=f'{2*side*side:,} faces')
    axes[1].bar_label(bars, labels=[f"{r['observed']}/32" for r in rows], padding=3, fontsize=9)
axes[0].set(xlabel='First reach of dual-graph distance at least r',
            ylabel='Conditional mean time (attempts / rooted operator)',
            title='Clock-averaged passage times; no exponent fitted', xticks=radii)
axes[0].text(0.03, 0.97, 'Small/large means agree exactly at r = 2, 4, 6\nthrough matched microscopic path prefixes.',
             transform=axes[0].transAxes, va='top', fontsize=9)
axes[1].set(xticks=np.arange(5), xticklabels=radii, ylim=(0, 37),
            xlabel='Target dual-graph distance r', ylabel='Trials reaching r before boundary censoring',
            title='Unresolved small-mesh trials are not discarded')
for ax in axes:
    ax.grid(axis='y', alpha=0.2)
    ax.legend(fontsize=9, loc='lower right' if ax == axes[1] else 'center left')
fig.suptitle('Exact event clock and paired pre-boundary spreading measurements', fontsize=17)
fig.supxlabel('32 independent trial seeds within each size; matching seeds across sizes are paired, not independent replicates.\n'
              'Stop when any enabled primitive touches the periodic seam. Gray-green curves show individual larger-mesh jump paths.', fontsize=10)
output = Path('docs/images/triangle-first-passage.png')
fig.savefig(output, dpi=160)
print(output)
