#!/usr/bin/env python3
"""Plot recorded conserved-sector trajectories alongside isolated controls."""
import argparse
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument('--input', type=Path, default=Path('data/d4-braid-rule-screen.json'))
parser.add_argument('--output', type=Path, default=Path('docs/images/derived-charge-collisions.png'))
args = parser.parse_args()
data = json.loads(args.input.read_text())
fig, axes = plt.subplots(1, 2, figsize=(10, 4.6), constrained_layout=True)
for axis, identifier, title in zip(axes, (11229, 17581), ('Mobile sectors reflect', 'Mobile sector displaces a stationary sector')):
    rule = next(r for r in data['selected_rules'] if r['solution_id'] == identifier)
    collision = next(c for c in rule['collisions'] if c['seeds'] == [1, 2])
    rows = collision['trajectory']
    for index, color in ((1, '#1768a4'), (2, '#c35e27')):
        axis.plot([r[index+2]-224 for r in rows], [r[0] for r in rows], '--', color=color, alpha=.4)
        axis.plot([r[index]-224 for r in rows], [r[0] for r in rows], color=color,
                  label=f'Sector {collision["seeds"][index-1]}', linewidth=2)
    axis.set(title=f'{title}\nFinite-table rule {identifier}', xlabel='Cell offset from first seed', ylabel='Update layers')
    axis.grid(alpha=.2)
    axis.legend(frameon=False)
fig.suptitle('Conserved defects in rules selected by algebraic constraints', fontsize=14)
fig.supxlabel('Solid: two-defect evolution. Dashed: isolated controls. Fixed chain and specified schedule; no physical energy or time calibration.', fontsize=9)
args.output.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(args.output, dpi=180)
print(args.output)
