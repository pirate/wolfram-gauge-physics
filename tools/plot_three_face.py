#!/usr/bin/env python3
"""Render the actual fixed-exterior witness, including commonly based holonomies."""
import json
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from face_energy_obstruction import derive_fiber_group
from run_three_face import LinkOracle


def main():
    data = json.loads(Path('data/d4-three-face.json').read_text())
    group = derive_fiber_group()
    oracle = LinkOracle(data['side'], group)
    witness = data['analysis']['fixed_exterior_witness']
    p = witness['prefix_patch']
    faces, paths = oracle.specs[p], oracle.patches[p]
    root = faces[0][0]
    nodes = sorted({v for face in faces for v in face})
    positions = {v: (v % data['side'], v//data['side']) for v in nodes}
    colors = {0: '#f5f6f7', 1: '#e79c52', 5: '#b297d0'}
    labels = {0: '1', 1: 'r', 4: 'rz', 5: 'z'}
    fig, axes = plt.subplots(2, 2, figsize=(10, 8))
    for row, prefix in enumerate(('initial', 'after')):
        for column, side in enumerate(('left', 'right')):
            ax = axes[row, column]
            values = witness[prefix+'_'+side]
            weights = []
            for face, path in zip(faces, paths):
                h = oracle.transport(values, path)
                q = data['analysis']['positive_weight_per_element'][h]
                weights.append(q)
                points = [positions[v] for v in face[:-1]]
                ax.add_patch(Polygon(points, facecolor=colors[group.sectors[h]], edgecolor='#354657', linewidth=1.8))
                x, y = (sum(point[i] for point in points)/3 for i in (0, 1))
                ax.text(x, y, f'{labels[h]}\nw={q}', ha='center', va='center', fontsize=13)
            for v, (x, y) in positions.items():
                ax.scatter(x, y, s=40 if v == root else 15, c='#12263a', zorder=3)
            x, y = positions[root]
            ax.annotate('common root', (x, y), xytext=(-45, 24), textcoords='offset points',
                        ha='right', va='bottom', fontsize=10,
                        arrowprops={'arrowstyle': '-', 'color': '#46596a', 'lw': .8})
            ax.set_title(('Before' if row == 0 else 'After the same update')+f' · connection {column+1}\nfan weight = {sum(weights)}', fontsize=13)
            ax.set_aspect('equal'); ax.autoscale_view(); ax.margins(.2); ax.axis('off')
    fig.suptitle('Same face classes and exterior links; different conserved-density response', fontsize=15, y=.98)
    fig.text(.5, .065, 'Orange: the same reflection class (r and rz). Purple: central curvature z. White: flat.\n'
             'Letter labels are holonomies in one common frame; colors and weights are gauge invariant.', ha='center', fontsize=10)
    fig.text(.5, .018, 'Actual selected fan from a supplied triangulation. The weight is derived from the rule, not a calibrated physical energy.',
             ha='center', fontsize=9)
    fig.tight_layout(rect=(0, .12, 1, .94))
    fig.savefig('docs/images/three-face-feedback.png', dpi=170)


if __name__ == '__main__':
    main()
