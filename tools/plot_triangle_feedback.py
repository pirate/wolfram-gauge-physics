#!/usr/bin/env python3
"""Draw actual before/after fan faces and measured full-generator reaction counts."""
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
import numpy as np

from triangle_feedback import FeedbackExperiment


def main():
    data = json.loads(Path('data/triangle-feedback.json').read_text())
    r = data['readout']; e = FeedbackExperiment(12, data['bank'])
    active, inactive = r['global_rate_witness']['states']
    before = r['preparations'][active]['checked']['final_links']
    after = e.compiled_bank([before], r['readout_schedule'][:1], 1)['runs'][1]['final_links']
    changed = {edge for edge, a, b in zip(e.geometry.edges, before, after) if a != b}
    points = np.array([(v % 12, v//12) for v in range(e.geometry.size)])
    faces = r['encounter_faces']; vertices = sorted({v for f in faces for v in e.geometry.faces[f][:3]})
    edges = {tuple(sorted(pair)) for f in faces for pair in zip(e.geometry.faces[f], e.geometry.faces[f][1:])}
    root = e.factor.oracle.fan.specs[r['encounter_patch']][0][0]
    colors = {0: '#f0f1f3', 1: '#9c85bf', 3: '#e7a141'}
    names = {0: 'vacancy\nq = 0', 1: 'reflection\nq = 1', 3: '3-cycle\nq = 2'}
    fig, axes = plt.subplots(1, 3, figsize=(13, 5.2), gridspec_kw={'width_ratios': [1, 1, 1.25]})
    fig.subplots_adjust(top=.79, bottom=.32, left=.03, right=.975, wspace=.3)
    for ax, links, title in zip(axes[:2], (before, after), ('Before the local reaction', 'After one primitive two-spoke update')):
        classes = e.geometry.sectors(links)
        for face in faces:
            xy = points[list(e.geometry.faces[face][:3])]
            ax.add_patch(Polygon(xy, facecolor=colors[classes[face]], alpha=.8, edgecolor='none'))
            center = xy.mean(axis=0)
            ax.text(*center, names[classes[face]], ha='center', va='center', fontsize=8)
        for edge in edges:
            xy = points[list(edge)]
            highlight = links is after and edge in changed
            ax.plot(xy[:, 0], xy[:, 1], color='#c45513' if highlight else '#53606a', lw=3 if highlight else 1.3)
        ax.scatter(points[vertices, 0], points[vertices, 1], s=23, c='#28353e', zorder=3)
        ax.scatter([points[root, 0]], [points[root, 1]], s=95, facecolor='white', edgecolor='#28353e', zorder=4)
        ax.set(aspect='equal', xlim=(4.75, 7.25), ylim=(4.75, 7.25), title=title)
        ax.axis('off')
        ax.text(.5, -.08, 'Actual supplied mesh patch; local charge = 3', transform=ax.transAxes, ha='center', fontsize=8)
    ax = axes[2]
    counts = [row['class_changing_triple_operators'] for row in r['global_class_rates']]
    bars = ax.bar(range(4), counts, color=['#8c70b0' if x else '#bbc2c8' for x in counts], width=.65)
    for bar, count in zip(bars, counts):
        ax.text(bar.get_x()+bar.get_width()/2, count+.35, str(count), ha='center', fontsize=11)
    ax.set(xticks=range(4), ylim=(0, 15), yticks=(0, 4, 8, 12), xlabel='Stored braid-memory state',
           ylabel='Enabled class-changing triple operators', title='Same face classes; different reaction rates')
    ax.spines[['right', 'top']].set_visible(False)
    ax.text(.5, -.28, f'States {active} and {inactive} also have identical\nencounter-boundary holonomy', transform=ax.transAxes, ha='center', fontsize=9)
    fig.suptitle('Classical braid memory can gate a local reaction', fontsize=16, y=.98)
    fig.text(.5, .87, 'Triangle adjacency → derived automorphisms → exhaustive rule census → actual link evolution', ha='center', fontsize=11)
    fig.text(.5, .075, 'All 11,232 rooted operator choices checked in C++ and Python; orange spokes are the changed links.', ha='center', fontsize=10)
    fig.text(.5, .025, 'Positive conserved q is a derived combinatorial charge—not physical energy, particle mass, or a quantum amplitude.', ha='center', fontsize=9, color='#444')
    fig.savefig('docs/images/triangle-feedback.png', dpi=180)


if __name__ == '__main__':
    main()
