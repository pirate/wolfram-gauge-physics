#!/usr/bin/env python3
"""Actual face-class snapshots on the supplied periodic triangulation."""
import json
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, Patch
from face_energy_obstruction import derive_fiber_group
from run_shared_edge import mesh_geometry


def main():
    data = json.loads(Path('data/d4-shared-edge.json').read_text())
    group = derive_fiber_group()
    side = data['side']
    edges, faces, _ = mesh_geometry(side)
    colors = {0: '#f6f7f8', 1: '#db7530', 2: '#148b7a', 3: '#327bc0', 5: '#c6a6d5'}
    names = {0: 'Flat', 1: 'Reflection class 1', 2: 'Reflection class 2', 3: 'Rotation class', 5: 'Central half-turn'}
    fig, axes = plt.subplots(2, 2, figsize=(10, 10), facecolor='white')
    for row, condition in enumerate(('reflection', 'rotation')):
        run = next(r for r in data['runs'] if r['condition'] == condition and r['trial'] == 0)
        for column, key in enumerate(('initial_links', 'final_links')):
            values = run[key]
            sectors = []
            for face in faces:
                h = group.identity
                for u, v in zip(face, face[1:]):
                    value = values[edges.index(tuple(sorted((u, v))))]
                    h = group.mul[group.inv[value] if u > v else value][h]
                sectors.append(group.sectors[h])
            ax = axes[row, column]
            for y in range(side):
                for x in range(side):
                    triangles = (((x, y), (x+1, y), (x+1, y+1)),
                                 ((x, y), (x+1, y+1), (x, y+1)))
                    for half, triangle in enumerate(triangles):
                        sector = sectors[2*(y*side+x)+half]
                        ax.add_patch(Polygon(triangle, facecolor=colors[sector], edgecolor='#ced3d8', linewidth=.45))
            tick = 0 if column == 0 else data['layers']
            defect_count = sum(s not in (0, 5) for s in sectors)
            if defect_count != sum(run['derived_charge_totals']):
                raise ValueError('snapshot charge disagrees with the measured invariant')
            ax.set(title=f'{condition.capitalize()} seed · layer {tick}\n'
                         f'{defect_count} conserved defects; {sum(s != 0 for s in sectors)} nonflat faces',
                   xlim=(0, side), ylim=(0, side), aspect='equal', xticks=[], yticks=[])
    fig.suptitle('Conserved defects, different curvature trails', fontsize=18, y=.985)
    fig.legend(handles=[Patch(facecolor=colors[s], edgecolor='#aaa', label=names[s]) for s in (0, 1, 2, 3, 5)],
               loc='lower center', ncol=3, bbox_to_anchor=(.5, .04), frameon=False)
    fig.text(.5, .012, f'Actual link-simulation output on a supplied periodic {side}×{side} triangulation. Colors are gauge classes, not energies.',
             ha='center', fontsize=9)
    fig.tight_layout(rect=(0, .11, 1, .96))
    fig.savefig('docs/images/shared-edge-trails.png', dpi=160)


if __name__ == '__main__':
    main()
