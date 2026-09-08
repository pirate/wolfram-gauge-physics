#!/usr/bin/env python3
"""Plot measured raw-link paths and the exhaustively observed gauge-state action."""
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Polygon
import numpy as np

from fiber_transport import FiberTransport


def main():
    data = json.loads(Path('data/fiber-transport-order.json').read_text())
    r = data['fibers'][0]
    e = FiberTransport(data['side'], r['fiber_vertices'], r['fiber_edges'])
    colors = ('#237cad', '#ce6120')
    fig, axes = plt.subplots(1, 2, figsize=(12, 6.5))
    fig.subplots_adjust(left=.065, right=.97, bottom=.28, top=.83, wspace=.22)
    ax = axes[0]
    vertices = np.array([(v % data['side'], v//data['side']) for v in range(e.geometry.size)])
    for u, v in e.geometry.edges:
        a, b = vertices[u], vertices[v]
        if np.max(np.abs(a-b)) <= 1:
            ax.plot([a[0], b[0]], [a[1], b[1]], color='#e1e5e9', lw=.6, zorder=0)
    for k, face in enumerate(r['preparation']['positions']):
        ax.add_patch(Polygon(vertices[list(e.geometry.faces[face][:3])], color='#564278', alpha=.8))
        p = np.array(e.point(face))/3
        ax.annotate(('moving defect', 'A encloses this', 'B encloses this', 'spectator')[k],
                    xy=p, xytext=(p[0]+.25, p[1]+.4), fontsize=9,
                    bbox={'facecolor': 'white', 'edgecolor': 'none', 'alpha': .8})
    for label, color, circuit in zip(('A', 'B'), colors, r['circuits']):
        lead, length = circuit['lead_steps'], circuit['cycle_steps']
        walk = np.array([e.point(f) for f in circuit['walk']])/3
        ax.plot(walk[:lead+1, 0], walk[:lead+1, 1], '--', color=color, lw=1.7)
        loop = walk[lead:lead+length+1]
        ax.plot(loop[:, 0], loop[:, 1], '.-', color=color, lw=2, label=f'{label}: 32 primitive updates')
        for k in (1, 5, 9):
            ax.annotate('', xy=loop[k+1], xytext=loop[k],
                        arrowprops={'arrowstyle': '->', 'color': color, 'lw': 2})
    ax.legend(loc='upper left', fontsize=8)
    ax.set(xlim=(.5, 10.5), ylim=(.5, 10.5), aspect='equal',
           xlabel='Supplied mesh x index', ylabel='Supplied mesh y index',
           title='Actual closed paths; every defect returns')

    ax = axes[1]
    orbit = r['controlled_circuit_orbit']
    # Layout of the measured state-transition graph, not a spatial embedding.
    positions = np.array([[-1, 1], [1, 1], [-1, -1], [1, -1]])
    for action, color in zip(orbit['generator_permutations'], colors):
        for source, target in enumerate(action):
            if source == target:
                continue
            ax.add_patch(FancyArrowPatch(positions[source], positions[target],
                         connectionstyle='arc3,rad=.13', arrowstyle='-|>',
                         mutation_scale=18, shrinkA=31, shrinkB=31, lw=2.2, color=color, zorder=1))
    for state, (x, y) in enumerate(positions):
        ax.scatter([x], [y], s=1650, color='#f4f0fa', edgecolor='#564278', lw=1.7, zorder=2)
        p = r['loop_action_audit']['projective_points_by_state'][state]
        ax.text(x, y, f'{state}\n[{p[0]} : {p[1]}]', ha='center', va='center', fontsize=10, zorder=3)
    fixed = [[i for i, j in enumerate(p) if i == j] for p in orbit['generator_permutations']]
    ax.text(.5, -.03, f'A fixes {fixed[0][0]}; B fixes {fixed[1][0]} (self-loops omitted)',
            transform=ax.transAxes, ha='center', fontsize=9)
    ax.set(xlim=(-1.6, 1.6), ylim=(-1.55, 1.55), aspect='equal',
           title=r'Four complete gauge states; induced action $A_4$')
    ax.axis('off')
    ax.plot([], [], color=colors[0], label='Circuit A')
    ax.plot([], [], color=colors[1], label='Circuit B')
    ax.legend(loc='center', fontsize=9, framealpha=.95)
    fig.suptitle('Order-dependent holonomy memory from a triangle fiber', fontsize=16, y=.97)
    fig.text(.5, .9, r'Derived microscopic group: $\mathrm{Aut}(C_3)=S_3$'
             '  |  Controlled classical transport, not quantum amplitudes', ha='center', fontsize=10)
    fig.text(.5, .105, r'A then B and B then A: identical face classes, inequivalent connections'
             '\nFull 432-vertex graph: first changed Laplacian trace at power 12; difference = 360',
             ha='center', fontsize=10, linespacing=1.6)
    fig.text(.5, .025, 'Controls: square fiber → 2 states / commuting action; commuting triangle seed → 1 state / trivial action',
             ha='center', fontsize=9, color='#444')
    fig.savefig('docs/images/fiber-transport-order.png', dpi=180)


if __name__ == '__main__':
    main()
