#!/usr/bin/env python3
"""Local raw-link witness and exact disk-cylinder reference probabilities.

These calculations concern the existing operator bank. A cylinder is an
event in its stationary reference, not a boundary condition imposed on
evolution. Root-fixed gauge orbits retain all relative boundary alignment.
"""
import json
import math
import random
from collections import deque
from fractions import Fraction
from functools import lru_cache

from triangle_gauge_current_corrector import GaugeCurrent
from triangle_reference import polynomial_power, subgroup


@lru_cache(maxsize=None)
def character_coefficients(faces):
    return tuple(tuple(polynomial_power(base, faces))
                 for base in ((1, 3, 2), (1, -3, 2), (1, 0, -1)))


def disk_orbit_probability(total_faces, total_charge, disk_faces,
                           disk_charge, boundary_standard_character):
    """Exact all-sector probability of one root-fixed disk gauge orbit.

    Sum over the complement's two handle elements, keeping the disk's
    boundary product. The characters are 1, sign, and fixed-points minus 1.
    The free vertex-frame factor cancels against the disk orbit size.
    """
    remaining = total_faces-disk_faces
    assert remaining >= 1 and total_charge % 2 == 0
    q = total_charge-disk_charge
    if not 0 <= q <= 2*remaining:
        return Fraction()
    outside = character_coefficients(remaining)
    whole = character_coefficients(total_faces)
    numerator = (outside[0][q]+(-1)**disk_charge*outside[1][q]
                 +Fraction(boundary_standard_character, 2)*outside[2][q])
    denominator = sum(p[total_charge] for p in whole)
    return numerator/denominator


def local_witness():
    side = 8
    g = GaugeCurrent(side, calculate_moments=False)
    x, oracle = g.x, g.x.e.factor.oracle
    geometry, group = x.e.geometry, x.e.geometry.group
    patch = next(p for p, pair in enumerate(oracle.pairs)
                 if geometry.edges[pair[0][0][0]] == (27, 28))
    reads = {e for path in oracle.pairs[patch] for e, _ in path}
    for fan in g.affected_fans[patch]:
        reads.update(oracle.fan.reads[fan])
        faces = set(x.p.supports[1][fan])
        faces.update(f for assignment in g.kernels[fan] for f, _ in assignment)
        reads.update(e for f in faces for e, _ in oracle.fan.face_paths[f])
    vertices = {v for e in reads for v in geometry.edges[e]}
    left, right = min(v % side for v in vertices), max(v % side for v in vertices)
    bottom, top = min(v//side for v in vertices), max(v//side for v in vertices)
    disk_faces = [2*(y*side+x)+t for y in range(bottom, top)
                  for x in range(left, right) for t in (0, 1)]
    disk_edges = {e for f in disk_faces for e, _ in oracle.fan.face_paths[f]}
    disk_vertices = {v for e in disk_edges for v in geometry.edges[e]}
    assert reads <= disk_edges
    assert len(disk_vertices)-len(disk_edges)+len(disk_faces) == 1
    rng = random.Random(150918001)
    for draw in range(256):
        raw = [rng.randrange(group.n) for _ in geometry.edges]
        before = g.activities(raw)
        weights = g.weights(x.p.charge(raw))
        moved = raw[:]
        old, new = x.apply(moved, x.supports+patch)
        after = g.activities(moved)
        delta = sum(k*(b-a) for k, a, b in zip(weights, before, after))
        if delta:
            break
    else:
        raise ValueError('No local current witness found')
    assert g.full_generator_current_drift(moved)-g.full_generator_current_drift(raw) == delta
    # Rooted parallel transports on the disk, without identifying other patches.
    adjacent = {v: [] for v in disk_vertices}
    for edge in disk_edges:
        a, b = geometry.edges[edge]
        adjacent[a].append((b, raw[edge]))
        adjacent[b].append((a, group.inv[raw[edge]]))
    root = min(disk_vertices)
    frames, queue = {root: group.identity}, deque([root])
    while queue:
        a = queue.popleft()
        for b, transport in sorted(adjacent[a]):
            if b not in frames:
                frames[b] = group.mul[transport][frames[a]]
                queue.append(b)
    based = [group.mul[group.inv[frames[b]]][group.mul[raw[e]][frames[a]]]
             for e in disk_edges for a, b in [geometry.edges[e]]]
    charges = x.p.charge(raw)
    assert len(subgroup(group, based)) == 6
    assert any(charges[f] == 1 for f in disk_faces)
    boundary = ([bottom*side+x for x in range(left, right)]
                +[y*side+right for y in range(bottom, top)]
                +[top*side+x for x in range(right, left, -1)]
                +[y*side+left for y in range(top, bottom, -1)])
    paths = [(geometry.ids[tuple(sorted((a, b)))], a > b)
             for a, b in zip(boundary, boundary[1:]+boundary[:1])]
    holonomy = oracle.fan.transport(raw, paths)
    standard = sum(i == v for i, v in enumerate(group.elements[holonomy]))-1
    charge = sum(charges[f] for f in disk_faces)
    z = 1/math.sqrt(2)
    p0 = 1+3*z+2*z*z
    limiting_probability = z**charge/p0**len(disk_faces)
    return {'witness_seed': 150918001, 'draw_index': draw, 'side': side,
            'elastic_patch': patch, 'elastic_rule': 1, 'pair_before': old, 'pair_after': new,
            'delta_L_b3': delta, 'read_edges': len(reads),
            'rectangle': [left, right, bottom, top], 'disk_faces': len(disk_faces),
            'disk_vertices': len(disk_vertices), 'disk_edges': len(disk_edges),
            'disk_charge': charge, 'disk_loop_group_order': 6,
            'boundary_holonomy': holonomy, 'boundary_standard_character': standard,
            'disk_raw_links': [[list(geometry.edges[e]), raw[e]] for e in sorted(disk_edges)],
            'density_one_limiting_orbit_probability': limiting_probability,
            'density_one_S_per_face_liminf_lower_bound': delta**2*limiting_probability/4,
            'finite_all_sector_orbit_probabilities': [
                {'faces': f, 'probability': float(disk_orbit_probability(
                    f, f, len(disk_faces), charge, standard))}
                for f in (128, 288, 512, 1152)],
            'scope': 'Positive-probability local witness for a bulk variational correction, not its useful numerical magnitude or a hydrodynamic limit.'}


if __name__ == '__main__':
    print(json.dumps(local_witness()), flush=True)
