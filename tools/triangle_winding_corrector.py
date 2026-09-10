#!/usr/bin/env python3
"""Exact one-function variational correction to conserved winding transport.

Current uses three times the supplied cell-coordinate displacement to keep
the entire moment census integral. No diffusion law is used in the dynamics.
"""
import itertools
import json
from collections import Counter
from fractions import Fraction
from pathlib import Path

import numpy as np

from triangle_charge_response import ChargeMarginal, coefficients
from triangle_elastic_scattering import ElasticExperiment


def add_term(poly, faces, values, coefficient):
    if coefficient:
        poly[tuple(sorted(zip(faces, values)))] += int(coefficient)


def squared_counts(poly, fixed=()):
    """Group an indicator polynomial's squared expectation by charge counts."""
    terms = [(key, value) for key, value in poly.items() if value]
    result = Counter()
    for i, (a, x) in enumerate(terms):
        for b, y in terms[i:]:
            merged = dict(fixed)
            compatible = True
            for face, value in (*a, *b):
                if face in merged and merged[face] != value:
                    compatible = False
                    break
                merged[face] = value
            if compatible:
                counts = tuple(sum(q == k for q in merged.values()) for k in range(3))
                result[counts] += x*y*(1 if a == b else 2)
    return Counter({key: value for key, value in result.items() if value})


def difference(poly, faces, before, after):
    old, new = dict(zip(faces, before)), dict(zip(faces, after))
    result = Counter()
    for assignment, coefficient in poly.items():
        if not coefficient or not any(f in old for f, _ in assignment):
            continue
        outside = tuple((f, q) for f, q in assignment if f not in old)
        if all(f not in old or old[f] == q for f, q in assignment):
            result[outside] -= coefficient
        if all(f not in new or new[f] == q for f, q in assignment):
            result[outside] += coefficient
    return Counter({key: value for key, value in result.items() if value})


def census(side, calculate_moments=True):
    bank = json.loads(Path('data/triangle-feedback.json').read_text())['bank']
    elastic = json.loads(Path('data/triangle-elastic-scattering.json').read_text())['pair_census']['elastic_tables']
    x = ElasticExperiment(side, bank, elastic)
    f = x.geometry.faces
    coordinates = np.array([[3*(i//2 % side)+(2 if i % 2 == 0 else 1),
                             3*(i//2//side)+(1 if i % 2 == 0 else 2)] for i in range(f)])
    period = 3*side
    displacements = np.array([((coordinates[b]-coordinates[a]+period//2) % period)-period//2
                              for a, b in x.p.dual])
    assert set(map(tuple, displacements)) <= {(a, b) for a, b in
        ((-1, 1), (-1, -2), (2, 1), (1, -1), (1, 2), (-2, -1))}
    poly = Counter()
    shot_counts = Counter()
    for arity in (0, 1):
        for faces, path in zip(x.p.supports[arity], x.p.paths[arity]):
            steps = [int(displacements[edge, 0])*sign for edge, sign in path]
            values = ((0, 1), (1, 0), (0, 2), (2, 0)) if arity == 0 else itertools.permutations((0, 1, 2))
            for q in values:
                current = (steps[0]*(q[0]-q[1]) if arity == 0 else
                           4*(steps[0]*(q[0]-1)+steps[1]*(1-q[2])))
                add_term(poly, faces, q, current)
                counts = tuple(q.count(k) for k in range(3))
                shot_counts[counts] += current**2 if arity == 0 else current**2//2
    assert shot_counts == {(1, 1, 0): 12*f, (1, 0, 1): 48*f, (1, 1, 1): 480*f}
    poly = Counter({key: value for key, value in poly.items() if value})
    if not calculate_moments:
        return x, poly, None, None
    variance = squared_counts(poly)
    twice_energy = Counter()
    for arity in (0, 1):
        for faces in x.p.supports[arity]:
            values = ((0, 1), (1, 0), (0, 2), (2, 0)) if arity == 0 else itertools.permutations((0, 1, 2))
            for old in values:
                new = old[::-1] if arity == 0 else (1, 1, 1)
                delta = difference(poly, faces, old, new)
                # Reversibility pairs every gauge-dependent birth with its
                # charge-only reverse death. This includes the full S3-sector
                # multiplicities without assuming a conditional gate average.
                weight = 1 if arity == 0 else 8
                for counts, coefficient in squared_counts(delta, tuple(zip(faces, old))).items():
                    twice_energy[counts] += weight*coefficient
    assert all(value % f == 0 for value in (*variance.values(), *twice_energy.values()))
    return x, poly, variance, twice_energy


def evaluate(x, variance, twice_energy, charge):
    f = x.geometry.faces
    marginal = ChargeMarginal(f, charge, 'nonabelian_reflections')
    v = marginal.expectation(variance.items())
    energy = marginal.expectation(twice_energy.items())/2
    a, b, chi, _, _ = coefficients(marginal)
    shot = 6*f*(a+10*b)
    bare = shot/(810*f)
    correction = 2*v*v/(810*f*energy) if energy else Fraction()
    assert v >= 0 and energy >= 0 and bare >= correction
    rational = lambda value: [value.numerator, value.denominator]
    return {'faces': f, 'charge': charge,
            'current_drift_variance': rational(v), 'current_drift_dirichlet_energy': rational(energy),
            'bare_xx_mobility': float(bare),
            'one_function_xx_mobility_upper_bound': float(bare-correction),
            'minimum_relative_reduction': float(correction/bare),
            'optimal_trial_coefficient': rational(v/energy) if energy else [0, 1],
            'static_susceptibility': float(chi),
            'scope': 'Finite-volume stationary winding mobility bound, not an evaluated exact mobility or a proved hydrodynamic diffusion coefficient.'}


def bulk_bound(z):
    """Analytic fixed-density bound; z is a fugacity, not a model parameter."""
    partition = (1+z)*(1+2*z)
    quartic = 68*z**4-108*z**3+1013*z**2-54*z+17
    cubic = 16*z**3+270*z**2+17*z+3
    eighth = (544*z**8+1560*z**7+11712*z**6+53070*z**5+65314*z**4
               +27459*z**3+3495*z**2+306*z+34)
    variance = 144*z**3*quartic/partition**5
    energy = 1152*z**3*eighth/partition**7
    bare = 2*z*cubic/(135*partition**3)
    reduction = 3*z*z*quartic*quartic/(cubic*eighth)
    assert 0 < reduction < 1
    return {'density': float(z*(4*z+3)/partition), 'bare_xx_mobility': float(bare),
            'xx_mobility_limsup_upper_bound': float(bare*(1-reduction)),
            'minimum_relative_reduction': float(reduction),
            'current_drift_variance_per_face': float(variance),
            'current_drift_dirichlet_energy_per_face': float(energy)}


if __name__ == '__main__':
    for side in (3, 6):
        x, poly, variance, twice_energy = census(side)
        print(json.dumps({'side': side, 'drift_indicator_terms': len(poly),
                          'variance_coefficients_per_face': [[list(c), value//x.geometry.faces] for c, value in sorted(variance.items())],
                          'twice_dirichlet_coefficients_per_face': [[list(c), value//x.geometry.faces] for c, value in sorted(twice_energy.items()) if value],
                          'results': [evaluate(x, variance, twice_energy, charge) for charge in (2, 4, x.geometry.faces)]}), flush=True)
    print(json.dumps({'fixed_density_bounds': [bulk_bound(Fraction(z, 4)) for z in (1, 2, 4, 8, 16)]}), flush=True)
