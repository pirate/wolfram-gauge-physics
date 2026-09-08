#!/usr/bin/env python3
"""Exact stationary charge projection and its two-step closure defect.

No heat equation is supplied to evolution. These are response observables of
the original uniform-attempt Markov operator on actual finite-group connections.
"""
import argparse
import itertools
import json
from collections import Counter
from fractions import Fraction
from functools import lru_cache
from pathlib import Path

import numpy as np

from triangle_reference import reference, fraction
from triangle_charge_current import CurrentProbe


def falling(n, k):
    value = 1
    for i in range(k):
        value *= n-i
    return value if n >= k else 0


class ChargeMarginal:
    def __init__(self, faces, charge, sector):
        self.reference = reference(faces, charge, sector)
        self.faces = faces

    @lru_cache(maxsize=None)
    def probability(self, counts):
        """Probability of ONE ordered assignment with these charge counts."""
        n = sum(counts)
        r = self.reference
        numerator = sum(row['canonical_connections']*
                        falling(row['populations'][0], counts[0])*
                        falling(row['populations'][1], counts[1])*
                        falling(row['populations'][2], counts[2]) for row in r['population_counts'])
        return Fraction(numerator, r['canonical_connections']*falling(self.faces, n))

    def expectation(self, polynomial):
        return sum(coefficient*self.probability(counts) for counts, coefficient in polynomial)


def term_value(term, q):
    indices, coefficient = term
    values = [q[i] for i in indices]
    if len(values) == 2:
        a, b = values
        return coefficient*(b-a) if (a == 0) != (b == 0) else 0
    return coefficient*(1-values[0]) if sorted(values) == [0, 1, 2] else 0


@lru_cache(maxsize=None)
def product_polynomial(first, second):
    """Integer coefficients grouped by charge populations on the union support."""
    size = max(*first[0], *second[0])+1
    result = Counter()
    for q in itertools.product(range(3), repeat=size):
        value = term_value(first, q)*term_value(second, q)
        if value:
            result[tuple(q.count(i) for i in range(3))] += value
    return tuple(sorted((key, value) for key, value in result.items() if value))


def canonical_terms(first, second):
    order = list(dict.fromkeys(first[0]+second[0]))
    ids = {face: i for i, face in enumerate(order)}
    return tuple((tuple(ids[f] for f in term[0]), term[1]) for term in (first, second))


class ResponseGeometry:
    def __init__(self, side, bank):
        self.p = CurrentProbe(side, bank)
        self.e = self.p.e
        self.side, self.faces = side, self.p.faces
        n = self.faces
        adjacency = np.zeros((n, n), dtype=np.int64)
        for a, b in self.p.dual:
            adjacency[a, b] = adjacency[b, a] = 1
        if not np.all(adjacency.sum(axis=1) == 3):
            raise ValueError('response requires the regular honeycomb dual')
        self.l1 = 3*np.eye(n, dtype=np.int64)-adjacency
        a2 = adjacency@adjacency-3*np.eye(n, dtype=np.int64)
        if not np.all((a2 == 0) | (a2 == 1)) or not np.all(a2.sum(axis=1) == 6):
            raise ValueError('distance-two dual relation is not simple six-regular')
        self.l2 = 6*np.eye(n, dtype=np.int64)-a2
        if not np.array_equal(self.l2, 6*self.l1-self.l1@self.l1):
            raise ValueError('dual distance-two Laplacian identity fails')
        self.squares = (self.l1@self.l1, self.l1@self.l2, self.l2@self.l2)
        fan_sum = np.zeros((n, n), dtype=np.int64)
        self.terms = []
        for f in range(n):
            terms = Counter()
            for a, b in self.e.factor.pairs:
                if f == a or f == b:
                    other = b if f == a else a
                    terms[(f, other)] += 1
            for faces in self.e.fan_faces:
                if f in faces:
                    terms[(f, *sorted(set(faces)-{f}))] += 4
            self.terms.append(tuple(sorted(terms.items())))
        for faces in self.e.fan_faces:
            for i, a in enumerate(faces):
                for j, b in enumerate(faces):
                    fan_sum[a, b] += 3*int(i == j)-1
        if not np.array_equal(fan_sum, 4*self.l1+self.l2):
            raise ValueError('actual fan incidences fail the projected stencil identity')
        # Translation prototypes for the two face orientations. Disjoint terms
        # have exactly zero covariance by exchangeability and zero conditional
        # drift averaged over the permutations of each local charge population.
        self.polynomials = []
        for root in (0, 1):
            row = []
            for other in range(n):
                total = Counter()
                for first in self.terms[root]:
                    for second in self.terms[other]:
                        if set(first[0]).isdisjoint(second[0]):
                            continue
                        for key, value in product_polynomial(*canonical_terms(first, second)):
                            total[key] += value
                row.append(tuple(sorted((key, value) for key, value in total.items() if value)))
            self.polynomials.append(row)

    def drift(self, q):
        return [sum(term_value(term, q) for term in terms) for terms in self.terms]

    def covariance(self, marginal):
        prototypes = [[marginal.expectation(poly) for poly in row] for row in self.polynomials]
        matrix = [[Fraction() for _ in range(self.faces)] for _ in range(self.faces)]
        for y in range(self.side):
            for x in range(self.side):
                for parity in (0, 1):
                    f = 2*(y*self.side+x)+parity
                    for other, value in enumerate(prototypes[parity]):
                        cell, orientation = divmod(other, 2)
                        yy, xx = divmod(cell, self.side)
                        target = 2*(((yy+y) % self.side)*self.side+(xx+x) % self.side)+orientation
                        matrix[f][target] = value
        if any(matrix[i][j] != matrix[j][i] for i in range(self.faces) for j in range(self.faces)) or any(sum(row) for row in matrix):
            raise ValueError('exact drift covariance fails symmetry or charge conservation')
        return matrix


def coefficients(marginal):
    f = marginal.faces
    r = marginal.reference
    # Ordered pair probabilities; both vacancy positions are included.
    a = 2*(marginal.probability((1, 1, 0))+4*marginal.probability((1, 0, 1)))
    p012 = 6*marginal.probability((1, 1, 1))
    b = Fraction(4, 3)*p012
    chi = Fraction(*r['single_face_charge_variance'])*f/(f-1)
    if chi <= 0:
        raise ValueError('charge response requires nonzero reference susceptibility')
    return a, b, chi, (a+4*b)/chi, b/chi


def sparse_formula(f):
    """Closed Q=4 formulas, independently compared with population sums."""
    if type(f) is not int or f < 18:
        raise ValueError('sparse response formula requires at least 18 faces')
    p = 20*f*f-59*f+24
    common = f*(f-1)*(10*f-3)
    return {'chi': Fraction(2*p, common),
            'a': Fraction(4*(f-3)*(20*f+1), common),
            'b': Fraction(432*(f-3), common*(f-2)),
            'h': Fraction(48*(40*f**3-198*f*f+3659*f-7575), common*(f-2)),
            'residual': Fraction(1296*(820*f**4-4479*f**3-16993*f*f+146682*f-229896), common*(f-2)**2*p)}


def sparse_identity():
    # Integer polynomial arithmetic clears every denominator. No interpolation
    # or floating fit is used to certify the rational response formula.
    from numpy.polynomial import Polynomial
    x = Polynomial(np.array([0, 1], dtype=object))
    a = 4*(x-3)*(20*x+1)*(x-2)
    b = 432*(x-3)
    c = 2*(20*x*x-59*x+24)*(x-2)
    h = 48*(40*x**3-198*x*x+3659*x-7575)
    numerator = 820*x**4-4479*x**3-16993*x*x+146682*x-229896
    residual = h*c-12*a*a-132*a*b-378*b*b-2592*numerator
    if any(residual.coef):
        raise ValueError('cleared sparse closure-defect identity is not zero')
    shifted = numerator(x+18)
    if not shifted.coef[0] > 0 or any(v < 0 for v in shifted.coef):
        raise ValueError('sparse defect positivity is not certified above 18 faces')
    denominator = x*(x-1)*(x-2)*(10*x-3)
    def leading_ratio(top, bottom, power):
        top, bottom = top.trim(), bottom.trim()
        if bottom.degree()-top.degree() != power:
            raise ValueError('claimed sparse asymptotic power disagrees with polynomial degrees')
        return fraction(Fraction(int(top.coef[-1]), int(bottom.coef[-1])))
    return {'cleared_identity_coefficients': [int(v) for v in residual.coef],
            'defect_numerator_in_F_minus_18': [int(v) for v in shifted.coef],
            'limit_F_cubed_times_defect': leading_ratio(2592*numerator, denominator*c, 3),
            'limit_F_squared_times_unresolved_fraction': leading_ratio(2592*numerator, c*h, 2),
            'limit_F_squared_times_nearest_rate_minus_2': leading_ratio(a+4*b-2*c, c, 2),
            'limit_F_squared_times_next_nearest_rate': leading_ratio(b, c, 2)}


def result(geometry, charge, sector='nonabelian_reflections'):
    marginal = ChargeMarginal(geometry.faces, charge, sector)
    a, b, chi, k1, k2 = coefficients(marginal)
    # The commuting integer stencil products avoid dense rational matrix
    # multiplication. Translation leaves only two orientation rows to store.
    h = [[marginal.expectation(poly) for poly in row] for row in geometry.polynomials]
    weights = (chi*k1*k1, 2*chi*k1*k2, chi*k2*k2)
    residual = [[h[i][j]-sum(w*int(matrix[i, j]) for w, matrix in zip(weights, geometry.squares))
                 for j in range(geometry.faces)] for i in (0, 1)]
    if h[0][0] != h[1][1] or residual[0][0] != residual[1][1] or any(sum(row) for row in h+residual):
        raise ValueError('response prototypes fail orientation symmetry or conservation')
    defect = residual[0][0]
    if defect < 0:
        raise ValueError('projection residual has negative squared norm')
    m = 39*geometry.faces
    if charge == 4 and sector == 'nonabelian_reflections':
        formula = sparse_formula(geometry.faces)
        if (chi, a, b, h[0][0], defect) != tuple(formula[k] for k in ('chi', 'a', 'b', 'h', 'residual')):
            raise ValueError('closed sparse formulas disagree with exact population/stencil sums')
    sparse = lambda rows: [[[j, *fraction(value)] for j, value in enumerate(row) if value] for row in rows]
    return {'side': geometry.side, 'faces': geometry.faces, 'charge': charge, 'sector': sector,
            'susceptibility': fraction(chi), 'pair_dirichlet_weight': fraction(a),
            'fan_dirichlet_weight': fraction(b), 'nearest_rate': fraction(k1), 'next_nearest_rate': fraction(k2),
            'long_wavelength_quadratic_coefficient': fraction((k1+6*k2)/3),
            'drift_covariance_orientation_rows': sparse(h[:2]),
            'closure_residual_covariance_orientation_rows': sparse(residual[:2]),
            'drift_squared_norm_per_face': fraction(h[0][0]),
            'closure_residual_squared_norm_per_face': fraction(defect),
            'unresolved_drift_variance_fraction': fraction(defect/h[0][0]) if h[0][0] else None,
            'two_step_trace_excess_per_face': fraction(defect/(m*m)),
            'exact_linear_mean_closure': defect == 0}


def exhaustive_charge_fields(geometry):
    """Independent census of every ordered field in the 18-face Q=4 subset."""
    if geometry.faces != 18:
        raise ValueError('complete charge-field census is bounded to 18 faces')
    r = reference(18, 4, 'nonabelian_reflections')
    total, fields = r['canonical_connections'], 0
    covariance = np.zeros((18, 18), dtype=np.int64)
    cross = np.zeros((18, 18), dtype=np.int64)
    second = np.zeros((18, 18), dtype=np.int64)
    for row in r['population_counts']:
        n0, n1, n2 = row['populations']
        for twos in itertools.combinations(range(18), n2):
            remaining = sorted(set(range(18))-set(twos))
            for ones in itertools.combinations(remaining, n1):
                q = np.zeros(18, dtype=np.int64)
                q[list(twos)] = 2; q[list(ones)] = 1
                # Independent full-root-support drift, not the compiled terms.
                d = np.zeros(18, dtype=np.int64)
                for f, g in geometry.e.factor.pairs:
                    if (q[f] == 0) != (q[g] == 0):
                        change = q[g]-q[f]
                        d[f] += change; d[g] -= change
                for fs in geometry.e.fan_faces:
                    if sorted(q[f] for f in fs) == [0, 1, 2]:
                        for f in fs:
                            d[f] += 4*(1-q[f])
                x = 18*q-4
                weight = row['per_field_canonical_connections']
                covariance += weight*np.outer(x, x)
                cross -= weight*np.outer(d, x)
                second += weight*np.outer(d, d)
                fields += 1
    a, b, chi, k1, k2 = coefficients(ChargeMarginal(18, 4, 'nonabelian_reflections'))
    h = geometry.covariance(ChargeMarginal(18, 4, 'nonabelian_reflections'))
    for i in range(18):
        for j in range(18):
            if Fraction(int(covariance[i, j]), total*18**2) != chi*(int(i == j)-Fraction(1, 18)):
                raise ValueError('full field census disagrees with susceptibility')
            if Fraction(int(cross[i, j]), total*18) != (a+4*b)*int(geometry.l1[i, j])+b*int(geometry.l2[i, j]):
                raise ValueError('full field census disagrees with projected first moment')
            if Fraction(int(second[i, j]), total) != h[i][j]:
                raise ValueError('full field census disagrees with local covariance polynomials')
    return {'ordered_charge_fields': fields, 'canonical_connection_weight': total,
            'covariance_integer_numerator': covariance.tolist(), 'covariance_denominator': total*18**2,
            'negative_drift_charge_cross_integer_numerator': cross.tolist(), 'cross_denominator': total*18,
            'drift_second_integer_numerator': second.tolist(), 'drift_second_denominator': total}


def audit():
    bank = json.loads(Path('data/triangle-feedback.json').read_text())['bank']
    records, exhaustive = [], None
    for side in (3, 6, 12):
        geometry = ResponseGeometry(side, bank)
        for q in (2, 4, geometry.faces):
            record = result(geometry, q)
            records.append(record)
            print('response', side, q, 'closure defect', float(Fraction(*record['closure_residual_squared_norm_per_face'])), flush=True)
        if side == 3:
            exhaustive = exhaustive_charge_fields(geometry)
    return {'records': records, 'exhaustive_18_face_charge4_census': exhaustive,
            'sparse_charge4_polynomial_certificate': sparse_identity(),
            'scope': 'Exact stationary linear charge projection and two-step nonclosure on supplied geometry. The Q=2 reflection subset has an exact closed mean heat equation; reactive sectors have nonzero projection residual. No continuum diffusion coefficient, quantum amplitude, or particle binding is asserted.'}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, default=Path('out/triangle-charge-response.json'))
    args = parser.parse_args()
    data = audit()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(data, separators=(',', ':'))+'\n')


if __name__ == '__main__':
    main()
