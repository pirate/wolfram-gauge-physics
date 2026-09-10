#!/usr/bin/env python3
"""Exact disk Wilson characters in the existing finite-group torus reference.

The formal variable z counts the already conserved charge. It is not a
Boltzmann parameter supplied to the update rules. All character and fusion
data are extracted from the actual derived triangle automorphisms.
"""
import itertools
import json
import math
from collections import Counter, defaultdict
from fractions import Fraction
from functools import lru_cache
from pathlib import Path

from triangle_charge_current import CurrentProbe
from triangle_reference import polynomial_power, reference, representation_audit, subgroup


def add(*terms):
    out = [Fraction(0)]*max(len(p) for _, p in terms)
    for scale, p in terms:
        for k, a in enumerate(p):
            out[k] += scale*a
    return tuple(out)


def multiply(a, b):
    out = [0]*(len(a)+len(b)-1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i+j] += x*y
    return tuple(out)


@lru_cache(None)
def power(p, n):
    return tuple(polynomial_power(p, n))


class WilsonReference:
    def __init__(self, group, charges):
        self.g, self.q = group, charges
        audit = representation_audit(group, charges)
        self.characters = audit['characters_from_permutation_action']
        self.dimensions = audit['dimensions']
        self.p = tuple(tuple(p) for p in audit['charge_polynomials'])
        self.no_reflection_p = tuple(tuple(sum(Fraction(c[h], d) for h in range(group.n)
                                                    if charges[h] == q and q != 1)
                                              for q in range(3))
                                     for c, d in zip(self.characters, self.dimensions))
        self.fusion = [[[sum(self.characters[r][g]*self.characters[s][g]*self.characters[t][g]
                            for g in range(group.n))//group.n for t in range(3)]
                        for s in range(3)] for r in range(3)]
        assert all(sum(self.fusion[r][s][t]*self.characters[t][g] for t in range(3))
                   == self.characters[r][g]*self.characters[s][g]
                   for r, s, g in itertools.product(range(3), range(3), range(group.n)))
        reflection = next(g for g in range(group.n) if charges[g] == 1)
        self.reflection_characters = [Fraction(c[reflection], d) for c, d in zip(self.characters, self.dimensions)]

    @lru_cache(None)
    def unfiltered_numerator(self, faces, area, representation, exclude_reflection_faces=False):
        p = self.no_reflection_p if exclude_reflection_faces else self.p
        terms = []
        for r, s in itertools.product(range(3), repeat=2):
            coefficient = Fraction(self.g.n*self.dimensions[r]*self.fusion[r][s][representation],
                                   self.dimensions[s]*self.dimensions[representation])
            if coefficient:
                terms.append((coefficient, multiply(power(p[r], area), power(p[s], faces-area))))
        return add(*terms)

    @lru_cache(None)
    def numerator(self, faces, area, representation, sector='nonabelian_reflections'):
        assert 0 <= area <= faces and representation in range(3)
        full = self.unfiltered_numerator(faces, area, representation)
        if sector == 'all':
            return full
        assert sector == 'nonabelian_reflections'
        no_reflections = self.unfiltered_numerator(faces, area, representation, True)
        ua, va = power((1, 1), area), power((1, -1), area)
        ub, vb = power((1, 1), faces-area), power((1, -1), faces-area)
        # Four handle pairs in each C2. Remove its reflection-free constant
        # term before subtracting the three C2 subgroups; these disjoint
        # reflected sectors then have no inclusion-exclusion overlap.
        c2 = add((1, multiply(add((1, ua), (1, va)), add((1, ub), (1, vb)))),
                 (self.reflection_characters[representation],
                  multiply(add((1, ua), (-1, va)), add((1, ub), (-1, vb)))))
        return add((1, full), (-1, no_reflections), (-3, c2), (12, (1,)))

    def mean(self, faces, charge, area, representation, sector='nonabelian_reflections'):
        total = reference(faces, charge, sector)['canonical_connections']
        assert self.numerator(faces, area, 0, sector)[charge] == total
        return Fraction(self.numerator(faces, area, representation, sector)[charge], total)

    def bulk(self, density):
        z = (math.sqrt((3-3*density)**2+4*(4-2*density)*density)-(3-3*density))/(2*(4-2*density))
        values = [sum(float(c)*z**i for i, c in enumerate(p)) for p in self.p]
        ratios = [p/values[0] for p in values]
        return {'density': density, 'counting_fugacity': z, 'character_area_factors': ratios,
                'absolute_area_decay_rates': [-math.log(abs(p)) if p else None for p in ratios]}

    @staticmethod
    @lru_cache(None)
    def conditional_standard(inside, outside):
        """Exact E[normalized standard disk character | entire charge field].

        In the full-S3, reflected sector, the answer depends only on the
        inside/outside populations, but conditioning retains the whole field.
        The genus belongs to the exterior; swapping inside and outside is
        therefore not a symmetry of this formula.
        """
        n1, n2 = inside[1]+outside[1], inside[2]+outside[2]
        assert n1 > 0 and n1 % 2 == 0 and min(*inside, *outside) >= 0
        def eigenvalues(n):
            a = 3**n[1]*2**n[2]
            return (a, (-1)**n[1]*a, 0 if n[1] else (-1)**n[2])
        a0, a1, a2 = eigenvalues(inside)
        b0, b1, b2 = eigenvalues(outside)
        numerator = (6*a2*(b0+b1)+Fraction(3, 2)*(a0+a1)*b2+3*a2*b2
                     -12*int(n2 == 0 and inside[1] % 2 == 0))
        denominator = 12*(3**n1*2**n2-int(n2 == 0))
        return Fraction(numerator, denominator)

    def direct_small_torus_count(self, faces=4):
        """Independent sum over actual group-valued faces and handle pairs."""
        g, chars = self.g, self.characters
        handles = defaultdict(list)
        for a, b in itertools.product(range(g.n), repeat=2):
            comm = g.mul[g.inv[b]][g.mul[g.inv[a]][g.mul[b][a]]]
            handles[comm].append((a, b))
        closure_size = {mask: len(subgroup(g, [h for h in range(g.n) if mask >> h & 1]))
                        for mask in range(1 << g.n)}
        counts, numerators = Counter(), Counter()
        for word in itertools.product(range(g.n), repeat=faces):
            prefixes = [g.identity]
            for h in word:
                prefixes.append(g.mul[prefixes[-1]][h])
            charge = sum(self.q[h] for h in word)
            mask = sum(1 << h for h in set(word))
            for a, b in handles[prefixes[-1]]:
                sectors = ['all']
                if any(self.q[h] == 1 for h in word) and closure_size[mask | (1 << a) | (1 << b)] == g.n:
                    sectors.append('nonabelian_reflections')
                for sector in sectors:
                    counts[sector, charge] += 1
                    for area, h in enumerate(prefixes):
                        for t in range(3):
                            numerators[sector, charge, area, t] += chars[t][h]
        comparisons = 0
        for (sector, charge), count in counts.items():
            assert reference(faces, charge, sector)['canonical_connections'] == count
            for area, t in itertools.product(range(faces+1), range(3)):
                direct = Fraction(numerators[sector, charge, area, t], count*self.dimensions[t])
                assert self.mean(faces, charge, area, t, sector) == direct
                comparisons += 1
        return {'faces': faces, 'canonical_connections_summed': sum(n for (s, _), n in counts.items() if s == 'all'),
                'exact_character_comparisons': comparisons}


def calculation():
    bank = json.loads(Path('data/triangle-feedback.json').read_text())['bank']
    e = CurrentProbe(3, bank).e
    theory = WilsonReference(e.geometry.group, e.charges)
    rows = []
    for faces, charge in ((72, 18), (128, 32), (72, 36)):
        rows.append({'faces': faces, 'charge': charge, 'bulk': theory.bulk(charge/faces),
                     'loops': [{'area': area, 'exact_means': [str(theory.mean(faces, charge, area, t)) for t in (1, 2)],
                                'means': [float(theory.mean(faces, charge, area, t)) for t in (1, 2)]}
                               for area in (1, 2, 4, 6, 8)]})
    return {'charge_character_polynomials': theory.p, 'fusion': theory.fusion,
            'direct_count': theory.direct_small_torus_count(), 'predictions': rows,
            'scope': 'Uniform raw-connection stationary reference on supplied triangular geometry. Area dependence is not a temporal Wilson loop or a confinement/binding energy.'}


if __name__ == '__main__':
    print(json.dumps(calculation()), flush=True)
