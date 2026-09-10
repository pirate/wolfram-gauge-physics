#!/usr/bin/env python3
"""Exact finite-group reference measures and uniform raw-link initial sampling.

Reference sampling is NOT a replacement update rule or an ergodicity claim.
"""
import argparse
import itertools
import json
import math
import random
from collections import Counter
from fractions import Fraction
from pathlib import Path

from triangle_charge_current import CurrentProbe


def fraction(x):
    return [x.numerator, x.denominator]


def subgroup(g, generators):
    found = {g.identity, *generators}
    while True:
        expanded = found | {g.mul[a][b] for a in found for b in found}
        if expanded == found:
            return tuple(sorted(found))
        found = expanded


def polynomial_power(base, exponent):
    coefficients = [1]
    for _ in range(exponent):
        out = [0]*(len(coefficients)+len(base)-1)
        for i, a in enumerate(coefficients):
            for j, b in enumerate(base):
                out[i+j] += a*b
        coefficients = out
    return coefficients


def representation_audit(g, charges):
    characters = [[1]*g.n,
                  [(-1)**sum(p[i] > p[j] for i in range(3) for j in range(i+1, 3)) for p in g.elements],
                  [sum(i == value for i, value in enumerate(p))-1 for p in g.elements]]
    dimensions = [row[g.identity] for row in characters]
    if sum(d*d for d in dimensions) != g.n or any(
            sum(a*b for a, b in zip(x, y)) != g.n*int(i == j)
            for i, x in enumerate(characters) for j, y in enumerate(characters)):
        raise ValueError('derived permutation characters are not a complete irreducible set')
    polynomials = []
    for character, d in zip(characters, dimensions):
        p = [sum(Fraction(character[x], d) for x in range(g.n) if charges[x] == q) for q in range(3)]
        if any(c.denominator != 1 for c in p):
            raise ValueError('reference character polynomial is not integral')
        polynomials.append([int(c) for c in p])
    commutators = Counter()
    for a, b in itertools.product(range(g.n), repeat=2):
        value = g.mul[g.inv[b]][g.mul[g.inv[a]][g.mul[b][a]]]
        commutators[value] += 1
    expected = [g.n*sum(Fraction(row[x], d) for row, d in zip(characters, dimensions)) for x in range(g.n)]
    if [commutators[x] for x in range(g.n)] != expected:
        raise ValueError('direct commutator counts disagree with character formula')
    return {'characters_from_permutation_action': characters, 'dimensions': dimensions,
            'charge_polynomials': polynomials, 'direct_commutator_counts': [commutators[x] for x in range(g.n)]}


def field_weight(n0, n1, n2, sector='all'):
    """Canonical tree-gauge connections for ONE ordered face-charge field."""
    if min(n0, n1, n2) < 0:
        return 0
    if sector not in ('all', 'nonabelian_reflections'):
        raise ValueError('unknown reference sector')
    if n1 % 2:
        return 0
    # On the torus, summing the two handle holonomies gives |G| times the
    # product of central class-sum eigenvalues, summed over irreps. For S3 the
    # (E,R,Z) eigenvalues are (1,3,2), (1,-3,2), (1,0,-1), giving this formula.
    # Only populations enter: this reference is spatially exchangeable, not
    # evidence that a single reachable component equilibrates or binds defects.
    total = 6*(2*3**n1*2**n2 + ((-1)**n2 if n1 == 0 else 0))
    if sector == 'all':
        return total
    if n1 == 0:
        return 0
    # Remove the three single-reflection subgroups. Their intersections have
    # no reflections and were already excluded above.
    return total-(12 if n2 == 0 else 0)


def reference(faces, charge, sector='all'):
    if type(faces) is not int or not 1 <= faces <= 1152 or type(charge) is not int or not 0 <= charge <= 2*faces:
        raise ValueError('invalid bounded reference dimensions')
    rows = []
    for n2 in range(faces+1):
        n1, n0 = charge-2*n2, faces-charge+n2
        weight = field_weight(n0, n1, n2, sector)
        if weight:
            arrangements = math.comb(faces, n2)*math.comb(faces-n2, n1)
            rows.append({'populations': [n0, n1, n2], 'per_field_canonical_connections': weight,
                         'ordered_fields': arrangements, 'canonical_connections': arrangements*weight})
    total = sum(r['canonical_connections'] for r in rows)
    if not total:
        raise ValueError('reference sector is empty')
    mean = Fraction(sum(r['populations'][2]*r['canonical_connections'] for r in rows), total)
    second = Fraction(sum(r['populations'][2]**2*r['canonical_connections'] for r in rows), total)
    q_second = Fraction(sum((r['populations'][1]+4*r['populations'][2])*r['canonical_connections'] for r in rows), faces*total)
    variance = q_second-Fraction(charge, faces)**2
    return {'faces': faces, 'charge': charge, 'sector': sector, 'population_counts': rows,
            'canonical_connections': total, 'mean_rotation_count': fraction(mean),
            'rotation_count_variance': fraction(second-mean**2),
            'single_face_charge_variance': fraction(variance),
            'distinct_face_charge_covariance': fraction(-variance/(faces-1)) if faces > 1 else None,
            'scope': 'Uniform raw connections in the stated invariant sector. Defined by raw multiplicities, not an assumed gauge-orbit weighting. Not a proven reachable component or a physical temperature ensemble.'}


class AxialReference:
    def __init__(self, side, bank, charge):
        if type(side) is not int or not 3 <= side <= 12:
            raise ValueError('exact initial sampler supports sides 3..12')
        self.probe = CurrentProbe(side, bank)
        self.e, self.g = self.probe.e, self.probe.e.geometry.group
        self.side, self.faces, self.charge = side, 2*side*side, charge
        if type(charge) is not int or not 0 <= charge <= 2*self.faces:
            raise ValueError('invalid sample charge')
        # Product order from the explicit axial reconstruction, not sorted face ids.
        self.face_order = [2*(y*side+x)+i for y in reversed(range(side)) for x in range(side) for i in (1, 0)]
        self.dp = [[[int(p == self.g.identity) for p in range(self.g.n)]]]
        for length in range(1, self.faces+1):
            table = [[0]*self.g.n for _ in range(min(charge, 2*length)+1)]
            for q, values in enumerate(self.dp[-1]):
                for p, count in enumerate(values):
                    if not count:
                        continue
                    for h, cost in enumerate(self.e.charges):
                        if q+cost < len(table):
                            table[q+cost][self.g.mul[p][h]] += count
            self.dp.append(table)
        self.handles = [(a, b) for a, b in itertools.product(range(self.g.n), repeat=2)]
        self.handle_weights = [self.dp[-1][charge][self.commutator(a, b)] for a, b in self.handles]
        if not sum(self.handle_weights):
            raise ValueError('total charge has no torus connections')
        if sum(self.handle_weights) != reference(self.faces, charge)['canonical_connections']:
            raise ValueError('group convolution and character population counts disagree')

    def commutator(self, a, b):
        g = self.g
        return g.mul[g.inv[b]][g.mul[g.inv[a]][g.mul[b][a]]]

    @staticmethod
    def choose(rng, weights):
        total = sum(weights)
        if not total or any(w < 0 for w in weights):
            raise ValueError('invalid integer conditional sampling weights')
        index = rng.randrange(total)
        for i, weight in enumerate(weights):
            if index < weight:
                return i
            index -= weight
        raise AssertionError('integer weighted draw missed its interval')

    def directed(self, links, x, y, dx, dy):
        s = self.side
        u, v = (y % s)*s+x % s, ((y+dy) % s)*s+(x+dx) % s
        value = links[self.e.geometry.ids[tuple(sorted((u, v)))]]
        return value if u < v else self.g.inv[value]

    def put(self, links, x, y, dx, dy, value):
        s = self.side
        u, v = (y % s)*s+x % s, ((y+dy) % s)*s+(x+dx) % s
        links[self.e.geometry.ids[tuple(sorted((u, v)))]] = value if u < v else self.g.inv[value]

    def cell_holonomies(self, links):
        # The engine cyclically rotates face loops to their smallest vertex id.
        # Axial formulas instead base both triangles at the cell's lower-left
        # corner; periodic wrap faces therefore need explicit based products.
        result, g = [], self.g
        for y in range(self.side):
            for x in range(self.side):
                h = self.directed(links, x, y, 1, 0)
                v = self.directed(links, x, y, 0, 1)
                d = self.directed(links, x, y, 1, 1)
                right = self.directed(links, x+1, y, 0, 1)
                top = self.directed(links, x, y+1, 1, 0)
                result.extend((g.mul[g.inv[d]][g.mul[right][h]],
                               g.mul[g.inv[v]][g.mul[g.inv[top]][d]]))
        if [self.e.charges[x] for x in result] != self.probe.charge(links):
            raise ValueError('cell-based and native face charges differ')
        return result

    def reconstruct(self, a, b, holonomies):
        s, g = self.side, self.g
        if any(type(x) is not int or not 0 <= x < g.n for x in (a, b)):
            raise ValueError('invalid handle elements')
        if len(holonomies) != self.faces or any(type(x) is not int or not 0 <= x < g.n for x in holonomies):
            raise ValueError('invalid prescribed face holonomies')
        product = g.identity
        for f in self.face_order:
            product = g.mul[product][holonomies[f]]
        if product != self.commutator(a, b):
            raise ValueError('prescribed face product violates the handle commutator')
        squares = [[g.mul[holonomies[2*(y*s+x)+1]][holonomies[2*(y*s+x)]] for x in range(s)] for y in range(s)]
        vertical, wraps = [], [a]
        for y in range(s):
            row = [b if y == s-1 else g.identity]
            for x in range(s-1):
                row.append(g.mul[row[-1]][squares[y][x]])
            vertical.append(row)
            product = g.identity
            for value in squares[y]:
                product = g.mul[product][value]
            if y < s-1:
                wraps.append(g.mul[wraps[-1]][g.inv[product]])
            elif g.mul[b][g.mul[wraps[-1]][g.mul[g.inv[product]][g.inv[b]]]] != a:
                raise ValueError('axial wrap recurrence does not close')
        links = [g.identity]*len(self.e.geometry.edges)
        for y in range(s):
            for x in range(s):
                h = wraps[y] if x == s-1 else g.identity
                self.put(links, x, y, 1, 0, h)
                self.put(links, x, y, 0, 1, vertical[y][x])
                diagonal = g.mul[vertical[y][(x+1) % s]][g.mul[h][g.inv[holonomies[2*(y*s+x)]]]]
                self.put(links, x, y, 1, 1, diagonal)
        if self.cell_holonomies(links) != list(holonomies):
            raise ValueError('reconstructed actual links do not realize every prescribed holonomy')
        return links

    def transform(self, links, frames):
        g = self.g
        return [g.mul[frames[v]][g.mul[x][g.inv[frames[u]]]] for (u, v), x in zip(self.e.geometry.edges, links)]

    def canonicalize(self, links):
        s, g = self.side, self.g
        paths = [g.identity]*(s*s)
        for y in range(s):
            if y:
                paths[y*s] = g.mul[self.directed(links, 0, y-1, 0, 1)][paths[(y-1)*s]]
            for x in range(1, s):
                paths[y*s+x] = g.mul[self.directed(links, x-1, y, 1, 0)][paths[y*s+x-1]]
        canonical = self.transform(links, [g.inv[x] for x in paths])
        a, b = self.directed(canonical, s-1, 0, 1, 0), self.directed(canonical, 0, s-1, 0, 1)
        holonomies = self.cell_holonomies(canonical)
        if self.reconstruct(a, b, holonomies) != canonical:
            raise ValueError('axial extraction/reconstruction is not a raw-connection bijection')
        return canonical, a, b, holonomies, paths

    def sample(self, rng, sector='all'):
        # Rejection is only conditioning an exact reference draw, never a physical update.
        reference(self.faces, self.charge, sector)
        for attempt in range(1, 10001):
            a, b = self.handles[self.choose(rng, self.handle_weights)]
            q, product, values = self.charge, self.commutator(a, b), []
            for remaining in reversed(range(self.faces)):
                weights = []
                for h, cost in enumerate(self.e.charges):
                    residual = self.g.mul[self.g.inv[h]][product]
                    weights.append(self.dp[remaining][q-cost][residual] if 0 <= q-cost < len(self.dp[remaining]) else 0)
                if sum(weights) != self.dp[remaining+1][q][product]:
                    raise ValueError('conditional draw weights do not exhaust the completion count')
                h = self.choose(rng, weights)
                values.append(h); q -= self.e.charges[h]
                product = self.g.mul[self.g.inv[h]][product]
            if q or product != self.g.identity:
                raise ValueError('conditional draw missed the exact product/charge constraint')
            holonomies = [0]*self.faces
            for f, value in zip(self.face_order, values):
                holonomies[f] = value
            canonical = self.reconstruct(a, b, holonomies)
            image = subgroup(self.g, canonical)
            if sector == 'nonabelian_reflections' and (len(image) != 6 or not any(self.e.charges[h] == 1 for h in values)):
                continue
            frames = [self.g.identity]+[rng.randrange(self.g.n) for _ in range(self.side*self.side-1)]
            links = self.transform(canonical, frames)
            extracted, aa, bb, hh, recovered = self.canonicalize(links)
            if (extracted, aa, bb, hh, recovered) != (canonical, a, b, holonomies, frames):
                raise ValueError('reference sample does not round-trip its rooted gauge frames')
            observed_group = subgroup(self.g, self.e.forest.based_loops(links)[0])
            if len(observed_group) != len(image) or sum(self.probe.charge(links)) != self.charge:
                raise ValueError('reference draw has the wrong global charge or holonomy image')
            return {'handles': [a, b], 'cell_based_holonomies': holonomies,
                    'root_fixed_frames': frames, 'links': links, 'holonomy_group_order': len(image),
                    'reference_draws_until_acceptance': attempt}
        raise RuntimeError('bounded reference rejection exhausted; no sample returned')


def sector_audit(bank):
    probe = CurrentProbe(3, bank)
    e, g = probe.e, probe.e.geometry.group
    subgroups = sorted({subgroup(g, values) for count in range(4) for values in itertools.combinations(range(g.n), count)})
    checked = 0
    for rule, table in enumerate(e.tables):
        arity = 3 if rule else 2
        for code, target in enumerate(table):
            before = [(code//g.n**i) % g.n for i in reversed(range(arity))]
            after = [(target//g.n**i) % g.n for i in reversed(range(arity))]
            if subgroup(g, before) != subgroup(g, after):
                raise ValueError('primitive changes its locally generated subgroup')
            if any(e.charges[x] == 1 for x in before) != any(e.charges[x] == 1 for x in after):
                raise ValueError('primitive creates or removes the last reflection')
            checked += 1
    # Independent direct enumeration of a small punctured-torus presentation.
    handles = [(a, b, g.mul[g.inv[b]][g.mul[g.inv[a]][g.mul[b][a]]]) for a, b in itertools.product(range(g.n), repeat=2)]
    counts = [Counter(), Counter()]
    for values in itertools.product(range(g.n), repeat=4):
        product = g.identity
        for value in values:
            product = g.mul[product][value]
        field = tuple(e.charges[x] for x in values)
        for a, b, commutator in handles:
            if product != commutator:
                continue
            counts[0][field] += 1
            if 1 in field and len(subgroup(g, [a, b, *values])) == 6:
                counts[1][field] += 1
    for field in itertools.product(range(3), repeat=4):
        for i, sector in enumerate(('all', 'nonabelian_reflections')):
            if counts[i][field] != field_weight(*(field.count(q) for q in range(3)), sector):
                raise ValueError('subgroup-filtered character weight fails direct tuple enumeration')
    return {'all_subgroups': subgroups, 'local_rule_inputs_checked': checked,
            'four_face_presentation_counts': [sorted(c.items()) for c in counts],
            'invariants': 'The global based-loop subgroup up to conjugation and presence of any reflection face are preserved. Boundary-fixed word reconstruction preserves subgroup membership; inverse updates give equality.'}


def constant_connection(probe, directions):
    side, g = probe.e.geometry.side, probe.e.geometry.group
    links = [g.identity]*len(probe.e.geometry.edges)
    for y in range(side):
        for x in range(side):
            u = side*y+x
            for (dx, dy), value in zip(((1, 0), (0, 1), (1, 1)), directions):
                v = side*((y+dy) % side)+(x+dx) % side
                links[probe.e.geometry.ids[tuple(sorted((u, v)))]] = value if u < v else g.inv[value]
    return links


class ActivityProbe(CurrentProbe):
    """Exact O(E) activity count from based star loops, not a trajectory timeout."""
    def __init__(self, side, bank):
        super().__init__(side, bank)
        fan = self.e.factor.oracle.fan
        rooted, actual_fans = {}, {}
        for spec, paths in zip(fan.specs, fan.patches):
            u, a, b, _ = spec[0]
            rooted.setdefault(u, {})[a] = (b, paths[0])
            actual_fans.setdefault(u, set()).add(tuple(tuple(p) for p in paths))
        self.stars = []
        for u, successors in sorted(rooted.items()):
            first = current = min(successors)
            paths, vertices = [], []
            for _ in range(6):
                vertices.append(current)
                current, path = successors[current]
                paths.append(path)
            if len(set(vertices)) != 6 or current != first or len(successors) != 6:
                raise ValueError('activity criterion requires a six-face cyclic star')
            # Verify that these are exactly the actual rooted fan supports.
            expected = {tuple(tuple(paths[(i+j) % 6]) for j in range(3)) for i in range(6)}
            if actual_fans[u] != expected:
                raise ValueError('cyclic star does not reproduce every primitive fan')
            self.stars.append(paths)

    @staticmethod
    def star_activity(values, charges):
        equality = [values[i] == values[(i+1) % 6] for i in range(6)]
        q = [charges[h] for h in values]
        active = sum(all(q[(i+j) % 6] == 1 for j in range(3)) and equality[i] != equality[(i+1) % 6]
                     for i in range(6))
        reverse = sum(sorted(q[(i+j) % 6] for j in range(3)) == [0, 1, 2] for i in range(6))
        return active, reverse

    def measure(self, links):
        q, e = self.charge(links), self.e
        vacancy = sum((q[a] == 0) != (q[b] == 0) for a, b in e.factor.pairs)
        forward = reverse = 0
        star_counts = Counter()
        for paths in self.stars:
            values = [e.factor.oracle.fan.transport(links, path) for path in paths]
            a, b = self.star_activity(values, e.charges)
            forward += a; reverse += b
            star_counts[a, b] += 1
        changing = vacancy+12*forward+4*reverse
        frozen_criterion = not any(q) or (all(q) and forward == 0)
        if frozen_criterion != (changing == 0):
            raise ValueError('connected-mesh frozen criterion disagrees with exact activity')
        if 0 < sum(q) < self.faces and changing == 0:
            raise ValueError('positive subunit charge density cannot be frozen')
        return {'charge': sum(q), 'populations': [q.count(i) for i in range(3)],
                'vacancy_operators': vacancy, 'active_reflection_fans': forward,
                'reverse_conversion_fans': reverse, 'changing_operators': changing,
                'attempted_operators': len(e.tables)*len(e.factor.pairs),
                'star_activity_histogram': [[a, b, n] for (a, b), n in sorted(star_counts.items())],
                'frozen': frozen_criterion}


def activity_audit(bank):
    p = ActivityProbe(3, bank)
    e, g = p.e, p.e.geometry.group
    # Raw-code activity, not just charge differences, rules out hidden link-only moves.
    for code in range(g.n**2):
        a, b = divmod(code, g.n)
        if (e.tables[0][code] != code) != ((e.charges[a] == 0) != (e.charges[b] == 0)):
            raise ValueError('vacancy activity formula fails raw table')
    for values in itertools.product(range(g.n), repeat=3):
        code = (values[0]*g.n+values[1])*g.n+values[2]
        q = [e.charges[h] for h in values]
        expected = 4 if sorted(q) == [0, 1, 2] else 12*int(q == [1]*3 and ((values[0] == values[1]) != (values[1] == values[2])))
        if sum(table[code] != code for table in e.tables[1:]) != expected:
            raise ValueError('reaction activity formula fails raw tables')
    reflections = [h for h in range(g.n) if e.charges[h] == 1]
    star_counts = Counter()
    for values in itertools.product(reflections, repeat=6):
        active, reverse = p.star_activity(values, e.charges)
        equal = [values[i] == values[(i+1) % 6] for i in range(6)]
        if reverse or active % 2 or (active == 0) != (all(equal) or not any(equal)):
            raise ValueError('homogeneous star wall criterion fails')
        # Independent lookup of all six actual reversed fan triples.
        actual = 0
        for i in range(6):
            a, b, c = [values[(i+j) % 6] for j in (2, 1, 0)]
            code = (a*g.n+b)*g.n+c
            actual += sum(table[code] != code for table in e.tables[1:])
        if actual != 12*active:
            raise ValueError('star equality walls fail the complete bank')
        star_counts[active] += 1
    # Every operator is followed by its inverse (the same involution), so the
    # compiled event count measures degree at the initial state, not along a walk.
    rng = random.Random(501001)
    inputs = [[0]*len(e.geometry.edges), constant_connection(p, (1, 1, 1)),
              constant_connection(p, (1, 2, 5)), constant_connection(p, (3, 3, 3))]
    inputs += [[rng.randrange(g.n) for _ in e.geometry.edges] for _ in range(4)]
    inputs += [[rng.choice(reflections) for _ in e.geometry.edges] for _ in range(4)]
    inputs += [AxialReference(3, bank, q).sample(rng)['links'] for q in (2, 4, 8, 16)]
    schedule = [op for op in range(len(e.tables)*len(e.factor.pairs)) for _ in range(2)]
    checked = e.compiled_bank(inputs, schedule, len(schedule), capture_links=True)
    measures = [p.measure(links) for links in inputs]
    for links, m in zip(inputs[8:12], measures[8:12]):
        if m['populations'] != [0, p.faces, 0] or m['frozen'] or len(subgroup(g, e.forest.based_loops(links)[0])) != g.n:
            raise ValueError('homogeneous reflection controls must be active with full nonabelian holonomy')
    for run in checked['runs']:
        expected = measures[run['condition']]
        degree = expected['changing_operators'] if run['mode'] == 'combined' else expected['vacancy_operators']
        if run['final_links'] != inputs[run['condition']] or len(run['events']) != 2*degree:
            raise ValueError('compiled operator/inverse census disagrees with star activity')
    return {'local_pair_inputs_checked': g.n**2, 'local_triple_inputs_checked': g.n**3,
            'reflection_star_inputs_checked': len(reflections)**6,
            'reflection_star_active_fan_histogram': sorted(star_counts.items()),
            'compiled_initial_states': [{'links': links, 'activity': m} for links, m in zip(inputs, measures)],
            'scope': 'Exact frozen-state criterion on the supplied connected triangular torus. No irreducibility, mixing-time, phase-transition, or emergent-geometry claim.'}


def frozen_audit(bank):
    probe = CurrentProbe(3, bank)
    e, g = probe.e, probe.e.geometry.group
    rows, frozen = [], []
    for directions in itertools.product(range(g.n), repeat=3):
        links = constant_connection(probe, directions)
        if probe.charge(links) != [1]*probe.faces:
            continue
        changing = sum(table[e.factor.oracle.code(links, 1, patch)] != e.factor.oracle.code(links, 1, patch)
                       for patch in range(len(e.factor.pairs)) for table in e.tables[1:])
        image = subgroup(g, e.forest.based_loops(links)[0])
        rows.append({'direction_links': directions, 'changing_operators': changing, 'holonomy_group_order': len(image)})
        if not changing:
            frozen.append(links)
    schedule = list(range(len(e.tables)*len(e.factor.pairs)))
    for start in range(0, len(frozen), 16):
        batch = frozen[start:start+16]
        checked = e.compiled_bank(batch, schedule, len(schedule))
        if any(run['events'] or run['final_links'] != batch[run['condition']] for run in checked['runs']):
            raise ValueError('claimed frozen connection moves under a compiled operator')
    witnesses = []
    for side in (3, 6, 12):
        p = CurrentProbe(side, bank)
        ee = p.e
        links = constant_connection(p, (1, 2, 5))
        full_scan = list(range(len(ee.tables)*len(ee.factor.pairs)))
        checked = ee.compiled_bank([links], full_scan, len(full_scan))
        if any(run['events'] or run['final_links'] != links for run in checked['runs']):
            raise ValueError('nonabelian frozen witness moves')
        loops = ee.forest.based_loops(links)[0]
        if len(subgroup(ee.geometry.group, loops)) != 6 or p.charge(links) != [1]*p.faces:
            raise ValueError('frozen witness does not have full nonabelian reflection holonomy')
        root_patches = [i for i, spec in enumerate(ee.factor.oracle.fan.specs) if spec[0][0] == 0]
        triples = []
        for patch in root_patches:
            code = ee.factor.oracle.code(links, 1, patch)
            triples.append([code//36, code//6 % 6, code % 6])
        witnesses.append({'side': side, 'links': links, 'global_based_loops': loops,
                          'root_star_fan_triples': triples, 'all_compiled_operators_fixed': True})
    return {'constant_direction_homogeneous_charge_census': rows, 'frozen_count': len(frozen),
            'nonabelian_frozen_count': sum(r['changing_operators'] == 0 and r['holonomy_group_order'] == 6 for r in rows),
            'nonabelian_frozen_witnesses': witnesses,
            'scope': 'Complete census only of the 216 constant-direction connections, not all frozen states. Nonabelian reflection sector is not a single reachable component.'}


def audit():
    source = json.loads(Path('data/triangle-feedback.json').read_text())
    probe = CurrentProbe(3, source['bank'])
    rep = representation_audit(probe.e.geometry.group, probe.e.charges)
    references, samples = [], []
    for side, charge in ((3, 4), (3, 18), (6, 4), (6, 72), (12, 288)):
        sampler = AxialReference(side, source['bank'], charge)
        for sector in ('all', 'nonabelian_reflections'):
            references.append(reference(sampler.faces, charge, sector))
            rng = random.Random(180100+side*100+charge)
            draws = [sampler.sample(rng, sector) for _ in range(8)]
            samples.append({'side': side, 'charge': charge, 'sector': sector, 'draws': draws})
        print('reference', side, charge, 'exact sampler round-trips passed', flush=True)
    stationary_checks = []
    for sample in samples:
        if sample['sector'] != 'nonabelian_reflections' or (sample['side'], sample['charge']) not in ((6, 72), (12, 288)):
            continue
        p = CurrentProbe(sample['side'], source['bank'])
        ee = p.e
        rng = random.Random(92100+sample['side'])
        schedule = [rng.randrange(len(ee.tables)*len(ee.factor.pairs)) for _ in range(2000)]
        checked = ee.compiled_bank([draw['links'] for draw in sample['draws']], schedule, 200)
        for run in checked['runs']:
            if len(subgroup(ee.geometry.group, ee.forest.based_loops(run['final_links'])[0])) != 6 or 1 not in p.charge(run['final_links']):
                raise ValueError('compiled evolution leaves the conditioned reference sector')
        stationary_checks.append({'side': sample['side'], 'schedule': schedule, 'checked': checked,
                                  'scope': 'A fixed composition of the original involutions preserves the exact reference law. These few realized samples do not establish mixing.'})
    return {'representations': rep, 'references': references, 'samples': samples,
            'sector_checks': sector_audit(source['bank']), 'frozen': frozen_audit(source['bank']),
            'activity': activity_audit(source['bank']),
            'stationary_evolution_checks': stationary_checks,
            'scope': 'Exact uniform raw-connection reference measures and independently drawn initial states. Initial sampling is not physical evolution, a mixing proof, or emergent thermodynamics.'}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, default=Path('out/triangle-reference.json'))
    args = parser.parse_args()
    data = audit()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(data, separators=(',', ':'))+'\n')


if __name__ == '__main__':
    main()
