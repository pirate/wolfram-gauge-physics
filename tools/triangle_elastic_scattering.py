#!/usr/bin/env python3
"""Exhaustive charge-preserving pair laws and a separate inverse-paired bank.

The fiber group and charge come from the existing derivation. Selection below
is algebraic, not based on whether a trajectory resembles a desired particle.
"""
import argparse
import hashlib
import itertools
import json
import math
import random
import subprocess
from collections import Counter
from pathlib import Path

from triangle_charge_response import ResponseGeometry
from triangle_patch_observer import PatchUnion, WordObserver
from triangle_reference import AxialReference, constant_connection, subgroup
from triangle_response_memory import activity_mask, raw_targets, second_drift


def inverse(table):
    if sorted(table) != list(range(len(table))):
        raise ValueError('expected a permutation table')
    result = [0]*len(table)
    for x, y in enumerate(table):
        result[y] = x
    return result


def order(table):
    inverse(table)
    seen, result = set(), 1
    for x in range(len(table)):
        length = 0
        while x not in seen:
            seen.add(x); length += 1; x = table[x]
        if length:
            result = math.lcm(result, length)
    return result


def pair_census(g, charges):
    """All equivariant product/charge-preserving bijections, not just involutions.

    Fix ab=p. This fiber is parameterized by a, with b=a^-1 p. Its permitted
    permutations must commute with C(p); transporting one such permutation to
    every conjugate p gives a well-defined full pair map. Conversely every
    equivariant full map restricts to these centralizer-equivariant choices.
    """
    n = g.n
    reps = sorted(set(g.sectors))
    choices, fibers = {}, []
    for p in reps:
        centralizer = [h for h in range(n) if g.conj[h][p] == p]
        cost = [charges[a]+charges[g.mul[g.inv[a]][p]] for a in range(n)]
        choices[p] = [list(f) for f in itertools.permutations(range(n))
                      if all(cost[a] == cost[f[a]] for a in range(n))
                      and all(f[g.conj[h][a]] == g.conj[h][f[a]]
                              for h in centralizer for a in range(n))]
        fibers.append({'product': p, 'centralizer': centralizer,
                       'factorization_charges': cost, 'permutation_count': len(choices[p])})
    transports = [(g.sectors[p], next(h for h in range(n) if g.conj[h][g.sectors[p]] == p))
                  for p in range(n)]
    tables = []
    for selected in itertools.product(*(choices[p] for p in reps)):
        functions, table = dict(zip(reps, selected)), []
        for a, b in itertools.product(range(n), repeat=2):
            p = g.mul[a][b]
            representative, h = transports[p]
            x = g.conj[h][functions[representative][g.conj[g.inv[h]][a]]]
            table.append(n*x+g.mul[g.inv[x]][p])
        tables.append(table)
    tables.sort()
    if len({tuple(t) for t in tables}) != len(tables):
        raise ValueError('centralizer construction duplicated a law')
    j = [n*g.inv[b]+g.inv[a] for a, b in itertools.product(range(n), repeat=2)]
    commute, reversal, orders = Counter(), Counter(), Counter()
    for table in tables:
        g.validate(table)
        o, inv = order(table), inverse(table)
        orders[o] += 1
        if all(table[j[x]] == j[table[x]] for x in range(n*n)):
            commute[o] += 1
        if all(j[table[j[x]]] == inv[x] for x in range(n*n)):
            reversal[o] += 1
    elastic = [t for t in tables if order(t) == 3]
    # Identify the exhaustively derived candidates with the known restricted
    # Hurwitz formula only AFTER enumeration; it did not prune the search.
    hurwitz = []
    for a, b in itertools.product(range(n), repeat=2):
        x, y = ((g.mul[a][g.mul[b][g.inv[a]]], a)
                if charges[a] == charges[b] == 1 and a != b else (a, b))
        hurwitz.append(n*x+y)
    if sorted(elastic) != sorted([hurwitz, inverse(hurwitz)]):
        raise ValueError('order-three census is not the two restricted Hurwitz laws')
    minimum = min(sum(i != x for i, x in enumerate(t)) for t in tables if order(t) > 2)
    minimal = [t for t in tables if order(t) > 2 and sum(i != x for i, x in enumerate(t)) == minimum]
    if sorted(minimal) != sorted(elastic):
        raise ValueError('minimal-support non-involutions differ from elastic candidates')
    def adjacent(values, position, table):
        moved = list(values)
        moved[position:position+2] = divmod(table[n*values[position]+values[position+1]], n)
        return tuple(moved)
    braid_failures = []
    for table in (hurwitz, inverse(hurwitz)):
        count = 0
        for values in itertools.product(range(n), repeat=3):
            left = adjacent(adjacent(adjacent(values, 0, table), 1, table), 0, table)
            right = adjacent(adjacent(adjacent(values, 1, table), 0, table), 1, table)
            count += left != right
        braid_failures.append(count)
    summary = lambda c: {str(k): v for k, v in sorted(c.items())}
    return {'product_fibers': fibers, 'bijections': len(tables),
            'order_counts': summary(orders), 'commuting_orientation_order_counts': summary(commute),
            'reversing_orientation_order_counts': summary(reversal),
            'minimum_noninvolution_support': minimum,
            'elastic_braid_failures_on_all_triples': braid_failures,
            'elastic_tables': [hurwitz, inverse(hurwitz)], 'tables': tables}


class ElasticExperiment:
    """Explicit alternative bank: vacancy, H, H^-1, then original 12 reactions."""
    def __init__(self, side, bank, elastic):
        if len(elastic) != 2 or inverse(elastic[0]) != elastic[1]:
            raise ValueError('expected two inverse-paired elastic laws')
        self.geometry = ResponseGeometry(side, bank)
        self.e, self.p = self.geometry.e, self.geometry.p
        self.tables = [self.e.tables[0]]+elastic+self.e.tables[1:]
        self.inverses = [inverse(t) for t in self.tables]
        self.supports = len(self.e.factor.pairs)
        self.operator_count = len(self.tables)*self.supports
        self.inverse_slots = [self.tables.index(t) for t in self.inverses]

    def inverse_operator(self, encoded):
        if type(encoded) is not int or not 0 <= encoded < self.operator_count:
            raise ValueError('invalid elastic-bank operator')
        rule, patch = divmod(encoded, self.supports)
        return self.inverse_slots[rule]*self.supports+patch

    def apply(self, links, encoded, reverse=False):
        if type(encoded) is not int or not 0 <= encoded < self.operator_count:
            raise ValueError('invalid elastic-bank operator')
        rule, patch = divmod(encoded, self.supports)
        return self.e.factor.oracle.update(links, int(rule >= 3), patch,
                                          (self.inverses if reverse else self.tables)[rule])

    def compiled(self, inputs, schedule, stride=1):
        if not inputs or len(inputs) > 16 or not schedule or len(schedule) > 1000000:
            raise ValueError('invalid bounded elastic run dimensions')
        if type(stride) is not int or stride <= 0 or len(schedule) % stride:
            raise ValueError('invalid elastic sampling stride')
        if any(type(op) is not int or not 0 <= op < self.operator_count for op in schedule):
            raise ValueError('invalid elastic schedule')
        for links in inputs:
            self.e.geometry.holonomies(links)
        protocol = [self.e.geometry.side, len(self.tables)-3, len(inputs), len(schedule), stride, 3]
        protocol += [x for t in self.tables for x in t]
        protocol += [x for links in inputs for x in links]+schedule
        result = json.loads(subprocess.run(
            ['build/wgphysics_mixed_bank_experiments', '--cycle', '3', '--pair-bank', '--raw-events'],
            input=' '.join(map(str, protocol))+'\n', text=True, capture_output=True,
            check=True, timeout=120).stdout)
        if len(result['runs']) != 2*len(inputs):
            raise ValueError('compiled elastic run count differs')
        for row, (condition, combined) in zip(result['runs'], itertools.product(range(len(inputs)), (False, True))):
            links, events = inputs[condition][:], []
            history = [self.e.factor.oracle.histogram(links)]
            charge = sum(self.p.charge(links))
            if row['condition'] != condition or row['mode'] != ('combined' if combined else 'transport'):
                raise ValueError('compiled elastic condition differs')
            for tick, op in enumerate(schedule, 1):
                rule, patch = divmod(op, self.supports)
                if combined or rule < 3:
                    old, before = links[:], self.e.geometry.holonomies(links)
                    code, target = self.apply(links, op)
                    after = self.e.geometry.holonomies(links)
                    faces = set(self.e.fan_faces[patch] if rule >= 3 else self.e.factor.pairs[patch])
                    writes = (self.e.factor.oracle.fan.writes[patch] if rule >= 3 else
                              {self.e.factor.oracle.pairs[patch][0][0][0]})
                    if any(a != b for i, (a, b) in enumerate(zip(old, links)) if i not in writes):
                        raise ValueError('elastic primitive changed a boundary link')
                    if any(a != b for i, (a, b) in enumerate(zip(before, after)) if i not in faces):
                        raise ValueError('elastic primitive changed a spectator holonomy')
                    if sum(self.p.charge(links)) != charge:
                        raise ValueError('elastic bank violated conserved charge')
                    if rule in (1, 2) and self.p.charge(old) != self.p.charge(links):
                        raise ValueError('elastic primitive changed a face charge')
                    if target != code:
                        events.append([tick, rule, patch, code, target, links[:]])
                    elif old != links:
                        raise ValueError('identity local target changed raw links')
                if tick % stride == 0:
                    history.append(self.e.factor.oracle.histogram(links))
            if events != row['events'] or links != row['final_links'] or history != row['histograms']:
                raise ValueError('C++ elastic trajectory differs from independent raw replay')
            for op in reversed(schedule):
                if combined or op//self.supports < 3:
                    self.apply(links, op, reverse=True)
            if links != inputs[condition] or not row['exact_link_inverse'] or not row['reverse_histogram_echo']:
                raise ValueError('elastic trajectory did not reverse exactly')
        return result


def encounter_witness(bank, elastic):
    """Search exact Q=4 initializations for an observable effect of one collision.

    This is an existence witness, not an unbiased estimate of collision rates.
    Subsequent readout uses the UNCHANGED baseline bank, isolating preparation.
    """
    experiment = ElasticExperiment(3, bank, elastic)
    geometry, p, e = experiment.geometry, experiment.p, experiment.e
    reference, rng = AxialReference(3, bank, 4), random.Random(907331)
    words = WordObserver(e.geometry.group, e.charges)
    for sample in range(512):
        initial = reference.sample(rng, 'nonabelian_reflections')['links']
        q = p.charge(initial)
        active = activity_mask(geometry, initial, q)
        for patch in range(experiment.supports):
            moved = initial[:]
            code, target = experiment.apply(moved, experiment.supports+patch)
            if code == target:
                continue
            after = activity_mask(geometry, moved, q)
            if active == after:
                continue
            old_rates, new_rates = raw_targets(geometry, initial, compiled=True), raw_targets(geometry, moved, compiled=True)
            if old_rates == new_rates:
                continue
            m = 39*geometry.faces
            old_rates[tuple(q)] += m-sum(old_rates.values())
            new_rates[tuple(q)] += m-sum(new_rates.values())
            delta = [b-a for a, b in zip(second_drift(geometry, q, active), second_drift(geometry, q, after))]
            if not any(delta):
                continue
            before_loops = words.canonical(e.forest.based_loops(initial)[0])
            after_loops = words.canonical(e.forest.based_loops(moved)[0])
            if before_loops == after_loops or p.charge(moved) != q:
                raise ValueError('collision witness lost its gauge/charge distinction')
            scan = [experiment.supports+patch, 2*experiment.supports+patch]
            experiment.compiled([initial], scan)
            return {'seed': 907331, 'sample_index': sample, 'side': 3, 'total_charge': 4,
                    'initial_links': initial, 'scattered_links': moved, 'patch': patch,
                    'local_codes': [code, target], 'charge_field': q,
                    'active_reaction_patches': [list(active), list(after)],
                    'global_gauge_orbits': [list(before_loops), list(after_loops)],
                    'readout_operators': m,
                    'one_step_charge_law_l1_numerator': sum(abs(old_rates[k]-new_rates[k]) for k in old_rates.keys() | new_rates.keys()),
                    'two_step_mean_difference_numerator': delta,
                    'two_step_mean_denominator': m*m,
                    'cpp_scatter_inverse_and_readout_verified': True}
    raise ValueError('bounded search found no gauge-mediated charge-readout witness')


def automatic_response(bank, elastic, witness):
    """Exact third-step effect of adding elastic rules, with matched clock.

    Let L0 be the original unnormalized generator and E the two elastic sums.
    E q = E L0 q = 0, since both q and its baseline drift depend only on face
    charges. Thus (P1^3-P0^3)q = E L0^2 q / M^3, M=45 F. P0 pads the baseline
    bank with two identity pair slots; it is not the old 39 F attempt clock.
    """
    experiment = ElasticExperiment(witness['side'], bank, elastic)
    geometry, p, e = experiment.geometry, experiment.p, experiment.e
    initial = witness['initial_links']

    def microscopic_ld(links):
        q, result = p.charge(links), [0]*geometry.faces
        d = geometry.drift(q)
        # Every baseline operator is independently run and undone in C++.
        # Only its measured output charges enter this verification.
        for target, count in raw_targets(geometry, links, compiled=True).items():
            result = [x+count*(a-b) for x, a, b in zip(result, geometry.drift(target), d)]
        if result != second_drift(geometry, q, activity_mask(geometry, links, q)):
            raise ValueError('raw-operator second drift disagrees with observer formula')
        return result

    base = microscopic_ld(initial)
    targets, scan = Counter(), []
    for rule in (1, 2):
        for patch in range(experiment.supports):
            moved = initial[:]
            op = rule*experiment.supports+patch
            code, target = experiment.apply(moved, op)
            if p.charge(moved) != p.charge(initial):
                raise ValueError('elastic generator does not annihilate charge')
            targets[tuple(moved)] += 1
            scan += [op, experiment.inverse_operator(op)]
    experiment.compiled([initial], scan, stride=len(scan))
    delta, rows = [0]*geometry.faces, []
    for links, multiplicity in sorted(targets.items()):
        if list(links) == initial:
            continue
        ld = microscopic_ld(list(links))
        contribution = [multiplicity*(a-b) for a, b in zip(ld, base)]
        delta = [a+b for a, b in zip(delta, contribution)]
        rows.append({'links': list(links), 'operator_multiplicity': multiplicity,
                     'third_step_numerator_contribution': contribution})
    if not any(delta) or sum(delta):
        raise ValueError('automatic scattering response vanished or violated conservation')
    return {'control': 'Original 13 rules plus two identity pair slots; same 15-slot uniform attempted clock.',
            'attempted_operators': experiment.operator_count,
            'mean_difference_zero_through_attempt': 2,
            'third_step_mean_difference_numerator': delta,
            'third_step_mean_denominator': experiment.operator_count**3,
            'elastic_target_contributions': rows,
            'cpp_all_scatter_and_baseline_readout_scans_verified': True}


def patch_activation(bank, elastic):
    """Complete three-face read-graph kernel, not an autonomous mesh quotient."""
    experiment = ElasticExperiment(3, bank, elastic)
    union = PatchUnion(experiment.e, (0,))
    old_channels, baseline, _ = union.kernel()
    channels = [(patch, rule+2 if rule else 0) for patch, rule in old_channels]
    channels += [(patch, rule) for patch, old_rule in old_channels if not old_rule for rule in (1, 2)]
    tables, scan = [], []
    for patch, rule in channels:
        table = []
        op = rule*experiment.supports+patch
        inverse_op = experiment.inverse_operator(op)
        scan += [op, inverse_op]
        for rep in union.reps:
            raw = union.lift(rep)
            experiment.apply(raw, op)
            table.append(union.observe(raw))
        inverse(table)
        tables.append(table)
    if tables[:len(baseline)] != baseline:
        raise ValueError('elastic kernel changed the baseline channels')
    if any(inverse(t) not in tables for t in tables):
        raise ValueError('patch bank is not inverse-closed')
    for start in range(0, len(union.reps), 16):
        reps = union.reps[start:start+16]
        runs = experiment.compiled([union.lift(rep) for rep in reps], scan, stride=len(scan))['runs'][1::2]
        for index, run in enumerate(runs, start):
            events = {row[0]: row[-1] for row in run['events']}
            for k, table in enumerate(tables):
                actual = union.observe(events[2*k+1]) if 2*k+1 in events else index
                if actual != table[index]:
                    raise ValueError('C++ patch kernel differs from gauge quotient')

    def components(rows):
        seen, result = set(), []
        for start in range(len(union.reps)):
            if start in seen:
                continue
            component = [start]; seen.add(start)
            for i in component:
                for row in rows:
                    j = row[i]
                    if j not in seen:
                        component.append(j); seen.add(j)
            result.append(sorted(component))
        return result

    old, new = components(baseline), components(tables)
    merges = []
    for component in new:
        parts = [part for part in old if part[0] in component]
        if len(parts) <= 1:
            continue
        frozen = [part[0] for part in parts if len(part) == 1]
        rows = [Counter(t[i] for t in tables) for i in frozen]
        exits = [sum(count for j, count in row.items() if j not in frozen) for row in rows]
        if not exits or len(set(exits)) != 1 or not exits[0]:
            raise ValueError('merged frozen states lack a common positive escape probability')
        merges.append({'previous_components': parts, 'combined_component': component,
                       'previously_frozen_states': frozen,
                       'frozen_transition_counts': [[list(x) for x in sorted(row.items())] for row in rows],
                       'escape_probability': [exits[0], len(tables)],
                       'mean_first_escape_attempts': [len(tables), exits[0]]})
    transition_counts = [Counter(t[i] for t in tables) for i in range(len(union.reps))]
    if any(transition_counts[i][j] != transition_counts[j][i]
           for i in range(len(union.reps)) for j in range(len(union.reps))):
        raise ValueError('inverse-paired patch bank failed exact detailed balance')
    control = baseline+[list(range(len(union.reps))) for _ in range(len(tables)-len(baseline))]
    q = [[experiment.p.charge(union.lift(rep))[f] for f in experiment.e.fan_faces[0]] for rep in union.reps]
    evolve = lambda values, rows: [[sum(values[t[i]][f] for t in rows) for f in range(3)] for i in range(len(union.reps))]
    values = q
    # Cayley-Hamilton: annihilating the first N powers of the N-state control
    # matrix proves annihilation of its entire observable Krylov space.
    for _ in range(len(union.reps)):
        a, b = evolve(values, control), evolve(values, tables)
        if a != b:
            raise ValueError('isolated-patch means are not an all-time negative control')
        values = a
    squares = [[x*x for x in row] for row in q]
    old_second = evolve(evolve(squares, control), control)
    new_second = evolve(evolve(squares, tables), tables)
    second_difference = [[b-a for a, b in zip(u, v)] for u, v in zip(old_second, new_second)]
    if not any(any(row) for row in second_difference):
        raise ValueError('isolated-patch fluctuation signal unexpectedly vanished')
    return {'scope': 'Boundary-contained three-face read graph; patch clock differs from full-mesh clock.',
            'gauge_states': len(union.reps), 'baseline_operators': len(baseline),
            'candidate_operators': len(tables), 'channels': [list(c) for c in channels],
            'transition_tables': tables, 'baseline_components': old, 'candidate_components': new,
            'component_mergers': merges, 'cpp_all_kernel_entries_verified': True,
            'symmetric_transition_counts': True,
            'matched_clock_mean_equal_all_times': True,
            'mean_krylov_powers_checked': len(union.reps),
            'two_step_square_charge_difference_numerator': second_difference,
            'two_step_square_charge_denominator': len(tables)**2}


def homogeneous_state(experiment, links):
    """All-q=1 freezing is exactly reduction to a parallel reflection section."""
    e, g = experiment.e, experiment.e.geometry.group
    if experiment.p.charge(links) != [1]*experiment.geometry.faces:
        raise ValueError('parallel-reflection criterion requires every face charge to be one')
    original = len(activity_mask(experiment.geometry, links))*12
    elastic = 2*sum(experiment.tables[1][e.factor.oracle.code(links, 0, p)] !=
                    e.factor.oracle.code(links, 0, p) for p in range(experiment.supports))
    stars = [set() for _ in range(e.geometry.size)]
    for face, path in zip(e.geometry.faces, e.geometry.paths):
        for i, vertex in enumerate(face[:-1]):
            stars[vertex].add(e.factor.oracle.fan.transport(links, path[i:]+path[:i]))
    section = [next(iter(star)) for star in stars] if all(len(star) == 1 for star in stars) else None
    if section is not None and any(g.conj[value][section[u]] != section[v]
                                   for (u, v), value in zip(e.geometry.edges, links)):
        raise ValueError('constant vertex-star holonomies failed edge-parallel transport')
    image = subgroup(g, e.forest.based_loops(links)[0])
    frozen = original+elastic == 0
    if frozen != (section is not None) or frozen != (len(image) == 2):
        raise ValueError('homogeneous frozen-state reduction theorem failed')
    return {'original_changing_operators': original, 'elastic_changing_operators': elastic,
            'holonomy_group_order': len(image), 'candidate_frozen': frozen,
            'parallel_reflection_section': section}


def homogeneous_census(bank, elastic):
    records = []
    for side in (3, 4, 6):
        experiment = ElasticExperiment(side, bank, elastic)
        rows, counts = [], Counter()
        for directions in itertools.product(range(6), repeat=3):
            links = constant_connection(experiment.p, directions)
            if experiment.p.charge(links) != [1]*experiment.geometry.faces:
                continue
            row = homogeneous_state(experiment, links)
            rows.append({'direction_links': list(directions), **row})
            counts[row['holonomy_group_order'], not row['original_changing_operators'], row['candidate_frozen']] += 1
        # Original full-S3 frozen control and reducible frozen control. Scan
        # every new-bank operator followed by its actual inverse in C++.
        inputs = [constant_connection(experiment.p, directions) for directions in ((1, 2, 5), (1, 1, 1))]
        scan = [encoded for op in range(experiment.operator_count)
                for encoded in (op, experiment.inverse_operator(op))]
        runs = experiment.compiled(inputs, scan, stride=len(scan))['runs'][1::2]
        for raw, run in zip(inputs, runs):
            row = homogeneous_state(experiment, raw)
            if len(run['events']) != 2*(row['original_changing_operators']+row['elastic_changing_operators']):
                raise ValueError('compiled homogeneous activity differs from theorem diagnostics')
        records.append({'side': side, 'constant_direction_assignments': 216,
                        'all_charge_one_assignments': len(rows),
                        'classification': [{'holonomy_group_order': h, 'original_frozen': old,
                                            'candidate_frozen': new, 'count': count}
                                           for (h, old, new), count in sorted(counts.items())],
                        'rows': rows, 'cpp_full_operator_scans_verified': True})
    return {'scope': 'On the supplied connected triangular torus, every face q=1: frozen iff the full based-loop image has order two. Other charge profiles are not covered.',
            'records': records}


def audit():
    bank = json.loads(Path('data/triangle-feedback.json').read_text())['bank']
    e = ResponseGeometry(3, bank).e
    census = pair_census(e.geometry.group, e.charges)
    witness = encounter_witness(bank, census['elastic_tables'])
    records = []
    for side in (3, 6, 12):
        experiment = ElasticExperiment(side, bank, census['elastic_tables'])
        reference, rng = AxialReference(side, bank, 4), random.Random(99203+side)
        initial = reference.sample(rng, 'nonabelian_reflections')['links']
        dense = [rng.randrange(6) for _ in experiment.e.geometry.edges]
        schedule = [rng.randrange(experiment.operator_count) for _ in range(4000)]
        formerly_frozen = constant_connection(experiment.p, (1, 2, 5))
        result = experiment.compiled([initial, dense, formerly_frozen], schedule, stride=400)
        persistence = []
        for run in result['runs'][4:]:
            # An exact certificate along the trace, not a timeout-based claim
            # of irreducibility. The proof in the documentation covers states
            # beyond this finite trajectory when Q=F and the image is S3.
            for raw in [formerly_frozen]+[event[-1] for event in run['events']]:
                q = experiment.p.charge(raw)
                if sum(q) != experiment.geometry.faces or len(subgroup(experiment.e.geometry.group, experiment.e.forest.based_loops(raw)[0])) != 6:
                    raise ValueError('full-S3 Q=F trajectory left its invariant sector')
                if 0 in q:
                    enabled = any((q[a] == 0) != (q[b] == 0) for a, b in experiment.e.factor.pairs)
                else:
                    enabled = any(experiment.tables[1][experiment.e.factor.oracle.code(raw, 0, p)] !=
                                  experiment.e.factor.oracle.code(raw, 0, p) for p in range(experiment.supports))
                if not enabled:
                    raise ValueError('claimed nonabsorbing state has no certified pair update')
            persistence.append({'mode': run['mode'], 'states_checked': 1+len(run['events']),
                                'full_S3_Q_equals_F_and_enabled_pair_throughout': True})
        records.append({'side': side, 'attempts': len(schedule),
                        'nonabsorption_checks': persistence,
                        'runs': [{'condition': ('Q4_reference', 'uniform_raw_links', 'formerly_frozen_full_S3')[r['condition']],
                                  'mode': r['mode'], 'events_by_rule': dict(sorted(Counter(str(v[1]) for v in r['events']).items())),
                                  'raw_events_sha256': hashlib.sha256(json.dumps(r['events'], separators=(',', ':')).encode()).hexdigest(),
                                  'final_links': r['final_links'], 'exact_link_inverse': r['exact_link_inverse']}
                                 for r in result['runs']]})
    return {'scope': 'Separate inverse-paired elastic candidate bank; fixed supplied triangular torus; no quantum amplitudes or fitted forces.',
            'pair_census': census, 'encounter_witness': witness,
            'automatic_matched_clock_response': automatic_response(bank, census['elastic_tables'], witness),
            'patch_activation': patch_activation(bank, census['elastic_tables']),
            'homogeneous_frozen_reduction': homogeneous_census(bank, census['elastic_tables']),
            'uniform_attempt_runs': records}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=Path('out/triangle-elastic-scattering.json'))
    args = parser.parse_args()
    result = audit()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k: v for k, v in result['pair_census'].items() if k not in ('tables', 'elastic_tables')}, indent=2))
    print('Wrote', args.output)


if __name__ == '__main__':
    main()
