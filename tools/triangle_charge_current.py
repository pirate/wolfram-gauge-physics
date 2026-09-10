#!/usr/bin/env python3
"""Derived scalar-charge continuity, exact drift, and gauge-dependent fluctuations."""
import argparse
import itertools
import json
import random
from collections import Counter
from fractions import Fraction
from pathlib import Path

from triangle_feedback import FeedbackExperiment


def outer(vector):
    return [[a*b for b in vector] for a in vector]


def reaction_moments(values, charges):
    """Sums over all twelve rules, derived independently of their lookup tables."""
    q = [charges[x] for x in values]
    if sorted(q) == [0, 1, 2]:
        # Four inverse channels have the same increment d; the other eight idle.
        d = [1-x for x in q]
        return [4*x for x in d], [[4*x for x in row] for row in outer(d)]
    active = q == [1, 1, 1] and ((values[0] == values[1]) != (values[1] == values[2]))
    # At an active RRR input, each permutation of (-1,0,1) occurs twice.
    # Their sum is zero and their outer-product sum is 4(3I-J). Cancellation
    # needs equal channel rates: charge-only drift is not charge-process closure.
    return [0]*3, [[4*(3*int(i == j)-1) if active else 0 for j in range(3)] for i in range(3)]


def local_census(bank):
    e = FeedbackExperiment(3, bank)
    n, rows = e.geometry.group.n, []
    for values in itertools.product(range(n), repeat=3):
        code = (values[0]*n+values[1])*n+values[2]
        q = [e.charges[x] for x in values]
        deltas, currents = [], []
        for table in e.tables[1:]:
            target = table[code]
            after = [e.charges[target//n**2], e.charges[target//n % n], e.charges[target % n]]
            d = [b-a for a, b in zip(q, after)]
            j = [-d[0], d[2]]
            if sum(d) or [-j[0], j[0]-j[1], j[1]] != d:
                raise ValueError('local fan does not satisfy exact continuity')
            deltas.append(d); currents.append(j)
        drift = [sum(d[i] for d in deltas) for i in range(3)]
        second = [[sum(d[i]*d[j] for d in deltas) for j in range(3)] for i in range(3)]
        if (drift, second) != reaction_moments(values, e.charges):
            raise ValueError('derived reaction moment formula disagrees with the full bank')
        rows.append({'input': values, 'charges': q, 'charge_drift_sum': drift,
                     'charge_raw_second_sum': second,
                     'current_drift_sum': [sum(v[i] for v in currents) for i in range(2)],
                     'current_raw_second_sum': [[sum(v[i]*v[j] for v in currents) for j in range(2)] for i in range(2)]})
    return rows


class CurrentProbe:
    def __init__(self, side, bank):
        self.e = FeedbackExperiment(side, bank)
        self.faces = len(self.e.geometry.faces)
        self.dual = sorted({tuple(sorted(p)) for p in self.e.factor.pairs})
        self.edge_ids = {edge: i for i, edge in enumerate(self.dual)}
        face_edges = [{self.e.geometry.ids[tuple(sorted((a, b)))] for a, b in zip(f, f[1:])}
                      for f in self.e.geometry.faces]
        shared = [face_edges[a] & face_edges[b] for a, b in self.dual]
        if any(len(edges) != 1 for edges in shared):
            raise ValueError('dual current edge does not cross exactly one primal link')
        self.dual_primal_edges = [next(iter(edges)) for edges in shared]
        self.supports = [self.e.factor.pairs, [f[::-1] for f in self.e.fan_faces]]
        self.paths = []
        for supports in self.supports:
            paths = []
            for face_ids in supports:
                paths.append([(self.edge_ids[tuple(sorted((a, b)))], 1 if a < b else -1)
                              for a, b in zip(face_ids, face_ids[1:])])
            self.paths.append(paths)
        for arity, paths in enumerate(self.paths):
            for patch, path in enumerate(paths):
                actual = {self.dual_primal_edges[edge] for edge, _ in path}
                writes = (self.e.factor.oracle.fan.writes[patch] if arity else
                          {self.e.factor.oracle.pairs[patch][0][0][0]})
                if actual != writes:
                    raise ValueError('current path does not coincide with primitive written links')

    def charge(self, links):
        return [self.e.charges[x] for x in self.e.geometry.holonomies(links)]

    def empty(self):
        return [Counter(), Counter(), Counter(), Counter()]

    def add(self, accumulator, arity, patch, drift, second):
        qd, qs, jd, js = accumulator
        faces, path = self.supports[arity][patch], self.paths[arity][patch]
        for i, f in enumerate(faces):
            qd[f] += drift[i]
            for j, g in enumerate(faces):
                if second[i][j]:
                    qs[f, g] += second[i][j]
        # For a dual path, continuity uniquely fixes J=(-dq_0, dq_2).
        # A pair has only J=-dq_0. No circulation or routing choice is inserted.
        rows = [(0, -1)] if arity == 0 else [(0, -1), (2, 1)]
        for i, ((edge, sign), (a, s)) in enumerate(zip(path, rows)):
            jd[edge] += sign*s*drift[a]
            for (other, orient), (b, t) in zip(path, rows):
                value = sign*s*orient*t*second[a][b]
                if value:
                    js[edge, other] += value

    def finish(self, accumulator):
        qd, qs, jd, js = accumulator
        projected_drift, projected_second = Counter(), Counter()
        for edge, value in jd.items():
            a, b = self.dual[edge]
            projected_drift[a] -= value; projected_drift[b] += value
        for (edge, other), value in js.items():
            for a, s in zip(self.dual[edge], (-1, 1)):
                for b, t in zip(self.dual[other], (-1, 1)):
                    projected_second[a, b] += s*t*value
        clean = lambda c: {key: value for key, value in c.items() if value}
        if clean(qd) != clean(projected_drift) or clean(qs) != clean(projected_second):
            raise ValueError('global current moments violate continuity')
        row_sums = Counter()
        for (a, _), value in qs.items():
            row_sums[a] += value
        if sum(qd.values()) or any(row_sums.values()):
            raise ValueError('global charge drift or fluctuations fail conservation')
        sparse = lambda c: [list(key)+[value] for key, value in sorted(c.items()) if value]
        return {'charge_drift_sum': [qd[f] for f in range(self.faces)],
                'charge_raw_second_sum': sparse(qs),
                'current_drift_sum': [jd[e] for e in range(len(self.dual))],
                'current_raw_second_sum': sparse(js),
                'attempted_operators': len(self.e.tables)*len(self.e.factor.pairs),
                'exact_continuity_checked': True}

    def analytic(self, links):
        e = self.e
        q = self.charge(links)
        result = self.empty()
        for patch, faces in enumerate(self.supports[0]):
            a, b = [q[f] for f in faces]
            delta = [b-a, a-b] if (a == 0) != (b == 0) else [0, 0]
            self.add(result, 0, patch, delta, outer(delta))
        for patch in range(len(self.supports[1])):
            code = e.factor.oracle.code(links, 1, patch)
            n = e.geometry.group.n
            values = [code//n**2, code//n % n, code % n]
            drift, second = reaction_moments(values, e.charges)
            self.add(result, 1, patch, drift, second)
        return self.finish(result)

    def enumerate(self, links, compiled=True):
        e = self.e
        q = self.charge(links)
        result, targets, changing = self.empty(), {}, 0
        supports = len(e.factor.pairs)
        for rule, table in enumerate(e.tables):
            arity = int(bool(rule))
            for patch, faces in enumerate(self.supports[arity]):
                code = e.factor.oracle.code(links, rule, patch)
                if code == table[code]:
                    continue
                moved = links[:]
                e.factor.oracle.update(moved, rule, patch, table)
                after = self.charge(moved)
                delta = [after[f]-q[f] for f in faces]
                if sum(delta) or any(a != b for f, (a, b) in enumerate(zip(q, after)) if f not in faces):
                    raise ValueError('actual primitive has charge flow outside its dual path')
                self.add(result, arity, patch, delta, outer(delta))
                targets.setdefault(tuple(moved), rule*supports+patch)
                changing += 1
        if compiled:
            scan = [encoded for encoded in range(len(e.tables)*supports) for _ in range(2)]
            run = e.compiled_bank([links], scan, len(scan), capture_links=True)['runs'][1]
            if len(run['events']) != 2*changing:
                raise ValueError('full current operator/inverse scan differs')
            if {tuple(state) for state in run['raw_event_links'][::2]} != set(targets):
                raise ValueError('C++ current targets differ from Python raw updates')
        return self.finish(result)

    def projected(self, links, weights):
        if len(weights) != self.faces or any(type(w) is not int for w in weights):
            raise ValueError('projection needs an integer weight for every face')
        moments = self.analytic(links)
        drift = sum(w*d for w, d in zip(weights, moments['charge_drift_sum']))
        second = sum(weights[a]*weights[b]*v for a, b, v in moments['charge_raw_second_sum'])
        denominator = moments['attempted_operators']
        if denominator*second < drift*drift:
            raise ValueError('conditional projected variance is negative')
        return drift, second, denominator


def readout_audit(source):
    probe = CurrentProbe(12, source['bank'])
    if json.loads(json.dumps(probe.e.readout())) != source['readout']:
        raise ValueError('raw braid-memory preparations do not reproduce')
    links = [r['checked']['final_links'] for r in source['readout']['preparations']]
    records = []
    for i, state in enumerate(links):
        predicted = probe.analytic(state)
        if predicted != probe.enumerate(state):
            raise ValueError('global analytic current moments differ from raw generator enumeration')
        records.append({'state': i, 'links': state, 'moments': predicted,
                        'all_unique_raw_targets_cpp_checked': True, 'all_operator_inverse_scan_checked': True})
    if len({tuple(probe.charge(s)) for s in links}) != 1:
        raise ValueError('readout does not hold the full charge field fixed')
    if len({tuple(r['moments']['charge_drift_sum']) for r in records}) != 1:
        raise ValueError('readout has hidden-state-dependent mean charge drift')
    a, b = source['readout']['global_rate_witness']['states']
    matrices = [{(u, v): n for u, v, n in records[i]['moments']['charge_raw_second_sum']} for i in (a, b)]
    difference = [[u, v, matrices[0].get((u, v), 0)-matrices[1].get((u, v), 0)]
                  for u, v in sorted(set(matrices[0]) | set(matrices[1]))
                  if matrices[0].get((u, v), 0) != matrices[1].get((u, v), 0)]
    if not difference:
        raise ValueError('prepared gauge-memory states have no current-noise witness')
    g, e = probe.e.geometry.group, probe.e
    gauge = []
    for i in (a, b):
        rng = random.Random(34071+i)
        frames = [rng.randrange(g.n) for _ in range(e.geometry.size)]
        transformed = [g.mul[frames[v]][g.mul[x][g.inv[frames[u]]]] for (u, v), x in zip(e.geometry.edges, links[i])]
        if probe.analytic(transformed) != records[i]['moments'] or probe.enumerate(transformed) != records[i]['moments']:
            raise ValueError('local gauge frames change the current moments')
        gauge.append({'state': i, 'frames': frames, 'exact_same_moments': True})
    return {'records': records, 'dual_edges': probe.dual, 'same_charge_field': probe.charge(links[0]),
            'same_boundary_witness_states': [a, b],
            'same_boundary_product': source['readout']['encounter_products'][a],
            'raw_second_sum_difference': difference, 'gauge_checks': gauge}


def evolution_audit(source):
    probe = CurrentProbe(6, source['bank'])
    e = probe.e
    # A specified half-torus cut, not an inferred Euclidean coordinate system.
    weights = [int((f//2) % 6 < 3) for f in range(probe.faces)]
    initial = source['evolution']['inputs'][2]
    records = []
    for source_run in source['evolution']['runs']:
        checked = e.compiled_bank([initial], source_run['schedule'], 1000)['runs'][1]
        if checked != dict(source_run['checked']['runs'][-1], condition=0):
            raise ValueError('source current history does not reproduce')
        state, previous = initial[:], 0
        x0 = sum(w*q for w, q in zip(weights, probe.charge(state)))
        drift_sum = variance_sum = residual_square_sum = raw_square = raw_square_expected_sum = 0
        events = checked['events']+[[source['evolution']['attempts'], None, None, None, None]]
        intervals, total_flux = [], Counter()
        for tick, rule, patch, code, target in events:
            if tick == previous:
                continue
            drift, second, denominator = probe.projected(state, weights)
            old_q = probe.charge(state)
            x = sum(w*q for w, q in zip(weights, old_q))
            duration = tick-previous
            if rule is not None:
                if e.factor.oracle.update(state, rule, patch, e.tables[rule]) != (code, target):
                    raise ValueError('current history misses an actual primitive')
                after = probe.charge(state)
                arity = int(bool(rule))
                faces, path = probe.supports[arity][patch], probe.paths[arity][patch]
                d = [after[f]-old_q[f] for f in faces]
                js = [-d[0]] if not arity else [-d[0], d[2]]
                global_delta = Counter()
                for (edge, sign), j in zip(path, js):
                    total_flux[edge] += sign*j
                    a, b = probe.dual[edge]
                    global_delta[a] -= sign*j; global_delta[b] += sign*j
                if [global_delta[f] for f in range(probe.faces)] != [b-a for a, b in zip(old_q, after)]:
                    raise ValueError('realized microscopic continuity fails')
                dx = sum(w*(b-a) for w, a, b in zip(weights, old_q, after))
            else:
                dx = 0
            drift_sum += duration*drift
            variance_sum += duration*(denominator*second-drift**2)
            residual_square_sum += duration*drift**2+denominator**2*dx**2-2*denominator*dx*drift
            raw_square += dx**2
            raw_square_expected_sum += duration*second
            intervals.append([previous, tick, x, dx, drift, second])
            previous = tick
        final_q = probe.charge(state)
        initial_q = probe.charge(initial)
        projected = Counter()
        for edge, flux in total_flux.items():
            a, b = probe.dual[edge]
            projected[a] -= flux; projected[b] += flux
        if state != checked['final_links'] or [projected[f] for f in range(probe.faces)] != [b-a for a, b in zip(initial_q, final_q)]:
            raise ValueError('integrated current does not recover the final charge field')
        final_x = sum(w*q for w, q in zip(weights, final_q))
        records.append({'seed': source_run['seed'], 'initial_region_charge': x0,
                        'final_region_charge': final_x, 'attempts': source['evolution']['attempts'],
                        'denominator': denominator, 'drift_compensator_numerator': drift_sum,
                        'martingale_residual_numerator': denominator*(final_x-x0)-drift_sum,
                        'predictable_variance_numerator_over_denominator_squared': variance_sum,
                        'realized_centered_square_numerator_over_denominator_squared': residual_square_sum,
                        'realized_raw_squared_changes': raw_square,
                        'raw_square_compensator_numerator': raw_square_expected_sum,
                        'constant_state_intervals_start_end_X_deltaX_driftSum_secondSum': intervals,
                        'integrated_dual_currents': [[edge, n] for edge, n in sorted(total_flux.items()) if n],
                        'raw_link_replay_and_integrated_continuity_checked': True})
        print('current run', source_run['seed'], 'residual numerator', records[-1]['martingale_residual_numerator'], flush=True)
    return {'weights': weights, 'dual_edges': probe.dual, 'runs': records,
            'scope': 'Exact predictable drift/variance and realized currents for a specified half-torus observable. Four pilot histories are not equilibrium evidence, a continuum limit, or a calibrated statistical discovery test.'}


def homogeneous_activation(source):
    """Same uniform charge field, different microscopic reflection connections."""
    probe = CurrentProbe(3, source['bank'])
    e, g = probe.e, probe.e.geometry.group
    reflections = [x for x, q in enumerate(e.charges) if q == 1]
    patch = 0
    support = sorted({edge for path in e.factor.oracle.fan.patches[patch] for edge, _ in path})
    if len(support) != 7:
        raise ValueError('homogeneous activation census expects seven fan links')
    triples = Counter()
    for assignment in itertools.product(reflections, repeat=len(support)):
        state = [reflections[0]]*len(e.geometry.edges)
        for edge, value in zip(support, assignment):
            state[edge] = value
        triples[e.factor.oracle.code(state, 1, patch)] += 1
    if len(triples) != 27 or set(triples.values()) != {81}:
        raise ValueError('independent reflection links do not induce uniform based triples')
    active_assignments = sum(n for code, n in triples.items()
                             if reaction_moments([code//g.n**2, code//g.n % g.n, code % g.n], e.charges)[1][0][0])
    probability = Fraction(active_assignments, 3**7)*Fraction(len(e.tables)-1, len(e.tables))
    contrast = Fraction(sum(count*sum(reaction_moments(
        [code//g.n**2, code//g.n % g.n, code % g.n], e.charges)[1][i][i] for i in range(3))
        for code, count in triples.items()), 3**7*len(e.tables))
    runs = []
    for side, seed in [(6, seed) for seed in range(4)]+[(12, 0)]:
        probe = CurrentProbe(side, source['bank'])
        e = probe.e
        frozen = [reflections[0]]*len(e.geometry.edges)
        rng = random.Random(44100+seed)
        active = [rng.choice(reflections) for _ in e.geometry.edges]
        if probe.charge(frozen) != [1]*probe.faces or probe.charge(active) != [1]*probe.faces:
            raise ValueError('initial reflection connections do not share uniform charge')
        initial_moments = probe.analytic(active)
        frozen_moments = probe.analytic(frozen)
        if any(initial_moments['charge_drift_sum']) or any(initial_moments['current_drift_sum']):
            raise ValueError('homogeneous reflection field has nonzero mean drift')
        if frozen_moments['charge_raw_second_sum'] or not initial_moments['charge_raw_second_sum']:
            raise ValueError('activation controls do not distinguish zero from nonzero noise')
        if initial_moments != probe.enumerate(active) or frozen_moments != probe.enumerate(frozen):
            raise ValueError('homogeneous initial generator moments do not reproduce')
        attempts = 10000*(side//6)**2  # Match attempted updates per rooted operator.
        rng = random.Random(88000+seed)
        schedule = [rng.randrange(len(e.tables)*len(e.factor.pairs)) for _ in range(attempts)]
        checked = e.compiled_bank([frozen, active], schedule, 1000)
        for control in checked['runs'][:3]:
            original = frozen if control['condition'] == 0 else active
            if control['events'] or control['final_links'] != original:
                raise ValueError('frozen or reaction-disabled control has changed')
        run = checked['runs'][3]
        if not run['events'] or run['events'][0][1] == 0:
            raise ValueError('reaction did not initiate the active evolution')
        state, event_index, flux = active[:], 0, Counter()
        snapshots = []
        for tick in (0, 20, 100, attempts//10, attempts):
            while event_index < len(run['events']) and run['events'][event_index][0] <= tick:
                _, rule, patch, code, target = run['events'][event_index]
                before = probe.charge(state)
                if e.factor.oracle.update(state, rule, patch, e.tables[rule]) != (code, target):
                    raise ValueError('homogeneous raw-link snapshot replay disagrees')
                after = probe.charge(state)
                arity = int(bool(rule))
                ids = probe.supports[arity][patch]
                d = [after[f]-before[f] for f in ids]
                current = [-d[0]] if not arity else [-d[0], d[2]]
                for (edge, sign), j in zip(probe.paths[arity][patch], current):
                    flux[edge] += sign*j
                event_index += 1
            q = probe.charge(state)
            divergence = Counter()
            for edge, value in flux.items():
                a, b = probe.dual[edge]
                divergence[a] -= value; divergence[b] += value
            if q != [1+divergence[f] for f in range(probe.faces)]:
                raise ValueError('activation current does not reconstruct the charge field')
            population = [q.count(i) for i in range(3)]
            if population[0] != population[2] or sum(q) != probe.faces:
                raise ValueError('homogeneous activation violates its conserved charge')
            snapshots.append({'tick': tick, 'links': state[:], 'charge_field': q,
                              'populations_0_1_2': population,
                              'charge_contrast_norm_squared': sum((x-1)**2 for x in q),
                              'neighbor_charge_contrast_product_sum': sum((q[a]-1)*(q[b]-1) for a, b in probe.dual),
                              'integrated_dual_currents': [[edge, n] for edge, n in sorted(flux.items()) if n]})
        if state != run['final_links']:
            raise ValueError('activation snapshots miss the C++ final connection')
        runs.append({'side': side, 'seed': seed, 'attempts': attempts,
                     'initial_frozen_links': frozen, 'initial_active_links': active,
                     'initial_active_moments': initial_moments,
                     'frozen_and_transport_only_controls_exactly_fixed': True,
                     'schedule': schedule, 'checked': checked, 'snapshots': snapshots})
        print('activation', side, seed, 'final populations', snapshots[-1]['populations_0_1_2'], flush=True)
    return {'local_fan_links': support, 'local_reflection_assignments': 3**7,
            'based_triple_multiplicities': sorted(triples.items()), 'active_assignments': active_assignments,
            'initial_ensemble_changing_probability_per_attempt': [probability.numerator, probability.denominator],
            'initial_ensemble_expected_charge_contrast_increment': [contrast.numerator, contrast.denominator],
            'runs': runs,
            'scope': 'Uniform initial scalar charge with specified microscopic gauge disorder. Classical fluctuation-activated dynamics, not spontaneous symmetry breaking from identical raw states, an absorbing phase transition, equilibrium sampling, or emergent matter.'}


def audit():
    source = json.loads(Path('data/triangle-feedback.json').read_text())
    return {'scope': 'Derived conserved scalar-charge currents and classical fluctuations; not electromagnetic charge, a quantum wave law, emergent geometry, or matter.',
            'local_census': local_census(source['bank']),
            'readout': readout_audit(source), 'evolution': evolution_audit(source),
            'homogeneous_activation': homogeneous_activation(source)}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, default=Path('out/triangle-charge-current.json'))
    args = parser.parse_args()
    result = audit()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, separators=(',', ':'))+'\n')


if __name__ == '__main__':
    main()
