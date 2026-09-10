#!/usr/bin/env python3
"""Twelve explicit relational reactions and their stationary shared-link activity.

The family reconstructs the existing triangle bank; it does not independently
select a physical matching law. No charge value enters rule construction.
"""
import itertools
import json
import math
import random
import statistics
from collections import Counter
from fractions import Fraction
from pathlib import Path

from cycle_charge_compatibility import cycle_group
from run_three_face import LinkOracle
from screen_braid_rules import nullspace


class RelationalRule:
    def __init__(self, group, layout, orientation):
        self.g, self.layout, self.orientation = group, tuple(layout), orientation
        self.kind = ['E' if h == group.identity else 'R' if group.inv[h] == h else 'Z'
                     for h in range(group.n)]

    def reverse(self, w):
        return tuple(self.g.inv[h] for h in w[::-1])

    def forward_left(self, r, s):
        g = self.g
        z = g.mul[r][s] if self.orientation == 1 else g.mul[s][r]
        h = g.mul[g.inv[z]][s] if self.layout.index('Z') < self.layout.index('R') else g.mul[s][g.inv[z]]
        entries = {'E': g.identity, 'R': h, 'Z': z}
        return tuple(entries[k] for k in self.layout)

    def inverse_left(self, w):
        g = self.g
        s = g.mul[g.mul[w[0]][w[1]]][w[2]]
        z = w[self.layout.index('Z')]
        r = g.mul[z][s] if self.orientation == 1 else g.mul[s][z]
        return r, r, s

    def apply(self, w):
        kinds = tuple(self.kind[h] for h in w)
        a, b, c = w
        if kinds == ('R', 'R', 'R'):
            if a == b and b != c:
                return self.forward_left(a, c)
            if b == c and a != b:
                return self.reverse(self.forward_left(c, a))
        elif kinds == self.layout:
            return self.inverse_left(w)
        elif kinds == self.layout[::-1]:
            return self.reverse(self.inverse_left(self.reverse(w)))
        return w

    def __getitem__(self, code):
        n = self.g.n
        w = (code//(n*n), (code//n) % n, code % n)
        a, b, c = self.apply(w)
        return (a*n+b)*n+c


class TransportRule:
    def __init__(self, group, kind):
        self.g, self.kind = group, kind

    def apply(self, w):
        a, b, c = w
        g = self.g
        if self.kind == 0:
            return (b, a, c) if g.identity in (a, b) else w
        if self.kind == 1:
            return g.mul[g.mul[a][b]][g.inv[a]], a, c
        return b, g.mul[g.mul[g.inv[b]][a]][b], c

    def __getitem__(self, code):
        n = self.g.n
        a, b, c = self.apply((code//(n*n), (code//n) % n, code % n))
        return (a*n+b)*n+c


def family(g):
    return [RelationalRule(g, layout, orientation)
            for layout in itertools.permutations(('E', 'R', 'Z')) for orientation in (1, -1)]


def exact_family(n):
    g = cycle_group(n)
    rules = family(g)
    states = list(itertools.product(range(g.n), repeat=3))
    ids = {w: i for i, w in enumerate(states)}
    labels = sorted(set(g.sectors)-{g.identity})
    q = [0 if h == g.identity else 1 if g.inv[h] == h else 2 for h in range(g.n)]
    equations, tables, forward, backward = set(), [], 0, 0
    for rule in rules:
        table, moved = [], 0
        for w in states:
            target = rule.apply(w)
            table.append(ids[target])
            assert rule.apply(target) == w
            assert rule.apply(rule.reverse(w)) == rule.reverse(target)
            assert g.mul[g.mul[w[0]][w[1]]][w[2]] == g.mul[g.mul[target[0]][target[1]]][target[2]]
            if target != w:
                moved += 1
                equations.add(tuple(sum(g.sectors[h] == s for h in w)-sum(g.sectors[h] == s for h in target) for s in labels))
                if all(q[h] == 1 for h in w):
                    forward += 1
                else:
                    backward += 1
        assert moved == 4*n*(n-1)
        tables.append(table)
    assert len({tuple(table) for table in tables}) == 12
    basis = nullspace(sorted(equations), len(labels))
    assert basis == [[q[h] for h in labels]]
    for w in states:
        active = sum(rule.apply(w) != w for rule in rules)
        expected = 12 if all(q[h] == 1 for h in w) and ((w[0] == w[1]) != (w[1] == w[2])) else 4 if sorted(q[h] for h in w) == [0, 1, 2] else 0
        assert active == expected
    result = {'cycle_vertices': n, 'rule_count': 12, 'moved_tuples_per_rule': 4*n*(n-1),
              'derived_charge_basis': basis, 'nonidentity_class_representatives': labels,
              'uniform_input_forward_probability_per_reaction_proposal': str(Fraction(forward, 12*g.n**3)),
              'uniform_input_backward_probability_per_reaction_proposal': str(Fraction(backward, 12*g.n**3))}
    if n == 3:
        bank = json.loads(Path('data/triangle-feedback.json').read_text())['bank']
        assert [list(p) for p in g.elements] == bank['group_automorphisms']
        result['existing_triangle_rule_mapping'] = [
            {'layout': rule.layout, 'rotation_orientation': rule.orientation,
             'existing_rule_id': bank['selected_rule_ids'][bank['tables'].index(table)]}
            for rule, table in zip(rules, tables)]
        assert set(map(tuple, tables)) == set(map(tuple, bank['tables']))
    return result


def embedding(source_n, target_n):
    a, b = cycle_group(source_n), cycle_group(target_n)
    ra, rb = family(a), family(b)
    d = target_n//source_n
    mapping = []
    for p in a.elements:
        sign = 1 if (p[1]-p[0]) % source_n == 1 else -1
        mapping.append(b.elements.index(tuple((sign*j+d*p[0]) % target_n for j in range(target_n))))
    assert len(set(mapping)) == a.n
    for w in itertools.product(range(a.n), repeat=3):
        lifted = tuple(mapping[h] for h in w)
        for x, y in zip(ra, rb):
            assert tuple(mapping[h] for h in x.apply(w)) == y.apply(lifted)
    return {'source_cycle': source_n, 'target_cycle': target_n, 'all_rule_tuple_comparisons': 12*a.n**3}


def mean_se(values):
    return {'mean': statistics.mean(values), 'standard_error_across_independent_trajectories': statistics.stdev(values)/len(values)**0.5}


def shared_link_stationary(n, trials=64, attempts=4000, side=4):
    g = cycle_group(n)
    oracle = LinkOracle(side, g)
    rules = [TransportRule(g, k) for k in range(3)]+family(g)
    q = [0 if h == g.identity else 1 if g.inv[h] == h else 2 for h in range(g.n)]
    populations = lambda raw: Counter(q[oracle.transport(raw, face)] for face in oracle.face_paths)
    rows, witness = [], None
    for trial in range(trials):
        rng = random.Random(252218000+n*1000+trial)
        initial = [rng.randrange(g.n) for _ in oracle.edges]
        raw = initial[:]
        before = populations(raw)
        counts = Counter()
        for step in range(attempts):
            channel, patch = rng.randrange(15), rng.randrange(len(oracle.patches))
            counts['reaction_proposals' if channel >= 3 else 'transport_proposals'] += 1
            rule = rules[channel]
            word = tuple(oracle.transport(raw, path) for path in oracle.patches[patch][::-1])
            target = rule.apply(word)
            if target == word:
                continue
            save = raw[:] if channel >= 3 and witness is None else None
            oracle.update(raw, patch, rule)
            if channel < 3:
                counts['changing_transport_events'] += 1
                continue
            direction = 'created_rotation' if all(q[h] == 1 for h in word) else 'consumed_rotation'
            counts[direction] += 1
            if save is not None:
                assert all(save[e] == raw[e] for e in range(len(raw)) if e not in oracle.writes[patch])
                replay = raw[:]
                oracle.update(replay, patch, rule)
                assert replay == save
                witness = {'trial': trial, 'attempt': step, 'channel': channel, 'patch': patch,
                           'before_tuple': word, 'after_tuple': target,
                           'before_links': save, 'after_links': raw[:], 'exterior_fixed': True, 'raw_inverse_exact': True}
        after = populations(raw)
        assert before[1]+2*before[2] == after[1]+2*after[2]
        assert after[2]-before[2] == counts['created_rotation']-counts['consumed_rotation']
        rows.append({'trial': trial, 'initial_face_populations': dict(before), 'final_face_populations': dict(after), **counts})
    forward = mean_se([r.get('created_rotation', 0)/attempts for r in rows])
    backward = mean_se([r.get('consumed_rotation', 0)/attempts for r in rows])
    total = mean_se([(r.get('created_rotation', 0)+r.get('consumed_rotation', 0))/attempts for r in rows])
    prediction = Fraction(2*(n-1), 5*n*n)
    return {'cycle_vertices': n, 'side': side, 'faces': len(oracle.faces), 'links': len(oracle.edges),
            'independent_trajectories': trials, 'attempts_per_trajectory': attempts,
            'initial_measure': 'Independent uniform raw links; exactly stationary mixture of conserved sectors, not an equilibration claim.',
            'scheduler': 'Uniformly choose one fan and one of 15 channels: vacancy/Hurwitz/inverse-Hurwitz on the first two based loops, or one of twelve relational reactions.',
            'clock_scope': 'All proposals, including idle ones. Channel layout is this generic fan experiment, not claimed identical to the compiled triangle pair-support scheduler.',
            'forward_activity_per_all_proposal': forward, 'backward_activity_per_all_proposal': backward,
            'total_reaction_activity_per_all_proposal': total, 'exact_total_reaction_activity': str(prediction),
            'total_activity_z_score': (total['mean']-float(prediction))/total['standard_error_across_independent_trajectories'],
            'independent_trajectory_counts': rows, 'raw_reaction_witness': witness}


def coupled_refinement(source_n=3, target_n=27, trials=8, attempts=4000):
    a, b = cycle_group(source_n), cycle_group(target_n)
    oa, ob = LinkOracle(4, a), LinkOracle(4, b)
    ra = [TransportRule(a, k) for k in range(3)]+family(a)
    rb = [TransportRule(b, k) for k in range(3)]+family(b)
    d = target_n//source_n
    mapping = []
    for p in a.elements:
        sign = 1 if (p[1]-p[0]) % source_n == 1 else -1
        mapping.append(b.elements.index(tuple((sign*j+d*p[0]) % target_n for j in range(target_n))))
    records = []
    for trial in range(trials):
        rng = random.Random(252219000+trial)
        x = [rng.randrange(a.n) for _ in oa.edges]
        y = [mapping[h] for h in x]
        initial_x, initial_y = x[:], y[:]
        counts = Counter()
        for _ in range(attempts):
            channel, patch = rng.randrange(15), rng.randrange(len(oa.patches))
            wx = tuple(oa.transport(x, p) for p in oa.patches[patch][::-1])
            wy = tuple(ob.transport(y, p) for p in ob.patches[patch][::-1])
            assert wy == tuple(mapping[h] for h in wx)
            tx, ty = ra[channel].apply(wx), rb[channel].apply(wy)
            assert ty == tuple(mapping[h] for h in tx)
            if tx != wx:
                oa.update(x, patch, ra[channel])
                ob.update(y, patch, rb[channel])
                counts['reactions' if channel >= 3 else 'transport'] += 1
            assert y == [mapping[h] for h in x]
        records.append({'trial': trial, 'event_counts_identical_at_both_resolutions': dict(counts),
                        'initial_coarse_links': initial_x, 'initial_fine_links': initial_y,
                        'final_coarse_links': x, 'final_fine_links': y})
    return {'source_cycle': source_n, 'target_cycle': target_n, 'attempts_per_trajectory': attempts,
            'trajectories': trials, 'every_raw_state_intertwines': True, 'embedding_element_ids': mapping,
            'scope': 'Fine initial links lie in the embedded coarse subgroup, not the full-fiber uniform ensemble.',
            'records': records}


def product_reference(n, rho):
    a, b = (2-rho)*(n-1), (1-rho)*n
    root = math.sqrt(b*b+4*a*rho)
    z = 2*rho/(b+root) if b >= 0 else (-b+root)/(2*a)
    partition = 1+n*z+(n-1)*z*z
    p0, pr, pz = 1/partition, n*z/partition, (n-1)*z*z/partition
    forward, backward = 2*(n-1)/(n*n)*pr**3, 2*p0*pr*pz
    assert abs(forward-backward) < 1e-14 and abs(pr+2*pz-rho) < 1e-14
    return {'cycle_vertices': n, 'mean_charge': rho, 'fugacity_reference_only': z,
            'probabilities_e_R_Z': [p0, pr, pz], 'activity_per_reaction_proposal': forward+backward,
            'n_scaled_activity': n*(forward+backward), 'limiting_n_scaled_activity': 4*min(rho, 2-rho)**3}


if __name__ == '__main__':
    report = {'scope': 'Explicit reconstruction and refinement of a selected relational law; no claim of unique rule selection, physical energy, or emergent matter.',
              'exact_families': [exact_family(n) for n in (3, 5, 9)],
              'embeddings': [embedding(a, b) for a, b in ((3, 9), (3, 27), (9, 27))],
              'coupled_raw_refinement': coupled_refinement(),
              'unconditioned_product_reference_density_scaling': [product_reference(n, rho)
                  for rho in (0.25, 0.5, 1.0, 1.5) for n in (3, 9, 27, 81, 243, 729, 6561)],
              'stationary_shared_link_runs': []}
    for n in (3, 5, 9, 27, 81):
        row = shared_link_stationary(n)
        report['stationary_shared_link_runs'].append(row)
        print(json.dumps({'cycle': n, 'activity': row['total_reaction_activity_per_all_proposal'],
                          'exact': row['exact_total_reaction_activity'], 'z_score': row['total_activity_z_score']}), flush=True)
    output = Path('data/cycle-relational-dynamics.json')
    output.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({'output': str(output)}), flush=True)
