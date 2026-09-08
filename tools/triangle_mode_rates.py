#!/usr/bin/env python3
"""Exact frozen-state spectral transition kernels of the derived primitive bank.

These are one-attempt probabilities conditional on a complete raw connection,
not particle decay constants or a Markov assumption for spectral summaries.
"""
import argparse
import itertools
import json
import random
from collections import Counter
from fractions import Fraction
from pathlib import Path

from fiber_mode_probe import lifted_graph_laplacian
from triangle_feedback import FeedbackExperiment
from triangle_modes import exact_band_count, reaction_audit


def saturation_audit(bank):
    e = FeedbackExperiment(3, bank)
    g = e.geometry.group
    rows = []
    for code, target in enumerate(e.tables[0]):
        if code == target:
            continue
        before, after = divmod(code, g.n), divmod(target, g.n)
        charges = [e.charges[x] for x in before]
        moved = sum(charges)
        if charges.count(0) != 1 or [e.charges[x] for x in after] != charges[::-1]:
            raise ValueError('vacancy rule does not move one nonzero charge')
        rows.append({'input': before, 'output': after, 'moved_charge': moved})
    # Nonnegative charges make every local triple at global Q<=2 one of these.
    low_charge = []
    for values in itertools.product(range(g.n), repeat=3):
        if sum(e.charges[x] for x in values) <= 2:
            code = (values[0]*g.n+values[1])*g.n+values[2]
            if any(table[code] != code for table in e.tables[1:]):
                raise ValueError('a reaction acts at total charge at most two')
            low_charge.append(code)
    return {'vacancy_charge_transfers': rows, 'reaction_fixed_input_codes_at_charge_at_most_two': low_charge,
            'saturation_bound': 'If n_+(L-12I)=Q/2, both orientation charges equal Q/2. Every changing vacancy move of charge q leaves at most Q/2-q above-band modes.',
            'charge_two_consequence': 'At global Q=2 every reaction is fixed. If one above-band mode exists, every changing primitive move destroys it. This concerns above-band count, not all localized structures.'}


class ModeRates:
    def __init__(self, side, bank):
        self.e = FeedbackExperiment(side, bank)
        self.cache = {}

    def observable(self, links):
        key = tuple(links)
        if key not in self.cache:
            e = self.e
            charges = [e.charges[x] for x in e.geometry.holonomies(links)]
            matrix = lifted_graph_laplacian(range(e.geometry.size), e.geometry.edges,
                                          links, e.geometry.group, e.geometry.fiber_edges)
            rank = exact_band_count(matrix)['positive']
            up, down = (sum(charges[c::2]) for c in range(2))
            if rank > min(up, down):
                raise ValueError('spectral count exceeds orientation capacity')
            self.cache[key] = (rank, up, down, charges.count(1), charges.count(2))
        return self.cache[key]

    def kernel(self, links, compiled=True):
        e = self.e
        initial = self.observable(links)
        supports = len(e.factor.pairs)
        counters = [Counter(), Counter()]
        changing, raw_targets, witnesses = [0, 0], {}, []
        saturated_checks = []
        for rule, table in enumerate(e.tables):
            for patch in range(supports):
                code = e.factor.oracle.code(links, rule, patch)
                if table[code] == code:
                    counters[bool(rule)][initial] += 1
                    continue
                target = links[:]
                actual = e.factor.oracle.update(target, rule, patch, table)
                if actual != (code, table[code]):
                    raise ValueError('frozen-state update disagrees with its table')
                obs = self.observable(target)
                if sum(obs[1:3]) != sum(initial[1:3]):
                    raise ValueError('frozen-state update does not conserve charge')
                counters[bool(rule)][obs] += 1
                changing[bool(rule)] += 1
                raw_targets.setdefault(tuple(target), rule*supports+patch)
                if obs[0] != initial[0]:
                    witnesses.append([rule*supports+patch, obs[0]])
                if rule == 0 and initial[0] == initial[1] == initial[2]:
                    moved = sum(e.charges[x] for x in divmod(code, e.geometry.group.n))
                    bound = initial[0]-moved
                    if obs[0] > bound:
                        raise ValueError('saturated state evades the vacancy-loss bound')
                    saturated_checks.append([patch, moved, obs[0], bound])
        # Independent C++ final raw connections, not only target local codes or
        # histograms. Check each distinct changed target with its representative
        # operator; the full scan below checks every operator's codes and inverse.
        if compiled:
            for target, encoded in raw_targets.items():
                run = e.compiled_bank([links], [encoded], 1)['runs'][1]
                if run['final_links'] != list(target):
                    raise ValueError('independent C++ frozen target disagrees')
            scan = [encoded for encoded in range(len(e.tables)*supports) for _ in range(2)]
            run = e.compiled_bank([links], scan, len(scan))['runs'][1]
            if len(run['events']) != 2*sum(changing) or run['final_links'] != links:
                raise ValueError('full operator/inverse scan disagrees')
        total = counters[0]+counters[1]
        denominator = len(e.tables)*supports
        if sum(total.values()) != denominator:
            raise ValueError('frozen kernel omitted attempted operators')
        rank_counts = Counter()
        for state, n in total.items():
            rank_counts[state[0]] += n
        losses = [sum(n for obs, n in c.items() if obs[0] < initial[0]) for c in counters]
        gains = [sum(n for obs, n in c.items() if obs[0] > initial[0]) for c in counters]
        return {'observable': list(initial), 'attempted_operators': denominator,
                'changing_operators_transport_reaction': changing,
                'loss_operators_transport_reaction': losses,
                'gain_operators_transport_reaction': gains,
                'rank_counts': sorted(rank_counts.items()),
                'outcomes_transport_reaction': [sorted(c.items()) for c in counters],
                'mode_changing_operators': witnesses,
                'saturated_vacancy_checks_operator_charge_count_bound': saturated_checks,
                'unique_changed_raw_targets': len(raw_targets),
                'all_unique_targets_compiled_checked': compiled,
                'full_operator_inverse_scan_checked': compiled}


def candidates(source, modes):
    """Every recorded reaction endpoint, plus all specified seed controls."""
    if source['evolution']['side'] != 6:
        raise ValueError('kernel endpoint audit requires the recorded side-six experiment')
    e = FeedbackExperiment(6, source['bank'])
    rows = [{'name': name, 'links': e.geometry.paired_seed(*pair)[0]}
            for name, pair in [('flat', (0, 0)), ('single_reflection', (1, 0)),
                               ('commuting_pair', (1, 1)), ('noncommuting_pair', (1, 2))]]
    if rows[3]['links'] != source['evolution']['inputs'][2]:
        raise ValueError('kernel seed differs from the common source history input')
    for index, reaction in enumerate(modes['dynamics']['reactions']):
        for endpoint in ('before', 'after'):
            rows.append({'name': f"seed_{reaction['seed']}_tick_{reaction['event'][0]}_{endpoint}",
                         'reaction_index': index, 'endpoint': endpoint,
                         'links': reaction[f'{endpoint}_links']})
    return rows


def two_step_escape(probe, initial, initial_kernel):
    """Exhaustive first-loss probability for a one-step-stable control, not a fit."""
    e = probe.e
    rank = initial_kernel['observable'][0]
    denominator = initial_kernel['attempted_operators']
    if dict(initial_kernel['rank_counts']) != {rank: denominator}:
        raise ValueError('two-step control must preserve its count under every first operator')
    supports = len(e.factor.pairs)
    targets = {}
    for rule, table in enumerate(e.tables):
        for patch in range(supports):
            code = e.factor.oracle.code(initial, rule, patch)
            if table[code] != code:
                target = initial[:]
                e.factor.oracle.update(target, rule, patch, table)
                targets.setdefault(tuple(target), []).append(rule*supports+patch)
    first_changing = sum(map(len, targets.values()))
    if not first_changing:
        raise ValueError('two-step control needs a changing operator')
    rows, numerator, witness = [], 0, None
    changing_probability = Fraction(0)
    for target, operators in sorted(targets.items()):
        kernel = probe.kernel(list(target))
        losses = sum(kernel['loss_operators_transport_reaction'])
        second_changing = sum(kernel['changing_operators_transport_reaction'])
        if not second_changing:
            raise ValueError('two-step control has an absorbing intermediate state')
        numerator += len(operators)*losses
        changing_probability += Fraction(len(operators), first_changing)*Fraction(losses, second_changing)
        rows.append({'first_operators': operators, 'intermediate_links': target,
                     'second_step_kernel': kernel})
        if witness is None and losses:
            second = next(encoded for encoded, count in kernel['mode_changing_operators'] if count < rank)
            schedule = [operators[0], second]
            checked = e.compiled_bank([initial], schedule, 1)['runs'][1]
            final = probe.observable(checked['final_links'])
            if final[0] >= rank:
                raise ValueError('compiled two-step escape did not lower the mode count')
            witness = {'schedule': schedule, 'final_links': checked['final_links'],
                       'final_observable': final, 'compiled_exact_inverse': checked['exact_link_inverse']}
    # No-op first steps contribute zero because the initial kernel has no loss.
    probability = Fraction(numerator, denominator**2)
    return {'candidate': 'commuting_pair', 'changed_first_step_targets': rows,
            'losing_ordered_operator_pairs': numerator, 'all_ordered_operator_pairs': denominator**2,
            'first_loss_by_two_attempts': [probability.numerator, probability.denominator],
            'first_loss_by_two_changing_updates': [changing_probability.numerator, changing_probability.denominator],
            'shortest_escape_length': 2 if witness else None, 'escape_witness': witness}


def audit():
    source = json.loads(Path('data/triangle-feedback.json').read_text())
    modes = json.loads(Path('data/triangle-modes.json').read_text())
    if json.loads(json.dumps(reaction_audit())) != modes['dynamics']:
        raise ValueError('recorded endpoint provenance does not reproduce')
    probe = ModeRates(6, source['bank'])
    rows, kernels = candidates(source, modes), {}
    for row in rows:
        key = tuple(row['links'])
        if key not in kernels:
            kernels[key] = probe.kernel(row['links'])
        row['kernel'] = kernels[key]
        print(row['name'], row['kernel']['observable'], 'loss', row['kernel']['loss_operators_transport_reaction'], flush=True)
    # Same observed macrostate but different next-rank marginal is already a
    # sufficient obstruction to strong lumpability of the richer macrostate.
    # Exclude the three separate control seeds. All remaining states are reached
    # from candidate 3 by the independently replayed primitive histories.
    witness = next(([i, j] for i, j in itertools.combinations(range(3, len(rows)), 2)
                    if rows[i]['kernel']['observable'] == rows[j]['kernel']['observable']
                    and rows[i]['kernel']['rank_counts'] != rows[j]['kernel']['rank_counts']), None)
    # A deterministic nonconstant local gauge-frame change checks that a rate
    # difference is not an artifact of representing fiber vertices by labels.
    gauge_rows = []
    e, g = probe.e, probe.e.geometry.group
    for index in (witness or [3]):
        rng = random.Random(73019+index)
        frames = [rng.randrange(g.n) for _ in range(e.geometry.size)]
        original = rows[index]['links']
        transformed = [g.mul[frames[v]][g.mul[x][g.inv[frames[u]]]]
                       for (u, v), x in zip(e.geometry.edges, original)]
        observed = probe.kernel(transformed)
        if observed != rows[index]['kernel']:
            raise ValueError('frozen-state mode kernel changes under local gauge frames')
        gauge_rows.append({'candidate_index': index, 'frames': frames,
                           'transformed_links': transformed, 'exact_same_kernel': True})
    return {'scope': 'Frozen-state one-attempt spectral transition counts for all rooted rules, including no-ops. Supplied scheduler; not physical decay rates or particle identity. Endpoint set is the complete recorded reaction set, not an equilibrium sample.',
            'side': 6, 'observable_fields': ['above_band_rank', 'Q_up', 'Q_down', 'N_reflection', 'N_rotation'],
            'source': ['data/triangle-feedback.json', 'data/triangle-modes.json'],
            'selection': 'Four specified seed controls and both endpoints of every recorded reaction in all four runs.',
            'candidates': rows, 'macrostate_nonclosure_witness': witness, 'gauge_checks': gauge_rows,
            'saturation_obstruction': saturation_audit(source['bank']),
            'one_step_stability_is_not_protection': two_step_escape(probe, rows[2]['links'], rows[2]['kernel'])}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, default=Path('out/triangle-mode-rates.json'))
    args = parser.parse_args()
    data = audit()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(data, separators=(',', ':'))+'\n')
    print('macrostate nonclosure witness', data['macrostate_nonclosure_witness'])


if __name__ == '__main__':
    main()
