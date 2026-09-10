#!/usr/bin/env python3
"""Derived triangle-fiber reactions and a raw-link readout of stored braid memory."""
import argparse
import hashlib
import itertools
import json
import random
import subprocess
from collections import Counter, deque
from pathlib import Path

from fiber_transport import FiberTransport
from triple_channels import positive_certificate
from triple_rule_search import algebra, assess, make_table, minimal_closures
from screen_braid_rules import nullspace
from triangle_spectral_capacity import capacity_certificate


def derive_bank(group):
    states, products, actions, _ = algebra(group)
    closures = minimal_closures(group)
    compiled = json.loads(subprocess.run(['build/wgphysics_triple_census', '--cycle', '3'],
                          text=True, capture_output=True, check=True, timeout=120).stdout)
    if compiled['group_order'] != group.n or [r['id'] for r in compiled['minimal_rules']] != list(range(len(closures))) or [tuple(map(tuple, r['transpositions'])) for r in compiled['minimal_rules']] != closures:
        raise ValueError('C++ breadth-first rule closures disagree with independent finite-action orbits')
    labels = sorted(set(group.sectors)-{group.identity})
    populations = [tuple(sum(group.sectors[a] == s for a in t) for s in labels) for t in states]
    counts, records, chosen, equations = Counter(), [], [], []
    for index, pairs in enumerate(closures):
        table = make_table(pairs, group.n**3)
        audit = assess(group, table)
        rows = [tuple(a-b for a, b in zip(populations[x], populations[y])) for x, y in pairs]
        positive = positive_certificate(rows, len(labels))
        changes_count = any(sum(populations[x]) != sum(populations[y]) for x, y in pairs)
        feedback = audit['same_classes_and_boundary_feedback'] is not None
        counts[positive['exists'], changes_count, feedback] += 1
        record = {'id': index, 'transpositions': pairs, 'positive_certificate': positive,
                  'changes_occupied_count': changes_count, **audit}
        records.append(record)
        # Positivity is a selection hypothesis, not an emergent prediction.
        # The subsequent nullspace derives the common charge of this selected
        # bank; adjacency/equivariance alone do not select it or its clock.
        if positive['exists'] and changes_count:
            chosen.append(index); equations += rows
    basis = nullspace(equations, len(labels))
    if len(basis) != 1 or basis[0] != [1, 2]:
        raise ValueError('triangle bank did not derive the expected single positive charge')
    charges = [0 if s == group.identity else basis[0][labels.index(s)] for s in group.sectors]
    # Independently obtain graph-generated word length, using every transposition
    # across an edge of the triangle. It agrees with the derived charge; it did
    # not select the rules or enter their update tables.
    generators = [i for i, p in enumerate(group.elements) if sum(a != b for a, b in enumerate(p)) == 2]
    distances, queue = {group.identity: 0}, deque([group.identity])
    while queue:
        a = queue.popleft()
        for b in generators:
            c = group.mul[a][b]
            if c not in distances:
                distances[c] = distances[a]+1; queue.append(c)
    if [distances[x] for x in range(group.n)] != charges:
        raise ValueError('derived charge differs from triangle-edge transposition length')
    tables = [make_table(closures[i], group.n**3) for i in chosen]
    sectors = [tuple(group.sectors[a] for a in t) for t in states]
    rate_rows, seen, witness = [], {}, None
    for code, state in enumerate(states):
        rates = tuple(sorted(Counter(sectors[t[code]] for t in tables).items()))
        key = sectors[code], products[code]
        if key in seen and seen[key][1] != rates and witness is None:
            old, previous = seen[key]
            witness = {'input_codes': [old, code], 'input_triples': [states[old], state],
                       'input_classes': sectors[code], 'same_boundary_product': products[code],
                       'outgoing_class_counts': [previous, rates]}
        seen.setdefault(key, (code, rates))
        rate_rows.append({'input': code, 'outgoing_class_counts': rates})
    if witness is None:
        raise ValueError('uniform triangle bank unexpectedly closes on classes and boundary product')
    # The local reaction law is exhaustively checked rather than inferred from
    # one witness. Equal adjacent reflections, but not three equal ones, react.
    reflections = [x for x, q in enumerate(charges) if q == 1]
    local_law = []
    for code, values in enumerate(states):
        count = sum(table[code] != code for table in tables)
        if all(x in reflections for x in values):
            a, b, c = values
            expected = len(tables) if (a == b) != (b == c) else 0
            if count != expected:
                raise ValueError('three-reflection gate differs from adjacent-equality law')
            local_law.append({'input': code, 'changing_rules': count})
        elif sorted(charges[x] for x in values) == [0, 1, 2]:
            if count != 4:
                raise ValueError('reverse reaction does not have four of twelve enabled rules')
        elif count:
            raise ValueError('bank has an unexpected reaction channel')
    return {'group_automorphisms': group.elements, 'sector_labels': labels,
            'minimal_rule_count': len(records), 'minimal_rules': records,
            'classification': [{'positive': a, 'changes_occupied_count': b, 'same_class_boundary_feedback': c, 'rules': n}
                               for (a, b, c), n in sorted(counts.items())],
            'selection': 'All minimal gauge/orientation-equivariant product-preserving involutions with a strictly positive additive class charge and a changing occupied-face count. No feedback or observed trajectory enters selection.',
            'selected_rule_ids': chosen, 'tables': tables, 'charge_basis': basis,
            'element_charges': charges, 'triangle_edge_transposition_lengths': [distances[x] for x in range(group.n)],
            'full_bundle_spectral_capacity': capacity_certificate(group, charges),
            'uniform_bank_class_boundary_rate_witness': witness, 'uniform_bank_rate_rows': rate_rows,
            'three_reflection_gate': local_law}


class FeedbackExperiment(FiberTransport):
    def __init__(self, side, bank):
        super().__init__(side, 3, [(0, 1), (1, 2), (2, 0)])
        if bank['group_automorphisms'] != self.geometry.group.elements and bank['group_automorphisms'] != [list(p) for p in self.geometry.group.elements]:
            raise ValueError('feedback bank and derived fiber group disagree')
        self.tables = self.factor.tables+bank['tables']
        self.charges = bank['element_charges']
        canonical = lambda f: min(tuple(f[k:]+f[:k]) for k in range(3))
        ids = {canonical(list(f[:3])): i for i, f in enumerate(self.geometry.faces)}
        self.fan_faces = [[ids[canonical(list(f[:3]))] for f in spec] for spec in self.factor.oracle.fan.specs]

    def compiled_bank(self, inputs, schedule, stride, capture_links=False):
        side, n = self.geometry.side, self.geometry.group.n
        if not inputs or len(inputs) > 16 or not schedule or len(schedule) > 1000000 or stride <= 0 or len(schedule) % stride:
            raise ValueError('invalid bounded feedback run dimensions')
        supports = len(self.factor.pairs)
        if any(not isinstance(x, int) or not 0 <= x < len(self.tables)*supports for x in schedule):
            raise ValueError('invalid mixed-rule schedule')
        for links in inputs:
            self.geometry.holonomies(links)
        protocol = [side, len(self.tables)-1, len(inputs), len(schedule), stride]
        protocol += [x for table in self.tables for x in table]
        protocol += [x for links in inputs for x in links]+schedule
        command = ['build/wgphysics_mixed_bank_experiments', '--cycle', '3']
        if capture_links:
            command.append('--raw-events')
        result = json.loads(subprocess.run(command,
                            input=' '.join(map(str, protocol))+'\n', text=True, capture_output=True,
                            check=True, timeout=120).stdout)
        if len(result['runs']) != 2*len(inputs):
            raise ValueError('unexpected compiled bank run count')
        for row, (condition, combined) in zip(result['runs'], itertools.product(range(len(inputs)), (False, True))):
            links = inputs[condition][:]
            if capture_links:
                if any(len(event) != 6 for event in row['events']):
                    raise ValueError('compiled raw event snapshot is missing')
                row['raw_event_links'] = [event.pop() for event in row['events']]
            if row['condition'] != condition or row['mode'] != ('combined' if combined else 'transport'):
                raise ValueError('compiled bank condition or mode differs')
            events, history, conversions = [], [self.factor.oracle.histogram(links)], Counter()
            charge = sum(self.charges[x] for x in self.geometry.holonomies(links))
            for tick, encoded in enumerate(schedule, 1):
                rule, patch = divmod(encoded, supports)
                if combined or not rule:
                    code = self.factor.oracle.code(links, rule, patch)
                    target = self.tables[rule][code]
                    if code != target:
                        before, old = self.geometry.holonomies(links), links[:]
                        actual = self.factor.oracle.update(links, rule, patch, self.tables[rule])
                        if actual != (code, target):
                            raise ValueError('independent feedback table application disagrees')
                        after = self.geometry.holonomies(links)
                        affected = set(self.fan_faces[patch] if rule else self.factor.pairs[patch])
                        writes = self.factor.oracle.fan.writes[patch] if rule else {self.factor.oracle.pairs[patch][0][0][0]}
                        if any(a != b for i, (a, b) in enumerate(zip(before, after)) if i not in affected):
                            raise ValueError('feedback changed a spectator holonomy')
                        if any(a != b for i, (a, b) in enumerate(zip(old, links)) if i not in writes):
                            raise ValueError('feedback changed a boundary link')
                        if sum(self.charges[x] for x in after) != charge:
                            raise ValueError('primitive feedback violated the derived additive charge')
                        events.append([tick, rule, patch, code, target])
                        if capture_links and row['raw_event_links'][len(events)-1] != links:
                            raise ValueError('compiled intermediate raw connection differs')
                        if rule:
                            conversions[tuple(sum(self.charges[x] == q for x in after)-sum(self.charges[x] == q for x in before) for q in (1, 2))] += 1
                if tick % stride == 0:
                    history.append(self.factor.oracle.histogram(links))
            if row['events'] != events or row['final_links'] != links or row['histograms'] != history:
                raise ValueError('compiled feedback history differs from independent raw-link replay')
            for _, rule, patch, _, _ in reversed(events):
                self.factor.oracle.update(links, rule, patch, self.tables[rule])
            if links != inputs[condition] or not row['exact_link_inverse'] or not row['reverse_histogram_echo']:
                raise ValueError('independent feedback inverse failed')
            row.update(independent_raw_link_replay=True, boundary_links_unchanged=True,
                       spectator_holonomies_unchanged=True, conserved_charge=charge,
                       conversions=[{'population_change': list(delta), 'count': count} for delta, count in sorted(conversions.items())])
        return result

    def global_class_rates(self, links):
        """All rooted supports and all bank rules, including no-ops, equally weighted."""
        original = tuple(self.geometry.sectors(links))
        rates, changing_by_arity = Counter(), [0, 0]
        for rule, table in enumerate(self.tables):
            for patch in range(len(self.factor.pairs)):
                code = self.factor.oracle.code(links, rule, patch)
                if table[code] == code:
                    rates[original] += 1
                    continue
                moved = links[:]
                self.factor.oracle.update(moved, rule, patch, table)
                classes = tuple(self.geometry.sectors(moved))
                rates[classes] += 1
                changing_by_arity[bool(rule)] += int(classes != original)
        encoded = json.dumps(sorted(rates.items()), separators=(',', ':'))
        return rates, {'attempted_operators': len(self.tables)*len(self.factor.pairs),
                       'class_changing_pair_operators': changing_by_arity[0],
                       'class_changing_triple_operators': changing_by_arity[1],
                       'full_class_rate_sha256': hashlib.sha256(encoded.encode()).hexdigest()}

    def rates_from_compiled_scan(self, initial, schedule, run):
        """Read class outcomes from C++ event targets, independently of Python updates."""
        initial_classes = self.geometry.sectors(initial)
        events = {row[0]: row for row in run['events']}
        rates = Counter()
        for index in range(0, len(schedule), 2):
            if schedule[index] != schedule[index+1]:
                raise ValueError('rate scan must undo each operator before the next')
            row, inverse = events.get(index+1), events.get(index+2)
            classes = initial_classes[:]
            if row is not None:
                _, rule, patch, code, target = row
                if inverse != [index+2, rule, patch, target, code]:
                    raise ValueError('compiled rate scan did not invert a changing event')
                faces = self.fan_faces[patch][::-1] if rule else self.factor.pairs[patch]
                values = []
                for _ in faces:
                    target, value = divmod(target, self.geometry.group.n)
                    values.insert(0, value)
                for face, value in zip(faces, values):
                    classes[face] = self.geometry.group.sectors[value]
            elif inverse is not None:
                raise ValueError('compiled inverse acted after an identity event')
            rates[tuple(classes)] += 1
        if run['final_links'] != initial:
            raise ValueError('compiled rate scan did not restore its input')
        return rates

    def readout(self):
        source = json.loads(Path('data/fiber-transport-order.json').read_text())
        if source['side'] != self.geometry.side:
            raise ValueError('readout needs the saved braid experiment mesh')
        source = source['fibers'][0]
        fan = self.factor.oracle.fan
        root = 6*self.geometry.side+6
        patch = next(i for i, spec in enumerate(fan.specs) if spec[0][0] == root)
        targets = self.fan_faces[patch]
        prepared, states = [], []
        for initial in source['controlled_circuit_orbit']['raw_representatives']:
            links, events = initial[:], []
            for origin, target in zip(source['preparation']['positions'][:3], targets):
                self.move(links, origin, target, events)
            schedule = [event[2] for event in events]
            checked = self.compiled(initial, schedule)
            if checked['final_links'] != links:
                raise ValueError('braid-memory encounter preparation replay disagrees')
            prepared.append({'schedule': schedule, 'checked': checked})
            states.append(links)
        if any(p['schedule'] != prepared[0]['schedule'] for p in prepared):
            raise ValueError('encounter used state-dependent routing')
        if len({tuple(self.geometry.sectors(s)) for s in states}) != 1:
            raise ValueError('readout does not start with identical full face-class fields')
        distributions, rates = zip(*(self.global_class_rates(links) for links in states))
        triples = [tuple(fan.transport(s, path) for path in fan.patches[patch][::-1]) for s in states]
        g = self.geometry.group
        products = [g.mul[g.mul[a][b]][c] for a, b, c in triples]
        differing = next((a, b, output) for a, b in itertools.combinations(range(len(states)), 2)
                         if products[a] == products[b]
                         for output in sorted(set(distributions[a]) | set(distributions[b]))
                         if distributions[a][output] != distributions[b][output])
        a, b, output = differing
        supports = len(self.factor.pairs)
        scan = [rule*supports+patch for rule in range(len(self.tables)) for patch in range(supports) for _ in range(2)]
        scan_checked = self.compiled_bank(states, scan, len(scan))
        for i, state in enumerate(states):
            if self.rates_from_compiled_scan(state, scan, scan_checked['runs'][2*i+1]) != distributions[i]:
                raise ValueError('full compiled global generator rates disagree with independent Python enumeration')
        # Every rule at the same encounter is independently applied from each
        # original state. Append its inverse to return before the next readout.
        schedule = [r*supports+patch for r in range(1, len(self.tables)) for _ in range(2)]
        checked = self.compiled_bank(states, schedule, 1)
        return {'source': 'data/fiber-transport-order.json: first fiber controlled gauge component',
                'encounter_patch': patch, 'encounter_faces': targets,
                'preparations': prepared,
                'encounter_triples': triples, 'encounter_products': products,
                'global_class_rates': rates,
                'global_rate_scan_schedule': scan, 'checked_global_rate_scan': scan_checked,
                'global_rate_witness': {'states': [a, b], 'output_face_classes': output,
                                        'same_encounter_boundary_product': products[a],
                                        'counts': [distributions[a][output], distributions[b][output]]},
                'readout_schedule': schedule, 'checked_readout_and_inverse': checked}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--attempts', type=int, default=100000)
    parser.add_argument('--output', type=Path, default=Path('out/triangle-feedback.json'))
    args = parser.parse_args()
    if not 1 <= args.attempts <= 1000000:
        raise ValueError('attempts must be in 1..1000000')
    bank = derive_bank(FiberTransport(3, 3, [(0, 1), (1, 2), (2, 0)]).geometry.group)
    readout = FeedbackExperiment(12, bank).readout()
    e = FeedbackExperiment(6, bank)
    inputs = [[0]*len(e.geometry.edges), e.geometry.paired_seed(1, 1)[0], e.geometry.paired_seed(1, 2)[0]]
    runs = []
    for seed in range(4):
        rng = random.Random(seed)
        schedule = [rng.randrange(len(e.tables)*len(e.factor.pairs)) for _ in range(args.attempts)]
        checked = e.compiled_bank(inputs, schedule, max(1, args.attempts//100) if args.attempts % 100 == 0 else args.attempts)
        runs.append({'seed': seed, 'schedule': schedule, 'checked': checked})
        print('seed', seed, 'combined noncommuting conversions', checked['runs'][-1]['conversions'], flush=True)
    result = {'scope': 'Classical derived-fiber reaction feedback and controlled braid-memory readout. Positive combinatorial charge is not physical energy; no quantum amplitudes, emergent geometry, or binding asserted.',
              'bank': bank, 'readout': readout,
              'evolution': {'side': 6, 'conditions': ['flat', 'commuting_reflections', 'noncommuting_reflections'],
                            'inputs': inputs, 'attempts': args.attempts, 'runs': runs}}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, separators=(',', ':'))+'\n')
    print('minimal rules', bank['minimal_rule_count'], 'selected', len(bank['tables']),
          'global class-changing reaction counts', [r['class_changing_triple_operators'] for r in readout['global_class_rates']])


if __name__ == '__main__':
    main()
