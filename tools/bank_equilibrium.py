#!/usr/bin/env python3
"""Reachable canonical states of a uniformly sampled nonabelian rule bank.

The class graph is an exact stochastic factor, not a deterministic per-rule
replacement for link dynamics. Positive paths are lifted back to actual links.
"""
import argparse
import hashlib
import itertools
import json
import math
import random
import statistics
import struct
import subprocess
from collections import Counter
from fractions import Fraction
from pathlib import Path

from face_energy_obstruction import derive_fiber_group
from run_channel_banks import MixedOracle, seeds, verify_run
from run_shared_edge import mesh_geometry
from triple_rule_search import assess, make_table


def encode(state):
    result = 0
    for value in state:
        result = 8*result+value
    return result


def decode(code, arity):
    return tuple((code//8**i) % 8 for i in reversed(range(arity)))


def fraction(value):
    return [value.numerator, value.denominator]


class BankFactor:
    def __init__(self, side, channels, family_id=0):
        self.group = derive_fiber_group()
        self.oracle = MixedOracle(side, self.group)
        self.family = channels['families'][family_id]
        self.tables = [channels['shared_edge_vacancy_transport']]+[
            make_table(channels['rules'][i]['transpositions']) for i in self.family['rule_ids']]
        if len(self.tables) != 49 or len({tuple(t) for t in self.tables[1:]}) != 48:
            raise ValueError('expected a complete distinct 48-rule family')
        for table in self.tables[1:]:
            assess(self.group, table)
        for a, b in itertools.product(range(8), repeat=2):
            expected = (b, a) if (a == 0) != (b == 0) else (a, b)
            if self.tables[0][encode((a, b))] != encode(expected):
                raise ValueError('pair table is not vacancy transport')
        species = self.family['reaction_species']
        self.h = species['heavy_class']
        self.a, self.b = species['light_classes']
        if {self.h, self.a, self.b} != {1, 2, 3}:
            raise ValueError('invalid class assignment')

        # Derive the averaged kernel from all raw representatives, not just the
        # artifact's coarse rows. A positive coarse edge then lifts from EVERY
        # link state with that observed input, even after earlier path choices.
        self.rates = {}
        for code in range(512):
            before = tuple(self.group.sectors[x] for x in decode(code, 3))
            counts = Counter(tuple(self.group.sectors[x] for x in decode(t[code], 3))
                             for t in self.tables[1:])
            if before in self.rates and self.rates[before] != counts:
                raise ValueError('uniform bank fails strong class lumpability')
            self.rates[before] = counts
        saved = {tuple(r['input']): Counter({tuple(s): n for s, n in r['outputs']})
                 for r in self.family['uniform_bank_class_rates']}
        if self.rates != saved:
            raise ValueError('saved class rates differ from raw rule bank')
        for before, outputs in self.rates.items():
            for after, count in outputs.items():
                if self.signature(before) != self.signature(after):
                    raise ValueError('local bank charge drift')
                if self.weight(before)*count != self.weight(after)*self.rates[after][before]:
                    raise ValueError('local degeneracy balance fails')
        self.targets = {encode(before): tuple(encode(after) for after in sorted(outputs) if after != before)
                        for before, outputs in self.rates.items()}
        canonical = lambda f: min(f[k:]+f[:k] for k in range(len(f)))
        ids = {canonical(tuple(f[:-1])): i for i, f in enumerate(self.oracle.fan.faces)}
        self.fans = [tuple(ids[canonical(tuple(f[:-1]))] for f in spec[::-1])
                     for spec in self.oracle.fan.specs]
        _, _, pair_specs = mesh_geometry(side)
        self.pairs = [tuple(spec[2:]) for spec in pair_specs]
        # Both roots of a shared edge give the same CLASS swap. Deduplicate only
        # for reachability, retaining an actual rooted patch for path lifting.
        distinct = {}
        for p, faces in enumerate(self.pairs):
            distinct.setdefault(tuple(sorted(faces)), p)
        self.supports = [(0, p, self.pairs[p]) for p in distinct.values()]
        self.supports += [(1, p, faces) for p, faces in enumerate(self.fans)]
        self.side = side

    def project(self, links):
        return bytes(self.oracle.fan.sectors(links))

    def signature(self, state):
        return (state.count(self.a)+state.count(self.h),
                state.count(self.b)+state.count(self.h), state.count(5))

    @staticmethod
    def weight(state):
        return 2**sum(x in (1, 2, 3) for x in state)

    def successors(self, state):
        for support, (kind, patch, faces) in enumerate(self.supports):
            if kind:
                code = encode(state[f] for f in faces)
                targets = self.targets[code]
            else:
                a, b = (state[f] for f in faces)
                targets = (encode((b, a)),) if (a == 0) != (b == 0) else ()
            for target in targets:
                out = bytearray(state)
                for f, value in zip(faces, decode(target, len(faces))):
                    out[f] = value
                yield bytes(out), support, target

    def lift_path(self, initial, path):
        links = initial[:]
        state = self.project(links)
        schedule = []
        count = len(self.fans)
        for support, target in path:
            kind, patch, faces = self.supports[support]
            code = self.oracle.code(links, kind, patch)
            choices = range(1, len(self.tables)) if kind else (0,)
            rule = next((r for r in choices if encode(self.group.sectors[x] for x in
                        decode(self.tables[r][code], len(faces))) == target), None)
            if rule is None:
                raise ValueError('positive coarse path did not lift from raw representative')
            self.oracle.update(links, rule, patch, self.tables[rule])
            expected = bytearray(state)
            for f, value in zip(faces, decode(target, len(faces))):
                expected[f] = value
            state = bytes(expected)
            if self.project(links) != state:
                raise ValueError('raw lift changed a spectator or missed a class target')
            schedule.append(rule*count+patch)
        final = links[:]
        for encoded in reversed(schedule):
            rule, patch = divmod(encoded, count)
            self.oracle.update(links, rule, patch, self.tables[rule])
        if links != initial:
            raise ValueError('lifted path failed exact raw inverse')
        return schedule, final


def canonical_counts(faces, charges):
    """All formal class configurations at the given three additive charges."""
    qa, qb, z = charges
    if faces < 0 or any(q < 0 for q in charges):
        raise ValueError('negative dimensions or charges')
    result = []
    for h in range(min(qa, qb)+1):
        a, b = qa-h, qb-h
        empty = faces-a-b-h-z
        if empty < 0:
            continue
        counts = [empty, a, b, h, z]
        configurations = math.factorial(faces)//math.prod(math.factorial(n) for n in counts)
        weight = 2**(a+b+h)
        result.append({'heavy': h, 'populations_empty_a_b_h_z': counts,
                       'configurations': configurations, 'weight_per_configuration': weight,
                       'stationary_weight': configurations*weight})
    total = sum(r['stationary_weight'] for r in result)
    for row in result:
        row['probability'] = fraction(Fraction(row['stationary_weight'], total))
    return result


def canonical_rates(factor, rows):
    """Stationary conditional reaction rates, NOT a closed population dynamics."""
    rates = []
    overall = Fraction(0)
    for row in rows:
        counts = dict(zip((0, factor.a, factor.b, factor.h, 5), row['populations_empty_a_b_h_z']))
        faces = sum(counts.values())
        if faces < 3:
            raise ValueError('three-face rates require at least three faces')
        changes = Counter()
        for before, outputs in factor.rates.items():
            remaining = counts.copy()
            multiplicity = 1
            for value in before:
                multiplicity *= remaining[value]
                remaining[value] -= 1
            if not multiplicity:
                continue
            # All rooted fans contain three distinct faces. Canonical conditional
            # exchangeability supplies this hypergeometric marginal, not an iid
            # approximation. Each conversion TABLE has attempt weight 1/49.
            probability = Fraction(multiplicity, faces*(faces-1)*(faces-2)*49)
            for after, count in outputs.items():
                delta = after.count(factor.h)-before.count(factor.h)
                if delta:
                    changes[delta] += probability*count
        rates.append({'heavy': row['heavy'], 'changes': [[d, fraction(p)] for d, p in sorted(changes.items())]})
        overall += Fraction(*row['probability'])*sum(changes.values())
    lookup = {r['heavy']: {d: Fraction(*p) for d, p in r['changes']} for r in rates}
    probabilities = {r['heavy']: Fraction(*r['probability']) for r in rows}
    for h, changes in lookup.items():
        for d, p in changes.items():
            if probabilities[h]*p != probabilities[h+d]*lookup[h+d][-d]:
                raise ValueError('canonical population flux balance fails')
    return {'conditional_rates': rates, 'conversions_per_attempt': fraction(overall),
            'interpretation': 'stationary conditional expectations; populations alone are not asserted Markov'}


def spatial_reference(factor, rows):
    """Class a-b correlations follow geometry alone in this complete component."""
    faces = len(factor.oracle.fan.faces)
    adjacency = [set() for _ in range(faces)]
    for a, b in factor.pairs:
        adjacency[a].add(b)
        adjacency[b].add(a)
    distances = []
    for root in range(faces):
        row, queue = [-1]*faces, [root]
        row[root] = 0
        for current in queue:
            for other in sorted(adjacency[current]):
                if row[other] < 0:
                    row[other] = row[current]+1
                    queue.append(other)
        if -1 in row:
            raise ValueError('spatial reference requires a connected dual graph')
        distances.append(row)
    shells = Counter(distances[a][b] for a in range(faces) for b in range(faces) if a != b)
    pairs = sum(Fraction(*r['probability'])*r['populations_empty_a_b_h_z'][1]*r['populations_empty_a_b_h_z'][2]
                for r in rows)
    return {'distance': 'shortest-path distance on the supplied face-dual graph, not causal or branchial distance',
            'distances': distances, 'ordered_face_pairs_by_distance': sorted(shells.items()),
            'expected_total_a_b_pairs': fraction(pairs),
            'expected_a_b_pairs_by_distance': [[d, fraction(pairs*n/(faces*(faces-1)))] for d, n in sorted(shells.items())],
            'expected_a_b_per_ordered_face_pair': fraction(pairs/(faces*(faces-1))),
            'scope': 'stationary face-class pairs only; does not exclude transient pairing or correlations of larger loop observables'}


def reachable(factor, initial, budget, progress=False):
    if budget < 1:
        raise ValueError('state budget must be positive')
    start = factor.project(initial)
    signature = factor.signature(start)
    states, parents, lookup = [start], [None], {start: 0}
    witnesses = {start.count(factor.h): 0}
    exhausted = True
    for index, state in enumerate(states):
        for after, support, target in factor.successors(state):
            if after in lookup:
                continue
            if len(states) == budget:
                exhausted = False
                break
            if factor.signature(after) != signature:
                raise ValueError('reachable state charge drift')
            lookup[after] = len(states)
            witnesses.setdefault(after.count(factor.h), len(states))
            states.append(after)
            parents.append((index, support, target))
        if not exhausted:
            break
        if progress and index and index % 50000 == 0:
            print('Expanded', index, 'discovered', len(states), flush=True)
    counts = Counter(s.count(factor.h) for s in states)
    formal = canonical_counts(len(start), signature)
    expected = {r['heavy']: r['configurations'] for r in formal}
    full = exhausted and counts == expected
    component_weights = {h: n*2**(signature[0]+signature[1]-h) for h, n in counts.items()}
    component_total = sum(component_weights.values())
    # Full equality proves realizability by induction from the raw seed, using
    # strong lumpability on every BFS edge. No independent-face ansatz is used.
    records = []
    for h, index in sorted(witnesses.items()):
        path, cursor = [], index
        while parents[cursor] is not None:
            previous, support, target = parents[cursor]
            path.append((support, target))
            cursor = previous
        path.reverse()
        schedule, links = factor.lift_path(initial, path)
        if factor.project(links) != states[index]:
            raise ValueError('witness lift reached the wrong state')
        records.append({'heavy': h, 'changing_events': len(path), 'path': path,
                        'schedule': schedule, 'final_links': links,
                        'final_classes': list(states[index]), 'independent_link_inverse': True})
    return {'side': factor.side, 'family': factor.family['id'], 'charges': signature,
            'initial_links': initial, 'initial_classes': list(start), 'state_budget': budget,
            'exhausted': exhausted, 'reachable_states': len(states),
            'reachable_counts_by_heavy': sorted(counts.items()),
            'component_stationary_probabilities': [[h, fraction(Fraction(w, component_total))]
                                                   for h, w in sorted(component_weights.items())] if exhausted else None,
            'formal_canonical_counts': formal, 'all_formal_states_reachable': full,
            'stationary_prediction_validated_for_this_seed_component': full,
            'minimum_idle_triple_rules': min(outputs[before] for before, outputs in factor.rates.items()),
            'canonical_reaction_rates': canonical_rates(factor, formal),
            'formal_canonical_spatial_reference': spatial_reference(factor, formal),
            'witnesses': records}


def compiled_run(factor, initial, schedule, stride):
    inputs = [factor.side, len(factor.tables)-1, 1, len(schedule), stride]
    inputs += [x for t in factor.tables for x in t]+initial+schedule
    process = subprocess.run(['build/wgphysics_mixed_bank_experiments'],
                             input=' '.join(map(str, inputs))+'\n', text=True,
                             capture_output=True, check=True, timeout=1800)
    runs = json.loads(process.stdout)['runs']
    if [(r['condition'], r['mode']) for r in runs] != [(0, 'transport'), (0, 'combined')]:
        raise ValueError('compiled runner omitted controls')
    for run in runs:
        verify_run(factor.oracle, initial, schedule, factor.tables, factor.family, run, stride)
    return runs


def residence(factor, initial, events, attempts, burn, block):
    """Exact post-attempt occupation times, with self loops and no subsampling."""
    if not 0 <= burn < attempts or block <= 0 or (attempts-burn) % block:
        raise ValueError('invalid residence clock')
    max_h = min(factor.signature(factor.project(initial))[:2])
    blocks = [[0]*(max_h+1) for _ in range((attempts-burn)//block)]
    heavy = factor.project(initial).count(factor.h)
    previous, conversions = 0, 0
    def accumulate(start, stop, h):
        start, stop = max(start, burn+1), min(stop, attempts+1)
        while start < stop:
            index = (start-burn-1)//block
            end = min(stop, burn+1+(index+1)*block)
            blocks[index][h] += end-start
            start = end
    for tick, rule, patch, code, target in events:
        if not previous < tick <= attempts:
            raise ValueError('event clock is not strictly increasing')
        accumulate(previous, tick, heavy)
        arity = 3 if rule else 2
        heavy += sum(factor.group.sectors[x] == factor.h for x in decode(target, arity))
        heavy -= sum(factor.group.sectors[x] == factor.h for x in decode(code, arity))
        if not 0 <= heavy <= max_h:
            raise ValueError('event population leaves charge sector')
        conversions += bool(rule and tick > burn)
        previous = tick
    accumulate(previous, attempts+1, heavy)
    if any(sum(row) != block for row in blocks):
        raise ValueError('residence clock lost or double-counted attempts')
    totals = [sum(row[h] for row in blocks) for h in range(max_h+1)]
    return {'occupation_attempts': totals, 'block_occupation_attempts': blocks,
            'post_burn_conversions': conversions,
            'occupation_probabilities': [x/(attempts-burn) for x in totals],
            'conversion_rate_per_attempt': conversions/(attempts-burn)}


def spatial_residence(factor, initial, run, attempts, burn, distances):
    """Measure pair-distance occupation from already raw-link-verified events."""
    if not 0 <= burn < attempts:
        raise ValueError('invalid spatial residence clock')
    state = factor.project(initial)
    diameter = max(max(row) for row in distances)
    totals, previous = [0]*(diameter+1), 0
    def accumulate(start, stop):
        duration = max(0, min(stop, attempts+1)-max(start, burn+1))
        if duration:
            for a, value in enumerate(state):
                if value == factor.a:
                    for b, other in enumerate(state):
                        if other == factor.b:
                            totals[distances[a][b]] += duration
    for tick, rule, patch, code, target in run['events']:
        if not previous < tick <= attempts:
            raise ValueError('spatial event clock is not strictly increasing')
        accumulate(previous, tick)
        faces = factor.fans[patch] if rule else factor.pairs[patch]
        before = tuple(factor.group.sectors[x] for x in decode(code, len(faces)))
        if before != tuple(state[f] for f in faces):
            raise ValueError('spatial observer differs from raw event input')
        after = bytearray(state)
        for f, x in zip(faces, decode(target, len(faces))):
            after[f] = factor.group.sectors[x]
        state = bytes(after)
        previous = tick
    accumulate(previous, attempts+1)
    if state != factor.project(run['final_links']):
        raise ValueError('spatial observer differs from final raw links')
    return {'a_b_pair_attempts_by_distance': totals,
            'mean_a_b_pairs_by_distance': [n/(attempts-burn) for n in totals]}


def relaxation(factor, component, attempts, burn, block, trials, seed):
    if not component['all_formal_states_reachable']:
        raise ValueError('canonical trajectory comparison requires exhaustive component validation')
    if not 0 < attempts <= 1000000 or not 0 <= burn < attempts or block <= 0 or (attempts-burn) % block or not 1 <= trials <= 32:
        raise ValueError('invalid bounded trajectory dimensions')
    runs = []
    for witness in component['witnesses']:
        for trial in range(trials):
            run_seed = seed+100003*trial+10000019*witness['heavy']
            rng = random.Random(run_seed)
            schedule = [rng.randrange(49*len(factor.fans)) for _ in range(attempts)]
            initial = witness['final_links']
            control, run = compiled_run(factor, initial, schedule, max(1, math.gcd(attempts, block)))
            measurement = residence(factor, initial, run['events'], attempts, burn, block)
            measurement.update(spatial_residence(factor, initial, run, attempts, burn,
                                                component['formal_canonical_spatial_reference']['distances']))
            null = residence(factor, initial, control['events'], attempts, burn, block)
            if null['post_burn_conversions'] or null['occupation_attempts'][witness['heavy']] != attempts-burn:
                raise ValueError('transport control changed class populations')
            run.update(measurement, initial_heavy=witness['heavy'], trial=trial, seed=run_seed,
                       initial_links=initial, schedule_sha256=hashlib.sha256(
                           b''.join(struct.pack('<I', s) for s in schedule)).hexdigest(),
                       transport_control={'final_links': control['final_links'],
                                          'independent_link_replay': control['independent_link_replay'],
                                          'exact_link_inverse': control['exact_link_inverse'], **null})
            runs.append(run)
            print('Initial heavy', witness['heavy'], 'trial', trial, 'occupations',
                  measurement['occupation_probabilities'], 'rate', measurement['conversion_rate_per_attempt'], flush=True)
    def summary(selected):
        def estimate(values):
            return {'mean': statistics.mean(values), 'independent_run_standard_error':
                    statistics.stdev(values)/math.sqrt(len(values)) if len(values) > 1 else None}
        changes = []
        for run in selected:
            blocks = run['block_occupation_attempts']
            mid = len(blocks)//2
            if not mid:
                break
            earlier = [sum(row[h] for row in blocks[:mid])/(mid*block) for h in range(len(blocks[0]))]
            later = [sum(row[h] for row in blocks[mid:])/((len(blocks)-mid)*block) for h in range(len(blocks[0]))]
            changes.append([a-b for a, b in zip(later, earlier)])
        return {'runs': len(selected), 'occupations': [estimate([r['occupation_probabilities'][h] for r in selected])
                                                     for h in range(len(runs[0]['occupation_attempts']))],
                'late_minus_early_occupation': [estimate([row[h] for row in changes])
                                               for h in range(len(changes[0]))] if changes else None,
                'a_b_pairs_by_distance': [estimate([r['mean_a_b_pairs_by_distance'][d] for r in selected])
                                         for d in range(len(selected[0]['mean_a_b_pairs_by_distance']))],
                'conversion_rate': estimate([r['conversion_rate_per_attempt'] for r in selected])}
    return {'attempts_per_run': attempts, 'burn_attempts': burn, 'block_attempts': block,
            'trials_per_initial_population': trials, 'seed': seed,
            'schedule_convention': 'Python random.Random(run.seed).randrange(49 * rooted_fan_count), iid attempts; sha256 uses little-endian uint32',
            'measurement': 'all post-attempt states after burn, including idle attempts; uncertainty across independently seeded run means, not iid ticks',
            'pooled': summary(runs),
            'by_initial_heavy': [{'initial_heavy': h, **summary([r for r in runs if r['initial_heavy'] == h])}
                                 for h in sorted({r['initial_heavy'] for r in runs})],
            'runs': runs}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--side', type=int, default=3)
    parser.add_argument('--family', type=int, choices=range(3), default=0)
    parser.add_argument('--state-budget', type=int, default=500000)
    parser.add_argument('--component', type=Path, help='replay a saved exhaustive result by re-running its search')
    parser.add_argument('--trials', type=int, default=0, help='raw trajectories per initial heavy population; zero skips')
    parser.add_argument('--attempts', type=int, default=1000000)
    parser.add_argument('--burn', type=int, default=200000)
    parser.add_argument('--block', type=int, default=50000)
    parser.add_argument('--seed', type=int, default=8320917)
    parser.add_argument('--output', type=Path, default=Path('out/bank-equilibrium.json'))
    args = parser.parse_args()
    if not 3 <= args.side <= 24 or not 1 <= args.state_budget <= 1000000:
        parser.error('invalid bounded reachability dimensions')
    channels = json.loads(Path('data/d4-triple-channels.json').read_text())
    factor = BankFactor(args.side, channels, args.family)
    # Pick the actual two-link seed h+a, independent of family class labels.
    initial = next(x for x in seeds(factor.oracle, args.side)
                   if factor.signature(factor.project(x)) == (4, 2, 0))
    result = reachable(factor, initial, args.state_budget, progress=True)
    for witness in result['witnesses']:
        if witness['schedule']:
            run = compiled_run(factor, initial, witness['schedule'], 1)[1]
            if run['final_links'] != witness['final_links']:
                raise ValueError('compiled path lift disagrees with independent links')
            witness['compiled_link_inverse'] = run['exact_link_inverse']
    if args.component:
        expected = json.loads(args.component.read_text())
        if any(json.loads(json.dumps(value)) != expected.get(key) for key, value in result.items()):
            raise ValueError('saved component differs from repeated exact search')
    if args.trials:
        result['relaxation'] = relaxation(factor, result, args.attempts, args.burn, args.block, args.trials, args.seed)
    result['scope'] = ('Finite classical gauge-link dynamics on a supplied torus under iid uniform '
                       'rule/placement sampling; coarse stationary weights are not full gauge-orbit counts.')
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, separators=(',', ':'))+'\n')
    print('Reachability exhausted:', result['exhausted'], 'all formal states reachable:',
          result['all_formal_states_reachable'], 'states:', result['reachable_states'])


if __name__ == '__main__':
    main()
