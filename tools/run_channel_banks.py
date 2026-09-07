#!/usr/bin/env python3
"""Run every compatible nonabelian channel bank on shared raw links, with controls."""
import argparse
import json
import random
import subprocess
from pathlib import Path

from face_energy_obstruction import derive_fiber_group
from run_three_face import LinkOracle
from run_shared_edge import mesh_geometry
from triple_rule_search import make_table
from triple_channels import analyze as analyze_channels, subgroup

CONDITIONS = ('flat', 'r', 's', 'rotation', 'r_s', 'r_rotation', 's_rotation')


class MixedOracle:
    def __init__(self, side, group):
        self.fan = LinkOracle(side, group)
        self.group = group
        _, _, specs = mesh_geometry(side)
        ids = {edge: i for i, edge in enumerate(self.fan.edges)}
        def path(loop): return [(ids[tuple(sorted((u, v)))], u > v) for u, v in zip(loop, loop[1:])]
        self.pairs = [(path(a), path(b)) for a, b, _, _ in specs]
        if len(self.pairs) != len(self.fan.patches):
            raise ValueError('rooted support counts disagree')

    def code(self, links, rule, patch):
        paths = self.fan.patches[patch][::-1] if rule else self.pairs[patch]
        result = 0
        for path in paths:
            result = 8*result+self.fan.transport(links, path)
        return result

    def update(self, links, rule, patch, table):
        if rule:
            return self.fan.update(links, patch, table)
        g = self.group
        a, b = self.pairs[patch]
        code = self.code(links, 0, patch)
        target = table[code]
        x, y = divmod(target, 8)
        q = self.fan.transport(links, a[1:])
        new = g.mul[g.inv[q]][x]
        edge, reverse = a[0]
        links[edge] = g.inv[new] if reverse else new
        if (self.fan.transport(links, a), self.fan.transport(links, b)) != (x, y):
            raise ValueError('independent shared-edge reconstruction missed targets')
        return code, target

    def histogram(self, links):
        sectors = self.fan.sectors(links)
        return [sectors.count(a) for a in range(8)]


def seeds(oracle, side):
    first = (side//2)*side+side//2
    other = (side//2)*side+(side//2+max(1, side//3)) % side
    second_end = (other//side)*side+(other+1) % side
    ids = {edge: i for i, edge in enumerate(oracle.fan.edges)}
    result = []
    values = ((), (1,), (2,), (3,), (1, 2), (1, 3), (2, 3))
    for generators in values:
        links = [0]*len(ids)
        for (u, v), a in zip(((first, first+1), (other, second_end)), generators):
            links[ids[tuple(sorted((u, v)))]] = a if u < v else oracle.group.inv[a]
        result.append(links)
    return result


def verify_run(oracle, initial, schedule, tables, family, run, stride):
    g = oracle.group
    labels = sorted(set(g.sectors)-{0})
    charges = [[0 if s == 0 else q[labels.index(s)] for s in g.sectors] for q in family['charge_basis']]
    weights = [0 if s == 0 else family['positive_certificate']['weights'][labels.index(s)] for s in g.sectors]
    total = lambda hist: [sum(q[a]*hist[a] for a in range(8)) for q in charges]
    hist = oracle.histogram(initial)
    conserved = total(hist)
    positive_weight = sum(weights[a]*hist[a] for a in range(8))
    if run['histograms'][0] != hist or len(run['histograms']) != len(schedule)//stride+1:
        raise ValueError('bank initial histogram or sampling clock disagrees')
    links, expected_events = initial[:], []
    count = len(oracle.fan.patches)
    max_nonflat = len(oracle.fan.faces)-hist[0]
    for tick, encoded in enumerate(schedule, 1):
        rule, patch = divmod(encoded, count)
        if run['mode'] == 'combined' or not rule:
            code = oracle.code(links, rule, patch)
            target = tables[rule][code]
            # Identity targets are proven to fix all written links in both lifts.
            if target != code:
                before, after = oracle.update(links, rule, patch, tables[rule])
                if (before, after) != (code, target): raise ValueError('bank local target disagreement')
                expected_events.append([tick, rule, patch, code, target])
                decode = lambda x: (x//64, (x//8) % 8, x % 8) if rule else divmod(x, 8)
                old, new = decode(code), decode(target)
                if rule and (len(subgroup(g, old)) != 8 or len(subgroup(g, new)) != 8):
                    raise ValueError('claimed full-group conversion lies in a proper subgroup')
                for a in old: hist[g.sectors[a]] -= 1
                for a in new: hist[g.sectors[a]] += 1
                if total(hist) != conserved: raise ValueError('derived bank charges drifted at an event')
                max_nonflat = max(max_nonflat, len(oracle.fan.faces)-hist[0])
        if tick % stride == 0:
            if hist != oracle.histogram(links) or hist != run['histograms'][tick//stride]:
                raise ValueError('bank histogram differs from link replay')
            if total(hist) != conserved: raise ValueError('derived bank charges drifted')
            max_nonflat = max(max_nonflat, len(oracle.fan.faces)-hist[0])
    if run['events'] != expected_events or run['final_links'] != links:
        raise ValueError('bank event log or final raw links differs from independent replay')
    if max_nonflat > positive_weight: raise ValueError('positive charge failed to bound nonflat support')
    for _, rule, patch, _, _ in reversed(expected_events):
        oracle.update(links, rule, patch, tables[rule])
    if links != initial or not run['exact_link_inverse'] or not run['reverse_histogram_echo']:
        raise ValueError('bank inverse echo failed')
    run.update(derived_charges=conserved, positive_weight=positive_weight,
               max_nonflat=max_nonflat, independent_link_replay=True,
               conversion_events=sum(bool(e[1]) for e in expected_events))


def experiment(channels, side, attempts, stride, trials, seed):
    g = derive_fiber_group()
    oracle = MixedOracle(side, g)
    initial = seeds(oracle, side)
    count = len(oracle.fan.patches)
    families = channels['families']
    if len({len(f['rule_ids']) for f in families}) != 1:
        raise ValueError('matched family experiment requires equal generator counts')
    schedules = []
    for trial in range(trials):
        rng = random.Random(seed+100003*trial)
        schedules.append([rng.randrange(count*(len(families[0]['rule_ids'])+1)) for _ in range(attempts)])
    runs = []
    for family in families:
        tables = [channels['shared_edge_vacancy_transport']]+[make_table(channels['rules'][i]['transpositions']) for i in family['rule_ids']]
        for trial, schedule in enumerate(schedules):
            inputs = [side, len(tables)-1, len(initial), attempts, stride]
            inputs += [x for table in tables for x in table]
            inputs += [x for links in initial for x in links]+schedule
            process = subprocess.run(['build/wgphysics_mixed_bank_experiments'], input=' '.join(map(str, inputs))+'\n',
                                     text=True, capture_output=True, check=True, timeout=1800)
            raw = json.loads(process.stdout)['runs']
            if [(r['condition'], r['mode']) for r in raw] != [(c, mode) for c in range(len(initial)) for mode in ('transport', 'combined')]:
                raise ValueError('bank runner drops, duplicates, or reorders controls')
            for run in raw:
                verify_run(oracle, initial[run['condition']], schedule, tables, family, run, stride)
                run.update(family=family['id'], trial=trial, condition=CONDITIONS[run['condition']])
            for control, together in zip(raw[::2], raw[1::2]):
                together['control_face_class_difference'] = sum(a != b for a, b in zip(oracle.fan.sectors(control['final_links']), oracle.fan.sectors(together['final_links'])))
                if together['conversion_events'] == 0 and together['final_links'] != control['final_links']:
                    raise ValueError('no-conversion run differs from the matched control')
            runs += raw
            print('Family', family['id'], 'trial', trial, 'conversions:', [r['conversion_events'] for r in raw if r['mode'] == 'combined'], flush=True)
    return {'schema': 1, 'side': side, 'attempts': attempts, 'stride': stride, 'trials': trials, 'seed': seed,
            'conditions': CONDITIONS, 'initial_links': initial, 'schedules': schedules, 'families': families,
            'rule_id_convention': '0: shared-edge vacancy transport; j>0: family.rule_ids[j-1] on a three-face fan',
            'schedule': 'uniform random generator-placement attempts; every one of 48 reaction rules and the transport rule has equal attempt weight',
            'scope': 'finite classical multirule gauge dynamics on a supplied torus; no calibrated physical time, mass, or quantum amplitudes',
            'runs': runs}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--side', type=int, default=6)
    parser.add_argument('--attempts', type=int, default=50000)
    parser.add_argument('--stride', type=int, default=1000)
    parser.add_argument('--trials', type=int, default=2)
    parser.add_argument('--seed', type=int, default=6291871)
    parser.add_argument('--output', type=Path, default=Path('out/channel-banks.json'))
    args = parser.parse_args()
    if not 3 <= args.side <= 24 or not 0 < args.attempts <= 1000000 or args.stride <= 0 or args.attempts % args.stride or not 0 < args.trials <= 32:
        parser.error('invalid bounded experiment dimensions')
    channels = json.loads(Path('data/d4-triple-channels.json').read_text())
    expected = analyze_channels(json.loads(Path('data/d4-triple-rule-search.json').read_text()))
    if channels != json.loads(json.dumps(expected)):
        raise ValueError('channel artifact differs from the independent symmetry and charge derivation')
    result = experiment(channels, args.side, args.attempts, args.stride, args.trials, args.seed)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, separators=(',', ':'))+'\n')
    print('Verified bank runs:', len(result['runs']), 'attempted clock ticks:', len(result['runs'])*args.attempts)


if __name__ == '__main__':
    main()
