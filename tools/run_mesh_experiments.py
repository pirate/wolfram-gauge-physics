#!/usr/bin/env python3
"""Evolve the same microscopic gauge table on shared-link meshes and audit local closure."""
import argparse
import itertools
import json
import subprocess
import math
import statistics
from collections import defaultdict
from pathlib import Path
from charge_quotient import analyze
from screen_braid_rules import nullspace


def analyze_result(result, table):
    group = [tuple(p) for p in result['group_permutations']]
    n = len(group)
    identity = group.index(tuple(range(len(group[0]))))
    multiply = [[group.index(tuple(a[b[i]] for i in range(len(a)))) for b in group] for a in group]
    inverse = [next(b for b in range(n) if multiply[a][b] == identity) for a in range(n)]
    def conjugate(g, a):
        return multiply[g][multiply[a][inverse[g]]]
    sectors = [min(conjugate(g, a) for g in range(n)) for a in range(n)]
    labels = sorted(set(sectors)-{identity})
    targets, representatives = defaultdict(set), {}
    equations = set()
    witness = None
    rows = result['local_census']['transitions']
    if len(rows) != n**4 or {tuple(row[:4]) for row in rows} != set(itertools.product(range(n), repeat=4)):
        raise ValueError('local census is not exhaustive in the four independent gauge-fixed holonomies')
    for row in rows:
        a, b, c, d = row[:4]
        aa, bb = divmod(table[a*n+conjugate(inverse[a], b)], n)
        bnew = conjugate(a, bb)
        link = multiply[aa][inverse[a]]
        predicted = [multiply[a][link], bnew, multiply[inverse[link]][c],
                     multiply[inverse[bnew]][multiply[b][d]]]
        if row[4:] != predicted:
            raise ValueError('microscopic four-face output disagrees with independent based-loop algebra')
        before, after = tuple(sectors[x] for x in row[:4]), tuple(sectors[x] for x in row[4:])
        targets[before].add(after)
        if before in representatives and representatives[before][1] != after and witness is None:
            prior, prior_after = representatives[before]
            gauge_equivalent = any(tuple(conjugate(g, x) for x in prior[:4]) == tuple(row[:4]) for g in range(n))
            if gauge_equivalent:
                raise ValueError('nonclosure witness is gauge equivalent: covariance or comparison is wrong')
            witness = {'initial_face_sectors': before, 'first': prior, 'second': row,
                       'first_output_sectors': prior_after, 'second_output_sectors': after,
                       'initial_states_gauge_equivalent': False}
        representatives.setdefault(before, (row, after))
        equations.add(tuple(before.count(label)-after.count(label) for label in labels))
    basis = nullspace(sorted(equations), len(labels))
    result['local_census'].update({
        'independent_algebra_checks': len(rows), 'element_sectors': sectors,
        'nonidentity_sector_labels': labels, 'incoming_sector_tuples': len(targets),
        'nonclosed_incoming_sector_tuples': sum(len(values)>1 for values in targets.values()),
        'maximum_output_sector_tuples': max(map(len, targets.values())),
        'nonclosure_witness': witness, 'charge_convention': 'rational one-face additive class functions, q(identity)=0',
        'charge_equations': sorted(equations), 'additive_charge_basis': basis})
    for run in result['runs']:
        for row in run['trajectory']:
            if sum(row[5:]) != result['faces'] or row[1] != result['faces']-row[5+identity]:
                raise ValueError('face histogram disagrees with active-face count')
            if row[2] > 2*row[0] or (run['condition'] == 'flat' and row[1]):
                raise ValueError('vacuum or local propagation check failed')
        if not run['exact_inverse_replay'] or not run['reverse_histogram_echo']:
            raise ValueError('microscopic inverse replay failed')
    reference = {}
    seed_elements = result['seed_group_elements']
    for condition in ('flat', 'reflection', 'rotation', 'noncommuting_pair', 'random_links'):
        generators = ([seed_elements[condition]] if condition in seed_elements else
                      list(seed_elements.values()) if condition == 'noncommuting_pair' else list(range(n)))
        if condition == 'flat':
            subgroup = [identity]  # Flat fixing is stronger than the stabilizer bound.
        else:
            centralizer = [g for g in range(n) if all(multiply[g][a] == multiply[a][g] for a in generators)]
            subgroup = [a for a in range(n) if all(multiply[g][a] == multiply[a][g] for g in centralizer)]
        if any(table[a*n+b]//n not in subgroup or table[a*n+b]%n not in subgroup
               for a in subgroup for b in subgroup):
            raise ValueError('proposed invariant subgroup is not closed under the actual table')
        probability = [sum(sectors[a] == label for a in subgroup)/len(subgroup) for label in range(n)]
        runs = [r for r in result['runs'] if r['condition'] == condition]
        for run in runs:
            if any(count and probability[label] == 0 for row in run['trajectory'] for label, count in enumerate(row[5:])):
                raise ValueError('face sectors escaped the proven invariant subgroup')
        fractions = [r['trajectory'][-1][1]/result['faces'] for r in runs]
        def entropy(row):
            return -sum((count/result['faces'])*math.log(count/result['faces']) for count in row[5:] if count)
        reference[condition] = {
            'invariant_subgroup_elements': subgroup,
            'stationary_sector_probabilities': probability,
            'stationary_active_fraction': 1-probability[identity],
            'stationary_sector_entropy_nats': -sum(p*math.log(p) for p in probability if p),
            'final_mean_active_fraction': statistics.fmean(fractions) if fractions else None,
            'final_active_fraction_se_across_trials': statistics.stdev(fractions)/math.sqrt(len(fractions)) if len(fractions)>1 else None,
            'final_mean_sector_entropy_nats': statistics.fmean(entropy(r['trajectory'][-1]) for r in runs) if runs else None}
    result['stationary_reference'] = reference
    result['reference_caveat'] = ('Uniform independent links in the proven invariant subgroup form an exact stationary labeled-state measure. '
                                  'Agreement of local histograms is not proof of ergodicity, thermalization, or uniform physical gauge orbits.')
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--binary', type=Path, default=Path('build/wgphysics_mesh_experiments'))
    parser.add_argument('--side', type=int, default=12)
    parser.add_argument('--layers', type=int, default=256)
    parser.add_argument('--trials', type=int, default=8)
    parser.add_argument('--seed', type=int, default=819031)
    parser.add_argument('--output', type=Path, default=Path('out/shared-mesh.json'))
    parser.add_argument('--reanalyze', type=Path)
    args = parser.parse_args()
    census = json.loads(Path('data/d4-involution-braid-search.json').read_text())
    screen = json.loads(Path('data/d4-braid-rule-screen.json').read_text())
    quotient = analyze(census, screen)
    selected = next(q for q in quotient['rules'] if q['exclusion_distinguished_symbol'] == 0)
    table = census['solutions'][selected['solution_id']]['table']
    if args.reanalyze:
        result = json.loads(args.reanalyze.read_text())
        if result['solution_id'] != selected['solution_id']:
            raise ValueError('recorded experiment uses a different microscopic law')
    else:
        completed = subprocess.run([str(args.binary), '--side', str(args.side), '--layers', str(args.layers),
                                    '--trials', str(args.trials), '--seed', str(args.seed)],
                                   input=' '.join(map(str, table))+'\n', text=True,
                                   stdout=subprocess.PIPE, check=True, timeout=1800)
        result = json.loads(completed.stdout)
    result.update({'solution_id': selected['solution_id'],
                   'geometry': 'supplied periodic triangulated two-dimensional lattice; no embedding used by updates',
                   'clock': 'one conflict-free link read/write matching layer, not calibrated physical time',
                   'schedule': 'all ordered adjacent-face closing-link choices; seeded ordering, greedy conflict coloring; fixed cyclic replay',
                   'randomness': 'schedule choice per trial; initial links only for random_links; no random update law',
                   'boundary_control': 'local four-face census writes internal links only; all exterior links fixed',
                   'global_caveat': 'no inferred energy, electric charge, physical amplitudes, or emergent spatial dimension'})
    analyze_result(result, table)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2)+'\n')
    local = result['local_census']
    print(f'{local["nonclosed_incoming_sector_tuples"]}/{local["incoming_sector_tuples"]} incoming sector tuples have multiple outputs; '
          f'additive class-charge rank {len(local["additive_charge_basis"])}')
    for condition in sorted({r['condition'] for r in result['runs']}):
        runs = [r for r in result['runs'] if r['condition'] == condition]
        print(condition, 'initial active', [r['trajectory'][0][1] for r in runs],
              'final active', [r['trajectory'][-1][1] for r in runs])


if __name__ == '__main__':
    main()
