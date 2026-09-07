#!/usr/bin/env python3
"""Run the microscopic gauge-table experiment and estimate transport with block CIs."""
import argparse
import json
import math
import random
import statistics
import subprocess
from pathlib import Path
from charge_quotient import analyze
from exact_tag_distribution import tag_distribution


def moments(blocks, block_size):
    count = len(blocks)*block_size
    results = []
    for t in range(len(blocks[0])):
        raw = [sum(block[t][k] for block in blocks)/count for k in range(4)]
        mean, second, third, fourth = raw
        variance = max(0, second-mean*mean)
        central_fourth = fourth-4*mean*third+6*mean*mean*second-3*mean**4
        results.append({'mean': mean, 'msd': second, 'variance': variance,
                        'kurtosis': central_fourth/(variance*variance) if variance else None})
    return results


def slope(x, y):
    mx, my = statistics.fmean(x), statistics.fmean(y)
    return sum((a-mx)*(b-my) for a, b in zip(x, y))/sum((a-mx)**2 for a in x)


def fits(times, measured):
    chosen = list(range(len(times)-3, len(times)))
    values = [measured[i]['variance'] for i in chosen]
    windows = [times[i] for i in chosen]
    return (slope(windows, values)/2,
            slope([math.log(t) for t in windows], [math.log(v) for v in values]) if min(values)>0 else None)


def interval(values):
    values = sorted(v for v in values if v is not None)
    if not values: return None
    def quantile(p):
        i = (len(values)-1)*p
        return values[int(i)]*(1-(i % 1))+values[math.ceil(i)]*(i % 1)
    return [quantile(.025), quantile(.975)]


def validate_recording(result):
    """Cross-check independently accumulated histograms and block raw moments."""
    layers, cells = result['layers'], result['cells']
    if layers < 16 or layers % 2 or cells % 4 or cells < 4*layers+16:
        raise ValueError('recording violates clock or causal-boundary guard')
    for ensemble in result['ensembles']:
        times, blocks = ensemble['times'], ensemble['block_raw_moment_sums']
        histograms = ensemble['displacement_histograms']
        if (times != sorted(set(times)) or times[0] != 0 or times[-1] != layers
                or any(t % 2 for t in times) or len(blocks) != 32
                or len(histograms) != len(times)
                or any(len(block) != len(times) for block in blocks)
                or ensemble['trials'] != len(blocks)*ensemble['block_size']):
            raise ValueError('invalid sampling metadata')
        for index, (tick, histogram) in enumerate(zip(times, histograms)):
            if (len(histogram) != 2*layers+1 or sum(histogram) != ensemble['trials']
                    or any(count < 0 or int(count) != count for count in histogram)
                    or any(count and abs(i-layers) > tick for i, count in enumerate(histogram))):
                raise ValueError('invalid histogram or superluminal tag in layer/cell units')
            for k in range(4):
                measured = sum(block[index][k] for block in blocks)
                independent = sum(count*(i-layers)**(k+1) for i, count in enumerate(histogram))
                if not math.isclose(measured, independent, rel_tol=1e-12, abs_tol=1e-8):
                    raise ValueError('histogram and raw trajectory moments disagree')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--binary', type=Path, default=Path('build/wgphysics_diffusion'))
    parser.add_argument('--output', type=Path, default=Path('out/diffusion.json'))
    parser.add_argument('--trials', type=int, default=4096)
    parser.add_argument('--cells', type=int, default=1040)
    parser.add_argument('--layers', type=int, default=256)
    parser.add_argument('--seed', type=int, default=731291)
    parser.add_argument('--reanalyze', type=Path, help='reuse recorded raw moments and histograms; do not rerun evolution')
    args = parser.parse_args()
    census = json.loads(Path('data/d4-involution-braid-search.json').read_text())
    screen = json.loads(Path('data/d4-braid-rule-screen.json').read_text())
    quotient = analyze(census, screen)
    selected = next(q for q in quotient['rules'] if q['exclusion_distinguished_symbol'] == 0)
    table = census['solutions'][selected['solution_id']]['table']
    payload = ' '.join(map(str, table+selected['element_symbols']))+'\n'
    if args.reanalyze:
        result = json.loads(args.reanalyze.read_text())
        if result['solution_id'] != selected['solution_id'] or result['element_symbols'] != selected['element_symbols']:
            raise ValueError('recorded experiment uses a different microscopic rule or charge projection')
    else:
        completed = subprocess.run([str(args.binary), '--trials', str(args.trials), '--cells', str(args.cells),
                                    '--layers', str(args.layers), '--seed', str(args.seed)],
                                   input=payload, text=True, stdout=subprocess.PIPE, check=True, timeout=1800)
        result = json.loads(completed.stdout)
    validate_recording(result)
    result.update({'solution_id': selected['solution_id'], 'element_symbols': selected['element_symbols'],
                   'initial_ensemble': 'independent occupied/vacant charge labels; uniform representatives within each label; one forced tag',
                   'sampling': 'randomness only in initial conditions; deterministic alternating matching updates; phases balanced',
                   'verification': 'every microscopic update compared with exact quotient; initial/final transported total checked',
                   'comparison': 'D=(1/rho-1)/2 in single-layer units; asymptotic benchmark from the exactly matched exclusion model'})
    generator = random.Random(114733)
    for ensemble in result['ensembles']:
        measured = moments(ensemble['block_raw_moment_sums'], ensemble['block_size'])
        diffusion, exponent = fits(ensemble['times'], measured)
        bootstrap = []
        blocks = ensemble['block_raw_moment_sums']
        for _ in range(1000):
            resampled = [blocks[generator.randrange(len(blocks))] for _ in blocks]
            bootstrap.append(fits(ensemble['times'], moments(resampled, ensemble['block_size'])))
        ensemble['observables'] = measured
        rho = ensemble['density']
        ensemble['late_fit'] = {'times': ensemble['times'][-3:],
                                'variance_slope_over_two': diffusion,
                                'diffusion': diffusion if 0 < rho < 1 else None,
                                'variance_slope_over_two_ci95': interval([x[0] for x in bootstrap]),
                                'diffusion_ci95': interval([x[0] for x in bootstrap]) if 0 < rho < 1 else None,
                                'variance_exponent': exponent, 'variance_exponent_ci95': interval([x[1] for x in bootstrap]),
                                'predicted_diffusion': (1/rho-1)/2 if rho else None,
                                'control': 'ballistic_empty' if rho == 0 else 'blocked_full' if rho == 1 else 'diffusion_test',
                                'ci_method': '1000 bootstrap resamples of 32 independent trajectory blocks'}
        ensemble['microscopic_updates'] = ensemble['trials']*result['layers']*(result['cells']-1)//2
        exact = tag_distribution(result['layers'], rho)
        finite = [tag_distribution(t, rho) for t in ensemble['times']]
        finite_slope, finite_exponent = fits(ensemble['times'], finite)
        ensemble['finite_time_prediction'] = {
            'variances': [row['variance'] for row in finite],
            'late_variance_slope_over_two': finite_slope,
            'late_variance_exponent': finite_exponent}
        observed_cdf, total = [], 0
        for count in ensemble['displacement_histograms'][-1]:
            total += count
            observed_cdf.append(total/ensemble['trials'])
        ensemble['exact_late_prediction'] = {
            'mean': exact['mean'], 'variance': exact['variance'], 'pmf': exact['pmf'],
            'cdf_max_error': max(abs(a-b) for a, b in zip(observed_cdf, exact['cdf']))}
        if 0 < rho < 1:
            predicted = (1/rho-1)/2
            cumulative = 0
            error = 0
            for i, count in enumerate(ensemble['displacement_histograms'][-1]):
                cumulative += count
                edge = i-result['layers']+.5
                expected = (1+math.erf(edge/math.sqrt(4*predicted*result['layers'])))/2
                error = max(error, abs(cumulative/ensemble['trials']-expected))
            ensemble['late_fit']['gaussian_cdf_max_error'] = error
        print(f'rho={rho}: half variance slope={diffusion:.5g}, alpha={exponent}, '
              f'diffusion CI={ensemble["late_fit"]["diffusion_ci95"]}')
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2)+'\n')


if __name__ == '__main__':
    main()
