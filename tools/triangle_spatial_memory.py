#!/usr/bin/env python3
"""Translation-resolved gauge memory, whole-charge conditioning, and spectral bounds.

Fourier phases are diagnostic coordinates on the supplied periodic mesh, not
physical quantum amplitudes. No microscopic rule or update clock is changed.
"""
import argparse
import json
import math
import random
from fractions import Fraction
from pathlib import Path

import numpy as np

from triangle_charge_response import ChargeMarginal, coefficients
from triangle_elastic_scattering import ElasticExperiment
from triangle_reference import AxialReference, field_weight, fraction, reference
from triangle_relational_memory import RelationalObserver, conditional_marginals, replay_samples


def conditioned_coefficients(n1, n2):
    """P(distinct RR | entire q field), P(active RRR | entire q field).

    For each specified charge field with at least one reflection, its uniform
    full-S3 weight is 12*3**n1*2**n2 - 12*[n2=0]. The handle character convolution
    gives the two numerators below; no particular placement of charges enters.
    """
    if type(n1) is not int or type(n2) is not int or n1 < 2 or n1 % 2 or n2 < 0:
        raise ValueError('requires a nonempty reflection-present torus charge field')
    denominator = field_weight(0, n1, n2, 'nonabelian_reflections')
    pair_rotation = 12*3**(n1-2)*2**n2-3*int(n1 == 2)*(-1)**n2
    beta = Fraction(6*pair_rotation, denominator)
    alpha = Fraction(144*3**(n1-3)*2**n2, denominator) if n1 >= 4 else Fraction()
    return beta, alpha


def conditioned_variances(faces, charge):
    ref = reference(faces, charge, 'nonabelian_reflections')
    variances = [Fraction(*ref['single_face_charge_variance']), Fraction(), Fraction()]
    for row in ref['population_counts']:
        _, n1, n2 = row['populations']
        weight = Fraction(row['canonical_connections'], ref['canonical_connections'])
        beta, alpha = conditioned_coefficients(n1, n2)
        for index, k, mean in ((1, 2, beta), (2, 3, alpha)):
            probability = Fraction(math.prod(n1-i for i in range(k)), math.prod(faces-i for i in range(k))) if n1 >= k else Fraction()
            variances[index] += weight*probability*mean*(1-mean)
    return variances


class Translations:
    """Free translation orbits constructed from actual oriented support specs."""
    def __init__(self, experiment):
        self.x, self.side = experiment, experiment.e.geometry.side
        e = experiment.e
        self.specs = [[(tuple(face),) for face in e.geometry.faces],
                      [(tuple(a), tuple(b)) for a, b, _, _ in e.geometry.pair_specs],
                      [tuple(map(tuple, fan)) for fan in e.factor.oracle.fan.specs]]
        self.layouts = []
        for kind, specs in enumerate(self.specs):
            keys = {self.key(spec, kind): i for i, spec in enumerate(specs)}
            if len(keys) != len(specs):
                raise ValueError('support keys are not unique')
            seen, layouts = set(), []
            for seed, spec in enumerate(specs):
                if seed in seen:
                    continue
                layout = np.array([[keys[self.key(self.shift_spec(spec, dx, dy), kind)]
                                    for dx in range(self.side)] for dy in range(self.side)])
                orbit = set(layout.flat)
                if len(orbit) != self.side**2 or orbit & seen:
                    raise ValueError('support translation action is not free')
                seen.update(orbit); layouts.append(layout)
            self.layouts.append(np.array(layouts))
        if [len(a) for a in self.layouts] != [2, 6, 6]:
            raise ValueError('unexpected translation-channel multiplicities')

    @staticmethod
    def key(spec, kind):
        if kind:
            return spec
        face = spec[0][:-1]
        return (min(face[i:]+face[:i] for i in range(len(face))),)

    def vertex(self, v, dx, dy):
        s = self.side
        return s*((v//s+dy) % s)+(v % s+dx) % s

    def shift_spec(self, spec, dx, dy):
        return tuple(tuple(self.vertex(v, dx, dy) for v in path) for path in spec)

    def links(self, raw, dx, dy):
        e, g = self.x.e, self.x.e.geometry.group
        result = [0]*len(raw)
        for (u, v), h in zip(e.geometry.edges, raw):
            a, b = self.vertex(u, dx, dy), self.vertex(v, dx, dy)
            edge = e.geometry.ids[tuple(sorted((a, b)))]
            result[edge] = h if a < b else g.inv[h]
        return result

    def charge_blocks(self):
        """Fourier symbol of the actual dual Laplacian, retaining both phases."""
        layout, s = self.layouts[0], self.side
        positions = {int(face): (a, y, x) for a in range(2) for y in range(s) for x in range(s) for face in [layout[a, y, x]]}
        neighbors = [set() for _ in range(2*s*s)]
        for a, b in self.x.p.dual:
            neighbors[a].add(b); neighbors[b].add(a)
        blocks = np.zeros((s, s, 2, 2), dtype=complex)
        for ky in range(s):
            for kx in range(s):
                for a in range(2):
                    blocks[ky, kx, a, a] = 3
                    for face in neighbors[int(layout[a, 0, 0])]:
                        b, dy, dx = positions[face]
                        blocks[ky, kx, a, b] -= np.exp(2j*np.pi*(kx*dx+ky*dy)/s)
                if not np.allclose(blocks[ky, kx], blocks[ky, kx].conj().T, atol=1e-12):
                    raise ValueError('dual Fourier block is not Hermitian')
        eigenvalues, eigenvectors = np.linalg.eigh(blocks)
        return blocks, eigenvalues, eigenvectors


class SpatialObserver(RelationalObserver):
    def __init__(self, experiment):
        faces = experiment.geometry.faces
        super().__init__(experiment, conditional_marginals(experiment.e.geometry.group, experiment.e.charges, faces, faces))
        self.variances = conditioned_variances(faces, faces)
        self.translations = Translations(experiment)
        self.blocks, self.eigenvalues, self.eigenvectors = self.translations.charge_blocks()
        levels = []
        for value in sorted(self.eigenvalues[:, :, 0].flat):
            if value > 1e-10 and not any(abs(value-v) < 1e-10 for v in levels):
                levels.append(float(value))
        self.shells = [[list(map(int, pair)) for pair in np.argwhere(abs(self.eigenvalues[:, :, 0]-v) < 1e-10)] for v in levels[:2]]
        self.chi = float(self.variances[0]*faces/(faces-1))
        _, _, _, k1, k2 = coefficients(ChargeMarginal(faces, faces, 'nonabelian_reflections'))
        self.rates = [float(k1)*v+float(k2)*(6*v-v*v) for v in [*levels[:2], 6.0]]

    def measure(self, links):
        raw = super().measure(links)
        n1, n2 = int(np.sum(raw[0] == 1)), int(np.sum(raw[0] == 2))
        beta, alpha = conditioned_coefficients(n1, n2)
        values = [raw[0]-1, raw[1][:, 0]-float(beta)*raw[1][:, 1], raw[2][:, 0]-float(alpha)*raw[2][:, 1]]
        fields = [values[layout] for values, layout in zip(values, self.translations.layouts)]
        fourier = [np.fft.fft2(field, axes=(-2, -1)) for field in fields]
        return {'populations': [len(raw[0])-n1-n2, n1, n2], 'fields': fields, 'fourier': fourier}

    def correlate(self, first, later):
        spatial, spectra = [], []
        for a, b, variance in zip(first['fourier'], later['fourier'], self.variances):
            spectrum = np.sum(a.conj()*b, axis=0)/(a.size*float(variance))
            spectra.append(spectrum)
            spatial.append(np.fft.ifft2(spectrum).real)
        # ifft(conj(fft(a))*fft(b)) is the SUM over positions, not its mean.
        # The division by N is already included in spectrum above.
        spatial = np.array(spatial)
        q0, qt = [np.einsum('yxab,ayx->byx', self.eigenvectors.conj(), sample['fourier'][0]) for sample in (first, later)]
        qproducts = q0.conj()*qt/(self.x.e.geometry.size*self.chi)
        shell_values = [np.mean([qproducts[0, y, x] for y, x in shell]) for shell in self.shells]
        shell_values.append(qproducts[1, 0, 0])
        for spectrum in spectra[1:]:
            shell_values.extend(np.mean([spectrum[y, x] for y, x in shell]) for shell in self.shells)
        y, x = self.shells[0][0]
        imaginary = [qproducts[0, y, x].imag, spectra[1][y, x].imag, spectra[2][y, x].imag]
        # The conserved total charge is the acoustic k=0 mode. Keeping both
        # face phases prevents confusing the nonconserved optical mode with it.
        if abs(qproducts[0, 0, 0]) > 1e-20:
            raise ValueError('conserved charge zero mode did not vanish')
        return spatial, [float(v.real) for v in shell_values], list(map(float, imaginary))


class Moments:
    def __init__(self):
        self.count = 0
        self.mean = None
        self.m2 = None

    def add(self, values):
        values = np.asarray(values)
        self.count += 1
        if self.mean is None:
            self.mean, self.m2 = values.copy(), np.zeros_like(values)
        else:
            delta = values-self.mean
            self.mean += delta/self.count
            self.m2 += delta*(values-self.mean)

    def result(self):
        if self.count < 2:
            raise ValueError('at least two independent replicates required')
        return {'mean': self.mean.tolist(), 'standard_error': np.sqrt(self.m2/(self.count-1)/self.count).tolist()}


def trial(observer, control, sampler, index, times):
    x, side = observer.x, observer.x.e.geometry.side
    initial_seed, schedule_seed = 18300000+side*1000+index, 28400000+side*1000+index
    initial = sampler.sample(random.Random(initial_seed), 'nonabelian_reflections')['links']
    ticks = [t*x.geometry.faces for t in times]
    rng = random.Random(schedule_seed)
    schedule = [rng.randrange(x.operator_count) for _ in range(ticks[-1])]
    records, maps = {}, {}
    for name, experiment in (('baseline_padded', control), ('elastic', x)):
        samples, check = replay_samples(experiment, initial, schedule, ticks, observer)
        measured = [observer.correlate(samples[0], sample) for sample in samples]
        maps[name] = np.array([r[0] for r in measured])
        records[name] = {'shell_correlations': [r[1] for r in measured], 'oriented_mode_imaginary_parts': [r[2] for r in measured],
                         'populations': [s['populations'] for s in samples], **check}
    return {'index': index, 'initial_seed': initial_seed, 'schedule_seed': schedule_seed, 'runs': records}, maps


def whole_mesh(side, bank, elastic, trials, times):
    x = ElasticExperiment(side, bank, elastic)
    control = ElasticExperiment(side, bank, [list(range(36))]*2)
    observer, sampler = SpatialObserver(x), AxialReference(side, bank, x.geometry.faces)
    spatial = {name: Moments() for name in ('baseline_padded', 'elastic', 'paired_elastic_minus_baseline')}
    rows = []
    for i in range(trials):
        record, maps = trial(observer, control, sampler, i, times)
        rows.append(record)
        maps['paired_elastic_minus_baseline'] = maps['elastic']-maps['baseline_padded']
        for name, values in maps.items():
            spatial[name].add(values)
    def estimates(key):
        arrays = {name: np.array([row['runs'][name][key] for row in rows]) for name in ('baseline_padded', 'elastic')}
        arrays['paired_elastic_minus_baseline'] = arrays['elastic']-arrays['baseline_padded']
        return {name: {'mean': a.mean(axis=0).tolist(), 'standard_error': (a.std(axis=0, ddof=1)/math.sqrt(trials)).tolist()} for name, a in arrays.items()}
    m = x.operator_count
    bounds = [[(1-rate/m)**(t*x.geometry.faces) for rate in observer.rates] for t in times]
    return {'side': side, 'faces': x.geometry.faces, 'trials': trials, 'times_in_attempts_per_face': list(times),
            'translation_channel_counts': [len(a) for a in observer.translations.layouts],
            'field_variances': [fraction(v) for v in observer.variances],
            'shell_wavevectors_yx': observer.shells,
            'acoustic_shell_laplacian_eigenvalues': [float(observer.eigenvalues[tuple(shell[0])][0]) for shell in observer.shells],
            'shell_observable_order': ['charge_acoustic_1', 'charge_acoustic_2', 'charge_optical_k0', 'pair_shell_1', 'pair_shell_2', 'fan_shell_1', 'fan_shell_2'],
            'charge_initial_generator_rates': observer.rates, 'charge_jensen_lower_bounds': bounds,
            'replicates': rows, 'spatial_correlations': {name: moment.result() for name, moment in spatial.items()},
            'shell_estimates': estimates('shell_correlations'), 'imaginary_estimates': estimates('oriented_mode_imaginary_parts'),
            'scope': 'Stationary, full-S3 reflection-present Q=F sector. Fourier channels are diagnostic translations, not emergent Euclidean coordinates or quantum amplitudes. Errors count independent trajectory pairs only.'}


def audit(trials=64, sides=(6, 12), times=(0, 1, 2, 3, 4, 8, 16, 32, 64, 128)):
    if type(trials) is not int or trials < 2:
        raise ValueError('at least two independent trials required')
    bank = json.loads(Path('data/triangle-feedback.json').read_text())['bank']
    elastic = json.loads(Path('data/triangle-elastic-scattering.json').read_text())['pair_census']['elastic_tables']
    return {'scope': 'New observers of the existing two banks; no changes to primitive physics.',
            'whole_mesh': [whole_mesh(side, bank, elastic, trials, times) for side in sides]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--trials', type=int, default=64)
    parser.add_argument('--sides', type=int, nargs='+', default=[6, 12])
    parser.add_argument('--output', type=Path, default=Path('out/triangle-spatial-memory.json'))
    args = parser.parse_args()
    result = audit(args.trials, args.sides)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2)+'\n')
    for row in result['whole_mesh']:
        print('side', row['side'], 'charge eigenvalues', row['acoustic_shell_laplacian_eigenvalues'])
        print('elastic shell means', row['shell_estimates']['elastic']['mean'])
    print('Wrote', args.output)


if __name__ == '__main__':
    main()
