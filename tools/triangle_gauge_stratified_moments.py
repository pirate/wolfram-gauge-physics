#!/usr/bin/env python3
"""Stationary gauge-current energies resolved by microscopic move family.

Batching evaluates the same finite-group path products, not an approximate
evolution law. Stratification changes only the estimator of a stationary
Dirichlet sum. All reference states are independent exact canonical draws.
"""
import json
import random
from pathlib import Path

import numpy as np

from triangle_charge_response import ChargeMarginal, coefficients
from triangle_gauge_current_corrector import GaugeCurrent
from triangle_reference import AxialReference


class HiddenBatch:
    def __init__(self, current):
        self.g = current
        oracle = current.x.e.factor.oracle
        group = oracle.group
        self.mul, self.inv = np.array(group.mul), np.array(group.inv)
        self.pair_paths = np.array(oracle.pairs)
        self.fan_paths = np.array([p[::-1] for p in oracle.fan.patches])
        self.face_paths = np.array(oracle.fan.face_paths)
        self.tables = [np.array(table) for table in current.x.tables]
        self.charges = np.array(current.x.e.charges)
        counts = [len(fs) for fs in current.affected_fans]
        assert len(set(counts)) == 1
        self.fans_per_pair = counts[0]
        self.pair_ids = np.repeat(np.arange(current.x.supports), counts)
        self.fan_ids = np.array([f for fs in current.affected_fans for f in fs])
        self.changed_paths = self.fan_paths[self.fan_ids]
        self.write_edges = self.pair_paths[:, 0, 0, 0]

    def transport(self, raw, paths, replacement=None):
        result = np.zeros((len(raw), *paths.shape[:-2]), dtype=np.int64)
        for step in range(paths.shape[-2]):
            edges, reverse = paths[..., step, 0], paths[..., step, 1]
            values = raw[:, edges]
            if replacement is not None:
                written, new = replacement
                values = np.where(edges == written, new, values)
            values = np.where(reverse, self.inv[values], values)
            result = self.mul[values, result]
        return result

    def codes(self, raw):
        pair = self.transport(raw, self.pair_paths)
        fan = self.transport(raw, self.fan_paths)
        return 6*pair[:, :, 0]+pair[:, :, 1], 36*fan[:, :, 0]+6*fan[:, :, 1]+fan[:, :, 2]

    @staticmethod
    def activity(values):
        return (values[:, :, 0] == values[:, :, 1]) != (values[:, :, 1] == values[:, :, 2])

    def evaluate(self, states):
        raw = np.asarray(states, dtype=np.int64)
        q = self.charges[self.transport(raw, self.face_paths)]
        weights = np.array([self.g.weights(row) for row in q])
        initial = self.activity(self.transport(raw, self.fan_paths))
        pair = self.transport(raw, self.pair_paths)
        codes = 6*pair[:, :, 0]+pair[:, :, 1]
        prefix = self.transport(raw, self.pair_paths[:, 0, 1:])
        total = np.zeros(len(raw), dtype=np.int64)
        squares = total.copy()
        for rule in (1, 2):
            target_first = self.tables[rule][codes]//6
            new = self.mul[self.inv[prefix], target_first]
            new = np.where(self.pair_paths[:, 0, 0, 1], self.inv[new], new)
            changed = self.transport(raw, self.changed_paths,
                                     (self.write_edges[self.pair_ids, None],
                                      new[:, self.pair_ids, None]))
            difference = weights[:, self.fan_ids]*(
                initial[:, self.fan_ids].astype(np.int64)-self.activity(changed))
            delta = difference.reshape(len(raw), self.g.x.supports, self.fans_per_pair).sum(axis=2)
            total += delta.sum(axis=1)
            squares += (delta*delta).sum(axis=1)
        assert not np.any(squares % 2)
        return total, squares//2

    def active_families(self, raw):
        pair, fan = self.codes(np.array([raw]))
        pair, fan = pair[0], fan[0]
        p = self.g.x.supports
        vacancy = np.flatnonzero(self.tables[0][pair] != pair).tolist()
        elastic = [r*p+i for r in (1, 2) for i in np.flatnonzero(self.tables[r][pair] != pair)]
        birth, death = [], []
        fan_q = self.charges[np.array([fan//36, fan//6 % 6, fan % 6]).T]
        all_reflections = np.all(fan_q == 1, axis=1)
        for rule in range(3, 15):
            for i in np.flatnonzero(self.tables[rule][fan] != fan):
                (birth if all_reflections[i] else death).append(rule*p+int(i))
        return vacancy, elastic, birth, death


def measure(side=6, draws=4096, seed=158329001, samples_per_family=1):
    current = GaugeCurrent(side)
    batch = HiddenBatch(current)
    bank = json.loads(Path('data/triangle-feedback.json').read_text())['bank']
    f = current.x.geometry.faces
    reference = AxialReference(side, bank, f)
    marginal = ChargeMarginal(f, f, 'nonabelian_reflections')
    v = marginal.expectation(current.variance_counts.items())
    w = marginal.expectation(current.twice_energy_counts.items())/2
    a, b, _, _, _ = coefficients(marginal)
    bare = float((a+10*b)/135)
    base = bare-float(2*v*v/(810*f*w))
    c = float(v/w)
    # Chosen before this sample: half the coefficient from the old 128-draw
    # pilot. Neither this coefficient nor c is optimized on these observations.
    d = -2.9208964015355687e-5/2
    rng = random.Random(seed)
    rows = []
    for start in range(0, draws, 32):
        states, moved, selected = [], [], []
        for _ in range(min(32, draws-start)):
            raw = reference.sample(rng, 'nonabelian_reflections')['links']
            states.append(raw)
            groups = []
            for family in batch.active_families(raw):
                indices = []
                for _ in range(samples_per_family if family else 0):
                    after = raw[:]
                    current.x.apply(after, int(rng.choice(family)))
                    indices.append(len(moved))
                    moved.append(after)
                groups.append((len(family), indices))
            selected.append(groups)
        h, s = batch.evaluate(states)
        other = np.concatenate([batch.evaluate(moved[i:i+64])[0]
                                for i in range(0, len(moved), 64)]) if moved else np.array([])
        for i, groups in enumerate(selected):
            energies = [count*float(np.mean((other[indices]-h[i])**2))/(2*f)
                        if indices else 0.0 for count, indices in groups]
            rows.append([float(s[i])/f, *energies, float(int(h[i])**2)/f])
        if len(rows) % 256 == 0 or len(rows) == draws:
            print(json.dumps({'draws': len(rows), 'mean_S_family_T_and_h_squared_per_face':
                              np.mean(rows, axis=0).tolist()}), flush=True)
    observed = np.array(rows)
    total_energy = observed[:, 1:5].sum(axis=1)
    changes = (4*c*d*observed[:, 0]+2*d*d*total_energy)/810
    birth_death = observed[:, 3]-observed[:, 4]
    se = lambda values: float(np.std(values, ddof=1)/np.sqrt(draws))
    return {'side': side, 'faces': f, 'draws': draws, 'seed': seed,
            'samples_per_family_per_state': samples_per_family,
            'families': ['vacancy', 'elastic', 'birth', 'death'],
            'S_per_face': float(observed[:, 0].mean()), 'S_standard_error': se(observed[:, 0]),
            'family_T_per_face': observed[:, 1:5].mean(axis=0).tolist(),
            'family_T_standard_errors': [se(observed[:, k]) for k in range(1, 5)],
            'T_per_face': float(total_energy.mean()),
            'T_standard_error': se(total_energy),
            'h_squared_per_face': float(observed[:, 5].mean()),
            'h_squared_standard_error': se(observed[:, 5]),
            'h_initial_decay_timescale_attempts_per_face': float(45*observed[:, 5].mean()/total_energy.mean()),
            'birth_minus_death': float(birth_death.mean()), 'birth_minus_death_standard_error': se(birth_death),
            'exact_one_function_upper_bound': base, 'charge_coefficient': c, 'fixed_hidden_coefficient': d,
            'estimated_trial_change': float(changes.mean()), 'trial_change_standard_error': se(changes),
            'estimated_trial_functional': base+float(changes.mean()),
            'extra_reduction_fraction_of_bare': -float(changes.mean())/bare,
            'scope': 'Independent exact stationary states; move-family importance sampling only estimates the Dirichlet functional. Empirical standard errors, not certified confidence bounds, no new physical update law or proved bulk mobility.'}


if __name__ == '__main__':
    print(json.dumps(measure()), flush=True)
