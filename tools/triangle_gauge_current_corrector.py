#!/usr/bin/env python3
"""An explicitly gauge-hidden corrector E L b for winding transport.

E is the positive elastic generator. It annihilates every function of the
complete charge field, so E L b is conditionally centered without guessing
an average activity or quotienting unrelated local frames.
"""
import itertools
import json
import random
from collections import Counter
import numpy as np

from triangle_reference import AxialReference
from triangle_winding_corrector import census, difference
from pathlib import Path


def value(poly, q):
    return sum(coefficient for assignment, coefficient in poly.items()
               if all(q[f] == a for f, a in assignment))


class GaugeCurrent:
    def __init__(self, side=3, calculate_moments=True):
        self.x, self.poly, self.variance_counts, self.twice_energy_counts = census(side, calculate_moments)
        self.kernels = []
        for faces in self.x.p.supports[1]:
            kernel = Counter()
            for target in itertools.permutations((0, 1, 2)):
                for assignment, coefficient in difference(self.poly, faces, (1, 1, 1), target).items():
                    kernel[assignment] -= 2*coefficient
            self.kernels.append({a: c for a, c in kernel.items() if c})
        terms = []
        for fan, kernel in enumerate(self.kernels):
            for assignment, coefficient in kernel.items():
                assert len(assignment) <= 2
                padded = list(assignment)+[(-1, -1)]*(2-len(assignment))
                terms.append((fan, *padded[0], *padded[1], coefficient))
        self.kernel_terms = np.array(terms, dtype=np.int64).T
        oracle = self.x.e.factor.oracle
        self.affected_fans = [
            [g for g, reads in enumerate(oracle.fan.reads) if pair[0][0][0] in reads]
            for pair in oracle.pairs]

    def weights(self, q):
        fan, f, a, g, b, coefficient = self.kernel_terms
        q = np.asarray(q)
        selected = ((f < 0) | (q[f] == a)) & ((g < 0) | (q[g] == b))
        weights = np.zeros(self.x.supports, dtype=np.int64)
        np.add.at(weights, fan, coefficient*selected)
        weights *= np.all(q[np.array(self.x.p.supports[1])] == 1, axis=1)
        return weights.tolist()

    def activities(self, raw):
        result = []
        for patch in range(self.x.supports):
            code = self.x.e.factor.oracle.code(raw, 1, patch)
            a, b, c = code//36, code//6 % 6, code % 6
            result.append(int((a == b) != (b == c)))
        return result

    def elastic_moments(self, raw):
        weights = self.weights(self.x.p.charge(raw))
        initial = self.activities(raw)
        total = squares = 0
        for rule in (1, 2):
            for patch in range(self.x.supports):
                moved = raw[:]
                a, b = self.x.apply(moved, rule*self.x.supports+patch)
                if a != b:
                    delta = 0
                    for fan in self.affected_fans[patch]:
                        if not weights[fan]:
                            continue
                        code = self.x.e.factor.oracle.code(moved, 1, fan)
                        a, b, c = code//36, code//6 % 6, code % 6
                        delta += weights[fan]*(initial[fan]-int((a == b) != (b == c)))
                    total += delta
                    squares += delta*delta
        assert squares % 2 == 0
        return total, squares//2

    def hidden_corrector(self, raw):
        return self.elastic_moments(raw)[0]

    def locality_bound(self):
        oracle = self.x.e.factor.oracle
        edge_patches = {}
        for patch, pair in enumerate(oracle.pairs):
            edge_patches.setdefault(pair[0][0][0], []).append(patch)
        reads = []
        for edge, patches in edge_patches.items():
            assert len(patches) == 2
            support = {e for patch in patches for path in oracle.pairs[patch] for e, _ in path}
            fans = set().union(*(set(self.affected_fans[patch]) for patch in patches))
            for fan in fans:
                support.update(oracle.fan.reads[fan])
                faces = set(self.x.p.supports[1][fan])
                faces.update(f for assignment in self.kernels[fan] for f, _ in assignment)
                support.update(e for face in faces for e, _ in oracle.fan.face_paths[face])
            reads.append(support)
        writes = [({pair[0][0][0]}, 3) for pair in oracle.pairs]
        writes += [(write, 12) for write in oracle.fan.writes]
        c1 = max(sum(bool(write & support) for support in reads) for write, _ in writes)
        c2 = max(sum(weight for write, weight in writes if write & support) for support in reads)
        return {'side': self.x.geometry.side, 'faces': self.x.geometry.faces,
                'maximum_local_corrector_read_edges': max(map(len, reads)),
                'maximum_terms_changed_by_one_slot': c1,
                'maximum_slots_touching_one_term': c2,
                'T_over_S_upper_bound': 12*c1*c2,
                'scope': 'Conservative read-support bound, not the measured T/S ratio.'}

    def full_generator_current_drift(self, raw):
        before = value(self.poly, self.x.p.charge(raw))
        total = 0
        for op in range(self.x.operator_count):
            moved = raw[:]
            a, b = self.x.apply(moved, op)
            if a != b:
                total += before-value(self.poly, self.x.p.charge(moved))
        return total

    def witness(self):
        bank = json.loads(Path('data/triangle-feedback.json').read_text())['bank']
        reference = AxialReference(self.x.geometry.side, bank, self.x.geometry.faces)
        rng = random.Random(114837001)
        for trial in range(64):
            raw = reference.sample(rng, 'nonabelian_reflections')['links']
            q = self.x.p.charge(raw)
            weights, activity = self.weights(q), self.activities(raw)
            for patch in range(self.x.supports):
                moved = raw[:]
                old, new = self.x.apply(moved, self.x.supports+patch)
                if old == new:
                    continue
                after = self.activities(moved)
                delta = sum(k*(b-a) for k, a, b in zip(weights, activity, after))
                if not delta:
                    continue
                assert self.x.p.charge(moved) == q
                u = self.full_generator_current_drift(raw)
                v = self.full_generator_current_drift(moved)
                assert v-u == delta
                h_before, h_after = self.hidden_corrector(raw), self.hidden_corrector(moved)
                assert h_before or h_after
                return {'side': self.x.geometry.side, 'faces': self.x.geometry.faces,
                        'reference_seed': 114837001, 'sample_index': trial,
                        'elastic_patch': patch, 'elastic_rule': 1,
                        'initial_links': raw, 'final_links': moved,
                        'complete_charge_field': q, 'current_drift_b3': value(self.poly, q),
                        'L_b3_before': u, 'L_b3_after': v, 'difference_L_b3': delta,
                        'E_L_b3_before': h_before, 'E_L_b3_after': h_after,
                        'changed_weighted_activities': [
                            {'fan': p, 'weight': k, 'before': a, 'after': b,
                             'contribution_to_difference': k*(b-a)}
                            for p, (k, a, b) in enumerate(zip(weights, activity, after)) if k and a != b],
                        'scope': 'Exact current-specific gauge-hidden witness; not a numerical value of the two-function mobility bound.'}
        raise RuntimeError('No current-specific witness in the bounded reference search')


if __name__ == '__main__':
    print(json.dumps(GaugeCurrent().witness()), flush=True)
