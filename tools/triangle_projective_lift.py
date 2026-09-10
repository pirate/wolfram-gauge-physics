#!/usr/bin/env python3
"""Auxiliary spin-lift diagnostic for the derived boundary-odd action.

No spinor variable, phase, Hamiltonian, or change of dynamics is added.
Clifford signs are computed separately from actual classical permutations.
"""
import json
from collections import Counter

import numpy as np

from triangle_local_control_group import calculation, compose, binary_rank
from triangle_elastic_scattering import inverse


def sign_exponent(v, w):
    """e_v e_w = (-1)^c(v,w) e_(v xor w), with e_i^2=+1."""
    return sum((w & ((1 << i)-1)).bit_count() for i in range(10) if v >> i & 1) % 2


def projective_calculation():
    local = calculation()
    reps, tau = local['representative_state_ids'], local['frame_partner']
    tables, quotient = local['framed_tables'], local['quotient_tables']
    base = {i: k for k, rep in enumerate(reps) for i in (rep, tau[rep])}
    records = []
    for channel, table in zip(local['channels'], tables):
        matrix = np.zeros((10, 10), dtype=np.int64)
        for i, rep in enumerate(reps):
            target = table[rep]
            matrix[base[target], i] = 1 if target == reps[base[target]] else -1
        assert np.array_equal(matrix.T@matrix, np.eye(10, dtype=np.int64))
        involution = np.array_equal(matrix@matrix, np.eye(10, dtype=np.int64))
        negative = (10-int(np.trace(matrix)))//2 if involution else None
        records.append({'channel': channel, 'odd_trace': int(np.trace(matrix)),
                        'negative_eigenvalues_if_involution': negative,
                        'spin_lift_square_sign_if_involution': (-1)**(negative//2) if involution else None})
    even = [v for v in range(1 << 10) if v.bit_count() % 2 == 0]
    basis = [(1 << i) | (1 << 9) for i in range(9)]
    gram_rows = [sum(((v & w).bit_count() % 2) << j for j, w in enumerate(basis)) for v in basis]
    radical = [v for v in even if all((v & w).bit_count() % 2 == 0 for w in basis)]
    stabilizer = [v for v in even if not v & 1]
    stabilizer_radical = [v for v in stabilizer if all((v & w).bit_count() % 2 == 0 for w in stabilizer)]
    assert all(sign_exponent(v, v) == (v.bit_count()//2) % 2 for v in even)
    assert all((sign_exponent(v, w)+sign_exponent(w, v)) % 2 == (v & w).bit_count() % 2
               for v in even for w in even)
    # Construct actual primitive words for two commuting sign flips that fix
    # both raw states of pair zero, rather than only naming abstract elements.
    twist_word = local['central_twist_channel_indices']
    twist = tuple(range(20))
    for k in twist_word:
        twist = compose(tables[k], twist)
    start = sum(1 << i for i, rep in enumerate(reps) if twist[rep] == tau[rep])
    orbit, queue = {start: []}, [start]
    for old in queue:
        for k, permutation in enumerate(quotient):
            new = sum(1 << permutation[i] for i in range(10) if old >> i & 1)
            if new not in orbit:
                orbit[new] = orbit[old]+[k]; queue.append(new)
    def four_flip(mask):
        w = orbit[mask]
        return [tables.index(tuple(inverse(tables[k]))) for k in reversed(w)]+twist_word+w
    def pair_flip(a, b):
        filler = sum(1 << i for i in (4, 5, 6))
        word = four_flip(filler | (1 << a))+four_flip(filler | (1 << b))
        permutation = tuple(range(20))
        for k in word:
            permutation = compose(tables[k], permutation)
        mask = (1 << a) | (1 << b)
        assert all(permutation[i] == (tau[i] if mask >> base[i] & 1 else i) for i in range(20))
        return word, permutation
    aw, ap = pair_flip(1, 2)
    bw, bp = pair_flip(2, 3)
    assert compose(ap, bp) == compose(bp, ap)
    assert ap[reps[0]] == bp[reps[0]] == reps[0]
    v, w = (1 << 1) | (1 << 2), (1 << 2) | (1 << 3)
    commutator_sign = (-1)**((sign_exponent(v, w)+sign_exponent(w, v)) % 2)
    assert commutator_sign == -1
    return {'primitive_odd_actions': records,
            'even_flip_vectors': len(even), 'commutator_form_rank': binary_rank(gram_rows),
            'radical_masks': radical, 'volume_element_square': (-1)**sign_exponent(1023, 1023),
            'classical_state_stabilizer_flip_vectors': len(stabilizer),
            'stabilizer_radical_masks': stabilizer_radical,
            'stabilizer_quadratic_form_counts': dict(Counter(sign_exponent(v, v) for v in stabilizer)),
            'minimal_complex_dimension_for_this_cocycle': 2**(binary_rank(gram_rows)//2),
            'spin_pullback_group_order': 2*local['framed_group_order'],
            'fixed_classical_state_id': reps[0], 'fixed_classical_based_word': local['framed_states'][reps[0]],
            'commuting_classical_words': [aw, bw], 'word_lengths': [len(aw), len(bw)],
            'projective_commutator_sign': commutator_sign,
            'scope': 'One nontrivial auxiliary projective class of the derived odd-observable action. Not the full Schur multiplier, a physical spinor encoding, quantum interference, or an altered state-update rule.'}


if __name__ == '__main__':
    print(json.dumps(projective_calculation()), flush=True)
