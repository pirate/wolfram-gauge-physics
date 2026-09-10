#!/usr/bin/env python3
"""Adjacency-assisted curvature transport, raw subgroup growth, angular diffusion.

New candidate primitives, explicitly using ambient fiber adjacency. No
Hamiltonian, spectral-energy filter, or continuum differential equation is
used to update links. The equal-rate scheduler remains a modeling choice.
"""
import itertools
import json
import math
import random
import statistics
from collections import deque
from pathlib import Path

import numpy as np

from cycle_charge_compatibility import cycle_group
from cycle_relational_dynamics import RelationalRule
from run_three_face import LinkOracle
from screen_braid_rules import nullspace
from triangle_reference import subgroup


class AngularFiber:
    def __init__(self, n):
        assert n >= 3 and n % 2
        self.n, self.m, self.g = n, (n-1)//2, cycle_group(n)
        g = self.g
        adjacency = [{(v-1) % n, (v+1) % n} for v in range(n)]
        distances, queue = {0: 0}, deque([0])
        while queue:
            v = queue.popleft()
            for w in adjacency[v]-distances.keys():
                distances[w] = distances[v]+1
                queue.append(w)
        self.steps = [h for h, p in enumerate(g.elements) if all(p[v] in adjacency[v] for v in range(n))]
        assert len(self.steps) == 2 and g.inv[self.steps[0]] == self.steps[1]
        self.rotations = [h for h in range(g.n) if h != g.identity and g.inv[h] != h]
        self.reflections = [h for h in range(g.n) if h != g.identity and g.inv[h] == h]
        self.powers = {}
        for u in self.steps:
            powers = [g.identity]
            for _ in range(1, n):
                powers.append(g.mul[u][powers[-1]])
            assert len(set(powers)) == n
            self.powers[u] = powers
        self.distance, self.unit = {}, {}
        for z in self.rotations:
            d = distances[g.elements[z][0]]
            units = [u for u in self.steps if self.powers[u][d] == z]
            assert len(units) == 1
            self.distance[z], self.unit[z] = d, units[0]
        self.maps = []
        for j in range(1, self.m+1):
            mapping = list(range(g.n))
            for z in self.rotations:
                d = self.distance[z]
                if j < self.m and d in (j, j+1):
                    mapping[z] = self.powers[self.unit[z]][2*j+1-d]
                elif j == self.m and d == self.m:
                    mapping[z] = g.inv[z]
            assert all(mapping[mapping[h]] == h for h in range(g.n))
            assert all(g.mul[g.inv[z]][mapping[z]] in self.steps for z in self.rotations if mapping[z] != z)
            self.maps.append(mapping)
        self.rules = [AngularRule(self, j) for j in range(self.m)]


class AngularRule:
    def __init__(self, fiber, index):
        self.fiber, self.index = fiber, index

    def apply_pair(self, w):
        f, g = self.fiber, self.fiber.g
        a, b = w
        product = g.mul[a][b]
        if a in f.reflections and b in f.distance:
            z = f.maps[self.index][b]
            return g.mul[product][g.inv[z]], z
        if a in f.distance and b in f.reflections:
            z = f.maps[self.index][a]
            return z, g.mul[g.inv[z]][product]
        return w

    def apply(self, w):
        return (*self.apply_pair(w[:2]), w[2])

    def __getitem__(self, code):
        n = self.fiber.g.n
        a, b, c = self.apply((code//(n*n), (code//n) % n, code % n))
        return (a*n+b)*n+c


def algebra(f):
    g = f.g
    labels = sorted(set(g.sectors)-{g.identity})
    rows, counts = set(), []
    reverse = lambda w: tuple(g.inv[h] for h in w[::-1])
    for rule in f.rules:
        moved = 0
        for w in itertools.product(range(g.n), repeat=2):
            target = rule.apply_pair(w)
            assert rule.apply_pair(target) == w
            assert g.mul[w[0]][w[1]] == g.mul[target[0]][target[1]]
            assert rule.apply_pair(reverse(w)) == reverse(target)
            if target != w:
                moved += 1
                assert all(rule.apply_pair(tuple(row[h] for h in w)) == tuple(row[h] for h in target) for row in g.conj)
                rows.add(tuple(sum(g.sectors[h] == s for h in w)-sum(g.sectors[h] == s for h in target) for s in labels))
        counts.append(moved)
    basis = nullspace(sorted(rows), len(labels))
    assert len(basis) == 2
    # Add the previously reconstructed RRR <-> ERZ relation only after
    # defining the adjacency rules; no charge enters their construction.
    r, s = f.reflections[:2]
    source, target = (r, r, s), (g.identity, r, g.mul[r][s])
    rows.add(tuple(sum(g.sectors[h] == z for h in source)-sum(g.sectors[h] == z for h in target) for z in labels))
    combined = nullspace(sorted(rows), len(labels))
    assert combined == [[1 if h in f.reflections else 2 for h in labels]]
    return {'cycle_vertices': f.n, 'adjacency_step_permutations': [g.elements[u] for u in f.steps],
            'moved_ordered_pairs_by_channel': counts, 'angular_class_charge_basis': basis,
            'combined_with_reaction_charge_basis': combined, 'nonidentity_class_representatives': labels}


def seed_pair(f, rotation):
    g, oracle, patch = f.g, LinkOracle(4, f.g), 0
    # Prepare, not evolve, a connection with fixed boundary reflection and
    # desired based pair. The third based face is identity throughout.
    p = f.reflections[0]
    paths = oracle.patches[patch][::-1]
    other = {e for path in paths[1:] for e, _ in path}
    rim = next(e for e, _ in paths[0] if e not in other)
    raw = [g.identity]*len(oracle.edges)
    raw[rim] = p
    assert tuple(oracle.transport(raw, path) for path in paths) == (p, g.identity, g.identity)
    target = (g.mul[p][g.inv[rotation]], rotation, g.identity)
    code = lambda w: (w[0]*g.n+w[1])*g.n+w[2]
    oracle.update(raw, patch, {code((p, g.identity, g.identity)): code(target)})
    return oracle, raw, p


def holonomy_image(oracle, raw):
    g = oracle.group
    adjacent = {}
    for edge, (a, b) in enumerate(oracle.edges):
        adjacent.setdefault(a, []).append((b, edge, False))
        adjacent.setdefault(b, []).append((a, edge, True))
    transports, queue = {0: g.identity}, [0]
    for a in queue:
        for b, edge, inverse in adjacent[a]:
            if b not in transports:
                transports[b] = g.mul[g.inv[raw[edge]] if inverse else raw[edge]][transports[a]]
                queue.append(b)
    based = [g.mul[g.inv[transports[b]]][g.mul[raw[e]][transports[a]]] for e, (a, b) in enumerate(oracle.edges)]
    return subgroup(g, based)


def growth_witness(f):
    g, k = f.g, f.n//3
    z = f.powers[f.steps[0]][k]
    oracle, raw, p = seed_pair(f, z)
    before = raw[:]
    h_before = holonomy_image(oracle, raw)
    q = lambda h: 0 if h == g.identity else 1 if h in f.reflections else 2
    charge = lambda values: sum(q(oracle.transport(values, path)) for path in oracle.face_paths)
    before_q = charge(raw)
    rule = f.rules[k-2]  # The edge between magnitudes k-1 and k.
    oracle.update(raw, 0, rule)
    h_after = holonomy_image(oracle, raw)
    assert len(h_before) == 6 and len(h_after) == 2*f.n
    assert charge(raw) == before_q
    changed = [e for e in range(len(raw)) if raw[e] != before[e]]
    assert len(changed) == 1 and set(changed) <= oracle.writes[0]
    reverse = raw[:]
    oracle.update(reverse, 0, rule)
    assert reverse == before
    return {'cycle_vertices': f.n, 'initial_rotation_step_count': k, 'final_rotation_step_count': k-1,
            'initial_holonomy_image_order': len(h_before), 'final_holonomy_image_order': len(h_after),
            'whole_mesh_charge': before_q, 'changed_edge': changed[0], 'boundary_reflection': p,
            'initial_links': before, 'final_links': raw, 'raw_inverse_exact': True}


def path_spectrum(f):
    rotations = f.powers[f.steps[0]][1:]
    ids = {z: i for i, z in enumerate(rotations)}
    N = len(rotations)
    laplacian = np.zeros((N, N))
    for mapping in f.maps:
        for i, z in enumerate(rotations):
            j = ids[mapping[z]]
            if i != j:
                laplacian[i, i] += 1
                laplacian[i, j] -= 1
    expected = np.diag([1]+[2]*(N-2)+[1])-np.eye(N, k=1)-np.eye(N, k=-1)
    assert np.array_equal(laplacian, expected)
    eigenvalues = [2-2*math.cos(math.pi*j/N) for j in range(N)]
    assert np.max(np.abs(np.linalg.eigvalsh(laplacian)-eigenvalues)) < 1e-12
    h = 2*math.pi/f.n
    return {'cycle_vertices': f.n, 'oriented_rotation_states': N,
            'generator_exactly_path_laplacian': True,
            'lowest_frame_odd_decay_rate': eigenvalues[1], 'lowest_gauge_even_decay_rate': eigenvalues[2],
            'scaled_frame_odd_rate': eigenvalues[1]/h**2, 'scaled_gauge_even_rate': eigenvalues[2]/h**2}


def heat_relaxation(f, trials=128):
    n, m, g = f.n, f.m, f.g
    h, k0 = 2*math.pi/n, n//3
    checkpoints = sorted(set(round(t*m/h**2) for t in (0, 0.25, 1, 2)))
    rotations = f.powers[f.steps[0]]
    exponent = {z: k for k, z in enumerate(rotations)}
    oracle, initial, _ = seed_pair(f, rotations[k0])
    paths = oracle.patches[0][::-1]
    samples = {tick: [] for tick in checkpoints}
    first_changes, examples = [], []
    for trial in range(trials):
        rng = random.Random(252220000+1000*n+trial)
        raw, first = initial[:], None
        for tick in range(checkpoints[-1]+1):
            word = tuple(oracle.transport(raw, p) for p in paths)
            k = exponent[word[1]]
            if tick in samples:
                samples[tick].append([math.cos(math.pi*j*(k-0.5)/(n-1)) for j in (1, 2)])
            if tick == checkpoints[-1]:
                break
            rule = f.rules[rng.randrange(m)]
            target = rule.apply(word)
            if target != word:
                oracle.update(raw, 0, rule)
                if first is None:
                    first = tick+1
                    # For n=3^a, either neighbor of n/3 is coprime to n.
                    assert len(holonomy_image(oracle, raw)) == 2*n
        first_changes.append(first)
        if trial == 0:
            examples.append({'initial_links': initial, 'final_links': raw})
    rows = []
    for tick in checkpoints:
        modes = []
        for i, j in enumerate((1, 2)):
            values = [sample[i] for sample in samples[tick]]
            eigenvalue = 2-2*math.cos(math.pi*j/(n-1))
            prediction = math.cos(math.pi*j*(k0-0.5)/(n-1))*(1-eigenvalue/m)**tick
            modes.append({'mode': 'boundary_frame_odd' if j == 1 else 'gauge_even',
                          'mean': statistics.mean(values), 'standard_error': statistics.stdev(values)/trials**0.5,
                          'exact_discrete_clock_mean': prediction})
        rows.append({'proposals': tick, 'scaled_time': tick*h*h/m, 'modes': modes})
    observed = [t for t in first_changes if t is not None]
    return {'cycle_vertices': n, 'independent_trajectories': trials, 'proposals_per_trajectory': checkpoints[-1],
            'scheduler': 'Uniform proposal among m involutions on one fixed shared pair; idle proposals count.',
            'initial_condition': 'Embedded C3 holonomy with rotation exponent n/3, fixed exterior links.',
            'mean_first_subgroup_growth_proposal': statistics.mean(observed),
            'standard_error_first_growth': statistics.stdev(observed)/len(observed)**0.5,
            'expected_first_growth_proposal': m/2, 'right_censored_trajectories': first_changes.count(None),
            'checkpoints': rows, 'raw_example': examples}


class AngularLastPair:
    def __init__(self, rule):
        self.rule = rule

    def apply(self, w):
        return (w[0], *self.rule.apply_pair(w[1:]))

    def __getitem__(self, code):
        n = self.rule.fiber.g.n
        a, b, c = self.apply((code//(n*n), (code//n) % n, code % n))
        return (a*n+b)*n+c


def reactive_activation(f, trials=128):
    g, n, m = f.g, f.n, f.m
    z0 = f.powers[f.steps[0]][n//3]
    oracle, raw, p = seed_pair(f, z0)
    r0 = g.mul[z0][p]
    source = (r0, r0, p)
    code = lambda w: (w[0]*g.n+w[1])*g.n+w[2]
    paths = oracle.patches[0][::-1]
    current = tuple(oracle.transport(raw, path) for path in paths)
    oracle.update(raw, 0, {code(current): code(source)})  # Initial-state preparation only.
    initial = raw[:]
    assert len(holonomy_image(oracle, initial)) == 6
    reaction = RelationalRule(g, ('E', 'R', 'Z'), 1)
    rules = [reaction]+[AngularLastPair(rule) for rule in f.rules]

    # Recover the entire invariant component from the actual primitive maps.
    rotations = f.powers[f.steps[0]][1:]
    states = [(g.mul[z][p], g.mul[z][p], p) for z in rotations]+[(g.identity, g.mul[z][p], z) for z in rotations]
    ids, N = {w: i for i, w in enumerate(states)}, len(rotations)
    laplacian = np.zeros((2*N, 2*N))
    for i, w in enumerate(states):
        for rule in rules:
            j = ids[rule.apply(w)]
            if i != j:
                laplacian[i, i] += 1
                laplacian[i, j] -= 1
    path = np.diag([1]+[2]*(N-2)+[1])-np.eye(N, k=1)-np.eye(N, k=-1)
    expected = np.block([[np.eye(N), -np.eye(N)], [-np.eye(N), np.eye(N)+path]])
    assert np.array_equal(laplacian, expected)
    eigenvalues = []
    for j in range(N):
        lam = 2-2*math.cos(math.pi*j/N)
        eigenvalues.extend([(lam+2-math.sqrt(lam*lam+4))/2, (lam+2+math.sqrt(lam*lam+4))/2])
    assert np.max(np.abs(np.linalg.eigvalsh(laplacian)-sorted(eigenvalues))) < 1e-12
    waiting, example = [], None
    for trial in range(trials):
        rng = random.Random(252221000+n*1000+trial)
        raw, event_history = initial[:], []
        for tick in range(1, 100001):
            channel = rng.randrange(m+1)
            rule = rules[channel]
            word = tuple(oracle.transport(raw, path) for path in paths)
            target = rule.apply(word)
            if target == word:
                continue
            oracle.update(raw, 0, rule)
            event_history.append({'proposal': tick, 'channel': channel, 'before_tuple': word, 'after_tuple': target})
            if channel != 0:
                assert len(holonomy_image(oracle, raw)) == 2*n
                waiting.append(tick)
                if example is None:
                    example = {'initial_links': initial, 'final_links': raw, 'changing_events': event_history}
                break
        else:
            raise ValueError('Activation observation ended without a first angular event')
    lam = 2-2*math.cos(2*math.pi/N)
    slow = (lam+2-math.sqrt(lam*lam+4))/2
    return {'cycle_vertices': n, 'initial_condition': 'Three reflections in the embedded C3 subgroup; no initial rotation face in the active triple.',
            'scheduler': 'Uniformly choose the existing rule 27 or one of m angular gates on the last two loops.',
            'trials': trials, 'first_full_subgroup_proposal_mean': statistics.mean(waiting),
            'first_full_subgroup_proposal_standard_error': statistics.stdev(waiting)/trials**0.5,
            'exact_mean_first_full_subgroup_proposal': 2*(m+1), 'first_growth_proposals': waiting,
            'component_states': 2*N, 'actual_generator_is_comb_graph': True,
            'slow_gauge_even_rate_per_gate_clock': slow, 'scaled_slow_gauge_even_rate': slow/(2*math.pi/n)**2,
            'raw_activation_example': example}


if __name__ == '__main__':
    fibers = {n: AngularFiber(n) for n in (9, 27, 81)}
    report = {'scope': 'New adjacency-assisted candidate on fixed cycle fibers and supplied base mesh; no derived physical energy, amplitudes, or evolving spacetime.',
              'exact_pair_algebra': [algebra(fibers[n]) for n in (9, 27)],
              'raw_subgroup_growth': [growth_witness(f) for f in fibers.values()],
              'derived_path_spectra': [path_spectrum(f) for f in fibers.values()],
              'autonomous_reaction_then_angular_activation': [reactive_activation(f) for f in fibers.values()],
              'nonstationary_shared_link_relaxation': []}
    for f in fibers.values():
        row = heat_relaxation(f)
        report['nonstationary_shared_link_relaxation'].append(row)
        print(json.dumps({'cycle': f.n, 'growth_mean': row['mean_first_subgroup_growth_proposal'],
                          'growth_expected': row['expected_first_growth_proposal'],
                          'last_checkpoint': row['checkpoints'][-1]}), flush=True)
    output = Path('data/fiber-angular-transport.json')
    output.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({'output': str(output)}), flush=True)
