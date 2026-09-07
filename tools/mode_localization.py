#!/usr/bin/env python3
"""Gauge-invariant localization diagnostics, not a physical wave evolution law."""
import argparse
import base64
import hashlib
import json
import random
import zlib
from pathlib import Path

import numpy as np
import scipy
from scipy.linalg import eigh
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import eigsh

from bank_equilibrium import BankFactor, compiled_run
from d4_loop_observer import LoopForest
from face_energy_obstruction import derive_fiber_group
from fiber_mode_probe import representation
from run_shared_edge import mesh_geometry

BAND_EDGE = 9.0
TOLERANCE = 1e-8


class ModeGeometry:
    def __init__(self, side):
        if not 3 <= side <= 192:
            raise ValueError('invalid bounded mode geometry')
        self.side, self.size = side, side*side
        self.group = derive_fiber_group()
        self.matrices = np.asarray(representation(self.group)['matrices'], dtype=np.int64)
        self.edges, self.faces, _ = mesh_geometry(side)
        self.ids = {e: i for i, e in enumerate(self.edges)}
        self.paths = [[(self.ids[tuple(sorted((u, v)))], u > v) for u, v in zip(f, f[1:])]
                      for f in self.faces]
        self.adjacency = [set() for _ in range(self.size)]
        for u, v in self.edges:
            self.adjacency[u].add(v)
            self.adjacency[v].add(u)

    def seed(self, reflection_only=False):
        y = x = self.side//2
        v, w, t = [y*self.side+(x+i) % self.side for i in range(3)]
        links = [0]*len(self.edges)
        for u, w, value in [(v, w, 1)]+([] if reflection_only else [(w, t, 2)]):
            links[self.ids[tuple(sorted((u, w)))]] = value if u < w else self.group.inv[value]
        return links

    def operator(self, links):
        if len(links) != len(self.edges) or any(not isinstance(x, int) or not 0 <= x < 8 for x in links):
            raise ValueError('invalid mode link vector')
        rows, columns, values = [], [], []
        degrees = [0]*self.size
        for (u, v), x in zip(self.edges, links):
            degrees[u] += 1
            degrees[v] += 1
            for a in range(2):
                for b in range(2):
                    value = -int(self.matrices[x, a, b])
                    if value:
                        rows.extend((2*v+a, 2*u+b))
                        columns.extend((2*u+b, 2*v+a))
                        values.extend((value, value))
        for v, d in enumerate(degrees):
            for a in range(2):
                rows.append(2*v+a); columns.append(2*v+a); values.append(d)
        result = coo_matrix((values, (rows, columns)), shape=(2*self.size, 2*self.size), dtype=float).tocsr()
        if (result-result.T).nnz:
            raise ValueError('mode operator is not symmetric')
        return result

    def sectors(self, links):
        result = []
        for path in self.paths:
            h = 0
            for e, reverse in path:
                x = self.group.inv[links[e]] if reverse else links[e]
                h = self.group.mul[x][h]
            result.append(self.group.sectors[h])
        return result

    def core_distances(self, links):
        core = sorted({v for face, c in zip(self.faces, self.sectors(links)) if c for v in face[:-1]})
        distances, queue = [-1]*self.size, core[:]
        for v in core:
            distances[v] = 0
        for v in queue:
            for w in sorted(self.adjacency[v]):
                if distances[w] < 0:
                    distances[w] = distances[v]+1
                    queue.append(w)
        return distances

    def parallel_dimension(self, links):
        forest = LoopForest(range(self.size), self.edges, self.group)
        first = None
        for values in forest.based_loops(links):
            for x in values:
                for row in self.matrices[x]-np.eye(2, dtype=np.int64):
                    a, b = map(int, row)
                    if not (a or b):
                        continue
                    if first is None:
                        first = a, b
                    elif first[0]*b-first[1]*a:
                        return 0
        return 2 if first is None else 1


def density(vectors):
    if vectors.shape[0] % 2 or not vectors.shape[1]:
        raise ValueError('density requires a nonempty two-component subspace')
    return np.sum(np.abs(vectors.reshape(-1, 2, vectors.shape[1]))**2, axis=(1, 2))/vectors.shape[1]


def band_report(operator, values, vectors, distances):
    if not len(values):
        return {'rank': 0, 'eigenvalues': [], 'density': None, 'ipr': None, 'effective_vertices': None}
    residuals = np.linalg.norm(operator@vectors-vectors*values, axis=0)
    orthogonal_error = float(np.max(np.abs(vectors.T@vectors-np.eye(len(values)))))
    if float(max(residuals)) > 1e-8 or orthogonal_error > 1e-8:
        raise ValueError('mode residual or orthogonality validation failed')
    p = density(vectors)
    if abs(float(p.sum())-1) > 1e-9:
        raise ValueError('projector density is not normalized')
    ipr = float(p@p)
    masses = []
    for radius in (0, 1, 2, 4, 8):
        mask = np.asarray([0 <= d <= radius for d in distances])
        masses.append({'radius': radius, 'mass': float(p[mask].sum()), 'uniform_mass': float(mask.mean())})
    return {'rank': len(values), 'eigenvalues': values.tolist(), 'density': p.tolist(),
            'ipr': ipr, 'effective_vertices': 1/ipr,
            'max_eigenpair_residual': float(max(residuals)), 'orthogonality_error': orthogonal_error,
            'mass_near_curvature_vertices': masses}


def full_profile(geometry, links):
    operator = geometry.operator(links)
    values, vectors = eigh(operator.toarray(), driver='evr')
    if len(values) != operator.shape[0] or not np.isfinite(values).all() or values[0] < -TOLERANCE:
        raise ValueError('invalid complete spectrum')
    if (abs(float(values.sum()-operator.diagonal().sum())) > 1e-7 or
            abs(float(values@values-operator.multiply(operator).sum())) > 1e-6):
        raise ValueError('full spectrum trace checks failed')
    low = values <= values[0]+TOLERANCE
    high = values > BAND_EDGE+TOLERANCE
    distances = geometry.core_distances(links)
    nullity = geometry.parallel_dimension(links)
    if int(np.sum(np.abs(values) <= TOLERANCE)) != nullity:
        raise ValueError('numerical zero modes disagree with exact holonomy-fixed sections')
    return {'method': 'complete dense symmetric eigendecomposition; no multiplicity truncation',
            'side': geometry.side, 'vertices': geometry.size,
            'parallel_section_dimension': nullity, 'flat_reference_band': [0, BAND_EDGE],
            'spectral_tolerance': TOLERANCE, 'near_upper_band_edge_count': int(np.sum(np.abs(values-BAND_EDGE) <= TOLERANCE)),
            'ground_cluster_gap': float(values[np.sum(low)]-values[0]) if np.sum(low) < len(values) else None,
            'ground': band_report(operator, values[low], vectors[:, low], distances),
            'above_flat_band': band_report(operator, values[high], vectors[:, high], distances)}


def positive_rank_bound(geometry, links, background=None):
    """Above-9 count bound relative to a gauge-trivial flat background, not arbitrary curvature."""
    background = [0]*len(links) if background is None else background
    if len(background) != len(links) or len(links) != len(geometry.edges):
        raise ValueError('background link dimensions disagree')
    if any(not isinstance(x, int) or not 0 <= x < 8 for x in list(links)+list(background)):
        raise ValueError('invalid rank-bound links')
    # Flat faces alone do not exclude nontrivial torus cycles. Require trivial
    # holonomy on every fundamental loop, hence global equivalence to identity.
    if any(background):
        forest = LoopForest(range(geometry.size), geometry.edges, geometry.group)
        if any(x for component in forest.based_loops(background) for x in component):
            raise ValueError('rank-bound background must be gauge-trivial')
    total = 0
    for x, y in zip(links, background):
        a, b = (geometry.matrices[x]-geometry.matrices[y]).tolist()
        total += 2 if a[0]*b[1]-a[1]*b[0] else int(any(a+b))
    return total


def exact_positive_witness(operator, vectors):
    """Exact dyadic Rayleigh matrix for L-9I; positive minors certify a subspace above 9."""
    if (vectors.ndim != 2 or vectors.shape[0] != operator.shape[0] or
            not np.isfinite(vectors).all() or vectors.shape[1] not in (1, 2)):
        raise ValueError('bounded exact certificate supports one or two directions')
    ratios = [[float(x).as_integer_ratio() for x in vectors[:, j]] for j in range(vectors.shape[1])]
    power = max(den.bit_length()-1 for column in ratios for _, den in column)
    integer = [[num << (power-(den.bit_length()-1)) for num, den in column] for column in ratios]
    shifted = operator.copy()
    shifted.setdiag(shifted.diagonal()-9)
    shifted.eliminate_zeros()
    rows = shifted.tocsr()
    image = []
    for column in integer:
        target = []
        for i in range(rows.shape[0]):
            value = 0
            for position in range(rows.indptr[i], rows.indptr[i+1]):
                coefficient = float(rows.data[position])
                if not coefficient.is_integer():
                    raise ValueError('exact certificate requires integer graph operator')
                value += int(coefficient)*column[int(rows.indices[position])]
            target.append(value)
        image.append(target)
    gram = [[sum(a*b for a, b in zip(left, right)) for right in image] for left in integer]
    minors = [gram[0][0]]
    if len(gram) == 2:
        if gram[0][1] != gram[1][0]:
            raise ValueError('exact Rayleigh matrix is not symmetric')
        minors.append(gram[0][0]*gram[1][1]-gram[0][1]*gram[1][0])
    if any(x <= 0 for x in minors):
        raise ValueError('candidate subspace is not strictly above the flat band')
    raw = np.ascontiguousarray(vectors, dtype='<f8').tobytes()
    return {'directions': vectors.shape[1], 'dyadic_scale_power': power,
            'scaled_integer_rayleigh_matrix': [[str(x) for x in row] for row in gram],
            'leading_principal_minors_positive': [True]*len(minors),
            'witness_shape': list(vectors.shape), 'witness_float64_le_zlib_base64': base64.b64encode(zlib.compress(raw)).decode(),
            'witness_sha256': hashlib.sha256(raw).hexdigest()}


def verify_positive_witness(operator, certificate):
    """Recompute the integer proof from saved dyadic vectors, without an eigensolver."""
    shape = certificate['witness_shape']
    if shape not in ([operator.shape[0], 1], [operator.shape[0], 2]):
        raise ValueError('invalid saved witness dimensions')
    raw = zlib.decompress(base64.b64decode(certificate['witness_float64_le_zlib_base64'], validate=True))
    if len(raw) != 8*shape[0]*shape[1] or hashlib.sha256(raw).hexdigest() != certificate['witness_sha256']:
        raise ValueError('saved witness bytes or hash disagree')
    vectors = np.frombuffer(raw, dtype='<f8').reshape(shape)
    actual = exact_positive_witness(operator, vectors)
    # Compression may vary across zlib versions; the raw bytes and math may not.
    for key in actual.keys()-{'witness_float64_le_zlib_base64'}:
        if actual[key] != certificate[key]:
            raise ValueError('saved witness arithmetic disagrees: '+key)
    return actual['directions']


def sparse_certified_profile(geometry, links):
    operator = geometry.operator(links)
    bound = positive_rank_bound(geometry, links)
    if bound not in (1, 2):
        raise ValueError('large sparse audit requires a one- or two-direction bound')
    rng = np.random.default_rng(620917)
    values, vectors = eigsh(operator, k=6, which='LA', v0=rng.normal(size=operator.shape[0]), tol=1e-12)
    order = np.argsort(values); values, vectors = values[order], vectors[:, order]
    selected = values > BAND_EDGE+TOLERANCE
    if int(np.sum(selected)) != bound:
        raise ValueError('sparse modes do not saturate the independently derived bound; completeness unproven')
    certificate = exact_positive_witness(operator, vectors[:, selected])
    return {'method': 'sparse eigenpairs with exact positive-subspace witness saturating a rank-perturbation upper bound',
            'side': geometry.side, 'vertices': geometry.size, 'above_band_rank_upper_bound': bound,
            'flat_reference_band': [0, BAND_EDGE], 'spectral_tolerance': TOLERANCE,
            'ground': None, 'exact_count_certificate': certificate,
            'above_flat_band': band_report(operator, values[selected], vectors[:, selected], geometry.core_distances(links))}


def converted_seeds(side):
    geometry = ModeGeometry(side)
    initial = geometry.seed()
    factor = BankFactor(side, json.loads(Path('data/d4-triple-channels.json').read_text()))
    patch = next(p for p in range(len(factor.fans)) if factor.oracle.code(initial, 1, p) == 66)
    result = {'flat': [0]*len(initial), 'single_reflection': geometry.seed(True), 'compact_noncommuting_pair': initial}
    for rule, label in ((33, 'converted_equal'), (34, 'converted_opposite')):
        links = initial[:]
        factor.oracle.update(links, rule, patch, factor.tables[rule])
        run = compiled_run(factor, initial, [rule*len(factor.fans)+patch], 1)[1]
        if links != run['final_links']:
            raise ValueError('scaled seed conversion differs from compiled raw-link replay')
        result[label] = links
    if geometry.sectors(result['converted_equal']) != geometry.sectors(result['converted_opposite']):
        raise ValueError('scaled hidden-bit witness has different face classes')
    return geometry, result


def static_study(sides, large_sides):
    records = []
    for side in sides:
        geometry, inputs = converted_seeds(side)
        for name, links in inputs.items():
            result = full_profile(geometry, links)
            result.update(condition=name, initial_links=links)
            records.append(result)
            print('Static', side, name, 'above-band', result['above_flat_band']['rank'],
                  'IPR', result['above_flat_band']['ipr'], flush=True)
    for side in large_sides:
        geometry = ModeGeometry(side)
        links = geometry.seed()
        result = sparse_certified_profile(geometry, links)
        result.update(condition='compact_noncommuting_pair', initial_links=links)
        records.append(result)
        print('Certified sparse', side, result['above_flat_band']['eigenvalues'], 'IPR', result['above_flat_band']['ipr'], flush=True)
    return records


def evolution(side, attempts, stride, trials, seed):
    if (not 3 <= side <= 48 or not 1 <= attempts <= 1000000 or
            stride < 1 or attempts % stride or not 1 <= trials <= 8):
        raise ValueError('invalid bounded evolution measurement clock')
    geometry = ModeGeometry(side)
    factor = BankFactor(side, json.loads(Path('data/d4-triple-channels.json').read_text()))
    initial = geometry.seed()
    output = []
    for trial in range(trials):
        rng = random.Random(seed+100003*trial)
        schedule = [rng.randrange(49*len(factor.fans)) for _ in range(attempts)]
        for run in compiled_run(factor, initial, schedule, stride):
            links, cursor, snapshots = initial[:], 0, []
            for tick in range(0, attempts+1, stride):
                while cursor < len(run['events']) and run['events'][cursor][0] <= tick:
                    _, rule, patch, code, target = run['events'][cursor]
                    if factor.oracle.update(links, rule, patch, factor.tables[rule]) != (code, target):
                        raise ValueError('spectral snapshot replay differs from the verified event')
                    cursor += 1
                profile = full_profile(geometry, links)
                profile.update(tick=tick, links=links[:], nonflat_faces=sum(bool(c) for c in geometry.sectors(links)))
                current = profile['above_flat_band']['density']
                previous = snapshots[-1]['above_flat_band']['density'] if snapshots else None
                profile['previous_density_affinity'] = float(np.sqrt(np.asarray(current)*np.asarray(previous)).sum()) if current is not None and previous is not None else None
                snapshots.append(profile)
            if links != run['final_links']:
                raise ValueError('mode snapshots missed final raw state')
            output.append({'trial': trial, 'mode': run['mode'], 'seed': seed+100003*trial,
                           'initial_links': initial, 'events': run['events'], 'snapshots': snapshots,
                           'conversion_events': run['conversion_events'], 'independent_link_replay': run['independent_link_replay'],
                           'exact_link_inverse': run['exact_link_inverse']})
            print('Evolution', trial, run['mode'], 'conversions', run['conversion_events'],
                  'band ranks', [s['above_flat_band']['rank'] for s in snapshots], flush=True)
    return {'side': side, 'attempts': attempts, 'stride': stride, 'trials': trials, 'seed': seed,
            'clock': 'iid uniform rule-placement attempts, including no-ops; not physical time',
            'affinity': 'sum_v sqrt(p_v(t) p_v(previous sample)); a density similarity, not quantum fidelity or mode identity',
            'runs': output}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--sides', type=int, nargs='+', default=[6, 12, 24, 48])
    parser.add_argument('--large-sides', type=int, nargs='*', default=[96, 192])
    parser.add_argument('--evolution-side', type=int, default=12)
    parser.add_argument('--attempts', type=int, default=200000)
    parser.add_argument('--stride', type=int, default=20000)
    parser.add_argument('--trials', type=int, default=4)
    parser.add_argument('--seed', type=int, default=5820391)
    parser.add_argument('--output', type=Path, default=Path('out/mode-localization.json'))
    args = parser.parse_args()
    if any(not 3 <= s <= 48 for s in args.sides) or any(not 3 <= s <= 192 for s in args.large_sides):
        parser.error('invalid bounded study sizes')
    if (not 0 <= args.trials <= 8 or not 3 <= args.evolution_side <= 48 or
            not 1 <= args.attempts <= 1000000 or args.stride < 1 or args.attempts % args.stride):
        parser.error('invalid bounded evolution measurement clock')
    result = {'schema': 1, 'numpy': np.__version__, 'scipy': scipy.__version__,
              'static': static_study(args.sides, args.large_sides),
              'scope': 'Graph-native two-component spectral diagnostics of classical gauge-link states. No physical Hamiltonian, wave evolution, particle identity, or binding energy is asserted.'}
    if args.trials:
        result['evolution'] = evolution(args.evolution_side, args.attempts, args.stride, args.trials, args.seed)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, separators=(',', ':'))+'\n')


if __name__ == '__main__':
    main()
