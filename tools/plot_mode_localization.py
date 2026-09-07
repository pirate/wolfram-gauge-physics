#!/usr/bin/env python3
"""Plot measured spectral densities; coordinates are supplied mesh indices, not emergent space."""
import base64
import hashlib
import json
import zlib
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.linalg import eigh

from mode_localization import ModeGeometry, band_report, verify_positive_witness


def main():
    source = Path('data/d4-mode-localization.json')
    study = json.loads(source.read_text())
    transition = json.loads(Path('data/d4-spectral-transition.json').read_text())
    pairs = [r for r in study['static'] if r['condition'] == 'compact_noncommuting_pair']
    separated = []
    for r in pairs:
        if 'exact_count_certificate' not in r:
            continue
        g = ModeGeometry(r['side'])
        operator = g.operator(r['initial_links'])
        certificate = r['exact_count_certificate']
        verify_positive_witness(operator, certificate)
        vectors = np.frombuffer(zlib.decompress(base64.b64decode(certificate['witness_float64_le_zlib_base64'])), dtype='<f8').reshape(certificate['witness_shape'])
        values = np.asarray(r['above_flat_band']['eigenvalues'])
        if len(values) != 2 or values[1]-values[0] < 0.1:
            raise ValueError('individual density audit requires resolved nondegenerate eigenvalues')
        modes = [band_report(operator, values[j:j+1], vectors[:, j:j+1], g.core_distances(r['initial_links'])) for j in range(2)]
        separated.append({'side': r['side'], 'vertices': g.size,
                          'modes': [{k: v for k, v in m.items() if k != 'density'} for m in modes]})
    audit = {'source_sha256': hashlib.sha256(source.read_bytes()).hexdigest(),
             'resolved_modes': separated,
             'scope': 'Individual rank-one projector densities are valid here because both eigenvalues are separated; no mode identity across time is asserted.'}
    Path('data/d4-localization-density-audit.json').write_text(json.dumps(audit, separators=(',', ':'))+'\n')

    fig, axes = plt.subplots(2, 2, figsize=(13, 8), constrained_layout=True)
    ax = axes[0, 0]
    ax.loglog([r['vertices'] for r in pairs], [r['above_flat_band']['effective_vertices'] for r in pairs], 'o-', label='Two-mode band density')
    low = [r for r in pairs if r['ground'] is not None]
    ax.loglog([r['vertices'] for r in low], [r['ground']['effective_vertices'] for r in low], 's-', label='Lowest-mode density')
    ax.loglog([r['vertices'] for r in pairs], [r['vertices'] for r in pairs], ':', color='gray', label='Uniform density')
    ax.set(xlabel='Base vertices', ylabel=r'Effective vertices $1/\sum_v p_v^2$', title='Lowest mode spreads; above-band density stays compact')
    ax.legend(fontsize=8)
    ax.grid(alpha=.2)

    ax = axes[0, 1]
    r = pairs[-1]; side = r['side']; center = side//2; radius = 16
    p = np.asarray(r['above_flat_band']['density']).reshape(side, side)
    cropped = p[center-radius:center+radius+1, center-radius:center+radius+1]
    im = ax.imshow(np.log10(np.maximum(cropped, 1e-12)), origin='lower', extent=(-radius-.5, radius+.5, -radius-.5, radius+.5), cmap='magma', vmin=-8, vmax=-1)
    ax.plot([0, 1], [0, 0], color='cyan', linewidth=2)
    ax.plot([1, 2], [0, 0], color='lime', linewidth=2)
    ax.set(xlabel='Supplied mesh x index (relative to defect)', ylabel='Supplied mesh y index', title=f'192 × 192 torus: {cropped.sum():.1%} band mass in this crop')
    fig.colorbar(im, ax=ax, label=r'$\log_{10} p_v$')

    ax = axes[1, 0]
    for run in study['evolution']['runs']:
        combined = run['mode'] == 'combined'
        ax.plot([s['tick']/1000 for s in run['snapshots']], [s['above_flat_band']['rank'] for s in run['snapshots']],
                '-' if combined else '--', color='#148f89' if combined else '#888888', alpha=.6,
                marker='.' if combined else 'x', label=run['mode'] if run['trial'] == 0 else None)
    ax.set(xlabel='Rule-placement attempts (thousands; not physical time)', ylabel='Eigenvalues above flat band', title='Four matched trials: sampled count is not conserved', yticks=range(4))
    ax.legend(fontsize=8)
    ax.grid(alpha=.2)

    ax = axes[1, 1]
    g = ModeGeometry(transition['side'])
    for i, e in enumerate(transition['endpoints']):
        values = eigh(g.operator(e['links']).toarray(), eigvals_only=True)
        ax.scatter(np.full(8, i), values[-8:], color=('#4072ba', '#d05a49')[i], s=35)
    ax.axhline(9, color='black', linestyle=':', label='Flat upper band edge')
    ax.set(xlim=(-.4, 1.4), xticks=[0, 1], xticklabels=['Before', 'After'], ylabel='Eight largest eigenvalues',
           title='One edge update: exact above-band count 1 → 0')
    ax.legend(fontsize=8)
    fig.suptitle('Graph-native fiber modes — spectral diagnostics, not a physical wave law', fontsize=15)
    fig.savefig('docs/images/mode-localization.png', dpi=180)


if __name__ == '__main__':
    main()
