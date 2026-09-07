# Localized fiber modes: an infinite-lattice witness and a primitive disappearance

The square-fiber graph diagnostic now has **two certified above-band modes for a
specified two-link defect on the infinite triangular lattice**, within its derived
two-component sector. A compact integer
witness proves their existence; it is not an extrapolation from a finite box. The
same mode count is **not conserved by the actual primitive dynamics**: an independently
replayed one-edge update changes it from one to zero without changing the curvature-class
histogram.

This is a spectral result about the [graph-derived two-component sector](complete-loop-observer.md),
not an atomic simulation. The base lattice, square fiber, and initial defects are supplied.
We have not derived geometry, quantum amplitudes, a physical wave law, energy, or stable
matter. No potential or localization target enters the update rules. These are applications
of established spectral methods to this explicit finite-fiber model, not a claim to have
invented spectral localization or advanced continuum gauge theory in general.

![Measured size scaling, actual projector density, sampled evolution, and the exact single-event disappearance](images/mode-localization.png)

The heatmap uses supplied mesh indices, not an emergent spatial embedding. The bottom-left
lines connect samples; they do not track individual modes or establish behavior between samples.

## 1. A flat reference band derived from graph adjacency

Use the undirected steps $(1,0),(0,1),(1,1)$ on $\mathbb Z^2$ or its periodic quotient.
The horizontal connection Laplacian satisfies

$$\langle\psi,L_U\psi\rangle
=\sum_{u\to v}\|\psi_v-R_{U_{uv}}\psi_u\|^2\ge0.$$

For identity links its Fourier symbol, repeated for both fiber components, is

$$\lambda(k_x,k_y)=6-2\bigl(\cos k_x+\cos k_y+\cos(k_x+k_y)\bigr).$$

Since $9-\lambda=|1+e^{ik_x}+e^{-ik_y}|^2$, the infinite flat spectrum is $[0,9]$.
The same bound holds on the finite tori; the upper endpoint is attained when the side
is divisible by three. Thus eigenvalues above nine are outside the flat reference band,
**not low-energy binding levels**. The total bundle-graph restriction adds the previously
derived vertical shift $2I$; this does not alter mode densities.

The band qualification is essential. The **full** flat bundle graph also has the square
fiber's eigenvalue-zero and eigenvalue-four sectors, so its overall spectrum is $[0,13]$.
Our modes shift to approximately $11.026$ and $11.281$ in that full graph: they are
square-summable eigenstates **embedded in the full essential spectrum**, not above it.
The exact invariant two-component fiber sector makes the restricted spectral argument
valid. No stability under changes that mix fiber sectors has been proved.

## 2. An integer witness on 169 vertices proves the infinite result

Change just the links $(0,0)\to(1,0)$ and $(1,0)\to(2,0)$ to the two derived,
noncommuting reflections $r$ and $s$. Write $L_U=L_0+\Delta$.

Each changed edge contributes an off-diagonal symmetric block whose positive inertia
is $\operatorname{rank}(R_g-I)=1$. Consequently $n_+(\Delta)\le2$, and $L_0\le9I$
implies at most two eigenvalues of $L_U$ above nine, counting multiplicity.

For a lower bound, restrict trial vectors to the $13\times13$ square $[-6,6]^2$ and
extend them by zero outside. **Keep diagonal degree six on the boundary**: this is a
principal restriction of the infinite operator, not the Laplacian of the induced patch.
Numerical eigensolving proposes two trial columns, rounded to small integers. Their
authoritative check uses only exact integer arithmetic. For the saved $338\times2$
integer matrix $Q$, the checked matrices are

$$G=Q^T(L_U-9I)Q=
\begin{pmatrix}82247&-517\\-517&4701533\end{pmatrix},$$

$$M=Q^TQ=
\begin{pmatrix}16780151&-2544\\-2544&16775929\end{pmatrix}.$$

Both leading principal minors of $250G-M$ are strictly positive. Every nonzero
vector in this two-dimensional trial space therefore has Rayleigh quotient greater
than $9+1/250$. Finite-rank perturbations preserve essential spectrum; the variational
principle and the upper bound together prove **exactly two discrete eigenvalues above
nine, both strictly above $9.004$**, on the infinite lattice. The relevant general results
are the min-max principle and compact-perturbation invariance in
[Teschl, sections 4.4 and 6.4](https://www.mat.univie.ac.at/~gerald/ftp/book-schroe/schroe2.pdf).

The certificate stores all integer trial columns, the two matrices, and a separately
recomputable dyadic witness. Verification does not run an eigensolver. A regression checks
the independently constructed restriction against the corresponding torus principal block.
The existence/count proof does not depend on the numerical eigenvalue approximations.

There is also an elementary localization bound. For a normalized eigenvector with
$\lambda>9$, set $\eta=\Delta\psi$, supported on the three defect endpoints $S$.
The norm-convergent series

$$\psi=(\lambda-L_0)^{-1}\eta
=\frac1\lambda\sum_{n\ge0}(L_0/\lambda)^n\eta$$

and finite propagation of $L_0$ imply, for graph-distance balls $B_R(S)$,

$$\|\mathbf1_{B_R(S)^c}\psi\|
\le\frac{\|\eta\|}{\lambda-9}(9/\lambda)^{R+1}
\le\frac4{\lambda-9}(9/\lambda)^{R+1}.$$

Here $\|\Delta\|\le4$ by the two edge contributions. This is a conservative exponential
tail bound, not a fitted localization length and not physical propagation in time.

## 3. Measured densities, with degeneracies handled explicitly

For a complete spectral subspace of dimension $d$, use its normalized projector density

$$p_v=\frac1d\operatorname{tr}P_{vv},\qquad
\operatorname{IPR}=\sum_vp_v^2,\qquad V_{\rm eff}=1/\operatorname{IPR}.$$

This is invariant under local orthogonal frames and changes of basis inside the subspace.
It is a graph-mode participation measure, not a Born probability postulate. A compact
band average alone would not prove that every constituent mode is compact, so the two
well-separated eigenvalues are additionally audited as individual rank-one projectors.

The torus study covers sides $6,12,24,48,96,192$. At the largest two sizes:

| Side | Lower above-band eigenvalue | Its effective vertices | Upper eigenvalue | Its effective vertices |
| ---: | ---: | ---: | ---: | ---: |
| 96 | 9.0259839056 | 56.580974 | 9.2811659668 | 8.187171 |
| 192 | 9.0259839042 | 56.580923 | 9.2811659668 | 8.187171 |

These are numerical estimates, not certified decimal enclosures of the infinite eigenvalues.
The combined band has about 18.081 effective vertices; that number is **not the support
volume of either individual mode**. Meanwhile the lowest mode spreads approximately
uniformly over the growing graph in the complete decompositions through side 48.

Controls include flat links, a single reflection, and the two existing class-indistinguishable
conversion outcomes. At side 48 their above-band counts are respectively zero, one, two,
and one. The raw conversions are independently checked against the C++ engine. Full dense
spectra are used through side 48; the two larger pair cases use sparse eigenpairs plus
exact positive-subspace certificates saturating the independent rank upper bound. Failure
to saturate the bound is reported as unproven completeness, never as absence of modes.

## 4. Curvature permits these modes but does not protect them

For each oriented triangle let $S_f$ be its $6\times6$ block matrix with identity diagonal
blocks and transported adjacency off the diagonal. On this degree-six closed mesh,

$$9I-L_U=\tfrac12\sum_f\iota_f S_f\iota_f^T.$$

Every edge belongs to two faces and every vertex to six, giving the identity exactly.
After local frame changes, a face matrix depends only on its holonomy $H_f$. For each
fiber holonomy eigenphase $\theta$, its three eigenvalues are
$1+2\cos((\theta+2\pi j)/3)$, $j=0,1,2$. Thus its negative inertia is zero for identity,
one for either reflection class, and two for the quarter-turn and central-half-turn classes.
The code checks this using exact inertia for all 512 raw triangle connections.

Subadditivity of negative inertia gives the gauge-invariant bound

$$n_+(L_U-9I)\le N_r+N_s+2N_\rho+2N_z.$$

It is an **upper bound**, not a particle count. In particular, flat faces prohibit
above-band modes even with nontrivial global torus holonomy. Nonzero curvature does not
guarantee any such mode. The class-indistinguishable conversion pair also shows that
face labels cannot determine their exact number.

## 5. One primitive update removes the last above-band mode

Four independently seeded, matched transport-only/combined experiments run 200,000
rule-placement attempts each, including no-ops. Each of the eight raw histories is
independently replayed and reversed. Spectra are sampled every 20,000 attempts; these
samples alone do not establish uninterrupted persistence or individual mode identities.

In transport trial zero, bisection within the first sampled positive-to-zero bracket
finds one crossing at attempt **39,304**, rule 0, rooted patch 610. This search does not
claim the first crossing in the full history. Raw edge 305 is the only changed edge.
Both states contain the same four nonflat face classes: two of each reflection class.
The exact inertia of $L_U-9I$, written (positive, negative, zero), is

$$ (1,286,1)\quad\longrightarrow\quad(0,287,1). $$

Before the update the above-band eigenvalue is approximately $9.041604956$. Afterward
there are exactly none; there is still one eigenvalue exactly at the band edge. Exact
symmetric integer congruence elimination verifies these counts without tolerance decisions.
An independent FLINT audit computes the exact integer characteristic polynomials and counts
positive roots by sign variations; symmetry makes all roots real, so the count is exact.
The curvature upper bound remains four on both sides. This is a spectral rearrangement
under a reversible gauge update, **not particle annihilation or a failure of physical
energy conservation**: neither particle number nor physical energy has been established.

## Reproduce and challenge the result

```bash
cmake --build build -j 4
export OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
uv run --with numpy --with scipy python tests/mode_localization.py
uv run --with numpy --with scipy python tools/compact_mode_certificate.py \
  --verify data/d4-compact-mode-certificate.json
uv run --with numpy --with scipy python tools/mode_localization.py \
  --output out/mode-localization.json
uv run --with numpy --with scipy python tools/spectral_transition.py \
  --output out/spectral-transition.json
uv run --with numpy --with scipy --with python-flint python tools/spectral_transition.py \
  --verify-with-flint data/d4-spectral-transition.json
uv run --with numpy --with scipy --with matplotlib python tools/plot_mode_localization.py
```

The plot script reads checked-in measurements and regenerates the density audit and image.
Cross-platform numerical comparisons use tolerances and invariant densities, not eigenvector
signs or byte-identical floating-point JSON. Integer certificates are recomputed exactly.

Next: identify which microscopic observables predict persistence, test isolated defects
and encounter controls under the actual updates, and seek genuinely long-lived localized
structures without selecting update rules to preserve this diagnostic. A static bound mode
does not supply the missing physical evolution law or a route to molecules by itself.
