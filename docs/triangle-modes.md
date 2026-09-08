# Full-graph localization and an orientation-resolved persistence obstruction

The specified noncommuting two-link defect has **exactly two isolated, square-summable
modes above the full infinite bundle graph's flat band**, both strictly above $12.05$.
An integer compact-support certificate proves this; a finite-box eigensolver is only
used to propose trial vectors and measure profiles. The commuting two-link control has
exactly one such mode, despite the same total conserved charge $Q=4$.

These modes are not automatically persistent. A sharper charge bound explains how
actual reactions can force their loss:

$$n_+(L_{\mathrm{bundle}}-12I)\le\min(Q_\uparrow,Q_\downarrow)
\le\lfloor Q/2\rfloor.$$

The face-orientation charges can redistribute even though their sum is conserved.
Putting all charge on one orientation prohibits **every** above-band mode. The saved
raw dynamics contains exactly replayed examples of that event.

These are graph spectral statements, not atomic energies or a physical wave equation.
The supplied base is still a two-dimensional triangular lattice, the fiber is supplied
triangle adjacency, and the initial link defects are specified. Above-band localization
is not low-energy binding, particle identity, emergent geometry, or quantum matter.

![Measured full-graph mode densities, scaling, and exact dynamics](images/triangle-modes.png)

## 1. A certificate for the entire infinite graph

For the flat base, the horizontal Laplacian symbol is

$$\lambda_\triangle(k_x,k_y)
=6-2\bigl(\cos k_x+\cos k_y+\cos(k_x+k_y)\bigr)\in[0,9].$$

The fiber $C_3$ has Laplacian eigenvalues $0,3,3$, so the full flat graph has spectrum
$[0,12]$. Unlike the earlier square-fiber modes, the new modes lie outside the **entire**
flat spectrum rather than being embedded in another fiber sector's continuum.

Change just the links $(0,0)\to(1,0)$ and $(1,0)\to(2,0)$ to the two derived noncommuting
reflections (group enumeration indices one and two). Write $L=L_0+\Delta$.
The difference $\Delta$ is supported on three base vertices, hence represented by an
exact $9\times9$ integer matrix. Congruence elimination gives positive index two,
negative index two, and nullity five. Since $L_0\le12I$, any subspace above twelve is
positive for $\Delta$; its dimension is therefore at most two.

For the lower bound, use a $13\times13$ base patch with all three fiber vertices per
site: 507 full graph vertices. **Retain diagonal degree eight at its boundary.** This is
the principal restriction of the infinite operator, not the Laplacian of the induced
finite subgraph. Extend its trial columns by zero outside the patch.

For the saved two integer columns $X$, exact arithmetic gives

$$G=X^T(L-12I)X=\begin{pmatrix}958554&490\\490&4204372\end{pmatrix},$$
$$M=X^TX=\begin{pmatrix}16778960&1955\\1955&16772604\end{pmatrix}.$$

Both leading principal minors of $20G-M$ are positive. Every nonzero vector in this
two-dimensional trial space consequently has Rayleigh quotient greater than
$12+1/20$. Finite-rank perturbations preserve the essential spectrum, so these are
two discrete eigenvalues, not a shifted continuum. Together with the upper bound this
proves the exact count. The underlying general spectral results are the min-max principle
and Weyl's theorem; see [Teschl, sections 4.4 and 6.4](https://www.mat.univie.ac.at/~gerald/ftp/book-schroe/schroe2.pdf).
Our contribution is the explicit full-graph construction and verifiable integer witness,
not a claim to have invented these spectral methods.

The verifier reconstructs the operator and perturbation, recomputes both forms, and
checks every inequality without calling an eigensolver. Signed-integer multiplication
uses an explicit absolute-intermediate bound; oversized trial integers are rejected.
All trial columns also have exactly zero sum within each fiber. This is a derived
invariant subspace, not a separately supplied two-component wave law: permutation
transport preserves the constant fiber vector, whose scalar sector stays below nine.

The same procedure supplies complete controls:

| Specified links | Exact number above the full flat band | Strict bound for each existing mode |
| :--- | ---: | :--- |
| Flat | 0 | — |
| One reflection link | 1 | $>12.05$ |
| Two equal reflection links | 1 | $>12.05$ |
| Two noncommuting reflection links | 2 | $>12.05$ |

The equal-reflection perturbation has positive index one, so the single-mode count is
not an assertion based on failing to find a second numerical eigenpair.

## 2. Localization follows from a resolvent bound

Let $\psi$ be a normalized eigenvector with $\lambda>12$, and put
$\eta=\Delta\psi$, supported on the three defect endpoints $S$. Since
$\|L_0-6I\|=6$, the norm-convergent expansion

$$\psi=\frac1{\lambda-6}\sum_{n\ge0}
\left(\frac{L_0-6I}{\lambda-6}\right)^n\eta$$

and finite propagation in the base graph imply

$$\|\mathbf1_{B_R(S)^c}\psi\|
\le\frac4{\lambda-12}
\left(\frac6{\lambda-6}\right)^{R+1}.$$

Each changed link contributes an operator of norm at most two, giving
$\|\eta\|\le\|\Delta\|\le4$. The certificate's gap gives the uniform, conservative
bound $80(120/121)^{R+1}$. This proves exponential spatial tails; it is not a fitted
localization length, a physical propagation velocity, or a statement about time evolution.

## 3. Finite-size profiles are measurements, not the proof

Embed the same integer certificate unchanged in tori of sides $16,24,48,96$.
Its two quadratic forms agree exactly with the infinite principal restriction in every
case. The perturbation upper bound then proves the finite mode counts as well.

Separately estimate eigenpairs and the normalized base-site projector density
$p_v=\sum_{a=0}^2|\psi_{v,a}|^2$. The participation volume
$V_{\mathrm{eff}}=1/\sum_vp_v^2$ measures concentration. It is invariant under local fiber
permutations, but is not a Born-probability assumption. For the noncommuting pair:

| Side | Lower eigenvalue estimate | Its participation volume | Upper eigenvalue estimate | Its participation volume |
| ---: | ---: | ---: | ---: | ---: |
| 48 | 12.0672501495 | 25.135744 | 12.2519080026 | 9.032702 |
| 96 | 12.0672500278 | 25.135135 | 12.2519080026 | 9.032702 |

The largest graph has 27,648 vertices. These decimals and densities are numerical
estimates with residual checks, not certified enclosures of infinite eigenvalues.
The exact existence/count claim comes from the compact integer witness instead.
The heatmaps show the side-24 measurements in supplied coordinates, with the two
changed base edges marked. They are not an emergent spatial embedding.

## 4. A stronger bound from the two face orientations

Each base edge belongs to exactly one upward and one downward triangle. Each base
vertex belongs to three triangles of each orientation. Consequently **either orientation
alone** gives the exact full-graph identity

$$12I-L_{\mathrm{bundle}}=
\sum_{f\in\uparrow}\iota_f S_f\iota_f^T+\bigoplus_v J_3
=\sum_{f\in\downarrow}\iota_f S_f\iota_f^T+\bigoplus_v J_3.$$

The [face inertia calculation](triangle-feedback.md#5-the-same-charge-bounds-full-graph-spectral-capacity)
gives $n_-(S_f)=q(H_f)$, where $q$ is the charge derived from the reaction rules.
Intersecting the nonnegative subspaces of the face forms separately for each orientation
therefore yields

$$n_+(L_{\mathrm{bundle}}-12I)\le\min(Q_\uparrow,Q_\downarrow),
\qquad Q_\uparrow+Q_\downarrow=Q.$$

This strengthens the earlier total-charge bound by at least a factor of two. Both
orientation identities are checked against an independently assembled integer graph
matrix. Only the sum $Q$ is conserved by the mixed reaction/transport dynamics; a
vacancy move can transfer charge between orientations, and a reaction can redistribute it.
Positive capacity is only permission for modes to exist, not a guarantee.

## 5. Every changing event is counted in a recorded evolution

Replay the [unrouted reaction experiment](triangle-feedback.md#4-unrouted-evolution-and-negative-controls)
with $Q=4$ on its supplied $6\times6$ torus. For scheduler seed zero, compute the exact
characteristic polynomial of the full 108-vertex graph after **every changing event**,
not merely every thousandth tick. Real symmetry makes sign variations of the shifted
characteristic polynomial an exact count above twelve. Initial and final endpoints
complete the history; no-op attempts leave the graph unchanged.

Also audit both endpoints of every reaction in all four recorded seeds. Independent
fraction-free integer congruence agrees with the characteristic-polynomial counts at
all 24 endpoints. Every reaction receives an isolated C++ replay and exact inverse,
and the full C++ histories are revalidated before extracting the endpoints.

The orientation bound forces every observed downward reaction crossing. For example:

$$\begin{array}{c|c|c}
\text{tick}&(Q_\uparrow,Q_\downarrow)&\text{above-band count}\\\hline
294&(1,3)\to(2,2)&1\to2\\
545&(2,2)\to(3,1)&2\to1\\
12579&(1,3)\to(0,4)&1\to0\\
18123&(4,0)\to(3,1)&0\to1\\
18296&(3,1)\to(4,0)&1\to0
\end{array}$$

The upward crossings are measured, not guaranteed by the bound. Other recorded reactions
raise capacity without creating another mode. The complete exported audit includes all
twelve reactions, not only the examples above, and all 860 changing events of seed zero.

Over its 100,000 attempted-update intervals, seed zero spends 26,940 intervals with no
above-band modes, 64,757 with one, and 8,303 with two. The longest zero-mode interval is
$[62917,64393)$, of length 1,476. These are exact residence counts for this particular
history, not equilibrium probabilities or physical lifetimes. They also do not track the
identity of an individual eigenmode. There are zero-mode states with positive orientation
capacity, so the upper bound is not an exact mode-count formula.

This rejects identifying this model's above-band mode count with a generally conserved
particle number. It does not rule out metastable structures under other initial
conditions or other derived rules. Establishing stable matter, calibrated interactions,
physical time, and molecular behavior remains unfinished.

The [full-generator follow-up](triangle-mode-rates.md) proves that *every* vacancy move
from a capacity-saturated state forces mode loss, and computes exact conditional kernels
instead of inferring a decay law from one history.

## Reproduce

```bash
cmake --build build -j 4
export OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
uv run --with numpy --with scipy --with python-flint python tools/triangle_modes.py --verify-only
uv run --with numpy --with scipy --with python-flint python tests/triangle_modes.py
uv run --with numpy --with scipy --with python-flint python tools/triangle_modes.py \
  --output out/triangle-modes.json
uv run --with numpy --with scipy --with matplotlib python tools/plot_triangle_modes.py
```

The compact certificate is verified exactly, not regenerated in CI. Re-proposing its
integer vectors uses `--generate-certificate out/triangle-compact-modes.json`; floating
eigensolver choices can change those proposed vectors without invalidating the saved
certificate. Likewise, numerical profile decimals need tolerance-based comparisons.
The dynamical counts, raw connections, and compact inequalities are exact.
