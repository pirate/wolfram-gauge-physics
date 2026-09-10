# Boundary alignment survives the isolated gauge quotient

The 40-state projective reflection chart is insufficient when a region
interacts with its environment. We now have an explicit boundary-framed
lift, a closed internal circuit that changes exterior-relative alignment,
and a primitive interaction that reads that alignment through charge
fluctuations. All use the existing shared-edge rules.

This is classical gauge memory and a controlled interaction experiment.
It is not a quantum phase, a fermion, or evidence of a bound particle.

## Keep the boundary and the external reference

Use the same five-face fan as in
[the overlap calculation](triangle-overlap-algebra.md), with ordered face
IDs $(70,71,10,1,0)$. Hold all seven actual outer links fixed, with ordered
boundary holonomy $p=1$, the first derived reflection. The region has four
internal links. Solving its native face holonomy equations reconstructs
those links uniquely for each compatible ordered tuple.

At $Q=5$, full $S_3$ image, and this fixed boundary connection, there are
**560 framed configurations**, including 80 all-reflection configurations.
An actual exterior triangle, face 61, has no internal edge of the region;
its independent outer edge allows its holonomy to be fixed as a reference.
Initially choose the derived rotation $z=3$ there. Every contained update
leaves that reference and all seven boundary links untouched.

The residual centralizer of the reflection boundary is

$$C(p)=\{e,p\}\cong\mathbb Z_2.$$

It acts on internal words by

$$\tau(a_0,\ldots,a_4)=(pa_0p^{-1},\ldots,pa_4p^{-1}).$$

There are no fixed points in the full-$S_3$ sector, so forgetting the
boundary-relative alignment identifies pairs and gives the previous 280
states, including 40 reflection states. Every one of the framed primitive
transitions projects exactly to the previously constructed unframed
transition. The 560-state framed transition graph is connected as well.

Keeping $z$ fixed while conjugating the interior by $p$ is **not** a global
gauge transformation. The joint stabilizer of $z$ and $p$ is trivial. In
the reflection chart, this distinguishes the vectors $v$ and $-v$ that
the isolated projective chart identifies.

For example, the observables

$$
O_i=\frac{q(zpa_i)-q(z^{-1}pa_i)}2
$$

are common-root, jointly gauge-invariant relative-alignment measurements.
They change sign under $\tau$ with $z$ held fixed. They are ordinary real
observables on classical configurations, not probability amplitudes.

## A closed projective circuit with a nontrivial framed result

Apply the elastic Hurwitz rule on pair supports $(10,2,4,6)$, in that order,
and repeat this four-update sequence five times. The resulting 20-update
circuit acts on every one of the 80 reflection configurations as $\tau$.
In the four-coordinate description its matrix is $-I$.

The internal projective state therefore returns exactly, but its alignment
with the exterior reference reverses. One actual example is

$$
(1,1,1,2,2)\longmapsto(1,1,1,5,5),
$$

with boundary product $p=1$ and exterior rotation $z=3$ unchanged. The
five-component observable $O$ changes from $(0,0,0,-1,-1)$ to
$(0,0,0,1,1)$. The raw-link circuit was checked for all 80 reflection inputs,
not just this example.

This is a concrete configuration-space holonomy. Repeating the circuit
twice gives identity on this sector because $p^2=e$. It is not a spatial
$2\pi$ rotation and does not establish spin-half behavior.

For any sequence of contained operators, conjugate initial words $a$ and
$\tau a$ have identical individual face-charge histories: each operator
commutes with $\tau$, and charges are conjugacy invariant. Their internal
unframed gauge histories likewise agree. An external comparison can still
distinguish them.

## An actual interaction reads the alignment

In a separate readout experiment, fix the exterior face 61 to reflection
$s=2$ instead of rotation $z=3$. This reflection does not commute with the
held boundary reflection $p=1$. Use the same interior input and the same
20-update circuit. Before the readout, every individual face charge is
identical in the two histories, and their internal projective states agree.

Now apply a boundary-crossing reaction on fan support 4. Its three face
IDs, in table order, are $(1,0,61)$:

- Before the twist, its based tuple is $(2,2,2)$, which is nonreactive.
- After the twist, its based tuple is $(5,5,2)$, which is reactive.

This is not dependent on picking a special rule. With uniform choice
among all twelve reaction rules, the first input stays at charge $(1,1,1)$
with probability one. The second produces each permutation of $(0,1,2)$
with probability $1/6$.

Consequently

$$
\sum_i q_i=3\quad\text{in both experiments},\qquad
\mathbb E[q_i]=1\quad\text{for each of the three faces},
$$

but

$$
\sum_i q_i^2=
\begin{cases}
3 &\text{before the twist},\\
5 &\text{after the twist},
\end{cases}
$$

with certainty. Relative boundary alignment thus has an exact charge-
fluctuation readout even though the mean charge readout misses it. This
does not violate gauge invariance: the exterior reference participates
in the interaction, and was never conjugated along with the interior.

## Exact return operator with both boundary characters

The framed transition matrix commutes with $\tau$. Choose one representative
$x_i$ of each pair $\{x_i,\tau x_i\}$, and split the operator into its two
centralizer-character blocks:

$$
P^+_{ij}=P(x_i,x_j)+P(x_i,\tau x_j),\qquad
P^-_{ij}=P(x_i,x_j)-P(x_i,\tau x_j).
$$

$P^+$ is the previous 280-state stochastic matrix. $P^-$ propagates the
boundary-odd observables. Its signed entries are differences of classical
probabilities, not negative probabilities or quantum amplitudes. Both
matrices are symmetric in the uniform framed reference.

Split each block into its 40 reflection representatives $R$ and 240 mixed
representatives $E$. Keeping the first-return clock gives

$$
F^\pm(z)=zP^\pm_{RR}
+z^2P^\pm_{RE}(I-zP^\pm_{EE})^{-1}P^\pm_{ER}.
$$

The two framed return probabilities are reconstructed by
$(F^++F^-)/2$ and $(F^+-F^-)/2$. Elimination at $z=1$ was done over exact
rationals; every reconstructed entry is nonnegative. The boundary-odd
mixed-excursion operator is nonzero, so retaining only the even block
really would omit reaction-return information.

For a uniform initial reflection configuration, the rounded probabilities
of returning to the **same projective point** are:

| First-return event | Same alignment | Reversed alignment |
| --- | ---: | ---: |
| Any first positive return | 0.584366896 | 0.006650747 |
| First return after departing into the mixed sector | 0.034366896 | 0.006650747 |

The second row is **not** normalized by the probability of departure.
Conditional on both a mixed excursion and return to the same projective
point, the reversal probability is approximately $0.162143575$.
These are exact finite-chain probabilities rounded for display, not Monte
Carlo estimates. Returns to other projective points are outside this table.

The largest nonconstant eigenvalue of the even one-attempt block is
approximately $0.990513796$; the largest eigenvalue of the odd block is
$0.990494925$ (numerical diagonalization). Boundary-relative information
can therefore relax on a timescale comparable to the slow internal modes.
Neither value is a whole-mesh dispersion law or a molecular lifetime.

Return-alignment labels between *different* projective points depend on
the chosen representatives. The diagonal return probabilities above and
the spectra do not. No path-independent identification of remote frames
has been assumed.

## A derived gluing coordinate, with limits

For two full-$S_3$ interiors matched across a reflection boundary, their
interior stabilizers are trivial and the available relative alignments
come from $C(p)=\mathbb Z_2$. Thus the extra binary gluing coordinate is
derived from the boundary centralizer, not declared as a new fundamental
gauge group. The even/odd operator split is its exact character decomposition.

This does not justify putting independent binary variables on an arbitrary
coarse lattice. Multiple overlaps impose compatibility constraints, and
boundary charge changes can change the centralizer itself. Those must be
derived from the actual shared graph when coupling evolving regions.

The present result establishes a controlled, locally readable classical
alignment memory and gives its exact return dynamics. It does not establish
spontaneous coherent structures, quantum interference, binding, or emergent
spacetime. The next mathematical task is to compose these framed operators
across actual evolving overlaps without discarding their compatibility or
elapsed-time information.

The first released-neighbor calculation and its derived short-time
variance/current law are now in
[Released overlap](triangle-variance-current.md).

## Reproduce

```sh
OPENBLAS_NUM_THREADS=1 uv run --with numpy --with scipy --with python-flint python tools/triangle_framed_overlap.py
```

This reconstructs all 560 boundary-held raw configurations, applies every
contained primitive, verifies projection to the old quotient, runs all 80
twist circuits, finds the boundary-crossing readout, and computes both exact
return blocks. It uses the unchanged reaction and elastic tables.
