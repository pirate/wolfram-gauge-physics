# Exact charge memory: 42 activities, 39 coupling channels

The variance-directed current does not close the charge dynamics. An exact
projection of the released six-face region retains the missing gauge memory
and shows precisely how it enters the charge evolution.

There are **121 complete charge patterns** in the 2,943-state region. Their
instantaneous hidden-gauge coupling has exact rank **39**. It factors through
42 independent conditional reaction activities, with exactly three output
relations arising from disjoint, commuting reaction supports.

The reduced calculation reproduces the prepared alignment response. It is
not a new Markov process with 39 states: the 39 channels carry a time-dependent
memory kernel and a term determined by the initial hidden state.

## Conditional projection, with the right measure

Use the unchanged 108-slot operator $P$ on the connected, boundary-held,
full-$S_3$, $Q=6$ region from
[the released-overlap calculation](triangle-variance-current.md). Its
stationary measure is uniform on 2,943 gauge-fixed configurations.

The 121 charge fibers have sizes 243 (one fiber), 54 (30 fibers), and 12
(90 fibers). Let

$$S_{xq}=\frac{\mathbf1_{q(x)=q}}{\sqrt{n_q}},\qquad
\Pi=SS^T,\qquad Q=I-\Pi.$$

In these Euclidean coordinates, equivalent to the uniform stationary
$L^2$ coordinates up to a common scale, $S^TS=I$. Thus $\Pi$ is conditional
expectation given the **entire six-face charge field**, not just total
charge or a selected Fourier component.

Define

$$A=S^TPS,\qquad F=QPS,\qquad D=QPQ.$$

For a state-distribution perturbation, write its charge component as $x_t$
and its hidden component as $y_t$. Since $P$ is symmetric,

$$x_{t+1}=Ax_t+F^Ty_t,\qquad y_{t+1}=Fx_t+Dy_t.$$

Eliminating $y$ gives the exact discrete memory equation

$$
x_{t+1}=Ax_t+\sum_{j=0}^{t-1}K_jx_{t-1-j}+F^TD^ty_0,
\qquad K_j=F^TD^jF.
$$

The last term is essential for the twisted/untwisted preparations: their
initial charge difference is zero, but their hidden difference is not.
Discarding that term would discard the very alignment response being studied.

This is a finite-state implementation of the established projection-memory
approach, not a new projection formalism. For the statistical-mechanical
framework see [Mori, *Transport, Collective Motion, and Brownian Motion*
(1965)](https://academic.oup.com/ptp/article/33/3/423/1925580). The operator
identity used here follows directly from the two block equations above.

## Where the rank 39 comes from

A gate can affect a one-step charge transition only when its three faces
currently have charge $(1,1,1)$. For each compatible charge pattern $q$ and
fan $g$, introduce the conditional activity observable

$$h_{qg}(X)=\mathbf1_{q(X)=q}
\left[a_g(X)-\mathbb E[a_g\mid q]\right].$$

There are 42 such observables in this region: six in the all-reflection
fiber and 36 in fibers with one rotation and one vacancy. Their centered
function space has exact rank 42.

For an active gate, its charge-output count vector is

$$v_{qg}=2\sum_{\sigma\in S_3}e_{q^{g,\sigma}}-12e_q,$$

where $q^{g,\sigma}$ replaces the three charges on $g$ by the permutation
$\sigma(0,1,2)$. The entire hidden part of the one-step charge transition
counts factors as the conditional activity matrix times the matrix of
these 42 vectors. Elastic updates and the other charge transitions depend
only on $q$ and disappear under this conditional centering.

### Three relations from commuting reactions

The disjoint fan pairs are $(0,1)$, $(2,4)$, and $(3,5)$. Their actual raw
write sets miss each other's read sets, so the corresponding primitive
updates commute before any gauge fixing.

For one such pair $(g,h)$, with $q_0=(1,1,1,1,1,1)$, the charge-output
vectors obey

$$
6(v_{q_0g}-v_{q_0h})
+\sum_{\sigma\in S_3}v_{q_0^{g,\sigma},h}
-\sum_{\sigma\in S_3}v_{q_0^{h,\sigma},g}=0.
$$

To see this, expand each vector. Both two-reaction terms contain the same
sum over all 36 final charge assignments, because the supports are disjoint.
Their intermediate-state terms cancel against the first term. This is a
relation between ordinary path counts, not interference of amplitudes.

The three relations are independent and exhaust the output matrix's left
nullspace. Its exact rank is $42-3=39$, and the rank of $F$ is also 39.
These are integer-rank calculations, not small singular values discarded
at a numerical tolerance.

Eliminating one all-reflection port from each relation gives an explicit
rational factorization

$$F=UW,\qquad U\in\mathbb R^{2943\times39},\quad
W\in\mathbb R^{39\times121},$$

and hence

$$K_j=W^T\underbrace{(U^TD^jU)}_{39\times39}W.$$

The implementation uses those 39 channels to propagate the memory. Since
$K_0=F^TF$ has rank 39, no smaller linear channel space can factor this
complete instantaneous coupling exactly.

This does **not** remove all but 39 hidden states. A combination invisible
to $F^T$ now can become visible after evolution under $D$. The memory
propagator retains that effect.

## Individual mean charges are initially blind

Every vector representing an individual charge $q_i$ lies in the nullspace
of $F$. Thus a purely hidden perturbation cannot change an individual mean
charge in one attempt. It can change the full charge distribution, including
its variance, in one attempt.

The reduced memory equation, including its initial-hidden-state term,
reproduces the two-step mean difference

$$\frac1{108^2}(-16,8,-8,16,-8,8)$$

for the prepared alignment pair. The complete projected propagator and the
prepared response are checked against the microscopic matrix through nine
attempts. The arithmetic ranks and factorization identities are exact;
the normalized memory propagation uses floating-point arithmetic.

## Why this feedback is not an unstable force law

For the current inverse-paired uniform scheduler, $P$ is a reversible Markov
contraction in its stationary $L^2$ space. In particular,

$$-I\leq P\leq I,\qquad \|P^tf\|\leq\|f\|.$$

Eliminating hidden variables cannot introduce an exponentially growing
linear mode that the full operator lacks. More explicitly, for $s>1$,

$$sI-A-F^T(sI-D)^{-1}F$$

is positive definite: it is the Schur complement of $sI-P$. Therefore the
exact charge-memory equation has no real growth pole $s>1$. The full
self-adjoint contraction also rules out complex eigenvalues outside the
unit disk.

There is a useful bound on transient transfer. At every even lag,
$T=P^{2n}$ satisfies $0\leq T\leq I$. Hence

$$\|\Pi P^{2n}Q\|
=\|\Pi(T-\tfrac12I)Q\|\leq\tfrac12.$$

This bounds transfer from a unit stationary-$L^2$ hidden perturbation into
the **complete charge component**, not a raw charge magnitude or a reaction
probability. In the six-face calculation, the transfer norms are approximately
0.21075 after two attempts, 0.22291 after four, and 0.03940 after 64.

The exact transfer norm can be evaluated without a full hidden basis:

$$\|\Pi P^tQ\|^2
=\lambda_{\max}(C_{2t}-C_t^2),\qquad C_t=S^TP^tS.$$

Thus the prepared charge response can grow transiently while the complete
perturbation contracts. Treating the positive variance-gradient term as
an independent self-amplifying force would omit the memory and damping
that enforce this constraint.

These are applications of standard reversible-operator facts. They do
not exclude metastable configurations, a gap closing with increasing
system size, or thermodynamic ordering. They also do not establish
binding, quantum dynamics, or a continuum limit. They rule out an
unstable linear feedback mode of the exact current probability operator.

## Scope and next mathematical question

The result is an exact, geometry-derived reduction of how hidden gauge
activity affects a complete charge field. It keeps the initial alignment
source and delayed feedback; it supplies no new update law.

The next issue is the lifetime and spatial extent of those memory channels
as regions grow and exchange boundary data. A long-lived collective mode
would need evidence of that scaling or metastability, rather than growth
manufactured by truncating the short-time current equation.

```sh
OPENBLAS_NUM_THREADS=1 uv run --with numpy --with scipy --with python-flint python tools/triangle_charge_memory_kernel.py
```
