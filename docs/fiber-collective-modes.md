# Collective subgroup memory in the full local bank

Follow-up: [which of these internal dynamics can charge actually
observe?](fiber-charge-coupling.md) Exact blindness on a single fan and
an additional prime-fiber obstruction qualify the physical interpretation
of the slow modes below.

Restoring all twelve selected reactions and both neighboring pair supports
changes the isolated-pair result. The complete fixed-boundary,
charge-three component has

$$\boxed{D_n=(n-1)(n+7)\text{ framed states}.}$$

Only a fraction $6/(n+7)$ contain a rotation face; the others store the
charge in three reflections. Their elastic dynamics cannot be replaced
by a passive waiting state.

With adjacency-assisted angular moves the component is connected, but
develops slow **gauge-invariant** modes increasingly associated with the
old holonomy-subgroup sectors. At $C_{243}$, sector labels retain 95.1%
of the slowest mode's squared norm, yet their averaged transition rates
overestimate its relaxation rate by a factor 8.62. We also obtain

$$\boxed{\limsup_{n=3^a\to\infty}n\gamma_n\le1458/73}$$

for the continuous-time spectral gap. A relaxation rate therefore vanishes
with fiber refinement. This is not a proven scaling exponent, a spatial
phase transition, or a massless particle: the base region stays at three
faces while its internal fiber resolution increases.

## The actual dynamics and state count

On a three-face fan, include all twelve relational reactions; vacancy,
Hurwitz, and inverse Hurwitz on each adjacent pair; and all
$m=(n-1)/2$ angular involutions on each pair. There are $n+17$ channels.
Reported continuous-time rates assign each channel a rate-one clock.
Uniformly proposing one channel divides the generator by $n+17$;
idle proposals still count. Exterior links and the boundary product stay
fixed. This is the full **chosen** bank on one fan, not every finite law
or all interactions of an extended mesh.

Fix the boundary product to $P=r_0$ and charge to $Q=3$. The possible
face inventories are $RRR$ and $ERZ$ in any order. Reflection tuples are

$$(r_a,r_b,r_{-a+b}),$$

giving $n^2$ possibilities. Exclude the frozen parallel tuple $(P,P,P)$.
An $ERZ$ tuple has one of $n-1$ rotation entries and six layouts; the
reflection is fixed by its product. Hence

$$D_n=(n^2-1)+6(n-1)=(n-1)(n+7).$$

The two-element centralizer of $P$ acts freely on these nonparallel
tuples. The gauge quotient has $D_n/2$ states. All reported gaps are
computed on this quotient, not from slow changes of a chosen frame.
The framed sizes at $n=9,27,81,243$ are $128,884,7040,60500$; the largest
eigenproblem has 30,250 gauge states.

## Classify and connect the old subgroup sectors

Here $q$ denotes the order of the generated **rotation subgroup**, not
the face charge. The full holonomy image has order $2q$. For an $RRR$
tuple it is

$$q=n/\gcd(n,a,b);$$

for an $ERZ$ tuple with rotation $t_k$, it is $n/\gcd(n,k)$. With $P$
fixed this identifies the actual rooted subgroup.

The old transports and reactions preserve $q$. On reflection labels,
the adjacent Hurwitz moves are

$$M_1=\begin{pmatrix}2&-1\\1&0\end{pmatrix},\qquad
M_2=\begin{pmatrix}1&0\\1&1\end{pmatrix}.$$

The conjugate $M_2^{-1}M_1M_2$ is the elementary upper shear, while $M_2$
is the elementary lower shear. They generate
$\operatorname{SL}(2,\mathbb Z_{3^a})$: after removing a common power of
three, a primitive vector has a unit coordinate, and elementary row
operations reduce and rescale it to a coordinate vector. Thus the
elastic orbits are exactly the common-divisor classes. Every $ERZ$ state
attaches to its reflection orbit through an inverse reaction.

For $q\mid n$, $q>1$, the sector size is consequently

$$\boxed{D_q^{\rm sector}=J_2(q)+6\varphi(q),}$$

where $J_2(q)$ counts primitive ordered pairs and $\varphi(q)$ counts
primitive rotation steps. For $q=3^b$ these are $8q^2/9$ and $2q/3$.

An angular move from any proper sector changes a rotation step divisible
by three to a neighboring step coprime to $n$. It goes directly to the
full sector $q=n$. The inverse connects back. The subgroup connectivity
graph is a star and the augmented component is connected. Connectivity
does not imply that motion on the star is Markovian.

## Exact fluxes do not give a closed sector process

Equal inverse-channel rates make the stationary measure uniform. In a
proper sector there are $\varphi(q)$ possible rotation entries and four
layouts where $R$ and $Z$ are adjacent. Each has two angular exits. Thus

$$F_{q,n}=F_{n,q}=8\varphi(q),\qquad
r_q=\frac{8\varphi(q)}{D_q^{\rm sector}}=\frac{12}{2q+9}.$$

The mean $r_q$ is not the rate at each microscopic state. Exit rate is
zero in every $RRR$ state and separated mixed layout, but two in each
adjacent mixed layout. Its conditional variance is

$$\operatorname{Var}(r\mid q)=2r_q-r_q^2>0.$$

Two states with the same subgroup label have different immediate exit
rates. This proves nonclosure, without a trajectory fit. Averaging
preserves stationary flux but loses information determining future exits.

Slot counting also gives a rotation-containing fraction $6/(n+7)$,
total angular activity $8(n-2)/[(n-1)(n+7)]$, and reaction activity
$48/(n+7)$ in the per-channel clock. Multiplying an isolated angular
diffusion coefficient by the rotation-containing fraction would still be
unjustified: the reflection configurations have nontrivial elastic motion.

## The full gauge-invariant modes

Let $L$ be the positive generator, with evolution $e^{-tL}$. Assemble
it from the actual primitive maps, then quotient by the remaining
fixed-boundary conjugation. Direct sparse diagonalization gives:

| Fiber | Full gauge-invariant gap | Sector-averaged gap | Slow-mode norm retained by sectors |
|---|---:|---:|---:|
| 9 | 0.307335 | 0.948148 | 43.8% |
| 27 | 0.100933 | 0.503704 | 67.4% |
| 81 | 0.0296522 | 0.213009 | 87.0% |
| 243 | 0.00905163 | 0.0780452 | 95.1% |

The averaged rates are exactly $F_{q,r}/D_q^{\rm sector}$, not fitted
coefficients. Their eigenvalues are Ritz values obtained by restricting
to sector-constant observables: upper bounds on corresponding full
eigenvalues, not proof of Markovian sector dynamics.

At $n=243$, the averaged gap is 8.62 times too large. Projecting the
actual slow eigenvector onto sector constants gives Rayleigh quotient
$0.0797269$, versus its true $0.00905163$. Most mode variance is retained,
but the omitted component contains important transition information.

Rotation-containing and reaction-enabled reflection states occupy only
$8/(n+7)=3.2\%$ of this component, yet account for about 80% of the
omitted squared norm. The true slow mode's Dirichlet energy is
$0.00756007$ in other within-sector moves, $0.000462703$ in within-sector
reactions, and $0.00102886$ across subgroups. Internal gradients greatly
reduce inter-sector flux. Variance-only compression misses this effect.

Numerical eigenpair residual norms are below $3\times10^{-12}$. This
checks the finite eigenproblem, not a continuum exponent or the physical
validity of the chosen microscopic law.

## A general vanishing-gap bound

For a proper sector $A$, the gauge-invariant indicator has Rayleigh
quotient $r_A/(1-\pi_A)$, where $\pi_A=D_A/D_n$. Choose $q=n/3$:

$$D_A=8n^2/81+4n/3,\qquad r_A=36/(2n+27).$$

Therefore

$$\gamma_n\le\frac{36}{2n+27}
\left(1-\frac{8n^2/81+4n/3}{n^2+6n-7}\right)^{-1},$$

which yields $\limsup n\gamma_n\le1458/73$. In the uniform proposal
clock the rate is divided by $n+17$, giving an $O(n^{-2})$ upper bound.
These are upper bounds on rates, or lower bounds on relaxation times.
They do not prove matching lower bounds, an asymptotic power law, or
spatial criticality.

## First-passage memory measured on shared links

Kill the process when it leaves a proper sector $A$. Its principal
generator submatrix $K=L_{A,A}$ is positive definite. The microscopic
mean first-exit times solve

$$Ku=\mathbf1,\qquad
\mathbb E_A T=|A|^{-1}\mathbf1^T K^{-1}\mathbf1$$

for an initially uniform state inside $A$. This is generally not $1/r_A$.
At $C_{243}$, numerical solutions of these exact backward equations give:

| Rotation subgroup order | Inverse mean exit rate | Mean first exit from uniform sector |
|---|---:|---:|
| 3 | 1.25 | 1.48077 |
| 9 | 2.25 | 5.08761 |
| 27 | 5.25 | 23.0192 |
| 81 | 14.25 | 85.9154 |

The same proper-sector exit problem occurs in any larger power-of-three
fiber: its internal legacy dynamics depends on $q$, and angular killing
has rate two on the same mixed configurations.

Raw-link trajectories were initialized uniformly in each fiber's largest
proper sector, then evolved with all $n+17$ channels until first exit.
There were 64 independent trajectories per fiber and 1,731,216 total
proposals. Times are proposals divided by $n+17$:

- $C_9$, $q=3$: $1.694\pm0.204$, predicted $1.481$.
- $C_{27}$, $q=9$: $4.466\pm0.588$, predicted $5.088$.
- $C_{81}$, $q=27$: $19.325\pm2.242$, predicted $23.019$.
- $C_{243}$, $q=81$: $95.830\pm12.048$, predicted $85.915$.

Errors are standard errors across independent trajectories. Initial and
exit holonomy images are computed from all rooted edge transports, not
merely inferred from an abstract sector label. State preparation is
separate from evolution; no creation from a flat vacuum is asserted.

The survival function is

$$S_A(t)=|A|^{-1}\mathbf1^T e^{-tK}\mathbf1.$$

It is a positive mixture of exponential decays. Its conditional escape
hazard decreases as faster components are depleted, even though the
parent process has a stationary equilibrium. The initial hazard is
$r_A$ and its initial derivative is $-(2r_A-r_A^2)$. Thus a decreasing
hazard alone is not evidence of nonequilibrium aging or a glass transition:
unresolved configurations already produce it here.

## A kinetic observer improves on the subgroup label

Extend $u=K^{-1}\mathbf1$ by zero outside $A$ and subtract its global
mean. It accounts for the microscopic expected time to escape rather
than assigning every state in $A$ the same value. Its Rayleigh quotient is

$$\mathcal R[u]=\frac{\mathbb E_A u}
{\mathbb E_A u^2-\pi_A(\mathbb E_A u)^2},\qquad
\gamma_n\le\mathcal R[u].$$

At $C_{243}$, $q=81$, this improves the indicator bound from $0.0781248$
to $0.0121703$, much closer to the actual gap $0.00905163$. The inequality
is the variational principle; the reported value is a numerical evaluation
using the finite backward-equation solution. This supplies a concrete
observer that respects kinetics rather than just explaining variance.

## Interpretation and open work

The isolated reaction-plus-angular chain's factor-of-two diffusion
correction must not be promoted to a property of the full bank. The
larger reflection configuration space adds subgroup memory and slower
gauge-invariant modes. We now have an exact state count, exact stationary
flows, a nonclosure witness, a general gap bound, and raw first-exit data.

These are properties of a classical gauge-covariant stochastic model on
a fixed local base region. They do not establish propagating fields,
bound states, masses, quantum interference, or evolving spacetime. A
vanishing Markov relaxation rate is not being identified with a physical
particle mass. Extended spatial interactions and a controlled joint
spatial/fiber limit remain open.

The underlying tools are established: see Levin and Peres,
[Markov Chains and Mixing Times](https://pages.uoregon.edu/dlevin/MARKOV/),
and [Quantitative coarse-graining of Markov chains](https://arxiv.org/abs/2201.10256).
No independent field-wide novelty review of this application has been made.

Calculation: [fiber_collective_modes.py](../tools/fiber_collective_modes.py).
Spectra, flows, backward equations, and raw first-exit observations:
[fiber-collective-modes.json](../data/fiber-collective-modes.json).
