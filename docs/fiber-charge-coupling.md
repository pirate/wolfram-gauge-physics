# Which internal gauge dynamics can charge actually observe?

Follow-up: [relative holonomy becomes observable as charge noise](fiber-relative-angle.md)
constructs the first non-blind prime-fiber sectors, including a five-face
identity-boundary example, and derives their exact drift and covariance.

The slow modes in [the full local bank](fiber-collective-modes.md) do not
automatically describe charge carriers. Here are two exact obstructions and
a small but nonzero coupling once patches overlap. These are results about
our selected classical model, not a derivation of a physical force.

## Model and clocks

The fiber is the supplied cycle $C_n$, with derived automorphism group
$D_n$. Write $E$ for identity, $R$ for a reflection, and $Z$ for a
nonidentity rotation. The conserved class charge is

$$q(E)=0,\qquad q(R)=1,\qquad q(Z)=2.$$

Take an ordered fan of $F$ triangles, all holonomies based at its root,
with ordered product fixed to a reflection $P$ and total charge three.
Include vacancy exchange and both Hurwitz orientations on every adjacent
pair, and the twelve selected relational reactions on every consecutive
triple. Their positive generator is $A$. The adjacency-assisted angular
involutions on every adjacent pair give a second positive generator $B$.
Each channel has rate one; adding $B$ does not slow the old clocks.
Both operators are symmetric under the same uniform measure.

Exclude the invariant parallel configurations whose three reflections
are all $P$. The remaining framed state count is

$$N_{n,F}=\binom F3(n^2-1)+F(F-1)(n-1).$$

The first term consists of three reflections and identities elsewhere;
the second has one reflection and one rotation. In the latter inventory,
the rotation and the boundary product determine the reflection.

No independent-face update is being substituted for shared-link geometry:
fixed-boundary fan holonomies can be lifted by solving successive radial
links. The recorded five-face witnesses are realized on the existing
32-triangle periodic mesh. Each angular step changes one shared link,
leaves the exterior fixed, and preserves every face charge. The outer
region contributes one extra unit, making the whole-mesh charge four.
Only the designated fan channels are evolved in the finite calculation;
this is not unrestricted evolution of that whole mesh.

## Theorem 1: a single fan's linear charges are exactly blind

For $F=3$, consider functions that vanish on $RRR$ and depend only on
the six $ERZ$ layouts, with their six values summing to zero. This is a
five-dimensional invariant subspace of $A$:

- An active $RRR$ state has twelve reaction targets, each layout twice.
  Their contributions cancel. An inactive $RRR$ state stays in $RRR$.
- From an $ERZ$ state four reactions lead to $RRR$, contributing $4f$ to
  the positive generator.
- Pair transports just permute layouts: a swap involving $E$ has weight
  three, while an $R,Z$ swap has weight two. These weights are independent
  of the actual group elements and preserve the zero-sum subspace.

Every angular move preserves the layout, so $B$ annihilates this subspace.
Every centered linear face-charge observable belongs to it. Consequently,
for every odd $n\ge3$, every $\alpha\ge0$, and all $t\ge0$,

$$e^{-t(A+\alpha B)}f=e^{-tA}f.$$

This proves equality of their mean responses from any initial state and
of their equilibrium two-time correlations. It does not claim equality
of arbitrary nonlinear charge histories for composite $n$.

In layout order $ERZ,EZR,REZ,RZE,ZER,ZRE$, the restriction is the
following matrix, acting on zero-sum vectors:

$$K=\begin{pmatrix}
9&-2&-3&0&0&0\\
-2&9&0&0&-3&0\\
-3&0&10&-3&0&0\\
0&0&-3&9&0&-2\\
0&-3&0&0&10&-3\\
0&0&0&-2&-3&9
\end{pmatrix}.$$

For $f=3q_0-3$ and $d=2(q_2-q_0)$, exact rational inversion gives

$$\tau_f=\frac{\langle f,K^{-1}f\rangle}{\langle f,f\rangle}
=\frac{24009}{184184},\qquad
\tau_d=\frac{71}{483}.$$

These normalized integrated correlation times are independent of fiber
size and angular rate, despite the changing full gauge spectrum.

## Theorem 2: prime fibers hide angular motion from entire charge histories

There is a stronger obstruction for prime $n=p$ at total charge three,
on a fan of any length. It applies to every function of the charge field,
not just its linear moments.

Set $P=r_0$. The abstract group automorphisms

$$\phi_u(t_k)=t_{uk},\qquad \phi_u(r_a)=r_{ua},
\qquad u\in(\mathbb Z/n\mathbb Z)^\times$$

fix $P$. They commute with the old rules: these use multiplication,
inversion, equality, identity, and reflection/rotation predicates, all
preserved by $\phi_u$. Thus $A$ preserves the space $\mathcal H$ of
functions invariant under all these automorphisms. Charge functions lie
in $\mathcal H$.

An angular move at charge three can only act on the unique $R,Z$ pair;
every other face is identity. For prime $p$, any two nonzero rotation
labels $k,k'$ are related by the unit $u=k'k^{-1}$. Their corresponding
reflections are carried along because the product stays $P$. Therefore
each angular move stays inside an orbit of these automorphisms, and

$$B|_{\mathcal H}=0,\qquad A\mathcal H\subseteq\mathcal H.$$

The two semigroups agree on $\mathcal H$. Multiplication by any charge
observable also preserves $\mathcal H$. Applying the Markov property
successively proves equality of **all finite-time joint charge
distributions**, with angular moves on or off, for any common initial
state or distribution. This is an analytic result, not an extrapolation
from a finite list of vanishing derivatives.

These extra automorphisms are not generally geometric symmetries of the
cycle fiber. They are an accidental symmetry of the old reaction and
transport laws as seen through charge. Angular dynamics can break that
symmetry in its internal state without making the breaking observable
in this charge-three sector.

For composite $n$, units are transitive only within a fixed
$\gcd(k,n)$, equivalently a fixed rotation order. Any angular edge that
keeps that order remains invisible to $\mathcal H$; order-changing edges
can couple to charge. Thus the coupling below specifically probes
arithmetic order classes, not merely a change in angular distance.

The same argument gives a useful **flat-boundary corollary**: for prime
fibers at total charge four and boundary product $E$, angular moves are
also invisible to entire charge histories. The inventories are $RRRR$,
$RRZ$, and $ZZ$, with identities elsewhere. Only $RRZ$ admits an angular
move. Its adjacent $R,Z$ pair has product equal to the remaining
reflection $S$, since the total product is identity. A group automorphism
centered at $S$ fixes $S$ and rescales any nonzero rotation to any other;
the paired reflection follows from the fixed pair product. Thus every
active angular edge again stays inside an abstract-automorphism orbit.
Merely adding a fourth charge does not remove the obstruction in a
flat-boundary region.

## Overlap breaks the linear blindness for composite fibers

For a self-adjoint positive pair $A,B$, let

$$r=\min\{j\ge0:\langle A^jf,BA^jf\rangle>0\}.$$

Positivity makes all earlier zero energies equivalent to $BA^jf=0$.
In the power-series expansion of
$C_\alpha(t)=\langle f,e^{-t(A+\alpha B)}f\rangle$, any word containing
$B$ must have at least $r$ copies of $A$ to its left and right to survive.
At the first possible order the only surviving word is $A^rBA^r$.
Hence, for $\alpha>0$,

$$\boxed{C_\alpha(t)-C_0(t)
=-\frac{\alpha t^{2r+1}}{(2r+1)!}
\langle A^rf,BA^rf\rangle+O(t^{2r+2}).}$$

This provides a finite, exact criterion for whether an internal primitive
can influence a chosen observable. The inner product here is uniform
with normalization $1/N$; $C$ is not normalized by the variance.

Integer Krylov powers and integer sums of squared differences give:

| Fiber | Fan faces | Framed states | First response depth $r$ | First changed correlation derivative |
| --- | ---: | ---: | ---: | ---: |
| $C_9$ | 4 | 416 | 6 | 13 |
| $C_9$ | 5 | 960 | 6 | 13 |
| $C_{27}$ | 5 | 7,800 | 6 | 13 |
| $C_{81}$ | 5 | 67,200 | 6 | 13 |
| $C_{25}$ | 4 | 2,784 | 7 | 15 |
| $C_{49}$ | 4 | 10,176 | 8 | 17 |

Both $f=Fq_0-3$ and $d=\sum_i(2i-F+1)q_i$ have the stated depth.
These are exact finite calculations, not a proven formula for other
fiber sizes or charges. Prime $C_5,C_7$ five-face calculations agree with
the analytic blindness theorem.

For example, the five-face $C_9$ endpoint coefficient is

$$\langle A^6f,BA^6f\rangle=327680/3.$$

A recorded pair has the same charge field $(0,1,2,0,0)$ but
$A^6f=-36318200$ versus $-36315640$. It is connected by one actual
angular link move. Its mean endpoint responses under $A$ first differ
at order six; their equilibrium angular-on/off correlation first differs
at order thirteen. Neither number is a lower bound on when individual
sample trajectories might separate.

The integrated effect is small. On five faces, angular moves reduce
the endpoint integrated correlation time by approximately 2.59, 0.376,
and 0.0452 parts per million for $n=9,27,81$, respectively. Dipole
reductions are 3.72, 0.530, and 0.0631 ppm. Numerical Poisson residuals
are below $10^{-12}$; the visibility depths themselves use exact integers.

The direction of the integrated change follows from the variational
formula for $\langle f,L^+f\rangle$: adding the positive form $\alpha B$
can only reduce its supremum. Here $f$ is orthogonal to every old
stationary component. This does **not** imply pointwise ordering of the
two autocorrelation functions at every time.

## Consequence for the research direction

Making the internal state space connected, refining the fiber, or finding
a vanishing internal gap is insufficient to establish a spatial force.
At low charge, the present rules can hide even genuinely changing
angular states exactly, or expose them only through weak order-sensitive
responses. For prime fibers the next sectors not excluded by these
arguments are charge five with reflection boundary, or charge six with
identity boundary. Their angular-active inventories can contain extra
holonomies that cannot all be fixed by the same rescaling. A nonidentity
rotation boundary is another possible anchor, but is a different boundary
condition and must be treated explicitly. Do their charge-observable
algebras still hide angular information, or does a relative angular
quantity survive? This is the next mathematical target, without adding
a desired potential, particle shape, or continuum equation.

The apparatus is [fiber_charge_coupling.py](../tools/fiber_charge_coupling.py),
with exact coefficients, numerical integrals, and raw shared-link
witnesses in [the calculation output](../data/fiber-charge-coupling.json).
These results concern supplied geometry and classical stochastic clocks;
no quantum amplitude, physical time scale, emergent dimension, or molecule
has been derived.
