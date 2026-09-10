# A sharp interaction-arity bound for a cycle-derived spectral observable

For a cycle fiber of odd prime size $p$, consider the quadratic displacement
of its lowest nonconstant Laplacian mode. If an update is a bijection of
ordered loop tuples, is gauge covariant, preserves their ordered product,
and exactly conserves the sum of this observable, then changing any
curvature-class population requires at least

$$\boxed{2p+1\text{ loop variables}.}$$

The bound is sharp. We construct an orientation-covariant involution on
$2p+1$ variables attaining it. This gives arities $7,11,15,23,35$ for
$p=3,5,7,11,17$, with explicit finite certificates.

This is a conditional mathematical result, **not a new energy law installed
in the simulation**. The choice to conserve this spectral observable is
not implied by the original rule census. It probes whether one natural
smooth-angle candidate is compatible with fixed-size reactive primitives.

## Derive the observable from the fiber adjacency

Let $L$ be the graph Laplacian of $C_p$. Writing
$\zeta=e^{2\pi i/p}$, the functions $j\mapsto\zeta^{\pm j}$ span its
lowest nonconstant real eigenspace $V$, of dimension two. This follows
directly from

$$L\zeta^j=(2-\zeta-\zeta^{-1})\zeta^j.$$

Every graph automorphism commutes with $L$ and acts orthogonally on $V$.
Let $U(g)$ be that action. A basis-independent displacement observable is

$$\varepsilon(g)=\frac12\|U(g)-I\|_F^2=2-\operatorname{tr}U(g).$$

The adjacency-derived automorphisms give

$$\varepsilon(e)=0,\qquad\varepsilon(r_a)=2,\qquad
\varepsilon(t_k)=2-\zeta^k-\zeta^{-k}
=2-2\cos(2\pi k/p).$$

This has quadratic small-angle behavior. It is a real observable of the
derived fiber action, not a complex amplitude, a continuum group supplied
to the dynamics, or a uniquely selected physical Hamiltonian. Selecting
the lowest mode rather than a different graph observable remains a
modeling choice whose consequences are being examined.

For the pentagon, its values on $(R,Z_1,Z_2)$ are

$$\left(2,\frac{5-\sqrt5}{2},\frac{5+\sqrt5}{2}\right).$$

None of the eight reactive triple-rule equation families in the
[complete pentagon census](cycle-charge-compatibility.md) conserves this
vector. The calculation uses exact coefficients of $1$ and $\sqrt5$,
not numerical tolerances.

Even at $p=3$, these spectral values are $(2,3)$, not the existing
triangle bank's derived charge $(1,2)$. The theorem therefore does not
contradict the active three-face reactions, which conserve a different
observable.

## Exact spectral conservation has arithmetic consequences

There are $m=(p-1)/2$ nonidentity rotation classes. Suppose a transition
changes the reflection population by $x$ and the rotation populations by
$y_1,\ldots,y_m$. Exact spectral conservation gives

$$2x+2\sum_k y_k-\sum_k y_k(\zeta^k+\zeta^{-k})=0.$$

The associated polynomial has degree at most $p-1$. For prime $p$, the
minimal polynomial of $\zeta$ is
$\Phi_p(z)=1+z+\cdots+z^{p-1}$; its irreducibility follows by applying
Eisenstein's criterion to $\Phi_p(z+1)$. Thus every coefficient of our
vanishing polynomial must be equal. In particular all $y_k$ are equal,
say $y_k=t$, and

$$2x+pt=0.$$

Since $p$ is odd and the population changes are integers, write
$x=-p\ell$, $t=2\ell$. Ordered-product preservation adds another
constraint: the determinant of the low-mode action is $-1$ for a
reflection and $+1$ for a rotation, so the reflection-count change is
even. Hence $\ell$ is even, giving

$$\boxed{\Delta N_R=-2p h,\qquad
\Delta N_{Z_k}=4h\ (1\leq k\leq m),\qquad
\Delta N_e=2h,\quad h\in\mathbb Z.}$$

Every nonzero population change consumes at least $2p$ reflections in
one direction. This proves an arity bound of $2p$ before using gauge
covariance and invertibility.

It also exposes extra exact integer invariants. Every such update
conserves

$$N_{Z_k}-N_{Z_1}\quad(k=2,\ldots,m),\qquad
2N_R+pN_{Z_1}.$$

Thus conserving one algebraic-valued spectral sum enforces $m$ independent
integer count combinations. For a pentagon, for example, an energy-neutral
conversion must change both rotation populations by the same amount;
it cannot freely exchange one rotation class for another.

## Gauge covariance strengthens the bound to $2p+1$

For a conjugation-equivariant bijection $T$, the stabilizer of a tuple is
preserved exactly:

$$\operatorname{Stab}(Tx)=\operatorname{Stab}(x).$$

Equivariance gives one inclusion, and the equivariant inverse gives the
other. This is the formal version of not discarding frame information
when an interaction changes curvature classes.

At arity $2p$, a nonzero conversion would have one endpoint consisting
entirely of $2p$ reflections. Its other endpoint has four members of
each rotation class and two identities, with no reflections.

The rotation-only endpoint has stabilizer $C_p$. A reflection-containing
tuple in the odd dihedral group has stabilizer either the order-two
centralizer of one reflection, or the identity if two distinct
reflections occur. Neither equals $C_p$. An equivariant bijection
cannot connect these endpoints. Therefore

$$\operatorname{arity}\geq2p+1.$$

This proof applies to arbitrary equivariant bijections with the stated
conservation properties, not only to the minimal-involution census.

## A matching sparse reaction

Use the following ordered source tuple of length $2p+1$:

$$x=(\underbrace{r_0,\ldots,r_0}_{2p-1},r_1,r_1).$$

For the target, retain $r_0$, then for every $k=1,\ldots,m$ insert two
inverse rotation pairs, and finally two identities:

$$y=(r_0,
t_1,t_{-1},t_1,t_{-1},\ldots,
t_m,t_{-m},t_m,t_{-m},e,e).$$

Both ordered products are $r_0$. Both tuples generate the full derived
dihedral group, and each has trivial conjugation stabilizer. Their
spectral sums are both $4p+2$, since
$\sum_{k=1}^m\varepsilon(t_k)=p$.

Close the exchange $x\leftrightarrow y$ under simultaneous conjugation
and reversed orientation. Each endpoint has an orbit of size $4p$ under
these operations. For the source, reversal changes the positions of its
two unequal-label blocks; for the target, it moves the two identities
from the end to the beginning. Neither reversed tuple is a conjugate of
the original. Source and target orbits are disjoint because their class
populations differ.

The resulting **$4p$ disjoint transpositions** define an involution:
swap those pairs and fix every other tuple. It moves only $8p$ states
out of the ambient $(2p)^{2p+1}$ possibilities. It fixes the flat tuple,
preserves ordered product, is gauge and orientation covariant, and
conserves the spectral sum exactly. This attains the lower bound.

The sparse construction is checked algebraically for
$p=3,5,7,11,17$. All spectral identities are reduced modulo $\Phi_p$;
there is no floating-point acceptance window. These higher-arity rules
are **not** installed in the shared-mesh engine, and the construction
does not establish a physical scheduling law or an emergent locality
scale for them. Existence of a sparse reaction also does not establish an
appreciable activation rate or a useful continuum dynamics.

## Consequence for a continuum-oriented refinement

Along odd-prime cycle sizes, any fixed arity eventually becomes too
small to change curvature-class populations under the stated exact
spectral conservation. Transport and class-preserving rearrangements
can continue, but the number of exact population constraints grows.
An overall rescaling of $\varepsilon$ does not change this obstruction.

Therefore refining the fiber while retaining a fixed-size reversible
reaction bank and exact additive spectral conservation cannot by itself
give unrestricted reactive angle exchange. This is a precise conditional
obstacle, not a proof against every continuum limit or every cycle-size
sequence. The [prime-power analysis](cycle-prime-power-reactions.md) now
shows that for $n=p^a$, $a\ge2$, the sharp minimum arity drops to $p$,
independent of $a$, while every population-changing reaction has spectral
sum at least $2p$. Thus the prime-only arity obstruction cannot be
extrapolated to arbitrary refinement sequences.

Possible research directions include larger-support conserved observables,
additional dynamically derived energy-carrying variables, or a different
microscopic selection principle. None has been assumed here. Simply
loosening an energy tolerance or inserting a desired field Hamiltonian
would not establish their physical derivation.

Finite gauge digitization has other known obstructions, such as the
weak-coupling freezing studied by Hartung et al. in
[Digitising SU(2) Gauge Fields and the Freezing Transition](https://arxiv.org/abs/2201.09625).
That Monte Carlo/digitization issue is not the theorem above: our bound
comes from exact additive conservation, cyclotomic arithmetic, and
stabilizer preservation in local tuple updates. We do not claim this
calculation resolves the digitization problem or establishes a novel
general theorem about gauge theory.

Calculation: [cycle_spectral_reactions.py](../tools/cycle_spectral_reactions.py).
Sparse exact witnesses are in
[cycle-spectral-reactions.json](../data/cycle-spectral-reactions.json).
