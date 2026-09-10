# An explicit reaction bank and its refinement-dependent kinetics

All twelve existing triangle reaction tables can be reconstructed from
relational matching and short transport formulas. The construction extends
to every odd cycle fiber, without running a positive-charge census at the
new fiber size. It preserves the derived reflection-length charge and
commutes with the natural injective embeddings of these fiber groups.

There are two different refinement results:

- **Same embedded state:** the entire shared-link trajectory is preserved,
  including every reaction and idle step.
- **Uniformly populated larger fiber:** stationary reaction activity falls
  as $1/n$ for a cycle with $n$ vertices.

Both are consequences of the same microscopic law. The difference is
whether new internal directions are occupied in the initial ensemble.
The experiments below use actual shared links and reproduce the predicted
rate dependence. They establish classical constrained reaction kinetics,
not physical energy, continuum gauge dynamics, matter, or molecules.

## Twelve relational rules, not twelve opaque tables

The fiber graph determines $G=\operatorname{Aut}(C_n)$, for odd $n\ge3$.
Write $E$ for identity, $R$ for nonidentity involutions (reflections), and
$Z$ for other nonidentity elements (rotations). These types are determined
by multiplication and equality, not supplied charge values.

An active all-reflection input has one of the forms

$$(r,r,s),\qquad(s,r,r),\qquad r\ne s.$$

Choose a layout $\pi$, one of the six permutations of $(E,R,Z)$, and
a sign $\eta\in\{+,-\}$. This gives twelve rules. On the left-pair input
$(r,r,s)$, define

$$z=\begin{cases}rs,&\eta=+,\\sr,&\eta=-.\end{cases}$$

Its boundary product is $s$. To place a reflection $h$ and rotation $z$
in the chosen layout while preserving that product, set

$$h=\begin{cases}
z^{-1}s,&Z\text{ precedes }R,\\
sz^{-1},&R\text{ precedes }Z.
\end{cases}$$

Output $(e,h,z)$ in layout $\pi$. The right-pair input is handled by
reversed inversion $J(A,B,C)=(C^{-1},B^{-1},A^{-1})$:

$$T(s,r,r)=J\,T(r,r,s).$$

There is an explicit inverse, so no root choice or search is needed.
For an output $w$ of layout $\pi$, let $s$ be its ordered product and
$z$ its rotation entry. Recover

$$r=\begin{cases}zs,&\eta=+,\\sz,&\eta=-,
\end{cases}\qquad w\longmapsto(r,r,s).$$

For layout $\pi$ reversed, conjugate this prescription by $J$.
Fix every remaining tuple. The two target layouts are disjoint because
their endpoint types differ. Nontrivial rotations are not involutions
when $n$ is odd, so recognition is unambiguous.

These formulas prove that each rule is an involution, fixes the flat
tuple, preserves ordered product, and commutes with simultaneous
conjugation and reversed inversion. Both source and target generate the
same subgroup: all target entries are words in $r,s$, while $s$ is the
target product and $r$ is recovered by the displayed formula.

Each rule moves $4n(n-1)$ tuples. Every reflection input with exactly
one equal adjacent pair enables all twelve rules. Every tuple containing
one $E$, one $R$, and one $Z$ enables four. All other inputs enable none.

The complete maps on $C_3$ match the existing twelve tables exactly,
including their inverse branches, not just selected trajectories.
The saved data records the layout-to-existing-rule correspondence.
The $(E,R,Z),+$ case is the previously reconstructed rule 27.

For larger fibers this is an explicit extension of the triangle bank,
**not** a claim to exhaust all equivariant finite laws or reproduce the
old positive-charge selection criterion. Other triple configurations can
support additional finite rules at larger $n$.

## Conservation derived after specifying the rules

Normalize a class charge by $q(e)=0$, and let its reflection value be $c$.
The exchange $(r,r,s)\leftrightarrow(e,r,rs)$ implies $q(rs)=2c$.
Every nonidentity rotation is a product of two distinct reflections, so
the full common additive class-charge space is exactly

$$q(e)=0,\qquad q(R)=c,\qquad q(Z)=2c.$$

These charges are sufficient because every reaction exchanges the types
$RRR$ and $ERZ$. With $c=1$, a forward event changes populations by

$$(\Delta N_E,\Delta N_R,\Delta N_Z)=(1,-2,1).$$

The conserved quantity is reflection length, not the lowest-mode spectral
sum. Its nonzero rotation value remains discontinuous as rotation angles
approach zero under refinement.

The construction was reconstructed from a bank historically selected by
a positive-charge criterion. A charge-free formula for those rules does
not retroactively make their original selection an unbiased derivation
from adjacency. The matching law remains a modeling choice; what is new
here is its explicit structure, extension, and consequences.

## Exact drift, before any statistical closure

For the generic shared-fan experiment, choose a fan uniformly and one of
fifteen channels uniformly: vacancy exchange, Hurwitz, inverse Hurwitz,
or one of the twelve reactions. The three transport channels act on the
first two entries in the fan's ordered based tuple. They preserve class
populations. This support convention is explicit; it is not claimed to
duplicate the compiled triangle engine's separate pair-support scheduler.

Let $X_p$ indicate a fan containing three reflections with exactly one
equal adjacent pair, and $Y_p$ indicate a fan with types $ERZ$ in any
order. For $M$ possible fans the exact one-proposal drift is

$$\boxed{\mathcal L N_Z=
\frac1M\sum_p\left(\frac45 X_p-\frac4{15}Y_p\right).}$$

This follows from twelve enabled forward channels versus four reverse
channels. It is a raw-state identity, not a mean-field equation. In
particular the forward term depends on actual reflection equalities, not
just three scalar face charges. No independence assumption is needed.

## Stationary raw-link activity and the $1/n$ bottleneck

Initialize every raw link independently and uniformly in $G$, on the
supplied periodic triangular mesh. Each shared-link update is bijective:
the exterior is fixed and the internal spokes are reconstructed uniquely
from the reversible target tuple. Therefore the uniform raw-link measure
is stationary for each channel and their state-independent mixture.

This is a mixture of conserved sectors, not an assertion that one
trajectory explores the whole ensemble. No equilibration is assumed.

Within any three-face fan, condition on all links except its three rim
edges. Those rim links are independent uniform group elements, each
entering one of the based loops. Hence the three based holonomies are
independent uniform elements of $G$, even on the closed mesh. Thus

$$p_E=\frac1{2n},\qquad p_R=\frac12,\qquad
p_Z=\frac{n-1}{2n}.$$

Given three independent reflections, exactly one adjacent equality has
probability $2(n-1)/n^2$. Consequently

$$\Pr(X)=\frac{n-1}{4n^2},\qquad
\Pr(Y)=\frac{3(n-1)}{4n^2}.$$

The expected forward and reverse rates agree, as required by stationarity.
Counting all proposals, including idle ones, the total reaction activity is

$$\boxed{\lambda_n=\frac{2(n-1)}{5n^2}\sim\frac{2}{5n}.}$$

Conditional on having selected a reaction channel, the activity is
$(n-1)/(2n^2)$. The slowdown comes from two matching probabilities:
coincident reflection axes become rare, and the identity occupies only
one of the increasingly many group elements. It is not a floating-point
acceptance threshold or a fitted kinetic parameter.

The scientific run used 64 independent trajectories per fiber, each with
4,000 proposals on a side-four mesh (32 faces, 48 links). Across five
fibers this is 1,280,000 proposals. Each event was evaluated from based
loop transports and realized by the existing shared-link reconstruction.
Errors below are standard errors across independent trajectories, not
across correlated events or spatial translations.

| Fiber size $n$ | Measured reaction activity | Exact stationary value |
|---|---:|---:|
| 3 | $0.08743\pm0.00167$ | $0.08889$ |
| 5 | $0.06811\pm0.00226$ | $0.06400$ |
| 9 | $0.04152\pm0.00192$ | $0.03951$ |
| 27 | $0.01332\pm0.00082$ | $0.01427$ |
| 81 | $0.00552\pm0.00055$ | $0.00488$ |

Each mean is within two estimated standard errors of its prediction.
This agreement measures stationary activity, not a relaxation exponent,
an ergodic theorem, or a universality class. The mean charge of this
particular ensemble is $3/2-1/n$, so these are not fixed-density runs.
Raw reaction witnesses and independent trajectory counts are saved;
whole-mesh charge and the event-count population balance hold exactly.

## A fixed-density reference, with its assumptions exposed

For independent based loops, the charge-weighted reference

$$\mu_z(g)=\frac{z^{q(g)}}{D},\qquad D=1+nz+(n-1)z^2$$

is invariant under every tuple gate, because each gate is a bijection
preserving the summed charge. Here $z>0$ is a reference fugacity, not a
coupling or temperature installed in the dynamics. The corresponding
class probabilities satisfy

$$p_E=\frac1D,\quad p_R=\frac{nz}D,\quad
p_Z=\frac{(n-1)z^2}D,\qquad
\boxed{\frac{p_Ep_Z}{p_R^2}=\frac{n-1}{n^2}.}$$

Independent reflections are uniform over their $n$ axes. The forward and
reverse rates per reaction-channel proposal are therefore

$$r_+=\frac{2(n-1)}{n^2}p_R^3,
\qquad r_-=2p_Ep_Rp_Z,$$

which coincide in this reference. This is an equilibrium counting
identity, not a closed kinetic equation for arbitrary graph states.
Nonuniform charge-weighted measures on a finite closed mesh can have
boundary and global correlations, so the product formulas must not be
silently substituted for those measures. The $z=1$ shared-link result
above has its own exact independence argument.

At fixed reference mean charge $\rho=p_R+2p_Z\in(0,2)$, as $n\to\infty$,

$$
(p_E,p_R,p_Z)\longrightarrow
\begin{cases}
(1-\rho,\rho,0),&\rho<1,\\
(0,2-\rho,\rho-1),&\rho>1.
\end{cases}$$

At $\rho=1$, $z=1/\sqrt{n-1}$ and
$p_E=p_Z\sim n^{-1/2}$, while $p_R\to1$. In particular the scaled total
reference activity satisfies

$$n(r_++r_-)\longrightarrow4\min(\rho,2-\rho)^3.$$

Below unit charge density, the rotation population is suppressed as
$p_Z\sim\rho^2/[n(1-\rho)]$. Above it, vacancies instead satisfy
$p_E\sim(2-\rho)^2/[n(\rho-1)]$. These are reference-ensemble scaling
laws, not a demonstrated phase transition or a simulated low-energy limit.

## Coherent refinement does not slow the same state

For odd $n,d$, the embedding

$$\phi:\operatorname{Aut}(C_n)\hookrightarrow\operatorname{Aut}(C_{dn}),
\qquad t_k\mapsto t_{dk},\quad r_k\mapsto r_{dk}$$

preserves and reflects the rule's equality, identity, and involution
predicates. The target formulas are group words. Therefore every one of
the twelve rules commutes with this embedding; so do the transport rules.

Shared-link reconstruction also uses only group words in the input links
and target holonomies. Embed every initial raw link and replay the same
schedule at both resolutions. Inductively,

$$\boxed{U^{\mathrm{fine}}_e(t)=\phi(U^{\mathrm{coarse}}_e(t))
\quad\text{for every edge and every proposal}.}$$

This preserves changing events, idle events, and their order exactly.
Eight coupled $C_3\to C_{27}$ trajectories, 4,000 proposals each, satisfy
this identity on every raw state. Full tuple calculations also establish
the rule correspondence for $C_3\to C_9$, $C_3\to C_{27}$, and
$C_9\to C_{27}$. Fine initial links in these runs occupy the embedded
coarse subgroup, not the uniform distribution over the full fine group.

Thus the observed $1/n$ activity law does **not** say that adding unused
fiber resolution slows an existing configuration. It describes a changed
ensemble that populates all of the larger fiber's internal directions.

## The remaining obstruction is dynamical access to new fiber structure

The holonomy-image subgroup is preserved, up to the usual change of root
frame. To see this, choose a gauge in which every link lies in that
subgroup, with spanning-tree links set to identity. Every new link is a
word in the old links on the chosen branch, so the new holonomy image is
a subgroup of the old one. Apply the inverse update to obtain the reverse
inclusion. Gauge covariance makes the conclusion independent of the
chosen representative. Explicit connectors remain part of this argument;
unrelated loop frames are not identified.

Consequently an embedded coarse state cannot spontaneously populate the
new group directions under this bank. Enlarging the fiber or accelerating
the same operations does not supply that missing dynamics. A candidate
time rescaling by $n$ may retain finite stationary activity in full-fiber
ensembles, but does not prove a kinetic limit or remove subgroup trapping.

Possible next mechanisms must make their additional structure explicit:
fiber-changing rewrites, relational rules that can grow holonomy image,
or additional microscopic variables. Their effects cannot be inferred
from the present fixed-fiber reaction bank alone.

The [adjacency-assisted angular construction](fiber-angular-transport.md)
now supplies one such extra mechanism: it preserves this bank's coarse
charge but uses the ambient fiber graph to move into new holonomy
directions. Its isolated shared-pair diffusion and reaction-mediated
waiting dynamics are derived and measured separately; it does not yet
derive a physical Hamiltonian or changing fiber geometry.

Kinetic constraints and slow relaxation are established topics; see
Hartarsky and Toninelli's [Kinetically constrained models](https://arxiv.org/abs/2412.13634).
The present result is the explicit algebra and measured refinement behavior
of this gauge-covariant shared-link rule family. No general novelty,
glassy phase, or emergent fundamental-physics claim is made.

Calculation: [cycle_relational_dynamics.py](../tools/cycle_relational_dynamics.py).
Exact correspondence, raw-link trajectories, and reference values:
[cycle-relational-dynamics.json](../data/cycle-relational-dynamics.json).
