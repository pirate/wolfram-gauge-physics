# An exact integer carrier for a rare reaction's memory

The [reaction-burst calculation](fiber-reaction-bursts.md) identified local
flatness as the source of correlated reactions. We can now follow that
correlation after it leaves a recognizable triangle or quadrilateral.
At fixed volume and fixed angular rate, the event-conditioned fine-fiber
limit closes on **raw link signs plus one primitive integer affine
constraint**. It retains the original local rules and needs no frozen
reference holonomies.

This is a limit for the response to one typical rare stationary reaction,
not a replacement for the entire finite-fiber model. It cannot describe
the arrival and interaction of multiple independent rare correlations on
times growing with fiber size. Its main result is a precise distinction
between information remaining present and information remaining locally
reactive.

## State and origin of the constraint

Use the derived affine coordinates of the actual odd-cycle automorphisms,

$$U_e(v)=s_ev+a_e\pmod n,\qquad s_e\in\{-1,1\}.$$

Let $m$ be the number of raw links. A typical marked post-reaction state
has either an identity face or one equal adjacent reflection pair, with
equal probabilities. Each imposes one independent linear relation on
the raw shifts. Its other nonzero requirements exclude sets of relative
measure $O(1/n)$.

Represent the inherited relation by

$$\boxed{(s,\ell,b),\qquad
\ell\in\mathbb Z^m\text{ primitive},\quad b\in\mathbb Z,\quad
\ell\cdot a=b.}$$

In a finite fiber this equation is read modulo $n$. In the limiting
description it is a formal integer relation: independent symbolic shifts
represent the remaining free coordinates, and an integer constant
represents an actual number of elementary cycle steps. In particular,
a constant one is not set equal to zero just because $1/n\to0$ in a
normalized angular plot. That distinction matters to exact-match rules.

Normalize by making the first nonzero entry of $\ell$ positive. Initially
$b=0$, and $\ell$ is the boundary covector of either a triangle or a
two-face region. This is a description of the generic event-conditioned
ensemble, not a claim that a numerical finite state has only one possible
algebraic relation.

For any local affine expression $v\cdot a+c$, its vanishing in this
ensemble is determined exactly:

$$v\cdot a+c\equiv0
\quad\Longleftrightarrow\quad
v=t\ell\text{ and }c+tb=0\text{ for some }t\in\mathbb Q.$$

The structurally zero expression is included with $t=0$. Every guard in
the existing rules is such a zero test, together with link/face signs:
identity holonomy, equality of two reflections, or whether an angular
step would reach identity. Thus all guard outcomes are determined by
$(s,\ell,b)$ in the generic limit.

## Updates come from the existing group words

Within a chosen branch of an update, actual raw shifts transform by an
integer affine map $a'=Ma+d$. Pair transport has an integer inverse,
because its inverse Hurwitz word is already in the bank. Its matrix is
therefore unimodular. Angular motion is a unit translation of one raw
shift, after the existing nonidentity guard has been evaluated.

For either kind of move,

$$\boxed{\ell'=M^{-T}\ell,\qquad b'=b+\ell'\cdot d.}$$

Only two internal spokes can be written. Let $W$ denote those edges and
$S$ the local read set. With covectors written as rows for this formula,

$$\ell'_W=\ell_W(M_{WW})^{-1},\qquad
\ell'_j=\ell_j-\ell'_WM_{Wj}\quad(j\in S\setminus W).$$

Everything outside the read set stays fixed. The implemented matrices
are obtained by composing the original group words on a local affine
basis and using the existing spoke reconstruction; they are not chosen
to generate a desired spreading law.

For a reaction or identity exchange, the applicable branch maps the
input constraint plane bijectively onto its target constraint plane.
The inverse rule supplies the inverse map on that plane. The target
relation is read from the actual new signs and based loops: either the
new identity face or the new equal-reflection pair. It has $b'=0$.
The number of independent inherited relations stays one.

An angular move changing $a_e$ by the signed unit $\delta$ has the
particularly transparent effect

$$\ell'=\ell,\qquad b'=b+\ell_e\delta.$$

It can move the correlation away from zero holonomy while preserving
every face charge. It cannot create an identity face: the original
forbidden-zero guard remains part of the update.

## Why this is a finite-time limit, rather than a new approximation

Fix the mesh, a finite $\alpha$, and an observation time $T$. The exact
two-direction representation of the angular generator makes the total
proposal rate $3F(15+2\alpha)$ independent of $n$. There are almost surely
finitely many proposals before $T$.

For any bounded proposal prefix, each branch is an integer affine word
in the independent initial shifts. If a guard is not identically zero on
the inherited plane, it becomes an additional nonzero integer affine
condition on those free shifts. Its probability modulo $n$ tends to zero.
For a nonzero coefficient $c$, it is bounded by $\gcd(c,n)/n\le|c|/n$;
a nonzero constant cannot vanish once $n$ exceeds its absolute value.
This argument works along all growing odd $n$, not only primes.

There are finitely many words for a bounded prefix, so a union bound
controls all accidental extra equalities. First cut off the Poisson
proposal count, then take $n\to\infty$, then remove the cutoff. This gives
finite-time coupling of the face charges, signs, and reaction history
to the integer constraint process with probability tending to one.

The inverse integer words preserve the generic one-plane description
after each event; excluding previously avoided extra equalities changes
the finite-prefix law by a vanishing amount. A finite-time statement does
not justify an infinite-time limit, a time of order $n$, or an angular
rate diverging with $n$.

## Gauge-invariant measures of where the information goes

The carrier is not a gauge-frame drawing. For an oriented canonical edge
$e:u\to v$, a translation of the vertex frames changes

$$a_e\mapsto a_e+\beta_v-s_e\beta_u.$$

Every inherited covector satisfies the exact twisted closure relation

$$\boxed{(C_s\ell)_v
=\sum_{e:\,\mathrm{head}(e)=v}\ell_e
-\sum_{e:\,\mathrm{tail}(e)=v}s_e\ell_e=0.}$$

Hence its relation does not depend on the translation frames. A general
vertex frame with sign $\tau_v$ changes $\ell_e$ to $\tau_v\ell_e$, before
the common normalization sign is chosen. Twisted closure cancels the
translation terms again. Consequently the support of $\ell$, its absolute
coefficient norms, and $|b|$ are gauge-invariant diagnostics.

This supplies two rigorously distinct reasons for inactivity:

- **Nonzero offset:** if $b\ne0$, no independent nonzero local homogeneous
  holonomy expression can vanish. There is no identity face or active
  reflection-equality reaction anywhere. Hurwitz moves preserve $b$ up to
  the harmless common sign; angular motion can change it again.
- **Spatially extended relation:** even with $b=0$, the covector may no
  longer match a triangle or an eligible two-face boundary. The relation
  is then present but not recognized by a local reaction guard.

There can also be a locally represented two-face relation with no eligible
third reflection. Inactivity alone is therefore not a sufficient measure
of spatial delocalization; the saved support, coefficients, signs, and
offset distinguish these cases.

All nonidle moves have available inverses in the constraint process.
An inactive state reached in finite time is not an absorbing state:
there is a positive-probability reverse path back to activity. Whether
the infinite-state reversible walk is recurrent or transient is a further
mathematical question, not settled by observing finite-time spreading.

## One correlation cannot form an independent branching cascade

Two distinct identity faces would impose two independent triangle
relations. An identity face and a separate equal-reflection pair would
likewise impose independent triangle and quadrilateral relations. Their
edge supports cannot represent the same primitive covector.

Therefore in the one-constraint limit:

$$N_E\in\{0,1\},\qquad
N_E=1\ \Longrightarrow\ \text{no forward }RRR\to ERZ\text{ reaction}.$$

When $N_E=0$, no backward $ERZ\to RRR$ reaction is possible. Every reaction
changes $N_E$ by one, so the reaction directions strictly alternate:

$$\boxed{\left|N_{RRR\to ERZ}(t)-N_{ERZ\to RRR}(t)\right|\le1.}$$

Multiple available sites are alternative uses of the same correlation,
not independently multiplying defects. Finite fibers can acquire extra
independent equalities and leave this regime; that is precisely the
multiple-constraint physics this limit does not replace.

## The signs define an auxiliary branched surface

There is a topological interpretation of $C_s$, not a claim that a new
physical space has been simulated. Make two copies of each base vertex.
A positive-sign link stays on its sheet; a negative-sign link changes
sheets. An even face lifts to two triangular disks. A reflection face
lifts to one disk whose boundary traverses both sheets, with a branch
point in its interior.

Let $R=N_R>0$. The resulting connected oriented surface has cell counts

$$\widetilde V=2V,\qquad\widetilde E=2E_g,\qquad
\widetilde F=2F-R.$$

The supplied base is a torus, so its Euler characteristic is zero. Thus
the auxiliary double cover has

$$\widetilde\chi=-R,\qquad\boxed{g=1+R/2.}$$

The anti-invariant edge chains under sheet exchange have twisted boundary
matrix $C_s$ after a choice of orientation signs. An odd face contributes
no anti-invariant face boundary; the even faces contribute their actual
holonomy covectors. Over $\mathbb Q$, therefore,

$$H_1^-(\widetilde\Sigma;\mathbb Q)
\cong\ker C_s\,/\operatorname{span}\{\text{even-face boundary covectors}\}.$$

For $R>0$, the sign cover is connected and $C_s$ has rank $V$. The
$F-R$ even-face boundary rows are independent: any dependence propagates
through the dual graph and vanishes at a missing odd-face row. Therefore

$$\boxed{\dim H_1^-=E_g-V-(F-R)=R.}$$

This quotient is a diagnostic, not a replacement state space for the
simulator. A nonidentity even face can carry nonzero curvature; adding its
boundary can change a carrier's actual holonomy condition and reactivity.
The dynamics retain the full integer covector, not just its homology class.

The connection between reflection braids and double-cover homology is
established mathematics; compare
[Brendle and Margalit, section 2.1](https://arxiv.org/html/1410.7416v3#S2.SS1)
and the repository's [earlier finite Burau calculation](triangle-overlap-algebra.md).
Here the construction is tied to the current raw sign field and its exact
integer carrier. No claim about a new braid representation is made.

The derived conserved charge gives a further identity:

$$Q=2F-R-2N_E
\quad\Longrightarrow\quad
\boxed{g+N_E=1+F-Q/2.}$$

This is the existing conserved charge expressed in cover language, not
an additional independent conservation law.

A forward reaction reduces the auxiliary genus by one and creates an
identity face; its inverse raises the genus by one and removes that face.
These are counts in a constructed cover. They are not creation or
destruction of handles in physical spacetime.

An identity-face carrier is an even-face boundary and hence trivial in
this homology quotient. In the reflection-pair phase of the present Palm
sector there are at least four reflection faces. The lifted two-face disk
is an annulus, while its complement has other branch points and is
connected; its core is nonseparating. The corresponding carrier has a
nonzero anti-invariant homology class. This also explains why losing a
short local representative need not erase the correlation.

## Observations with an independent finite-fiber shadow

The integer process uses 64 independent initial states per 32-face
condition, with $\alpha=0,1,8$, and 32 initial states on 128 faces with
$\alpha=1$. Starts are matched across the three 32-face conditions.
Each trajectory is followed to model time two.

Every trajectory also evolves actual $D_{1,000,000,007}$ raw links under
the same proposal schedule. The finite initialization is sampled from
the corresponding Palm constraint plane, with the required nonidentity
conditions retained. Across 1,191,840 proposals there were no detected
disagreements in local targets; signs, the modular carrier relation, and
sampled face charges also agreed. This comparison supports the apparatus;
the finite-time limit follows from the affine-word argument above.

The carrier spreads while reactions die down:

| Faces / angular rate | Mean support at $t=2$ | Median $\|\ell\|_1$ at $t=2$ | Mean reactions through $t=2$ |
| --- | ---: | ---: | ---: |
| 32 / 0 | $45.97\pm0.55$ of 48 edges | 5,883 | $5.875\pm0.840$ |
| 32 / 1 | $46.00\pm0.34$ of 48 edges | 4,597 | $3.938\pm0.509$ |
| 32 / 8 | $46.48\pm0.32$ of 48 edges | 5,474 | $2.938\pm0.388$ |
| 128 / 1 | $135.56\pm4.79$ of 192 edges | 19,946 | $3.156\pm0.514$ |

At $\alpha=0$, the offset remains exactly zero, yet all 64 samples are
reaction-inactive at $t=2$: coefficient spreading alone can hide the
initial local relation. At positive angular rates, all saved samples at
$t=2$ also have nonzero offset. No sampled reaction sequence violates the
proved alternation or the bound of one identity face.

The absence of activity at that checkpoint is not a zero-probability
theorem. With 64 independent samples, observing none still permits a
one-sided 95% upper probability of about 4.6%. The coefficient means have
large tails; medians are reported above, and the full observations remain
available. This is not a fitted Lyapunov exponent or an infinite-volume
spreading law.

The topological calculation uses exact integer matrix ranks on sixteen
saved initial/final states. For example, a 32-face state with $Q=50$
starts just after a forward reaction with $R=12$, $N_E=1$, $g=7$, and a
homologically trivial carrier. Its final state has $R=14$, $N_E=0$,
$g=8$, and a nontrivial carrier. The invariant $g+N_E=8$ and the carrier's
twisted closure hold exactly. All saved final carriers are nontrivial in
the quotient even though no local reaction is active there.

## Consequence and next mathematical problem

The observed single-collision memory is not a growing collection of
independent particles. It is one relation whose integer coefficients and
support can spread through the gauge field. Its information survives,
but local reactivity need not. The auxiliary homology detects that
survival without relying on a spatial layout of the raw graph.

This does not yet give a persistent localized object or molecular binding.
The next problems are to determine return versus escape in this reversible
integer process, and to extend the algebra to two independent relations
so their encounters can be distinguished from repeated uses of one
relation. Neither task is answered by imposing an attractive potential
or interpreting the auxiliary cover as emergent physical spacetime.

The follow-up [joint two-constraint calculation](fiber-two-constraints.md)
now establishes non-additive contact rates and identifies the minimal
angular-sensitive sector in an all-reflection reference preparation.
The [exact dilute encounter law](fiber-dilute-encounter.md) distinguishes
volume-dependent pair lifetimes and slowing decay from positional binding.

Apparatus and observations:
[integer carrier dynamics](../tools/fiber_constraint_dynamics.py),
[shared-schedule observations](../data/fiber-constraint-dynamics.json),
[exact cover-cycle ranks](../tools/fiber_constraint_topology.py),
[topological states](../data/fiber-constraint-topology.json).
