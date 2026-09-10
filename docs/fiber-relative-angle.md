# Relative holonomy becomes observable as charge noise

Follow-up: [angular control of conserved-charge transport](fiber-bulk-transport.md)
derives a rate-comparison theorem, proves that the winding-current response
has a nonzero bulk short-time coefficient, and evolves the full bank on
32-, 128-, and 288-face periodic meshes.

Independent holonomies remove the [low-charge blindness](fiber-charge-coupling.md).
There is an explicit, gauge-invariant angular move that switches a neighboring
reaction from twelve active channels to none. It changes no face charge at
that instant, but changes the instantaneous charge fluctuation covariance.
Subsequent mean charge responses then separate.

The key distinction is exact: **for the existing equal-rate bank, relative
holonomy controls charge noise, not instantaneous mean charge drift.** This
holds beyond the finite examples below. It is a microscopic stochastic
coupling, not a derived Coulomb force, quantum phase, or molecular bond.

The [earlier triangle-fiber current calculation](triangle-charge-current.md)
already established the drift/noise distinction for that fiber. The new
work here extends its formulas to the odd-cycle bank and identifies how
the newer adjacency-assisted angular primitive reaches that noise channel:
an explicit relative-holonomy switch, minimal prime-fiber sectors, exact
response orders, and full-bank reaction-time changes. The general stochastic
moment identities are not being claimed as new mathematics.

## An explicit relative-holonomy switch

Use the group derived from the supplied odd cycle $C_n$. Label its actual
permutations $t_k(v)=v+k$ and $r_a(v)=a-v$, with indices modulo $n$.
All face loops are based at the same root. For every odd $n\ge3$,

$$w_1=(r_1,r_0,r_0,t_1),\qquad
w_2=(r_1,r_0,r_1,t_2).$$

Both have ordered product $r_0$ and charge field $(1,1,1,2)$.
The first existing angular channel on the last pair takes $w_1$ to $w_2$:
it moves the rotation by one fiber edge and adjusts the paired reflection
to preserve their product. The first three entries of $w_1$ have exactly
one equal adjacent reflection pair. All twelve relational reactions are
enabled. Those of $w_2$ have neither adjacent pair equal; none is enabled.

This is not a special rotation-order effect: it works for prime cycles,
including $C_3$, where all nonidentity rotations have the same order.

More generally, consider

$$w_k=(r_1,r_0,r_{k-1},t_k),\qquad 1\le k\le n-1.$$

The anchor $h=r_1r_0=t_1$ and the mobile rotation $z=t_k$ obey $z=h^k$.
The relative exponent $k$ is unchanged by simultaneous conjugation, or
even by any abstract group automorphism applied to both $h$ and $z$.
The first triple reacts precisely when $z=h$, equivalently $k=1$.
Changing $z$ while keeping the two anchor reflections fixed is therefore
observable in a way that was impossible with only one reflection and one
rotation. These seeds specify initial states, not frozen anchors: all
the usual local channels run afterward.

Append $r_0$ to both four-tuples and their boundary product becomes identity.
This gives the same switch in a five-face region with charge six and no
boundary curvature. The saved shared-link realization changes one radial
link, fixes the exterior, and leaves the whole mesh charge at six.
The reflection-boundary four-face realization has charge five inside
and one unit outside. Neither construction assigns independent fake face
values in place of a shared connection.

## Derive the drift and noise separately

Let $G$ be the backward Markov generator, with each existing channel at
rate one. For a pair with charges $q_i,q_j$, vacancy exchange and the two
Hurwitz orientations give

$$b_i^{\rm pair}=c_{ij}(q_j-q_i),\qquad b_j^{\rm pair}=-b_i^{\rm pair},
\qquad c_{ij}=2+\mathbf1_{\{q_i=0\text{ or }q_j=0\}}.$$

For a three-face reaction support, the entire twelve-channel contribution
to $b_i=Gq_i$ is

$$b_i^{\rm reaction}=4\,\mathbf1_{ERZ}(1-q_i),$$

where $\mathbf1_{ERZ}$ means the charges are a permutation of $(0,1,2)$.
An active $RRR$ state has twelve targets, each of the six charge layouts
twice. Their mean charge at each position is one, so its drift is zero.
An $ERZ$ state has four targets with charges $(1,1,1)$. All other cases
have no reaction contribution. Angular moves preserve every charge and
contribute no direct drift.

Thus $b(q)$ is a function of charges alone. Summing local contributions
extends this statement to overlapping patches on the base graph. It
uses the equal-rate twelve-channel bank; it is not a theorem about
arbitrary asymmetrically weighted reaction laws.

Now calculate the infinitesimal covariance

$$\Gamma_{ij}(w)=\sum_{w'}c(w,w')
[q_i(w')-q_i(w)][q_j(w')-q_j(w)].$$

It is the conditional charge increment covariance per unit time; subtracting
the drift product changes only terms of order $dt^2$. A pair contributes

$$\Gamma^{\rm pair}=c_{ij}(q_j-q_i)^2
(e_i-e_j)(e_i-e_j)^\mathsf T.$$

On a triple, write $a(w)=1$ if all three entries are reflections with
exactly one adjacent equality, and zero otherwise. The reaction part is

$$\boxed{\Gamma^{\rm reaction}
=4a(w)(3I-\mathbf1\mathbf1^\mathsf T)
+4\mathbf1_{ERZ}(\mathbf1-q)(\mathbf1-q)^\mathsf T.}$$

Unlike the drift, this depends on the actual holonomies through $a$.
Angular moves themselves still contribute zero charge covariance, but
can change $a$ for the next reaction. There is no externally injected
Gaussian noise in this calculation: the covariance comes from the exact
microscopic jump increments.

For $w_1,w_2$, the common drift is $(0,0,2,-2)$, while

$$\Gamma(w_1)-\Gamma(w_2)=
\begin{pmatrix}
8&-4&-4&0\\
-4&8&-4&0\\
-4&-4&8&0\\
0&0&0&0
\end{pmatrix}.$$

Let $Z$ count rotation faces. Since
$Z=(\sum_iq_i^2-Q)/2$, it directly measures this fluctuation channel.
More generally, the exact full-bank identity is

$$GZ=\sum_{\text{triples}}(12a-4\mathbf1_{ERZ}).$$

Consequently $GZ(w_1)=12$ and $GZ(w_2)=0$. The angular step controls the
rate at which the conserved charge is reorganized into reflections and
rotations, without inserting an angle-dependent reaction rate.

## Noise feeds back into later mean motion

The charge-only drift is nonlinear in the charge field, so equal initial
drifts need not produce equal later mean responses. For the four-face
witness, direct summation of the unchanged channels gives

$$G^2q(w_1)-G^2q(w_2)=(0,8,16,-24).$$

One can obtain this without group enumeration: sum the charge-only drift
over the six reaction target layouts, each twice. All other channel
contributions cancel between the two states. Hence, under the old bank,

$$\mathbb E_{w_1}[q(t)]-\mathbb E_{w_2}[q(t)]
=\tfrac12t^2(0,8,16,-24)+O(t^3).$$

The same leading difference holds when angular channels are included,
because they annihilate both $q$ and its charge-only drift. These are
mean responses from two initial states linked by one angular move, not
individual deterministic trajectories.

For equilibrium comparisons write $A=-G_{\rm old}$ and let $B$ be the
positive angular generator. For every linear charge observable $f$,
$Bf=BAf=0$. The witnesses above, and their reflected ordering, show that
$BA^2f$ need not vanish. Thus the first angular-on/off autocorrelation
change can occur at order five, rather than the order thirteen found in
the lower-charge composite examples. Rotation count is visible one
level earlier and changes at order three.

For example, on the $C_5$ five-face identity-boundary component, with
$f=5q_0-6$ and uniform inner product on its 4,900 framed states,

$$\begin{aligned}
C_{Z,\alpha}(t)-C_{Z,0}(t)&=-\alpha\frac{144}{49}t^3+O(t^4),\\
C_{f,\alpha}(t)-C_{f,0}(t)&=-\alpha\frac{185}{49}t^5+O(t^6).
\end{aligned}$$

Here $C_{v,\alpha}(t)=\langle v,e^{-t(A+\alpha B)}v\rangle$; centering $Z$
does not change the displayed difference. The coefficients follow from
exact integer sums, not small-time curve fits.

For the identity-boundary endpoint response, an explicit alternative
pair is $(r_{-1},r_0,r_0,r_0,t_1)$ and
$(r_{-1},r_0,r_0,r_1,t_2)$. The angular step activates the second triple
without deactivating the first. The first component of the difference
of $G^2q$ is four. This supplies the endpoint witness for all odd $n\ge3$,
not just the enumerated $C_5$ example.

## Full-bank reaction times and equilibrium responses

The calculation uses all adjacent transports and all consecutive
three-face reactions, with angular channels added on every adjacent
pair. It does not freeze the anchor reflections, suppress elastic
motion, or replace evolution by a one-coordinate diffusion.

The nonzero-reflection sectors enumerated here have sizes

$$\begin{aligned}
N_{Q=5,F=4,P=R}&=4n^2(n-1)+12(n-1)^2,\\
N_{Q=6,F=5,P=E}&=5n^3(n-1)+30n(n-1)^2.
\end{aligned}$$

For the second line, the separate all-rotation sector is excluded; it
cannot exchange states with the reflection-containing sector under these
rules. The computed components for the listed prime fibers are connected
with angular motion both off and on.

For first reaction times, restrict the positive generator to the
one-rotation inventory and solve $Ku=\mathbf1$. Leaving that inventory
is the first reaction. This is an absorbing observation of the full
process, not deletion of any channel before the event.

For a uniform initial distribution on this inventory, the initial hazard
is $12(n-1)/n^2$ for four faces and $144(n-1)/(5n^2)$ for the flat-boundary
five-face region. Angular motion does not change this initial hazard.
Nevertheless the subsequent mean waiting time changes:

| Identity-boundary five-face fiber | Mean first reaction, angular off | Angular on | Reduction of integrated $Z$ correlation time | Reduction of integrated charge-dipole correlation time |
| --- | ---: | ---: | ---: | ---: |
| $C_3$ | 0.378385 | 0.356901 | 4.74% | 0.248% |
| $C_5$ | 0.716771 | 0.635521 | 10.28% | 0.548% |
| $C_7$ | 1.141271 | 0.973779 | 13.92% | 0.826% |

The correlation columns use the full equilibrium component, not the
initial distribution used for first reaction times. All times are in
per-channel rate-one units. No physical seconds or continuum transport
coefficient is assigned.

There is also genuine initial-state dependence. On the flat $C_5$ region,
adding angular motion increases the mean first reaction time from the
initially active $w_1$ seed, $0.266611\to0.284007$, but decreases it from
the inactive $w_2$ seed, $0.980755\to0.803612$. It can move away from
reaction alignment as well as help reach it. A single rate multiplier
does not describe the effect. Linear-solve relative residuals in these
calculations are below $10^{-12}$.

## Charges plus reaction flags still do not close

For $n\ge5$, $w_2$ and $w_3$ have identical charge fields and identical
activity flags for every triple: all are inactive. Consider transitions
to the same charge field with the first triple active. The total rates
are

$$w_2:\quad 2\text{ (Hurwitz)}+1\text{ (angular)},\qquad
w_3:\quad 0+0.$$

For $w_2=(r_1,r_0,r_1,t_2)$, Hurwitz on the first pair gives
$(r_2,r_1,r_1,t_2)$, and inverse Hurwitz on the middle pair gives
$(r_1,r_1,r_2,t_2)$. The angular step gives $w_1$. For $w_3$ none of
the charge-preserving channels makes either adjacent reflection pair
equal. Appending $r_0$ preserves this counterexample with identity boundary.

Therefore a Markov model retaining only charges and binary reaction flags
is not exact. Relative holonomy affects the evolution of the flags
themselves. Projecting it away requires a history-dependent description,
or a justified limiting approximation, rather than an assumed independent
reaction clock.

## What this advances, and what it does not

For prime fibers, these are minimal charge/support combinations within
the stated boundary sectors: lower odd charges with reflection boundary
are blind; a three-face charge-five region has no applicable reaction.
Lower even charges with identity boundary are blind, and a charge-six
region of at most four faces also has no applicable reaction. The first
observable switches occur at charge five on four faces, or charge six
on five faces with identity boundary. This minimality does not apply
to the lower-charge composite-fiber arithmetic effect.

We now have an explicit microscopic route from a gauge-invariant relative
holonomy to charge fluctuations, reaction timing, and subsequent mean
motion, with no fitted force law. The mechanism is conditional activation
in a classical reversible stochastic system. It does not establish a
physical particle, energy, attraction, quantum interference, or emergent
spacetime geometry. Its persistence under spatial growth, fiber refinement,
and changes of initial density remains to be determined before assigning
it a continuum physical interpretation.

The next mathematical problem is to retain the relative-holonomy variables
needed by the exact current/noise equations and determine their memory
and scaling on an expanding base region. The covariance formula above
provides a concrete observable to follow, rather than identifying slow
internal modes with particles by appearance.

Apparatus: [fiber_relative_angle.py](../tools/fiber_relative_angle.py).
Results, exact moments, and actual link states:
[fiber-relative-angle.json](../data/fiber-relative-angle.json).
