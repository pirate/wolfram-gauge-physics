# Intersection topology distinguishes constraint overlap from proximity

The [two-constraint experiment](fiber-two-constraints.md) compared a pair
of overlapping reflection relations with two separated relations. Their
initial reaction rates really are $72$ and $96$, respectively, but the
comparison changes more than distance: the inherited cycles have
different intersection rank on the auxiliary sign cover.

This supplies both a correction to the interpretation of that comparison
and a new diagnostic. Elastic motion preserves the intersection form;
reactions need not. A nonzero pairing obstructs separating the relations
into disjoint cycle representatives during elastic evolution, but does
not imply that their joint support stays small.

## Construct the form from oriented incidence, not a drawing

The raw signs define the branched double cover described in
[the one-constraint note](fiber-constraint-dynamics.md). A raw covector
$\ell$ obeying $C_s\ell=0$ lifts to an anti-invariant cycle: on the lift
of $e:u\to v$ ending at sheet $t\in\{+1,-1\}$, give the oriented edge
coefficient $t\ell_e$. That edge starts at sheet $s_et$ over $u$.

At the positive-sheet lift of a vertex, its outward dart coefficient is

$$d_{ve}(\ell)=
\begin{cases}
s_e\ell_e,&v=\operatorname{tail}(e),\\
-\ell_e,&v=\operatorname{head}(e).
\end{cases}$$

Twisted closure says $\sum_{e\ni v}d_{ve}=0$. The oriented mesh triangles
give a cyclic ordering of the incident darts. No Euclidean projection
or renderer is used. With $i<j$ referring to a cut of that cyclic order,
the lifted intersection pairing is

$$\boxed{I_s(\ell,m)=\sum_v\sum_{i<j}
\left[d_{vi}(\ell)d_{vj}(m)-d_{vi}(m)d_{vj}(\ell)\right].}$$

To obtain the formula, thicken the lifted graph to an oriented ribbon
surface and separate coincident cycle segments in the edge ribbons.
The signed crossing count in each vertex disk is half the displayed
cyclic-order wedge. The negative sheet negates both sets of dart
coefficients, so its contribution is equal to the positive sheet's.
The two halves add to the formula above. Changing the cut of the cyclic
order changes the expression by multiples of the zero dart sums, so
the result is cut-independent. Reversing the surface orientation changes
only its overall sign.

The form vanishes between every even-face boundary and every closed
cycle. It therefore descends to the anti-invariant homology quotient.
For a connected sign cover with $R>0$ reflection faces, its rank on the
twisted cycle space is $R$, with the $F-R$ even-face boundaries as its
radical. This agrees with the previously derived dimension of
$H_1^-$ and the nondegenerate surface intersection form.

If $N$ has a basis of the inherited constraint space as its rows, the
observable we use is

$$\sigma(C)=\operatorname{rank}(N J_s N^T),\qquad
I_s(\ell,m)=\ell J_sm^T.$$

Changing the constraint basis acts by congruence and leaves this rank
unchanged. Gauge changes merely relabel sheets and transport the cycles,
so the rank is gauge invariant too. Individual matrix entries depend on
the chosen integral representatives; their absolute magnitude is not
being treated as a normalized invariant of an arbitrary rational basis.

## What the update rules preserve

The Hurwitz rule exchanges two puncture holonomies by the usual braid
word. On the corresponding branched cover it induces an
orientation-preserving homeomorphism, hence preserves intersection of
transported cycles. This is the established topological interpretation
of the integral Burau action, not a newly discovered representation; see
[Brendle and Margalit, section 2.1](https://arxiv.org/html/1410.7416v3#S2.SS1).
Here the form is evaluated on the actual mesh's shared-link constraints.

There is also a direct raw-coordinate identity. If an elastic update is
$a'=Ma$, its transported covectors are $\ell'=\ell M^{-1}$. For a matrix
$Z$ whose columns span $\ker C_s$,

$$Z^T\left(M^{-1}J_{s'}M^{-T}-J_s\right)Z=0.$$

The calculation derives $M$ from the original group words. It verifies
this whole-cycle identity with exact rational arithmetic for both
Hurwitz orientations and all 128 sign assignments in the seven-edge
read patch, with the exterior signs positive. This finite certificate
is not an enumeration of all global sign fields; the general explanation
is the lifted homeomorphism.

An identity exchange agrees with the Hurwitz map on its applicable
plane. Angular steps change offsets but neither the linear constraint
space nor the signs. Consequently all nonreaction channels preserve
$\sigma(C)$.

Reactions are different: $RRR\to ERZ$ removes two branch points and
lowers the auxiliary cover genus by one. The new identity boundary is
in the radical of the new intersection form. More generally, $k$
independent identity-face boundaries in a rank-$r$ constraint space imply

$$\boxed{\sigma(C)\le2\left\lfloor\frac{r-k}{2}\right\rfloor.}$$

Thus at rank two, one identity face already forces $\sigma=0$.
The original inverse reaction can restore a nonzero pairing. It is not
a conserved quantity of the full reaction process.

## The preparations are not on the same elastic orbit

For the same 32-face sign field, the exact initial matrices are

$$I_{\rm overlap}=\begin{pmatrix}0&-4\\4&0\end{pmatrix},\qquad
I_{\rm separated}=\begin{pmatrix}0&0\\0&0\end{pmatrix}.$$

The factor four reflects the particular lifted boundary representatives;
the invariant distinction used here is rank two versus rank zero.
No sequence of elastic updates can identify these two constraint planes.
Thus the old $72$ versus $96$ comparison is not a distance-only experiment.
It remains a valid exact demonstration of shared-relation inhibition.

A nonzero pairing means no independent basis of the two-dimensional
constraint space can have disjoint lifted cycle representatives: disjoint
curves have zero intersection. This is a topological obstruction to
splitting the relations, not a proof of binding. Both cycles can spread
over an arbitrarily large region while still intersecting. Once a
reaction occurs, the obstruction can disappear.

We added a third preparation with nearby but face-disjoint pairs,
$(6,1)$ and $(30,31)$. The minimum dual distance between these pairs is
one. It has $\sigma=0$, like the distant preparation, and initial reaction
rate $96$. Its joint edge support has size seven, the same as the
overlapping preparation: support size alone misses the distinction.

All changing single steps from the three initial planes were enumerated:

| Preparation | Elastic steps and pairing rank | Reaction steps and pairing rank |
| --- | --- | --- |
| Overlapping | 184, $2\to2$ | 72, $2\to0$ |
| Nearby disjoint | 184, $0\to0$ | 96, $0\to0$ |
| Distant disjoint | 184, $0\to0$ | 96, $0\to0$ |

Counts include operator multiplicity. Matching this intersection rank is
not a complete classification of all global topological invariants or
dynamical components. It removes the specific confound identified above.

## A proximity-dependent angular interaction survives that correction

Let $N(t)$ count reactions and $\lambda$ be the instantaneous total
reaction rate. Write the original generator as $G_\alpha=G_0+\alpha A$,
where $A$ contains the two angular directions at unit rate. From a fixed
initial state $x$,

$$\mathbb E_xN(t)=\lambda(x)t+\frac12(G_\alpha\lambda)(x)t^2
+\frac16(G_\alpha^2\lambda)(x)t^3+O(t^4).$$

At an all-reflection state, every angular channel is idle, so $Af(x)=0$
for any observable $f$. The first angular dependence is consequently

$$\boxed{\mathbb E_xN(t;\alpha_1)-\mathbb E_xN(t;\alpha_0)
=\frac{\alpha_1-\alpha_0}{6}(G_0A\lambda)(x)t^3+O(t^4).}$$

Only a first reaction can contribute to $G_0A\lambda(x)$: elastic motion
leaves every face a reflection, so angular motion remains unavailable.
Thus this coefficient is a finite exact sum over one original reaction,
then one original angular update, followed by a rate readout. It is not
obtained by fitting a time series. At fixed mesh and finite rates the
proposal rate and $\lambda$ are bounded, justifying the expansion even
though the integer carrier state space is infinite.

The exact coefficients for the three 32-face preparations are:

| Preparation | $\lambda(x)$ | $(G_0\lambda)(x)$ | $(G_0A\lambda)(x)$ |
| --- | ---: | ---: | ---: |
| Overlapping | 72 | $-4608$ | $-960$ |
| Nearby disjoint | 96 | $-6720$ | $-3456$ |
| Distant disjoint | 96 | $-4992$ | 0 |

The [volume calculation](../data/fiber-constraint-proximity.json)
repeats all three preparations on 128 and 512 faces. Every displayed
coefficient is unchanged. The larger calculation uses the exact read
supports: an elastic write outside the inherited constraint support
leaves the plane and all-reflection sign field unchanged; an angular
channel requires an actual $R/Z$ pair. These exclusions follow from the
rules and do not discard rare eligible channels. The support-local
calculation also reproduces the complete 32-face enumeration above.

For nearby disjoint constraints, 24 first reaction channels each lead to
angular rate drift $-144$, totaling $-3456$. For example, reaction fan 0,
channel 11 enables four angular directions that each reduce reaction
intensity by $36$. No first reaction in the distant preparation allows
an angular change of reaction intensity at the next step.

Writing $\Delta_\alpha$ for the difference between angular rates
$\alpha_1$ and $\alpha_0$, the nearby-minus-distant comparison gives

$$\boxed{\Delta_\alpha\left(
\mathbb E N_{\rm nearby}(t)-\mathbb E N_{\rm distant}(t)\right)
=-576(\alpha_1-\alpha_0)t^3+O(t^4).}$$

The two states have identical face charges, equal initial reaction rate,
and equal intersection rank. Their geometric arrangement nevertheless
changes the leading angular feedback. Independently of angular rate,
their expected counts first separate by $-864t^2+O(t^3)$.

This is a derived local kinetic interaction between joint constraints,
not an attractive or repulsive positional force. It concerns early-time
reaction counts, not a potential, an equilibrium bond, or persistent
localization. The zero cubic coefficient for distant pairs does not
exclude a later interaction after transport.

### Why the leading angular correction cannot be positive at rank two

In a homogeneous zero-offset rank-two preparation, the first reaction
leaves one identity boundary and at most one other independent relation.
No all-reflection fan can then have three equal based reflections: that
would require two independent local relations in addition to the
identity boundary. Equivalently, the two pair cycles of an all-equal
fan have nonzero intersection, contradicting $\sigma=0$ when an identity
boundary already lies in the rank-two space. Its angularly accessible forward-reactive fans have
at most one equality.

Immediately after this first reaction all affine constants are still
zero. An angular unit translation cannot make a nonzero linear form
identically zero; it can only destroy an existing zero equality or leave
it unchanged. There is no all-equal fan whose inhibition it could remove.
Backward reaction rates depend only on the unchanged face charges.
Therefore $A\lambda\le0$ after every first reaction, and

$$\boxed{(G_0A\lambda)(x)\le0}$$

for this class of rank-two starts. This sign restriction is not a
long-time monotonicity theorem. Once offsets have changed, later angular
steps can restore equalities. At rank three, two residual relations can
remain after creating an identity face, so the above exclusion of an
all-equal fan no longer follows; whether enhancement actually occurs
must be calculated rather than assumed.

The [activation follow-up](fiber-angular-activation.md) resolves this:
three constraints still cannot produce a net angular intensity increase
immediately after the first reaction, because openings are offset by
other local closures. An exact flat-edge boundary formula proves the
restriction, and four-constraint constructions attain a positive
conditional response. Averaging over all first reactions from their
precursors nevertheless remains suppressive.

Apparatus and exact states:
[intersection calculation](../tools/fiber_constraint_intersection.py),
[matrices and transition records](../data/fiber-constraint-intersection.json).
