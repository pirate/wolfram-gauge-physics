# Two joint constraints: contact inhibition and minimal angular feedback

The [one-constraint process](fiber-constraint-dynamics.md) follows one
inherited holonomy relation exactly. To study encounters, we must retain
two relations **jointly**. Evolving two separately labeled carriers would
discard their linear combinations and can give wrong reaction rates.

Two results follow directly from the existing rules. Overlapping flat
reflection pairs inhibit one another's initial reaction channels. Also,
in the all-reflection reference sector, two independent constraints are
the minimum needed for angular motion to affect the limiting reaction
history. Neither result establishes attraction or binding.

## Exact joint affine state

Write the derived cycle automorphisms as $U_e(v)=s_ev+a_e$. Parameterize
the raw shifts together:

$$a=Bz+c,\qquad B\in\mathbb Z^{m\times(m-r)},\qquad
\operatorname{rank}B=m-r.$$

Here $z$ contains formally independent shifts, $c$ counts elementary
cycle steps, and $r$ is the number of inherited independent relations.
All local group words act on these integer affine forms. A zero guard
holds only when every free coefficient and the constant vanish. Thus
the calculation retains every consequence of the joint relations.

The constraint row space and its offset functional are

$$C=\ker B^T,\qquad b(\ell)=\ell\cdot c.$$

A homogeneous raw holonomy covector is identically zero precisely when
it belongs to $C_0=\ker(b|_C)$. Therefore

$$\boxed{\dim C_0=r-\mathbf1_{b\ne0}.}$$

In particular, two independent affine constraints always retain at least
one homogeneous zero relation. That relation may be spatially extended
and invisible to every local reaction guard; this is not guaranteed
local activity.

For each applicable branch, the original inverse rule gives an inverse
on the affine plane. Its dimension is preserved. For these initial
planes the pivot minors have determinant $\pm1$, so free integer
coordinates also parameterize their finite-fiber counterparts without
a divisibility restriction. The bounded-affine-word and Poisson-cutoff
argument in the one-constraint note consequently extends to fixed $r$,
fixed volume, fixed finite angular rate, and fixed observation time.
It does not establish a long-time or finite-density limit.

We measure joint support, the rank of locally visible zero holonomies,
and the rank of relations actually enabling reactions. These satisfy

$$r_{\rm active}\le r_{\rm visible}\le\dim C_0\le r.$$

They do not depend on a choice of constraint basis. To measure coefficient
complexity for $r=2$, form the minors of the extended normals:

$$P_{ij}=\ell^{\rm ext}_{1i}\ell^{\rm ext}_{2j}
-\ell^{\rm ext}_{1j}\ell^{\rm ext}_{2i},\qquad
\ell^{\rm ext}_k=(\ell_k,-b_k).$$

Divide all minors by their common greatest common divisor. Their absolute
values are invariant under rational changes of basis: such a change
multiplies every minor by the same determinant. Vertex reflection frames
only change column signs, while twisted closure cancels vertex translation
frames. Joint support and the primitive absolute-minor norm are therefore
also gauge invariant. Norms at different exterior degrees are not a
common physical size or energy.

## Controlled encounter preparation

Start with every face a reflection, $q_f=1$. The conserved charge is

$$Q=F+N_Z-N_E=F,\qquad N_Z=N_E.$$

This is a homogeneous reference, not a derived vacuum or electromagnetic
neutrality. Sample raw signs uniformly conditional on those face signs,
retaining the global sign sectors. Impose one or two zero holonomies on
adjacent reflection pairs, leaving all remaining shifts generic.

On the 32-face torus, the contact preparation uses pairs $(6,1)$ and
$(1,0)$; the separated preparation uses $(6,1)$ and $(26,27)$. The latter
is selected by dual-graph separation, not by measured reaction activity.
The sign samples and initial charge fields are matched across controls.
These are deliberately prepared encounters, not a stationary ensemble
conditioned on two spontaneous events.

Each isolated flat pair lies in four rooted three-face fans. In the
all-reflection background, each such fan has twelve reaction channels:
the initial rate is $48$. Two sufficiently separated constraints give
$96$. At contact, their common fan has three equal based reflections;
the original exactly-one-equality guard makes all its reaction channels
idle. The other six fans remain active, giving

$$\boxed{\lambda_{\rm contact}=72<96=\lambda_{\rm separated}.}$$

Independent addition would wrongly assign rate $24$ to the common fan.
This is a classical rule-level inhibition, not Pauli exclusion.

The subsequent [intersection calculation](fiber-constraint-intersection.md)
finds an important confound: the contact and separated constraint planes
also have different intersection rank on the auxiliary sign cover. They
are not related by elastic motion alone. Thus this is not a distance-only
comparison. A nearby, face-disjoint preparation has the same zero
intersection rank and initial rate $96$ as the distant preparation.
The same follow-up derives a proximity-dependent angular contribution
to reaction counts at order $t^3$, with identical exact coefficients on
32-, 128-, and 512-face meshes, and proves its nonpositive leading sign
for the homogeneous zero-offset rank-two preparation.

## Why angular feedback first appears at rank two in this sector

Distinct identity faces give independent relations in the present
geometries, so $N_E\le r$. At rank one and $Q=F$ there are only two phases:

- If $N_E=N_Z=0$, all faces are reflections and no angular $R/Z$ move
  is available.
- If $N_E=N_Z=1$, the sole relation is the identity-face boundary.
  Angular moves write an $R/Z$ edge, which is not incident to that
  identity face. Thus $\ell_e=0$ on every angular write, and
  $b'=b+\ell_e\delta=b$.

Angular moves leave the entire rank-one carrier state unchanged. Other
generator terms therefore give the same sign, carrier, charge, and
reaction-history law for every fixed finite $\alpha$. Free raw angles
may move, but no generic limiting guard notices them.

More generally, if $N_E=N_Z=r$, the independent identity-face boundaries
span the whole constraint space, and angular writes again leave every
relation unchanged. Angular feedback requires

$$\boxed{0<N_E=N_Z<r.}$$

Hence $r=2$, during its one-pair phase, is the first possible neutral
sector with angular feedback. This statement concerns the $Q=F$
preparation; it does not contradict angular effects in the earlier
one-constraint experiments with a mixed-charge background.

### An explicit two-step witness

Starting from the contact plane with seed 914407001, apply the existing
reaction at fan 3, channel 9, followed by the existing positive angular
generator direction at fan 0, channel 15. Immediately before the angular
step there is one identity face (31) and one nonidentity rotation face
(6). The angular step leaves **every face charge unchanged**, but changes

$$\lambda_{\rm reaction}:52\longrightarrow16,\qquad
r_{\rm active}=r_{\rm visible}:2\longrightarrow1.$$

Three forward-reactive fans each lose rate twelve. The backward rate
sixteen remains. Thus the instantaneous identity-count drift changes
from $36-16=20$ to $0-16=-16$. This is a nonlinear observable: it does
not imply that the instantaneous drift of an individual linear face
charge changes. The full reaction bank cancels those linear forward
contributions.

The joint support remains seven edges. The lost channels result from
moving one independent relation off zero, not spatial separation or an
imposed interaction potential. Actual $D_{1,000,000,007}$ raw-link
states and the unchanged charge field are retained in the
[complete witness](../data/fiber-two-constraint-witness.json).

## Released trajectories and what they do not establish

The released experiment follows 32 independent sign preparations per
condition to model time one. All joint forms evolve under the original
rules; neither relations nor reference faces are frozen.

At $\alpha=1$, the mean accumulated reaction count is $5.469\pm0.954$
for contact and $8.250\pm1.390$ for separation. The paired difference is
$2.781\pm1.524$ (one standard error), insufficient to establish a
long-time suppression effect. Mean joint support reaches 42 of 48 edges
for contact and 36 for separation. These observations do not demonstrate
a compact bound object.

Across all six conditions, 302,396 proposals were compared with actual
finite-fiber raw-link evolution on the same schedule, with no detected
disagreement. That numerical comparison supports the calculation; the
rank statements, initial contact rates, and angular-visibility result
come from the algebra above.

The [dilute encounter calculation](fiber-dilute-encounter.md) supplies an
exact nonbinding baseline. Long-lived pairs and slowing annihilation can
occur through symmetric diffusion and survival selection alone. A
two-constraint binding claim must establish more than either effect.

Apparatus and observations:
[joint affine evolution](../tools/fiber_two_constraints.py),
[released observations](../data/fiber-two-constraints.json),
[witness construction](../tools/fiber_two_constraint_witness.py).
