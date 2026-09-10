# Four constraints are the minimum for net angular activation after a first reaction

The [rank-two calculation](fiber-constraint-intersection.md) established
angular suppression of early reaction counts. Extending to rank three
does not reverse that sign in the first-born-pair regime. A local counting
identity explains why and gives a constructive threshold: **four
independent constraints suffice, and fewer cannot give a positive net
angular change of total reaction intensity immediately after the first
reaction**, under the preparation assumptions below.

The constructed four-constraint states exhibit both $76\to88$ and
$52\to76$ intensity increases with no face-charge change. Summing over
all angular channels remains positive. These are conditional activation
mechanisms, not evidence of binding or a self-sustaining active phase.
In fact, averaging over every possible first reaction from their
all-reflection precursors still gives a negative leading angular effect.

## Scope of the minimum-rank statement

Start on the supplied triangular torus, side at least four, with all
faces reflections and $r$ generic independent **homogeneous zero-offset**
relations. Use the existing equal-rate reaction bank and exact angular
directions. The fine-fiber interpretation is the fixed-volume,
fixed-time limit described in the earlier affine-carrier notes.

Consider the state immediately after one forward reaction. There is
one identity face $E$, one nonidentity rotation face $Z$, and otherwise
only reflections. Their dual distance is one or two because they were
born in the same three-face fan. All affine constants are still zero.
No intervening elastic spreading or angular offset has yet occurred.

The threshold below applies to this state, not every possible rank-three
state at arbitrary time or in a mixed-charge background.

## Reaction activity is a boundary count on the flat-edge graph

Take the subgraph of the honeycomb dual induced by reflection faces.
Mark an edge flat when the two adjacent reflection holonomies, transported
to their common basepoint, are equal. Equivalently their two-face boundary
holonomy is identity. This mark is gauge invariant and independent of
which endpoint of the common primal edge supplies the basepoint.

For a reflection face $v$, write $d_v$ for its number of reflection
neighbors and $k_v$ for the number of incident flat edges. Each rooted
all-reflection fan is a length-two dual path. Its twelve reaction
channels are enabled precisely when one of its two edges is flat and
the other is not. Therefore

$$\boxed{\lambda_+=12\sum_{v\in R}k_v(d_v-k_v).}$$

This is also a boundary count in the line graph of the reflection
subgraph: marked and unmarked edges form the two sides of the boundary.
No spatial potential has been assigned. The backward intensity is the
separate charge-only count $4N_{ERZ}$.

An angular update on an $R/Z$ edge changes one reflection face $v$.
Because the state has zero affine offsets, a nonzero free linear form
cannot become identically zero by adding a unit constant. All initially
flat $R/R$ edges incident to $v$ lose flatness; other flat-edge marks
stay unchanged. Let $S$ be those $k_v$ neighboring reflection faces.
Removing this flat star from the counting formula gives

$$\boxed{\frac{\Delta\lambda}{12}
=-k_v(d_v-k_v)+\sum_{u\in S}(2k_u-d_u-1).}$$

Backward intensity is unchanged because every face charge is unchanged.
The formula counts all newly opened and closed forward sites, not just
a selected favorable site.

## Why three constraints still cannot give a net increase

Since $v$ neighbors $Z$, it has $d_v\le2$. Each reflection neighbor $u$
of $v$ has $d_u=3$ in this first-born-pair geometry. It cannot also
neighbor $Z$, because the honeycomb has no triangles. If it neighbored
$E$, the path $Z-v-u-E$ together with the birth path of length one or
two would produce a forbidden short cycle or violate bipartite parity.

Three flat edges incident to a reflection face require three independent
relations on the four-face tree. The identity boundary is an additional
independent relation. Independence can be seen from private leaf boundary
edges; the identity face cannot close a short dual cycle around that
tree. Thus for $r\le3$, no such neighbor has $k_u=3$.

- If $k_v=0$, the intensity does not change.
- If $k_v=1$, then $k_u\le2$ and $d_u=3$, so every term in the displayed
  change formula is nonpositive.
- If $k_v=2$, the two incident flat edges already use both relations
  available beyond the identity boundary. Another flat edge at either
  neighbor would require a fourth independent relation. Hence both have
  $k_u=1$, and $\Delta\lambda=-48$.

This proves $\Delta\lambda\le0$ for every angular direction, and hence
$A\lambda\le0$, at ranks at most three in this regime. It also extends
the earlier nonpositive cubic coefficient $(G_0A\lambda)(x)\le0$ to
homogeneous all-reflection rank-three starts.

The completed finite census adds every distinct third adjacent-pair
constraint to the earlier overlapping pair on 32 faces: 46 independent
preparations, 5,376 first-reaction channels, and 3,576 subsequent angular
selections with nonzero intensity change. Every such change is negative.
The census supports the local argument; it is not a classification of
all rank-three affine planes.

## A minimal four-constraint activator

Choose a reflection hub $u$ with three reflection neighbors and make all
three incident edges flat. One leaf $v$ is adjacent to the newborn $Z$.
Impose the identity-face boundary as the fourth relation. Keep all other
raw coordinates formally free. The construction uses an integer affine
parameterization with a determinant-one pivot, not independently assigned
face holonomies.

Here $k_u=d_u=3$ and $k_v=1$. An angular step breaks the $u-v$ relation.
At the hub, two previously all-equal, inhibited fans become reactive.
There are two geometries:

| Birth geometry | $d_v$ | Forward rate | Backward rate | Total intensity change | $A\lambda$ |
| --- | ---: | --- | ---: | --- | ---: |
| $E,Z$ adjacent | 2 | $60\to72$ | 16 | $76\to88$ | 48 |
| $E,Z$ at distance two, with common neighbor $v$ | 1 | $48\to72$ | 4 | $52\to76$ | 96 |

For the adjacent case, opening rate $24$ at the hub costs rate $12$ at
the leaf, leaving a net $12$. For distance two there is no reactive fan
centered on the leaf to lose, so the net gain is $24$. Two roots and two
angular directions produce the four equal positive contributions to
$A\lambda$. All other angular directions have zero intensity change.

These states are reachable by an original reaction. Applying an available
inverse $ERZ\to RRR$ channel gives an all-reflection, homogeneous
rank-four precursor. Applying the same involution returns exactly to
the constructed state. On 32 faces the forward-then-angular witnesses
are $(\text{fan }0,\text{channel }5)$ then $(4,15)$ for adjacent birth,
and $(0,3)$ then $(4,15)$ for distance-two birth.

The same construction on 128 and 512 faces gives identical rates and
angular responses. Actual $D_{1,000,000,007}$ raw connections replay
both updates exactly, including agreement of the entire per-fan rate
array. Those finite raw states are saved with the symbolic constraints.

The reaction changes the inherited intersection rank from four to two.
The angular update preserves intersection rank two while reducing the
zero-offset rank from four to three. Thus activation here comes from
breaking excess local alignment, not adding a new correlation or
changing the charge field.

## A positive conditional response is not bulk amplification

Start the full original process from one of the constructed post-reaction
states $y$, with no anchors frozen. For the accumulated reaction count,

$$\mathbb E_yN(t;\alpha_1)-\mathbb E_yN(t;\alpha_0)
=\frac12(\alpha_1-\alpha_0)(A\lambda)(y)t^2+O(t^3).$$

The two witnesses therefore give respectively $24\Delta\alpha\,t^2$
and $48\Delta\alpha\,t^2$ as positive leading changes. This averages
over every future original channel, not a forced sequence after $y$.

But starting from their all-reflection precursors is a different ensemble.
Enumerating every possible first reaction on 32 faces gives

$$G_0A\lambda=-2928\quad\text{and}\quad-5088,$$

respectively. Their leading angular changes in total reaction counts
are consequently $-488\Delta\alpha\,t^3$ and
$-848\Delta\alpha\,t^3$. Positive-response descendants exist but are
outweighed by suppressive descendants in these initial ensembles.

The result is an exact minimum-rank activation mechanism with a clear
limitation, not an active phase, a bound object, an attractive force, or
molecular emergence. The next question is whether the original dynamics
can replenish such configurations often enough to sustain localized
activity, rather than merely visit them transiently. The two-point
flat-edge count determines an instantaneous rate, not a closed evolution
law: it does not replace the full joint constraints.

Apparatus and observations:
[rank-three census](../tools/fiber_three_constraints.py),
[complete census](../data/fiber-three-constraints.json),
[constructive activation calculation](../tools/fiber_angular_activation.py),
[six volume/geometry witnesses](../data/fiber-angular-activation.json).
