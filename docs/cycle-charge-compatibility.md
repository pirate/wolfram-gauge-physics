# Refining a fiber does not automatically preserve its derived charge

The triangle's rule-selection procedure does not extend unchanged to
larger odd cycle fibers. Individually admitting a positive additive charge
is not enough: combining every such rule can destroy every nonzero common
additive class charge.

For a pentagon fiber, an exact census finds 1,435 minimal closures, of
which 184 meet the triangle's individual positivity and occupied-count
change criteria. Their joint charge space is zero. There are thirteen
maximal jointly positive subsets of those 184 rules; the largest contains
96 and assigns the same charge to both nonidentity rotation classes.

This is an obstruction to a particular bottom-up selection policy, not
to continuous gauge theory in general. No new rule or charge has been
installed in the running triangle model.

## What is derived from the fiber graph

Start with the adjacency of $C_n$, not an input group. The image of one
vertex and one of its neighbors fixes the rest of a cycle automorphism.
Enumerating adjacency-preserving continuations produces $2n$ permutations.
Only afterward identify them in cycle coordinates as

$$t_k(j)=j+k,\qquad r_a(j)=a-j\pmod n.$$

For odd $n$, the reflections form one conjugacy class. The nonidentity
rotation classes are $Z_k=\{t_k,t_{-k}\}$ for
$1\leq k\leq(n-1)/2$. Denote their proposed additive charges by
$q_R,q_k$, with $q(e)=0$.

The existing minimal-rule criterion takes one exchanged pair of ordered
triples and closes it under simultaneous conjugation and reversed
orientation $(a,b,c)\mapsto(c^{-1},b^{-1},a^{-1})$. A closure is an
involution precisely when its exchanged pairs have no conflicting
endpoints. It must preserve ordered product and fix the flat triple.
The triangle selection keeps every such minimal rule that changes the
number of occupied faces and individually admits a strictly positive
additive class charge.

## An explicit obstruction for every odd cycle of size at least five

For each nonzero class $Z_k$, consider the exchanged pair

$$ (r_0,r_0,r_k)\longleftrightarrow(r_0,t_{-k},e). $$

Both products are $r_k$. The conjugation/reversal closure is valid: each
endpoint has trivial simultaneous-conjugation stabilizer, reversal moves
the repeated-reflection pattern or the identity position to a distinct
orbit, and the source and target orbits have different class populations.
There are $4n$ transpositions, moving $8n$ triples.

Each rule individually admits a positive charge and imposes

$$q_k=2q_R.$$

For odd $n\geq5$, there is also the valid closure of

$$ (t_1,t_1,e)\longleftrightarrow(t_2,e,e). $$

The product is $t_2$ on both sides. Both nonidentity rotations have the
same rotation-subgroup centralizer, and their four-element
conjugation/reversal tuple orbits have matching stabilizers and no
conflicting endpoints. Their classes are distinct because
$2\not\equiv\pm1\pmod n$.

This rule individually admits a positive charge but requires

$$q_2=2q_1.$$

The reflection conversions for $k=1,2$ give $q_1=q_2=2q_R$.
Together with the rotation fusion, they force $q_R=q_1=q_2=0$.
Including the reflection conversion for every other $k$ then forces all
remaining $q_k$ to vanish:

$$\boxed{\text{The full individually-positive selected bank has no
nonzero additive class charge for odd }n\geq5.}$$

For $C_3$, $t_2$ and $t_1$ are in the same rotation class. That fusion
would instead require $q_1=0$ and is excluded by individual positivity.
The reflection conversion leaves the triangle's known charge ratio
$(q_R,q_1)=(1,2)$. This explains why the triangle is exceptional under
this selection policy.

The explicit closures and joint nullspaces were computed at
$n=3,5,7,9$, giving dimensions $1,0,0,0$. The argument above, rather than
extrapolation from those four cases, establishes the odd-$n$ statement.

## The conflict exists on actual shared links

For $C_5$, all three witness closures were separately applied through
the existing generic shared-link fan reconstruction on a side-four base
mesh. Only the two internal spokes of the three-face disk change; all
exterior links are held fixed, and applying each involution again returns
the exact initial raw connection.

For the coarse assignment $(q_R,q_1,q_2)=(1,2,2)$, the two reflection
conversions preserve the whole-mesh charge, but the rotation fusion gives

$$Q_{\rm after}-Q_{\rm before}=-2.$$

Its actual based words are $(3,3,0)\to(5,0,0)$ in the saved group
enumeration. This is not an artifact of independently quotienting or
gluing patches. These are isolated mathematical witness replays; the
pentagon bank has not replaced the active triangle dynamics.

## Complete pentagon compatibility classification

In coordinates $(q_R,q_1,q_2)$, the 184 eligible rules give eight distinct
conservation-equation families:

| Charge equation | Minimal rules |
| --- | ---: |
| $q_1=3q_2$ | 1 |
| $q_1=2q_2$ | 40 |
| $2q_1=3q_2$ | 3 |
| $2q_1=q_2$ | 40 |
| $3q_1=2q_2$ | 3 |
| $3q_1=q_2$ | 1 |
| $2q_R=q_1$ | 48 |
| $2q_R=q_2$ | 48 |

Every positive equation plane meets another in a strictly positive ray.
Enumerating those pairwise intersections and collecting every equation
satisfied there gives all maximal compatible subsets. Two independent
equations fix a ray; a lone equation is not maximal because it can be
extended at a positive intersection. In primitive integer normalization:

| Positive charge ray | Compatible eligible rules |
| --- | ---: |
| $(1,2,2)$ | 96 |
| $(1,1,2),(1,2,1),(1,2,4),(1,4,2)$ | 88 each |
| $(1,2,3),(1,3,2),(3,4,6),(3,6,4)$ | 51 each |
| $(1,2,6),(1,6,2),(3,2,6),(3,6,2)$ | 49 each |

These are maximal subsets **within the 184 individually eligible minimal
occupied-changing rules**, not a classification of all composite update
laws. Vacancy transport and class-permuting elastic moves conserve every
one of these additive class charges and do not resolve the choice.
This calculation also does not impose a braid relation or establish
causal invariance for overlapping updates.

No physical principle in the current selection policy chooses one of
these thirteen possibilities. Maximizing the number of admitted rules
would pick $(1,2,2)$, but that is a further policy choice, not a derivation
of an angle-sensitive field energy.

## Why the largest bank is not a smooth curvature cost

Keeping all the reflection conversions enforces $q_R=1$ and $q_k=2$ for
every nonidentity rotation class. This equals word length with respect
to all reflections: a reflection needs one factor and every nontrivial
rotation needs two. The triangle's separate identification with edge
transpositions should not be extended blindly: individual adjacent
vertex swaps cease to be cycle automorphisms for $n\geq4$.

The low fiber mode distinguishes a one-step rotation from a two-step
rotation, but this conserved scalar charge does not. In particular
$t_1$ approaches the identity in the low-mode rotation representation
as $n\to\infty$, while $q(t_1)=2$ remains fixed.

No overall rescaling repairs that discontinuity into a nontrivial
continuous angle cost. If the rescaling tends to zero, the cost of
every nonidentity rotation tends to zero; otherwise an identity jump
remains or the cost diverges. This does not rule out larger-support
observables, different compatible banks, or other scaling constructions.

The follow-up [spectral reaction-arity theorem](cycle-spectral-reactions.md)
examines a quadratic observable derived from the lowest fiber mode,
without adopting it as a physical energy.

Calculation: [cycle_charge_compatibility.py](../tools/cycle_charge_compatibility.py).
Exact closures, charge certificates, raw-link witnesses, and the full
pentagon classification are in
[cycle-charge-compatibility.json](../data/cycle-charge-compatibility.json).
