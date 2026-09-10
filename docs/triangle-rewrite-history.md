# The auxiliary spin sign is not holonomy of the raw rewrite history

The [spin-cover calculation](triangle-projective-lift.md) found a genuine
projective class of the local transformation group. Examining the actual
rewrite walks now gives a stronger reason not to interpret that class as
emergent physical spin: **two reaction attempts can leave every raw link
unchanged in either order, while their auxiliary spin lifts anticommute.**

This is a statement about that particular proposed interpretation, not
a no-go theorem for quantum physics from rewriting. No dynamics has been
modified to produce these results.

## A witness using just two existing reaction rules

Use the same fixed-boundary three-face disk and exterior reference as in
the [local group calculation](triangle-local-control-group.md). Take
framed state 0, with based face word

$$x=(0,2,3).$$

Let $A$ and $B$ be the existing reaction rules 7 and 8 on fan support 0
(channel indices 16 and 17 in the local calculation). Direct composition
of the raw-link-derived tables gives

$$A^2=B^2=(AB)^2=I,\qquad AB=BA,\qquad Ax=Bx=x.$$

The two chronological histories are therefore

$$x\xrightarrow{A}x\xrightarrow{B}x,
\qquad
x\xrightarrow{B}x\xrightarrow{A}x.$$

Their complete raw configurations agree at every attempt. They have the
same elapsed attempted-update time and the same multiset of attempted
rules. They differ only in the order of externally named inactive rules.
The fixed based word uniquely determines the two internal links; all
remaining links are unchanged by the contained primitives. Thus this is
not merely equality of a coarse charge projection.

On the ten-dimensional boundary-odd observable representation,

$$\operatorname{tr}A=\operatorname{tr}B=\operatorname{tr}(AB)=6.$$

Each involution consequently has a two-dimensional negative eigenspace.
Since $A$ and $B$ commute, their negative projectors commute, and

$$\dim(E_-(A)\cap E_-(B))
=\operatorname{tr}\frac{(I-A)(I-B)}4
=\frac{10-6-6+6}{4}=1.$$

In common orthonormal coordinates they negate planes with one shared
axis. Their spin lifts are, up to signs, Clifford products $e_1e_2$ and
$e_2e_3$. Hence

$$\widetilde A\widetilde B=-\widetilde B\widetilde A,
\qquad[\widetilde A,\widetilde B]=-1.$$

The classical subgroup is a Klein four-group. Its spin preimage is the
quaternion group: the lifts of its three nonidentity elements all square
to $-1$. This is a group-theoretic fact about the observable representation,
not a quaternion variable present on the evolving links.

These two primitives share twelve fixed classical states. Across the 24
channel slots, the exact census finds twelve unordered primitive pairs
with commuting classical actions and anticommuting spin lifts. Some slots
have identical action on this sector, so twelve is a channel-pair count,
not a count of distinct abstract subgroups.

Any functional of the raw configuration history and elapsed attempt clock
must give the same result for the two displayed histories. It therefore
cannot reproduce their relative auxiliary spin sign. Allowing extra
variables to respond differently to ordered inactive proposals would
change the state/history description; that mechanism is not already
derived by taking a spin cover.

## Inverse cancellation is another independent obstruction

An inverse-compatible connection along actual changing state transitions
assigns a transport $U(y,x)$ satisfying

$$U(x,y)=U(y,x)^{-1}.$$

This definition makes immediate retracing have identity holonomy, in any
matrix dimension. For vacancy channel 0, corresponding to support 4 and
rule 0, the actual state path is

$$0\longrightarrow8\longrightarrow0.$$

Its inverse-compatible transport is $I$. But assigning the same global
spin lift of that involution to each attempted update gives
$\widetilde T^2=-I$. The lift of an operation on the entire configuration
set is therefore not, by itself, a connection along the visited state
edges. Making a reverse step use the inverse lift would require an
additional state/orientation-dependent prescription. That would no longer
be the original global-lift assignment.

This argument does not declare two forward attempts to take zero time.
It separates geometric parallel transport from elapsed dynamical time.
In particular it does not rule out a dynamical time phase or other memory
in a different, independently justified state description.

## The longer fixed-state witness also collapses under retracing

The previous projective calculation constructed two words that flip pair
masks $\{1,2\}$ and $\{2,3\}$ while fixing state 0. Their global lifts
anticommute. Following their actual state walks gives:

| Word | Attempts | Inactive attempts | Changing transitions | After inverse cancellation |
| --- | ---: | ---: | ---: | --- |
| First pair flip | 26 | 24 | 2 | Empty walk |
| Second pair flip | 28 | 24 | 4 | Empty walk |
| Their commutator | 108 | 96 | 12 | Empty walk |

Cancellation holds even when changing transitions retain their primitive
channel labels and inverse partners. It is not a consequence of merging
different rules that happen to connect the same endpoints. Thus an
inverse-compatible connection on these changing-event paths gives identity
holonomy even if its fibers have arbitrary matrix dimension.

More generally, inverse-compatible scalar edge phases give identity on a
commutator of two closed walks, because scalar holonomies commute. In this
particular witness the stronger statement holds: both loops already
reduce to the empty walk themselves.

## There are histories, but no selected phase law

The labeled changing-event graph has twenty vertices and 88 unoriented
edges, with inverse slots paired. It is connected, so its fundamental
group before attaching any two-cells is free of rank $88-20+1=69$.
Merging all parallel transitions to form raw-state adjacency leaves 68
edges and free rank 49. These are configuration-walk graphs, not spacetime.

These counts do not select a physical history quotient. We have neither
attached a two-cell to every same-endpoint word nor declared all overlapping
updates reorderable. Idle attempts are omitted only in this changing-walk
diagnostic and remain counted in the attempted-update clock. Retaining
idle loops or imposing additional path relations would give different
history objects, and needs an explicit physical justification.

Assigning phases to graph cycles is possible mathematics, but the raw
rules do not choose those phases. The auxiliary spin class in particular
does not pass the raw-history compatibility checks above.

## What survives this result

The independently derived boundary alignment $\tau$ remains valid: its
known twist circuit actually changes raw links and an exterior-relative
observable, even when the isolated gauge quotient returns. That is genuine
classical relational memory, with an existing interaction readout.

The auxiliary spin-cover deck sign is different. The present result
rules out treating it as an already observed phase of the same raw
history. Further work should distinguish effects measured directly from
evolving holonomies from representation-theoretic structures that require
additional physical state variables or amplitude rules.

Calculation: [triangle_rewrite_history.py](../tools/triangle_rewrite_history.py).
It uses the actual boundary-held primitive tables, exact traces, explicit
walks, and inverse-edge reduction. No complex amplitudes are evolved.
