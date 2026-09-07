# Interacting cell holonomies through explicit fiber transport

The fiber determines $G=\operatorname{Aut}(F)$. The new
`src/gauge_interactions.hpp` implements two reversible word maps on that derived group and
realizes them as connection updates on explicitly supplied cell loops. These are candidate
microscopic dynamics. The choice between them is an input to the experiment; fiber symmetry
alone does not select a law.

## Comparing different fibers' frames

Let $A$ be a cell holonomy based at $p$ and $B$ one based at $q$. Supply an actual connector
path $\gamma:q\to p$ and compute its transport $T$. Compare the two holonomies at $p$ using

$$\widetilde B=TBT^{-1},\qquad P=A\widetilde B.$$

Under independent local frames, $T\mapsto g_pTg_q^{-1}$, so both $A$ and $\widetilde B$
transform by conjugation with $g_p$. Any simultaneous-conjugation-equivariant pair rule
therefore defines a gauge-covariant update. Transport its second target back with
$B'=T^{-1}\widetilde B'T$.

The connector is part of the operator specification. Different paths can give different
results in a curved connection. There is no implicit identification of distant frames and no
assumption of path independence. Current homogeneous fibers are copies of the same finite
graph; interactions across non-isomorphic fibers still require further mathematics.

The implementation changes each loop's exclusive closing link by
$U'=(H'H^{-1})U$. Neither closing link may occur in the other loop or the connector. Simple
loops and existing paths are checked before mutation. This realizes the target holonomies
exactly, preserves connector transport, and leaves all other links unchanged. Other incident
faces containing a written link can nevertheless change; these must be included when a full
cell complex is evolved.

$P$ is a based ordered product. Calling it the holonomy of a geometric outer boundary requires
additional incidence and orientation assumptions, which this generic patch API does not make.

## A sector-changing braid map

The existing Hurwitz primitive is

$$H(A,B)=(ABA^{-1},A).$$

Its outputs have conjugacy classes $([B],[A])$. Any sequence of these moves conserves the
multiset of cell conjugacy classes. It can transport existing curvature but cannot change
that inventory. This is a group-theoretic restriction, independent of simulation size.

We also implement the candidate called `BoundaryShear` in the API:

$$S(A,B)=(ABA,A^{-1}),\qquad S^{-1}(C,D)=(D^{-1},DCD).$$

It is the dual of case (10) in Theorem 1.1 of
[Ito's classification of Wada-type representations](https://arxiv.org/abs/1105.2633).
The word map is established mathematics; our work is its transported link realization,
gauge checks, and finite-fiber comparison with Hurwitz dynamics.

Both maps fix $(1,1)$ and preserve the ordered product $AB$. Both commute with simultaneous
conjugation, so they descend to the gauge quotient. With orientation reversal
$J(A,B)=(B^{-1},A^{-1})$, both satisfy $J F J=F^{-1}$ at the holonomy level.

For a square fiber and a quarter turn $R$,

$$S(R,1)=(R^2,R^{-1}).$$

The half-turn sector differs from the quarter-turn sector and the identity. Thus a single
occupied cell can become two occupied cells while conserving the pair product. This is a
change of curvature sectors, not a measured particle-production process or an energy law.
In particular, $S(A,A^{-1})=(A,A^{-1})$: the entire identity-product pair sector is fixed.

The braid relation $S_{12}S_{23}S_{12}=S_{23}S_{12}S_{23}$ holds for every group: either
side reduces on $(A,B,C)$ to

$$\left(ABCBA,\ A^{-1}B^{-1}A^{-1},\ A\right).$$

Hurwitz satisfies the same braid relation. These identities extend to our link realization
when the loops have exclusive closing links and their connectors compose consistently through
a fixed tree. Tests exercise this with three different cell basepoints and nontrivial connector
transport. Arbitrary overlapping cell patches still require separate compatibility analysis.

## Exact square-fiber census

`data/d4-pair-interactions.json` enumerates all 64 ordered pairs in
$\operatorname{Aut}(C_4)^2$ and quotients by simultaneous conjugation. There are 28 physical
pair states in this **fixed, ordered loop basis**. This does not quotient cell permutations,
base-graph isomorphisms, or changes of loop basis.

| Candidate | Raw pairs changing the sector multiset | Cycles on the 28-state gauge quotient |
|---|---:|---|
| Hurwitz | 0 of 64 | 8 fixed points and 10 two-cycles |
| Boundary shear | 4 of 64 | 8 fixed points, 8 two-cycles, and 1 four-cycle |

Both pass all 512 based-holonomy braid triples. Connection tests additionally check 65,536
local-frame transformations, 1,024 transported three-cell braid identities, exact inversion,
unchanged off-support links, invalid-patch rejection, and commuting independent patches.

The loop observer now canonicalizes conjugacy **inside $\operatorname{Aut}(F)$**. Permutation
cycle type alone would merge the square's half-turn with edge reflections: both have cycle
type $(2,2)$ but belong to different gauge sectors. The normalized permutation character
remains useful as an observable, but is not a complete sector classifier.

Reproduce with:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build -L wgphysics --output-on-failure
./build/wgphysics_pair_census --output out/pair-interactions.json
diff -u data/d4-pair-interactions.json out/pair-interactions.json
```

## Dependencies and the next experiment

Each patch exposes its exact link read/write sets. Disjoint writes, including absence of
write/read overlaps, are a sufficient condition for two patches to commute. Shared read-only
connectors are allowed. These footprints can drive event dependency construction; they are
not yet integrated into the official rewriting engine's event allocator.

A dependency edge records what an update consumed. A spatial-locality test must additionally
bound patch extent and connector length in the graph, then measure propagation with scale and
schedule controls. Branchial relations concern alternative histories; an edge by itself does
not establish entanglement or signaling.

The [shared-face mesh experiment](shared-face-dynamics.md) now implements full incident-face
accounting, link-derived dependencies, flat controls, and exact reversal. It uses a table
selected by the later finite-law census and shows that the independent-cell sector factor
does not close over all mesh faces. Correlation, persistence, and multi-face invariant tests
remain necessary. The braid identities provide specific schedule equivalences to test; they
do not imply general causal invariance of shared-link updates, a continuum limit, or QED.
