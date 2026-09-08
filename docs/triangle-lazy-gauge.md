# Exact whole-mesh evolution in incrementally maintained tree gauge

The original primitive bank can now evolve complete root-fixed connection
coordinates on the **whole supplied mesh**, without rebuilding those coordinates
after every change. Every independent loop, including global torus information,
is retained. This complements the [overlapping-patch construction](triangle-patch-observer.md):
it does not treat separately canonicalized local patches as autonomous states.

Only connection coordinates change. The microscopic rule tables, supports,
uniform attempted-update scheduler, and conserved charge remain unchanged.
The geometry is still a fixed, supplied triangular torus; this does not establish
emergent space, quantum evolution, or particle binding.

## 1. The complete tree-gauge state

Choose a spanning tree $T$ with root $o$. If $P_v$ is transport along its root-to-$v$
path, a chord $e=(u,v)$ carries the common-root loop

$$H_e=P_v^{-1}U_eP_u.$$

Tree links become identity. The $E-V+1$ chord values retain the complete connection
modulo gauge transformations that fix the root frame. On our torus this gives
$F+1$ independent loops, not merely the $F$ individual face classes. Common-root
conjugation remains internally; a separate invariant observer can quotient it.

Maximal-tree loop descriptions are established gauge-theory tools; see, for
example, [Burbano and Bauer's general-graph treatment](https://arxiv.org/abs/2409.13812).
Here we implement exact incremental maintenance for the existing finite-fiber
primitive bank. We do not import that work's $SU(2)$ group, Hamiltonian, or quantum
degrees of freedom.

## 2. Why restoring a tree edge can look nonlocal

Suppose a primitive changes a parent-to-child tree link from identity to $x$.
Applying $x^{-1}$ to every vertex frame in the child's subtree restores that link
to identity. Internal tree links remain identity, and chord coordinates transform
at their endpoints. That may alter many stored coordinates even though the
physical primitive wrote only one or two local links.

This is a **gauge-coordinate change**, not a physical signal or a nonlocal force.
The new implementation stores its effect lazily instead of visiting every chord.

For a canonical oriented chord $e=(u,v)$, store a base value $B_e$ and vertex-frame
tags $G_v$ so that the current tree-gauge connection is

$$U_e=G_v B_eG_u^{-1},\qquad U_e=1\ \text{for }e\in T.$$

The $G_v$ tags are a computational cache, not extra physical variables. Different
cache histories may encode the same chord state; cache bytes must not be used as
a physical state key.

## 3. One primitive, with all nonabelian ordering retained

1. Decode only its actual read links. A fan involves at most five vertices; a pair
   at most four. Cache those frame queries for the duration of the update.
2. Apply the **existing raw-link oracle and table** through that local view. Buffer
   its writes so the oracle can verify the simultaneous local target.
3. Store each new chord value $U'_e$ in the current frame:

   $$B'_e=G_v^{-1}U'_eG_u.$$

4. Restore changed tree edges **deepest first**. For a parent-to-child value $x$,
   left-multiply all tags in that child's subtree by $x^{-1}$. Reverse-oriented
   stored edges use the corresponding inverse convention.

Deepest-first order matters in a general tree. Restoring a descendant cannot alter
an ancestor's pending edge; a subsequent common conjugation of the larger subtree
keeps the restored identities fixed. Tests use a second valid spanning tree with
nested primitive writes. Deliberately using ancestor-first restoration produces
a detectable nonabelian error. The normal breadth-first tree alone did not exercise
that case, so it was not treated as sufficient coverage.

A depth-first ordering of **vertices in the fixed tree** makes each subtree a
contiguous interval. A segment tree then supports interval left multiplication
and point queries in $O(\log V)$ time. Its lazy tags preserve chronological matrix
order: a newer ancestor tag multiplies older descendant tags on the left.

Consequently one primitive needs at most five point queries, two interval updates,
and a constant number of link operations. It does not scan the mesh. This argument
uses group associativity and inverses, not commutativity.

Full snapshots are separate: a traversal materializes all tags in $O(V)$, then all
chords in $O(E)$. Complete invariant signatures, rendering exports, and global
deduplication keys are not claimed to be logarithmic-time operations.

## 4. What has actually been checked

The checked dataset runs 1,000 original attempted updates at each of four mesh
sizes: 9, 36, 144, and 576 vertices. After **every attempt**, all 19, 73, 289, or
1,153 rooted loops agree exactly with the independent C++ raw connection. Each
complete reverse schedule restores the initial loop tuple.

Additional tests cover:

- All small pairs of noncommuting interval updates, plus interleaved lazy pushes
  and queries on uneven tree sizes.
- Actual primitive evolution in an alternate spanning tree with nested writes.
- Nonconstant local gauge transformations, including a nonidentity root frame.
  The trajectories differ only by the expected fixed common-root conjugation.
- Flat connections with identical local face charges but different global torus
  handles. Their global loop distinction survives; neither is collapsed into the
  other merely because all local charges vanish.
- A guard that forbids whole-mesh loop, snapshot, or frame materialization calls
  from the primitive update path.

The largest audited mesh needed at most five frame queries, two interval updates,
55 point-query tree-node visits, and 40 interval-update node visits in an attempt.
These are observed work counts, not a claim that every larger tree has the same
constant visit count; the tree-depth factor grows logarithmically.

## 5. Performance: useful for gauge maintenance, not faster than raw evolution

The benchmark uses identical schedules and verifies complete final rooted-loop
tuples and event counts for three Python adapters. Initialization and the common
final verification are outside the timing; eager in-loop gauge restoration is
included. Medians below are five runs of 20,000 attempts on an Apple M1 Max.

| Vertices | Raw links | Eager tree-gauge restoration | Lazy tree gauge |
| ---: | ---: | ---: | ---: |
| 9 | 19 ms | 36 ms | 165 ms |
| 36 | 20 ms | 73 ms | 199 ms |
| 144 | 20 ms | 220 ms | 242 ms |
| 576 | 22 ms | 827 ms | 306 ms |

At the largest tested size, lazy maintenance is about **2.7 times faster than eager
restoration**, but about **14 times slower than raw-link evolution**. The raw C++
engine remains the execution reference; this Python adapter is not a replacement
performance winner. Timings are not CI thresholds, and this is not a GPU benchmark.

The benefit is complete gauge-coordinate evolution without a mesh-wide restoration
per event. Storage remains $O(E+V)$ and includes redundant cache data; no memory
compression ratio is claimed. Changing the graph or coordinate tree requires
rebuilding and verifying the coordinate map.

## Reproduce

```bash
uv run --with numpy --with scipy --with python-flint python tools/triangle_lazy_gauge.py --output out/triangle-lazy-gauge.json
diff -u data/triangle-lazy-gauge.json out/triangle-lazy-gauge.json
uv run --with numpy --with scipy --with python-flint python tests/triangle_lazy_gauge.py
uv run --with numpy --with scipy --with python-flint python tools/bench_triangle_lazy_gauge.py
```

The next physics-facing use is to track complete relational observables through
longer encounters and test reduced descriptions against the same full trajectories.
Lazy coordinate changes themselves are not new physical behavior, and must not be
interpreted as propagation or evidence of entanglement.
