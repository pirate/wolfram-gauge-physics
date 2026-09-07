# A boundary-fixed shared-edge primitive restores transported charges

The [face-energy obstruction](face-energy-obstruction.md) applies to the **exclusive closing-link
lift**: two selected faces choose a move that also changes unread spectator faces. It does not
apply to every realization of the same pair table. A two-face disk admits a different, uniquely
determined one-edge realization that changes exactly the faces used to choose the update.

This is a new microscopic dynamics option, not an optimization of the old dynamics. Both lifts
remain implemented and separately labeled. No force law, energy function, or desired trajectory
was added to the pair table.

![Actual initial and final gauge-class states from the shared-edge simulation](images/shared-edge-trails.png)

The figure uses recorded link states from trial zero, not a fitted field or particle template.
The periodic triangulation and its coordinates are supplied geometry, not emergent spacetime.

## Derivation from the two-face boundary

Take consistently oriented faces sharing the edge $u\to v$. Base both loops at $u$, choosing
the first loop to start with that edge and the second to end with its reverse. Let $S=U_{uv}$.
If $Q$ and $P$ are the exterior paths from $v$ to $u$ and $u$ to $v$, their holonomies are

$$A=QS,\qquad B=S^{-1}P,\qquad AB=QP.$$

The common edge cancels in $AB$, which is exactly the exterior boundary transport. For any
product-preserving rule $F(A,B)=(X,Y)$, keep every exterior link fixed and set

$$\boxed{S'=Q^{-1}X=SA^{-1}X.}$$

Then $A'=QS'=X$ and $B'=(S')^{-1}P=X^{-1}AB=Y$. This is the unique shared-link value realizing
those targets with the boundary fixed. A nontrivial target requires at least one link change,
so this realization has minimum write support on this disk.

Under independent vertex frames, $S\mapsto g_vSg_u^{-1}$, $A\mapsto g_uAg_u^{-1}$, and
$X\mapsto g_uXg_u^{-1}$. Therefore $S'\mapsto g_vS'g_u^{-1}$. The construction is locally
gauge covariant. An inverse pair map recovers the original $S$ exactly because the exterior
path $Q$ is unchanged.

For a manifold interior edge, these are its only two incident faces. Thus no spectator face
changes. The implementation rejects a third incident face, incompatible orientations, and
face pairs that do not form a two-face disk. These restrictions matter: the result is not a
construction for arbitrary nonmanifold hypergraphs.

Both rooted orientations of every interior edge are enumerated. They are distinct candidate
operators; equivalence under exchanging their roots or reordering overlapping events is not
assumed. The scheduler uses actual link read/write conflicts and records actual read dependencies.

## Global conservation without imposing an energy

If a class function $q$ satisfies the pair equation

$$q(A)+q(B)=q(X)+q(Y),$$

then $\sum_f q(H_f)$ is conserved on **any consistently oriented mesh supporting these patches**:
only the two selected face terms change. Rotating a face basepoint only conjugates its holonomy.
Unlike the previous lift, this permits a density to move between faces while its global sum stays
fixed.

We reused table 11229 from the existing minimum-support strict braid screen. Solving all 4,096
four-face patch equations afresh gives three normalized charges, the indicators of the three
noncentral conjugacy classes of the derived square-fiber group. In the repository's element
ordering these classes have representatives $1,2,3$; the identity $0$ and central half-turn $5$
both have zero charge. The unchanged two spectator faces are included in every equation.

The old lift's charge rank is zero; the new lift's rank is three. The earlier obstruction remains
valid under its stated hypotheses. Changing the written edge changes those hypotheses and the
actual dynamics, rather than violating the theorem.

These are conserved combinatorial densities. They have not been identified with electric charge,
mass, or energy. In particular, the zero-charge central half-turn still has nontrivial curvature.

## An exact class-level explanation of the curvature trail

The full five-class update is autonomous for this table. Its zero-charge states can be written
as $v_0=[1]$ and $v_1=[z]$, where $z$ is the central half-turn. The three nonzero species are
$c_1,c_2,c_3$. Exhaustive projection of the microscopic table yields

$$
(v_b,c_i)\longmapsto(c_i,v_{b\oplus\eta_i}),\qquad
(c_i,v_b)\longmapsto(v_{b\oplus\eta_i},c_i),\qquad
(\eta_1,\eta_2,\eta_3)=(1,0,0).
$$

Pairs with two occupied faces or two zero-charge faces are unchanged at the class level.
After forgetting the vacancy bit, this is exactly the colored exclusion factor already derived
from the pair table, now realized on a genuine two-dimensional shared-link mesh. The bit records
an additional reversible memory carried by the zero-charge faces. Reflection-class-1 motion
toggles it; the other two species do not. The asymmetry is a property of this selected table,
not a claim that these combinatorial species are physical particles.

Consequently a reflection seed can leave many central-curvature faces while retaining only two
noncentral defects. A rotation seed in an initially flat background cannot create central faces
under this rule and keeps exactly two nonflat faces. The observed contrast follows from the
microscopic table, not a fitted phenomenological model.

### A topology-qualified parity law

If the dual face-adjacency graph is bipartite, let $\epsilon_f\in\{0,1\}$ be a fixed bipartition.
The same class rule conserves

$$I=\left(N_{v_1}+\sum_{f:H_f\in c_1}\epsilon_f\right)\bmod2.$$

Every $c_1$ hop flips both contributions; every other hop changes neither modulo two. The
analyzer checks all 50 class transitions with both bipartitions of an edge, reconstructs the
dual bipartition from adjacency, and verifies $I$ after every replayed layer. It rejects an odd
cycle rather than assigning a bipartition to it. Exchanging the two bipartition labels changes
$I$ by the already conserved number of $c_1$ defects modulo two.

This is a torsion-valued, sublattice-dependent invariant, not a real-valued energy or quantum
phase. Its independence from other topological constraints has not been established.

There is also a kinematic closed-surface check: projected to the abelianization
$G/[G,G]$, all oriented face holonomies multiply to identity because each edge occurs with
opposite orientations. For the square fiber, $[G,G]=\{1,z\}$. Every recorded layer satisfies
this constraint. It is not an independently derived electromagnetic Gauss law.

## Verification and measured evolution

- The C++ test enumerates all $8^5=32,768$ raw link assignments on a two-triangle disk, with
  both rooted orientations: **65,536 updates**. It checks the compiled $SA^{-1}X$ kernel against
  an independent exterior-path $Q^{-1}X$ realization, each face target, every exterior link,
  and exact raw-link inversion. The test uses the non-involutive BoundaryShear rule, so inverse
  correctness is not merely an involution check.
- A torus test verifies local-frame covariance, read dependencies, and commuting layers. Guards
  reject nonmanifold incidence and incompatible face orientations.
- The selected strict table is run on 288 faces for 256 layers with eight schedules, and on
  1,152 faces for 512 layers with four independently seeded schedules. Each schedule uses five
  initial conditions: flat, reflection seed, rotation seed, a noncommuting pair, and random links.
  All **60 runs** preserve the derived charges and reverse to every initial raw link, including
  an exact face-histogram echo at every intermediate layer.
- Independent Python path algebra replays **4,915,855 forward updates**, checks every layer's
  face histogram and the final raw links, and measures actual charge hops and vacancy flips.
  The propagation bound is at most one dual-face hop per conflict-free layer.

The 288-face reflection runs end with 5–93 nonflat faces, despite retaining exactly two charged
defects; the 1,152-face replication ends with 5–59. Every rotation run retains exactly two nonflat
faces while those defects travel. These ranges are schedule-sensitive observations, not estimates
of a limiting density, diffusion coefficient, or binding energy. One reflection schedule on the
smaller mesh stays within dual distance two; others reach distances fourteen or fifteen. Neither
mixing nor schedule-independent propagation follows from these runs.

## What is still missing

The current pair table still projects exactly to an autonomous exclusion process. Hidden
nonabelian connection information cannot feed back into those projected charges. Recovering
transported charges solves one problem, but does not by itself produce interacting matter,
complex amplitudes, a physical energy scale, or a bound state.

The next primitive search should include complete three-or-more-face neighborhoods, preserve
their exterior transports, and derive invariants and observer closure together. The discriminating
target is conserved transport **with feedback from relative loop data**, not another autonomous
colored-particle factor. Ordinary gauge covariance and pair-level braid identities are not enough
to select physical dynamics or guarantee consistency of overlapping mesh updates.

## Reproduce

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build -L wgphysics --output-on-failure
python3 tests/shared_edge.py
python3 tools/run_shared_edge.py --full-replay --output out/shared-edge.json
diff -u data/d4-shared-edge.json out/shared-edge.json
python3 tools/run_shared_edge.py --side 24 --layers 512 --trials 4 --seed 1819031 --full-replay --output out/shared-edge-replication.json
diff -u data/d4-shared-edge-replication.json out/shared-edge-replication.json
uv run --with matplotlib python tools/plot_shared_edge.py
```

The JSON includes initial/final raw links, complete schedules, layer histories, derived charge
bases, the exact class-memory rule, and independent replay counts. No molecular template enters
the construction or these measurements.
