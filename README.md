# Wolfram Physics Engine Exploration:<br/>Gauge Theory & the Standard Model

[![CI](https://github.com/pirate/wolfram-gauge-physics/actions/workflows/ci.yml/badge.svg)](https://github.com/pirate/wolfram-gauge-physics/actions/workflows/ci.yml)

Steven Wolfram embarked on the great challenge of discretizing the standard model into a properly computable simulation. Instead of simulating physics top-down by encoding the known equations from the standard model, navier stokes, relativity, etc. into a particle engine, he went the other direction and just a bottom-up graph-based system for modeling propagation of arbitrary information between causal edges (not even necessarily local ones!). Interestingly natural analogs for space and time arise in the evolution of the graph over time, you can slice across the hypergraph of all possible next-step graphs to get different possibilities within the future lightcone, and the number of hypergraph nodes that describe a given quantum state map to probabilistic certainty of outcomes. the same way (Amplitude)^2 of the wavefunction does in normal physics models. These emergent properties of the graph and others are surprising because the graph is a lossy, discretized approximation of a nearly infininitely complex reality, and yet the same rules of space and time and the speed of light fall out of these simple rules applied in a loop to a graph. (just like the game of life)

<table><tr>
<td>
<a href="https://www.youtube.com/watch?v=yAJTctpzp5w"><img src="https://github.com/user-attachments/assets/3695f717-4d87-499e-8d38-3770bdf57508"/><br/><small><code>Stephen Wolfram: Can space and time emerge from simple rules?</code></small></a>
</td>
<td>
<a href="https://www.wolframcloud.com/obj/wolframphysics/Tools/hands-on-introduction-to-the-wolfram-physics-project.nb"><img src="https://github.com/user-attachments/assets/b5bffd92-8464-4061-9711-52046b77b5be"/><br/><small><code>Hands-On Introduction to the Wolfram Physics Project
</code></small></a>
</td>
</tr></table>

This repo aims to progress the state of the art when it comes to modeling the actual laws of our reality in Wolframs hypergraph system.

(goal: see how much of gauge theory and QED, pops out as a result of some minimal rules + simulation time. not expecting to be able to model anything QCD-ish or sub-proton scale just yet, focusing more on getting a simple molecule or two atom system working without hardcoding *any* of the known laws of physics beyond some basic constants and transforms)

This repository asks a deliberately bottom-up question:

> Can gauge structure be computed from discrete fibers and rewrite symmetries before naming a
> continuum group, particle, force, lattice, or molecular geometry?

It combines the fiber/connection hierarchy developed by Wolfram Institute's
[InfraGaugeTheory](https://github.com/WolframInstitute/InfraGaugeTheory) with exact multiway
evolution from the
[HypergraphRewritingEngine](https://github.com/WolframInstitute/HypergraphRewritingEngine), then
adds rewrite-compatible gauge primitives that are currently missing between those two layers.

> [!IMPORTANT]
> This is experimental mathematical software, not a demonstrated derivation of the Standard
> Model, electromagnetism, particles, or continuum spacetime. Every claim below is scoped to an
> implemented finite combinatorial construction and its tests.

## What is implemented

| Capability | Construction | Evidence |
|---|---|---|
| Gauge group from a fiber | Exhaustively compute `Aut(F)` from the fiber graph itself | A square fiber derives the eight-element nonabelian `D4` group |
| Discrete connection | Fiber isomorphisms attached to oriented base edges | Reverse orientation gives the inverse map |
| Parallel transport | Ordered composition along base paths | Exact permutation result |
| Curvature / holonomy | Transport around closed paths | Gauge-invariant conjugacy signature |
| Local gauge transformations | Independent changes of fiber frame at every base vertex | Holonomy changes only by conjugation |
| Horizontal sections | Fiber assignments preserved by every edge transport | Fixed points agree with holonomy obstruction |
| Gauge-preserving rewrites | Factor old transport through each newly created fiber | Boundary transport is unchanged |
| Fresh-fiber gauge quotient | Prove all rewrite factorizations form one local gauge orbit | `|Aut(F)|` labelings become one physical extension |
| Quantum rewrite isometry | Normalized amplitude over that orbit | Norm one; distinct boundary holonomies are orthogonal |
| Exact connection quotient | Gauge-fix a spanning forest and canonicalize chord holonomies | Invariant under arbitrary local frame changes |
| GPU representation | Dense `uint16` multiplication, inverse, and action tables | Exact CPU/Metal differential check |
| Multiway base evolution | Pinned official C++ rewriting engine | Exact isomorphism-aware quotient exploration |
| Fiber inference | Exact canonical classes of induced vertex links | A wheel derives a square link and its order-eight symmetry |
| Gauge-aware rewrite identity | Base relabeling plus local-frame quotient | Eight subdivision labelings register as one physical child |
| Generalized fibers | Non-isomorphic fibers and partial induced embeddings | Missing horizontal lifts remain explicitly undefined |
| Dynamics census | Reversible maps preserving identity, inversion, and conjugation | Eight exact candidates for `Aut(C4)` |
| Unary-dynamics no-go | Exhaustively quotient every `D4` candidate by conjugacy | All eight maps fix all five physical curvature sectors |
| Two-cell curvature transport | Reversible Hurwitz exchange on adjacent oriented cells | Moves a nontrivial sector while conserving total holonomy |
| Explicit two-complex product | Oriented faces carried through official engine events and exact identity | Shared-edge subdivision updates every incident face |
| Causal curvature | Loop sectors compared with causal event reachability | Reports propagated and off-causal changes separately |
| Compression diagnostic | Exact Schmidt spectrum and discarded-norm bound | The `D4` subdivision orbit is rank eight and flat-spectrum |
| Engine × gauge product | Connection sectors propagated through actual raw engine events | 16 raw states and 15 events become three physical depth sectors |
| Automatic causal sectors | Before/after holonomy attached to every supported engine event | Six real causal edges and zero subdivision curvature violations |
| Rule × fiber census | Every order-four fiber crossed with the subdivision closure | 18 fibers, every conjugacy sector, exact orbit-norm checks |

The implementation is generic over finite fiber graphs. It does **not** declare `U(1)`, `SU(2)`,
or another Lie group as its microscopic starting point. Different fibers derive different finite
automorphism groups; a continuous effective group, if any, remains a coarse-graining question.

## Thirty-second demo

Requirements: CMake 3.20+ and a C++20 compiler.

```bash
git clone https://github.com/pirate/wolfram-gauge-physics.git
cd wolfram-gauge-physics
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build -L wgphysics --output-on-failure
./build/fiber_demo
./build/research_demo
./build/cell_dynamics_demo
```

Expected core output:

```text
fiber vertices:              4
derived group order:         8
derived group is nonabelian: true
fundamental cycle count:     1
rewrite factorizations:      8
physical fresh-fiber orbits: 1
quantum orbit support:       8
quantum orbit norm:          1
```

Run a Wolfram-model rule with exact state canonicalization:

```bash
./build/wgphysics_evolve \
  --rule '0,1;0,2->0,2;0,3;1,3;2,3' \
  --init '1,2;1,3' \
  --steps 3 \
  --output out/evolution.json
```

Generate the deterministic small-fiber census:

```bash
./build/wgphysics_census --output out/fiber-census.json
./build/wgphysics_product_census --output out/rule-fiber-census.json
./build/wgphysics_dynamics_census --output out/d4-cell-dynamics-census.json
```

Run the exact engine × gauge product evolution:

```bash
./build/wgphysics_product_evolve --steps 2 --output out/product-evolution.json
./build/wgphysics_product_evolve \
  --steps 2 --cell-dynamics-index 7 --output out/product-unary-rule-7.json
```

The second export deliberately still reports zero physical curvature changes: candidate 7 moves
six of eight based `D4` elements but, like every admissible unary candidate, fixes all conjugacy
sectors. This is a reproducible no-go check, not a propagation demo.

## Meaningful extensions relative to the inspected upstream state

The comparison is against `InfraGaugeTheory` commit
`edd9bdca46b7838d6b3e940e8ae8cde90b60ef2c` (paclet 1.0.7) and
`HypergraphRewritingEngine` commit `03fe60ddf338983060b6bb4b23e8b4b5d7ae7337`.

1. **The local group is derived, not supplied.** `Aut(F)` is enumerated from microscopic fiber
   adjacency. This gives a concrete bridge from local combinatorial symmetry to gauge freedom.
2. **Gauge transformations act on complete connections.** The code implements
   `U_xy -> g_y U_xy g_x^-1` and tests the conjugacy invariance of closed holonomy.
3. **Gauge redundancy is quotiented exactly before general graph isomorphism.** A deterministic
   spanning forest removes tree transports; fundamental-cycle holonomies are canonicalized under
   simultaneous conjugation at each component root.
4. **Base-space rewrites carry their fibers.** For `(x--y) -> (x--w--y)`, every factorization
   `U_xy = U_wy U_xw` is generated and proven to be one gauge orbit at the fresh vertex.
5. **The rewrite has a quantum map.** The fresh-frame orbit is a normalized superposition rather
   than a collection of falsely distinct multiway branches.
6. **The representation is ready for persistent GPU evolution.** Repeated transport is compiled
   to compact integer tables and verified against a CPU oracle.
7. **Physics-facing failures are explicit.** Tests distinguish flat and curved sectors, obstructed
   and unobstructed horizontal sections, gauge copies and physical sectors, and abelian versus
   nonabelian derived groups.

These are software and finite-model improvements. We do not claim priority over unpublished work
or claim that the construction is the unique microscopic interpretation of Wolfram gauge
symmetry. See [State of the art and contribution boundary](docs/STATE_OF_THE_ART.md).

## GPU result

The compact group-table kernels were tested on an Apple M1 Max 32-core GPU using a custom MLX
Metal kernel. At one million independent items:

- local frame transformations: **0.318 ms**, approximately **3.15 billion/s**;
- one million length-eight paths: **0.333 ms**, approximately **24.1 billion group
  multiplications/s**.

Every device result matched the CPU table oracle exactly. This is a kernel microbenchmark, not an
end-to-end simulation claim; multiway width and exact state canonicalization remain the expensive
parts. Reproduce it with:

```bash
uv run --with mlx --with numpy python bench/metal_group_transport.py
```

Raw measurements are in
[`bench/results/m1-max-metal-group-transport.json`](bench/results/m1-max-metal-group-transport.json).

## Research program: useful next contributions

### 1. Infer fibers from rewrites — first exact prescription implemented

The code now extracts induced open-link graphs, groups them by exact graph isomorphism, and derives
their automorphism groups. A wheel exposes a square link with order-eight symmetry. The stronger
Wolfram-model result is to compare this prescription with local rule automorphisms, branchlike
paths, and repeated neighborhoods in the multiway causal graph.

### 2. Put connection state inside multiway evolution — real event product implemented

Small uniform-fiber states now have a joint identity under arbitrary base relabeling and local frame
changes. The product runner carries connections through every raw state and actual consumed/produced
edge event from the official engine. At depth two, 16 raw states and 15 events quotient to three
physical product states. Unsupported morphologies fail closed. The next step is embedding the joint
identity in the official engine's arena and match keys so exact quotient exploration can happen
without first expanding raw provenance.

### 3. Generalize the bundle — partial non-isomorphic lifts implemented

The code supports different graph fibers at different base vertices and injective partial induced
embeddings as horizontal lifts; undefined lifts fail explicitly. Twisted bundles, directed/ribbon
fibers, hypergraph fibers, higher connections, and dynamically changing fiber type remain open.

### 4. Search for dynamics instead of inserting forces — finite census implemented

The first census exhaustively enumerates reversible maps on the derived finite group and keeps only
updates preserving identity, edge-orientation inversion, and conjugation equivariance. `Aut(C4)`
has eight such maps, but the exact product census now rejects all eight as physical unary dynamics:
each fixes all five conjugacy sectors despite sometimes changing a based group element. The first
nontrivial replacement is a two-cell Hurwitz exchange `(A,B) -> (ABA^-1,A)`, which is reversible,
gauge-covariant, conserves `AB`, and transports curvature between neighboring cells. Carrying
oriented face boundaries through engine rewrites is now exact; introducing independent gauge events
with causal dependencies derived from their face/link read and write sets is the next step.

### 5. Connect gauge propagation to causal propagation — automatic event attachment implemented

The product runner now attaches exact before/after holonomy sectors to every supported real engine
event and imports the engine's causal edges automatically. An analyzer measures changes reachable
from declared disturbances, maximum causal depth, off-causal changes, and causal alignment.
Transport-preserving subdivision is a zero-change control. A derived nontrivial update law,
propagation-speed convergence, and foliation comparisons remain next.

### 6. Compress amplitudes without silently pruning physics — first obstruction measured

The elementary `D4` subdivision orbit has an exact flat rank-eight Schmidt spectrum: retaining rank
four necessarily discards half its norm. This rules out a free low-rank shortcut for that primitive.
The next step is to evaluate multi-event contraction graphs and permit controlled tensor-network
compression only where spectra decay, always reporting discarded norm and observable error.

### 7. Integrate the GPU data path — physical-only CPU path implemented

Gauge-aware evolution now fixes the fresh frame and constructs one physical subdivision
representative directly, while retaining the raw orbit size as event metadata. It no longer
materializes `|Aut(F)|` copies before deduplication. Dense multiplication, inverse, and action table
buffers are exposed as the device ABI. The next step is to fuse match, affected-cycle update,
canonical signature, dedup, and queue insertion in the rewriting engine's persistent CUDA kernel.

### 8. Publish a rule-and-fiber census — first product closure published

The committed dataset covers all 18 unlabeled simple graph fibers with one through four vertices.
It records automorphism-group order and commutativity, triangle curvature-sector count, admissible
local-dynamics count when the exhaustive search is bounded, subdivision orbit reduction, exact
Schmidt rank, and half-rank discarded norm. A second dataset now crosses all 18 fibers and every
initial conjugacy sector with the real subdivision-engine closure, recording raw and physical state
counts, actual events and causal edges, curvature violations, and orbit-norm verification. The next
census should add nontrivial rule morphisms, effective dimension, persistence, and localization.

## Exact product-evolution result

The official engine is run with full raw provenance because quotienting the bare base first can
erase matches distinguished by an attached connection. The product layer reconstructs each
supported subdivision from immutable consumed/produced edge IDs, propagates its connection, and then
applies the joint base-isomorphism/local-gauge quotient.

For a `C4` fiber over an initial triangle at depth four:

- 436 raw states become **five physical product states**;
- all 435 actual events receive exact connection sectors;
- all 336 engine causal edges are retained;
- curvature changes remain zero, as required for transport-preserving subdivision;
- median end-to-end CPU time was **0.852 s** over three M1 Max repetitions.

See [Exact engine × gauge product evolution](docs/product-evolution.md) and the
[raw scaling result](bench/results/m1-max-product-scaling.json).

## Questions for the Wolfram community

1. Which object should be regarded as primary: a fiber extracted from rule automorphisms, a fiber
   over spatial vertices, or branchlike equivalence classes in the multiway causal graph?
2. Is `Aut(F)` the right microscopic local symmetry object, or should the relevant group be a
   stabilizer/quotient determined jointly by the fiber and allowed rewrites?
3. Which InfraGaugeTheory data model should a rewrite extension preserve for eventual upstream
   compatibility: total graph plus projection, connection subgraph, or explicit transport maps?
4. Which observables would constitute the strongest early evidence for an emergent gauge field:
   holonomy propagation, center sectors, confinement-like loop scaling, causal invariance, or
   something else?
5. Which small rule/fiber families are the best shared benchmark suite?

## Repository map

- [`src/infragauge.hpp`](src/infragauge.hpp) — fiber graphs, derived automorphisms, connections,
  holonomy, exact gauge quotient, rewrite extensions, and quantum orbit amplitudes.
- [`src/research.hpp`](src/research.hpp) — inferred link fibers, gauge-aware evolution identity,
  oriented cell complexes, partial lifts, dynamics census, causal-curvature analysis, and finite
  research diagnostics.
- [`src/main.cpp`](src/main.cpp) — exact multiway rewrite export harness.
- [`src/product_evolution.hpp`](src/product_evolution.hpp) — strict engine-event morphisms, joint
  product identity, automatic curvature sectors, causal attachment, and rule-by-fiber census.
- [`tests/infragauge.cpp`](tests/infragauge.cpp) — finite differential/invariance tests.
- [`docs/infragauge-foundations.md`](docs/infragauge-foundations.md) — mathematical construction.
- [`docs/gpu-execution.md`](docs/gpu-execution.md) — persistent-kernel integration plan.
- [`docs/product-evolution.md`](docs/product-evolution.md) — product-state semantics, supported
  morphism, exactness boundary, causal integration, and measured scaling.
- [`docs/causal-dynamics.md`](docs/causal-dynamics.md) — exact unary no-go census and reversible
  two-cell curvature transport.
- [`docs/STATE_OF_THE_ART.md`](docs/STATE_OF_THE_ART.md) — related work and contribution boundary.
- [`docs/research-milestones.md`](docs/research-milestones.md) — exact definitions, evidence, and
  remaining limitations for milestones 1–8.
- [`data/fiber-census-n4.json`](data/fiber-census-n4.json) — deterministic census of every unlabeled
  simple fiber graph through four vertices.
- [`data/rule-fiber-census-subdivision.json`](data/rule-fiber-census-subdivision.json) — all 18
  fibers and curvature sectors crossed with the real depth-two subdivision closure.
- [`data/d4-cell-dynamics-census.json`](data/d4-cell-dynamics-census.json) — every admissible unary
  `D4` map and its induced physical-sector and causal behavior.
- [`FORUM_POST.md`](FORUM_POST.md) — concise draft for a Wolfram Community introduction.

## Build without the rewrite-engine dependency

The fiber library, tests, and demo are header-only apart from the executable entrypoints:

```bash
cmake -S . -B build -DWGPHYSICS_BUILD_REWRITE_HARNESS=OFF
cmake --build build -j
ctest --test-dir build -L wgphysics --output-on-failure
```

## License and relationship to Wolfram Institute

MIT licensed. This is an independent experimental project and is not an official Wolfram
Institute or Wolfram Research repository. It depends on and cites their MIT-licensed research
software; no Wolfram Language source is vendored.
