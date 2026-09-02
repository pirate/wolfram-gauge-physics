# Wolfram Gauge Physics

[![CI](https://github.com/pirate/wolfram-gauge-physics/actions/workflows/ci.yml/badge.svg)](https://github.com/pirate/wolfram-gauge-physics/actions/workflows/ci.yml)

Experimental combinatorial gauge and fiber dynamics for Wolfram-model hypergraph rewriting.

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
| Causal curvature | Loop sectors compared with causal event reachability | Reports propagated and off-causal changes separately |
| Compression diagnostic | Exact Schmidt spectrum and discarded-norm bound | The `D4` subdivision orbit is rank eight and flat-spectrum |

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
```

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

### 2. Put connection state inside multiway evolution — exact subdivision prototype implemented

Small uniform-fiber states now have a joint identity under arbitrary base relabeling and local frame
changes. Gauge-aware subdivision evolution registers eight raw `D4` factorizations as one physical
child. The next step is embedding this identity in the official engine's arena, match keys, event
identity, causal relations, and CPU/GPU differential tests.

### 3. Generalize the bundle — partial non-isomorphic lifts implemented

The code supports different graph fibers at different base vertices and injective partial induced
embeddings as horizontal lifts; undefined lifts fail explicitly. Twisted bundles, directed/ribbon
fibers, hypergraph fibers, higher connections, and dynamically changing fiber type remain open.

### 4. Search for dynamics instead of inserting forces — finite census implemented

The first census exhaustively enumerates reversible maps on the derived finite group and keeps only
updates preserving identity, edge-orientation inversion, and conjugation equivariance. `Aut(C4)`
has eight such maps. They remain candidates: reject them unless coupled evolution produces
propagating curvature disturbances, stable localized sectors, reproducible scattering, and bound
composites without microscopic particle names or potentials.

### 5. Connect gauge propagation to causal propagation — analyzer implemented

An exact analyzer now measures changed holonomy sectors reachable from declared disturbances along
causal event edges, maximum causal depth, off-causal changes, and causal alignment. Automatic
event-sector attachment, propagation-speed convergence under refinement, and comparison of
gauge-equivalent observables across foliations remain next.

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

### 8. Publish a rule-and-fiber census — first exact dataset published

The committed dataset covers all 18 unlabeled simple graph fibers with one through four vertices.
It records automorphism-group order and commutativity, triangle curvature-sector count, admissible
local-dynamics count when the exhaustive search is bounded, subdivision orbit reduction, exact
Schmidt rank, and half-rank discarded norm. The next census should cross these fibers with small
hypergraph rules and add causal propagation, effective dimension, persistence, and computational
cost—including negative results.

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
  partial lifts, dynamics census, causal-curvature analysis, and finite research diagnostics.
- [`src/main.cpp`](src/main.cpp) — exact multiway rewrite export harness.
- [`tests/infragauge.cpp`](tests/infragauge.cpp) — finite differential/invariance tests.
- [`docs/infragauge-foundations.md`](docs/infragauge-foundations.md) — mathematical construction.
- [`docs/gpu-execution.md`](docs/gpu-execution.md) — persistent-kernel integration plan.
- [`docs/STATE_OF_THE_ART.md`](docs/STATE_OF_THE_ART.md) — related work and contribution boundary.
- [`docs/research-milestones.md`](docs/research-milestones.md) — exact definitions, evidence, and
  remaining limitations for milestones 1–8.
- [`data/fiber-census-n4.json`](data/fiber-census-n4.json) — deterministic census of every unlabeled
  simple fiber graph through four vertices.
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
