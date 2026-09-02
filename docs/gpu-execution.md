# GPU execution strategy

The official `HypergraphRewritingEngine` already supplies a serious CUDA backend: a persistent
device kernel, device-resident work queue and termination detection, structure-of-arrays state,
exact device canonicalization, quotient exploration, and CPU/GPU differential tests. We should
extend that backend rather than create a second graph runtime.

## What modern ML kernels teach us

The transferable lesson from FlashAttention is fusion and IO accounting, not attention itself.
Do not materialize connection factorizations, gauge copies, lifted paths, and holonomy matrices in
global memory. A rewrite worker should load one local base neighborhood and its compact group
tables, generate the physical extension, reduce its cycle transports, form a gauge-quotiented
signature, and emit only a surviving child.

The transferable lesson from sparse GNN runtimes is shape-aware scheduling. Graph degrees, rule
arities, fiber sizes, automorphism-group orders, and cycle lengths are ragged. Work should be
bucketed by `(rule, fiber type, state-size class, cycle-length class)` and mapped differently:

- one thread for tiny sequential transports;
- one subwarp or warp for medium cycle batches;
- one block for large canonicalization or high-degree matching;
- segmented scans/reductions for ragged path and cycle batches.

The transferable lesson from mixture-of-experts/grouped-GEMM kernels is device-side routing of
many heterogeneous small problems through a fixed set of resident workers. Here the “expert” is
a specialized rewrite kernel for a rule/fiber/size bucket. It performs integer permutation-table
operations, not GEMM.

For quantum amplitudes, the useful ML-adjacent architecture is a tensor network rather than a
transformer or learned surrogate. A local base rewrite and its fiber map can be represented as a
small tensor; the evolving state is their contraction graph. Exact contraction, or controlled
SVD/MPS compression when entanglement permits it, can avoid materializing every amplitude branch.
This must be evaluated against exact small closures because bond dimension can still grow
exponentially. It is a representation strategy, not permission to prune physics with a neural
score.

## Device representation

Per fiber template, resident read-only tables contain:

- `mul[g,h]`, `inverse[g]`, and `action[g,fiberVertex]`;
- conjugacy-class or canonical-conjugate IDs;
- valid local rewrite extension tables.

Per state, structure-of-arrays storage contains:

- base hyperedges and offsets;
- one `uint16` transport ID per oriented binary base edge;
- cycle/path offsets for the selected local fundamental cycles;
- amplitudes only for physical sectors, never fresh-frame gauge copies.

Parallel transport is ordered group multiplication. Holonomy over ragged cycles is an associative,
noncommutative segmented reduction. Gauge transformation is `mul[g_v, mul[U_uv, inv[g_u]]]`.
Subdivision is `first=k; second=mul[U, inv[k]]`. These are lookup-heavy kernels; tensor cores are
not the target.

## Integration pipeline

1. Match a base rewrite using the engine's existing signature-partitioned join.
2. Load the affected connection IDs and fiber-template table into registers/cache.
3. Generate connection extensions and quotient fresh-fiber gauge copies immediately.
4. Update only cycles intersecting changed edges; preserve unaffected holonomies.
5. Fuse the gauge signature into the engine's state identity/canonicalization input.
6. Deduplicate and enqueue on device without a host round trip.

The exact CPU implementation remains the oracle. Every CUDA cell must compare states, physical
connection signatures, holonomy sectors, amplitudes, and causal relations against it.

## Expected limits

GPU acceleration changes throughput, not the exponential width of an unpruned multiway system.
The existing engine's own RTX 4090 measurements show that exact graph canonicalization can consume
roughly 79% of device cycles on symmetric quotient workloads, while rewrite and matching can be
negligible. Gauge quotienting before general canonicalization is therefore more valuable than a
standalone “edge propagation” kernel.

The first benchmark gate should compare:

- CPU vs existing CUDA engine without fibers;
- CPU vs CUDA with compact connection tables;
- gauge quotient enabled vs disabled;
- full multiway vs canonical-state quotient exploration;
- tiny, medium, and saturation-width workloads.

No speed claim should be made until canonical state counts, holonomy signatures, and causal
relations agree exactly across all configurations.

On Apple Silicon, `bench/metal_group_transport.py` provides a deliberately narrow MLX custom
Metal-kernel probe for the compact table representation. It checks every device result against a
CPU table oracle before reporting throughput. Run it without modifying the project environment:

```bash
uv run --with mlx --with numpy python bench/metal_group_transport.py
```

This measures only local frame transformation and fixed-length path transport. It cannot be used
as an end-to-end multiway evolution speedup claim.

The first M1 Max run is recorded in `bench/results/m1-max-metal-group-transport.json`. It found an
approximately 0.26--0.28 ms dispatch/evaluation floor from 1,000 through 100,000 items. At one
million items, a 0.318 ms frame-transform kernel corresponds to 3.15 billion transforms/s; a
0.333 ms batch of length-eight paths corresponds to 24.1 billion group multiplications/s. These
figures establish that compact fiber propagation is viable and that small frontiers require
persistent/fused execution to amortize their floor.

## What not to borrow

- Do not use a GNN embedding as physical state identity. The exact engine already measured a
  Weisfeiler--Lehman prebucket followed by exact isomorphism checking and found it slower.
- Do not let a learned router discard rewrite branches. A learned cost model may select a kernel
  shape or rank work, but correctness still requires exact matching and quotienting.
- Do not target tensor cores for permutation-table propagation. The work is integer lookup and
  irregular reduction; tensor hardware becomes relevant only if the quantum state is represented
  as dense or block-sparse tensor contractions.
