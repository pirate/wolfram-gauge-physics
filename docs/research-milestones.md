# Research milestones 1–8

This document separates implemented finite constructions from the larger physical claims they are
intended to test. The executable evidence is `tests/research.cpp`; run it through CTest or use
`research_demo` for a short readable summary.

## 1. Candidate fibers inferred from local structure

`infer_link_fibers` constructs the graph induced by each base vertex's open neighborhood, computes
an exact canonical graph code, groups repeated isomorphism classes, and derives each representative's
automorphism order. A wheel graph therefore exposes a square link with `Aut(C4) = D4` without
supplying `D4` or a continuum group.

This is one falsifiable fiber-extraction prescription, not evidence that open links are uniquely the
correct Wolfram-model fibers. The next comparison should evaluate rule automorphisms, causal-graph
branch classes, and repeated-radius neighborhoods on the same rule census.

## 2. Gauge state participates in rewrite identity

`canonical_connection_state` jointly canonicalizes the base graph under every small-graph vertex
relabeling and the connection under every local fiber-frame change. `GaugeSubdivisionEvolution`
uses that identity while registering children. For a square fiber, eight raw factorizations of an
edge subdivision become exactly one physical child, while a different holonomy sector remains a
different state.

The product runner now propagates this state through every supported raw event emitted by the
official engine and performs the joint quotient afterward. It is not yet embedded in the engine's
internal arena, match keys, or CUDA storage. The unquotiented run remains the exact oracle until the
engine's own canonical identity includes connection data.

## 3. Non-isomorphic fibers and partial lifts

`PartialFiberMap` is an injective partial induced-graph embedding. `VariableFiberBundle` attaches a
different finite graph to each base vertex and composes partial horizontal lifts along a path. A
missing image returns an explicit failed lift; the code never completes it by choosing a convenient
vertex.

This handles a useful strict subset of the general InfraGaugeTheory data model. Twisted bundles,
hypergraph fibers, non-injective relations, and dynamically changing fiber type remain open.

## 4. Exact local-dynamics census

`search_local_gauge_dynamics` exhaustively enumerates bijections on the derived finite group and
retains only maps that:

- fix the identity;
- commute with inversion, so reversing an oriented edge remains consistent;
- commute with conjugation, so the update is independent of the local fiber frame.

For `Aut(C4)`, eight maps pass. This is a finite candidate-rule census, not a force law. Candidates
still have to be coupled locally across a base graph and rejected unless they show robust propagation,
localized sectors, scattering, or binding.

## 5. Curvature changes measured against causality

`curvature_sector` produces an exact conjugacy-canonical loop-holonomy sector.
`analyze_causal_curvature` consumes before/after sectors, causal event edges, and declared disturbance
sources. It reports reached changes, off-causal changes, maximum causal depth, and a causal-alignment
fraction. The test fixture deliberately contains one propagated and one off-causal change, proving
that the diagnostic exposes rather than hides the latter.

Connection sectors are now attached automatically to actual engine subdivision events, and the
engine's causal edges feed this analyzer. Propagation-speed convergence and foliation comparisons
require a nontrivial derived update law; transport-preserving subdivision is intentionally a
zero-change control.

## 6. Compression starts with an exact obstruction

The normalized fresh-fiber subdivision state is a scaled permutation matrix across its two edge
transports. Its Schmidt spectrum is therefore computed exactly from the stored amplitudes. For
`Aut(C4)`, all eight singular values equal `1/sqrt(8)`: the state has rank eight, and a rank-four
approximation necessarily discards one half of the norm.

This negative result prevents a misleading low-rank speedup claim. Future tensor-network work must
report discarded norm and observable error, and must first search for multi-event structures whose
spectrum actually decays.

## 7. Do not materialize gauge copies on the device path

`subdivide_edge_representative` fixes the new fiber's frame to identity and constructs one physical
child directly. `GaugeSubdivisionEvolution` records `|Aut(F)|` as the raw orbit size but only
registers that representative. Its test compares the result with every exhaustively generated orbit
member under the exact physical state identity. The dense multiplication, inverse, and action arrays
are also exposed as contiguous `uint16` buffers for a device backend.

This removes the most obvious gauge-copy allocation before GPU integration. The official engine's
persistent CUDA kernel still needs fiber-aware matching, affected-cycle updates, canonicalization,
deduplication, and differential tests before an end-to-end performance claim is possible.

## 8. Deterministic small-fiber census

`wgphysics_census` enumerates every labeled simple graph through a chosen small order, quotients them
by exact graph isomorphism, and emits one record per unlabeled fiber. The committed order-four census
contains the known `1 + 2 + 4 + 11 = 18` classes. Each row records:

- graph size and canonical code;
- automorphism-group order and whether it is nonabelian;
- conjugacy classes, equal to the connection sectors on a single triangle;
- the reversible equivariant dynamics count when group order is at most eight;
- raw versus physical subdivision counts;
- exact subdivision Schmidt rank and half-rank discarded norm.

The null dynamics entries for order-24 groups are an explicit computation bound, not zero results.
A companion dataset now crosses every entry and conjugacy sector with the real depth-two subdivision
closure. The next dataset should add other explicitly defined rule morphisms and record nontrivial
causal propagation and localization.
