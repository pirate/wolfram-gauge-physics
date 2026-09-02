# Draft Wolfram Community post

## Deriving finite gauge structure from fiber automorphisms and carrying it through Wolfram-model rewrites

I have been experimenting with a small C++ bridge between the Wolfram Institute's
`InfraGaugeTheory` and `HypergraphRewritingEngine` repositories.

The starting point is a finite graph fiber `F`, without declaring `U(1)`, `SU(2)`, or a particle
model. The code derives the local symmetry group as `Aut(F)`, represents connections as fiber
isomorphisms on base edges, and computes parallel transport, holonomy, local frame changes, flatness
obstructions, and horizontal sections.

The part I most wanted to test was what happens when the base graph itself rewrites. For the first
elementary move `(x--y) -> (x--w--y)`, the old transport is factored through the fresh fiber in
every possible way. There are `|Aut(F)|` factorizations, but the implementation verifies that they
form one gauge orbit under a frame change at `w`, so they should not become falsely distinct
physical multiway branches. The corresponding quantum map is a normalized superposition over that
orbit.

There is also an exact connection quotient: gauge-fix a spanning forest, retain fundamental-cycle
holonomies, then canonicalize their simultaneous conjugacy class. For a square fiber the code
derives the nonabelian eight-element dihedral group and tests flat versus curved sectors without a
named continuum gauge group.

A second pass now implements the eight listed research milestones as finite tests: exact fiber
candidate inference from vertex links, gauge-aware rewrite-state identity, non-isomorphic fibers
with partial lifts, an exhaustive census of reversible conjugation-equivariant local updates, and a
causal-curvature reachability score. The first compression result is negative but useful: the
elementary `D4` subdivision state has a flat rank-eight Schmidt spectrum, so halving the rank must
discard half the norm. The execution path now constructs one gauge-fixed physical child without
materializing the eight frame copies, and a reproducible dataset covers all 18 unlabeled simple
fiber graphs through four vertices.

The repeated group operations compile to small integer lookup tables. A differential MLX/Metal
probe on an M1 Max reached about 3.15 billion local-frame transforms/s and 24.1 billion group
multiplications/s for length-eight paths, though this is deliberately reported only as a kernel
microbenchmark; exact multiway state canonicalization remains the harder performance problem.

Repository: https://github.com/pirate/wolfram-gauge-physics

I would especially appreciate feedback on these questions:

1. Should the microscopic local group be `Aut(F)`, or a stabilizer/quotient defined jointly by the
   fiber and allowed rewrite rules?
2. What is the most natural route for extracting candidate fibers from rule automorphisms or
   branchlike structure, rather than supplying `F` explicitly?
3. Should a rewrite-aware extension align with InfraGaugeTheory's total-graph/projection model,
   its connection-subgraph model, or explicit transport maps?
4. Which early observable would be most informative: holonomy propagation, center sectors,
   loop-scaling behavior, or compatibility with multiway causal invariance?

This is not presented as a derivation of electromagnetism or the Standard Model. It is intended as
a small, reproducible set of combinatorial primitives and exact tests that may help make that
research question more computationally concrete.
