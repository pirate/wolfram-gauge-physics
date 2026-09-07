# State of the art and contribution boundary

This note records the public work used to position the repository. It is not a comprehensive
history of discrete gauge theory; its scope is gauge/fiber computation around the Wolfram Physics
Project.

## Wolfram-model motivation

The Wolfram Physics technical introduction proposes that a spatial hypergraph acts as a fiber
bundle base, while branchlike choices in the multiway causal graph act like gauge choices. It
suggests that local phenomena may reflect an effective Lie algebra and that causal consequences
encode gauge-field propagation:

- [Local Gauge Invariance](https://www.wolframphysics.org/technical-introduction/potential-relation-to-physics/local-gauge-invariance/)
- [Gauge Groups Meet Hypergraphs](https://writings.stephenwolfram.com/2020/07/a-burst-of-physics-progress-at-the-2020-wolfram-summer-school/#gauge-groups-meet-hypergraphs)

That conceptual proposal does not by itself specify a production data structure for connections
that evolve with hypergraph rules.

## Public computational foundations

### InfraGaugeTheory

[WolframInstitute/InfraGaugeTheory](https://github.com/WolframInstitute/InfraGaugeTheory) develops
fibered graphs, graph projections, bundle predicates, sections, connection subgraphs, horizontal
lifts, parallel transport, holonomy matrices, flatness, and horizontal leaves. It explicitly lists
obtaining fibered graphs from hypergraph rewriting and gauge dynamics in the Wolfram Physics
Project among its goals.

Relative to the [inspected InfraGaugeTheory revision](https://github.com/WolframInstitute/InfraGaugeTheory/tree/edd9bdca46b7838d6b3e940e8ae8cde90b60ef2c), this repository adds
an exact fiber-automorphism group, explicit local-frame action on connections, a spanning-forest
gauge quotient, rewrite transport factorization, quantum amplitudes over rewrite gauge orbits, a
compact device representation, and an exact connection-aware product evolution over real engine
event provenance.

### HypergraphRewritingEngine

[WolframInstitute/HypergraphRewritingEngine](https://github.com/WolframInstitute/HypergraphRewritingEngine)
provides exact canonicalization, multiway/causal/branchial evolution, quotient exploration,
incremental matching, and a persistent CUDA backend. This repository uses pinned commit
[this engine revision](https://github.com/WolframInstitute/HypergraphRewritingEngine/tree/03fe60ddf338983060b6bb4b23e8b4b5d7ae7337) as its base-space evolution engine.

The product runner now attaches connection sectors to every raw state and supported subdivision
event from that engine, then applies a joint base/gauge quotient. Connection data is not yet part of
the engine's internal match key, canonical state, event identity, or CUDA storage. The runner must
therefore retain full raw provenance for exactness. The integration design is in `gpu-execution.md`.

### Infrageometry projects

[Infrageometry](https://github.com/WolframInstitute/Infrageometry) and
[SyntheticInfrageometry](https://github.com/WolframInstitute/SyntheticInfrageometry) investigate
which geometric structures can be built or observed on discrete graphs at different information
ceilings. They motivate keeping the fiber graph and observer-accessible invariants explicit rather
than assuming continuum coordinates.

## Selected community investigations

- Graham Van Goffrier, [Full Discretization of Fiber Bundle Topology for Gauge Theory](https://community.wolfram.com/groups/-/m/t/2030337), Wolfram Summer School 2020.
- Chang Wu, [Exploring Gauge Symmetries in the Wolfram Model](https://community.wolfram.com/groups/-/m/t/2162318), Wolfram Winter School 2021.
- Matthew Maddock, [Gauge Field Theories in Terms of a Discrete Principal Fibre Bundle](https://community.wolfram.com/groups/-/m/t/2312018), Wolfram Summer School 2021.
- Omar Medina, [SU(2) Gauge Theory in the Wolfram Model](https://community.wolfram.com/groups/-/m/t/2163358), Wolfram Winter School 2021.
- Ioana-Alexandra Milea, [An Investigation of Discrete SU(2) Gauge Theory through the Hopf Fibration and Wilson Loops](https://community.wolfram.com/groups/-/m/t/3497643), Wolfram Summer School 2025.
- [The Fine-Structure Constant Challenge](https://community.wolfram.com/groups/-/m/t/2131169), whose discussion explicitly identifies constructing a specific gauge field from a Wolfram model as an open target.

These works explore important discrete gauge constructions and physical interpretations. This
repository's narrower contribution is an executable bridge between combinatorial fibers and exact
rewrite evolution, with strict treatment of gauge copies.

## New finite-model results

The [pair-gate classification](braid-quotient-classification.md) and
[exact diffusion bridge](emergent-diffusion.md) identify known exclusion-process dynamics
inside the independent-cell construction. These are finite classifications and validation,
not a claim to have discovered diffusion or generalized permutation gates.

The [shared-face follow-up](shared-face-dynamics.md) supplies a fixed-boundary counterexample
to extending that closed sector description to an actual mesh. It derives the failure of all
nontrivial rational one-face additive class charges for the selected local law, implements
full incident-face updates, and observes reversible coarse relaxation near an exact
stationary uniform-link reference. This is a concrete improvement in the fidelity and
testability of our gauge/fiber simulator, not evidence of a continuum gauge theory or matter.

The [one-face energy obstruction](face-energy-obstruction.md) generalizes the selected-rule
failure to a finite-group normal-subgroup criterion and excludes every nonconstant one-face
energy for all 5,730 nonidentity strict rules in our square-fiber census. Surviving quotient
observables in other controls are fixed locally, not transported. The accompanying
[observable audit](observable-memory-and-invariants.md) measures exact coarse-memory errors,
frozen sectors and schedule dependence, and eliminates a specified neighboring-loop density
ansatz with larger-mesh counterexamples and a universal constant-function certificate.
These are scoped mathematical results with executable checks; literature priority and a
connection to physical energy remain unestablished.

The [shared-edge construction](shared-edge-transport.md) gives a constructive escape specific
to a different link realization: the same pair table now preserves all exterior links, affects
only its two input faces, and transports the derived additive charges globally. Its two-state
class memory, bipartite parity, and 4.92 million independently replayed link updates are exact
finite-model results. The autonomous exclusion factor remains an explicit limitation; this is
not yet a microscopic derivation of interacting matter or physical energy.

## Contribution matrix

| Layer | Public upstream capability | Implemented here | Still open |
|---|---|---|---|
| Fiber | Total graph, projection, fiber predicates | Exact `Aut(F)` derivation, compact tables, and exact open-link inference | Compare rule, link, and multiway-derived fibers |
| Connection | Connection subgraph and horizontal lifts | Explicit transport maps, local-frame action, non-isomorphic fibers, and partial lifts | Twists, hypergraph fibers, dynamic fiber type |
| Curvature | Holonomy matrices and flatness | Conjugacy observables and exact gauge signature | Dynamical curvature action from rewrites |
| Rewriting | Exact bare-hypergraph multiway evolution | Strict event-provenance subdivision, explicit oriented faces, and joint physical state identity | Carry cell incidence in upstream arena and add proven morphisms |
| Quantum | Multiway structure; separate community work | Normalized fresh-fiber orbit isometry and exact Schmidt bound | Multi-event complex amplitudes and interference |
| Performance | Persistent CUDA rewrite engine | `uint16` group tables, physical-only subdivision, and verified Metal probe | End-to-end fiber-aware CUDA integration |
| Physics | Conceptual gauge emergence program | Exact independent-cell diffusion factor, shared-face multi-loop feedback, and reversible coarse relaxation | Derived conserved energy, correlation scaling, engine-integrated geometry changes, continuum limit, matter |
| Reproducibility | Public demonstrations and technical documents | Finite fiber/rule censuses, complete pair classification, exhaustive shared-face patch oracle, and controlled mesh runs | Broader rules, independent replication, multi-face invariants and localization census |

## Claim discipline

The repository does not claim:

- that `Aut(F)` is uniquely the gauge group of a Wolfram model;
- that a continuous Lie group has emerged;
- that branch counts alone define quantum amplitudes;
- that a gauge-invariant kinematics specifies a force law;
- that nuclei, atoms, molecules, photons, or fermions have been derived;
- that a microkernel throughput number is an end-to-end simulation speedup.

Each such statement is instead represented as a testable research milestone.
