# Research map and current boundaries

The project is organized around open questions in
[Wolfram's interview](https://www.youtube.com/watch?v=yAJTctpzp5w), not a
claim to have advanced every layer of gauge theory. The
[research direction](research-direction.md) specifies the questions and
what would count as an answer. The [README](../README.md) explains the
conceptual sequence.

## Rewrite-based work

[HypergraphRewritingEngine](https://github.com/WolframInstitute/HypergraphRewritingEngine)
is the upstream engine for raw rewrite histories and their causal and
multiway records. The local [product construction](product-evolution.md)
adds a supplied finite-fiber connection and jointly identifies equivalent
base/connection states. Its supported base-event morphology is edge
subdivision. It is not a general derivation of gauge interactions from
arbitrary hypergraph rewrites.

[Intrinsic geometry probes](phenomenon-detection.md) measure graph-ball
growth and random-walk returns. The bounded rule survey is exploratory;
its small graphs and short histories do not establish emergent physical
dimension. Projection, source sampling, and finite-size limitations are
part of the measurement.

These are available starting points for investigating geometry, local
rewrite alternatives, and event dependencies. A debugger or a faster
kernel is useful only insofar as it makes those experiments possible.

## The missing gauge connection

Wolfram's [local gauge-invariance proposal](https://www.wolframphysics.org/technical-introduction/potential-relation-to-physics/local-gauge-invariance/)
suggests investigating symmetry among local rewrite choices and the
causal consequences of selecting them.
[InfraGaugeTheory](https://github.com/WolframInstitute/InfraGaugeTheory)
provides combinatorial fiber and connection constructions and identifies
obtaining fibered graphs from rewriting as a goal.

Our finite-fiber machinery starts with a specified internal graph.
Its automorphisms and gauge quotients follow from that graph, but the
choice of graph has not been derived from the raw rewriting process.
The [foundations note](infragauge-foundations.md) states those inputs.

The next construction must specify a projection, determine its internal
alternatives, and recover consistent transport from actual events—or
identify why that construction fails. It must not silently replace
general partial relations with a desired finite group.

## Supporting finite-model results

### Equivalent descriptions and lost information

[Complete loop observables](complete-loop-observer.md) and
[triangle patch gluing](triangle-patch-observer.md) give exact descriptions
for specified finite gauge systems. The patch results show why separate
local summaries may lose relative alignment needed to predict their
interaction. This is a concrete reference for testing proposed reduced
descriptions, not a derivation of a physical observer.

### What local rules can preserve or change

The [unary census](causal-dynamics.md) separates apparent motion within a
gauge orbit from changes in invariant loop sectors. The
[shared-edge](shared-edge-transport.md) and
[three-face](three-face-feedback.md) constructions study particular
boundary-preserving interaction laws. The
[cycle reaction formulas](cycle-relational-dynamics.md) make the chosen
rule bank and its conserved weight explicit.

These results characterize selected finite dynamics. Reversibility,
gauge covariance, and a positive conserved weight do not uniquely select
a physical law. The original positive-charge criterion is a selection
assumption, not an unbiased prediction of charge.

### Internal relationships, encounters, and memory

[Relative-angle dependence](fiber-relative-angle.md),
[reaction bursts](fiber-reaction-bursts.md), and
[transported constraints](fiber-constraint-dynamics.md) study how internal
relations affect subsequent updates on a fixed mesh.
[Encounter-resolved measurements](triangle-encounter-memory.md) distinguish
untouched regions from states that change and return.

These are useful for asking whether candidate structure survives
activity. A finite-mesh memory effect is not automatically a particle,
a force, or a quantum correlation.

### Restrictions on persistence and binding

The [diffuse-refinement result](fiber-refinement-limit.md) bounds a
reduction of the actual charge process to known exclusion dynamics in a
specified regime. It identifies the disappearance of reactions in that
limit; it is not a discovery of diffusion or a general no-go theorem for
matter.

The [localization construction](triangle-modes.md) certifies modes for a
prepared defect and bounds how allowed dynamics can destroy them.
[Reachable equilibrium](reachable-gauge-equilibrium.md) excludes a
separation preference in a particular finite model. Together these
provide reasons not to identify static localization or clustering with
persistent bound matter.

## Boundaries that remain open

- A fiber or internal state structure obtained from underlying rewrites,
  rather than supplied at each base node.
- Consistent interacting geometry and internal dynamics beyond the
  restricted subdivision construction.
- Reproducible geometric scaling on growing evolved networks.
- Persistent structures that move and survive encounters under those
  rewrites.
- A construction of quantum amplitudes, interference, and probabilities.
  The existing normalized subdivision-orbit map is only a kinematic
  isometry, not a quantum dynamics.
- Identification of physical charge, energy, mass, and length/time scales.
- Predictions for bound systems that were not put into the microscopic
  rules or initial geometry.

No literature-priority claim follows from these local calculations.
Their value is the specific construction, counterexample, or restriction
they establish under stated assumptions. Scripts and datasets remain
available through the linked notes; they are not a list of completed
steps toward the Standard Model.
