# Draft Wolfram Community post

## From local rewrite choices to fibers: looking for a finite construction

This project began with
[Stephen Wolfram's World Science Festival interview](https://www.youtube.com/watch?v=yAJTctpzp5w),
particularly its unresolved questions about dimension, electric charge,
and persistent particle-like structures.

The gauge question I would like to make concrete is the route from
actual local rewrite alternatives to internal structure.
[Local Gauge Invariance](https://www.wolframphysics.org/technical-introduction/potential-relation-to-physics/local-gauge-invariance/)
suggests examining equivalent local choices and the causal consequences
of making them. [InfraGaugeTheory](https://github.com/WolframInstitute/InfraGaugeTheory)
provides relevant graph-fiber constructions and lists obtaining fibers
from rewriting as a goal.

The repository currently has two distinct kinds of work:

- Raw histories and causal/multiway records from
  [HypergraphRewritingEngine](https://github.com/WolframInstitute/HypergraphRewritingEngine),
  with intrinsic graph probes and a restricted connection-aware
  edge-subdivision extension.
- Exact finite-gauge laboratories that **supply** a fiber graph and
  interactions, then calculate their symmetry, invariant state
  descriptions, and dynamics.

The second layer does not derive the supplied fiber from the first.
Calculating a chosen graph's automorphisms does not resolve that gap.
The existing normalized subdivision-orbit construction is also only
kinematic; it is not a quantum evolution law.

Some finite results are useful controls. Local patch descriptions can
lose relative alignment needed to predict interactions. A diffuse
large-cycle limit loses reactions and reduces to known exchange
diffusion. A prepared defect can support localized graph modes that
allowed updates subsequently destroy. These results help reject
particular interpretations, but do not establish fields or particles.

The next target is a small, completely explored rewrite example in
which we can specify a projection, identify candidate internal
alternatives, and derive the relations between them from actual
events. Partial or inconsistent transport would be an informative
outcome rather than something to replace with a preferred group.

Questions for people working on this construction:

1. Which explicit rule and finite seed provide a useful first example
   of internal alternatives derived from rewriting, rather than
   attached as extra state?
2. What criterion separates gauge-equivalent descriptions from
   physically distinct branches that merely look the same to a
   coarse observer?
3. Can transport between those alternatives be obtained from event
   provenance, including where lifts are partial or nonunique?
4. What finite counterexample or consistency result would most help
   clarify the proposed connection to gauge-field propagation?

The [visual introduction](https://github.com/pirate/wolfram-gauge-physics#readme)
separates primitive rewrites, causal histories, supplied-fiber examples,
and the remaining research questions. The
[research direction](https://github.com/pirate/wolfram-gauge-physics/blob/main/docs/research-direction.md)
states what new experiments should resolve. Scripts and saved witnesses
are linked from the mathematical notes.

This is an independent experimental project. No physical charge,
particle model, quantum probabilities, or molecular dynamics have
been derived.
