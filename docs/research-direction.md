# Research direction: unresolved questions from the interview

The objective is to make mathematical progress on the program described in
[Stephen Wolfram's World Science Festival interview](https://www.youtube.com/watch?v=yAJTctpzp5w):
start from discrete relationships and local replacements, then determine
what geometric and physical structures follow. The software is an
experimental instrument for those questions, not the product being developed.

The interview motivates a program; it is not a proof of its proposed
physical correspondences. Results in this repository must identify the
microscopic assumptions and the scope of their conclusions separately.

## Questions, in working order

### 1. Can a rewrite system supply its own internal structure?

Wolfram's [local gauge-invariance proposal](https://www.wolframphysics.org/technical-introduction/potential-relation-to-physics/local-gauge-invariance/)
relates gauge freedom to equivalent local rewrite choices and relates
field-like effects to their subsequent causal consequences.
[InfraGaugeTheory](https://github.com/WolframInstitute/InfraGaugeTheory)
lists obtaining fibers from hypergraph rewriting and natural clustering
into fibers among its goals. This is the missing construction to address,
not an instruction to pick a finite group with desirable properties.

The first tractable target is a completely explored finite example:

1. Specify the raw rewrite rule, seed, event semantics, and exploration bound.
2. Define what a proposed observer retains. This projection is an explicit
   assumption unless independently derived; it cannot be selected afterward
   just to produce a desired group.
3. Identify detailed alternatives over the same retained description.
4. Derive which actual transitions relate them. Record where maps are
   absent, partial, nonunique, or incompatible across overlaps.
5. Separate arbitrary relabelings from alternatives with different future
   invariant consequences. A coarse observer's inability to distinguish
   states does not make those states gauge copies.

For a candidate projection, the set-theoretic fiber is

```math
F_x=p^{-1}(x).
```

That formula does not define fiber adjacency, a group action, or transport.
Each must be supplied by an explicit construction and checked against raw
histories. Candidate transports may be relations rather than bijections;
do not replace them with permutations merely to fit the existing library.

For a deterministic raw update, an exact reduced update exists only if
the proposed reduction is consistent:

```math
p(s)=p(s')\ \Longrightarrow\ p(Ts)=p(Ts').
```

For branching evolution, compare the projected successor structure under
the declared semantics instead. Preserve multiplicities when those are
part of the question. A failure can identify information the observer
must retain; it is not permission to invent a closure law.

**Useful outcome:** a finite construction with explicit projection,
transport, and equivalence conditions; or a counterexample showing why a
particular proposed reduction cannot support them. Neither outcome alone
establishes electromagnetism.

### 2. How does a local choice change later possibilities?

Start with two valid alternatives from the same raw state. Record their
consumed and produced relationships and any additional read supports.
Compare later admissible matches and invariant measurements, not node
positions in the viewer.

Separate four possibilities: pure renaming, independent updates that
can be exchanged, histories that merge, and persistent differences in
later behavior. Common final states alone do not prove equality of causal
histories or causal invariance of the rule.

**Useful outcome:** an explicit local-choice effect and its dependency
support, or a proof that a proposed effect disappears under the correct
equivalence. To interpret it as a field requires additional spatial,
scaling, and dynamical results.

### 3. Which rules sustain a geometric regime?

The interview identifies [dimension fluctuations](https://www.youtube.com/watch?v=yAJTctpzp5w&t=4277s)
as a direction to investigate. First establish whether a dimension can
be estimated consistently at all in the chosen evolution.

Extend bounded rule searches using multiple seeds, larger graphs, and
declared update schedules. Compare graph-ball growth and return probes
across scale. Test source sampling and the loss of information introduced
by projecting a hypergraph to a simple graph. Do not interpret a short
plateau or a three-coordinate drawing as a three-dimensional continuum.

**Useful outcome:** a growing range of stable geometric scaling, a
characterization of genuine spatial/temporal variation, or a restriction
showing that a rule family cannot sustain the proposed regime.
Matching a target dimension during search is selection, not a prediction;
any such selection must be stated.

### 4. Can a structure persist and move through rewriting?

The interview leaves a [concrete particle model](https://www.youtube.com/watch?v=yAJTctpzp5w&t=5138s)
unresolved. Search for invariant structure with a traceable history, not
a named particle or a familiar rendered shape.

Distinguish persistence through change from an untouched region. Check
whether a candidate survives replacement of its node IDs, moves relative
to surrounding structure, and remains identifiable after an encounter.
Measure lifetimes relative to background activity under specified clocks.

**Useful outcome:** a reproducible persistent or propagating structure,
its invariant descriptor and failure conditions; or a no-go result for
a specified candidate family. Localized eigenvectors, frozen defects,
and finite periodic orbits alone do not establish particles or binding.

## How existing finite-fiber work fits

The fixed-mesh studies choose the spatial mesh, fiber adjacency, rule
bank, and stochastic scheduler. Calculating consequences of those choices
can be useful, but does not derive the choices from bare rewriting.

- Exact gauge quotients and patch gluing help test whether a proposed
  description loses predictive information or double-counts equivalent states.
- Constraint and encounter calculations can test whether internal
  relationships persist through activity, rather than merely remaining untouched.
- Localization and equilibrium restrictions test particular candidates
  for persistence and binding.
- The diffuse-refinement limit excludes one proposed route to sustained
  reactions on the stated time scales. Its heat equation is a known
  limiting process, not the project's physical target.

These are supporting laboratories. Further work needs an explicit link
to one of the questions above or another named unresolved suggestion
from the interview. A finite-model theorem may be useful without being
a derivation of nature; that boundary belongs in its statement.

## What should drive a compute run

Before a substantial run, write down the question, assumptions, measured
quantity, and what outcome would change the next research decision.
Choose the budget to resolve that uncertainty: finite enumeration for a
bounded classification, long histories for persistence, or multiple
sizes and seeds for a scaling claim. Report incomplete exploration as
incomplete, not as absence of a phenomenon.

Keep raw witnesses and counterexamples. Use known physics to evaluate
an independently produced result, not to supply hidden forces or shapes.
Prepared initial conditions are legitimate controlled experiments when
identified as such; they do not demonstrate spontaneous formation.

Deprioritize repeated confirmations of a known diffusion limit, larger
fixed-fiber parameter sweeps without a discriminating question, cosmetic
viewer work, and standalone kernel speed records. Performance work is
in scope when it makes a named mathematical experiment feasible.

The next milestone need not be a molecule. A consistent rewrite-derived
fiber construction, a precise obstruction to one, a robust geometric
regime, or a persistent structure would each answer a nearer question.
Quantum probability, physical charge, mass, and unit calibration remain
separate obligations before molecular claims.
