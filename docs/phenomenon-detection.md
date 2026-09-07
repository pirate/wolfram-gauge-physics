# From finite gauge dynamics to measurable phenomena

The present code can test finite algebra, explore small rewrite histories, and estimate intrinsic
graph statistics. It has not recovered a physical space, time unit, particle, or molecular state.
The following measurements separate intermediate progress from those later claims.

## Geometry before a molecule template

The export projects each hypergraph to its undirected simple 2-section. This is an explicit
measurement choice: it forgets edge ordering and multiplicity and connects every pair of vertices
in each hyperedge. Report changes under other reasonable projections before claiming universality.

The current probes estimate ball volume $V(r)$ and lazy-walk return probability $P(t)$, then
compute discrete versions of

$$d_H(r)=\frac{d\log V(r)}{d\log r},\qquad
d_s(t)=-2\frac{d\log P(t)}{d\log t}.$$

Diffusion time here counts auxiliary random-walk steps, not physical simulation time. The use of
return scaling as an intrinsic dimension probe has precedent in
[Ambjørn, Jurkiewicz and Loll](https://arxiv.org/abs/hep-th/0505113).

Tests compare a cycle and a square torus without supplying their dimensions to the estimator.
The current plateau selector uses a relative-span heuristic; it supplies neither a statistical
confidence interval nor evidence of a continuum limit. Source selection is deterministic by label
order. Full-source averages are relabel invariant; a truncated sample can change under arbitrary
relabeling. Multi-seed sampling, uncertainty estimates, boundary controls, and increasing graph sizes
are necessary before using these estimates as physical evidence.

To compare eventually with ordinary three-dimensional geometry, infer dimension first and then
look for consistent regimes near three with scale windows that widen as graph size increases.
Measure anisotropy and inhomogeneity separately. A low-distortion 3-D picture of a small patch
does not provide those tests.

## Gauge observables and persistent excitations

The loop probe reports the exact conjugacy sector in $\operatorname{Aut}(F)$ and the normalized
character of its permutation representation,

$$W(C)=\frac{\operatorname{Tr}\operatorname{Hol}(C)}{|V(F)|}.$$

It also reports the fraction of supplied loops with nontrivial holonomy. Defining the diagnostic
weight $w_i=1-W(C_i)$ gives effective occupied support

$$N_{\mathrm{eff}}=\frac{(\sum_i w_i)^2}{\sum_i w_i^2}.$$

This is a participation statistic over the supplied loop collection, not a derived energy or
particle number. Loop choice, duplication, refinement, and representation affect it. Track exact
sectors alongside characters because a single character does not distinguish all sectors.

The new [pair-interaction census](pair-interactions.md) establishes sector changes and constrained
schedule identities. Next, track spatially connected concentrations of these invariant observables
through full multi-cell evolution. Require persistence across many updates, stability under
relabeling and gauge changes, and reproducibility across schedules and initial conditions.

For an invariant cell observable $O$, a candidate connected correlation statistic is

$$C_O(r)=\langle O(x)O(y)\rangle_{d(x,y)=r}-\langle O\rangle^2.$$

Define the state ensemble and distance sampling before interpreting its decay. A finite correlation
length alone does not establish a quantum mass gap. Loop area/perimeter scaling likewise requires
ensembles, explicit composite loops and surfaces, and finite-size checks; current length bins
alone are insufficient.

## Dependency, spatial locality, and quantum correlation

An event dependency records data read from an earlier event. The interaction API exposes read/write
sets, so off-support commutation can be tested exactly. Spatial locality additionally requires
bounded interaction extent in the spatial graph. Long connector paths must be recorded, not
silently treated as one short-range physical interaction.

A propagation experiment compares matched runs differing by one specified local connection
perturbation, tracks invariant changes against graph distance and causal depth, and repeats under
different compatible schedules. The clock and coarse-graining prescription must be specified.
Demonstrating stable signaling cones is a further result beyond observing dependency edges.

Branchial adjacency alone supplies neither amplitudes nor an entanglement measure. Quantitative
entanglement requires a specified state space, state, subsystem split or observable algebras, and
a suitable correlation or reduced-state calculation. Entanglement correlations do not by themselves
constitute a signaling channel.

## Admission criteria for the molecular stage

1. A reproducible large-scale geometric regime, tested independently of its display coordinates.
2. Persistent localized excitations and derived conserved charge sectors.
3. Measured interactions between excitations: response, dispersion, and any long-range potential
   must come from the update law rather than an inserted Coulomb term.
4. Stable neutral bound configurations, with a defined clock, excitation spectrum, and scale
   calibration. A periodic orbit of two finite group elements is insufficient.
5. Only then compare recovered distances and spectra with known molecular structure. Use those
   comparisons as held-out tests, not as a shape objective fed into the microscopic updates.

Null controls include flat connections, frame-only changes, transport-preserving subdivision,
relabelings, and disconnected independent patches. Predeclare detection windows and failure
criteria; keep unsuccessful runs, perform finite-size scaling, and verify on held-out rules.

`tools/novelty_sweep.py` currently explores a bounded 38-rule family by graph-statistic diversity.
The checked-in `data/novelty-sweep-v1.json` is a three-step census, with no target geometry supplied.
Its averages over scale-dependent dimension estimates are exploratory features, not inferred
continuum dimensions. It does not yet evolve the new gauge interactions or constitute an
unbounded autonomous research loop.

```bash
python3 tools/novelty_sweep.py --steps 3 --keep 8 --output out/novelty-sweep.json
```
