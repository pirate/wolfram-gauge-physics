# Sparse encounters, kinetic sectors, and an exact rate benchmark

The absence of sparse-seed reactions in the original repeating schedules was **not** a
no-interaction theorem. The same two-link seed and unchanged three-face rule admit encounters.
An exact reachable-state search finds shortest paths; the independent link oracle and compiled
two-spoke kernel both replay them and restore the original raw links in reverse.

The stronger finding is a kinetic restriction: the four-unit sector splits into a reactive
component and a central-only component on the two completely enumerated meshes. Additive
charges alone do not identify the accessible states. A stationary counting argument then gives
an exact, strongly dilute-suppressed conversion rate for a specified random scheduling ensemble.
Neither that ensemble nor the supplied triangulation is claimed to emerge from the model.

## 1. Derive the reduction from the links

The [selected rule](three-face-feedback.md) and the opposite-reflection two-link seed preserve

$$K=\{1,r,r',z\}\cong C_2\times C_2,\qquad r'=rz.$$

This follows from subgroup closure, the table's restriction to $K^3$, and the two-spoke
reconstruction formulas: every new link is a product of existing subgroup links and target
holonomies in the same subgroup. Because these transports commute, changing the root of an
oriented face does not change its holonomy. Consequently the three face holonomies form an
**exact closed factor** of these link trajectories. No force or substitute face evolution is
inserted. Tests exhaust all $4^7=16,384$ assignments on the seven read links of a fan, including
the surrounding faces, against the independent raw-link oracle.

The factor retains the distinction between $r$ and $r'$, which individual $D_4$ face classes
erase. It uses a chosen subgroup frame; its state counts are **not counts of full physical
gauge orbits**, nor a proof that general noncommuting link configurations admit this reduction.
The encounter path is separately checked after independent vertex-frame changes in the full
$D_4$ group.

Formal face fields are not accepted without a link realization. On this oriented torus, a dual
spanning-tree elimination constructs a $K$-valued link preimage for any face assignment whose
total group product is identity. It fixes one leaf-face target at a time by changing the edge
to its parent. The final root condition is exactly the total-product constraint. The tool
implements this construction and checks the resulting holonomies.

## 2. Exhaustive components within one sparse charge sector

Let each dual sublattice contain $m=L^2$ faces. The two-link seed has, on **each** part,

$$Q_r=2,\qquad \prod_f H_f=z,\qquad Q_2=Q_3=0.$$

The middle face never changes. On the outer pair, both weight and product are conserved; these
claims are checked on every subgroup triple. Thus each part contains either one $z$, or one $r$
and one $r'$ at different sites. There are $m+m(m-1)=m^2$ possibilities per part, or $m^4$ in total.
All satisfy the global identity-product condition and therefore have link preimages.

The $m^2$ states containing two $z$ faces and no reflections form a closed set: the rule can
transport them, but splitting a $z$ requires a neighboring reflection. Reversibility also
prevents a reflection-containing state from entering that set. The remaining $m^4-m^2$ states
have at most one central face. This separation holds at any supported size; connectivity
**within** those sets is established by complete search only at $L=3,4$:

| Supplied torus | Reactive component | Central-only component | Total compatible face states |
| --- | ---: | ---: | ---: |
| $L=3$, 18 faces | 6,480 | 81 | 6,561 |
| $L=4$, 32 faces | 65,280 | 256 | 65,536 |

Every changing update is explored; applying it twice must return the previous state. Each
new state is checked against its sublattice charge/product signature. Inert updates can be
omitted from reachability because unchanged target holonomies leave the raw links unchanged.
A sparse face-to-patch index skips only all-flat triples, which are proven inert. State-budget
termination is explicitly **inconclusive**, never an unreachability result.

From the specified opposite-reflection link pair, the shortest encounter takes three changing
events at $L=3,4$, and five at $L=12$. At $L=3$ one path is

$$
(1,1,r)\to(r,1,1),\qquad
(1,1,r)\to(r,1,1),\qquad
(r',r',r)\to(1,r',z),
$$

on successive, overlapping fans. These are triples at the selected event, not successive
values of one fixed fan. The first two events are transport; the last is the first conversion.
This is an existence and shortest-path result under arbitrary sequential scheduling, not the
probability of a scattering event or a naturally selected clock.

## 3. An exact stationary prediction, with an explicit scheduling assumption

Choose one of the $6m$ fan operators independently and uniformly at every attempt, retaining
idle attempts. Each operator is an involution, so the Markov transition matrix is symmetric:
the uniform measure on each closed component is stationary. Complete connectivity and positive
self-loop counts prove a unique stationary distribution in each of the two small-mesh
components. No convergence time is inferred for the larger meshes.

In the reflection-containing set, the counts at $N_z=0,1$ are respectively

$$m^2(m-1)^2,\qquad 2m^2(m-1).$$

Therefore, under the uniform measure on that set,

$$\langle N_z\rangle=\frac{2}{m+1}.$$

For each fan, each of its eight active collision triples leaves $m-1$ choices for the other
reflection on the middle face's sublattice. There are thus $48m(m-1)$ directed collision events
in the entire set. Dividing by the state count and the number of attempted operators gives

$$
p_{\rm conversion/attempt}=\frac{8}{m^2(m+1)},\qquad
\mathbb E[C_{\rm sweep}]=\frac{48}{m(m+1)}.
$$

Here a *sweep* means $6m$ independent attempted updates, not a complete permutation and not
physical time. Exhaustive graph counts give 3,456 and 11,520 directed collision events at
$L=3,4$, exactly matching the formula. The stationary measures and counts are well-defined at
larger sizes too, but convergence to that measure from the seed is not established there.

With fixed total weight four, conversion is suppressed as $m^{-2}$ per sweep. That is a
consequence of the local catalytic rule and available configurations, **not** a derivation of
a physical cross section. Omitting idle attempts from occupation statistics would instead
sample a degree-weighted jump chain and invalidate this comparison.

## 4. Random-schedule controls and limits

The checked runs use the same two-link initial condition as the reachability search, without
an inserted collision prefix. Every run has a matched transport-only control receiving the
same sequence of uniformly sampled attempted patches. The exact factor, independent raw-link
oracle, and compiled link kernel agree on all changing events and final links; both link
implementations verify exact inversion. Idle-event omission in the raw replays uses the
separately proven identity property. Recorded event times retain the attempted-update clock.
Together the forty combined/control pairs cover 43,872,000 attempted updates and 1,138,770
changing events; every changing event is independently replayed by both link implementations.

- $L=3$: eight runs of 1,000 sweeps, 486–571 conversions per run.
- $L=4$: eight runs of 1,000 sweeps, 146–180 conversions per run.
- $L=12$: sixteen runs of 500 sweeps; four runs react, with 2–4 conversions each.
- $L=24$: eight runs of 500 sweeps; one run has one conversion, the other seven none.

These finite runs do not prove stationarity, a Poisson reaction law, or a continuum limit.
At $L=24$, breadth-first search hit its 200,000-state budget without a collision witness;
the independently sampled trajectory nevertheless supplies one. This is a concrete example
of why bounded negative searches must not be reported as impossibility results.

The two-central-defect sector is not frozen: its defects can move. It is *reaction-inactive*.
Likewise, a central face that persists because a catalyst is absent is not evidence of a
bound molecule. There is no attractive force, two-body binding energy, quantum amplitude,
emergent three-dimensional geometry, or schedule-invariant spacetime here.

Restricted transitions producing nontrivial dynamics are well studied in
[kinetically constrained models](https://arxiv.org/abs/cond-mat/0210382); neither that principle
nor detailed balance under symmetric updates is claimed as new. The contribution here is
the explicit derivation from this finite gauge-link rule, its complete small-sector
classification, and independently replayable sparse encounters. Literature priority for the
particular construction remains unestablished.

The next rule comparison should measure activation requirements, additional kinetic sectors,
and dilute encounter rates across the symmetry-derived census. Selecting a force law or
declaring a persistent label to be a particle would bypass the unresolved physics.

## Reproduce

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
python3 tests/three_face_reachability.py
python3 tools/three_face_reachability.py --side 3 --exhaustive --random-trials 8 --sweeps 1000 --output out/reachability-3.json
diff -u data/d4-three-face-reachability-3.json out/reachability-3.json
python3 tools/three_face_reachability.py --side 4 --exhaustive --random-trials 8 --sweeps 1000 --output out/reachability-4.json
diff -u data/d4-three-face-reachability-4.json out/reachability-4.json
python3 tools/three_face_reachability.py --side 12 --max-states 200000 --random-trials 16 --sweeps 500 --output out/encounters-12.json
diff -u data/d4-three-face-encounters-12.json out/encounters-12.json
python3 tools/three_face_reachability.py --side 24 --max-states 200000 --random-trials 8 --sweeps 500 --seed 6930917 --output out/encounters-24.json
diff -u data/d4-three-face-encounters-24.json out/encounters-24.json
```
