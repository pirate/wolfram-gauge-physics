# A local gauge seed activates a frozen nonabelian background

Changing one link, while leaving every face charge equal to one, activates
spreading charge fluctuations under the existing primitive bank. In all 36
tested preparations, every face eventually departs from charge one, including
on the largest 1,152-face mesh. The untouched background and reaction-disabled
controls remain fixed. Each changing event has a checked dependency ancestry
back to the initial perturbation.

This is **externally seeded classical activity**, not spontaneous nucleation,
emergent spacetime, a quantum wave, or a localized particle. The experiment
supplies a periodic triangular mesh and independent random operator selections.
It does not insert an additional propagation equation or modify the update bank.

![Actual first-departure map, evolved charges, and boundary-marked spreading curves](images/triangle-spreading.png)

## 1. Preparation and exact first-step response

Use the [nonabelian frozen connection](triangle-reference.md#6-full-nonabelian-holonomy-still-permits-frozen-states):
the positive horizontal, vertical, and diagonal links carry the three different
reflections in the triangle-adjacency-derived group $S_3$. All face charges are
one, all primitive gates are inactive, and the global based-loop image is $S_3$.

At the central vertex, select one of its three positive-direction links and
replace its reflection by either other reflection. This gives six specified
single-link preparations, with no selection based on observed trajectories.
Every triangular face still has a reflection holonomy, hence

$$q_f=1\quad\text{for all }f,\qquad Q=F=2V.$$

The full holonomy image remains $S_3$. This change is **initial preparation**, not
an allowed move out of the frozen state: no primitive can move that state. It
changes microscopic information even though it changes neither total charge
nor the entire charge field.

The exact activity detector finds eight active all-reflection fans for each
preparation at each tested size. With twelve enabled rules per fan, there are
96 changing operator selections out of $M=78V$. For independent uniform
selection, the first-change time $T$, counting attempts from one, obeys

$$\Pr(T=t)=\left(1-\frac{96}{M}\right)^{t-1}\frac{96}{M},\qquad
\mathbb E[T]=\frac{M}{96}=\frac{13V}{16}.$$

Before the first change, every rejected attempt still sees the initial
connection. The script directly predicts the first event from that connection
and the saved schedule seed, then checks it against the compiled trajectory.
The exact mean waiting times are $117/4$, $117$, and $468$ attempts at the three
tested sizes. These are schedule-time predictions, not physical lifetimes.

Each enabled initial update changes $(1,1,1)$ to a permutation of $(0,1,2)$.
Let $C=\sum_f(q_f-1)^2$ be total charge contrast. Uniform bank averaging gives

$$\mathbb E[\Delta q_f\mid U_0]=0\quad\text{for every face},\qquad
\mathbb E[\Delta C\mid U_0]=\frac{192}{M}>0.$$

Exhaustive raw-table enumeration of all six 72-face preparations independently
checks both the zero first-moment vector and the integer contrast numerator 192.
This is a microscopic consequence of the derived reaction gates, not a fitted
Gaussian noise or diffusion term. It is an initial-state identity, not a closed
equation for the subsequent mean field.

## 2. Controlled protocol and actual outcomes

There are two independent schedules at each of sides $L=6,12,24$. All six seeds
at a given side/trial share that schedule. They are paired controls, **not six
independent statistical trials**. The horizons are respectively 12,000,
48,000, and 192,000 attempts, giving the same attempted updates per vertex and
per rooted operator. No-op attempts count toward the clock.

The datasets contain 21 charge/activity snapshots per seeded trajectory, every
face's first charge-departure time, final raw links, and canonical SHA-256
digests of schedules and changing-event streams. The schedule is reproduced
from its integer seed and documented uniform draw rule. Event digests record
identity; they are not a replacement for the actual replay verification.

The 36 trajectories contain 319,207 changing events in total. Every event is
independently replayed from raw links; the original validator also checks each
attempt, boundary links, spectator holonomies, total charge, and full inverse
recovery. The spreading observer then separately replays those verified events.

At the end of every trajectory, **every face has changed charge at least once**.
The maximum visited dual distances from the seed are 7, 15, and 31, equal to the
maximum distances on the respective finite tori. Final charge populations obey
$N_0=N_2$, as required by $Q=F$. On the largest mesh, $N_0=N_2$ ranges from 260
to 289 across the twelve runs, rather than remaining a few localized defects.

The unperturbed background remains fixed under the full bank. Every seeded
reaction-disabled run remains fixed too: initially there are no vacancies,
so vacancy transport alone cannot start. These controls receive the same
attempt schedules as their corresponding active runs.

This demonstrates spreading in the measured finite systems. It does **not**
prove an infinite-volume speed law, mixing to the stationary reference,
a phase transition, or the absence of all possible localized structures.
The particular response tested here does not stay spatially confined.

## 3. Checking ancestry from actual read/write supports

Couple the seeded connection $U_t$ to the fixed background $\bar U$ using the
same schedule, and define the current raw-link difference set

$$D_t=\{e:U_t(e)\ne\bar U(e)\}.$$

For an attempted primitive with read set $R$ and write set $W$, if
$R\cap D_t=\varnothing$, its input agrees with the fixed background. That
primitive must therefore be a no-op. Every changing event must satisfy

$$R\cap D_t\ne\varnothing,\qquad D_{t+1}\subseteq D_t\cup W.$$

The code checks the first condition for every actual event and updates the
difference set from the actual written links. A restored background value is
removed from the difference set; an old difference is not kept alive merely
because an edge was visited in the past. Full-state comparisons at snapshots
check that the bookkeeping matches the raw connections.

For a constructive ancestry certificate, the initial differing link has depth
zero. A currently differing output link receives one plus the largest depth of
the currently differing input links. A directed support graph has an arc from
every possible read link to every possible write link of a primitive. Its
shortest distance $d_{\rm dep}$ from the seed obeys

$$d_{\rm dep}(e)\le\ell(e)$$

for every currently differing link with assigned ancestry depth $\ell(e)$.
The script checks this bound on every write. This is conservative structural
dependency ancestry, not a minimal causal explanation for each output value.
Nor is this support-graph distance Euclidean distance or a physical light cone.
Locality of the supports is supplied by this model; it has not emerged here.

Raw-link equality between two arbitrary, independently reframed connections
would be ambiguous. Here the two runs are in a **common frame**. Applying the
same arbitrary local gauge transformation to both preserves their difference
set. A nonconstant-frame replay checks equal event locations, exactly transformed
final links, identical charge histories, identical first-passage times, and
identical ancestry and seam diagnostics.

## 4. Distances and finite-box safeguards

Spatial distances are shortest paths on the actual dual face graph, starting
from the two faces adjacent to the initially changed link. No layout coordinates
or fitted embedding enter these distances. For each shell, snapshots record its
size, current charge populations, and the count of faces that have ever departed
from charge one. The current contrast radius and cumulative visited radius are
different observables; the cumulative radius cannot diagnose persistence.

The coordinate cut of the supplied periodic mesh identifies seam edges. We
record the first changing event whose **read support** includes a seam edge.
Before that event, the seed's changing-event ancestry has not touched the chosen
periodic cut. Afterward, the plot changes the radius curve to dotted styling.
This is a conservative cut-dependent warning, not an exact time at which every
observable becomes contaminated by periodicity. It is not gauge-frame dependent.

All three sizes eventually fill their periodic boxes, so fitting a power law
to the full radius curves would confuse propagation with finite-size saturation.
The first-departure image is plotted in supplied mesh coordinates and includes
the seed link. It is actual simulation data, not an inferred spacetime embedding.

## 5. Relation to physics, and what still needs testing

Kinetic constraints can produce nontrivial dynamics without complicated
equilibrium interactions; this is an established statistical-physics subject,
reviewed by [Ritort and Sollich](https://arxiv.org/abs/cond-mat/0210382).
Spreading activity alone is therefore not a new general phenomenon. The
model-specific result here is its construction and exact audit through actual
nonabelian link updates, including a seed invisible to the initial charge field,
the exact first-step fluctuation rate, and common-frame dependency checks.

The next useful tests are pre-seam first-passage statistics at larger scales,
directional dependence, and persistence of relational loop or spectral structures
inside the active region. Independent schedule ensembles are needed before
claiming a propagation exponent or a universal macroscopic law. A charge pattern
that spreads is not by itself a molecule, and the current dynamics has no
established quantum amplitudes or electromagnetic interpretation.

## Reproduction

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j4
uv run --with numpy --with scipy --with python-flint python tools/triangle_spreading.py --output out/triangle-spreading.json
diff -u data/triangle-spreading.json out/triangle-spreading.json
uv run --with numpy --with scipy --with python-flint python tests/triangle_spreading.py
uv run --with numpy --with scipy --with matplotlib python tools/plot_triangle_spreading.py
```

A short control run uses `--sides 6 --trials 1 --attempts 1200`. CI reproduces the
full first 72-face case, including its gauge-transformed controls, and checks the
saved spatial observations at every size against first-passage arrays and final
raw connections. CI does not rerun the complete 36-trajectory large-size audit;
the full reproduction command above does.
