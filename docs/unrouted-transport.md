# Gauge-memory detection without prescribed circuits

The [controlled triangle-fiber experiment](noncommuting-fiber-transport.md) routes
defects along chosen loops to identify their exact action. This follow-up removes that
routing: at each attempted tick, a seeded generator samples a rooted adjacent-face patch
uniformly, independently of the state. Only the unchanged boundary-fixed vacancy rule
updates the links. No defect is pinned, no path target is supplied, and no observable
feeds back into the update selection.

A detector finds histories that return all four **tagged** defects to the same placement
while changing the complete gauge state. Tags are passive provenance: a vacancy move
transfers a tag to the newly occupied face. They never enter the rule or scheduler.
The result is unrouted classical gauge memory, not autonomous binding, a continuum
limit, or an emergent physical time law. The base torus and fiber remain supplied.

## Detection and complete accounting

Use the same triangle-derived $S_3$ group and a compact two-link noncommuting-reflection
seed. After each changing update, compare the tagged placement with its **most recent**
visit. If the full spanning-forest gauge signature differs, retain the interval and stop.
Gauge-equivalent returns replace the saved visit. This policy is explicit: it does not
search all possible pairs of observation times, and a budget miss does not prove absence.

For sides $3,4,6$ and scheduler seeds $0,1$, use a 10,000-attempt budget per run. The
first detected intervals are:

| Torus side | Scheduler seed | Attempt interval | Changing moves in interval | Commuting control gauge return |
| ---: | ---: | ---: | ---: | :--- |
| 3 | 0 | 401–586 | 68 | yes |
| 3 | 1 | 156–1085 | 332 | no |
| 4 | 0 | 1552–4457 | 636 | no |
| 4 | 1 | 4585–7031 | 564 | yes |
| 6 | 0 | none in 10,000 attempts | — | not evaluated |
| 6 | 1 | none in 10,000 attempts | — | not evaluated |

Every detected noncommuting interval changes its gauge state by definition; the
commuting comparisons are **not** part of the selection criterion. Both bounded negatives
are exported, including their entire schedules and C++ replay, rather than suppressed.
Attempt ticks count no-ops as well as changing moves. Two seeds are a reproducible pilot,
not an estimate of event rates or a finite-size scaling law.

The C++ engine independently replays each complete history, each selected prefix, and
each return interval. Python independently verifies boundary links, spectator holonomies,
the full gauge quotient, and exact reversal. A separate occupancy-only replay reconstructs
every tag trajectory and agrees with the changing-event clock. Random local-frame tests
find the same interval and appropriately transformed raw endpoints.

For each selected interval, replay the **entire same history** from a seed whose two link
reflections commute. Do not reset its state at the beginning of the interval. The vacancy
rule gives precisely the same occupied-face trajectories and changing-update times in
both conditions. This controls the routing and the observation clock.

## Periodic topology is measured, not ignored

All four detected intervals cross periodic seams and have nonzero individual winding
vectors. Consequently these are not isolated planar braids, and we do not identify their
action with the controlled four-state $A_4$ component.

Lift each dual-face step to its unambiguous short displacement in the covering mesh.
Integer coordinates are three times each face barycenter. When every tag returns, its
net displacement is $3n(w_x,w_y)$ for side $n$. Both displacements and integer winding
vectors are exported. Zero abelian winding alone would still not imply contractibility
in a punctured torus.

There is an exact prediction for the commuting control. All its links lie in
$\{1,h\}\cong\mathbb Z_2$, with $h^2=1$. Each changing vacancy move toggles its one shared
link by $h$. The endpoint difference

$$D_e=U_e^{\mathrm{after}}(U_e^{\mathrm{before}})^{-1}$$

is therefore the mod-two crossing cochain of the dual defect worldlines. Returning
all tags restores the raw face holonomies in this subgroup, so $D$ is flat on every
triangle. Its two primal periods are the intersection parities

$$p_x=\sum_{k=1}^4 w_{y,k}\pmod2,\qquad
p_y=\sum_{k=1}^4 w_{x,k}\pmod2.$$

A flat $\mathbb Z_2$ connection on this torus is a vertex-frame change exactly when both
periods vanish. To see sufficiency, integrate the flat cochain from a root on the simply
connected triangular covering mesh. Flatness makes integration independent of contractible
path changes, and vanishing periods make the resulting vertex frames periodic. Necessity
follows because vertex-frame factors cancel around every loop. Equivalently, the only
remaining obstruction after fixing a tree is the two torus periods.

The implementation constructs $D$ from the **actual endpoint links**, verifies its
flatness, measures periods along explicit horizontal and vertical primal loops, and checks
the predicted worldline parities against an independent full-connection frame witness.
The four observed period pairs are $(0,0),(1,0),(1,1),(0,0)$ in table order, explaining
both positive and negative commuting-control returns.

In the first and fourth intervals the noncommuting connection retains memory despite
the commuting control's zero periods. This distinguishes the response from this abelian
control, but does not eliminate nonabelian torus-handle effects or prove a purely local
braiding interpretation.

## Graph consequences and the remaining dynamical gap

Each selected noncommuting interval changes the exact characteristic polynomial of its
full unweighted bundle-graph Laplacian. The first distinct trace powers are $3,3,4,5$;
the signed differences (after minus before) are respectively $-54,-12,104,-120$.
The two gauge-returning commuting controls have identical exact spectra. The others
change, consistent with their nontrivial torus periods. These very small tori have short
noncontractible graph cycles; the low trace powers are not evidence of a scale-independent
local interaction or an atomic wave law.

This rule still has an important limitation: **the hidden gauge state cannot influence
where the defects move**. Its occupancy update is exactly

$$ (\eta_u,\eta_v)\mapsto(\eta_v,\eta_u),\qquad \eta_f=\mathbf1[H_f\ne1], $$

with equal-occupancy pairs unchanged. Thus a state-independent patch scheduler yields a
closed classical exclusion process for positions, regardless of the holonomy history.
The control's identical trajectories are not an accident; they follow from this factor.
Finding unrouted memory closes the prescribed-path gap, **not the back-reaction or binding
gap**. A next constructive test needs derived, gauge-covariant relative-holonomy feedback
in the primitive updates and must check whether its effects survive beyond supplied
finite periodic fixtures. An imposed potential or an assumed quantum Hamiltonian would
not establish that.

## Reproduce

```bash
cmake --build build -j 4
export OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
uv run --with numpy --with scipy --with python-flint python tools/unrouted_transport.py \
  --output out/unrouted-transport.json
diff -u data/unrouted-transport.json out/unrouted-transport.json
uv run --with numpy --with scipy --with python-flint python tests/unrouted_transport.py
```

The detector refuses budgets above 100,000 attempts rather than silently truncating.
Flat-vacuum and zero-attempt controls return no witness; a missing witness always records
the actual budget and consumed attempts. The saved data contain exact values only.
