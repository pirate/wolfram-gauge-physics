# Stored braid memory changes actual reaction rates

The triangle-fiber model now connects its previously measured holonomy memory to
**local reaction feedback on raw links**. Two states with identical face-curvature
classes everywhere and the same encounter-boundary holonomy have different outgoing
class-transition rates, even after averaging over the entire uniform update bank and
every rooted mesh support. One has twelve enabled class-changing reaction operators;
the other has none. Both C++ and Python enumerate all 11,232 choices independently.

The construction starts with the triangle adjacency and its derived automorphisms,
$\operatorname{Aut}(C_3)=S_3$. An exhaustive minimal-rule census derives reversible,
gauge-covariant three-face updates. A positive additive charge emerges from their
conservation equations; no potential, Hamiltonian, quantum amplitude, or binding target
is supplied. The new rules remove the vacancy model's exact absence of positional
feedback, but do **not** establish stable matter, physical energy, quantum statistics,
emergent geometry, or schedule-independent spacetime dynamics.

![Actual fan reaction and full-generator readout of braid memory](images/triangle-feedback.png)

The diagram uses the supplied mesh coordinates. The orange spokes are actual changed
links. The right panel is an exact generator count, not a fitted reaction probability.

## 1. Derive the rule bank, then derive its charge

Enumerate all minimal transposition closures on $S_3^3$ satisfying:

- preservation of the ordered boundary product $abc$;
- covariance under simultaneous conjugation;
- covariance under orientation reversal $(a,b,c)\mapsto(c^{-1},b^{-1},a^{-1})$;
- involutivity and a fixed flat vacuum $(1,1,1)$.

“Minimal” means the conjugation/reversal closure of one exchanged pair of distinct
triples, provided the resulting transpositions do not overlap inconsistently. This is
not a census of all composite rules and does not impose a triple braid relation.

The C++ breadth-first closure and independent Python finite-action orbit construction
agree on **144** minimal rules. Select **all** of those that change the occupied-face
count and admit a strictly positive additive class charge on nonidentity classes. This
gives twelve rules, each moving 24 of the 216 raw triples. Feedback behavior and observed
trajectories do not enter the selection criterion; no favored member is chosen afterward.

Write $R$ for the reflection class and $C$ for the three-cycle class. The full bank's
additive conservation equations have the one-dimensional solution

$$q(1)=0,\qquad q(R)=1,\qquad q(C)=2,$$
$$Q=\sum_f q(H_f)=N_R+2N_C.$$

The population equation is $2q(R)-q(C)=0$. Positivity is checked exactly, not by a bounded
search over guessed masses. Independently, breadth-first word length using the three
triangle-edge transpositions gives the same values $0,1,2$ on the derived automorphisms.
That agreement is an algebraic interpretation of the derived charge, not an inserted
energy law or a derivation of a physical mass ratio.

Combine these twelve three-face primitives with the existing one-edge vacancy transport.
Each attempted tick chooses one of the thirteen rules and one rooted support uniformly.
The scheduler does not inspect the charge or holonomy state.

## 2. The complete local rate law retains hidden holonomy information

All class conversions have the population form

$$3R\ \rightleftarrows\ R+C+1,$$

so a forward event has $(\Delta N_R,\Delta N_C)=(-2,+1)$ and a reverse event $(+2,-1)$.
The surviving reflection participates in the rule: it is not assumed to be an unchanged
spectator or a physical catalyst.

For three based reflections $(a,b,c)$, **all twelve** rules act precisely when

$$[a=b]\mathbin{\oplus}[b=c]=1.$$

Thus exactly one neighboring pair must be equal. Three equal reflections, equal endpoints
with a different middle reflection, and three distinct reflections are inert. An active
triple has six possible outgoing class placements, each realized by two rules. A tuple
whose charge values are a permutation of $(0,1,2)$ has exactly four of the twelve rules
enabled in the reverse direction. Every one of the 216 input triples is checked.

These equality tests are gauge invariant because the three holonomies share an explicit
basepoint. They are not comparisons of labels in unrelated local frames.

For an explicit same-boundary witness, the exported group elements satisfy

$$ (1,1,5)\quad\text{and}\quad(1,2,1). $$

Here the numbers are **group enumeration indices**, not charges or the identity symbol;
index zero is identity. Both triples have three reflection-class entries and ordered
product equal to group element five, but they are not simultaneous-conjugation copies.
The first has twelve outgoing conversions; the second has none. This nonclosure survives
uniform averaging over the whole compatible bank, unlike the earlier square-fiber bank's
[closed class-rate factor](reachable-gauge-equilibrium.md). The comparison is about these
specified banks, not a general impossibility theorem for square fibers.

## 3. Read out the previously generated braid memory on a real mesh

Take the four complete gauge states reached by the
[two closed triangle-fiber circuits](noncommuting-fiber-transport.md). Apply the **same
25 vacancy moves** to each state, bringing three defects into a common three-face fan
while leaving the fourth outside. The route depends only on occupancy. It contains no
reaction rule and no edit of the internal holonomy labels.

The resulting based triples, in the existing exported group enumeration, are

$$ (2,5,5),\quad(5,5,1),\quad(5,2,5),\quad(1,5,2). $$

All four connections have the same complete face-class field and remain mutually
gauge-inequivalent. States one and two also have exactly the same based encounter product.
Their class-changing reaction counts differ: twelve versus zero.

To avoid confusing a selected local rule with the full dynamics, enumerate all thirteen
rules on all 864 rooted supports of the $12\times12$ mesh. C++ applies each operator and
its exact inverse before the next test; Python independently computes the resulting
full face-class fields. The complete outgoing class-rate distributions agree between
implementations. The four class-changing **triple** counts are $(12,12,0,0)$.
Consequently the next attempted update produces a reaction with probability $1/936$ in
the first two states and zero in the other two under this particular uniform scheduler.
Vacancy moves still occur; their class-transition rates agree across the four states.

This is a concrete bridge from a stored classical holonomy relation to a change in
local curvature populations. It reads a partition of the four states, not all four
states in a single measurement. It is neither quantum measurement nor evidence that
the simulation autonomously manufactures the preceding controlled braid histories.

## 4. Unrouted evolution and negative controls

On a $6\times6$ torus, start from a flat connection, a commuting two-reflection-link seed,
or a noncommuting two-reflection-link seed. The latter two each have $Q=4$. Run four
fixed scheduler seeds for 100,000 attempted ticks each, comparing the full bank with
transport-only evolution on the **same clock**; triple-rule attempts are no-ops in the
transport control.

- Scheduler seed zero has five forward and five reverse conversions. Its first forward
  event is at tick 294.
- Seeds one and two have no conversions within the budget.
- Seed three has one forward and one reverse conversion.
- Every flat and commuting control has zero conversions, as do all transport-only runs.

All conditions, including bounded negatives, are exported. Four scheduler seeds at one
size and charge are a pilot, not a reaction-rate estimate, scaling law, or evidence of
typical long-time behavior. There is also an exact reason for the commuting controls:
every triple in a single reflection subgroup $\{1,h\}$ is fixed by every reaction rule,
and vacancy transport cannot leave that subgroup.

Independent C++ and Python replay verifies each event, every sampled histogram, all
final links, every boundary link and spectator holonomy, and conservation of $Q$ after
every changing update. Reversing the actual primitive history recovers the exact initial
connection, with matching intermediate histogram echoes. Local-frame-transformed
readout states evolve covariantly. These tests establish this classical model's link
dynamics, not a physical continuum interpretation.

## 5. The same charge bounds full-graph spectral capacity

The derived charge also has a graph-native spectral meaning. For the flat triangular
base, the horizontal Laplacian band is $[0,9]$. The triangle fiber has Laplacian
$L(C_3)=3I-J_3$ with eigenvalues $0,3,3$, so the **full** infinite flat bundle has band
$[0,12]$. This differs from selecting a two-component sector of the square fiber:
twelve is the upper edge of the whole flat triangle-fiber graph.

For an oriented base face, let $S_f$ be the $9\times9$ matrix consisting of identity
diagonal blocks and transported adjacency on its three edges, using the full three-vertex
permutation representation. On the supplied triangular tori, direct incidence counting gives

$$2(12I-L_{\mathrm{bundle}})=\sum_f\iota_f S_f\iota_f^T
 +2\bigoplus_v J_3.$$

Here $\iota_f$ inserts a face's three fibers into the full vertex space and
$J_3=\mathbf1\mathbf1^T\ge0$. Gauge fixing two edges of a face leaves only its holonomy
$H_f$. For a holonomy eigenphase $\theta$, the corresponding three eigenvalues of $S_f$
are $1+2\cos((\theta+2\pi j)/3)$, $j=0,1,2$. An identity eigenphase contributes no
negative eigenvalue, and each nonidentity eigenphase contributes exactly one. Thus

$$n_-(S_f)=\operatorname{rank}(I-P_{H_f})=q(H_f)\in\{0,1,2\}.$$

Exact integer congruence elimination checks this equality on all $6^3=216$ raw face
connections. A separately assembled full 27-vertex bundle matrix checks the incidence
identity entry by entry, with no floating-point decisions.

Each embedded face form is nonnegative on a subspace of codimension $q(H_f)$. The
intersection of these subspaces has codimension at most $Q$, and the $J_3$ term is
nonnegative everywhere. Consequently

$$\boxed{\quad n_+(L_{\mathrm{bundle}}-12I)\le Q=N_R+2N_C.\quad}$$

The primitive dynamics therefore conserves an **upper bound** on the number of modes
above the full flat reference band. It does not conserve the actual mode count, guarantee
such modes exist, establish their persistence, or assign them physical frequencies.
This connects a charge derived from the reaction rules to a spectral constraint on the
actual bundle graph without declaring the graph Laplacian to be a Hamiltonian.

## Remaining physical gap

Reversibility, positive conserved charge, and hidden-state feedback do not imply
attraction, long-lived localization, or binding. Uniform scheduling of reversible
microscopic permutations still admits a uniform stationary measure on a finite closed
component. Whether useful persistent structures occur, and whether any observed
correlation survives size, density, scheduler, and topology controls, remains to be
tested. The earlier obstruction to treating graph-localized modes as conserved particles
also remains. Do not rename $Q$ “energy” or prescribe a Schrödinger evolution to bypass
these missing derivations.

## Reproduce

```bash
cmake --build build -j 4
export OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
uv run --with numpy --with scipy --with python-flint python tools/triangle_feedback.py \
  --output out/triangle-feedback.json
diff -u data/triangle-feedback.json out/triangle-feedback.json
uv run --with numpy --with scipy --with python-flint python tests/triangle_feedback.py
uv run --with numpy --with scipy --with matplotlib python tools/plot_triangle_feedback.py
```

The census and mixed-bank runner accept `--cycle 3`; omitting the argument retains the
historical square-fiber protocol. The runner also supports cycle fibers with five and
six vertices, tested with independently replayed transport. The exhaustive census keeps
its separate bound of at most eight group elements and refuses larger groups rather than
claiming completeness. The exported dataset contains exact values only.
