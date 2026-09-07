# Reachable gauge equilibrium and a finite-sector binding restriction

The [noncommuting channel banks](nonabelian-channel-banks.md) now have an exact,
link-realizable finite-system equilibrium benchmark, not just a formal product measure.
On the 18-face torus, complete search reaches all **468,180** class configurations at
charges $(Q_a,Q_b,Q_z)=(4,2,0)$ from the two-link $h+a$ seed. Independent raw-link
trajectories approach the predicted occupation proportions and conversion frequency.

The same result supplies a limitation relevant to the molecule objective: the equilibrium
class distribution is invariant under permutations of face positions. It cannot prefer
nearby over distant arrangements beyond the supplied graph's pair-distance multiplicities.
Changing positive, state-independent sampling weights of these same reversible operators
cannot change that equilibrium, although it can change kinetics and observer memory.

![Raw-link residence times and conversion frequencies compared with exact predictions](images/bank-equilibrium.png)

## 1. Prove reachability from links before assigning equilibrium weights

The full 48-rule family has the checked **strong** class-lumpability property: every raw
triple with a given ordered class tuple has the same number of rules leading to each
output class tuple. Thus every positive coarse transition lifts from **every** raw
connection exhibiting its input classes. It is not just a transition between a conveniently
chosen pair of representatives. The boundary-fixed link lifts alter no spectator faces.

`tools/bank_equilibrium.py` independently reconstructs those counts from all $512$ raw
triples and all 48 tables, checks the saved kernel, charge conservation, and detailed
balance, then explores every changing class transition. Two equivalent rooted vacancy
swaps are deduplicated only for reachability, never for sampling rates. Starting with an
actual link connection, induction along every positive search path proves realizability
of every discovered class state. Representative shortest paths to $N_h=1$ and $N_h=0$
take one and two changing events; Python and C++ lift and reverse them on raw links,
including a test with independent local frame changes.

At $F=18$ faces and charges $(4,2,0)$, write $k=N_h$. Then

$$N_a=4-k,\qquad N_b=2-k,\qquad N_1=12+k,\qquad N_z=0.$$

The number of formal configurations is

$$A_k=\frac{18!}{(12+k)!(4-k)!(2-k)!k!},\qquad k\in\{0,1,2\}.$$

| $N_h$ | Formal configurations | Discovered configurations |
| ---: | ---: | ---: |
| 0 | 278,460 | 278,460 |
| 1 | 171,360 | 171,360 |
| 2 | 18,360 | 18,360 |

The queue is exhausted, not stopped at a budget. Equality in each population bin proves
there are no additional inaccessible class configurations in this sector. This does not
enumerate microscopic connection components or establish connectivity at larger sizes.
Tests deliberately exhaust a different, single-species seed: only its 153 two-heavy-face
placements are reachable, although its formal charge sector contains additional states.
An interrupted search explicitly withholds the canonical comparison.

## 2. Canonical equilibrium, not products of correlated marginal densities

Each noncentral class has two raw group elements, whereas the identity and central class
have one. Detailed balance therefore gives class-configuration weight $2^{6-k}$ here.
The full reachable sector has population probabilities

$$\Pr(N_h=k)=\frac{A_k2^{6-k}}{\sum_{j=0}^2 A_j2^{6-j}}
  =\frac1{241}(182,56,3)_k,$$

and hence $\mathbb E[N_h]=62/241$. Every raw triple has at least 36 identity actions among
the 48 rules. The specified 49-operator clock therefore has a self-loop probability of
at least $36/49$ at every state. Connectedness and these self loops give a unique limiting
class distribution by the standard finite-chain convergence theorem; see
[Levin and Peres, Theorem 4.9](https://pages.uoregon.edu/dlevin/MARKOV/mcmt2e.pdf).
This proves asymptotic convergence, **not a numerical mixing-time bound**.

For comparison, products of the finite canonical marginal densities give

$$\frac{\mathbb E[N_a]\mathbb E[N_b]}{\mathbb E[N_h]\mathbb E[N_1]}
 =\frac{13530}{6541}\simeq2.06849,$$

not the grand-canonical product-measure value two. No fitted correction or independence
assumption is needed: the conserved-charge correlations are counted exactly.

## 3. A conversion-frequency prediction in the specified clock

Conditioned on $N_h=k$ at equilibrium, class positions are uniform. An ordered fan's
three-face marginal is consequently hypergeometric, not three independent draws. If
$m_c(t)$ is the multiplicity of class $c$ in tuple $t$, its probability is

$$\Pr(t\mid k)=\frac{\prod_c (N_c(k))_{m_c(t)}}{(18)_3},$$

where $(n)_j=n(n-1)\cdots(n-j+1)$. Multiply this by the number of conversion rules
giving a particular change and by $1/49$, the attempt weight of each table. Summing
over all ordered input tuples yields stationary conditional rates

$$R_{0\to1}=\frac8{2499},\quad R_{1\to0}=\frac{26}{2499},\quad
  R_{1\to2}=\frac1{1666},\quad R_{2\to1}=\frac4{357}.$$

The exact population fluxes balance in both directions. Averaging over the canonical
population distribution gives

$$\mathbb E[C_{\rm attempt}]=\frac{440}{86037}\simeq0.00511408.$$

These are **stationary conditional expectations**, not a claim that $N_h$ alone is a
Markov process: spatial arrangements affect instantaneous collision opportunities.
The factor $1/49$ is an explicitly supplied scheduling choice, not a physical coupling
or an emergent unit of time. A five-face whole-configuration enumeration independently
checks the hypergeometric implementation.

## 4. Raw-link measurements, controls, and uncertainty

The first checked dataset uses four independently seeded runs from each of three actual
link states with initial $N_h=0,1,2$. Each run has one million attempts; the first 200,000
are excluded from the reported stationary comparison. Matched transport-only controls
receive the same schedule and skip the conversion slots. Their populations stay fixed.

| Observable | Exact prediction | Measured mean $\pm$ one run-based standard error |
| --- | ---: | ---: |
| $\Pr(N_h=0)$ | 0.755187 | $0.753483\pm0.002623$ |
| $\Pr(N_h=1)$ | 0.232365 | $0.233819\pm0.002259$ |
| $\Pr(N_h=2)$ | 0.012448 | $0.012698\pm0.000768$ |
| Conversions per attempt | 0.00511408 | $0.00513844\pm0.00005541$ |

Occupations integrate **all post-attempt states**, including self loops, directly from
the independently checked changing-event logs. They are not event-conditioned averages
or occasional histogram samples. C++ reverses every enabled attempt and intermediate
histogram; independent Python link algebra checks each enabled attempt and every changed
target, final connection, charge, and reverse echo. Schedule seeds and checksums, event
logs, 50,000-attempt occupation blocks, initial-state strata, and final links are retained.

Uncertainty is estimated across independently seeded run means, not by pretending that
millions of correlated ticks are independent observations. Initial-state strata and
early/late post-burn comparisons are retained. They support consistency of these measured
observables, not full-state equilibration at the chosen burn or a continuum limit.

A separately seeded replication keeps the same initial states, run lengths, burn, and
measurement rules. Its twelve run means give occupation proportions
$(0.756339,0.230869,0.012791)$ and conversion frequency
$0.00508177\pm0.00005953$ per attempt, again consistent with the predictions.
Together the datasets contain 24 combined runs and 24 matched controls, 48 million
attempted clock ticks including skipped control slots, and 122,261 conversions.
The stationary measurements cover 19.2 million combined-run post-burn attempts and
98,114 conversions. Nothing was fitted to the replication outcomes.

## 5. Why this sector cannot supply equilibrium class binding

For distinct face positions $f,g$, exchangeability gives

$$\Pr(c_f=a,c_g=b)=\frac{\mathbb E[N_aN_b]}{18\cdot17}
  =\frac{812}{36873},$$

independent of their separation. Here $\mathbb E[N_aN_b]=1624/241$, which is not the
product of mean populations. The dual graph has 54, 108, 108, and 36 ordered distinct
face pairs at distances 1 through 4. Multiplying these counts by $812/36873$ gives the
exact expected class-pair distance histogram. This uses **spatial dual-graph distance**,
not causal or branchial adjacency. The tool also integrates that histogram from verified
raw-link events, checking its input classes at every event and its final classes against
the final raw connection.

The first dataset has a nearest-neighbor pair deficit of 2.63 estimated run-based standard
errors: $1.173942\pm0.005787$, versus the prediction $1.189163$. We retain that discrepancy.
The independent replication instead gives $1.200579\pm0.008379$, so the deficit does not
repeat. Pooling all 24 independently seeded runs gives

| Dual distance | Predicted mean $a$-$b$ pair count | Measured mean $\pm$ one standard error |
| ---: | ---: | ---: |
| 1 | 1.189163 | $1.187260\pm0.005702$ |
| 2 | 2.378326 | $2.376998\pm0.005530$ |
| 3 | 2.378326 | $2.378620\pm0.007639$ |
| 4 | 0.792775 | $0.793445\pm0.004499$ |

This supports the spatial prediction within the reported sampling uncertainty; it is not
a rigorous finite-time mixing certificate or a significance test for every possible observable.

This is an equilibrium restriction for these local class observables in the proven
sector. It does not exclude long-lived transient correlations, information in larger
loops, other kinetic components, or different microscopic constructions. It does mean
that a stationary molecular-looking density cluster in this observer would contradict
the present model and needs investigation, not celebration.

### Positive schedule reweighting cannot evade the restriction

There is a stronger elementary consequence of reversibility. Let $C$ be any connected
component of the **raw-link** graph under the same supported generators. Each local
update $T_i$ is an involutive permutation. For positive state-independent probabilities
$p_i$, the transition matrix

$$P=\sum_i p_i T_i$$

is symmetric and preserves the uniform distribution on $C$. Changing the positive
weights leaves both the supported raw graph and this stationary distribution unchanged.
For uniform weights, strong lumpability and the proven connected class sector imply
that the projection of that uniform raw distribution is exactly the canonical class
distribution above. Therefore its projection remains the same for any other positive,
state-independent weights on these same generator-placement pairs—even when the new
class observer is no longer Markovian. Self-loop support also survives positive reweighting.

This argument does not require counting $C$ or declaring different gauge frames physical.
It concerns raw states in one component and their gauge-invariant class projection.
Tests independently enumerate every raw **local triple** component: projected state
counts have the predicted class-degeneracy ratios, and weights $1,2,\ldots,48$ leave
uniform raw stationarity intact while breaking class-rate closure. The global conclusion
follows from the permutation argument, not extrapolation from that local enumeration.

Thus adjusting constant reaction frequencies can alter relaxation and memory, but cannot
create equilibrium attraction in this sector. To seek more physics, the open questions
are genuinely richer observables, differently constrained reachable components, and
justified geometry/fiber/state constructions—not a hand-written attractive potential.

## Scope and reproduction

This is a finite classical model with $D_4$ derived from the chosen square fiber, a supplied
2D torus, and an externally specified random update clock. It establishes an exact
microscopic-to-statistical benchmark and a scoped obstruction; it is not a new general
theorem of statistical mechanics, a field-priority claim, quantum matter, or a molecule.

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j 4
python3 tests/bank_equilibrium.py
python3 tools/bank_equilibrium.py --trials 4 --component data/d4-bank-equilibrium.json --output out/bank-equilibrium.json
diff -u data/d4-bank-equilibrium.json out/bank-equilibrium.json
python3 tools/bank_equilibrium.py --trials 4 --seed 10472917 --component data/d4-bank-equilibrium.json --output out/bank-equilibrium-replication.json
diff -u data/d4-bank-equilibrium-replication.json out/bank-equilibrium-replication.json
uv run --with matplotlib python tools/plot_bank_equilibrium.py
```

The component is re-enumerated on reproduction, not trusted from its saved Boolean flag.
The three family class kernels agree after a checked renaming of $h,a,b$; the measurements
reported here use family 0 only. No general large-mesh component or mixing claim is made.
