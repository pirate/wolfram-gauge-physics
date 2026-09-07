# Noncommuting conversion channels and their exact statistical reduction

The next construction moves beyond the commuting-subgroup collision rule. An exhaustive
positivity/activation audit of the same 945 minimal triple rules finds **144** that convert
face-class populations from two noncentral input holonomies generating all of $D_4$.
They fall into three families of 48 with compatible derived charges. All three families
are implemented and tested on shared links, with one-step transport and matched controls.

There is also a crucial observer limit: under uniform sampling of each complete family,
the face-class transition probabilities close exactly. The resulting observable process is
a classical spatial conversion-and-vacancy-transport model. Noncommuting microscopic links
do not, by themselves, establish quantum matter or hidden gauge feedback in that observer.

## 1. Positivity without choosing an energy function

Let the four nonidentity conjugacy classes be the two reflection classes $[r],[s]$, the
quarter-turn class $[\rho]$, and the central half-turn $[z]$. For each minimal symmetry
closure, derive its population-change row $d$ and solve

$$d\cdot q=0,\qquad q(1)=0.$$

Here conjugation and orientation reversal preserve class populations, so all exchanges in a
minimal closure have the same row up to sign. The normalized charge space therefore has
dimension three or four. A strictly positive weight on every nonflat class exists precisely
when $d=0$ or $d$ has both positive and negative entries. A nonzero one-signed row is an
impossibility certificate, not a failure to find small integer weights.

For a mixed-sign row, set $P=\sum_{d_i>0}d_i$ and $N=-\sum_{d_i<0}d_i$. An explicit
positive solution is $q_i=N$ where $d_i>0$, $q_i=P$ where $d_i<0$, and $q_i=1$ otherwise.
This is a feasibility certificate, not a preferred mass assignment.

Of the 945 rules, 755 admit such a weight: 377 preserve all class populations, and 378
change them. Among those 378, exactly 144 have a class-changing channel whose input contains
two nonflat, noncentral holonomies generating the full group. All 144 move sixteen raw
triples, so none are discarded by the least-support tie. The census retains every rule,
its exact charge equations, positivity certificate or obstruction, and activation channels.

## 2. Three compatible charge families

The 144 selected rules partition by their derived charge equation:

| Family | Equation | Rules |
| --- | --- | ---: |
| 0 | $q_r=q_s+q_\rho$ | 48 |
| 1 | $q_s=q_r+q_\rho$ | 48 |
| 2 | $q_\rho=q_r+q_s$ | 48 |

In each row, $q_z$ remains independent. Mixing rules from different rows forces at least
one noncentral class weight to zero: subtracting the first two equations, for example,
forces $q_\rho=0$. Thus these are three mutually incompatible families **within the selected
144-rule set** if a strictly positive one-face weight is required. This is not a classification
of every compatible composite law, larger-support energy, or rule in the full 945-rule census.

Write the class with the summed weight as $h$, and the other two as $a,b$. A positive basis
of conserved quantities is

$$Q_a=N_a+N_h,\qquad Q_b=N_b+N_h,\qquad Q_z=N_z.$$

Choosing unit weights on $a,b,z$ gives a convenient positive bound

$$W=N_a+N_b+2N_h+N_z,\qquad N_{\rm nonflat}\le W.$$

The labels “heavy” and “light” refer only to this normalization. There is no calibrated
physical mass, energy, temperature, or coupling. Schematically the conversion channels are

$$h+a+1\leftrightarrow a+a+b,\qquad h+b+1\leftrightarrow a+b+b,$$

with specific ordered placements and multiplicities in the exported tables. They preserve
the ordered boundary product and the exact frame stabilizer. Both input and output generate
the full nonabelian group. This is not the forbidden isolated $r+r'\to z$ fusion: noncentral,
noncommuting information remains on both sides.

## 3. Mixed-arity local primitives, not a table stitched across overlaps

A *bank* is the set of all 48 reversible triple operators in one family, not their arbitrary
union as a function or an assertion that overlapping operators commute. Each selected event
uses the proven boundary-fixed two-spoke lift. The engine now registers multiple exact tables,
keeps spatial support separate from rule identity, accepts explicit per-event rule choices,
and records the chosen rule alongside causal parents. Inversion reverses both support and
rule history. Arity mismatches fail before changing state or provenance.

Transport uses the smaller, two-face primitive

$$(1,g)\leftrightarrow(g,1),\qquad g\ne1,$$

and its one-shared-edge lift. It preserves every additive class charge and moves a vacancy
and defect across one dual edge. Unlike the earlier three-face transport control, it does
not restrict motion to two dual hops. This is a separately declared microscopic operator,
not inferred as a unique transport law.

Solving the full sublattice-dependent conservation equations gives four dimensions for each
48-rule triple bank: two global noncentral charges and separately conserved central counts
on the two parts. Adding shared-edge transport leaves exactly three dimensions, all ordinary
global charges. This removes those *additive* sublattice restrictions; it does not prove
ergodicity or absence of other kinetic sectors.

## 4. Exact closure of the uniformly sampled face-class observer

For a raw triple $x$, count how many of the 48 operators produce each output class tuple.
For all $8^3$ raw inputs in each family, the resulting count distribution depends only on
the input class tuple. These give 125 exact coarse rows per family, each totaling 48.

This is the usual strong-lumpability condition: every microscopic representative of a coarse
state has the same total probability of entering each coarse state. Consequently independent
uniform sampling of the rule bank yields a Markovian class process on an actual shared-link
mesh, not just independent face variables. Local-frame conjugations preserve classes, each
lift changes only its input faces, and the pair transport is already class-determined.
The proof does **not** replace any microscopic replay with sampled class states.

This does not say an individual rule is a deterministic class map. It is a result about
averaging the complete family with the specified weights. Nonuniform rule probabilities,
fixed rule schedules, or more informative loop observables require their own closure checks.
For background on Markov coarse graining, see [Geiger and Temmel](https://arxiv.org/abs/1212.4375).

## 5. A stationary occupation relation derived from internal degeneracies

Let $n(c)$ be the size of conjugacy class $c$: it is two for $h,a,b$ and one for $1,z$.
Because each microscopic operator is an involution, the coarse transition counts obey

$$D(c)K(c,d)=D(d)K(d,c),\qquad D(c_1,c_2,c_3)=\prod_i n(c_i).$$

The tool verifies every directed entry of all three coarse kernels. Thus the formal spatial
class process has stationary weight

$$\mu(\mathbf c)\propto\prod_f n(c_f),$$

which can be restricted to any closed reachable component or conditioned on conserved
charges. This assertion does not count the full gauge-orbit fibers, identify all realizable
coarse configurations, or prove mixing of a particular seed's component.

An explicit stationary product family is obtained with positive bookkeeping fugacities
$x,y,t$ and one-face weights

$$w_1=1,\qquad w_a=2x,\qquad w_b=2y,\qquad w_h=2xy,\qquad w_z=t.$$

Normalizing them gives the exact grand-canonical occupation relation

$$\frac{\rho_a\rho_b}{\rho_h\rho_1}=2.$$

The factor two follows from the internal class sizes: $(2\times2)/(2\times1)$. It was
not inserted as a reaction constant. Tests substitute positive fugacities into every coarse
transition and verify detailed balance directly. Fixed-charge finite runs have correlated
occupations, so their marginal densities must **not** be substituted into this product-measure
formula without a justified factorization limit. We do not claim those runs have equilibrated.

Product-form stationary measures in reaction models are established mathematics; for example,
[Anderson, Craciun, and Kurtz](https://arxiv.org/abs/0803.3042) study a different, mass-action
setting. Our finite-exclusion spatial process is not their well-mixed Poisson construction;
the stationarity claim here follows directly from the checked degeneracy balance above.

## 6. Controlled raw-link evolution and retained null results

Each attempted event samples uniformly among the 49 registered operators and their rooted
placements: 48 conversion tables on fans, plus one vacancy-transport table on shared edges.
The transport-only control receives the same schedule but skips conversion attempts.
This specifies an experimental clock and relative attempt weights; neither is physical time
or a derived coupling. Idle slots remain in the clock.

The two datasets cover three families, seven initial conditions, two rule-set variants,
and two or four independently seeded schedules: **252 runs, 37.8 million attempted clock
ticks, and 17,571 changing events**, including **311 noncommuting conversions**. C++ replays
the full reverse schedule and intermediate histogram echoes. Independent Python link algebra
checks every attempted enabled operator, every changing target, every event's derived charges,
sampled histograms, final links, and exact reversal. No-op targets are skipped only after
checking the current raw holonomies; their link-identity property follows from the unique lifts.

- On the 72-face mesh, the mixed pairs containing $h$ react in all twelve such combined runs,
  with 3–52 conversions. The six two-light-pair runs have no conversion in those histories.
- On the 288-face mesh, twelve of 24 heavy/light paired runs react, with 1–19 conversions;
  the other twelve, and all twelve two-light-pair runs, are null. The smaller CI mesh does
  contain two-light-pair reactions, so those larger-mesh nulls are not impossibility claims.
- Flat and single-species seeds have no conversions. Subgroup closure supplies a structural
  control: the initial cyclic subgroup remains closed and cannot supply an active full-group
  tuple. All no-conversion histories end at exactly the same raw links as their matched controls.
- The positive weight bounds support throughout. No persistent localized bound state,
  binding energy, force law, continuum limit, or quantum interference is demonstrated.

The next physics-facing test is the derived stationary/reaction-diffusion benchmark with
finite-size, component, and mixing controls. For gauge-dependent observables beyond this
closed class process, the uniform averaging assumption is itself an information ceiling.
Changing the quantum or geometrical state construction remains separate work, not something
that more classical sampling can establish.

## Reproduce

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build -L wgphysics --output-on-failure
python3 tools/triple_channels.py --output out/triple-channels.json
diff -u data/d4-triple-channels.json out/triple-channels.json
python3 tests/channel_banks.py
python3 tools/run_channel_banks.py --output out/channel-banks.json
diff -u data/d4-channel-banks.json out/channel-banks.json
python3 tools/run_channel_banks.py --side 12 --attempts 200000 --stride 2000 --trials 4 --seed 7193901 --output out/channel-banks-replication.json
diff -u data/d4-channel-banks-replication.json out/channel-banks-replication.json
```
