# Releasing the boundary reverses the gauge-dependent delay

The compact-region threefold quiet-state delay is **not** an intrinsic
lifetime of a localized object. Releasing its boundary into the full
periodic primitive dynamics both shortens the delay and, in an aligned
dense background, reverses the ordering of the preparations.

Gauge alignment remains dynamically important: with the same complete
initial charge field, the nine preparations have a fourfold range of
initial global rotation-production rates. What changed was the available
cross-boundary interaction, not the microscopic rules.

## Preparations with shared exteriors, then no held links

Use the nine compact quiet configurations with charge field
$(0,1,1,1,1,2,1)$ from [the boundary-held calculation](triangle-strip-memory.md).
Preparation 0 belongs to the two-state quiet component; the other eight
are isolated quiet states under the contained rules. Their held-boundary
mean first exits were 22.5 and 7.5 attempts per face, respectively.

Embed all nine in an eight-by-eight periodic triangular mesh, with 128
faces and 5,760 uniformly proposed operator slots. Within each comparison,
every exterior link and **every initial face charge on the whole mesh**
is identical across the nine preparations. The evolution applies all
operators, including updates wholly outside the old region. No patch
boundary, local product, or exterior reference is held during evolution.

The observable remains first exit from the old region's quiet set: exactly
one vacancy and one rotation there, with every contained elastic and
reaction primitive inactive. Exiting that set does not mean the rotation
annihilated or the gauge information was erased.

Three explicit exterior preparations were used:

- **Vacuum exterior:** identity links outside the inserted region. Its
  boundary reflection requires a neighboring reflection; total mesh charge
  is 8. This is a sparse preparation, not a dense reference ensemble.
- **Aligned reflection exterior:** start from the constant reflection
  connection, and match its boundary-tree frame to the patch using an
  actual vertex gauge transformation. Every exterior face keeps charge 1;
  the inserted vacancy/rotation pair gives total charge 128.
- **Haar exterior:** independent uniform $S_3$ links outside the fixed
  initial patch links. The total charge varies with the exterior draw but
  is identical across its nine paired preparations. This is not sampling
  a fixed-total-charge equilibrium conditional distribution.

Boundary-frame matching matters in the second preparation. Overwriting
its outer links without first transforming the exterior creates a ring
of vacuum faces and changes the experiment. That unmatched preparation
was not used for the dense-background result below. The corrected gauge
map is applied to all edges incident to the transformed vertices; it
preserves exterior face charges and produces the required boundary links.

## First-exit results

| Exterior | Independent trials | Former two-state preparation | Mean of eight former isolated preparations | Paired difference |
| --- | --- | --- | --- | --- |
| Vacuum | 512 | $1.3634\pm0.0600$ | $1.2743\pm0.0563$ | $0.0891\pm0.0220$ |
| Aligned reflections | 2048 | $0.2318\pm0.0050$ | $0.5070\pm0.0094$ | $-0.2752\pm0.0079$ |
| Haar links | 512 | $0.9109\pm0.0446$ | $0.8709\pm0.0358$ | $0.0399\pm0.0315$ |

Times are attempts per face; uncertainties are standard errors across
independent trials. All nine trajectories in a trial share the same
scheduler sequence and exterior. The eight alternatives are averaged
within a trial before computing its paired difference; they are not
counted as eight independent samples. No trajectory reached the censoring
cap of 32 attempts per face. Haar total charges ranged from 130 to 172
in this batch.

The dense-background reversal is large compared with its sampling error.
It demonstrates that the previously slow preparation is not protected
against coupling to this exterior. Neither these first-exit means nor
their reversal establishes loss of all later gauge memory.

## Exact initial reaction growth, with the same charge field

The aligned background is inactive before insertion: vacancy rules have
no identity faces to swap, equal-reflection elastic rules do nothing,
and equal-reflection triples are nonreactive. The local nonabelian seed
introduces gauge-dependent active fans. It is then left to evolve without
further intervention.

Let $N_2$ be the total number of rotation faces and let $B(X),D(X)$ count
the primitive slots that create or remove one rotation. With $F$ total
faces and the unchanged $45F$-slot scheduler,

$$\boxed{F\,\mathbb E[N_2(t+1)-N_2(t)\mid X_t=X]
=\frac{B(X)-D(X)}{45}.}$$

This is an exact one-attempt drift expressed in attempts-per-face units,
not an assumed continuum rate equation. Initially all nine states have
$N_2=N_0=1$ and no annihilating fan. Exhausting the full microscopic
operator bank gives:

| Initial growth class | Number among the nine preparations | Active birth slots | Growth in attempts-per-face units |
| --- | --- | --- | --- |
| Strong | 4 | 192 | $64/15$ |
| Intermediate | 4 | 120 | $8/3$ |
| Weak | 1 | 48 | $16/15$ |

The active all-reflection fan counts are 16, 10, and 4, respectively;
each contributes twelve birth slots. These rates differ by a factor of
four despite identical complete initial charge fields. Every birth also
creates a vacancy, preserving the total charge exactly.

Preparation 0, the formerly slow held-boundary state, is in the strong
growth class. The weak class is one of the formerly isolated states.
The measured ordering reversal therefore has a microscopic explanation:
alignment with the evolving exterior opens different sets of reaction
fans. A local charge-only trapping time omits these gauge-dependent
cross-boundary channels.

## Autonomous evolution: the seed activates its surroundings

Initial production does not establish a self-sustaining front or a bound
object. Following the strong preparation 0 and weak preparation 3 under
32 independent paired schedules on each of two periodic meshes gives:

| Faces | Time, attempts per face | Strong-class mean rotation count | Weak-class mean rotation count |
| --- | --- | --- | --- |
| 128 | 0 | 1 | 1 |
| 128 | 16 | $24.38\pm0.88$ | $18.19\pm1.13$ |
| 128 | 32 | $30.28\pm0.54$ | $29.81\pm0.51$ |
| 128 | 128 | $29.91\pm0.53$ | $31.88\pm0.51$ |
| 512 | 0 | 1 | 1 |
| 512 | 16 | $26.44\pm1.45$ | $17.84\pm1.09$ |
| 512 | 32 | $81.53\pm3.14$ | $64.13\pm2.78$ |
| 512 | 128 | $123.03\pm1.11$ | $124.22\pm1.21$ |

The same two local raw-link preparations are translated into the larger
mesh; every other link is the aligned reflection background. No subsequent
forcing, reinsertion, boundary control, or change of rules occurs. The
vacancy count equals the rotation count throughout because $Q=F$.
Uncertainties are standard errors over schedules; the two classes are
paired within each size, but microscopic prefixes are not coupled between
sizes in this calculation.

On 512 faces, the paired rotation-count difference is $17.41\pm3.88$ at
time 32 and $-1.19\pm1.89$ at time 128. Thus this measured population
response retains a sizable alignment-dependent transient, but does not
show a resolved late-time difference in that larger-mesh batch. This is
not a statement that all gauge observables have forgotten the preparation.

The independently counted uniform-sector reference means are 30.9288
rotations for 128 faces and 124.1032 for 512. Late-time populations are
near these values. Agreement of a population mean is not proof of
convergence to the full stationary distribution or irreducibility.

The nonreflection support also expands. At time 16 on the 512-face mesh,
its mean maximum dual distance from the initial seven-face region is
$11.28\pm0.42$ in the strong class and $10.19\pm0.45$ in the weak class.
These are distances in the supplied face-adjacency graph. Later measurements
are influenced by the periodic box; they are not used to fit a front speed,
propagation exponent, dimension, or continuum law.

[Earlier calculations already established seeded activity spreading](triangle-spreading.md)
under the reaction/vacancy bank. Spreading itself is not a new discovery
here. This experiment connects the current elastic/reaction bank's
**boundary-held trapping calculation** to its released evolution and
quantifies the contrasting transients of specific, charge-identical gauge
alignments. It does not provide a particle or a molecule.

## Consequence for the search

The held quiet-set lifetime cannot be used as a boundary-independent
particle diagnostic. For these preparations, releasing the surroundings
produces gauge-dependent activation and spreading, not evidence of a
self-bound configuration. A stronger candidate would have to retain
localized relational structure under exchange with an evolving environment;
an initially slow closed-patch mode or an early population-growth difference
does not establish that.

Calculation: [triangle_released_traps.py](../tools/triangle_released_traps.py).
Autonomous continuation:
[triangle_seeded_activity.py](../tools/triangle_seeded_activity.py).
