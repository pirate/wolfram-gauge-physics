# Exact observable memory, frozen sectors, and finite-size false invariants

The shared-link simulations are reversible, but coarse observables discard information.
This audit measures that loss exactly, separates update-schedule behavior from reachability,
and rejects neighboring-loop invariants that exist only on the smallest test geometry.

## Coarse randomness can be missing memory

Let $F$ be the selected involutive four-face update and $O(x)$ the ordered tuple of individual
face conjugacy classes. There are 4,096 rooted gauge-fixed configurations and 625 coarse labels.
Under the uniform measure on these rooted configurations, define

$$m_i=|O^{-1}(i)|,\qquad
M_{ij}=\frac{|\{x:O(x)=i,\ O(Fx)=j\}|}{m_i}.$$

The one-step coarse transition counts are symmetric, so $\pi_i=m_i/4096$ obeys exact detailed
balance. Nevertheless, propagating with $M^2$ is wrong: it implicitly resamples the hidden state
after the first step. The actual microscopic map satisfies $F^2=1$.

$$\Pr\bigl(O(F^2x)=O(x)\bigr)=1,\qquad
\sum_i\pi_i(M^2)_{ii}=\frac{467}{512}.$$

For 192 of the 625 coarse labels, $M^2$ invents spreading that does not occur; its smallest
return probability is $1/4$. The stationary-weighted total-variation error is $45/512$ and
the largest conditional error is $3/4$.

The augmented observation $(O(x),O(Fx))$ has 757 labels and is closed for this **single
involution**, since the next augmented observation is just the reversed pair. After one step,
this can be expressed through the previous/current observations. It does not show that one
step of memory suffices for multiple different mesh generators or arbitrary schedules.

The weighting is uniform over rooted holonomies, not uniform over physical gauge orbits.
This is classical hidden-state memory, not a derivation of quantum amplitudes or entanglement.

## Exact small-mesh dynamics and frozen states

For a fixed labeled tetrahedron, a spanning-tree gauge leaves three independent group elements:
$8^3=512$ rooted configurations. Simultaneous conjugation reduces them to **176 physical gauge
states**. We enumerate all 48 adjacent-triangle closing-link patches and all 8,448 transitions.
Independent Python path algebra checks the C++ output, including frame normalization, face
classes, stabilizers, and neighboring-pair features.

Allowing every generator gives 62 connected components:

| Component size | Number of components |
|---:|---:|
| 1 | 56 |
| 6 | 2 |
| 12 | 1 |
| 24 | 2 |
| 48 | 1 |

For this rule, the center/noncenter distinction gives an exact frozen-state criterion on a
mesh with connected face adjacency. The raw pair map is identity when both inputs are central
or both are noncentral, and exchanges those statuses when exactly one is central. Therefore:

- All-central or all-noncentral face configurations are fixed under every available patch.
- A mixed configuration has an adjacent mixed pair; updating it changes a gauge-invariant
  central status and therefore changes the physical state on the fixed labeled base.

Eight tetrahedron gauge states are all-central and 48 are all-noncentral. Flatness is thus
not the only source of frozen behavior. Reversibility also prevents a nonfixed state from
entering a state fixed by every generator.

Reachability is not the same as a deterministic trajectory. Applying all 48 generators in
their exported index order gives the identity permutation on all 176 states. A fixed alternative
ordering, sorting patch IDs $p$ by $17p\bmod48$, gives cycles of lengths 1, 2, 3, 4, and 6.
No claim of schedule-independent dynamics follows from pair-level braid identities.

Individual face labels happen to be a closed observer on the tetrahedron: refinement under
all generators remains at 149 classes. They still do not distinguish all 176 physical states.
The four-face open-patch counterexample shows that this closure is not a general mesh property.

## A neighboring-loop energy ansatz and its falsification

For each ordered adjacent-face patch, transport both loop holonomies to a common basepoint
and quotient by simultaneous conjugation. There are 28 pair classes. Let $n_i(x)$ count patches
in class $i$, and consider the same coefficient function at every patch:

$$Q_w(x)=\sum_{i=1}^{28} w_i n_i(x).$$

On all tetrahedron transitions, the exact rational equations
$w\cdot[n(F_px)-n(x)]=0$ have coefficient nullity 23. After removing zero functions and
state-independent constants, **12 independent nonconstant invariants** remain.

They are finite-size artifacts. On an octahedron, 128 seeded random link configurations with
16 elementary updates each provide 2,048 additional exact constraints. Every before/after
feature vector is checked by independent link algebra. Each of the twelve tetrahedron basis
functions has an explicit violation, and the combined nullspace has no nonconstant combination.

### Why this conclusion is stronger than "constant on our samples"

Reversing an ordered patch exchanges its commonly based holonomies up to simultaneous
conjugation. Let $r(i)$ be this class reversal. The complete patch family includes each reversal,
so **on every mesh using that family**, $n_i=n_{r(i)}$.

There are ten exchanged class pairs and eight self-reversing classes. Every vector in the
11-dimensional joint nullspace satisfies

$$\frac{w_i+w_{r(i)}}2=c\quad\text{independently of }i.$$

Its antisymmetric part cancels universally, leaving $Q_w=c\,N_{\mathrm{patches}}$. Thus all
surviving weights define constants on any such mesh, not merely on the sampled states. A
universal invariant must satisfy the finite checked constraints; their complete remaining
space is already universally trivial. This excludes nonconstant invariants in **this common
28-coefficient ansatz** on a mesh family including these tests.

It does not exclude longer loops, unequal positional/orientation coefficients, different
connector/path families, stroboscopic invariants, or nonlocal and topological quantities.
The [one-face normal-subgroup theorem](face-energy-obstruction.md) is a separate, more general
result for individual face densities across the full strict rule census.

## Consequence for the relaxation experiments

The uniform independent-link measure remains an exact stationary measure of the reversible
mesh dynamics. Agreement of seeded local histograms with that measure is not a proof of
ergodicity, thermalization, or an energy ensemble. The exact memory error, frozen sectors,
and schedule-sensitive cycles are discriminating controls that must accompany such plots.
The tetrahedron does not establish the large-mesh component structure; that remains unmeasured.

## Reproduce

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --target wgphysics_mesh_orbits -j
python3 tools/mesh_observable_audit.py --output out/mesh-observable-audit.json
diff -u data/d4-mesh-observable-audit.json out/mesh-observable-audit.json
python3 tests/observable_memory.py
```

The checked-in dataset contains all tetrahedron representatives, patches, successors and
features, every octahedron input and feature transition, exact memory probabilities, invariant
bases, and explicit larger-mesh violations. The test suite also uses a closed coarse-observer
control and rejects corrupted successor, feature, and trajectory records.
