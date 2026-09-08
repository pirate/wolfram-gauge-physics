# Gauge-complete patch observables and exact gluing

We can now represent the gauge state of an overlapping five-face region without
retaining its raw vertex-frame labels, and evolve that state exactly under every
existing primitive whose read support is contained there. The important correction
is **relative alignment**: two individually complete patch descriptions do not,
by themselves, determine how the patches interact.

The result is a lossless finite gauge-coordinate construction, a 1,393-state kernel,
and an explicit mean-charge witness for the missing alignment. It is not a new force,
a quantum amplitude, or an emergent spatial geometry.

![Actual overlapping patches and their alignment-dependent charge response](images/triangle-patch-observer.png)

## 1. Seven word charges determine a three-loop gauge orbit

The triangular fiber derives $G=\operatorname{Aut}(C_3)\cong S_3$. Its scalar charge
$q$ distinguishes the three conjugacy classes: identity, reflection, and nonidentity
rotation have charges $0,1,2$ respectively. For three loops $(a,b,c)$ at a **common
basepoint**, consider

$$\mathcal W(a,b,c)=
\big(q(a),q(b),q(c),q(ab),q(ac),q(bc),q(abc)\big).$$

Products use the existing multiplication convention, $xy=x\circ y$; algebraic word
order is not silently reversed into path traversal order. All these quantities are
invariant under simultaneous conjugation. Enumerating all 216 triples proves that
$\mathcal W$ distinguishes all **49 simultaneous-conjugacy orbits**.

Each of the four added product words is necessary **within this seven-word
catalogue, with the three individual charges retained**. For every omitted word,
the dataset contains two nonconjugate triples that agree on all six remaining
measurements. This is not a claim of minimality among arbitrary encodings.

For $r$ common-based loops, the charges of all increasing-index subwords of lengths
one through three give complete signatures in the exhaustive censuses at $r=3,4,5$:
7, 14, and 25 words distinguish 49, 251, and 1,393 orbits. Independent Burnside counts
using the actual centralizers give

$$|G^r/G|={6^r+3\,2^r+2\,3^r\over6}.$$

### Linear-size coordinates, rather than a cubic list of words

There is also a constructive invariant representation for arbitrary tuple length.
The anchors below are selected by first occurrence, so their positions are themselves
gauge invariant. The code checks reconstruction against direct conjugacy orbits.

- With rotations but no reflections, choose the first rotation $z$. Each other
  rotation is distinguished as $z$ or $z^{-1}$ by $q(zv)$.
- With reflections but no rotations, retain their equality partition, labeled in
  order of first occurrence. Conjugation permutes the three reflections as the full
  permutation group, so this partition is complete. At most two earlier reflection
  representatives are needed for comparisons.
- With both types, choose the first rotation $z$ and reflection $r$. Rotation labels
  again use $q(zv)$. For another reflection $v$, $q(rv)=0$ identifies $v=r$; otherwise
  $q(rvz)\in\{0,2\}$ distinguishes the two remaining possibilities.

Why is the mixed case complete? Conjugation acts transitively and freely on the six
ordered choices of a nonidentity rotation and a reflection. Once $(z,r)$ is fixed to
a reference pair, every remaining element is uniquely determined by the stated
labels. In the one-type cases the residual stabilizer acts trivially on the recovered
tuple. These are coordinates in the **derived finite group**, not supplied continuum
phases.

`anchored()` and `reconstruct()` therefore use $O(r)$ operations and storage. The
all-short-word list is useful for small-patch interpretation, but is not the intended
representation for a large tuple. Tests reconstruct 100,000-loop tuples, exercise all
six common frame changes, and exhaust every raw tuple through rank five. No numerical
tolerance enters these checks.

## 2. Local completeness does not imply completeness across an overlap

Take actual rooted fan patches 0 and 2 on the supplied side-six torus. They have one
face in common. Their union has seven vertices, eleven links, and five independent
cycles:

$$r=E-V+1=11-7+1=5.$$

There are $6^{11}=362{,}797{,}056$ raw link assignments, $6^5=7{,}776$ tuples after
fixing a spanning tree, and **1,393 gauge orbits** after quotienting the remaining
common frame. These are configuration counts for a small region, not a claim to
have simulated hundreds of millions of spatial nodes.

Describing the two fans by their separate complete 49-state orbits yields only
**845 compatible pairs**. The number of three-loop orbits with a specified shared
loop class is 11, 20, or 18, hence $11^2+20^2+18^2=845$. The forgotten relative
alignment is real information:

$$1393=432\cdot1+305\cdot2+99\cdot3+9\cdot6.$$

Thus 432 compatible pairs have one gluing, 305 have two, 99 have three, and nine have
six. Choosing an arbitrary canonical representative for each fan independently
would select one of these possibilities without justification.

## 3. Retain exactly the relative alignment: a double coset

Let $A,B$ be representatives of the two patch orbits. First align their shared
holonomy to the same element $h$. Define

$$C_A=\{g:gAg^{-1}=A\},\qquad
C_B=\{g:gBg^{-1}=B\},\qquad
C_h=\{g:ghg^{-1}=h\}.$$

The allowed relative alignments are classified by

$$\boxed{C_A\backslash C_h/C_B}.$$

To see this, after fixing $A$, the other patch may be rotated by an element of
$C_h$ without disturbing the shared loop. Its own stabilizer $C_B$ acts on the right;
the remaining freedom preserving $A$, namely $C_A$, acts on the left. Quotienting
both removes redundant choices, while preserving inequivalent alignments.

The implementation explicitly enumerates those double cosets for **every one of
the 845 compatible pairs**. Their sizes independently reproduce the fibers of the
1,393-state raw-connection census. The coordinate

$$(\text{orbit of }A,\ \text{orbit of }B,\ \text{double-coset index})$$

is a bijection onto the union's complete gauge quotient. No freely chosen physical
frame is introduced. For the displayed witness, $h$ is identity, $C_A$ is trivial,
and $C_B$ has order two, leaving three possible relative alignments.

The role of stabilizers and double cosets in combining symmetry data is established
mathematics; see [Baez's discussion of groupoidification and double cosets](https://math.ucr.edu/home/baez/groupoidification/).
The contribution here is its explicit, independently checked realization for these
microscopic gauge rules, not the invention of gluing theory.

### Different basepoints require actual transport

If the fan roots differ, take the explicit link $K$ from the first root to the second
along their shared triangle. Transport the second patch's loops back by

$$b\longmapsto K^{-1}bK.$$

Only then form the common-based gluing coordinates. This makes the shared face
holonomy exactly equal in both tuples and respects arbitrary local gauge frames.
The code records and recomputes that connector, including when a primitive changes
it. It does not assume distant frame identification or path independence.

Tests exhaust all 1,393 union states for each of the twelve one-face overlaps read
by the chosen actor's written links, plus a reverse-connector case. Additional
nonconstant-frame tests and a deliberately omitted-connector control check that
transport is doing necessary work. These tests certify this one-face-overlap
construction; larger intersections and arbitrary patch covers are not yet covered.

## 4. An exact boundary-contained transition kernel

The read-support union contains eight rooted vacancy operators and three rooted
reaction fans, each with twelve rules. **All 44 existing operators** are retained,
including their no-ops. Applying raw link surgery, observing the resulting union
state, and enumerating all initial gauge orbits gives 44 tables of length 1,393.

Every table is a permutation involution. The tables preserve union charge and orbit
stabilizer size. Orbit-stabilizer weights also reconstruct the complete $6^{11}$ raw
assignment count; unequal gauge-orbit sizes have not been silently discarded.

The C++ verification checks **all $44\times1393=61{,}292$ transitions**, using actual
intermediate raw links in operator/inverse scans. It checks the entire raw reset
after each inverse, not merely the observed charge or orbit. Nonconstant local frame
changes give the same quotient transitions. Spanning-tree construction reuses the
existing connection-forest observer.

This is an exact reusable kernel for that bounded region under its contained bank.
It does **not** make the union autonomous under the full-mesh scheduler: an outside
operator can read or modify its boundary. Applying these kernels throughout an
evolving mesh requires explicit overlapping context and consistent gluing data.

## 5. The missing alignment affects later mean charge

Two raw witnesses have identical full charge fields and identical complete orbits
for both selected fans. In the union's tree-based five-loop coordinates they are

$$(z,e,r,r,e),\qquad(z,e,r,s,e),$$

where $z$ is the encoded nonidentity rotation and $r\ne s$ are reflections. An actor
reaction gives different relative-reflection states in the second fan, despite the
two initial separate patch descriptions being identical. The distinction remains
after uniform averaging over the actor's twelve reactions; it is not solely a
particular rule-label effect.

For a controlled evolution that chooses uniformly among the **44 contained
operators**, exact integer powers of the kernel give identical mean charge through
attempts zero, one, and two. At attempt three, in face order $(0,1,10,70,71)$,

$$\mathbb E[q_3\mid U_1]-\mathbb E[q_3\mid U_2]
={1\over44^3}(24,-24,24,-48,24).$$

The numerator sums to zero. Independent enumeration of unquotiented raw paths gives
the same result and checks total path multiplicity $44^n$ at each step. The denominator
is the controlled bank's attempted-update clock, **not** the full torus clock; no
rescaling to a full-mesh physical prediction is asserted.

The word observables also resolve the previous
[charge-plus-activity counterexample](triangle-charge-response.md#6-current-activity-is-insufficient-an-exact-three-attempt-witness):
six patches already differ in products of two reflection holonomies before a
three-reflection encounter becomes active. This identifies pre-encounter relational
information behind that earlier failure, rather than inserting a fitted memory term.

## Reproduce and extend

```bash
uv run --with numpy --with scipy --with python-flint python tools/triangle_patch_observer.py --output out/triangle-patch-observer.json
diff -u data/triangle-patch-observer.json out/triangle-patch-observer.json
uv run --with numpy --with scipy --with python-flint python tests/triangle_patch_observer.py
uv run --with numpy --with matplotlib python tools/plot_triangle_patch_observer.py
```

The next implementation boundary is a consistent cover of an evolving region:
retain shared-frame stabilizers and overlap compatibility, and update the affected
read-support unions without losing boundary information. The present result supplies
an exact local kernel and counterexamples that any proposed whole-mesh reduction
must survive. It establishes neither quantum interference nor particle binding.
