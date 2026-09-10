# The exact three-face move group, and what chronology retains

The actual primitives on the smallest reactive disk generate

$$\boxed{G\cong V_{\rm even}\rtimes A_{10},\qquad
V_{\rm even}=\{v\in\mathbb F_2^{10}:\sum_i v_i=0\}.}$$

It has order $2^9\cdot10!/2=928\,972\,800$. Its action is on twenty
boundary-framed classical configurations. The quotient permutes ten
gauge classes, while the binary kernel changes their relative framing.
This identifies the move algebra of this model; it is not a new abstract
group, a quantum state space, or a physical law selecting a schedule.

Two consequences matter for further bottom-up work. The group is perfect,
so its one-dimensional complex phase characters are trivial. Also, making
observables insensitive to **arbitrary overlapping-circuit order** would
collapse this entire local sector to one observable state. These are
constraints on proposed constructions, not claims that quantum behavior
or Wolfram causal invariance is impossible in general.

## Actual links, fixed boundary, and an exterior reference

Use fan 0 on the side-6 mesh. Its three faces read seven edges and write
two internal spokes, edge IDs 2 and 3. Hold the five outer links fixed,
with boundary holonomy $p=1$, a derived reflection. The ordered based
tuple satisfies

$$a_0a_1a_2=p,\qquad\sum_iq(a_i)=3,\qquad
\langle a_0,a_1,a_2\rangle=S_3.$$

There are eight all-reflection tuples: nine with product $p$, minus the
parallel tuple. There are twelve mixed tuples: six orders of $(0,1,2)$
times two rotation orientations. Enumerating the $6^2$ actual internal
link assignments produces exactly these twenty states, without gluing
independently reduced patches.

Keep exterior face 61 at rotation $z=3$, based at the same actual root.
It reads no internal spoke and is unchanged by every contained move.
The joint centralizer of $p$ and $z$ is trivial. This exterior reference
matters: with an entirely commuting exterior, the remaining conjugation
could still be a global gauge transformation.

The boundary involution is

$$\tau(a_0,a_1,a_2)=(pa_0p^{-1},pa_1p^{-1},pa_2p^{-1}).$$

It has no fixed points here, giving ten pairs $\{x,\tau x\}$: four
reflection pairs and six mixed pairs. The jointly gauge-invariant
relative observables

$$O_i=\frac{q(zpa_i)-q(z^{-1}pa_i)}2$$

change sign under $\tau$ with the exterior held fixed. At least one is
nonzero on every state in this sector. Thus the retained sign has an
actual exterior-relative readout, not an artificially identified frame.

All complete read supports inside the disk are included: four rooted
pair patches, each with vacancy and two elastic rules, plus twelve
reaction rules on the fan. These are 24 primitive slots. Each transition
is computed on raw links, leaves the exterior fixed, and commutes with
$\tau$.

## First determine the ten-class quotient

On the twenty framed states, a vacancy or reaction generator consists of
four transpositions; on the ten pairs it consists of two transpositions.
An elastic generator consists of two three-cycles on framed states and
one three-cycle on the quotient. Thus both actions are even, and the
quotient group is a subgroup of $A_{10}$.

Conjugate one actual elastic three-cycle by the actual quotient generators.
The finite conjugacy orbit contains all 240 oriented three-cycles on ten
points: $2\binom{10}{3}$. These generate $A_{10}$. This proves equality
without enumerating all $10!/2$ quotient permutations or appealing to a
numerical group-order estimate.

The calculation exports conjugating words for each star generator
$(0\ 1\ i)$, $i=2,\ldots,9$, as well as every primitive transition.

## Then recover the framing kernel

Any permutation commuting with $\tau$ acts by permuting the ten pairs
and flipping selected pairs. Its sign on the twenty states is
$(-1)^{\text{number of pair flips}}$: exchanging two whole pairs is
even. Since all primitive raw permutations are even, the sign kernel is
contained in $V_{\rm even}$, a nine-dimensional binary space.

Take the adjacent Hurwitz moves on rooted pair supports 4 and 6, in
that chronological order, repeated three times. This six-update word
acts as $\tau$ on the eight reflection states and fixes the twelve
mixed states. It induces identity on the quotient, flipping exactly
four of its ten pairs.

Conjugating this word moves its support through all
$\binom{10}{4}=210$ four-pair subsets. Their flip vectors span precisely
$V_{\rm even}$. For example, two four-subsets sharing three entries
have a symmetric difference of size two, and two-entry vectors generate
the even-weight space. The direct binary rank is nine.

We therefore have the entire kernel, and the quotient is already
$A_{10}$. The group equals the full preimage of $A_{10}$ in the
even-flip permutation group, establishing the displayed semidirect
product. Pair labels and representative choices change coordinates,
not the group or the parity constraint.

This is a constraint on **whole transformations**, not an additional
state charge. It forbids a word that flips exactly one pair and fixes
every other framed state. It does not forbid flipping the alignment of
a particular prepared state: a permitted word may also flip another pair.

## Perfectness and the scalar-phase obstruction

The alternating factor is perfect. Explicitly, on five distinct labels,
the commutator of $(a\ b)(c\ d)$ and $(a\ b)(c\ e)$ is a three-cycle
on $c,d,e$ (orientation depends on the commutator convention). Hence its
commutators generate all three-cycles and all of $A_{10}$.

The binary part is also generated by commutators with that factor.
Choose $v=e_i+e_j$ and a three-cycle $\pi$ fixing $i$ and taking $j$
to $k$. The corresponding commutator has flip vector

$$v+\pi v=e_j+e_k.$$

Such vectors span $V_{\rm even}$. Consequently

$$\boxed{[G,G]=G,\qquad G_{\rm ab}=0.}$$

Every ordinary one-dimensional complex representation factors through
the abelianization, so any homomorphism $G\to U(1)$ is trivial. The
boundary-centralizer sign is therefore **not** a nontrivial scalar phase
character of the full move group.

Perfectness alone does not settle projective representations. The follow-up
[spin-lift calculation](triangle-projective-lift.md) constructs a nontrivial
projective class, but also proves a stabilizer obstruction to encoding
these classical states as rays carrying that class. The full Schur
multiplier remains undetermined. Neither result forbids higher-dimensional
representations or emergent continuous symmetries in a suitable limit;
merely assigning complex coefficients to classical states does not supply
a physical derivation.

## Arbitrary overlap-order erasure would erase the sector

Suppose a state observable $F$ obeyed

$$F(ABx)=F(BAx)\quad\text{for every }A,B\in G\text{ and every }x.$$

Substituting $x=B^{-1}A^{-1}y$ shows that $F$ is invariant under group
commutators. Perfectness then makes it invariant under all of $G$.
The action on these twenty framed states is transitive, so $F$ is
constant throughout the sector.

This hypothesis allows reordering **overlapping** circuit blocks
indiscriminately. It is much stronger than commuting independent updates
with disjoint dependency supports, and is not the definition of causal
invariance in Wolfram-model work. The result says that chronology here
cannot simply be declared gauge redundancy while preserving the local
charge and alignment information.

## A controlled cycle is not spontaneous coherence

An explicit nine-update word cycles three quotient classes and fixes
the other seven. In $(\text{patch},\text{bank rule})$ notation it is

$$[(4,1),(4,0),(0,13),(6,2),(4,1),(6,1),(0,13),(4,0),(4,2)].$$

Rule 0 is vacancy, rules 1 and 2 are inverse elastic partners, and rule
13 is an existing reaction-bank entry. The word is a conjugate of an
elastic generator, not a newly introduced reaction. One raw based-tuple
cycle is

$$(0,2,3)\longrightarrow(0,3,5)\longrightarrow(1,5,5)
\longrightarrow(0,2,3),$$

with charge patterns

$$(0,1,2)\longrightarrow(0,2,1)\longrightarrow(1,1,1)
\longrightarrow(0,1,2).$$

All twenty inputs were replayed through that word on shared links. Six
framed states move, fourteen are fixed, and the boundary and exterior
reference stay fixed. The displayed arrows each stand for nine
primitive updates. This bounded reaction cycle is not transport of a
conserved particle identity or net winding through the torus.

Uniformly preparing its three states and repeating the word produces
an oscillatory normalized autocorrelation $\cos(2\pi n/3)$ for the
centered first-face charge. This follows from an ordinary classical
three-cycle; it is not quantum interference. The oscillation is selected
by the externally chosen circuit and preparation. We have **not**
adopted that circuit as an autonomous physical schedule.

The contrast with the reversible random schedule is therefore precise:
ordered words retain enough control to produce cycles, while averaging
inverse-paired proposals gives the previously derived decay spectrum.
Neither a complex eigenvalue of a chosen permutation nor a periodic
charge trace answers the missing question of what selects physical
chronology or produces quantum probabilities.

Calculation: [triangle_local_control_group.py](../tools/triangle_local_control_group.py).
