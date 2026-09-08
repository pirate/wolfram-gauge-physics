# Minimal inverse-paired collisions and gauge-mediated activation

The triangle-fiber model now has a **separate, automatically scheduled elastic
candidate bank**. Its two new pair rules are derived by exhaustive algebraic
selection, not fitted to a desired trajectory. They preserve every individual
face charge during a collision, but can change relative gauge information that
enables or disables a neighboring reaction. An exact matched-clock comparison
detects their effect on subsequent whole-mesh mean charge.

This extends the model, not the established laws of nature. The fiber and triangular
torus are still supplied; the process is a classical finite-state Markov chain.
No quantum amplitudes, Hamiltonian, binding potential, continuum gauge group, or
molecular geometry have been inserted or derived. The original bank remains unchanged.

## 1. Exhaust all pair bijections before selecting an extension

Use the existing graph-derived $G=\operatorname{Aut}(C_3)\cong S_3$ and its previously
derived charge $q(e)=0$, $q(r)=1$ for reflections, $q(z)=2$ for nonidentity rotations.
Enumerate all bijections $T(a,b)=(a',b')$ satisfying

$$a'b'=ab,\qquad q(a')+q(b')=q(a)+q(b),$$

and equivariance under simultaneous conjugation. We do **not** initially require
$T^2=I$, an orientation condition, a braid relation, or a particular physical effect.

For a fixed product $p=ab$, parameterize the six factorizations by $a$, so
$b=a^{-1}p$. The permissible permutations of $a$ must preserve the factorization's
charge and commute with conjugation by the centralizer $C(p)$. Choose one such
permutation for each conjugacy class of $p$, then transport it to the other products.
Centralizer equivariance makes the choice of transporting frame irrelevant.
Conversely every admissible full pair map restricts to exactly these choices.

The identity, reflection, and rotation product representatives allow respectively
2, 16, and 6 choices. Thus there are **192** laws in this precise search class.
Their permutation orders have counts

$$N_1=1,\quad N_2=47,\quad N_3=2,\quad N_4=16,\quad N_6=94,\quad N_{12}=32.$$

For orientation reversal $J(a,b)=(b^{-1},a^{-1})$, distinguish two conditions:

- $JTJ=T$: just 32 laws, all involutions (including identity).
- $JTJ=T^{-1}$: 144 laws, including both order-three candidates. This is the
  orientation/time-reversal condition already enforced by the C++ engine.

The earlier minimal-involution search used the stronger restriction. The engine
itself already retained true inverse tables; it did not need a relaxed validator.

## 2. The unique smallest non-involutive extension

Among all 192 laws, the smallest non-involution support has six changed ordered
inputs. Exactly two laws attain it, and they are mutual inverses. Identifying
their exhaustively obtained tables afterward gives

$$H(a,b)=(aba^{-1},a),\qquad H^{-1}(a,b)=(b,b^{-1}ab)$$

on **distinct reflection pairs**, with identity action elsewhere. Both have order
three, conserve the ordered boundary product, fix the flat vacuum, and preserve
the charge of each face separately. They obey $JHJ=H^{-1}$, not $JHJ=H$.

The formula is a restriction of the established
[Hurwitz action on reflection factorizations](https://arxiv.org/abs/2001.08238),
not a newly discovered group operation. Restricting it matters: the resulting
table fails the braid relation on **36 of 216** arbitrary triples (as does its
inverse). We do not claim a Yang–Baxter solution or integrable scattering law.

Why did the involution restriction exclude this elastic channel? For a rotation
product, its two charge-two vacancy/rotation factorizations are fixed by $C(p)$,
while the three distinct-reflection factorizations form a free three-element orbit.
Equivariant bijections cannot mix these orbit types. The latter orbit admits cyclic
rotation, but no nonidentity equivariant involution. This also explains why simply
removing involutivity does **not** permit two reflections to convert into a vacancy
and a rotation within this pair-law class.

The new experimental bank is vacancy transport, $H$, $H^{-1}$, and the original
12 three-face reaction rules. All 15 slots and all rooted supports are sampled
uniformly, independently of state. Inverse partners have equal weight. The bank
therefore has detailed balance for uniform raw connections within an invariant
sector; this alone does not prove that sector is connected or mixing.

“Elastic” here means preservation of individual face charges. It does not mean
we have derived particle momenta, an energy-dependent cross section, or a quantum
scattering matrix. Minimal changed-table support is an explicit model-selection
criterion, not a proof that nature selects this bank.

## 3. A collision invisible inside the pair changes a neighboring reaction

The complete simultaneous-conjugacy class of the two colliding loops is unchanged
by either rule. A pair-only observer sees nothing happen. Their alignment relative
to the rest of the connection, however, need not be unchanged.

The checked whole-mesh witness has total charge four on an 18-face torus. One $H$
update on rooted pair 52 changes **one of 27 links**, leaves all face charges intact,
and removes the enabled reaction at fan 48. The complete whole-mesh loop signature
changes, so this is not merely relabeling vertex frames. Applying the inverse
restores the exact original links.

Reading the two preparations with the **unchanged original bank**, the one-attempt
charge distributions have total-variation distance

$$\frac{24}{2\cdot702}=\frac{2}{117}.$$

Their one-step means agree, but their two-step means differ. Every original
operator, including inverse restoration and no-op multiplicity, is checked against
the C++ raw-link engine. The saved witness is an existence example selected by a
bounded search, not an unbiased collision-frequency estimate.

## 4. Automatic evolution changes the response even with the clock matched

Let $L_0$ denote the unnormalized original generator, $E$ the sum of the two elastic
generators, and $D=L_0q$. There are $M=45F$ attempted operators in the candidate
bank. Compare

$$P_0=I+\frac{L_0}{M},\qquad P_1=I+\frac{L_0+E}{M}.$$

$P_0$ is the original bank **padded with two identity pair slots**. It is not the
original $39F$-operator attempt clock. This prevents an artificial speed difference
from being counted as physics.

The scalar field and its original first drift are unchanged by elastic updates:
$Eq=0$ and $ED=0$. The latter follows from the existing exact local moment law:
$D$ depends only on the full charge field, even though higher response does not.
Consequently, for every initial connection,

$$P_1q=P_0q,\qquad P_1^2q=P_0^2q,\qquad
(P_1^3-P_0^3)q=\frac{E L_0^2q}{M^3}.$$

For the witness above the last numerator is nonzero. On face indices
$(4,6,8,9,14,15,16,17)$ it equals

$$(-64,-64,-64,96,192,-128,-128,160),\qquad M^3=810^3.$$

It sums to zero, respecting charge conservation. All elastic first targets and
all baseline readout transitions needed to evaluate the numerator are independently
scanned in C++. This is an exact finite-time mean-response effect, not a Monte Carlo
fit, a long-time transport coefficient, or a continuum-limit claim.

In the two-charge reflection sector, elastic rules still leave the existing
pointwise drift $D=-2L_1q$ intact. At the matched clock, the exact mean heat law
survives. The extension does not manufacture a force between two defects by fiat.

## 5. Complete patch kernel: activation that the mean alone misses

On one actual three-face read graph, there are 49 complete gauge states. The original
boundary-contained bank has four rooted vacancy operators and 12 reaction operators.
Adding both elastic rules on the four rooted pairs gives a **24-operator** kernel.
This local experiment has its own clock; it is not an autonomous patch approximation
to the entire torus.

Every one of its $49\cdot24=1,176$ entries agrees with raw C++ evolution. The transition
count matrix is exactly symmetric. Two previously frozen states join an eight-state
active component, forming a ten-state component; the total component count falls
from 30 to 28. Many other components remain separate, including the all-equal
reflection state. This is not whole-sector ergodicity.

From either newly released state the probability per local attempt of entering the
old active component is $4/24=1/6$, irrespective of which of the two states the
process currently occupies. The first-entry time therefore obeys

$$\Pr(T>n)=(5/6)^n,\qquad \mathbb E[T]=6.$$

The patch's total charge is three; this is a patch condition with an exterior, not
a claim that the closed torus admits an odd total charge.

Importantly, the candidate and identity-padded baseline have **identical mean
charge at every time on this isolated patch**, for all 49 initial states. We verify
that their difference annihilates the first 49 powers of the control transition
matrix acting on charge; Cayley–Hamilton extends this to all powers. Nevertheless,
their mean squared charges differ already at two attempts. For a fixed initial
state the equal means make this a variance difference too.

Thus activation, first passage, and fluctuations are necessary diagnostics: an
unchanged mean is not evidence of unchanged dynamics. Neighboring mesh interactions
permit the nonzero mean response in the whole-mesh example above.

## 6. Frozen homogeneous states are exactly commuting reductions

There is also a whole-mesh result, not just a patch census. On the supplied connected
triangular torus, suppose every face has $q=1$. Under the candidate bank,

$$\boxed{\text{frozen}\ \Longleftrightarrow\ \text{a parallel reflection section exists}
\ \Longleftrightarrow\ |\operatorname{Hol}|=2.}$$

Here $\operatorname{Hol}$ is the subgroup generated by **all** common-root loops,
including handles. A parallel reflection section assigns a reflection $r_v$ to
each vertex with

$$r_v=U_{uv}r_uU_{uv}^{-1}$$

on every edge. This is a gauge-invariant reduction, not equality of labels in
unrelated vertex frames.

**Proof.** If the state is frozen, every rooted pair must have equal reflections:
distinct ones enable $H$ and $H^{-1}$. Adjacent faces around each vertex therefore
have the same vertex-based holonomy $r_v$. The vertex stars are connected cycles.
Cyclically rebasing an incident face across an edge gives the parallel-transport
identity above. Every closed transport then commutes with the root reflection.
Its centralizer in $S_3$ is precisely $\{e,r\}$. Some face holonomy is a reflection,
so the full holonomy image has order exactly two, not one.

Conversely, an order-two based-loop image determines a path-independent parallel
reflection section. Every face holonomy, transported to the root, is its unique
nonidentity element, because each face has charge one. All pairs therefore have
equal based reflections, all reaction triples have three equal reflections, and
there is no vacancy. Every primitive fixes the state. This proves both directions.

The constant-direction census exhausts all 216 assignments on each of sides 3, 4,
and 6. Each size has 108 all-charge-one assignments: 12 reducible frozen examples,
90 already active full-$S_3$ examples, and **six formerly frozen full-$S_3$ examples
that the new rules activate**. Nonconstant local frame tests verify actual parallel
sections, and the tests also check 100 nonhomogeneous reflection-link assignments.
The proof, rather than this bounded census, establishes the equivalence for every
state on the specified mesh class.

There is a useful extension to the entire $Q=F$ sector. A frozen positive-charge
state cannot contain both vacant and occupied faces: connectedness of the dual
graph would give a vacancy-transport update on their interface. Thus a frozen
$Q=F$ state must have no vacant faces, and its total charge forces every face to
have $q=1$. The theorem then applies. Consequently the **full-$S_3$, $Q=F$ sector
has no absorbing states** under the candidate bank, even when a state temporarily
contains both vacancies and charge-two faces.

The existing reactions preserve the subgroup generated by their based loops;
both new pair maps do too. Boundary-fixed reconstruction and its inverse therefore
preserve the full holonomy subgroup up to conjugation. Total charge is also conserved.
A trajectory starting in this sector cannot enter the reducible frozen sector.
With the finite uniform-attempt scheduler it almost surely keeps making updates
indefinitely. This rules out arrest, **not** multiple active components, slow
relaxation, recurrences, or a lack of spatial propagation.

The automatic traces now include the formerly frozen constant-direction connection
$(1,2,5)$ on three sizes. At every recorded changing state, complete-loop subgroup
and charge checks confirm $|\operatorname{Hol}|=6$, $Q=F$, and an enabled pair update.
This is a verified route to non-arresting finite-gauge evolution, not evidence of
a persistent localized object or a continuum phase transition.

## Reproduce and scope of verification

```bash
cmake --build build --target wgphysics_mixed_bank_experiments -j 4
uv run --with numpy --with scipy --with python-flint python tools/triangle_elastic_scattering.py
diff -u data/triangle-elastic-scattering.json out/triangle-elastic-scattering.json
uv run --with numpy --with scipy --with python-flint python tests/triangle_elastic_scattering.py
```

The optional C++ `--pair-bank` protocol inserts a pair-table count after the existing
five header fields, then reads those pair tables followed by the triple tables.
Schedule slots index pairs first, triples second. The legacy protocol and indices
remain unchanged. In the new protocol, the output mode `transport` means **all pair
rules only**, including elastic rules; it is not the old vacancy-only control.
Reverse replay uses the actual inverse table, never another forward application
of an order-three rule.

Tests also exhaust all $6^5=7,776$ assignments on the pair's five read links,
check nonconstant local gauge transformations, and compare C++/Python histories
on sides 3, 6, and 12 for sparse reference, dense random, and formerly frozen
full-$S_3$ initial connections.
Each trajectory has 4,000 uniform attempts, with exact reversal. The sparse larger
boxes see no elastic encounters in these short traces; these runs are implementation
checks, not evidence for a size-independent collision rate.

Next: extend unbiased encounter and residence measurements with this separately
identified bank, retain relative-loop information across encounters, and test for
persistent correlations. Neither a local activation event nor an exact relaxation
law establishes a bound particle or a molecule.
