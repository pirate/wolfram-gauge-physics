# Composite fibers: bounded reactions, but an exact spectral threshold

The [odd-prime arity obstruction](cycle-spectral-reactions.md) does **not**
extend unchanged to composite cycle fibers. For an odd prime $p$ and
$n=p^a$ with $a\ge2$, the corresponding sharp results are

$$\boxed{K_{\min}=p,\qquad E_{\min}=2p.}$$

Here $K$ is the number of ordered loop variables in an update and $E$ is
the sum of the same lowest-mode spectral observable

$$\varepsilon(g)=2-\operatorname{tr}U(g),\qquad
\varepsilon(r_j)=2,\qquad
\varepsilon(t_k)=2-2\cos(2\pi k/n).$$

The minimum is over updates that change conjugacy-class populations,
preserve ordered product, are bijective and conjugation equivariant, and
exactly conserve $E$. A sparse orientation-covariant involution attains
both minima simultaneously. The energy lower bound in fact needs neither
bijectivity nor covariance, only spectral conservation and product parity.

**This is conditional mathematics about a candidate diagnostic.** We
have not derived its conservation from the original primitive rule bank,
adopted it as physical energy, or installed these rules as the simulation's
new dynamics. Nor is an arity bound a derived spatial interaction radius.

## Exact population relations

Put $m=p^{a-1}$ and $\zeta=e^{2\pi i/n}$. As before, derive the
automorphisms from cycle adjacency; write their rotations and reflections
as $t_k(j)=j+k$ and $r_b(j)=b-j$ modulo $n$.

There is one reflection class $R$ and $(n-1)/2$ nonidentity rotation
classes $Z_k=\{t_k,t_{-k}\}$. Partition these rotation classes into:

- A special band $B_0=\{Z_{jm}:1\le j\le(p-1)/2\}$, of size $(p-1)/2$.
- For $1\le r\le(m-1)/2$, a generic band
  $B_r=\{Z_{r+jm}:0\le j<p\}$, of size $p$.

Indices of classes are identified modulo sign and $n$. These bands are
disjoint and cover all nonidentity rotation classes. Their spectral sums
are $p$ for $B_0$ and $2p$ for each $B_r$, $r>0$.

Let $x=\Delta N_R$ and $y_k=\Delta N_{Z_k}$. Conservation means

$$2x+2\sum_k y_k-\sum_k y_k(\zeta^k+\zeta^{-k})=0.$$

The minimal polynomial is

$$\Phi_{p^a}(z)=1+z^m+\cdots+z^{(p-1)m}.$$

Our polynomial has degree at most $n-1$, so it vanishes exactly when
its coefficients at $r,r+m,\ldots,r+(p-1)m$ agree, separately for
each residue $0\le r<m$. Indeed its quotient by $\Phi_{p^a}$ has
degree at most $m-1$, and these shifted copies have disjoint supports.

For the nonzero residues, this forces a common integer population change
$s_r$ throughout each generic band. At residue zero, write the common
special-band change as $s_0$. The constant coefficient equation gives

$$2x+p s_0+2p\sum_{r>0}s_r=0.$$

Since $p$ is odd, $s_0=2v$ for an integer $v$. Fixed tuple length and
the equation then give the complete integer energy kernel:

$$\boxed{
\begin{aligned}
\Delta N_e&=v,\\
\Delta N_R&=-p\left(v+\sum_{r>0}s_r\right),\\
\Delta N_Z&=2v &&(Z\in B_0),\\
\Delta N_Z&=s_r &&(Z\in B_r,\ r>0).
\end{aligned}}$$

Ordered-product preservation further requires

$$v+\sum_{r>0}s_r\equiv0\pmod2,$$

because the reflection count can only change by an even integer. This
parity condition is necessary, not a claim that every population vector
has an ordered-product-preserving realization at its smallest support.

The integer energy kernel has rank $(m+1)/2$. Equivalently, exact
conservation of this one algebraic-valued sum preserves
$\varphi(n)/2=(n-m)/2$ independent rational coefficient charges on the
nonidentity populations. Composite refinement permits more resonances
than prime refinement, but still introduces a growing family of exact
population constraints.

## Lower bounds on reaction size and spectral sum

If any $s_r\ne0$, one endpoint must contain at least one element of
each of the $p$ classes in $B_r$. Thus it has at least $p$ slots and
spectral sum at least $2p$. All other contributions are nonnegative.

If every $s_r=0$, a nonzero change instead has $v\ne0$. Product parity
forces $v$ even, so one endpoint contains at least $2p$ reflections,
with spectral sum at least $4p$.

Therefore every population-changing reaction has $K\ge p$ and
$E\ge2p$, regardless of additional spectator loops. Equality is possible.

In particular, **below total spectral sum $2p$, class populations cannot
change under any such update, even if arbitrarily many loops are read.**
This does not prohibit transport, positional rearrangements, or changes
of other observables. It is not a claim of total dynamical freezing.

## Construction attaining both bounds

Define the nonzero residue

$$s=\frac{p^2-1}{4}\pmod p,\qquad 1\le s<p,\qquad
r=s\frac{m}{p}.$$

The condition $a\ge2$ ensures $m/p$ is an integer. Choose the $p$-tuples

$$x=(t_r,t_{r+m},\ldots,t_{r+(p-1)m}),$$

$$y=(t_m,t_m,t_{2m},t_{2m},\ldots,
t_{(p-1)m/2},t_{(p-1)m/2},e).$$

Their spectral sums both equal $2p$. For the products, the sum of source
exponents is congruent to $pr=sm$ modulo $n$, and the sum of target
exponents is

$$2m\sum_{j=1}^{(p-1)/2}j=m\frac{p^2-1}{4}\equiv sm\pmod n.$$

Both tuples have conjugation stabilizer equal to the rotation subgroup
$C_n$: they contain nonidentity rotations and $n$ is odd. Each has two
conjugates, related by inversion of all rotations. Reversed inversion
doubles each orbit. For $x$, equality with its reversal-inverse would
require $2r\equiv m\pmod n$, impossible for $0<r<m$ and odd $m$;
equality with its conjugate inverse would require $m\equiv0\pmod n$.
For $y$, reversal moves the identity from the end to the beginning.

The two extended orbits are disjoint because only the target contains an
identity. Closing $x\leftrightarrow y$ under conjugation and reversed
inversion therefore produces **four disjoint transpositions**, moving
eight tuples and fixing all others. This is a flat-fixed, reversible,
product-preserving, gauge- and orientation-covariant rule.

For $C_9$, the concrete exchange is

$$\boxed{(t_2,t_5,t_8)\longleftrightarrow(t_3,t_3,e).}$$

Both products are $t_6$ and both spectral sums are $6$. This persists
on $C_{27}$ as $(t_6,t_{15},t_{24})\leftrightarrow(t_9,t_9,e)$ and
on $C_{81}$ as $(t_{18},t_{45},t_{72})\leftrightarrow(t_{27},t_{27},e)$.
The angles in these particular witnesses stay fixed under refinement;
they do not approach the identity.

## Exact computations and a shared-link realization

The calculation derives each automorphism group from adjacency and uses
integer polynomial remainders modulo $\Phi_{p^a}$ throughout. It constructs
the sharp exchanges on $C_9,C_{27},C_{81},C_{25},C_{49}$ and computes the
full energy-kernel dimensions as $2,5,14,3,4$, respectively.

For $C_9$, a complete search through all one-, two-, and three-loop
tuples, grouped by exact spectral sum and ordered product, gives:

- No reactive minimal exchange closures at arity one or two.
- Six reactive minimal exchange closures at arity three, all at $E=6$.

These six closures are the complete **minimal transposition-closure**
census at this arity, not a count of all composite bijections.

One exchange is reconstructed on the existing side-four triangular mesh
using the generic shared-link oracle. Only its two internal links change;
all exterior links remain fixed. The actual based face tuple changes as
claimed, the inverse returns every raw link exactly, and the whole-mesh
spectral sum remains $9$. The reaction patch contributes $6$; the other
faces contribute $3$. Thus the witness is not an independent-face update
that ignores neighboring faces or boundary transport. It remains an
isolated raw-link realization, not evidence of equilibration or continuum
physics under an installed evolution bank.

## What this changes about the research direction

The prime-only result could have suggested that refinement always forces
interaction size to grow. Prime-power fibers disprove that extrapolation:
fixed $p$ gives fixed reaction size for arbitrarily large $p^a$.

But the escape has a precise low-spectral-sum limitation. A fixed-size
tuple of rotations whose angles all tend to zero has $E\to0$, and
eventually cannot undergo any class-population conversion of this kind.
Multiplying the observable by a common scale factor changes the numerical
threshold by the same factor, not which exact reactions are permitted.
Large-support collective processes, different observables, additional
primitive degrees of freedom, and composite sizes with several prime
factors are not ruled out by this calculation.

Most importantly, this does not supply the missing physical principle
that would select exact conservation of $\varepsilon$. It tells us what
we would be committing to if we imposed it. We should not mistake the
existence of a smooth-looking spectral diagnostic for a derived dynamical
energy, or quietly impose it to obtain desired field behavior.

The [word-rule analysis](fiber-word-conservation.md) further distinguishes
these finite tables from multiplication-and-inversion primitives: the
sharp exchanges reduce the generated rotation subgroup from order $p^2$
to $p$, which no bijective word map on a finite group can do. Spectral
conservation for a word map on any odd cycle actually forces preservation
of the entire conjugacy-class inventory. Additional fiber structure or
state-dependent relational selection is therefore essential to these
particular reactions.

The roots-of-unity arithmetic has established mathematical background;
see Lam and Leung, [On vanishing sums for roots of unity](https://arxiv.org/abs/math/9511209).
Our derivation above is self-contained for prime powers. The application
to these tuple updates and its sharp constructions have not received an
independent novelty review; no field-wide priority claim is made.

Calculation: [cycle_prime_power_reactions.py](../tools/cycle_prime_power_reactions.py).
Exact tuples, closures, kernel bases, census, and raw links:
[cycle-prime-power-reactions.json](../data/cycle-prime-power-reactions.json).
