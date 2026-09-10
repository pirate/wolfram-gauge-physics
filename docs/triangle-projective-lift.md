# A projective sign exists, but it is not a spinor state encoding

The [twenty-state local move group](triangle-local-control-group.md)
admits a nontrivial projective class through its action on boundary-odd
observables. However, that same class obstructs an equivariant assignment
of spinor rays to the existing classical configurations. This separates
two questions that must not be conflated: whether a transformation algebra
admits spinors, and whether the actual state dynamics supplies them.

The subsequent [raw-history calculation](triangle-rewrite-history.md)
sharpens this distinction: two commuting reaction primitives can both be
inactive on the same raw state while their spin lifts anticommute. Thus
this auxiliary sign is not a phase already carried by the visited raw
configuration history.

Everything below concerns the fixed-boundary, three-face sector with an
explicit exterior reference. The triangular base geometry and triangle
fiber remain supplied inputs. No state-update rule, amplitude, physical
spin, or continuum gauge group is added by this calculation.

## The representation comes from the actual framed states

The derived group is

$$G=V_{\mathrm{even}}\rtimes A_{10},\qquad
V_{\mathrm{even}}=\{v\in\mathbb F_2^{10}:\textstyle\sum_i v_i=0\}.$$

Write the twenty states as ten pairs $\{x_i,\tau x_i\}$. In the real
permutation representation, the ten vectors

$$f_i=\frac{\delta_{x_i}-\delta_{\tau x_i}}{\sqrt2},\qquad 0\leq i<10,$$

span an invariant subspace. The actual primitive tables act on it by
signed permutation matrices. Both the underlying permutation and the
number of flips are even, so this faithful representation lands in
$SO(10)$. These ten coordinates describe observables, **not spatial
dimensions**. The embedding into $SO(10)$ is not a derivation of an
$SO(10)$ physical gauge theory.

Pulling back the standard double cover $\mathrm{Spin}(10)\to SO(10)$ gives

$$1\longrightarrow\{\pm1\}\longrightarrow\widetilde G
\longrightarrow G\longrightarrow1,\qquad
|\widetilde G|=1\,857\,945\,600.$$

This is an auxiliary mathematical construction. For the standard Clifford
construction and half-spin representations, see Peter Woit's
[Clifford Algebras and Spin Groups](https://www.math.columbia.edu/~woit/LieGroups-2012/cliffalgsandspingroups.pdf).
We use the convention $e_i^2=+1$ rather than the negative-square convention
used by default in those notes.

## An explicit nontrivial cocycle

Let $e_ie_j=-e_je_i$ for $i\ne j$, and choose the ordered lift of an
even-weight mask $v$:

$$s(v)=e_0^{v_0}\cdots e_9^{v_9}.$$

Directly moving generators past each other gives

$$s(v)s(w)=(-1)^{c(v,w)}s(v+w),\qquad
c(v,w)=\sum_{i>j}v_iw_j\pmod2.$$

Bilinearity proves the cocycle identity. For even-weight masks,

$$s(v)^2=(-1)^{|v|/2},\qquad
s(v)s(w)s(v)^{-1}s(w)^{-1}=(-1)^{v\cdot w}.$$

Take $v$ supported on $\{1,2\}$ and $w$ on $\{2,3\}$. Their classical
flips commute, but their lifts anticommute. Rephasing either lift by an
arbitrary scalar cannot change this commutator. Thus the cocycle is
nontrivial even with $U(1)$ coefficients; this is stronger than observing
that a particular involution has a lift whose square is $-1$.

On the nine-dimensional space $V_{\mathrm{even}}$, the commutator form
$b(v,w)=v\cdot w$ has rank eight and radical

$$\operatorname{rad}b=\{0,(1,1,\ldots,1)\}.$$

Consequently its lifted subgroup has 1024 elements and center
$\{\pm1,\pm\Omega\}\cong C_4$, where

$$\Omega=e_0\cdots e_9,\qquad\Omega^2=-1.$$

The element $\Omega$ lies over the classical global alignment flip
$\tau$, represented by $-I$ on the odd subspace. It is not the deck
transformation $-1$, which lies over the identity.

This establishes one nontrivial projective class, not the full Schur
multiplier of $G$.

## Primitive relations do not become endpoint amplitudes

Every vacancy or reaction involution has trace six on the ten-dimensional
odd subspace. It therefore has two negative eigenvalues and is a rotation
by $\pi$ in one plane. Both of its spin lifts square to $-1$, although
the raw transformation $T$ satisfies $T^2=I$.

An endpoint function satisfying $F(Tx)=\widetilde T F(x)$ would obey

$$F(x)=F(T^2x)=\widetilde T^{\,2}F(x)=-F(x),$$

so it must vanish. This argument concerns these spin lifts and linear
equivariance; by itself it does not rule out rays or arbitrary rephasings.
The stabilizer argument below supplies the stronger, rephasing-invariant
obstruction.

Elastic generators have trace seven and order three in the odd action.
Their two spin lifts have orders three and six. No physical choice of
lift or phase is selected here.

## The classical stabilizer has no invariant spinor ray

Fix $x_0$. Its stabilizer contains

$$W_0=\{v\in V_{\mathrm{even}}:v_0=0\}\cong(\mathbb Z_2)^8.$$

These transformations flip other pairs while leaving $x_0$ unchanged.
The form $b$ restricted to $W_0$ is nondegenerate: a vector orthogonal to
every even vector on the remaining nine coordinates must have those nine
coordinates all equal; the all-ones choice has odd weight and is excluded.

In particular the two masks $\{1,2\}$ and $\{2,3\}$ both fix $x_0$, but
their lifts anticommute. An equivariant ray assignment would require a
common eigenline of these lifts. If a nonzero vector on that line had
eigenvalues $\alpha,\beta$, anticommutation would imply
$\alpha\beta=-\alpha\beta$, an impossibility for invertible operators.

Therefore there is **no equivariant ray assignment from the twenty
classical states to a projective representation with this cocycle**.
Equivariance here is with respect to the full local move group, with the
deck element acting as scalar $-I$. This is not a no-go theorem for every
projective class, every coarse graining, or every possible dynamics.

The obstruction has an explicit primitive witness. The calculation builds
words of lengths 26 and 28 realizing those two masks, using only the
existing channels. Exact composition of the raw-link-derived transition
tables confirms that both fix state 0, whose based word is $(0,2,3)$,
and commute on all twenty states. Their spin commutator is $-1$.

## Minimal projective dimension, and why mixed states do not fix it

A nondegenerate rank-eight binary commutator form gives a twisted complex
group algebra

$$\mathbb C_c[W_0]\cong\operatorname{Mat}_{16}(\mathbb C).$$

One way to see this is to choose four symplectic pairs for $b$: each pair
generates a two-dimensional matrix factor, and different pairs commute.
Thus every module with this cocycle has dimension divisible by $2^4=16$.
The two complex half-spin representations of $\mathrm{Spin}(10)$ realize
dimension 16, so the lower bound is attained for the full pulled-back
group, not merely for the flip subgroup.

In such a minimal module, any covariant density encoding would satisfy

$$\rho(gx)=U_g\rho(x)U_g^\dagger.$$

At $x_0$ it must commute with all lifts of $W_0$. Their matrix algebra is
full, so its commutant is scalar. Normalization forces

$$\rho(x_0)=I_{16}/16.$$

Transitivity gives that same matrix for every classical state. Thus minimal
spinor density matrices cannot nontrivially encode this sector either.
Larger reducible modules are not excluded by this density-matrix argument.
The number 16 here is a representation-theoretic dimension, not evidence
for a particle generation or a molecular degree of freedom.

## A lower bound on any discrete state extension

There is a further constraint if we try to retain extra history labels.
Suppose a discrete $G$-set $Y$ has an equivariant surjection to the twenty
classical states and an equivariant assignment of rays carrying this
cocycle. For $y$ lying over $x_0$, set

$$K=\{v\in W_0:v y=y\}.$$

Every lift of $K$ preserves the assigned ray. The same anticommutation
argument forces $b|_{K\times K}=0$: $K$ is isotropic. Nondegeneracy of the
eight-dimensional form bounds $\dim K\leq4$. Orbit-stabilizer then gives

$$|W_0y|=\frac{2^8}{|K|}\geq16.$$

All these labels lie over the same classical state. Transitivity gives
at least sixteen labels over every state, hence

$$\boxed{|Y|\geq320.}$$

This is a necessary bound, not a construction or a claim that 320 is
attainable under the full stabilizer. In particular, merely appending one
binary history label per endpoint cannot meet these assumptions. The
bound applies to a discrete extension with an action of $G$ itself;
histories acted on by a different group or groupoid require a separate
analysis. It is not a count of already observed physical degrees of freedom.

## What remains missing

The primitive transformation group supplies a mathematically nontrivial
spin-cover class, but it does not supply physical amplitudes. Even its
classical endpoints cannot be relabeled as rays in that class. A proposed
history-based mechanism must explain which distinctions are retained,
which history relations are physically identified, and why the resulting
structure carries this cocycle rather than one inserted by hand. A
probability and interference law would still need a separate derivation.

Calculation: [triangle_projective_lift.py](../tools/triangle_projective_lift.py).
It constructs the odd actions from actual primitive tables and performs
exact binary/Clifford calculations; it does not evolve spinor states.
