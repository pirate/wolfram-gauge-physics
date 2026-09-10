# Longer gauge-memory tails do not yet give stronger collective response

The same primitive rules on three-, five-, and seven-face triangular strips
produce increasingly slow hidden modes. But the slowest mode couples more
weakly to the charge-square response. Summing the **entire** corresponding
memory spectrum gives a normalized memory duration of about two to three
attempts per face, not the 56 or 152 suggested by the longest lifetime alone.

This is a finite-region result, not a bulk limit, a particle, or evidence of
binding. It identifies a concrete distinction that further searches must
retain: lifetime, observable weight, and spatial response are separate facts.

## Identical primitives and a comparable clock

Each region is an odd, one-cell-wide triangular disk with the existing
derived $S_3$ fiber, total charge equal to its face count $F$, full $S_3$
holonomy image, and a held reflection boundary holonomy $p$. The updates
are exactly the existing vacancy, elastic Hurwitz/inverse, and twelve
reaction rules whose complete read supports lie in the region.

Every strip vertex is on its boundary. An outer tree is fixed to identity
and its closing edge to $p$. All $F-1$ internal edges are enumerated before
filtering by the conserved sector. Simultaneous conjugation by $p$ acts
freely on these framed states. Quotienting this residual pair is legitimate
for the internal charge dynamics considered here, but would lose alignment
information needed for coupling the strip to an exterior gauge reference.

Let $C$ be the integer transition-count matrix and $M=18F-30$ the number
of contained operator slots. Set

$$L=MI-C,\qquad P=I-\frac{L}{45F}.$$

Padding the omitted slots with no-ops keeps the same per-slot proposal
intensity as the full mesh. Time is $\tau=t/F$, in attempts per face.
Renormalizing each strip to only its surviving slots would change this clock.

## Weight the hidden lifetime by the actual observable

Let $\Pi$ be stationary conditional expectation given the **entire** strip
charge field, and $Q=I-\Pi$. Here $Q$ denotes a projector, not total charge.
On the hidden subspace define

$$H=(QLQ)|_{\operatorname{ran}Q},\qquad
D=I-\frac{H}{45F},\qquad b_i=QLf_i,\quad f_i(X)=q_i(X)^2.$$

The enumerated sector is connected, so $H$ is positive definite. Also
$L\leq2MI$ and $2M<45F$, hence $0<D<I$. The normalized trace of the
charge-square hidden-force memory is

$$c_j=\frac{\sum_i b_i^T D^j b_i}{\sum_i\|b_i\|^2}
=\sum_\alpha w_\alpha
\left(1-\frac{\kappa_\alpha}{45F}\right)^j,$$

$$w_\alpha=\frac{\sum_i|v_\alpha^Tb_i|^2}{\sum_i\|b_i\|^2},
\qquad Hv_\alpha=\kappa_\alpha v_\alpha.$$

These weights are nonnegative and sum to one over the complete hidden
eigenbasis. They are spectral weights of an observable, not probabilities
of particles or events. Charge squares probe second moments; they are not
themselves variances for a general evolving distribution.

The full, discrete left-sum memory area has a direct expression:

$$\boxed{\mathcal T=\frac1F\sum_{j=0}^\infty c_j
=45\frac{\sum_i b_i^TH^{-1}b_i}{\sum_i\|b_i\|^2}
=45\sum_\alpha\frac{w_\alpha}{\kappa_\alpha}.}$$

Solving $Hz_i=b_i$ includes **all** hidden modes without truncating the
spectrum. The exponential lifetime of one mode is instead

$$\tau_\alpha=-\frac1{F\log(1-\kappa_\alpha/(45F))}.$$

These two time summaries need not track one another.

| Faces | Framed / quotient states | Charge patterns | Slowest charge-coupled lifetime | Its charge-square force weight | Full memory area $\mathcal T$ |
| --- | --- | --- | --- | --- | --- |
| 3 | 20 / 10 | 7 | 3.0446 | 100% | 3.2143 |
| 5 | 560 / 280 | 51 | 55.8836 | 0.23928% | 2.1748 |
| 7 | 16212 / 8106 | 393 | 151.9288 | 0.09914% | 2.4338 |

Times use the padded attempts-per-face clock. For three faces the hidden
rates are $12,14,16$, but only rate $14$ couples to the charge field;
the slower rate $12$ is dark. Here $\mathcal T=45/14$ exactly.
The larger-region values are numerical evaluations of fully enumerated
finite operators, not exact symbolic spectral formulas. Linear-solve
relative residuals are below $10^{-11}$ and the reported eigen-residuals
below $1.2\times10^{-11}$.

The two slowest charge-coupled modes together contribute only about
0.1341 and 0.1514 to the full memory areas for five and seven faces,
respectively: about 6.2% in either case. The raw force norm squared per
state per face stays near 6.2--6.4, so the comparison is not solely hiding
a large change in the overall force amplitude.

In those same modes, more than 99.7% (five faces) and 99.94% (seven faces)
of the squared state-space eigenvector norm lies in the one-rotation-face
sector. Their charge-square response is strongest in the middle of the
strip; the seven-face response weights are approximately
$(0.0035,0.0653,0.2249,0.4127,0.2249,0.0653,0.0035)$.
That is a response profile, not a spatial probability density or proof
that the underlying gauge mode is localized. All vertices are still
boundary vertices, so central response does not establish bulk behavior.

## A bound on integrated feedback

For any charge-visible observable $f=\Pi f$, put $b=QLf$. Complete the
square in the full microscopic Dirichlet form:

$$\min_{h\in\operatorname{ran}Q}(f+h)^TL(f+h)
=f^TLf-b^TH^{-1}b\geq0.$$

The minimizer is $h_*=-H^{-1}b$. Consequently,

$$\boxed{0\leq b^TH^{-1}b\leq f^TLf.}$$

Thus hidden relaxation can reduce the zero-frequency visible dissipation,
but its integrated feedback cannot exceed the direct Dirichlet energy.
This is a block-operator consequence of the existing reversible rules,
not an added force law or a new general theorem. The resulting Schur
complement is a static response operator, not an autonomous Markov update
on charge patterns.

For the sum over the individual charge-square fields, the measured ratios

$$\frac{\sum_i b_i^TH^{-1}b_i}{\sum_i f_i^TLf_i}$$

are 3.6474%, 2.7579%, and 3.1507% for three, five, and seven faces.
These are observable-specific integrated-feedback fractions, not bounds
on every charge observable's response. They give no growing collective
enhancement in this small strip comparison.

## The slow pair comes from endpoint kinetic traps

Resolving the modes by their **configurations**, rather than their response
profiles, changes the interpretation. For seven faces, the slowest mode has
99.54% of its squared norm in configurations with a single rotation at
an endpoint. Its parallel-reflection configurations carry 98.09%, and
97.15% lies in just eight quiet configurations described below. The next
mode has almost the same support. These are state-space norm fractions,
not probabilities of observing those configurations in equilibrium.

All face holonomies are compared in the common boundary-tree frame:
the chosen tree links are identity, so the actual tree connectors identify
the native face basepoints. This is not a comparison of untransported
holonomies at unrelated frames.

Put the rotation at face 0, a vacancy at face $j$, and the same reflection
on every other face. For $j\geq3$:

- Elastic moves on equal reflections do nothing.
- Every all-reflection reaction fan is inactive.
- No three-face fan contains both the rotation and the vacancy.
- The only changing primitive moves swap the vacancy with a neighboring
  reflection, with two proposal slots per neighboring pair.

There is a reflected copy with the rotation at the other endpoint.
In the fixed boundary sector, each allowed vacancy position gives one
such state after residual conjugation. These configurations are therefore
two exact vacancy paths, not a fitted low-dimensional dynamics.

### Exact first-passage law

For one endpoint set $n=F-3$ and $\ell=j-2\in\{1,\ldots,n\}$.
Stop at the **first exit** from the quiet path, when the vacancy reaches
$j=2$. The positive killed generator is exactly

$$K_n=\begin{pmatrix}
4&-2&&\\
-2&4&\ddots&\\
&\ddots&\ddots&-2\\
&&-2&2
\end{pmatrix}.$$

The missing neighbor at $\ell=0$ is absorbing and the far end is reflecting.
No reaction or elastic update changes a quiet state; this also holds in the
actual primitive transition tables. The smallest killed rate is

$$\kappa_{\rm quiet}=4\left[1-\cos\frac{\pi}{2n+1}\right].$$

Solving the first-passage equation $K_n T=45\mathbf1$ gives, in the same
attempts-per-face clock,

$$\boxed{\mathbb E_\ell[\tau_{\rm exit}]
=\frac{45}{4}\ell(2n+1-\ell).}$$

Starting with the vacancy at the far end gives 67.5 attempts per face for
five faces and 225 for seven. This is a time to enter the reaction-accessible
region, **not** a time to annihilate the rotation or erase all gauge memory.
The quadratic waiting time is derived from the microscopic vacancy moves;
no diffusion equation or attraction was supplied.

### Why the trap survives projection onto hidden gauge information

A quiet-state indicator alone also identifies its charge pattern. To
obtain a genuinely hidden trial function it must be centered within that
complete charge fiber. The correction can be calculated exactly.

For any one-rotation/one-vacancy pattern in an odd $F$-face strip, the
unframed charge fiber has

$$d=3^{F-3}$$

states. One way to count this is to use the ordered boundary-tree loop
factorization. There are $F-2$ reflection factors; choose the nontrivial
rotation in two ways and $F-3$ reflection factors freely. The remaining
reflection is fixed by the boundary product. Dividing by the free residual
conjugation by $p$ leaves $d$ states. Having both a rotation and a reflection
already ensures full $S_3$ image.

Let $E$ inject indicators of the $2n$ quiet configurations, and set $U=QE$.
Their complete charge patterns are all different, so

$$U^TU=(1-d^{-1})I.$$

Let $g_\ell$ count the all-reflection fans in a quiet charge pattern.
Their conditional activity probability is $4/9$: of the 27 reflection
triples, 12 are active, with four for each fixed reflection product.
Fixing the product outside a fan therefore does not change this fraction.
Each active fan has twelve changing reaction slots. Hence the conditional
average reaction-exit count is $16g_\ell/3$.

Vacancy transitions between these charge patterns are independent of their
hidden gauge state. Reaction transitions leave the one-rotation sector,
and elastic transitions do not change charges. Expanding $QLQ$ therefore
gives the exact normalized compression, for each endpoint,

$$\boxed{\frac{U^TLU}{1-d^{-1}}
=K_n+\operatorname{diag}\left(\frac{16g_\ell}{3(d-1)}\right).}$$

Here the equation is blockwise: the full $2n$-column compression is two
identical endpoint blocks. Equivalently, with integer columns $W=dQE$,

$$W^TLW=d(d-1)(K_n\oplus K_n)
+\frac{16d}{3}\operatorname{diag}(g).$$

This identity was evaluated in integer arithmetic against the complete
five- and seven-face operators. It explains the cost of removing the
charge-visible part of the trap, rather than omitting that cost.

| Faces | Quiet-path rate | Minimum centered-trap Rayleigh rate | Actual first two hidden rates |
| --- | --- | --- | --- |
| 5 | 0.763932 | 1.225148 | 0.803806, 0.819510 |
| 7 | 0.241230 | 0.386062 | 0.296052, 0.296563 |

The centered trap is a variational subspace, not an invariant one: the
true modes lower their cost by extending into other configurations.
A centered lowest quiet-path trial has squared projection 92.33% (five
faces) and 97.83% (seven faces) onto the two computed slow hidden modes.
This strongly identifies their mechanism without claiming the restricted
path reproduces the whole dynamics.

### A slow-gap bound without fitting three sizes

For the natural extension of this same odd zigzag strip family, each
endpoint block supplies a hidden trial direction. Since $g_\ell\leq F-4$,
the variational principle gives the following upper bound on **both** of
the first two hidden-generator eigenvalues, ordered with multiplicity:

$$\kappa_2(H)\leq
4\left[1-\cos\frac{\pi}{2F-5}\right]
+\frac{16(F-4)}{3(3^{F-3}-1)}.$$

The right side is $O(F^{-2})$ with an exponentially small centering
correction. This is an analytic upper bound from two boundary traps,
not a measured bulk critical exponent or a matching lower bound.
It also does not establish irreducibility for arbitrary strip length.
It does not guarantee that these slow hidden modes have appreciable
coupling to charge observables. Their small measured charge-square
weights remain essential.

## A compact region changes the traps, not the primitives

As a geometric control, take the six triangles around an interior vertex
and add one adjacent triangle. The face order used in the calculation is
$(36,37,39,52,54,55,21)$ on the existing side-eight mesh. The dual face graph
is a six-cycle with one leaf attached. Its seven outer edges are held;
the center is now an interior vertex whose gauge frame is fixed using
one spoke. No internal physical link is held by that frame choice: after
each primitive update, only an interior gauge transformation restores
the coordinate convention.

This region has exactly the same $Q=7$, reflection boundary sector,
full $S_3$ image, 16,212 framed configurations, 8,106 residual-quotient
configurations, and 393 complete charge patterns as the seven-face strip.
Its different geometry admits 138 contained slots instead of 96. Both
use the same padded $45F$ denominator, so individual primitive slots
have identical proposal rates; the compact region genuinely has more
local interactions, rather than a renormalized scheduler.

| Seven-face geometry | First two hidden rates | Slowest charge-coupled lifetime | Full charge-square memory area |
| --- | --- | --- | --- |
| Strip | 0.296052, 0.296563 | 151.9288 | 2.4338 |
| Six-cycle plus leaf | 0.975137, 1.059176 | 46.0759 | 1.7041 |

The long tail is strongly geometry-sensitive. The second-moment
integrated-feedback fraction is 2.5652% in the compact region, versus
3.1507% in the strip. One interior vertex is not a bulk limit, and the
outer boundary is still held.

### Find the quiet configurations without guessing their gauge pattern

For this comparison, define a quiet configuration operationally: it has
exactly one rotation and **every** elastic and reaction operator fixes it.
Only vacancy slots can change it. Enumerating this condition finds 36
quiet states in the compact region and 92 in the strip. These need not
have the same globally based reflection holonomies.

The two dominant compact quiet components have positive killed generators

$$K_A=\begin{pmatrix}4&0&-2\\0&4&-2\\-2&-2&4\end{pmatrix},
\qquad K_B=\begin{pmatrix}6&-2\\-2&2\end{pmatrix}.$$

Both have smallest rate $4-2\sqrt2\simeq1.171573$, even though their
configuration graphs differ. Component $A$ has its rotation on the leaf;
its vacancy occupies the three far-side cycle faces. Component $B$ has
the rotation on the cycle face opposite the attachment, and its vacancy
alternates between the attachment face and the leaf.

The first slow hidden mode carries 92.44% of its squared norm on component
$A$; the second carries 93.76% on component $B$. Their smaller full hidden
rates include excursions outside the quiet components. The longest mean
first exit from either component is exactly 45 attempts per face, compared
with 225 for the strip endpoint paths. Thus the surviving compact slow
modes also have a concrete vacancy-escape mechanism.

### Identical charge fields, threefold different first-exit times

There is a particularly direct gauge-dependent transport comparison.
In the compact face order above, the charge field

$$q=(0,1,1,1,1,2,1)$$

occurs in nine quiet gauge configurations. One belongs to component $B$.
The other eight are isolated vertices of the quiet-state graph, each
with killed generator $(6)$: any changing vacancy update leaves the
quiet set. The complete charge field alone does not distinguish them.

For the first configuration, two of the six changing vacancy slots move
to the quiet leaf-vacancy state, from which two slots return. The other
four leave the quiet set. For the other eight configurations, all six
changing slots leave it. Therefore

$$\mathbb E[\tau_{\rm exit}\mid\text{component }B,\ q]=22.5,
\qquad
\mathbb E[\tau_{\rm exit}\mid\text{isolated quiet state},\ q]=7.5.$$

These are exact first-passage means from the actual primitive transition
counts. No extra barrier energy or modified hopping law is inserted.
The factor of three measures how gauge alignment changes whether a
vacancy move preserves a locally inactive configuration. It is not a
threefold change in the initial hopping rate, an annihilation lifetime,
or a molecular binding energy. It supplies a charge-blind, gauge-sensitive
kinetic mechanism that any coarse description would have to retain.

## What remains open

The slow tails are real features of the finite primitive dynamics. Their
weak charge-square weights and narrow geometry prevent interpreting them
as robust emergent structures. No exponent can be inferred from these
three sizes. A different complete-charge observable could couple more
strongly than the particular second-moment fields studied here.

The endpoint and compact calculations resolve kinetic bottlenecks underlying
the observed slow modes. The remaining question is whether a comparable
weighted response survives in growing regions with many interior vertices
and evolving boundary data. That question requires retaining
shared-link transport and boundary framing, not gluing independent
unframed strip states or inserting a continuum dynamics.

The [released-boundary follow-up](triangle-released-traps.md) now shows that
the threefold local delay reverses in an aligned dense exterior. Those
preparations activate and spread through the mesh rather than supplying
evidence of a boundary-independent bound object.

Calculation: [triangle_strip_memory.py](../tools/triangle_strip_memory.py).
Mechanism and exact compression:
[triangle_strip_mechanism.py](../tools/triangle_strip_mechanism.py).
The compact comparison uses `Strip(7, shape='compact').spectrum()` and
`kinetic_components(Strip(7, shape='compact'))` from those calculation files.
For the exact projection equation and its initial-hidden-state source,
see [the charge-memory reduction](triangle-charge-memory-kernel.md).
