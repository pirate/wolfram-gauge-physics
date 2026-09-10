# A density-dependent bound on conserved-charge transport

The primitive rates give an instantaneous current fluctuation, but that is
not generally the long-time transport rate. An exact variational calculation
now proves a nonzero correction in reactive sectors. At density one, the
large-volume winding mobility is bounded at least **5.315% below** its
instantaneous-current estimate.

This is a bound on transport of the derived scalar charge on the supplied
triangular mesh. It is not an identification with electrical conductivity,
a fitted diffusion coefficient, or a proof of a hydrodynamic limit.

## Use winding current, not bounded charge displacement

The [primitive current construction](triangle-charge-current.md) gives a
current on exactly the dual edges crossed by each update's written links.
Use the actual unwrapped displacement of each such dual edge in the
supplied cell coordinates. Face centers are

$$r_A=(x+2/3,y+1/3),\qquad r_B=(x+1/3,y+2/3).$$

Across the periodic seam, use the local lifted displacement, not the long
difference between coordinates in the chosen fundamental box. The resulting
edge cochain has nonzero circulation around the torus. Summing its current
increments produces an unbounded integrated winding current $J_N$.

This distinction is essential: in a held finite disk, the corresponding
position-weighted charge change is bounded and its long-time variance per
unit time vanishes. A torus winding current cannot be written as the change
of a single-valued bounded charge polarization.

For the $x$ component define the finite-volume stationary mobility

$$\sigma_F=\lim_{N\to\infty}\frac{\mathbb E[J_N^2]}{2N}
=\lim_{\tau\to\infty}\frac{\mathbb E[J_{F\tau}^2]}{2F\tau}.$$

Here $\tau=N/F$ is the existing attempts-per-face clock. The normalization
uses one unit of volume per face; coordinates are the stated cell coordinates.
Expectation is in the exact uniform raw-connection reference at fixed
total charge, full $S_3$ holonomy image, and at least one reflection face.
No ergodicity assumption identifies this reference with a single trajectory.

## The exact variational identity

Let $M=45F$ and

$$P=\frac1M\sum_aT_a,\qquad L=\sum_a(I-T_a)=M(I-P).$$

For an operator slot $a$ at configuration $X$, let $j_a(X)$ be its winding
current increment, and define

$$b(X)=\sum_a j_a(X),\qquad a_2(X)=\sum_a j_a(X)^2.$$

The inverse-paired bank obeys
$j_{a^{-1}}(T_aX)=-j_a(X)$. With uniform stationary weights this gives

$$\mathbb E\sum_a j_a\,[f(T_aX)-f(X)]=-2\langle b,f\rangle.$$

Completing the square in the corrected current yields

$$\boxed{\sigma_F=\frac1{90F}\inf_f
\left\langle\sum_a[j_a+f(T_aX)-f(X)]^2\right\rangle
=\frac{\langle a_2\rangle-2\langle b,L^+b\rangle}{90F}.}$$

The corrector solves $Lf=b$ on each component. Antisymmetry makes $b$
orthogonal to component constants, so the pseudoinverse is well-defined
even if the stationary reference includes multiple components. Equivalently,
inverse pairing gives the discrete increment correlation

$$\mathbb E[j_0j_\ell]=-\frac1{M^2}\langle b,P^{\ell-1}b\rangle,
\qquad \ell\geq1,$$

which leads to the same long-time variance formula. This does not require
replacing the attempted-update process by continuous-time dynamics.

Variational transport methods are established; see
[Arita, Krapivsky and Mallick, *Variational calculation of transport
coefficients in diffusive lattice gases*](https://arxiv.org/abs/1611.07719).
The result here is the specialization and exact moment calculation for
the derived gauge/fiber rules, not a new general Green–Kubo formalism.

## A strict correction without solving the whole state space

Restrict the trial corrector to $f=cb$. Put

$$V=\langle b^2\rangle,\qquad W=\langle b,Lb\rangle.$$

The minimizing coefficient is $c=V/W$, giving

$$\boxed{0\leq\sigma_F\leq
\sigma_F^{(0)}-\frac{2V^2}{90F W},\qquad
\sigma_F^{(0)}=\frac{\langle a_2\rangle}{90F}.}$$

This is an **upper** bound on mobility, or a **lower** bound on the
correction to the instantaneous estimate. It is not an evaluation of the
full memory correction. Whenever $V>0$, that correction is strictly positive.

The current drift $b$ is charge-field-only, but nonlinear. Elastic moves
leave it unchanged. Thus this first trial captures nonlinear charge memory,
not all hidden-gauge memory; adding the elastic bank does not change its
numerator or energy at the matched proposal clock. Higher correctors, such
as $Lb$, can depend on the actual gauge-dependent reaction activities.

## Exact moments from local charge patterns

The calculation uses $b_3=3b$ and $j_{3,a}=3j_a$ so every coefficient is
integer. Write $V_3=\langle b_3^2\rangle$ and
$W_3=\langle b_3,Lb_3\rangle$. The corresponding mobility correction is
$2V_3^2/(810F W_3)$.

Expand $b_3$ into indicators of ordered two- and three-face charge patterns.
Multiplying these indicators and merging compatible assignments groups
their expectations into the already exact probabilities

$$p_F(c_0,c_1,c_2)
=\frac{\mathbb E[(N_0)_{c_0}(N_1)_{c_1}(N_2)_{c_2}]}{(F)_{c_0+c_1+c_2}}.$$

For the Dirichlet energy, the gauge-dependent births need not be replaced
by an assumed average activity. Reversibility pairs their squared changes
with their reverse deaths. Each $(0,1,2)$ fan has four such death slots.
For any charge-only function $f$,

$$\langle f,Lf\rangle=
\frac12\sum_{\text{rooted vacancy slots}}
\mathbb E[\mathbf1_{\rm one\ vacancy}(\Delta f)^2]
+4\sum_{\text{fans}}\mathbb E[\mathbf1_{\mathrm{Perm}(012)}
(f(q^{111})-f(q))^2].$$

This identity uses the actual stationary raw-connection multiplicities.
It is not a Markov closure on charge fields. For $f=b_3$, only terms
overlapping the changed faces contribute; each is restricted to the old
or new local charges before squaring. The final moment census requires
at most seven distinct charge variables.

Disjoint local current contributions have zero covariance by charge-field
exchangeability and their zero permutation-averaged current. The remaining
variance and energy coefficients are local. For side lengths at least four,
their connected support unions contain at most six dual edges, too few
to make a new periodic identification on these tori. The per-face integer
censuses agree at sides four and six; the side-three energy correctly
retains its additional short periodic identifications.

Independently summing local current squares gives

$$\frac{\langle\sum_a j_{3,a}^2\rangle}{F}
=12p_F(1,1,0)+48p_F(1,0,1)+480p_F(1,1,1).$$

At $Q=2$, $V_3=W_3=0$: the exact dilute exclusion/heat-equation sector has
no current-drift correction. At density one the exact rational evaluations
give:

| Faces | Instantaneous estimate $\sigma_F^{(0)}$ | One-function upper bound | Minimum reduction |
| --- | --- | --- | --- |
| 18 | 0.02518572 | 0.02376033 | 5.6595% |
| 32 | 0.02430600 | 0.02297301 | 5.4842% |
| 72 | 0.02373010 | 0.02245236 | 5.3845% |

These are finite-state expectation calculations with exact integer/rational
arithmetic, not fitted simulation slopes or Monte Carlo estimates.

## An explicit fixed-density bound

For admissible even-charge sequences $Q/F\to\rho\in(0,2)$, the local
charge marginals of the counted reference converge to

$$p_0=Z^{-1},\quad p_1=3z/Z,\quad p_2=2z^2/Z,
\qquad Z=(1+z)(1+2z),\qquad
\rho=\frac{z(4z+3)}{(1+z)(1+2z)}.$$

This follows from the fixed-charge coefficient weights
$3^{N_1}2^{N_2}$; the excluded proper-group and reflection-free terms are
exponentially negligible at interior densities. The positive fugacity $z$
parameterizes the **reference density**, not a new update parameter.

Define the polynomials

$$A=68z^4-108z^3+1013z^2-54z+17,$$
$$B=16z^3+270z^2+17z+3,$$
$$C=544z^8+1560z^7+11712z^6+53070z^5+65314z^4
+27459z^3+3495z^2+306z+34.$$

The local moment census gives

$$\frac{V_3}{F}\longrightarrow\frac{144z^3A}{Z^5},\qquad
\frac{W_3}{F}\longrightarrow\frac{1152z^3C}{Z^7}.$$

Consequently the explicit thermodynamic **limsup bound** is

$$\boxed{\limsup_{F\to\infty}\sigma_F
\leq\frac{2zB}{135Z^3}
\left(1-\frac{3z^2A^2}{BC}\right).}$$

The relative reduction is strictly positive for every interior density.
$B,C$ have positive coefficients; $A>0$ follows, for example, by writing
it as $2z^2(34z^2-54z+27)+(959z^2-54z+17)$, two positive quadratics
with the indicated nonnegative prefactor.

At density one, $z=1/\sqrt2$. The instantaneous estimate tends to
0.023296916, whereas the one-function upper bound tends to 0.022058679:
a minimum reduction of 5.315025%. This is an analytic bound, not a
finite-size extrapolation of measured transport.

## What is still missing

The full mobility has not been computed, nor has its thermodynamic limit
or a hydrodynamic equation been proved. Dividing this bound by a static
susceptibility and calling it a measured diffusion constant would go beyond
the result. A systematic next correction must include gauge-dependent
functions such as $Lb$ and their actual boundary-transported activities.
An [explicit purely hidden-gauge corrector](triangle-gauge-current-corrector.md)
now supplies a provably stricter two-function finite-volume bound. Its
additional moment values, and therefore its numerical improvement, remain
to be evaluated.
Neither this transport result nor a nonzero diffusion coefficient would
by itself establish particles, binding, quantum dynamics, or spacetime.

Calculation: [triangle_winding_corrector.py](../tools/triangle_winding_corrector.py).
