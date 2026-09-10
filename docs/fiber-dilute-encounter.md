# An exact nonbinding baseline for pair lifetimes

In the one-constraint $Q=F$ sector, a forward reaction creates one identity
face $E$ and one nonidentity rotation face $Z$ in an otherwise reflection
background. Until their first annihilation, the original rules close
exactly on their two positions. This gives a solvable baseline for
distinguishing diffusive persistence from binding.

The result is an exact stopped phase of the fixed-volume, generic
fine-fiber limit. It is not a closure of the full finite-fiber model or
the [two-constraint process](fiber-two-constraints.md), where another
pair can appear and angular motion can alter reaction availability.
The triangular torus geometry and stochastic rule bank remain supplied.

## Position rates counted from the original rules

Let $i$ be the identity position and $j\ne i$ the rotation position on
the honeycomb dual graph. Until annihilation:

- An adjacent $E/R$ exchange has rate $6$: two roots, each supplying
  $H$, $H^{-1}$, and the identity exchange.
- An adjacent $Z/R$ exchange has rate $4$: two roots, each supplying
  $H$ and $H^{-1}$.
- An adjacent $E/Z$ exchange has rate $6$.

The particles cannot occupy the same face. These are derived rates,
not fitted mobilities. There is no additional forward reaction while
the sole inherited relation is the identity-face boundary. Angular
updates preserve the entire position generator in this sector.

Annihilation occurs at rate

$$\kappa(i,j)=4\,\#\{\text{rooted three-face fans containing }i,j\}.$$

On the tori considered here (side at least four), this is $16$ at dual
distance one, $4$ at distance two, and zero farther away. Every motion
transition has the same rate as its reverse. Without annihilation, the
ordered-pair process has the uniform invariant distribution on
$\Omega=\{(i,j):i\ne j\}$: there is no equilibrium positional attraction.

## Birth law and exact mean lifetime

Let $L$ be the positive motion Laplacian and

$$K=L+\operatorname{diag}\kappa.$$

Then $K$ is symmetric, $K\mathbf1=\kappa$, and the killed semigroup is
$e^{-tK}$. The finite connected position graph and nonempty reaction set
make $K$ positive definite. Its mean annihilation times solve

$$Ku=\mathbf1.$$

A forward-reactive fan has twelve equal-rate channels. Each of its six
ordered $E/Z$ placements occurs twice. If the parent fan is sampled
uniformly, the resulting position birth law is therefore

$$\nu(i,j)=\frac{\kappa(i,j)}{\sum_{x\in\Omega}\kappa(x)}.$$

Each face belongs to nine rooted fans, each with two choices for the
other position. Consequently

$$|\Omega|=F(F-1),\qquad
\sum_{x\in\Omega}\kappa(x)=72F.$$

Symmetry now gives the mean without solving the individual hitting times:

$$\boxed{\mathbb E_\nu T
=\frac{\kappa^T K^{-1}\mathbf1}{72F}
=\frac{\mathbf1^T\mathbf1}{72F}
=\frac{F-1}{72}.}$$

Translations and rotations make the parent fans equivalent on this
torus, so the mean also holds for a fixed parent fan with its six layouts
sampled uniformly. The full occupation identities below refer to the
fan-averaged birth law, not an arbitrary fixed initial pair.

The identity uses symmetry and the birth/killing relation, not a chosen
transport speed. It persists under other irreducible symmetric motion
rates with the same birth and killing laws. It is a finite-state balance
identity, not a proposed new force law.

## Stronger result: no lifetime-integrated positional enrichment

For any ordered state $x$, the expected occupation before annihilation is

$$\boxed{\mathbb E_\nu\int_0^T\mathbf1_{X_t=x}\,dt
=\left(\nu^TK^{-1}\right)_x=\frac1{72F}.}$$

Thus all allowed ordered positions receive exactly the same expected
total time. More generally, for a set $A\subseteq\Omega$,

$$\frac{\mathbb E_\nu\int_0^T\mathbf1_{X_t\in A}\,dt}
{\mathbb E_\nu T}=\frac{|A|}{F(F-1)}.$$

This is the ratio of ensemble means, or the occupation fraction obtained
by concatenating independently born excursions. It is not the mean of
the separately normalized occupation fraction of each excursion.

There are three distance-one and six distance-two choices around each
face. Hence the total expected time in the annihilation neighborhood is

$$\mathbb E_\nu T_{d=1}=\frac1{24},\qquad
\mathbb E_\nu T_{d=2}=\frac1{12},\qquad
\mathbb E_\nu T_{1\le d\le2}=\frac18.$$

The fraction of expected lifetime in that neighborhood is exactly
$9/(F-1)$, tending to zero as the supplied volume grows. The growing
mean lifetime is not growing integrated residence near the partner.

This rules out lifetime-integrated positional enrichment for this
specific birth ensemble. It does not claim that the distribution at
each time is uniform; births start close together and killing selects
which paths survive.

## Slowing decay arises without stabilization

Let $Kv_a=\lambda_a v_a$ be an orthonormal eigenbasis, with
$\lambda_a>0$, and write $c_a=\langle v_a,\mathbf1\rangle$. Survival
from the birth law is

$$S_\nu(t)=\nu^Te^{-tK}\mathbf1
=\sum_a w_a e^{-\lambda_a t},\qquad
w_a=\frac{\lambda_a c_a^2}{72F}\ge0,\qquad\sum_a w_a=1.$$

Its hazard $h(t)=-\partial_t\log S_\nu(t)$ is a weighted mean of the
eigenvalues. Differentiation yields

$$\boxed{h'(t)=-\operatorname{Var}_{w(t)}(\lambda)\le0,\qquad
w_a(t)=\frac{w_ae^{-\lambda_at}}{S_\nu(t)}.}$$

The initial hazard is $12$: two thirds of newborn placements are
distance one with rate $16$, and one third are distance two with rate
$4$. Later survivors are increasingly weighted toward the slowly
decaying modes. A falling annihilation rate is therefore expected here
without an attractive interaction or a dynamically stabilized object.
On each finite torus the eventual tail is exponential, governed by the
smallest eigenvalue. No infinite-volume tail law is asserted.

## Exact position calculation across volumes

Translation symmetry reduces the $F(F-1)$ ordered positions to
$2(F-1)$ relative states: two choices of identity sublattice and all
distinct relative rotation positions. This is an exact orbit reduction,
not a continuum or mean-field approximation.

| Faces $F$ | Relative states | Newborn mean lifetime $(F-1)/72$ | Mean from uniform positions |
| ---: | ---: | ---: | ---: |
| 32 | 62 | 0.430556 | 0.751203 |
| 128 | 254 | 1.763889 | 5.107786 |
| 512 | 1,022 | 7.097222 | 29.816739 |
| 2,048 | 4,094 | 28.430556 | 157.945468 |
| 8,192 | 16,382 | 113.763889 | 787.800517 |

The different means reflect different initial ensembles. Newborns start
within reaction range; uniformly placed pairs must generally search for
one another. These five volumes do not by themselves prove an
asymptotic scaling law for the uniform-start mean.

At the largest volume the unreduced position space has 67,100,672 states.
The reduced calculation satisfies the mean-time equation to a maximum
residual below $4\times10^{-11}$. The mean, occupation, and monotone
hazard formulas above are algebraic consequences of the generator, not
fits to these observations.

The next physical question is whether joint constraints create a
persistent, gauge-invariant localization mechanism **beyond** this
symmetric encounter baseline. A long lifetime, nearby birth, or slowing
decay alone is insufficient. No attractive potential has been inserted.

Apparatus and observations:
[exact position generator](../tools/fiber_dilute_encounter.py),
[volume calculations](../data/fiber-dilute-encounter.json).
