# Relative-holonomy reaction gates and their continuum limits

The [diffuse refinement limit](fiber-refinement-limit.md) loses finite-time
reactions. An exact conditional calculation identifies a longer-time
possibility and its missing ingredient: an alignment reference that lives
long enough. The gate reduces to a partially absorbing boundary for an
internal angular coordinate. It is not yet an effective equation for the
full moving mesh.

## What is held fixed

Use the already constructed four-face family

$$w_j=(r_1,r_0,r_{j-1},t_j),\qquad 1\le j\le n-1.$$

Its first two reflections supply a relative-angle reference. The last
$R/Z$ pair has fixed product $r_{-1}$ and forms a degree-one angular star.
Under its two rooted angular descriptions, $j$ steps to $j\pm1$ at rate
$2\alpha$, excluding zero; its state graph is the reflecting path with
$N=n-1$ vertices.

The first triple has exactly one adjacent reflection equality if and only
if $j=1$. All twelve marked reaction channels are active there and nowhere
else. Thus its total marked reaction hazard is $\kappa=12$ at the first
endpoint. This follows directly from the actual relational maps, not a
fitted potential or chosen distance threshold.

For this calculation only, the anchors are fixed, old transport and
other reactions are omitted, and the first marked reaction terminates
the angular process. This defines a **conditional killed-angular
operator**, an ingredient of an effective description. It is not the
unconditional first-reaction law of the original full dynamics. In
particular, freezing these anchors is not evidence of a stable composite.

## Exact finite-fiber waiting time

With $\mathcal L_N$ the positive path Laplacian, the killed operator is

$$K_N=2\alpha\mathcal L_N+\kappa e_1e_1^T.$$

For $\alpha>0$, the mean waiting times solve $K_Nu=\mathbf1$.
Summing the equations gives $\kappa u_1=N$. The reflecting last endpoint
and successive interior differences give

$$\boxed{u_j=\frac N\kappa+
\frac{(j-1)(2N-j)}{4\alpha},\qquad
\overline u=\frac N\kappa+
\frac{(N-1)(2N-1)}{12\alpha}.}$$

This separates reaction-limited and angle-search-limited waiting. The
uniform instantaneous hazard is $\kappa/N$, but its reciprocal is not the
mean first waiting time unless angular mixing is sufficiently fast.
For fixed $\alpha$, the mean grows as $N^2/(6\alpha)$, not merely $N$.

Even starting at the active endpoint gives $u_1=N/\kappa$ for every
$\alpha>0$. Rare escapes with long return times contribute to the mean.
This does not extend continuously to $\alpha=0$, when the endpoint wait
is $1/\kappa$ and all other states never react.

## A continuum reaction boundary without rescaling reaction clocks

Choose the explicit asymptotic regime $\alpha=cN$ and observe time
$t=N\tau$. The coefficient $c>0$ is a rate-control parameter, not a
derived physical coupling. Set $x_j=(j-1/2)/N$. With uniform normalized
inner product, the accelerated operator has quadratic form

$$\langle f,NK_Nf\rangle_N
=2cN\sum_{j=1}^{N-1}(f_{j+1}-f_j)^2+\kappa f_1^2.$$

For smooth sampled functions this tends to

$$2c\int_0^1|f'(x)|^2\,dx+\kappa|f(0)|^2.$$

The bulk generator is consequently $2c\partial_x^2$, with a reflecting
endpoint at one and a partially absorbing (Robin) endpoint at zero:

$$\partial_\tau v=2c\partial_x^2v,\qquad
2c\,\partial_xv(0,\tau)=\kappa v(0,\tau),\qquad
\partial_xv(1,\tau)=0.$$

Here $v$ is a survival probability, not a wave function. The boundary
coefficient comes from the original twelve unit-rate reaction clocks.
No reaction-rate multiplier or force term was added. The coordinate is
relative fiber angle; it is not particle separation in physical space.

The form limit also follows at the spectral level. Piecewise-linear
interpolation turns bounded discrete energy into an $H^1$ bound;
one-dimensional compactness and convergence of boundary traces give the
lower-bound direction. Sampling smooth functions gives recovery sequences.
The finite-dimensional min-max principle then identifies each fixed
limiting eigenvalue of this Sturm-Liouville problem.

Its principal decay rate is

$$\lambda_\mathrm{Robin}=2c k^2,\qquad
k\tan k=\frac\kappa{2c},\qquad 0<k<\frac\pi2.$$

The exact finite waiting-time formula independently gives

$$\frac{u_j}{N}\longrightarrow
\frac1\kappa+\frac{2x-x^2}{4c},\qquad
\frac{\overline u}{N}\longrightarrow\frac1\kappa+\frac1{6c}.$$

Thus three regimes are already distinguished by the finite formula:

- $\alpha/N\to0$: angular search dominates the uniform mean wait.
- $\alpha/N\to c\in(0,\infty)$: reaction and search survive together,
  giving the Robin problem on times of order $N$.
- $\alpha/N\to\infty$: the leading uniform mean wait is $N/\kappa$.

Only the critical-regime spectral limit is established above. The mean
asymptotics alone do not prove an exponential waiting-time law in the
last regime.

## Finite anchor separation gives an interior reactive point

There is a singularity in the endpoint example that must not be hidden:
$r_1$ and $r_0$ are only one cycle step apart, so their normalized relative
angle approaches zero. It is a boundary-layer preparation, not two
anchors at fixed distinct continuum angles.

The construction extends without that restriction. For any
$1\le d\le n-1$, use

$$w_j^{(d)}=(r_d,r_0,r_{j-d},t_j).$$

The last pair has fixed product $r_{-d}$ and the same angular path.
The marked triple reacts exactly at $j=d$, still at rate $\kappa=12$.
The condition is a relative-holonomy equality. Simultaneous conjugation
can reverse both the angular coordinate and the gate; it leaves matching
and the spectrum unchanged. No absolute gauge-frame angle is observable.

Now $K_N=2\alpha\mathcal L_N+\kappa e_de_d^T$. Its exact solution is

$$u_j=\frac N\kappa+\frac1{4\alpha}
\begin{cases}
d(d-1)-j(j-1),&j\le d,\\
(j-d)(2N-j-d+1),&j\ge d.
\end{cases}$$

Consequently

$$\overline u=\frac N\kappa+
\frac{d(d-1)(2d-1)+(N-d)(N-d+1)(2N-2d+1)}{12\alpha N}.$$

When $(d-1/2)/N\to a\in(0,1)$, the same critical scaling gives form

$$2c\int_0^1|f'|^2\,dx+\kappa|f(a)|^2.$$

The continuum process reflects at both ends, is continuous at $a$, and
has the reactive derivative jump

$$2c[f'(a^+)-f'(a^-)]=\kappa f(a).$$

Its principal rate is $2ck^2$, where

$$k[\tan(ka)+\tan(k(1-a))]=\frac\kappa{2c},\qquad
0<k<\frac\pi{2\max(a,1-a)}.$$

The limiting uniform mean is

$$\boxed{\lim\frac{\overline u}{N}
=\frac1\kappa+\frac{a^3+(1-a)^3}{6c}.}$$

In particular, central alignment has search contribution $1/(24c)$,
one quarter of the endpoint contribution. This describes classical
absorption at an internal angle-matching point, not a quantum delta
potential and not spatial attraction. Fixed finite-separation anchors
admit the mechanism; their lifetime in the full dynamics is still missing.

## Why anchor lifetime is the decisive unresolved variable

An elementary competing-clock calculation quantifies the limitation.
Suppose, as an explicitly auxiliary model, the reference is destroyed
at an independent exponential time of rate $\gamma>0$. Let

$$R_\gamma=(\gamma I+2\alpha\mathcal L_N)^{-1},\qquad
g_\gamma=(R_\gamma)_{dd}.$$

A rank-one inverse formula and
$\mathbf1^TR_\gamma=\mathbf1^T/\gamma$ give the probability of reacting
before that loss, starting uniformly:

$$\boxed{\Pr(\text{reaction before reference loss})
=\frac{\kappa}{N\gamma(1+\kappa g_\gamma)}
\le\min\left\{1,\frac{\kappa}{N\gamma}\right\}.}$$

This upper bound is independent of angular speed. A reference with an
order-one lifetime does not rescue a uniformly searched gate in the
refinement limit. In this auxiliary model, a nonvanishing probability
requires a reference lifetime at least of order $N$, or a nonuniform
initial concentration near alignment.

The actual full-model loss process has not been shown to be independent,
exponential, or even described by one scalar rate. That remains a real
gap, not permission to insert such a clock. What this calculation supplies
is the scale a full-model persistent-alignment mechanism must confront.

The bound follows one marked reference. It does not rule out an order-one
total reaction count on times of order $N$ from the continual turnover of
many short-lived references. Persistent references and repeated renewal
are distinct possible longer-time mechanisms; a full-model kinetic limit
must account for that distinction instead of inferring global inactivity
from the lifetime of one pair.

The [released full-dynamics calculation](fiber-reaction-bursts.md) now
observes this distinction directly: rare spontaneous reactions occur in
bursts carried by local loop-flatness constraints, including reactions
outside the original fan. No reference holonomies are frozen in those runs.

## Finite operator observations

Direct tridiagonal solutions use $n=9,33,129,513,2049$ and
$c=0.1,1,10$, always $\kappa=12$. For the endpoint gate at $n=2049$:

| $c$ | Exact scaled mean wait | Continuum mean | Scaled principal rate | Robin rate |
| ---: | ---: | ---: | ---: | ---: |
| 0.1 | 1.748779 | 1.750000 | 0.477666 | 0.477437 |
| 1 | 0.249878 | 0.250000 | 3.644047 | 3.642586 |
| 10 | 0.099988 | 0.100000 | 9.943546 | 9.942347 |

Additional sequences at $c=1$ place the limiting gate at $a=1/4$ and
$a=1/2$. For the central gate at $n=2049$, the scaled mean wait is
$0.12500002$ against limit $0.125$, and the scaled principal rate is
$7.8129543$ against the derived point-reaction value $7.8129580$.

The linear solve agrees with the closed-form mean to relative error below
$9\times10^{-12}$ across the saved cases. The competing-clock resolvent
and rank-one formula agree numerically. At $c=1$, $N=2048$, and
$\gamma=1$, the reaction-before-loss probability is $0.0049403$, below
the upper bound $0.0058594$.

This is a mathematically specified route by which a discrete exact-match
reaction becomes a continuum point or boundary reaction, conditional on a
surviving reference. It does not demonstrate that the original model supplies that
reference, nor does it establish attraction or a molecule.

Apparatus: [fiber_alignment_boundary.py](../tools/fiber_alignment_boundary.py).
Observations: [fiber-alignment-boundary.json](../data/fiber-alignment-boundary.json).
