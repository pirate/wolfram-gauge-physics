# Angular control of conserved-charge transport

Follow-up: [exact angular orbits and the fast-angular effective
generator](fiber-angular-stars.md) identifies the kernel $\ker B$,
derives a conditional mixing bound independent of base volume, and
computes exact orbit-averaged reaction rates.

The [relative-holonomy reaction switch](fiber-relative-angle.md) is not
confined to a held five-face region. It changes winding-current memory
on the actual periodic shared-link mesh, and its short-time contribution
per face stays nonzero as the volume increases.

There is also an exact nonperturbative comparison: at matched microscopic
clocks, the long-time winding mobility is a **nondecreasing, concave
function of the angular rate**. A separate criterion is needed to prove
that its increase is strict. Our local witness proves a finite-time
response, not that stricter zero-frequency statement.

The new bulk runs contain 313,122,734 changing raw-link updates across
32-, 128-, and 288-face tori. They resolve a small finite-cutoff current
memory reduction under the faster-angular control. The natural-rate
effect remains unresolved on the largest mesh. None of this identifies
the conserved scalar with electric charge, supplies a physical clock,
or derives the geometry on which the model runs.

## Same primitives, actual links, matched clocks

At every rooted three-face fan, run vacancy exchange and both Hurwitz
orientations on its first ordered pair, all twelve relational reactions
on the triple, and the $m=(n-1)/2$ angular channels on that pair. Sliding
the rooted fan covers each rooted pair once. This avoids counting its
first and last pair twice. All loop products use actual common-basepoint
transport and the existing internal-spoke reconstruction.

Write $A$ for the positive old generator and $B$ for the angular one.
The positive generator is $L_\alpha=A+\alpha B$. Old channel rates are
one and each angular channel has rate $\alpha$. On the triangular torus
there are $3F$ rooted fans, so the total proposal rate is
$3F(15+\alpha m)$, including idle slots. The simulations use Poisson
counts and uniform channel choices, implementing those continuous-time
clocks without reducing the old rates when angular channels are added.

The reference measure is independent uniform raw group elements on every
link. It is stationary because all moves are bijections with equally
rated inverses. It is a mixture over conserved charges and components,
not an assertion of ergodicity or a thermal physical ensemble. Its mean
charge density is $3/2-1/n$, equal to $1.3$ for $C_5$.

Winding current uses the dual edges crossed by the written primal links,
with local unwrapped displacements across periodic seams. Coordinates
are the supplied cell coordinates, not a derived Euclidean embedding.
The integrated current agrees with the change in bounded charge
polarization modulo the torus periods. Angular events have identically
zero charge current, even when the underlying link changes.

## A rate-comparison theorem that does not require coarse-graining

Fix one winding direction. Let $j_a(X)$ be the old channel's current
increment and define

$$b(X)=\sum_a j_a(X),\qquad s(X)=\sum_a j_a(X)^2.$$

Both are independent of $\alpha$. With inner product in the common
stationary measure, the established reversible-current variational
identity specializes to

$$\boxed{\sigma_F(\alpha)
=\frac{\langle s\rangle}{2F}
-\frac1F\langle b,L_\alpha^+b\rangle.}$$

Equivalently, it is the infimum over $f$ of

$$\frac1{2F}\left\langle
\sum_{a\in A}[j_a+f(T_aX)-f(X)]^2
+\alpha\sum_{a\in B}[f(T_aX)-f(X)]^2\right\rangle.$$

For each $f$ this is affine and nondecreasing in $\alpha$. Its infimum is
therefore nondecreasing and concave. This applies at every finite volume
and to every winding direction. No positive increase has been assumed.
The general variational method is not new; see the
[earlier derivation and primary reference](triangle-winding-transport.md#the-exact-variational-identity).
The specialization here uses the new angular channels' exactly vanishing
current and the unchanged stationary measure.

For $\alpha>0$, where the common kernel is fixed, let
$u_\alpha=L_\alpha^+b$. Differentiation gives

$$\sigma_F'(\alpha)=\frac1F\langle u_\alpha,Bu_\alpha\rangle\ge0,$$
$$\sigma_F''(\alpha)=-\frac2F
\langle Bu_\alpha,L_\alpha^+Bu_\alpha\rangle\le0.$$

At $\alpha=0$ extra stationary components can appear. One must then
minimize the angular energy over the otherwise arbitrary old component
constants; blindly differentiating a particular pseudoinverse representative
can give the wrong one-sided derivative. The variational monotonicity and
concavity do not have this ambiguity.

The precise obstruction to a strict increase is useful:

$$\boxed{\sigma_F(\alpha)=\sigma_F(0)\ (\alpha>0)
\quad\Longleftrightarrow\quad
\exists f:\ Af=b\ \text{and}\ Bf=0.}$$

Indeed, such a corrector minimizes both forms with zero added cost.
Conversely, equality of the minima forces an old minimizer with zero
angular Dirichlet energy. A nonzero $BAb$ does **not** by itself rule out
this invariant corrector. Thus observable current memory is not sufficient
evidence for a changed asymptotic conductivity.

In the finite-volume $\alpha\to\infty$ limit, the same variational problem
restricts the corrector to $\ker B$. This is an exact fast-angular limit,
not permission to average uniformly over an entire charge field: angular
orbits can be smaller than charge fibers.

## An exact current-specific angular witness

Let $b_3=3b$, keeping supplied displacement arithmetic integral. Since
$b_3$ depends only on the complete charge field, $Bb_3=0$. The hidden
part of $Ab_3$ comes from three-reflection reaction activities:

$$Ab_3=u_{\rm charge}(q)+\sum_g k_g(q)a_g,$$
$$k_g(q)=-2\sum_{\pi\in\mathrm{Perm}(012)}
[b_3(q^{g,\pi})-b_3(q)],$$

with $k_g=0$ unless the three current charges are all one. A single angular
step preserves $q$ and therefore changes this response by exactly

$$\Delta(Ab_3)=\sum_g k_g(q)\Delta a_g.$$

On a $C_5$ torus with 288 faces, one actual angular step changes one raw
edge, preserves the complete charge field and the common drift
$b_3=(34,26)$, but gives

$$Ab_3:(2372,268)\longrightarrow(2364,300).$$

The weighted activity difference $(-8,32)$ agrees with independently
summing the full old generator on both connections. This is a response
of the winding current, not just a local inventory count.

The difference reads only 67 specified links in an unwrapped region of
extent $6\times5$ cells; its read graph has no torus winding. Keeping that
local assignment fixes the same response difference regardless of the
outside links. Under the uniform raw reference its probability is
$(2n)^{-67}$. Translating the same chosen angular slot once per vertex
gives $F/2$ copies. Hence, for the $x$ current on sufficiently large tori,

$$\boxed{\frac{\langle Ab,BAb\rangle}{F}
\ge\frac{8^2}{36}\,10^{-67}>0.}$$

This conservative bound establishes positivity, not a useful magnitude.
More generally, each angular difference of $Ab$ is a finite-range local
function. Under independent uniform raw links its expectation has only
finitely many translated local types. Consequently
$\langle Ab,BAb\rangle/F$ is volume-independent once those supports embed
without periodic identifications. This statement is about a local moment,
not a hydrodynamic or infinite-time limit.

## What the bulk moment actually measures

Set $C_\alpha(t)=\langle b,e^{-tL_\alpha}b\rangle$. Because $Bb=0$, its
first angular-dependent Taylor term is

$$C_\alpha(t)-C_0(t)
=-\frac{\alpha t^3}{6}\langle Ab,BAb\rangle+O(t^4).$$

The stationary winding variance satisfies

$$\mathbb E[J(t)^2]=t\langle s\rangle
-2\int_0^t(t-v)C_\alpha(v)\,dv.$$

Therefore

$$\boxed{\frac{\mathbb E[J_\alpha(t)^2]-\mathbb E[J_0(t)^2]}F
=\frac{\alpha t^5}{60}
\frac{\langle Ab,BAb\rangle}{F}+O(t^6).}$$

The angular primitive thus produces a strictly positive fifth-order
winding-variance response whose coefficient survives volume growth.
It can do so without changing the asymptotic variance slope; the invariant
corrector criterion above remains an independent question.

Independent local sampling estimates the actual $x$-coefficient moment as
$991.1\pm33.6$ on 128 faces and $1040.9\pm36.5$ on 288 faces, using
50,000 independent raw draws and angular slots for each. Errors are one
standard error. Approximately 9% of those draws give a nonzero winding
response. These values concern $\langle Ab,BAb\rangle/F$, not mobility;
the large numerical coefficient is meaningful only with its $t^5/60$
factor and in the short-time expansion's range of validity.

## Actual autonomous bulk evolution

After an exploratory run, the larger calculation used fresh seeds,
512 stationary trajectories for each combination of sides $4,8,12$ and
angular rates $0,1,8$, each for 32 per-channel time units. Initial raw
connections are paired across rates. Each sample evolves shared links
with every channel enabled at its stated rate. The $\alpha=8$ arm is
an explicit rate-control experiment, not an emergent coupling constant.

The exact instantaneous mobility tensor in cell coordinates is

$$\sigma^{\rm bare}=\begin{pmatrix}1.22&0.61\\0.61&1.22\end{pmatrix},$$

independent of angular rate and mesh size. The current-memory subtraction
is much smaller, about $0.028$ in either diagonal component.

For comparison define the half-trace finite-cutoff memory

$$I_\alpha(T)=\frac1{2F}\int_0^T
[\langle b_x,e^{-tL_\alpha}b_x\rangle+
\langle b_y,e^{-tL_\alpha}b_y\rangle]dt.$$

At the recorded diagnostic cutoff $T=0.2$, paired differences are:

| Faces | $I_0-I_1$, in $10^{-4}$ | $I_0-I_8$, in $10^{-4}$ |
| ---: | ---: | ---: |
| 32 | $4.96\pm2.23$ | $7.82\pm2.25$ |
| 128 | $4.38\pm2.23$ | $6.27\pm2.10$ |
| 288 | $2.09\pm2.14$ | $6.45\pm2.11$ |

Errors are paired trajectory-cluster standard errors, not independent
time-window errors. All seven recorded cutoffs, including the noisier
longer ones, remain in the output. The rate-eight comparison is consistent
across these volumes; the rate-one result on 288 faces is not clearly
resolved. These are exploratory finite-cutoff comparisons, not multiple
independent discovery claims or a finite-size scaling fit.

Current correlations were sampled at $dt=0.001$ and integrated with
composite Simpson quadrature. The largest mean trapezoid–Simpson
difference was $6.70\times10^{-6}$, smaller than the reported sampling
errors. That comparison does not bound every quadrature error or the
unmeasured long-time tail. Direct winding-variance windows are also
recorded but are too noisy to resolve these small angular differences.

Since $C_\alpha(t)\ge0$ for a reversible positive generator, truncating
each memory integral omits a nonnegative tail. But the **difference** of
two omitted tails has no known sign here. A finite-cutoff reduction must
not be reported as a certified DC mobility increase.

## Remaining mathematical gap

The angular-to-current connection is now explicit on an expanding base,
and its local bulk coefficient is both provably positive and numerically
resolved. What remains is the low-frequency question: does the exact
current corrector necessarily vary along angular orbits, and does the
corresponding transport increase remain nonzero in the volume limit?
That requires controlling the corrector or the correlation tail, rather
than declaring a particle or a force from local response measurements.

Geometry, cycle fibers, the reaction bank, and stochastic scheduling are
still supplied ingredients. No physical energy, quantum amplitude,
emergent spatial dimension, or molecular bound state has been derived.

Apparatus and records:
[bulk evolution](../tools/fiber_bulk_transport.py),
[current witness and local moment](../tools/fiber_angular_current.py),
[full trajectory statistics](../data/fiber-bulk-transport-confirmation.json),
[paired summary](../data/fiber-bulk-transport-summary.json),
[exact link witness](../data/fiber-angular-current.json), and
[bulk local moments](../data/fiber-angular-current-moments.json).
