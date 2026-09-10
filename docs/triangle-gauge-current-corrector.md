# A purely hidden-gauge correction to winding transport

The first winding-mobility bound used only the nonlinear charge-current
drift. There is now an explicit **gauge-hidden** trial function that strictly
improves that one-function bound in the audited reactive sectors. It is
constructed from the actual elastic and reactive operators, without assuming
a conditional gate-activity average.

The [stratified stationary calculation](triangle-gauge-resolved-transport.md)
now resolves the correction's size on 72- and 128-face tori. A local
operator bound below separately supplies a coefficient that guarantees
improvement without estimating the noisy hidden Dirichlet moment. Neither
replaces the previous numerical 5.315% bulk bound with a sharper certified
bulk number.
The follow-up [disk-cylinder argument](triangle-gauge-bulk-correction.md)
does prove that a strictly positive correction survives the volume limit,
although its explicit guaranteed size is far below the reported precision.

## Construct a hidden function without estimating a conditional mean

Use the positive generators

$$L=L_0+E,\qquad E=\sum_{\text{elastic slots}}(I-T_a),$$

where $L_0$ contains vacancy and reaction slots. Both are self-adjoint and
positive in the exact stationary raw-connection reference. Let $b_3$ be
three times the unnormalized $x$-winding current drift from
[the transport calculation](triangle-winding-transport.md). It depends only
on the complete charge field.

Every elastic primitive preserves that field. Hence, for **every** function
$f(q)$, $Ef(q)=0$. Define

$$u=Lb_3,\qquad \boxed{h=Eu=ELb_3.}$$

For any charge-visible function,

$$\langle f(q),h\rangle=\langle Ef(q),Lb_3\rangle=0.$$

Therefore $\mathbb E[h\mid\text{entire charge field}]=0$ exactly. This
does not assume the elastic moves connect a whole charge fiber: summing
their stationary generator gives zero even when that fiber has many
elastic components.

Since $Eb_3=0$, the same function is $h=[E,L]b_3$. This is ordinary
noncommutation of classical update operators, not quantum interference
or a postulated commutator evolution law.

## Its local expression retains actual gauge transport

Only the all-reflection reaction activities can make $u=Lb_3$ depend on
gauge information beyond charges. For a fan $g$, set

$$k_g(q)=-2\sum_{\sigma\in\mathrm{Perm}(012)}
[b_3(q^{g,\sigma})-b_3(q)]$$

when its current charges are $(1,1,1)$, and zero otherwise. Then

$$u(X)=u_{\rm charge}(q(X))+\sum_g k_g(q(X))A_g(X).$$

The activity $A_g$ compares the three holonomies transported to the fan's
actual common basepoint. For an elastic step $T_a$, charges and the weights
$k_g$ stay fixed, so

$$u(T_aX)-u(X)=\sum_gk_g(q)[A_g(T_aX)-A_g(X)].$$

Summing the negative of these differences over the elastic slots produces
$h$. Thus the hidden corrector measures how elastic rearrangements change
reaction activity relevant to the global winding current. No independent
local gauge quotients are glued together.

## Exact current-specific witnesses

An exact reference draw at $Q=F$ supplies the following elastic pairs.
Each pair has identical complete charge fields and identical $b_3$.
Their positive-generator responses differ:

| Faces | Elastic patch | $b_3$ before/after | $Lb_3$ before | $Lb_3$ after | $ELb_3$ before/after |
| --- | --- | --- | --- | --- | --- |
| 18 | 12 | $-16$ | $-340$ | $-708$ | $2288,-2816$ |
| 72 | 64 | $38$ | $852$ | $860$ | $368,432$ |

For 18 faces, the changed weighted activities are exactly $-192$ and
$-176$, summing to $-368$. For 72 faces, one contribution is $+8$.
These differences also agree with independently summing $Lb_3$ over the
complete microscopic operator bank on both raw connections. The calculation
exports the initial and final link arrays and the reference seed.

The numbers $Lb_3$ are generator responses, not current increments. With
$M=45F$, the next expected drift is $b_3-Lb_3/M$. Thus the pairs are also
explicit witnesses that identical present charge fields and current drift
do not determine the next expected current drift.

## A strictly better two-function variational bound

Use $f=c_1b_3+c_2h$ as the scaled corrector. Define

$$V=\langle b_3^2\rangle,\quad W=\langle b_3,Lb_3\rangle,\quad
S=\langle u,Eu\rangle,\quad T=\langle h,Lh\rangle.$$

The source vector and Dirichlet Gram matrix are exactly

$$\begin{pmatrix}\langle b_3,b_3\rangle\\\langle b_3,h\rangle\end{pmatrix}
=\begin{pmatrix}V\\0\end{pmatrix},\qquad
G=\begin{pmatrix}W&S\\S&T\end{pmatrix}.$$

In particular, $\langle b_3,Lh\rangle=\langle Lb_3,ELb_3\rangle=S$.
The elastic Dirichlet identity gives

$$S=\frac12\sum_{\text{elastic slots}}
\mathbb E[(u(T_aX)-u(X))^2].$$

The exact witnesses prove $S>0$ in their finite-volume reference sectors:
their configurations have positive stationary weight. They also show
$h\ne0$. Since $h$ is in the range of $E$ and $L\succeq E$, $T>0$.
The current drift is orthogonal to component constants. Together with
$\langle b_3,h\rangle=0$, this makes $G$ positive definite on these two
nonzero independent trial functions, so $WT-S^2>0$.

Minimizing over the two coefficients gives

$$c_1=\frac{V}{W-S^2/T},\qquad c_2=-\frac STc_1,$$

$$\boxed{\sigma_F\leq\sigma_F^{(0)}
-\frac{2V^2}{810F\,(W-S^2/T)}
<\sigma_F^{(0)}-\frac{2V^2}{810F W}.}$$

This proves a strict improvement over the specific one-function bound,
using a function orthogonal to **all** charge-only observables. It does
not prove that this two-function ansatz beats the optimal unrestricted
charge-only ansatz. Nor does the finite witness alone quantify an additional
bulk correction by itself. The subsequent
[local disk witness](triangle-gauge-bulk-correction.md) supplies a positive
bulk lower bound, and the [stratified calculation](triangle-gauge-resolved-transport.md)
estimates a useful finite-volume size. A useful numerical bulk value remains open.

## Locality controls the hidden Dirichlet cost

The large configuration space need not appear in an upper bound on $T/S$.
Group the elastic generator by its written primal edge:

$$E=\sum_e E_e,\qquad h=\sum_e h_e,\qquad h_e=E_eu.$$

There are two rooted pair patches per edge. Each contributes
$2I-H-H^{-1}$. The restricted Hurwitz permutation has $H^3=I$, so this
positive operator has eigenvalues $0,3$. Its lift changes just one link:
after three applications the based pair returns, and the fixed surrounding
links uniquely determine that written link. Thus the same order statement
holds on raw connections, not merely on pair classes. Consequently

$$0\preceq E_e\preceq6I,\qquad
\sum_e\|h_e\|^2\leq6\sum_e\langle u,E_eu\rangle=6S.$$

Although $u$ is extensive, $h_e$ is local: the charge-only part cancels,
and only fan activities whose read paths contain $e$ can change. Its read
set $R_e$ includes the pair paths, those fan paths, and the face paths
needed for their charge-dependent weights. These are actual raw-link read
sets, including the common-basepoint connectors.

Let $C_1$ bound the number of these read sets touched by any one primitive's
writes, and $C_2$ bound the number of primitive slots touching any one read
set. Count all rule and rooted-patch multiplicities, including slots that
may happen to be no-ops in a particular state. With $\Delta_a=T_a-I$,
Cauchy--Schwarz and the unitarity of each stationary permutation give

$$\begin{aligned}
T&=\frac12\sum_a\left\|\sum_e\Delta_a h_e\right\|^2\\
&\leq\frac{C_1}{2}\sum_{a,e:\,a\text{ touches }R_e}
\|\Delta_a h_e\|^2\\
&\leq2C_1C_2\sum_e\|h_e\|^2
\leq12C_1C_2S.
\end{aligned}$$

The explicit support census on each of the side-6 and side-8 tori gives
$|R_e|\leq67$, $C_1=75$, and $C_2=2202$. Therefore on these geometries,

$$\boxed{T\leq C S,\qquad C=1\,981\,800.}$$

This is deliberately conservative, not a measured relaxation rate. The
same constants on two sizes do not alone prove that precise numerical
constant for every torus. Bounded support and bounded incidence give a
volume-independent bound on this fixed local lattice family; identifying
its sharp constant is a separate question. This argument is not a
locality claim for arbitrary growing hypergraphs.

### A guaranteed correction without estimating $T$

Keep the exact charge-trial coefficient $c=V/W$ and use

$$f=c\left(b_3-\frac{h}{C}\right).$$

If $U_1=\sigma_F^{(0)}-2V^2/(810FW)$ is the previous bound, its change in
variational cost is exactly

$$\Delta U=\frac{2c^2}{810F}
\left(-\frac{2S}{C}+\frac{T}{C^2}\right)
\leq-\frac{2c^2S}{810FC}<0.$$

No empirical choice of coefficient is needed. Optimizing just the scale
of $b_3-h/C$ also yields

$$\boxed{\sigma_F\leq\sigma_F^{(0)}
-\frac{2V^2}{810F\,(W-S/C)}<U_1.}$$

The denominator is positive: $WT>S^2$ and $T\leq CS$ imply $W>S/C$.
This replaces the difficult moment $T$ by a rigorous local bound and the
nonnegative elastic energy $S$. It still requires a quantitative estimate
or lower bound for $S/F$ to give a numerical improvement. A finite witness
alone proves $S>0$, but not $\liminf S/F>0$. The subsequent
[boundary-sensitive cylinder count](triangle-gauge-bulk-correction.md)
fills that gap with a positive local-event probability and a separate
volume-independent cost bound.

## Earlier unstratified estimates and their limitation

These exploratory estimates are retained to distinguish them from the
later [resolved fixed-coefficient evaluations](triangle-gauge-resolved-transport.md).

For the side-6 torus at $Q=F=72$, the exact charge-only trial gives
$U_1=0.022452355023951657$ and $c=0.02516796844742563$.
An exploratory coefficient $d=-2.9208964015355687\times10^{-5}$ in
$f=cb_3+dh$ was chosen from 128 independent exact stationary draws
(seed 125948001). This fitted coefficient is different from the
conservative analytic choice $d=-c/C$ above.

A first independent evaluation with 512 draws (seed 136059001) estimated
$\Delta U=-6.73\times10^{-5}$ with standard error
$2.57\times10^{-5}$. A larger fresh evaluation with 4096 draws
(seed 147160001), keeping the coefficient fixed, instead gave

$$\widehat{\Delta U}=-2.47\times10^{-5},\qquad
\operatorname{SE}(\widehat{\Delta U})=1.70\times10^{-5}.$$

The larger estimate is only about 1.46 standard errors below zero:
**the fitted trial's numerical improvement is not resolved**. In
particular, the initial pilot should not be reported as a replicated
effect. These standard errors are empirical, not rigorous coverage
guarantees.

Each draw sums every elastic slot for $S$, but samples only one uniform
full-bank slot for $T$. The larger evaluation estimated
$S/F\approx50\,965$ and $T/F\approx76\,106\,642$, with nonzero sampled
$T$ contributions in 556 of 4096 draws. The training estimate of $T/F$
was only about $44\,157\,420$. The noisy hidden Dirichlet cost was the
numerical bottleneck. The subsequent calculation conditions transition
sampling on the four changing-move families and uses a smaller coefficient
fixed before new draws. These are stationary variational
calculations, not imposed physical dynamics or evidence of a molecule.

## A sign constraint on adding elastic mixing

There is a separate exact comparison between the reaction/vacancy bank
and the current elastic extension. Keep the common attempted-update clock
and let

$$L_\lambda=L_0+\lambda E,\qquad 0\leq\lambda\leq1.$$

At $\lambda=0$, elastic proposal slots are no-ops; at $\lambda=1$, they
are the current bank. Intermediate values only describe an analytical
comparison of rates. No such rate change was added to the simulation.

Elastic slots carry no charge current, so the bare noise and source $b$
are unchanged. For each trial corrector, increasing $\lambda$ adds the
nonnegative cost $2\lambda\langle f,Ef\rangle$. The variational identity
therefore proves

$$\boxed{\sigma_F(1)\geq\sigma_F(0).}$$

Where the component kernel is fixed and the derivative is defined, the
physical-current corrector $\phi_\lambda=L_\lambda^+b$ gives

$$\frac{d\sigma_F}{d\lambda}
=\frac{\langle\phi_\lambda,E\phi_\lambda\rangle}{45F}\geq0.$$

Thus adding these reversible, charge-preserving elastic updates cannot
decrease the stationary winding mobility at matched primitive rates. This
is an application of variational monotonicity, not a new general transport
theorem. The witnesses above do not by themselves prove that this particular
comparison is **strict**, since that requires information about the optimal
corrector $\phi_\lambda$, not just $Lb_3$.

Calculation: [triangle_gauge_current_corrector.py](../tools/triangle_gauge_current_corrector.py).
Stationary moment estimates:
[triangle_gauge_transport_moments.py](../tools/triangle_gauge_transport_moments.py).
