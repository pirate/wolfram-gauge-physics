# Fiber refinement produces exclusion diffusion in the diffuse ensemble

For the present cycle-fiber rule bank, a generic finely populated fiber
does not retain finite-time charge reactions. Its charge process converges
to symmetric simple exclusion on the supplied dual mesh, uniformly in the
angular rate. This gives both a derived macroscopic diffusion law and a
specific obstruction to obtaining binding simply by increasing resolution.

The result concerns this selected classical stochastic model. It is not
a derivation of quantum matter, electromagnetic charge, or spacetime.
Exclusion-to-heat convergence is established mathematics; see
[Fathi and Simon](https://arxiv.org/abs/1507.06489).
The contribution here is the coupling from our actual shared-link rules,
its error bounds, and the resulting refinement constraint. No priority
claim about a new result in the literature is made.

## The microscopic process and its clock

Take an odd cycle fiber $C_n$. Every cycle automorphism is fixed by the
image of one vertex and the choice of which neighbor follows it. Thus its
actual automorphism group has the exact coordinates

$$U_e(v)=s_ev+a_e\pmod n,\qquad s_e\in\{-1,1\},\quad a_e\in\mathbb Z_n.$$

This is a faithful representation of the derived $D_n$ permutations, not
an imposed continuous gauge group. Shared links, actual based holonomies,
and the old two-spoke reconstruction are retained.

There are $F=2L^2$ faces and $3F$ rooted three-face fans. Each fan has the
vacancy exchange, the two Hurwitz orientations $H,H^{-1}$, and the twelve
[relational reactions](cycle-relational-dynamics.md), each at rate one.
Each angular channel has rate $\alpha\ge0$. Idle proposals do not advance
a separate normalized clock. All times below use these per-channel rates.

The derived conserved face weight is $q(E)=0$, $q(R)=1$, $q(Z)=2$.
Define $\eta_f=1$ when the face holonomy is even, including identity. Then

$$q_f=1+\eta_f-2\mathbf1_{\{E_f\}}.$$

Every $H$ or $H^{-1}$ exchanges the two face parity bits, regardless of
their angles. Angular moves preserve those bits and cannot create an
identity face. There are two rooted descriptions of each dual edge, each
with both Hurwitz orientations. Consequently their combined parity
exchange clock has rate **four per undirected dual edge**.

Couple an auxiliary exclusion field $\xi$ to the same Hurwitz clocks,
starting at $\xi_0=\eta_0$. This auxiliary process never drives the links;
it is an observable used to identify their limiting law. If there is no
initial identity and no reaction up to time $T$, then

$$q_t=1+\xi_t\quad\text{for every }0\le t\le T.$$

Indeed, Hurwitz moves only permute face conjugacy classes, vacancy moves
are idle in the absence of identities, and angular moves preserve $q$.
Reactions are the only remaining way to create an identity.

## Stationary path coupling, uniform in angular speed

Independent uniform raw links form a stationary measure for every
$\alpha$: reaction and angular channels are permutations, and the two
Hurwitz permutations have inverse partners. Under this reference,

$$p_E=\frac1{2n},\qquad p_R=\frac12,\qquad p_Z=\frac{n-1}{2n}.$$

The three based fan holonomies are independent uniform group elements.
One can expose their three distinct rim links last; each independently
sets one holonomy by a bijection. Each reaction moves $4n(n-1)$ out of
$8n^3$ tuples. The expected total changing-reaction intensity is therefore

$$\lambda_{\rm reaction}
=3F\times12\times\frac{4n(n-1)}{8n^3}
=\frac{18F(n-1)}{n^2}.$$

The union bound and the expected reaction count give a full-path coupling:

$$\boxed{
\Pr\!\left(\exists t\le T:q_t\ne1+\xi_t\right)
\le\min\!\left\{1,\frac{F}{2n}+\frac{18FT(n-1)}{n^2}\right\}.}
\tag{1}$$

This is also an upper bound on total variation between their path laws.
It is independent of $\alpha$, so it holds along any sequence of finite
rates $\alpha_n$, including rates growing arbitrarily fast with $n$.
No interchange with an unconstructed infinite-rate process is needed.

The initial exclusion law is uniform over parity configurations with
even total particle count. To see this, raw link signs are independent
fair bits, and face signs are their incidence image over $\mathbb F_2$.
The connected dual graph has incidence rank $F-1$; its image consists of
even-sum face signs. Since $F$ is even, the same constraint holds for
$\eta$. This is a mixture of conserved-number exclusion equilibria, not
a single prescribed-density sector.

There is also a volume-free fixed-time density bound. Common Hurwitz
swaps preserve $\|\eta-\xi\|_1$. Each reaction changes two actual parity
bits. A vacancy event changes parity only on an $E/R$ pair, whose
stationary probability per root is $1/(2n)$. Thus

$$\frac{\mathbb E\|\eta_T-\xi_T\|_1}{F}
\le\left(\frac{39}{n}-\frac{36}{n^2}\right)T.$$

Using $\mathbb E[2N_E/F]=1/n$ gives

$$\boxed{\frac{\mathbb E\|q_T-(1+\xi_T)\|_1}{F}
\le\frac1n+\left(\frac{39}{n}-\frac{36}{n^2}\right)T.}
\tag{2}$$

The apparatus counts all changing vacancy events, including $E/Z$
exchanges. Only $E/R$ events enter the parity estimate above.

## Nonuniform density profiles: a stopped reference-process argument

Uniform equilibrium alone cannot demonstrate relaxation of a prepared
density profile. A broader class of initial states admits a similar bound:
choose any distribution of raw signs, then choose the shifts independently
uniform conditional on those signs.

Consider a reference link process containing only $H,H^{-1}$ and angular
moves. It preserves the family of distributions with conditional-uniform
shifts. Hurwitz moves map one sign fiber bijectively to another, with a
sign update independent of the shifts. Angular moves preserve each sign
fiber and act by reversible shift permutations. This remains true for an
arbitrary initial distribution over sign fibers, not just equilibrium.

Conditional on signs, the three fan shifts are independent uniform. For
an $RRR$ fan, exactly one adjacent equality has probability
$2(n-1)/n^2$, giving changing-reaction rate $24(n-1)/n^2$. For a fan with
one reflection and two even holonomies, an $ERZ$ arrangement has probability
$2(n-1)/n^2$, giving rate $8(n-1)/n^2$. Other parity patterns cannot react.
Hence the reference expectation of reaction intensity is at most
$72F(n-1)/n^2$.

Couple the full link process to that reference until the first reaction.
If no initial identity is present, vacancy proposals are idle before this
time. Bounding the reaction clocks evaluated on the reference trajectory
therefore proves

$$\boxed{\Pr\!\left(\exists t\le T:q_t\ne1+\xi_t\right)
\le\min\!\left\{1,\frac Fn+\frac{72FT(n-1)}{n^2}\right\}.}
\tag{3}$$

Crucially, this proof does **not** assume that the full reacting process
preserves conditional-uniform shifts after a reaction. It uses that
property only for the reference trajectory and a stopped coupling.

## The macroscopic equation, derived without fitting

Write the two dual sublattices as $A(x,y)$ and $B(x,y)$. The three
neighbors of $A(x,y)$ are $B(x,y)$, $B(x,y-1)$, and $B(x+1,y)$.
The positive dual graph Laplacian has Bloch branches

$$\lambda_\pm(k)=3\pm|\phi(k)|,\qquad
\phi(k)=1+e^{-ik_y}+e^{ik_x}.$$

The exclusion mean closes exactly:

$$\partial_t\mathbb E\xi=-4\mathcal L_{\rm dual}\mathbb E\xi.$$

No independence approximation is needed: exchanging two occupation bits
changes either one by the difference of the two bits. Expanding the
acoustic branch,

$$\lambda_-(k)
=\frac{k_x^2+k_y^2+k_xk_y}{3}+O(|k|^4).$$

With $X=x/L$, $Y=y/L$ and microscopic time $t=L^2\tau$, the limiting
density satisfies

$$\boxed{\partial_\tau\rho=\nabla\cdot D\nabla\rho,\qquad
D=\begin{pmatrix}4/3&2/3\\2/3&4/3\end{pmatrix}.}$$

These are mesh cell coordinates. If the supplied triangular mesh is
embedded with unit primal edges and basis
$a_1=(1,0)$, $a_2=(-1/2,\sqrt3/2)$, then
$[a_1\ a_2]D[a_1\ a_2]^T=I$: the same equation is isotropic with
diffusion coefficient one in that embedding and the specified clock.
The two-dimensional geometry is supplied, not emergent.

This identification goes beyond the mean. For a real acoustic mode $h$
with $\sum h_f^2=F/2$, set $Z=F^{-1}\sum h_f\xi_f$. Its exact drift is
$-4\lambda_- Z$. Its martingale quadratic-variation rate is bounded by

$$\frac4{F^2}\sum_{\{i,j\}}(h_i-h_j)^2(\xi_i-\xi_j)^2
\le\frac{2\lambda_-}{F}.$$

Prepare independent Bernoulli occupations with smooth probabilities
bounded away from zero and one, conditioned on even total count. If
$p_{\rm even}$ is the probability of that conditioning event, the initial
mode variance is at most $1/(8Fp_{\rm even})$. The conditioning changes
one-site means by an exponentially small amount, since

$$\mathbb E[\xi_i\mid\text{even}]
=p_i\frac{1-\prod_{j\ne i}(1-2p_j)}{1+\prod_j(1-2p_j)}.$$

Solving the linear martingale equation gives the finite-volume bound

$$\operatorname{Var}Z_t
\le\frac{e^{-8\lambda_-t}}{8Fp_{\rm even}}
+\frac{1-e^{-8\lambda_-t}}{4F}.$$

Thus each fixed macroscopic Fourier mode concentrates as $F\to\infty$.
The exact Bloch mean and approximation of smooth functions by finitely
many Fourier modes give weak empirical-density convergence to the heat
equation at fixed macroscopic times. The optical branch decays on the
microscopic scale and does not supply a second conserved density.

Equation (3) transfers this conclusion to the actual raw-link charge
field under the sufficient scaling

$$n/L^4\longrightarrow\infty,$$

uniformly over angular rates. Indeed, $F=2L^2$ and $T=L^2\tau$ make its
error $O(L^4/n)$. This is a sufficient condition for coupling the entire
charge path, not a necessary resolution threshold. The normalized charge
density tends to $1+\rho$. Under the stationary reference alone, (2) gives
the weaker sufficient condition $n/L^2\to\infty$ for fixed-time normalized
$L^1$ approximation; that statement does not establish relaxation of a
nonuniform profile.

## Raw-link observations

The stationary experiment uses 4,096 independent trajectories per
$(n,\alpha)$, $F=32$, $T=0.02$, $n=9,27,81,243,729$, and
$\alpha=0,1,8$. All fifteen saved runs completed. At $\alpha=1$, the
observed probability of any charge-path mismatch falls from
$0.9072\pm0.0045$ at $n=9$ to $0.03027\pm0.00268$ at $n=729$.
The latter has rigorous upper bound $0.03773$. At $n=729$, the mean
reaction count is $0.01587\pm0.00224$, versus exact expectation $0.01578$.
Errors here and below are one independent-trajectory standard error.

For the nonuniform experiment, use the real lowest acoustic mode at
$k=(2\pi/L,0)$ and initial probability $p_f=1/2+0.2h_f$, conditioned on
even count. A dual spanning tree solves the raw-sign incidence equation;
non-tree signs remain independent fair bits. All raw shifts are uniform.
This samples the specified sign law without imposing future motion.

The fiber has $n=1,000,000,007$ vertices, represented implicitly by its
exact affine automorphisms. The simulation stores actual link group
elements, not a billion-vertex expanded copy of every fiber. All old
reactions remain enabled; none was artificially removed. Summing the
angular involutions gives precisely the two allowed $z\to z\pm1$
transitions, excluding zero, so idle slots can be thinned without changing
the generator or the physical clock.

At $\tau=0.02$ and $\alpha=1$:

| Side / faces | Trajectories | Predicted final mode | Observed raw charge mode |
| --- | ---: | ---: | ---: |
| 4 / 32 | 2,048 | 0.037613 | $0.035619\pm0.001394$ |
| 8 / 128 | 1,024 | 0.035537 | $0.035175\pm0.000964$ |
| 16 / 512 | 512 | 0.035055 | $0.034470\pm0.000724$ |

The initial expected mode is $0.1$. The scaled microscopic decay rates
are $48.89165$, $51.72924$, and $52.41201$, tending to the independently
derived continuum value $16\pi^2/3=52.63789$. These are spectral predictions,
not fitted exponents. At side 8 and $\alpha=8$, another 1,024 trajectories
give $0.035121\pm0.001006$ against the same prediction.

All 4,608 nonuniform trajectories had no initial identities, no changing
reactions or vacancy events, and exact charge-path agreement with their
coupled exclusion processes. This is consistent with the small bounds;
it does not imply those event probabilities are zero. The largest
per-trajectory bound is $1.893\times10^{-4}$ on the 512-face mesh.
Final face charges are independently reread from the evolved raw links.
Representative initial and final link vectors are saved for side 8.

## Consequence for the molecule question

Increasing angular speed cannot offset the disappearing reaction measure
in these ensembles: bounds (1) and (3) do not depend on that speed. The
limiting fixed-particle-number exclusion equilibrium is uniform over
placements, with no energetic preference for neighboring particles.
Its off-site covariance is only the number-constraint term
$-\rho(1-\rho)/(F-1)$. This limit supplies diffusion and exclusion, not
an attractive binding mechanism.

This is not a no-go theorem for every initial state or every time scale.
The [coherently embedded coarse states](cycle-relational-dynamics.md) and
[aligned relative-angle sectors](fiber-relative-angle.md) violate the
diffuse-shift assumption and can retain reactions. Also, an $O(1/n)$
reaction rate per face need not disappear on times of order $n$.
Neither long-time kinetics in those regimes nor an invariant dynamically
selected population of aligned states has been derived here.

The next interaction question is therefore specific: can the existing
rules create or preserve relative-holonomy alignment at a nonvanishing
rate as $n$ grows, without preparing the desired alignment by hand?
Simply rescaling reaction rates, adding attraction, or interpreting the
diffusion equation as quantum evolution would not answer it. Nor does
the finite-time result justify interchanging refinement and infinite-time
limits in the unresolved DC-mobility problem.

The [conditional alignment-boundary calculation](fiber-alignment-boundary.md)
now resolves one possible longer-time mechanism exactly: a fixed
relative-angle reference produces a partially absorbing continuum boundary
on times of order $n$. Whether the full dynamics preserve such references
for that long remains unresolved; the conditional construction does not
freeze them in the real simulation.

The [full-dynamics reaction-burst calculation](fiber-reaction-bursts.md)
now establishes another mechanism without fixed anchors: a typical rare
reaction leaves local loop-flatness constraints that support further,
spatially overlapping reactions. Their instantaneous total rate tends to
25 at fixed volume even as the ordinary stationary rate vanishes.

Apparatus and observations:
[stationary coupling](../tools/fiber_refinement_limit.py),
[stationary data](../data/fiber-refinement-limit.json),
[raw-link profile evolution](../tools/fiber_refinement_heat.py),
[profile data](../data/fiber-refinement-heat.json).
