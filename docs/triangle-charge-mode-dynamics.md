# Charge relaxation, transport, and a spectral constraint

The winding-current experiments suggest classical diffusive transport.
A complementary measurement follows the relaxation of nonconstant charge
modes under the same raw-link update bank. These modes are observables of
the supplied graph, not an imposed wave equation or quantum amplitudes.

## Intrinsic modes and an exact initial slope

Let $\phi_\alpha$ be orthonormal, nonconstant eigenvectors of the
honeycomb dual graph Laplacian $L_1$, with eigenvalue $\lambda$. Define

$$A_\alpha(X)=\sum_f\phi_\alpha(f)(q_f(X)-Q/F).$$

The stationary covariance is exactly
$\mathbb E[A_\alpha A_\beta]=\chi\delta_{\alpha\beta}$, where
$\chi=F\operatorname{Var}(q_f)/(F-1)$. For a band of degeneracy $d$,
measure

$$C_\lambda(\tau)=\frac1{d\chi}\sum_\alpha
\mathbb E[A_\alpha(X_{F\tau})A_\alpha(X_0)].$$

Each primitive's actual face-charge changes update these amplitudes.
Their endpoint values agree with projections of the final raw-link
charges. Correlations average time origins without subtracting a fitted
time mean; the exact ensemble mean is zero. Independent trajectories,
not time origins or degenerate modes, are the statistical clusters.

The previous exact charge projection gives

$$a=2(p_{110}+4p_{101}),\quad b=8p_{111},\qquad
\gamma_0(\lambda)=\frac{(a+10b)\lambda-b\lambda^2}{45\chi}.$$

Here $\gamma_0$ is the initial rate on the attempts-per-face clock for
the current 15-rule bank. Elastic moves preserve every charge, so they
do not change the numerator. The denominator is 45, not the old
13-rule bank's 39. Exactly one attempted update has

$$C_\lambda(1/F)=1-\gamma_0(\lambda)/F.$$

The expression $[1-\gamma_0/F]^{F\tau}$ extrapolates a memoryless
charge projection. It is not generally the true correlation curve.

## Memory slows this autocorrelation: an exact lower envelope

The complete proposal operator $P$ is self-adjoint in the stationary
reference. Its eigenvalues $p_j$ lie in $[-1,1]$. A normalized
autocorrelation has a nonnegative spectral measure:

$$C_\lambda(\tau)=\sum_j w_jp_j^{F\tau},\qquad
w_j\geq0,\quad\sum_jw_j=1.$$

All these tori have even $F$. For integer $\tau\geq1$, convexity of
$x\mapsto x^{F\tau}$ and the exact one-attempt slope imply

$$\boxed{C_\lambda(\tau)\geq
[1-\gamma_0(\lambda)/F]^{F\tau}.}$$

Thus the memoryless projection is a lower envelope of this correlation,
not a closure to assume in the evolution. Sampling estimates can fall
below it through noise; that does not reverse the exact inequality.

More strongly, setting $z_j=p_j^F\in[0,1]$ gives
$C_\lambda(\tau)=\sum_jw_jz_j^\tau$. Hence at these integer times it
is nonnegative, nonincreasing, and discretely completely monotone:

$$(-1)^k\Delta^k C_\lambda(\tau)
=\sum_jw_jz_j^\tau(1-z_j)^k\geq0.$$

It is also log-convex wherever positive. Effective logarithmic decay
rates over fixed-width intervals therefore cannot increase as the
interval moves later, for the exact correlation.

### Why this matters for the molecule goal

This reversible stochastic description cannot supply an undamped
oscillatory equilibrium charge mode merely through longer simulation or
a more elaborate charge-memory closure. Its exact spectral contributions
are decays or constants on this clock, not oscillatory wave poles.
Appropriate pointwise scaling limits retain the complete-monotonicity
constraint; averaging more such modes does not create coherence.

This is an application of the spectral theorem to our chosen dynamics,
not a new general no-go theorem. It does not rule out quantum behavior in
other Wolfram-model constructions, nonreversible or deterministic
schedules, or a physically justified extension of the state description.
It does prevent interpreting Fourier analysis of this model, or a noisy
wiggle in a relaxation trace, as a derived quantum amplitude.
Multiplying the generator by $i$ or analytically continuing the clock
would introduce an additional physical interpretation; it would not by
itself derive real-time quantum dynamics from these trajectories.

## The conditional bridge to winding mobility

For cell wavevector $k$, the lowest graph band has

$$\lambda_-(k)=3-|1+e^{ik_x}+e^{-ik_y}|
=\frac{k_x^2+k_y^2+k_xk_y}{3}+O(|k|^4).$$

If a diffusive fluctuation limit exists with equilibrium susceptibility
$\chi$ and isotropic winding-mobility tensor in these cell coordinates,
then its relaxation rate should satisfy

$$\gamma(k)\sim\frac{\sigma}{\chi}
(k_x^2+k_y^2+k_xk_y).$$

This is a **conditional comparison**, not a derived hydrodynamic limit.
The microscopic initial slope is consistent with it: as $k\to0$,
$\gamma_0\sim\sigma^{(0)}(k_x^2+k_y^2+k_xk_y)/\chi$ with the exact
bare estimate $\sigma^{(0)}=(a+10b)/135$.

For the lowest band on side $s$, set $K_s=(2\pi/s)^2$. Report the
finite-wavevector proxy $\chi\gamma_{\rm eff}/K_s$, and separately
$\chi\gamma_{\rm eff}/(3\lambda_-)$ to expose the known graph-dispersion
correction. Neither is an exact finite-volume mobility identity. Even
the initial-slope value of the latter is
$\sigma^{(0)}-b\lambda_-/135$, not $\sigma^{(0)}$.

At fixed nonzero wavevector, the integrated gradient current is a bounded
charge-mode difference; its zero-frequency variance rate vanishes.
Winding current is different: it carries a nontrivial torus circulation.
The comparison therefore needs a controlled long-wavelength/diffusive
limit, not an interchange of these distinct fixed-volume observables.

## Measurements at three wavelengths

Each size uses 128 independent exact stationary starts and independent
uniform schedules. Sides 6 and 8 run for 4096 attempts per face, and side
12 for 2048. Seed bases are 218984001, 230095001, and 241106001; initial
and schedule seeds are base $+2i$ and base $+2i+1$. The longest run at
side 12 stays within the existing engine's proposal limit.

The lowest nonconstant band has degeneracy six at all three sizes. Use
the fixed-interval effective rate

$$\gamma_{\rm eff}=-\frac{\log[C(t_2)/C(t_1)]}{t_2-t_1}.$$

Intervals were specified before collecting the data: $(8,32)$, $(16,64)$,
and $(32,128)$ attempts per face at sides 6, 8, and 12. These interval
slopes are diagnostics, not a claim that the entire curve is exponential.
The delta-method error retains covariance between the two endpoint
correlations using the independent trajectory clusters.

| Side; faces | $\gamma_0$ | Measured $\gamma_{\rm eff}$ | $s^2\gamma_{\rm eff}$ |
| --- | --- | --- | --- |
| 6; 72 | $0.050215$ | $0.045927\pm0.000401$ | $1.6534\pm0.0144$ |
| 8; 128 | $0.028838$ | $0.026582\pm0.000357$ | $1.7013\pm0.0228$ |
| 12; 288 | $0.013007$ | $0.012185\pm0.000302$ | $1.7547\pm0.0435$ |

The approximately inverse-square rate scaling is consistent with
diffusion. Its remaining size dependence is visible; this is not an
exact $s^{-2}$ law or a measured universal exponent. Both known lattice
dispersion and dynamical finite-wavelength corrections matter.

Relative to the exact initial rate, the later effective slopes are
slower by $8.54\pm0.80\%$, $7.82\pm1.24\%$, and $6.32\pm2.32\%$.
This resolves memory retardation on the shorter two lattices and is
consistent with it on the longest wavelength. It does not separate
nonlinear charge memory from hidden-gauge memory.

For example, at 72 faces and $\tau=32$, the measured correlation is
$0.22785\pm0.00294$, compared with the exact lower envelope $0.20040$.
At 128 faces and $\tau=64$, it is $0.17903\pm0.00395$, compared with
$0.15789$. The correlations therefore retain substantially more memory
than the unclosed charge projection predicts.

## Comparison with the independent winding measurement

The resulting mobility proxies are:

| Faces | $\chi\gamma_{\rm eff}/K_s$ | $\chi\gamma_{\rm eff}/(3\lambda_-)$ |
| --- | --- | --- |
| 72 | $0.020457\pm0.000179$ | $0.021109\pm0.000184$ |
| 128 | $0.020990\pm0.000282$ | $0.021358\pm0.000286$ |
| 288 | $0.021605\pm0.000536$ | $0.021771\pm0.000540$ |

The two columns expose the finite graph-dispersion correction; choosing
one does not eliminate the unknown dynamical correction. Their approach
toward the direct winding scale near $0.022$ is a useful cross-observable
consistency result. In particular, the longest-wavelength proxy is
consistent with that scale at the current precision.

The [winding experiments](triangle-winding-fluctuations.md) used separate
trajectories on 72 and 128 faces, and their finite-time biases remain
unquantified. They did not measure the 288-face thermodynamic mobility.
Thus these results support, but do not establish, a common limiting
transport coefficient or the Einstein relation in this model.

All $\pm$ values above denote one empirical standard error, not a
certified confidence interval. Exact normalization uses $C(0)=1$ in
expectation; sample zero-lag values were not forced to one. Small negative
late correlation estimates are compatible with sampling noise, not
evidence against the exact positive-spectrum constraint.

Calculation: [triangle_charge_mode_dynamics.py](../tools/triangle_charge_mode_dynamics.py).
Lowest-band numerical record:
[triangle-charge-mode-dynamics.json](../data/triangle-charge-mode-dynamics.json).
