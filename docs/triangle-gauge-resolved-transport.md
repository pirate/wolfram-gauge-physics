# A numerically resolved hidden-gauge transport correction

The explicit hidden-gauge function $h=ELb_3$ now yields a resolved
finite-volume reduction of the charge-drift trial functional on both
72- and 128-face tori. Its coefficient was fixed before these draws;
no physical dynamics or forces were fitted to obtain the result.

This complements the [strict bulk existence proof](triangle-gauge-bulk-correction.md).
The latter is rigorous but numerically tiny; the estimates here have useful
precision but empirical, not confidence-certified, error bars. Neither is
an evaluation of the full long-time mobility or a derivation of quantum
or molecular behavior.

## The quantity being estimated

With $L$ the complete positive generator, $E$ its elastic part, and
$b_3$ the integer winding-current drift, recall

$$h=ELb_3,\quad V=\langle b_3^2\rangle,\quad
W=\langle b_3,Lb_3\rangle,\quad
S=\langle Lb_3,ELb_3\rangle,\quad T=\langle h,Lh\rangle.$$

The exact conditional identity $\mathbb E[h\mid q]=0$ makes this a
genuinely hidden-gauge correction, not an additional charge-only fit.
Use the trial $f=cb_3+dh$, with $c=V/W$ obtained from exact charge counts.
Its change from the previous one-function functional is

$$\boxed{\Delta U=\frac{4cdS+2d^2T}{810F}.}$$

The same coefficient $d=-1.4604482007677843\times10^{-5}$ was used
at both sizes. It is half the coefficient from the old 128-draw pilot,
chosen **before** collecting either sample below. That pilot had
underestimated the expensive hidden energy $T$. Halving its coefficient
reduces the quadratic cost relative to the linear gain. This is a choice
of variational trial, not a rate in the update bank.

## Sum the easy energy, stratify the difficult one

For each independent exact stationary raw connection $X$, evaluate

$$s(X)=\frac12\sum_{a\in E}[Lb_3(T_aX)-Lb_3(X)]^2$$

over every elastic slot. Thus $\mathbb E[s]=S$, with no transition
sampling error in that part.

For $T$, separate changing primitive slots into vacancy, elastic,
reaction-birth, and reaction-death families. Let $K_r(X)$ be the number
of changing slots in family $r$. Draw one uniformly from each nonempty
family and use

$$\widehat t_r(X)=\frac{K_r(X)}2[h(T_{a_r}X)-h(X)]^2.$$

Then $\mathbb E[\widehat t_r\mid X]$ is exactly the family's Dirichlet
sum at $X$. Empty families contribute zero. Omitting no-op proposals
from this **estimator**, with their correct multiplicity weights, does
not change the attempted-update clock or the physical schedule.

Each transformed state's $h$ is calculated from the actual group
multiplications along based face paths and the reconstructed shared
link. Batching those products reproduces the scalar raw-connection
calculation; it does not replace link updates with charge transitions
or independently quotient overlapping patches.

The error estimate uses one whole stationary state as its independent
unit. In particular, the four sampled family energies from a shared
state are not treated as four independent observations.

## Two fresh, fixed-coefficient evaluations

Both references have $Q=F$, full based-loop image $S_3$, and reflection
faces. Each row uses 4096 independent exact reference draws.

| Faces | Seed | Estimated $\Delta U$ | Empirical standard error | Extra reduction relative to bare estimate |
| --- | --- | --- | --- | --- |
| 72 | 158329001 | $-5.53317\times10^{-5}$ | $1.03926\times10^{-6}$ | $0.23317\%$ |
| 128 | 169430001 | $-5.26945\times10^{-5}$ | $8.89899\times10^{-7}$ | $0.22387\%$ |

The corresponding trial-functional estimates are $0.0223970233$ and
$0.0222250410$. Their exact one-function baselines are $0.0224523550$
and $0.0222777355$, respectively. Sampling the functional of a valid
trial does **not** turn its sample mean into a certified numerical upper
bound: the exact expectation is the upper bound.

The effects are well separated from zero on the empirical error scale.
Two sizes support persistence of a useful correction but do not determine
its bulk value or justify extrapolating a precise thermodynamic number.
The old unstratified trial with twice this $d$ remains a separate,
inconclusive experiment; this result does not retroactively validate it.

The local energies are

$$\begin{array}{c|cc}
F&S/F&T/F\\\hline
72&51695.9\ \pm653.1&(73.1109\ \pm1.4366)\times10^6\\
128&49740.6\ \pm471.5&(71.5816\ \pm1.3831)\times10^6
\end{array}$$

where $\pm$ denotes one empirical standard error, not a rigorous interval.

## Which microscopic moves contribute to the hidden energy?

Each entry below is an estimate of $T_r/F$ in millions:

| Family | 72 faces | 128 faces |
| --- | --- | --- |
| Vacancy | $5.455\pm0.209$ | $4.852\pm0.187$ |
| Elastic | $21.770\pm0.853$ | $20.895\pm0.789$ |
| Birth | $23.438\pm0.821$ | $22.425\pm0.763$ |
| Death | $22.448\pm0.631$ | $23.409\pm0.713$ |

These are Dirichlet costs of this particular hidden function, not reaction
probabilities, energy levels, or separate contributions to a force.
Most of its variation occurs across the reaction and elastic channels.

Detailed balance gives an exact identity, applicable to any real function:

$$T_{\rm birth}=T_{\rm death}.$$

Indeed, $(X,a)\mapsto(T_aX,a^{-1})$ pairs each changing birth with a
death, preserves stationary weight, and preserves the squared increment.
The observed birth-minus-death differences are
$(0.990\pm1.015)\times10^6$ and $(-0.984\pm1.051)\times10^6$.
They are consistent with this independently derived identity. The
comparison can expose estimator or reconstruction errors; agreement
alone is not proof that every implementation detail is correct.

## An initial memory timescale, not a relaxation-time fit

The 128-face calculation also estimates
$H=\langle h^2\rangle$:

$$H/F=981210\pm26505.$$

For the actual discrete proposal operator $P=I-L/(45F)$,

$$\frac{\langle h,Ph\rangle}{H}=1-\frac{T}{45FH}.$$

Therefore the inverse initial normalized-correlation slope in attempted
updates per face is $45H/T\approx0.617$. This is a ratio of estimated
moments, not a fitted exponential decay time, a bound on slow tails,
or a physical spacetime unit. A relatively rapid initial hidden response
can still affect the long-time variational transport functional.

The subsequent [autonomous winding measurements](triangle-winding-fluctuations.md)
compare this analysis with actual finite-time currents at both volumes.
They resolve a variance crossover and short-window correlations, but
their long-window precision does not isolate this small additional
variational correction. Direct local bulk moments and a controlled
long-time residual-corrector estimate remain useful next steps. None of
these should be replaced by fitting a continuum law into the updates.

Calculation: [triangle_gauge_stratified_moments.py](../tools/triangle_gauge_stratified_moments.py).
