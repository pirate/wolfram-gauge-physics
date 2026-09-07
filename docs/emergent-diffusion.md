# From finite gauge updates to deterministic diffusion

We can now explain a measured transport phenomenon through an exact factor map, rather than
recognizing a suggestive shape. A charge projection of a derived square-fiber interaction
is exactly a colored exclusion process. Its microscopic deterministic evolution displays
ballistic motion, diffusive tagged-particle spreading, or blocking depending on the initial
background density. The geometry and alternating update schedule are still supplied.

This is a bridge to **known statistical physics**, not a claim to have discovered diffusion,
electromagnetism, quantum mechanics, or emergent spacetime. The vacancy-plus-two-color restriction
matches equation (1) of [Medenjak, Klobas and Prosen, *Diffusion in deterministic interacting
lattice systems*](https://arxiv.org/abs/1705.04636). Their full time step contains two of our
matching layers. Their collective charge-transport calculation and our tagged-particle
calculation are different observables; particle-order conservation supplies the latter below.

![Microscopic spreading, transport estimates, and exact displacement distribution](images/derived-diffusion.png)

## An exact observable map, not a fitted replacement

Let $F:G^2\to G^2$ be a table on $G=\operatorname{Aut}(C_4)$ and let $\pi$ send an element
to the vector of its independently derived additive class charges. We test all 64 inputs for
an autonomous map $f$ satisfying

$$ (\pi\times\pi)\circ F=f\circ(\pi\times\pi). $$

This means all microscopic representatives of the same incoming charge pair give the same
outgoing charge pair. It implies the identity after every composition of supported local
updates, not only on sampled trajectories. Conjugating a cell through its connector does not
change a class charge, so the statement also holds with nonidentity static connectors.

All fourteen minimum-support subgroup-growing strict braid rules pass. Each projects to

$$f(a,b)=\begin{cases}(b,a),&a=v\text{ or }b=v,\\(a,b),&\text{otherwise},\end{cases}$$

for some distinguished symbol $v$. Six designate the zero-charge symbol as $v$; the others
designate one of two charged symbols. This is a discovered label-exchange law, not a force
inserted into the microscopic table. Here zero charge merges the identity and central half-turn;
it does **not** mean microscopic flatness. Full one-cell conjugacy sectors also form a closed
factor. The [all-rule classification](braid-quotient-classification.md) extends the analysis
beyond these fourteen laws and identifies its precise limitations.

## Exact finite-time tagged-particle distribution

Use a zero-vacancy rule (table 11229 in the checked-in census), one unique tag at the origin,
and background occupancy $n_i\in\{0,1\}$ independently Bernoulli with density $\rho$.
Other occupied sites have a different conserved color. The tag occupies $n_0=1$ by construction.

Two facts follow directly from the effective rule:

1. Occupancy bits evolve as unconditional adjacent swaps: swapping equal bits changes nothing.
   After an even number $t$ of layers, the bit initially at $i$ is at
   $r_t(i)=i+(-1)^{i+\phi}t$, where $\phi\in\{0,1\}$ chooses the first matching.
2. Occupied colors cannot pass each other. Therefore the tag retains its rank in occupied-site
   order even though the occupancy bits themselves follow those free trajectories.

The rank condition gives a finite sum despite the infinite reference chain:

$$\Pr[X_t\le x\mid\phi]
=\Pr\!\left[\sum_i n_i\big(\mathbf1_{r_t(i)\le x}-\mathbf1_{i\le0}\big)\ge0\right].$$

Let $m_+$ and $m_-$ count non-origin coefficients $+1$ and $-1$, and let $c_0$ be the origin's
fixed coefficient. Then this is exactly

$$\Pr[U-V+c_0\ge0],\qquad
U\sim\operatorname{Bin}(m_+,\rho),\quad V\sim\operatorname{Bin}(m_-,\rho),$$

with independent $U,V$. Average the two phases and difference adjacent CDF values to obtain
the PMF. `tools/exact_tag_distribution.py` evaluates this finite sum with floating-point
binomial probabilities; the formula is exact, the numeric evaluation is not rational arithmetic.
Tests compare it against exhaustive two-layer initial states at three densities.

At long times and fixed $0<\rho<1$, occupancy flux has variance
$t\rho(1-\rho)+O(1)$. Rank preservation converts flux to tag displacement by density $\rho$:
the leading fluctuation is the flux divided by $\rho$. The Bernoulli central-limit scaling thus
gives the asymptotic prediction

$$\operatorname{Var}X_t\sim\frac{1-\rho}{\rho}t=2D_{\rm layer}t,
\qquad D_{\rm layer}=\frac{1-\rho}{2\rho}.$$

The exact finite-time CDF, rather than this asymptotic argument alone, is used as the
distribution-level comparator. At $\rho=0$ the phase-balanced tag has $X_t=\pm t$, so variance
is $t^2$; at $\rho=1$ it is blocked. The divergence of $D$ as $\rho\to0$ does not replace the
ballistic endpoint with a diffusion process.

## Microscopic experiment and uncertainties

The C++ runner evolves raw group elements with `CellChain`, checks **every pair update** against
the proved charge factor, and verifies initial/final transported total holonomy. Background
group representatives are sampled uniformly within each charge label. Tag representatives and
matching phases are balanced deterministically; factor closure makes their correlation irrelevant
to the tag trajectory. No stochastic collision, fitted force, PDE update, or molecular target
is used. Static identity connectors select a convenient frame on the open chain.

The first run uses 1,040 cells and 256 layers; the independent-seed replication uses 2,064 cells
and 512 layers. Each has 4,096 trajectories per intermediate density and 256 per endpoint
control: 8,462,336,000 microscopic pair updates in total. The boundary guard $N\ge4t+16$ keeps
the union of possible tagged-particle past cones away from the ends, not merely the final tag.

| Background $\rho$ | Asymptotic $D_{\rm layer}$ | 256-layer estimate [95% interval] | 512-layer estimate [95% interval] |
|---|---:|---:|---:|
| 0.25 | 1.5 | 1.5292 [1.4490, 1.6150] | 1.4289 [1.3480, 1.5132] |
| 0.50 | 0.5 | 0.4950 [0.4763, 0.5154] | 0.4880 [0.4607, 0.5152] |
| 0.75 | 1/6 | 0.1721 [0.1636, 0.1805] | 0.1565 [0.1495, 0.1647] |

These are half-slopes of variance over the last three checkpoints, with 1,000 bootstrap
resamples of 32 independent trajectory blocks. Time correlations are retained within each
block. The intervals are pointwise, exploratory intervals, not simultaneous acceptance gates.
**The high-density replication interval misses the prediction.** We retain it and do not
resample until a preferred result appears. Overall agreement is supportive, not perfect.

All six final empirical CDFs differ from the exact CDF by at most 0.0166. This is a descriptive
error, not a fitted tolerance or formal goodness-of-fit test. The low-density first run's
log-variance exponent interval excludes the asymptotic value 1, but includes the exact
finite-time fit value 1.01287. Finite-window exponents should be compared with that prediction,
not assumed to equal their limiting values. Raw displacement histograms and joint-time block
moments are retained so these statistical choices can be independently changed.

## Reproduce

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build -L wgphysics --output-on-failure
python3 tests/effective_dynamics.py
python3 tools/charge_quotient.py --output out/charge-quotient.json
diff -u data/d4-charge-quotient.json out/charge-quotient.json
# Full microscopic runs (billions of updates; not run on every CI build):
python3 tools/run_diffusion_experiment.py --output out/diffusion.json
python3 tools/run_diffusion_experiment.py --cells 2064 --layers 512 --seed 1731291 --output out/diffusion-replication.json
# Recalculate statistics from recorded microscopic observations, without evolving again:
python3 tools/run_diffusion_experiment.py --reanalyze data/d4-diffusion.json --output out/diffusion-reanalyzed.json
uv run --with matplotlib python tools/plot_diffusion.py
```

The contribution is an independently checkable gauge-to-effective-dynamics bridge and its
microscopic validation. It is not yet a derivation of a physical gauge field, a calibrated
Hamiltonian, quantum amplitudes, three-dimensional space, or a bound state.
