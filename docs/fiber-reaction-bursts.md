# Spontaneous reaction bursts carried by local flatness constraints

The unrestricted shared-link model retains a local interaction mechanism
that the [fixed-time diffuse limit](fiber-refinement-limit.md) does not
resolve: rare reactions occur in correlated bursts. A reaction exchanges
one local flatness constraint for another, leaving nearby reaction sites
active. No holonomies are held fixed in the dynamics studied here.

There are three model-specific results:

- The exact rate immediately after a marked reaction is 8 on the same
  fan and $1/3+16(n-1)/(3n^2)$ on a specified one-face-overlap neighbor.
- At fixed mesh volume, the total post-reaction rate tends to 25 as
  $n\to\infty$, despite the ordinary stationary rate tending to zero.
  The two reaction directions contribute limiting rates 20 and 30.
- Spontaneous stationary trajectories exhibit the predicted clustering
  and cross-region reactions. Conserved-charge sector mixing accounts for
  much, but not all, of their large count variance.

These are classical, reversible reaction-transport effects in this
selected rule bank. They do not establish molecular binding, quantum
amplitudes, or physical electric charge. Event conditioning and reversible
spectral identities are standard methods; the local rates, flatness
carriers, and bounds below are the calculation for this model. No claim of
priority over the literature is made.

## Observing a typical reaction, not inserting one

Use the actual $D_n$ cycle automorphisms on the supplied triangular torus,
with odd $n$, $F=2L^2$ faces, $L\ge4$, and the existing per-channel clock:
fifteen old channels per rooted fan at rate one, angular channels at rate
$\alpha$. Write $\mu$ for independent uniform raw links and $\lambda_p(x)$
for the total changing-reaction intensity on a marked rooted fan $p$.

$$\lambda_p=\begin{cases}
12,&RRR\text{ with exactly one adjacent equality},\\
4,&ERZ\text{ in any order},\\
0,&\text{otherwise}.
\end{cases}$$

Direct group counting gives

$$\bar\lambda_p=\frac{6(n-1)}{n^2},\qquad
\mathbb E_\mu[\lambda_p^2]=\frac{48(n-1)}{n^2}.$$

The law seen immediately after a typical marked stationary reaction is

$$\mu_p^+(x)=\frac{\mu(x)\lambda_p(x)}{\bar\lambda_p}.$$

Each reaction is an involution, so its incoming and outgoing stationary
weights agree; summing the incoming marked channels proves the formula.
It is the event-observer, or Palm, law, rather than the law at an arbitrary
time. For background on this distinction, see
[Coeurjolly, Møller, and Waagepetersen](https://arxiv.org/abs/1512.05871).

There is an exact sampler. With equal probabilities, choose a uniform
active $RRR$ tuple or a uniform $ERZ$ tuple. Draw all other raw coordinates
uniformly, then solve the three distinct rim links for that tuple. Each
rim link is a bijective coordinate for one of the actual based holonomies.
This samples the states a stationary event observer encounters; it does
not supply a new interaction or select a desired future trajectory.

All links are released immediately and evolve with the original rules.
Separate long runs start with completely unconditioned uniform links and
observe spontaneous reactions directly.

## Exact local and neighboring rates

The marked post-event rate is exactly

$$\boxed{\mathbb E_{\mu_p^+}\lambda_p
=\frac{\mathbb E\lambda_p^2}{\bar\lambda_p}=8.}$$

Half the post-event tuples are $RRR$ with intensity 12; half are $ERZ$
with intensity 4. The ordinary rate is asymptotic to $6/n$, so the
post-event enhancement grows as $4n/3$.

Choose a fan $q$ sharing exactly one face with $p$, with each of its two
exterior faces having a private rim edge outside the marked fan and the
other exterior face. This is an actual read-set property, not an
independent-face evolution assumption. It makes the two exterior based
holonomies independent uniform variables conditional on the marked data.
There are twelve such neighbors, four for each marked face.

Under $\mu_p^+$, the shared face has class probabilities

$$\Pr(E)=\frac16,\qquad\Pr(R)=\frac23,\qquad\Pr(Z)=\frac16.$$

The conditional expected target intensity for each shared class is

$$\mathbb E[\lambda_q\mid E]=\frac{2(n-1)}n,\quad
\mathbb E[\lambda_q\mid R]=\frac{8(n-1)}{n^2},\quad
\mathbb E[\lambda_q\mid Z]=\frac2n.$$

For example, with a shared identity, either ordering of independent
$R/Z$ exterior faces activates four channels. With a shared reflection,
the $ERZ$ contribution is $2(n-1)/n^2$ and the active-$RRR$ contribution
is $6(n-1)/n^2$. Therefore

$$\boxed{\mathbb E_{\mu_p^+}\lambda_q
=\frac13+\frac{16(n-1)}{3n^2}\longrightarrow\frac13.}$$

This is more than immediate reversal on the original fan. Suppose a
forward $RRR\to ERZ$ event puts its new identity on the shared face and
the target's exterior classes are $R/Z$. Before the parent event the
target was $RRZ$, hence inactive; afterwards it is $ERZ$, hence active.
The exterior face classes are unchanged by the parent rewrite.

The probability of this newly opened target gate under the marked Palm
law is $(n-1)/(12n)$. Across the twelve specified neighbors the expected
number of newly opened gates is exactly $(n-1)/n$, tending to one per
typical reaction, or two conditional on a forward reaction.

These gates compete for the **same** new identity face. Treating them as
independent offspring would give a false branching-process picture.
If a neighbor consumes that identity and an exterior $Z$, while the
parent's new $Z$ remains, the two reactions together move one unit of
conserved $q$ from the exterior $Z$ face to the parent's $Z$ face. This is
reaction-mediated charge transport, not a fitted force.

The saved $C_5$ witness realizes this on actual raw links, with parent
faces $(6,1,0)$ and target faces $(0,3,2)$:

| Face | Before parent reaction | Between reactions | After target reaction |
| --- | ---: | ---: | ---: |
| 3, exterior source | 2 | 2 | 1 |
| 0, shared face | 1 | 0 | 1 |
| 6, parent destination | 1 | 2 | 2 |

Every other face returns the same charge it had initially, and total
$Q=52$ is unchanged. The source lies outside the parent's support. The
target is inactive before the first reaction and active after it. The
three raw link vectors and native based holonomies are retained in the
carrier data; this is an allowed two-event witness, not a forced trajectory
used to estimate spontaneous probabilities.

## Two local flatness carriers explain the number 25

The adjacent equality $a=b$ of two based reflections is equivalent to
$ab=E$, since a reflection is its own inverse. Their internal common
spoke cancels in the product, leaving identity holonomy around the
four-edge boundary of the two faces. This is gauge-invariant flatness of
a two-face region, not equality of unrelated local coordinate labels.

The reaction $RRR\leftrightarrow ERZ$ interconverts:

1. A flat two-face boundary whose two faces are reflections.
2. A flat triangular face, adjacent in the reacting fan to $R$ and $Z$.

In the diffuse event-conditioned ensemble, either is one independent
linear equality in the raw affine shifts. Additional distinct local
flatness equalities cost another factor of order $1/n$. Distinct triangle
or two-face boundaries have different edge supports, so their shift
equations have independent unit-coefficient pivots. A change of basepoint
can reverse the sign of a rotational holonomy, but does not change its
zero set.

At fixed $F$, this identifies every leading post-event reaction site:

- **After a forward event:** the identity face belongs to nine rooted
  three-face fans. The original is certainly $ERZ$. In each of the other
  eight, the remaining parities are $R/Z$ with probability $1/2$ at leading
  order. Their total expected reaction rate is $4(1+8/2)=20$.
- **After a backward event:** the flat adjacent reflection pair belongs
  to four rooted three-face fans. The original has a distinct third
  reflection; each of the other three has a reflection as its third face
  with probability $1/2$. Their total expected rate is $12(1+3/2)=30$.

Unrelated sites contribute $O(F/n)$. The two directions have equal
stationary event flux, so

$$\boxed{\lim_{n\to\infty}\mathbb E_{\mu_p^+}\lambda_{\rm total}
=\tfrac12(20+30)=25,\qquad F\text{ fixed}.}$$

This is a statement about an event-conditioned neighborhood, not an
order-one unconditioned reaction density. The latter remains
$18F(n-1)/n^2$. Nor can the fixed-volume limit be used unchanged when
$F/n$ is appreciable and unrelated background events contribute.

At $n=1,000,000,007$, 65,536 independent Palm samples on 32 faces give
$25.0010\pm0.0379$ for the total rate. Separating directions gives
$19.9504\pm0.0300$ and $30.0680\pm0.0574$. Another 32,768 samples on
128 faces give $25.0436\pm0.0538$. These uncertainties are one independent
sample standard error; none of these constants was fitted to the samples.

## A full-dynamics recurrence bound

A fan reads seven edges. Exactly nineteen rooted fans can write at least
one of them: each edge has four writing roots, and nine roots write two
edges in this read set, giving $7\times4-9=19$.

The angular generator can be represented by its two directional
proposals per root, each at rate $\alpha$, with forbidden moves idle.
Thus the total proposal intensity of all potentially interfering roots is

$$C=19(15+2\alpha),$$

independent of $n$ and mesh volume. Consider the event that the first such
proposal is a changing marked reaction and occurs before $\Delta$.
No other interfering proposal has yet altered the marked input. Distant
updates are allowed throughout. Averaging the marked rate gives

$$\Pr_{\mu_p^+}(\text{another marked reaction by }\Delta)
\ge\frac8C(1-e^{-C\Delta}).$$

For a newly opened one-face-overlap target, the same argument gives

$$\boxed{\Pr_{\mu_p^+}(\text{newly opened target reacts by }\Delta)
\ge\frac{n-1}{3nC}(1-e^{-C\Delta}).}$$

These are events inside the original full process, not dynamics with
anchors frozen. For fixed finite $\alpha$, both bounds remain positive
as $n\to\infty$. If $\alpha$ grows with $n$, these particular integrated
bounds can vanish; the instantaneous Palm rates alone do not control that
joint limit.

## Counting fluctuations, without an arbitrary burst threshold

Let $N_p(T)$ count marked reactions. Since the marked jump kernel is
reversible, its two-event density at positive lag is
$\langle\lambda_p,e^{-t\mathcal L}\lambda_p\rangle_\mu$, where
$\mathcal L$ is the full positive generator. Hence

$$\operatorname{Var}N_p(T)=\bar\lambda_p T+
2\int_0^T(T-t)\,
\langle f,e^{-t\mathcal L}f\rangle_\mu\,dt,
\qquad f=\lambda_p-\bar\lambda_p.$$

The spectral covariance is nonnegative. Also, only the nineteen
interfering roots act on $f$, and their positive operator has norm at
most $2C$. Therefore $\langle f,\mathcal Lf\rangle\le2C\langle f,f\rangle$.
Jensen's inequality for the spectral measure gives

$$\langle f,e^{-t\mathcal L}f\rangle
\ge\operatorname{Var}(\lambda_p)e^{-2Ct}.$$

Since $\operatorname{Var}(\lambda_p)/\bar\lambda_p=8-\bar\lambda_p$,

$$\frac{\operatorname{Var}N_p(T)}{\mathbb E N_p(T)}
\ge1+2(8-\bar\lambda_p)
\left[\frac1{2C}-\frac{1-e^{-2CT}}{4C^2T}\right].$$

On times $T=n\tau$, $\tau>0$, at fixed $\alpha$ this implies a limiting
lower bound $1+8/C$ on the variance-to-mean ratio. At $\alpha=1$, that is
$1.02477$. This conservative local bound already excludes a
Poisson approximation at the level of these moments. It does not prove
convergence to a compound Poisson process or justify exchanging moment
and weak-distribution limits. Reducible conserved sectors can contribute
additional nondecaying covariance.

## Spontaneous trajectories and released event observers

Unconditioned runs use independent uniform raw links, 32 independent
trajectories per condition, 32 faces, and duration $T=n/4$. All original
channels remain enabled. The predicted mean count is
$144(n-1)/n$ per trajectory. A representative raw starting state, final
state, and every intervening reaction are retained for each condition.

For $\alpha=1$, consider a $\Delta=0.05$ window following each spontaneous
reaction. Events too close to the run endpoint are excluded. Uncertainty
uses independent whole trajectories, not an assumption of independent
events within a burst.

| Fiber size | Observed reactions | Another reaction in window | Independent-clock prediction | One-face-overlap reaction in window |
| ---: | ---: | ---: | ---: | ---: |
| 243 | 4,444 | $0.5757\pm0.0110$ | 0.1113 | $0.1841\pm0.0068$ |
| 729 | 4,630 | $0.5653\pm0.0147$ | 0.0387 | $0.1735\pm0.0081$ |
| 2,187 | 4,590 | $0.5671\pm0.0078$ | 0.0131 | $0.1804\pm0.0067$ |

The independent-clock comparison is $1-\exp(-\bar\lambda_{\rm total}\Delta)$,
not a fitted null model. Same-fan recurrence is only about $0.16$–$0.17$
in these windows. Neighboring-region activity is therefore not exhausted
by exact local undoing. All four saved window sizes, $0.01,0.05,0.2,1$,
are retained; the conclusion does not depend on defining a single burst
boundary.

At $n=729$, $\alpha=8$, another 4,622 spontaneous reactions give
$0.4994\pm0.0109$ for any following reaction in the same window, against
the unchanged independent-clock prediction $0.0387$.

Released Palm runs use 2,048 independent trajectories for each of
$n=27,243,2187$ at $\alpha=1$ and $n=1,000,000,007$ at $\alpha=1,8$.
At the finest fiber and $\alpha=1$, by time $0.2$:

- $0.2251\pm0.0092$ have reacted again on the original marked fan.
- $0.04395\pm0.00453$ have reacted on the specified one-face-overlap target.
- $0.6875\pm0.0102$ have had a reaction change charge outside the original
  fan's three faces.

The total conditional intensity decreases from about 25 initially to
$6.16\pm0.24$ at time $0.2$ and $0.066\pm0.030$ at time one. This is
evidence of dispersing interaction memory, not proof that its tail is
integrable or that the reference has a particular lifetime law.

## Conserved sectors explain much of the large raw count variance

The whole-run count variance is much larger than its mean, but a large
part has a static explanation. From

$$Q=F+N_Z-N_E,$$

a fixed sector with $F\le Q\le2F$ has, in its leading identity-free
ensemble, $K=2F-Q$ reflections and $Q-F$ rotations. The $K$ reflection
positions are uniform, with $K$ even. A marked triple is $RRR$ with
probability $(K)_3/(F)_3$, where $(x)_3=x(x-1)(x-2)$.

Its forward changing rate, conditional on being $RRR$, is
$24(n-1)/n^2$. Distinct additional identity constraints cost another
factor of $1/n$, so excluding initial identities does not change this
leading coefficient. Stationary backward flux equals forward flux within
each conserved $Q$ sector. Consequently, at fixed $F$,

$$\boxed{\bar\lambda_{\rm total}\mid Q
=\frac{144F}{n}\frac{(2F-Q)_3}{(F)_3}+O(n^{-2}).}$$

The rate is not a phenomenological parameter: the factor 144 comes from
three rooted fans per face, twelve reaction channels, the two possible
adjacent equalities, and equal forward/backward flux.

Under the original raw reference, $K$ tends to a binomial $(F,1/2)$ count
conditioned to be even. On $T=n/4$ with $F=32$, the variance of the
sector-dependent mean counts alone tends to $6212.44$, or $43.14$ times
the unconditional expected count 144. This can produce a very large
unstratified variance without implying a critical avalanche.

The observed whole-run variance-to-mean ratios at $\alpha=1$ are
$45.64,74.10,54.11$ for $n=243,729,2187$. Centering each trajectory by
the derived leading mean for its actual conserved charge reduces the
sample residual second moment divided by the mean to $7.20,7.24,6.46$.
At $n=2187$, the mean prediction error is $0.105\pm5.466$ events.
No rates were fitted to these counts.

Those residual ratios are not exact finite-$n$ conditional Fano factors:
the centering is asymptotic, there are only 32 independent long runs,
and other conserved or slowly varying data may remain. The local Palm
bounds and the observed short-window recurrence establish temporal
clustering separately from this sector-mixture issue.

## What remains to derive

The interaction memory is carried by actual loop-flatness relations,
not fixed anchors or an added attraction. Pair transport can move or
extend the loop carrying a relation; angular shifts can make it invisible
to a local equality check. The next mathematical task is to follow this
constraint through the full dynamics and determine whether its return
statistics give a finite collision memory or a longer-lived structure.

This turn does not establish a closed kinetic limit on times of order
$n$, a finite mean burst size, an invariant bound object, or a quantum
state space. In particular, multiple available reactions sharing one
identity must not be treated as independent branching particles. A
successful next description must preserve their competition, the conserved
charge, and the actual based holonomy relations.

The [integer-carrier derivation](fiber-constraint-dynamics.md) now follows
one inherited correlation under the full rules in the fixed-time Palm
limit. It proves alternating reaction directions, distinguishes spatial
spreading from nonzero angular offset, and identifies an auxiliary
branched-cover homology class that can survive after local reactivity is
lost. It does not yet close the multiple-constraint kinetic problem.

Apparatus and observations:
[released and spontaneous trajectories](../tools/fiber_reaction_bursts.py),
[trajectory data](../data/fiber-reaction-bursts.json),
[flatness-carrier census](../tools/fiber_reaction_carriers.py),
[carrier data](../data/fiber-reaction-carriers.json),
[conserved-sector prediction](../tools/fiber_reaction_sector_rates.py),
[sector comparison](../data/fiber-reaction-sector-rates.json).
