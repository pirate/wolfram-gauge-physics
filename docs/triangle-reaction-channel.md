# What a reaction actually remembers

The uniformly scheduled twelve-rule reaction bank has an exactly soluble
local channel. **A changing reaction retains the ordered boundary product,
but randomizes the remaining active local factorization.** Individual rules
remain invertible. Forgetting which rule occurred is essential to this result.

This is a result about our current finite classical gauge model, not a
derivation of quantum scattering, molecular binding, or emergent spacetime.
The uniform independent rule scheduler is an existing modeling choice, not a
law of physics derived here.

## Derivation from the primitive symmetry orbits

Use the automorphisms of the triangle, derived in the existing model, with
identity $e$, three reflections, and two nonidentity rotations. All three
holonomies are based at the same actual event root. Write $p=abc$ in the
table's ordered-product convention, and let

$$
\begin{aligned}
R_p&=\{(s,s,p),(p,s,s):s\text{ a reflection},\ s\ne p\},\\
E_p&=\{(a,b,c):abc=p,\ \{q(a),q(b),q(c)\}=\{0,1,2\}\}.
\end{aligned}
$$

For each reflection $p$, $|R_p|=4$ and $|E_p|=12$. The latter count is six
charge orders times two choices of rotation; the product then fixes the
remaining reflection. These are exactly the active inputs of the reactions.

Common conjugation and reversal $(a,b,c)\mapsto(c^{-1},b^{-1},a^{-1})$
form a twelve-element symmetry action. The twelve active reflection triples
form one free orbit: a triple contains two distinct reflections, so has no
nonidentity conjugation stabilizer, and reversal exchanges its adjacent
equality position. Each of the three mixed charge-order/reversed-order pairs
also forms a free twelve-element orbit.

Fix one active reflection triple $x$. Choosing any of the twelve mixed
triples $y$ with the same product defines a unique equivariant matching
$h x\leftrightarrow h y$ for every symmetry $h$. The matching is an
involution, extended by identity elsewhere. These twelve maps are exactly
the already-derived minimal reaction bank; this reconstructs the bank from
the orbit argument without selecting for any desired trajectory.

Consequently every edge of the complete bipartite graph
$R_p\leftrightarrow E_p$ appears in exactly one rule. This is the source of
the following exact channel, not a fit to simulation.

## Averaged channel and exact spectrum

For a uniform choice among the twelve reaction rules on this particular
support, in the order $(R_p,E_p)$,

$$
P_p=\frac1{12}
\begin{pmatrix}
0_{4\times4}&J_{4\times12}\\
J_{12\times4}&8I_{12}
\end{pmatrix},
$$

where $J$ is the all-ones matrix. Every $R_p$ state changes; an $E_p$ state
changes with probability $1/3$. Conditioned on a change,

$$
\Pr(Y=y\mid X=x,\text{change})=
\begin{cases}
1/12 &x\in R_p,\ y\in E_p,\\
1/4 &x\in E_p,\ y\in R_p.
\end{cases}
$$

For each fixed reaction direction, the input and output are conditionally
independent given $p$, for **any** distribution of the incoming active state.
This does not discard the direction/sector information itself.

The spectrum of $P_p$ is exactly

$$
1^{(1)},\quad (-1/3)^{(1)},\quad 0^{(3)},\quad (2/3)^{(11)}.
$$

The corresponding spaces are the constant vector; the vector equal to $3$
on $R_p$ and $-1$ on $E_p$; zero-sum contrasts within $R_p$; and zero-sum
contrasts within $E_p$. The last eigenvalue comes from unsuccessful reaction
attempts, not transport of an internal mode through a changing reaction.
The negative eigenvalue describes ordinary classical sector alternation.

At uniform equilibrium on these sixteen states, with sector means $\mu_R$
and $\mu_E$, the local Dirichlet form is

$$
\langle f,(I-P_p)f\rangle
=\frac14\left[\operatorname{Var}_{R_p}f+
\operatorname{Var}_{E_p}f+(\mu_R-\mu_E)^2\right].
$$

The complete three-loop gauge quotient has two active reflection orbits and
six active mixed orbits. Its active-block spectrum is
$1,-1/3,0,(2/3)^{(5)}$. The other 41 gauge orbits are fixed by reactions on
this support alone. In particular, this calculation gives **no global mixing
bound**: the inactive gates and overlapping supports matter.

Two independently selected changing reactions on the same isolated support,
$E_p\to R_p\to E_p$, return uniformly over all twelve $E_p$ states. The exact
return probability is $1/12$; the probability of restoring the charge order
is $1/6$; conditional on that charge order, the exact return probability is
$1/2$. Repeating the *same known involution* instead returns with probability
one. Those are different experiments.

## What survives in the surrounding gauge field?

The theorem does **not** say that a consumed rotation loses every physical
correlation. Its information can enter $p$. With any untouched external
holonomies transported to the same root, joint conjugation invariants such
as $q(sp)$ survive this reaction exactly. A raw reflection label $p$ alone
is not gauge invariant; its relationships to the environment are.

For example, in the repository's derived group indexing, the two inputs
$(e,1,3)$ and $(e,1,4)$ have boundary products $5$ and $2$. An external
reflection $s=2$ distinguishes them through $q(sp)=2$ versus $0$, respectively.
Every output reaction preserves that distinction. Conversely, at fixed $p$
the two mixed factorizations with any fixed charge order have identical
distributions of changing outputs. There is no residual channel distinguishing
those factorizations once the rule choice is unobserved.

External comparisons require actual connector paths unaffected by the event,
or explicit accounting for how those paths change. This local calculation
does not identify frames at distant vertices and does not replace the
shared-edge simulation by independent face variables.

The useful next object is therefore the **overlapping boundary-product
algebra**, including its relation to untouched loops, not a particle label
declared to survive every reversible reaction. Different supports do not
commute; the one-support randomization result must not be extrapolated to
an arbitrary intervening history.

## Separating retained pairs from apparent three-loop memory

For a rotation-only triple, define three gauge-invariant signs at its common
root:

$$
\phi_{ij}=\begin{cases}+1&z_i=z_j,\\-1&z_i\ne z_j.\end{cases}
$$

The four possible gauge states give orthogonal, zero-mean coordinates
$(\phi_{12},\phi_{13},\phi_{23})$ under the exact uniform conditional
reference already derived for rotation-only disks. For two such states,

$$
\sum_{i<j}\phi_{ij}(X)\phi_{ij}(Y)
=4\mathbf1_{[X]=[Y]}-1.
$$

Thus the complete three-dimensional rotation-memory statistic is precisely
the average of **three pairwise memories**. It is not an irreducible
three-body correlation.

The two-reaction calculation above makes this concrete without inserting a
randomization rule. Take two unchanged rotation spectators at the same root,
and a mixed input containing a third rotation. Apply the actual
$E_p\to R_p\to E_p$ channel and condition on restoration of the mixed charge
order. The returned rotation is equally likely to have either orientation;
the spectators are unchanged. Two pair contrasts therefore average to zero
and the spectator-pair contrast remains one: the normalized three-rotation
memory is **exactly $1/3$**. The calculation enumerates all 144 choices of
mixed input and spectator rotations, and every allowed two-rule history.
The spectator/connector assumptions describe this local experiment, not an
arbitrary intervening whole-mesh history.

An event-derived decomposition in `triangle_rotation_episodes.py` separates
these pair contributions using actual episode histories and actual measured
common-root holonomies. A reaction ending an episode does not assert the
destruction of all its information. A retained episode pair does not assume
its connector, or its measured relative orientation, remained unchanged.

## Whole-mesh check: the unexplained excess does not replicate

The earlier side-12 sample had 64 independent initializations. At
$\tau=4$ attempts per face, its reconfigured rotation-memory contribution
was $0.033390\pm0.007524$. The episode decomposition gives:

- Same original three episodes, same faces: $0.007806\pm0.003576$.
- Other initial episodes: $-0.000434\pm0.000434$.
- At least one reaction-created episode: $0.026018\pm0.007125$.

The last term splits into $0.016045\pm0.003016$ from the two surviving
original episodes at their original faces, and $0.009974\pm0.005148$ from
the remaining pairs. This is a **reanalysis of the earlier trajectories**,
not independent confirmation of their signal.

A fresh 256-initialization ensemble at the same side and the fixed time
$\tau=4$ gives:

| Contribution within patches containing a reaction-created episode | Mean | Standard error |
| --- | ---: | ---: |
| Total three-rotation memory | 0.015611 | 0.002917 |
| Pair of original episodes still at the same faces | 0.014961 | 0.001610 |
| Remaining pair contributions | 0.000650 | 0.002735 |

There were 144 qualifying rooted patches, of which 72 returned to the same
complete gauge orbit. Of their 432 pair contributions, 138 involved the
same two original episodes; all 138 retained the measured relative
orientation. The other 294 pair contributions had signed sum only $6$.
These overlapping patches are **not independent replicates**; all quoted
standard errors are calculated across independent initialized trajectories.

The fresh residual is consistent with zero. The earlier excess beyond
surviving pairs has not replicated. This does not prove the absence of all
reaction-mediated correlations, or establish a bound uniform in time or
system size. It does remove the current positive three-loop statistic as
evidence for a persistent collective object.

The independent run used initial seeds $79312000+i$, scheduler seeds
$89412000+i$, $i=0,\ldots,255$, the unchanged elastic bank, the exact
stationary full-$S_3$, reflection-present reference at $Q=F=288$, and
$4F=1152$ attempted updates per run. Every trajectory was replayed through
actual shared links against the compiled event sequence. Aggregate results
are in [the numerical record](../data/triangle-reaction-memory-summary.json).

## Consequence for the next mathematical detector

For $k$ rotations, write their two orientations as $\sigma_i=\pm1$ at a
common root. Gauge conjugation identifies simultaneous reversal of all
signs. The complete real function space therefore has the character basis

$$
\chi_A=\prod_{i\in A}\sigma_i,\qquad |A|\text{ even}.
$$

There are $2^{k-1}-1$ nonconstant characters. For $k=3$ they are all pairs;
the first degree-four character occurs at $k=4$. It is orthogonal to the
linear span of pair characters under the uniform conditional reference,
but its temporal correlation can still factor into independent pair
memories. A degree-four signal alone would not establish collective dynamics.

This identifies two concrete next questions: how conserved boundary products
carry correlations through overlapping reaction supports, and whether
four-loop return correlations have a component beyond separately measured
pair returns. Both require actual connector transport and a dynamical
baseline; neither requires inventing particle identities or force laws.

## Reproduce the exact algebra

From the repository root:

```sh
uv run --with numpy --with scipy --with python-flint python tools/triangle_reaction_channel.py
```

The calculation reconstructs all twelve primitive maps, resolves all three
fixed-product blocks, checks a full integer eigenbasis of $12P_p$, counts
the two-changing-reaction returns, and independently checks the exact gauge
quotient. It uses no numerical eigensolver and no fitted parameters.

To repeat the independent whole-mesh calculation:

```sh
OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 uv run --with numpy --with scipy --with python-flint python - <<'PY'
import json, sys
sys.path.insert(0, 'tools')
from triangle_rotation_episodes import ensemble
result = ensemble(trials=256, times=(0, 4), initial_seed=79312000, schedule_seed=89412000)
print(json.dumps(result['estimates'], indent=2))
PY
```
