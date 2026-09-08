# Exact mode-loss kernels and a saturation obstruction

The full primitive generator makes a sharper statement than the recorded trajectory:
**a state with the maximum $Q/2$ above-band modes must lose modes whenever a vacancy
move acts**. At $Q=2$, every reaction is inert, so the very next changing primitive
update destroys the single possible above-band mode. This is an obstruction to treating
these particular modes as persistent particles, not a no-go theorem for every localized
structure or every rule family.

We also computed exact one-step kernels at all recorded reaction endpoints. Even
mode count, both orientation charges, and both defect populations together do not close
the dynamics: two states reached from the same seed have different next-count distributions
despite agreeing on all five quantities.

## 1. Saturated spectral capacity is fragile under transport

The [orientation-resolved spectral theorem](triangle-modes.md#4-a-stronger-bound-from-the-two-face-orientations)
gives

$$r:=n_+(L-12I)\le\min(Q_\uparrow,Q_\downarrow)\le Q/2.$$

Suppose $r=Q/2$. Necessarily $Q_\uparrow=Q_\downarrow=Q/2$. The vacancy rule exchanges
one nonidentity face holonomy with an adjacent identity face. Adjacent faces have
opposite orientations, so transporting charge $q>0$ gives

$$(Q_\uparrow',Q_\downarrow')=(Q/2-q,Q/2+q)$$

or its reversal. Applying the same theorem after the update yields

$$\boxed{r'\le Q/2-q=r-q.}$$

No eigensolver or stochastic assumption enters this argument. Exhausting all 36 raw
pair inputs verifies that each of the ten changing inputs transports one charge,
and every saturated-state vacancy transition in the kernel dataset satisfies the bound
with an exact characteristic-polynomial count. The statement concerns the specified
triangle-fiber bundle and its vacancy rule, not arbitrary graph dynamics.

At global $Q=2$, positivity makes every local triple have charge at most two. The full
twelve-rule reaction bank fixes every such raw triple, checked exhaustively. If a
mode exists, both orientation charges are one, hence there are two reflection defects
on opposite orientations. Every changing primitive update is then a vacancy move of
charge one, forcing $r'=0$.

For the side-six single-reflection-link control, exactly eight of the 2,808 equally
weighted rooted operators change the graph. Each eliminates its one above-band mode;
all other attempts leave the complete state unchanged. Under the supplied independent
uniform scheduler, the **first-loss** waiting time is therefore exactly geometric:

$$\Pr(T>t)=(1-8/2808)^t,\qquad\mathbb E[T]=2808/8=351.$$

These are attempted-update ticks, not physical time. Later mode reappearance is possible;
this calculation neither tracks particle identity nor proves an exponential law for
other initial states. The denominator includes empty regions and inactive rules, so
comparisons across mesh sizes must account for this supplied clock.

## 2. Enumerate every next update, not only observed transitions

The dataset includes four specified seeds and both endpoints of all twelve reactions
from all four recorded runs: 28 labeled candidates, including repeated states. There is
no selection by a favorable spectral effect. The complete source histories are replayed
before extracting their endpoints.

For each candidate, enumerate all thirteen rules at all 216 rooted supports. Every
no-op contributes to the denominator. Deduplicate equal raw targets only to avoid
recomputing their exact integer characteristic polynomials; **operator multiplicities
remain in the counts**. Export the complete next-state summary distribution separately
for vacancy transport and reactions.

Each distinct changed raw target is independently obtained as a C++ final connection.
A second C++ scan applies every operator followed by its inverse, checking all table
targets and exact restoration. Local gauge-frame transformations of the nonclosure
witness reproduce the entire kernel, including its operator-resolved changes.

For example, consider the initial noncommuting seed and its state after the recorded
seed-zero update at tick 306. Both have

$$(r,Q_\uparrow,Q_\downarrow,N_R,N_C)=(2,2,2,4,0).$$

Their exact kernels differ:

| Complete input state | Next count 1 | Next count 2 | Loss through vacancy / reaction |
| :--- | ---: | ---: | :--- |
| Initial noncommuting seed | $16/2808$ | $2792/2808$ | 16 / 0 operators |
| Same run, after tick 306 | $28/2808$ | $2780/2808$ | 12 / 16 operators |

Thus there cannot be a single exact next-step kernel depending only on those five
summary values for all raw states in this reachable component. In Markov-chain terms,
this partition is not strongly lumpable. This does **not** by itself disprove a useful
approximate or stationary effective model. Nor does it isolate hidden holonomy as the
cause: these two inputs need not share their full spatial face-class field. The earlier
[braid-memory readout](triangle-feedback.md) is the separate controlled result that holds
that field fixed.

## 3. Zero one-step loss is not protection

The commuting two-link seed has one above-band mode and loses none under **any** of
the 2,808 first-step operators. That alone might look like protection. Exhausting the
second-step kernels of its eight distinct changed targets shows otherwise:

$$\Pr(T\le2\text{ attempted updates})=\frac{48}{2808^2}=\frac1{164268}.$$

There is a C++-replayed two-vacancy-move witness ending at $(r,Q_\uparrow,Q_\downarrow)
=(0,0,4)$, so two is the exact minimum number of primitive attempts needed for loss.
The second-step enumeration retains the multiplicities of all first-step operators;
no-op first steps contribute zero because the initial kernel has no loss.

On the clock that counts only **changing** updates, the same exact calculation gives

$$\Pr(T\le2\text{ changing updates})=
\sum_y\frac{m_{x\to y}}{16}\frac{\#\text{loss operators at }y}
{\#\text{changing operators at }y}=\frac{13}{80}.$$

Each of the eight changed targets has multiplicity two; its changing-operator count
is either 16 or 20. The first expression is about six in a million; the second is
16.25 percent. Neither is a physical lifetime. Their difference illustrates how many
empty or inactive attempts can make a fragile structure appear long-lived. No fitted
decay law is needed to diagnose this effect.

## What this changes about the next experiment

A single fitted decay constant for a mode count would discard state dependence that is
already visible exactly. Persistence measurements must retain spatial and holonomy
information, distinguish attempted ticks from changing events, and distinguish survival
of a structure from repeated destruction and recreation. Saturated above-band modes
cannot provide a transport-stable particle count for the present vacancy rule.

These results neither introduce a wave equation nor derive matter. They narrow the
search: investigate other structural observables and nonsaturated sectors, or derive
and justify genuinely different primitives, rather than hiding loss with an imposed
particle-conserving update law.

## Reproduce

```bash
export OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
uv run --with numpy --with scipy --with python-flint python tools/triangle_mode_rates.py \
  --output out/triangle-mode-rates.json
diff -u data/triangle-mode-rates.json out/triangle-mode-rates.json
uv run --with numpy --with scipy --with python-flint python tests/triangle_mode_rates.py
```

The exported kernel counts and certificate are exact integers. There are no numerical
eigenvalue tolerances or fitted rates in this dataset.
