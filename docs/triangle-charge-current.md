# Conserved currents, hidden noise, and homogeneous-field activation

The triangle-fiber primitives supply more than a conserved total: they define an exact
local scalar-charge current and a computable fluctuation law. A useful distinction is
now explicit: **the full mean charge/current drift factors through the charge field,
but its conditional covariance does not**. Stored gauge memory can change fluctuations
while leaving the entire charge field, mean drift, and encounter-boundary holonomy fixed.

The autonomous control experiment starts both connections with $q_f=1$ on every face.
One stays frozen under every rule. The other, with specified microscopic reflection
disorder, develops vacancies and charge-two faces through the derived reactions, then
supports vacancy transport. Neither has an initial mean charge current. All observed
charge changes are reconstructed exactly from integrated primitive currents.

This is classical gauge-dependent activation and conserved stochastic transport on
supplied geometry. It is not electromagnetic charge, quantum vacuum fluctuations,
spontaneous symmetry breaking from identical raw states, a demonstrated phase transition,
or emergent matter.

![Actual homogeneous-field evolution and exact gauge-noise witness](images/triangle-charge-current.png)

## 1. The current follows the actual written links

Orient each dual edge from the smaller face id to the larger one and let $B$ be its
incidence matrix, with entries $-1,+1$ at source and destination. A vacancy update
writes one primal link, which crosses one dual edge. If its two local charge increments
are $(\delta q_0,\delta q_1)$, its directed local current is $j=-\delta q_0$.

A three-face fan writes two internal primal links, giving a dual path through its three
faces. The path currents are uniquely fixed by continuity:

$$j_{01}=-\delta q_0,\qquad j_{12}=\delta q_2,$$
$$\delta\mathbf q=(-j_{01},j_{01}-j_{12},j_{12}).$$

The code checks that these dual edges cross **exactly** the primitive's written links.
It does not choose an arbitrary route through the whole mesh or add a circulation.
After orienting each local current consistently, every actual update obeys

$$\mathbf q_{t+1}-\mathbf q_t=B\mathbf j_t,\qquad
\mathbf q_T-\mathbf q_0=B\sum_{t<T}\mathbf j_t.$$

Local gauge-frame changes leave these scalar quantities invariant. Their interpretation
depends on the specified primitive support; this does not define an electromagnetic
current or physical time.

## 2. Exact first and second moments of all twelve reactions

Use the [derived charge](triangle-feedback.md) $q(1)=0$, $q(\text{reflection})=1$,
$q(\text{three-cycle})=2$. For a based local triple $h=(a,b,c)$, define the **unnormalized**
bank sums

$$D(h)=\sum_{r=1}^{12}\delta\mathbf q_r(h),\qquad
A(h)=\sum_{r=1}^{12}\delta\mathbf q_r(h)\delta\mathbf q_r(h)^T.$$

If its charge tuple is a permutation of $(0,1,2)$, put $d_i=1-q_i$. Exactly four rules
convert it to three reflections, giving

$$D=4d,\qquad A=4dd^T.$$

If all three charges are one, the gate is

$$g(h)=\mathbf1\{(a=b)\mathbin{\mathrm{xor}}(b=c)\}.$$

An active triple has six outgoing charge permutations, each with multiplicity two:

$$D=0,\qquad A=4g(h)(3I-J_3).$$

Every other input has $D=A=0$. All 216 raw triples are checked against every table
entry. In particular, among the 27 three-reflection triples, twelve are active.

Conditional on selecting the reaction bank at this one patch, divide these sums by
twelve. Thus an active three-reflection encounter has charge covariance
$I-J_3/3$, and its two path-current covariance is

$$\frac1{12}\begin{pmatrix}8&4\\4&8\end{pmatrix}.$$

The inactive encounter has the same charge tuple and zero covariance. For the reverse
$(0,1,2)$ encounter, the centered covariance is $2dd^T/9$, not its uncentered second
moment $dd^T/3$.

Vacancy transport has an equally explicit drift: $(q_1-q_0,q_0-q_1)$ if exactly one
face is vacant, zero otherwise. Consequently the full drift depends only on the local
charge field. The three-reflection noise retains information about relative based
holonomies that the charges discard.

This is a pointwise statement about the generator. It does **not** close an equation
for $\mathbb E[\mathbf q_t]$: averaging its nonlinear local functions requires charge
correlations, whose evolution can depend on hidden gauge information. No deterministic
diffusion equation, Gaussian noise, or white-noise limit has been assumed.

## 3. The controlled memory readout isolates the fluctuation effect

Reproduce all four [prepared braid-memory states](triangle-feedback.md#3-read-out-the-previously-generated-braid-memory-on-a-real-mesh)
and enumerate all 11,232 rooted operators from each, including no-ops. States one and
two, using zero-based numbering, have the same full charge field and the same encounter
boundary product. Their entire mean charge and current fields coincide.

On encounter faces $(130,133,156)$ their conditional **centered** charge covariance
difference is exactly

$$\Sigma_1-\Sigma_2=\frac1{11232}
\begin{pmatrix}8&-4&-4\\-4&8&-4\\-4&-4&8\end{pmatrix},$$

and vanishes elsewhere. The common drift makes the centering terms cancel. This is not
just a difference between spatially distinct configurations with equal global counts.
Nonconstant local gauge transformations reproduce the full moment arrays.

The analytic formulas agree with independent raw-link enumeration. A bounded C++ event
snapshot stream checks every changed intermediate connection in a complete operator/
inverse scan; each scan restores its original links exactly. The default runner output
is unchanged. The optional `--raw-events` diagnostic stops after eight million emitted
link values rather than allowing unbounded trajectory output.

## 4. A uniform charge field can be frozen or dynamically active

Assign the same reflection $h$ to every primal link. Every face holonomy is $h$ and
every based fan triple is $(h,h,h)$. There are no vacancies, and every reaction gate is
inactive. **Every primitive fixes the complete raw connection**, on any of the supplied
triangular tori.

Now instead assign each link independently one of the three reflections. Each triangle
still has reflection holonomy: the product of three odd permutations is odd, and every
odd element of $S_3$ is a reflection. Hence the whole observed charge field is still
exactly $q_f=1$, and the mean charge and current drift still vanish everywhere.

The based triples need not agree. Exhausting all $3^7=2187$ assignments on the seven
links of a fan gives all 27 reflection triples with multiplicity 81. Thus a fan is active
with probability $4/9$ in this **specified initial ensemble**. Since twelve of the thirteen
equally weighted rules are reactions, the initial ensemble has exact one-attempt rate

$$\Pr(\text{changing update})=\frac{12}{13}\frac49=\frac{16}{39}.$$

Every changing initial reaction creates one charge-zero and one charge-two face from
three charge-one faces. Its squared charge contrast increases by two, so

$$\mathbb E\!\left[\sum_f(q_f(1)-1)^2\right]=\frac{32}{39}.$$

These are initial-ensemble identities, not stationary rates or a fitted long-time law.
Overlapping fan gates are not independent.

Four side-six runs and one side-twelve run use the same attempted-update clock per
rooted operator: 10,000 and 40,000 attempts respectively. Every run includes the frozen
connection and a reaction-disabled control of the disordered connection. Both controls
remain exactly fixed. In the full bank, a reaction initiates evolution, followed by
charge transport and further reactions. Final $(N_0,N_1,N_2)$ populations are
$(15,42,15)$, $(17,38,17)$, $(19,34,19)$, $(18,36,18)$ on side six, and $(70,148,70)$
on side twelve. Exact conservation forces $N_0=N_2$ throughout, not just at the end.

All raw histories, selected snapshots, and integrated currents are exported. The image
uses supplied mesh coordinates and the final side-twelve snapshot. These finite runs
do not prove equilibrium, criticality, or binding. Because every primitive is an
involution, a nontrivial component cannot enter a globally fixed state: reversing such
a last step would contradict fixedness. The controls therefore demonstrate dynamically
disconnected behavior, not an observed transition into an absorbing state.

## 5. Predictable fluctuation accounting on autonomous histories

For an integer-weighted charge observable $X=w^T\mathbf q$, let $M$ be the number of
rooted operators, $b(s)=w^TD(s)$ and $a(s)=w^TA(s)w$, now summed over the full generator.
Under the supplied independent uniform scheduler,

$$\mu(s)=b(s)/M,\qquad v(s)=a(s)/M-b(s)^2/M^2.$$

The compensated process and its predictable variance are

$$Z_T=X_T-X_0-\sum_{t<T}\mu(s_t),\qquad V_T=\sum_{t<T}v(s_t).$$

Conditional expectation gives zero mean increments of $Z$, and $Z_T^2-V_T$ is also a
martingale. This is the standard martingale/Doob decomposition, not a new stochastic
calculus result; see [Varadhan, chapter 5](https://math.nyu.edu/~varadhan/course/PROB.ch5.pdf).
Here the local generator formulas supply the drift and variance from the actual gauge
primitives instead of fitting them to the observed path.

For a specified half-torus region in all four earlier 100,000-attempt histories, the
export records these quantities with exact integer numerators. Constant-state intervals
are weighted by **every attempted tick**, including no-ops and the changing tick at the
end. An independently expanded tick-by-tick test checks the interval accounting. Realized
currents reconstruct the full final charge field, and total-charge weights give exactly
zero drift and variance. The four observed sums of squared region-charge changes are
88, 77, 87, and 91; their predictable raw-square sums are respectively
$247058/2808$, $219246/2808$, $261262/2808$, and $259432/2808$.

Agreement in a few pilots is not a calibrated discovery test. These diagnostics are
inputs for future concentration bounds and controlled coarse-scale comparisons, not
p-values, equilibrium error bars, or proof of hydrodynamics.

## Reproduce

```bash
cmake --build build -j 4
export OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
uv run --with numpy --with scipy --with python-flint python tools/triangle_charge_current.py \
  --output out/triangle-charge-current.json
diff -u data/triangle-charge-current.json out/triangle-charge-current.json
uv run --with numpy --with scipy --with python-flint python tests/triangle_charge_current.py
uv run --with numpy --with scipy --with matplotlib python tools/plot_triangle_charge_current.py
```

All stored moment sums, histories, and continuity checks use exact integer arithmetic.
Only the presentation plot converts values to floating point.
