# Exact charge response: diffusion, its failure, and the hidden variable

The primitive bank gives an exact discrete heat equation in its simplest
two-defect sector. Reactions spoil that closure. We can measure the failure
without fitting a continuum equation, and exhibit two actual connections with
the same charge field and initial drift but different **two-step mean charge**.
Adding current reaction-activity flags still fails to close evolution: a second
pair agrees through two-step mean response but separates after three attempts.

These are finite-model results on the supplied triangular torus, not a derivation
of spacetime or electromagnetism. The scalar charge is derived from the finite
fiber symmetry; it has not been identified with electric charge. Uniform random
rule selection and the mesh are supplied. No wave law or force has been added.

![Exact projected response, dilute correction, and actual hidden-gauge response](images/triangle-charge-response.png)

## 1. Keep the actual discrete clock

There are $F=2s^2$ faces and $M=39F$ rooted operators: vacancy transport and twelve
reaction involutions, including all no-ops. Write their Markov operator as

$$P=\frac1M\sum_a T_a,\qquad \mathcal L=M(P-I),\qquad D(U)=\mathcal L q(U).$$

$\mathcal L$ is a scaled algebraic difference operator, not a new continuous-time
evolution. One physical simulation tick here is one attempted operator. The
[exact event clock](triangle-event-clock.md) preserves those attempts.

Use the [exact raw-connection reference](triangle-reference.md), conditioned on
fixed $Q=\sum_f q_f$, full $S_3$ based-loop image, and at least one reflection face.
That subset is invariant. Uniform weight on its raw connections is stationary and
reversible because each operator is an involution. It need not be one communicating
component; no mixing assumption enters the calculation.

Let $X=q-(Q/F)\mathbf1$. Exchangeability of the exact reference gives

$$S=\mathbb E[XX^T]=\chi(I-J/F),\qquad
\chi=\frac F{F-1}\operatorname{Var}(q_f).$$

The conserved constant mode is removed throughout; an inverse of the full singular
covariance matrix is never required.

## 2. The projected drift follows from actual local supports

For populations $(N_0,N_1,N_2)$, an ordered assignment of $k$ distinct faces with
counts $(c_0,c_1,c_2)$ has probability

$$p(c_0,c_1,c_2)=
\frac{\mathbb E[(N_0)_{c_0}(N_1)_{c_1}(N_2)_{c_2}]}{(F)_k},$$

where $(n)_k$ is the falling factorial. The expectation uses connection
multiplicities, **not uniform weighting of charge fields**.

The [pointwise drift](triangle-charge-current.md) is a sum of vacancy differences
and $4(1-q_i)$ on each rooted fan whose charges are a permutation of $(0,1,2)$.
Three-reflection fans have zero first drift, even when active. Define

$$a=2[p(1,1,0)+4p(1,0,1)],\qquad b=8p(1,1,1).$$

Here $L_1$ is the degree-three nearest-neighbor Laplacian of the honeycomb dual,
and $L_2$ its degree-six distance-two Laplacian. Direct enumeration of the actual
fan supports proves

$$\sum_{\text{rooted fans}}L_{K_3}=4L_1+L_2,\qquad L_2=6L_1-L_1^2.$$

Consequently

$$K=-\mathbb E[DX^T]=(a+4b)L_1+bL_2,\qquad
\Gamma=K/\chi=\kappa_1L_1+\kappa_2L_2,$$

$$\kappa_1=(a+4b)/\chi,\qquad\kappa_2=b/\chi.$$

$-\Gamma X$ is the best linear projection of the actual drift in this stationary
reference. It is not generally the drift itself. With cell wavevector $(k_x,k_y)$,
the two eigenvalues of $L_1$ and the corresponding projected decay symbols are

$$\lambda_\pm=3\pm|1+e^{ik_x}+e^{-ik_y}|,\qquad
g_\pm=\kappa_1\lambda_\pm+\kappa_2(6\lambda_\pm-\lambda_\pm^2).$$

In particular,

$$g_- =\frac{\kappa_1+6\kappa_2}{3}
(k_x^2+k_y^2+k_xk_y)+O(|k|^4).$$

This is a diffusion-type **initial projected response**, not a proven long-time
diffusion coefficient or a wave frequency. Fourier complex numbers are an analysis
device here, not microscopic quantum amplitudes.

## 3. An exact heat equation, then an exact obstruction

At $Q=2$ in this reflection-present subset there are exactly two charge-one faces
and no charge-two faces. No reaction is possible. Summing the two rooted vacancy
operators on each dual edge gives the pointwise identity

$$D=-2L_1q,\qquad
\mathbb E[q_{t+1}]-\mathbb E[q_t]=-{2\over M}L_1\mathbb E[q_t].$$

It holds for every initial law supported in that subset, not just equilibrium.
The familiar closure of mean occupation in symmetric exclusion is established
mathematics, not a new general discovery; see
[van Ginkel and Redig's exclusion/heat-equation analysis](https://link.springer.com/article/10.1007/s10955-019-02420-2).
The result here identifies that behavior and its exact clock in the derived bank.

For reactive sectors, define

$$H=\mathbb E[DD^T],\quad R=D+\Gamma X,\quad
B=\mathbb E[RR^T]=H-\chi\Gamma^2\succeq0.$$

Orthogonality $\mathbb E[RX^T]=0$ proves the last identity. For stationary lagged
covariance $C_n=\mathbb E[X_nX_0^T]$, reversibility gives

$$C_1=S-K/M,\qquad C_2=S-2K/M+H/M^2.$$

A closed linear model would instead predict

$$C_2^{\rm linear}=(I-\Gamma/M)^2S,\qquad
\boxed{C_2-C_2^{\rm linear}=B/M^2}.$$

The positive-semidefinite correction is an exact two-step obstruction to linear
closure. At $Q=2$, $B=0$. At $Q=4$ and $Q=F$, it is nonzero on all three audited
sizes. This first residual arises from **nonlinear charge configurations**;
$D$ itself is charge-field-only. It is not solely a hidden-gauge noise term.

Projection and memory are established coarse-graining ideas; see
[te Vrugt and Wittkowski's pedagogical account](https://arxiv.org/abs/2001.01572).
The equations above are elementary finite reversible-Markov identities specialized
to this bank, not a claim to have invented projection theory.

## 4. The dilute four-charge correction is algebraic, not fitted

For $Q=4$ the two allowed populations are $(F-4,4,0)$ and $(F-3,2,1)$, with exact
probabilities $10(F-3)/(10F-3)$ and $27/(10F-3)$. Put

$$p_F=20F^2-59F+24,\qquad d_F=F(F-1)(10F-3).$$

Then

$$\chi={2p_F\over d_F},\quad
a={4(F-3)(20F+1)\over d_F},\quad
b={432(F-3)\over d_F(F-2)},$$

$$H_{ff}={48(40F^3-198F^2+3659F-7575)\over d_F(F-2)},$$

$$B_{ff}={1296(820F^4-4479F^3-16993F^2+146682F-229896)
\over d_F(F-2)^2p_F}.$$

The numerator has strictly positive coefficients when expanded in $F-18$.
The code clears denominators in the identity
$B_{ff}=H_{ff}-(12a^2+132ab+378b^2)/\chi$ using integer polynomials.
Their leading coefficients give

$$F^3B_{ff}\longrightarrow {26568\over5},\qquad
F^2{B_{ff}\over H_{ff}}\longrightarrow {1107\over40},$$

$$F^2(\kappa_1-2)\longrightarrow {81\over2},\qquad
F^2\kappa_2\longrightarrow {54\over5}.$$

These limits run along supplied tori $F=2s^2$, $s\ge3$, with **fixed total charge**,
not fixed density. Vanishing instantaneous residual per face does not prove that
its accumulated long-time effect vanishes or establish a hydrodynamic limit.

## 5. Which hidden variable affects the next response?

For any function $f$ of the entire charge field, the exact microscopic operator is

$$\begin{aligned}
\mathcal L(f\circ q)(U)={}&
\sum_{p:\,\text{one vacancy}}[f(q^{p,\mathrm{swap}})-f(q)]\\
&+4\sum_{p:\,q_p\in\mathrm{Perm}(012)}[f(q^{p,111})-f(q)]\\
&+2\sum_{p:\,q_p=111}A_p(U)
\sum_{\sigma\in\mathrm{Perm}(012)}[f(q^{p,\sigma})-f(q)].
\end{aligned}$$

All sums retain **rooted** patch multiplicities. The activity flag uses the
three holonomies $(h_0,h_1,h_2)$ transported to the patch's common basepoint:

$$A_p=\mathbf1\{(h_0=h_1)\mathbin{\mathrm{xor}}(h_1=h_2)\},\qquad q_p=111.$$

No unrelated local frames are identified. Local gauge transformations conjugate
all three together and preserve the flag. All 216 local inputs reproduce the
transition multiset: an active $111$ has all six $012$ permutations with multiplicity
two; each $012$ has four returns to $111$. This formula determines the next charge
transition law from $(q,A)$. The joint $(q,A)$ process, however, is not closed, as
the next counterexample shows.

Applying it to $f=D$ gives $E(U)=\mathcal L D(U)$ and the pointwise two-attempt mean

$$\mathbb E[q_2\mid U]=q(U)+2D(U)/M+E(U)/M^2.$$

The existing [same-boundary braid-memory preparations](triangle-charge-current.md#3-the-controlled-memory-readout-isolates-the-fluctuation-effect)
one and two have identical $q$ and $D$. Their activity masks differ at rooted fan
468. In face order $(107,130,131,132,133,156,157,159)$,

$$E(U_1)-E(U_2)=(16,-24,16,32,-48,-24,16,16),$$

with zero difference on every other face. Here $M=11232$. The sum is zero, as charge
conservation requires. Thus even retaining both the **full charge field and its
first drift** fails to predict the next drift, and the mean charge difference is
visible after exactly two attempts. There is no stationarity assumption in this
witness. It strengthens the earlier one-step noise distinction to a mean-response
distinction, without claiming a propagating quantum phase.

## 6. Current activity is insufficient: an exact three-attempt witness

A bounded search using the exact raw-connection sampler finds a pair on 18 faces
with $Q=4$. Both have charge-one faces $(5,7,8,11)$, all other charges zero, and
**no active reflection fans anywhere**. Both are in the full $S_3$ reflection-present
subset. They have identical next-charge transition laws, identical $D$, and identical
$E=\mathcal L D$.

Nevertheless, the exact next-state laws of $(q,A)$ differ. For one charge target,
one preparation activates fan 41 while the other does not. For another charge
target, the difference is activation of fan 31; fan 32 is active in both. Each
transition difference has multiplicity two, and the charge marginals cancel exactly.
The total variation distance between the joint transition laws is $2/351$ with
$M=702$ rooted operators. All no-ops are included in that comparison.

The same data determine the first difference in mean charge. Expand the **actual**
discrete operator, without a continuous-time approximation:

$$P^3q=q+{3D\over M}+{3E\over M^2}+{\mathcal L E\over M^3}.$$

Since $E$ is a function of $(q,A)$, summing it against the joint-target count
difference gives, in face order $(2,4,5,6,7,12)$,

$$\mathbb E[q_3\mid U_1]-\mathbb E[q_3\mid U_2]
={1\over702^3}(-32,-32,112,32,-112,32),$$

and zero on the remaining faces. The full means agree at attempts zero, one, and
two, then differ at attempt three. Thus the distinction is not just an irrelevant
microscopic change in an auxiliary flag.

The checked-in witness includes both 27-link connections, the deterministic search
seed and sample indices, and all four differing joint-target counts. An independent
two-level enumeration through raw-link surgery reproduces the third-response
numerator **without using the activity-based generator**. Complete C++ operator/
inverse scans check the joint-target law, also after nonconstant local frame changes.

## Verification and next question

The covariance calculation groups products of local drift terms by their overlap
pattern and exact charge populations. Disjoint terms cancel by exchangeability;
overlaps use at most six faces. Translation leaves two orientation prototype rows.
This is a bounded exact diagnostic, not a million-node evolution implementation.

An independent census on 18 faces checks all 5,508 allowed four-charge fields with
their correct connection weights (3,060 fields of weight 960 and 2,448 of weight
216). It reproduces $S$, $K$, and $H$ exactly. Independent raw-link surgery and a C++
scan of every operator followed by its inverse verify charge-target multiplicities.
Nonconstant local gauge transformations preserve the two-step witness.

```bash
uv run --with numpy --with scipy --with python-flint python tools/triangle_charge_response.py --output out/triangle-charge-response.json
diff -u data/triangle-charge-response.json out/triangle-charge-response.json
uv run --with numpy --with scipy --with python-flint python tools/triangle_response_memory.py --output out/triangle-response-memory.json
diff -u data/triangle-response-memory.json out/triangle-response-memory.json
uv run --with numpy --with scipy --with python-flint python tests/triangle_charge_response.py
uv run --with numpy --with matplotlib python tools/plot_triangle_charge_response.py
```

The next closure test is concrete: retain transported relative-loop relations
before a three-reflection encounter becomes active, and test whether they predict
the activity flags' transition law. The three-attempt witness is now a mandatory
negative control for any proposed reduction. Do not replace measured memory with
a fitted force or promote the projected diffusion symbol into a Hamiltonian.
