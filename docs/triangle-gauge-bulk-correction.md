# The hidden-gauge transport correction survives the volume limit

The [hidden current corrector](triangle-gauge-current-corrector.md) now has
a strictly positive **bulk variational improvement**, not just a
finite-volume witness. The proof combines an exact boundary-sensitive disk
marginal with a local operator bound. It concerns the existing classical
finite-fiber model on its supplied triangular lattice, not emergent
spacetime, quantum dynamics, or molecular binding.

The explicit lower bound on the extra improvement is only about
$1.07\times10^{-57}$ at density one. Its purpose is to establish a
nonvanishing volume-limit contribution, **not** a practically significant
estimate of its size. A subsequent
[stratified stationary calculation](triangle-gauge-resolved-transport.md)
resolves a much larger finite-volume correction on 72 and 128 faces,
but does not establish a sharper numerical bulk value.

## Count a disk without discarding its boundary alignment

Let $D$ be a fixed contractible disk of $r$ triangular faces, $v_D$ vertices,
and $e_D=v_D+r-1$ edges, embedded without wrapping in a large torus. Fix a
raw-link assignment $x_D$, its total disk charge $q_D$, and its oriented
boundary holonomy $g_D$. Write

$$p_0=1+3z+2z^2,\qquad p_1=1-3z+2z^2,\qquad p_2=1-z^2.$$

The corresponding characters are $1$, permutation sign, and
$\chi_2(g)=\#\operatorname{Fix}(g)-1$, with dimensions $1,1,2$.
Summing over the complement's handles, while retaining its boundary
product, gives

$$\mu^{\rm all}_{F,Q}(U|_D=x_D)
=6^{1-v_D}\,
\frac{[z^{Q-q_D}]
\left(p_0^{F-r}+\operatorname{sgn}(g_D)p_1^{F-r}
+\frac{\chi_2(g_D)}2p_2^{F-r}\right)}
{[z^Q]\left(p_0^F+p_1^F+p_2^F\right)}.$$

Here $Q$ is even. This is a counting variable, not a physical temperature
or a modification of the update bank.

One way to obtain the complement factor is to start with the direct
commutator count on the two handles and convolve with the central face
weight $z^{q(h)}$ for each remaining face. Its character eigenvalues are
$p_\chi$; the handle count at prescribed product $g$ is
$6\sum_\chi\chi(g)/d_\chi$. Restoring outside vertex frames produces
$6^{V-v_D+1}$ before division by the full raw count $6^V\sum p_\chi^F$.
This is a specialization of the classical finite-group surface-counting
framework, not a new character-theory identity; see
[Snyder's account](https://arxiv.org/abs/math/0703073).

An independent integer recurrence, starting with all 36 actual handle
pairs, agrees with all 486 boundary/charge coefficients for zero through
eight complement faces. Summing the one-face disk probabilities over all
six holonomies also gives exactly one at $F=18,72,128$.

For a useful event, include every transformation of $x_D$ by vertex
frames with one root frame fixed. This action is free on a connected
disk, so its orbit contains exactly $6^{v_D-1}$ assignments. The
probability of that event, denoted $\mathcal O_D$, is the displayed
coefficient ratio **without** the prefactor $6^{1-v_D}$.

This orbit is used only to count a gauge-invariant event. We do not evolve
independently quotiented patches, glue unrelated boundary frames, or
assume independence of overlapping events. The boundary character remains
in the exact finite-volume count.

## Its fixed-density probability stays positive

For $Q/F\to\rho\in(0,2)$ along admissible even charges, choose $z_\rho>0$
by

$$\rho=\frac{z_\rho p_0'(z_\rho)}{p_0(z_\rho)}
=\frac{z_\rho(3+4z_\rho)}{1+3z_\rho+2z_\rho^2}.$$

The disk boundary obeys $\operatorname{sgn}(g_D)=(-1)^{q_D}$.
Since $p_1(z)=p_0(-z)$, the two leading coefficient contributions add,
rather than cancel, for even total $Q$. The absolute $p_2$ contribution
is bounded using $(1+z^2)^{F-r}$, which is exponentially smaller than the
$p_0$ contribution at any fixed interior density.

For completeness, the remaining coefficient ratio is a local limit
calculation for independent variables with probabilities
$\Pr(Y=j)=a_jz^j/p_0(z)$, $a=(1,3,2)$. With $z=z_{Q/F}$, their mean is
$Q/F$, their variance is positive, and their integer span is one.
The local central limit estimate for sums of $F$ and $F-r$ such variables
at arguments differing by a fixed amount gives

$$\frac{[z^{Q-q_D}]p_0(z)^{F-r}}{[z^Q]p_0(z)^F}
\longrightarrow\frac{z_\rho^{q_D}}{p_0(z_\rho)^r}.$$

These independent variables are counting devices, not independently
evolving physical faces. Thus

$$\boxed{\lim_{F\to\infty}\mu^{\rm all}_{F,Q}(\mathcal O_D)
=\frac{z_\rho^{q_D}}{p_0(z_\rho)^r}>0.}$$

## A local witness forces an extensive elastic energy

Recall $b_3$ is three times the unnormalized winding-current drift,
$u=Lb_3$, and $S=\langle u,Eu\rangle$. The explicit witness has:

- A $5\times6$ cell rectangle: 60 faces, 101 edges, 42 vertices.
- Disk charge 60, boundary holonomy a nonidentity rotation, and full
  based-loop group $S_3$ inside the disk, with reflection faces present.
- One restricted Hurwitz step on its interior horizontal edge, with
  $u(T_aX)-u(X)=-48$. Every individual face charge is unchanged.
- A complete 67-edge read set for this response difference, contained in
  the disk. Summing the whole raw operator bank independently gives the
  same $-48$.

The calculation exports the raw disk links, the rooted patch, and the
seed. Full local holonomy image and a reflection face imply
$\mathcal O_D$ is contained in the required global nonabelian-reflection
sector. Hence its probability under that conditioned reference is at
least its all-sector probability. No assumption about sector mixing is
required.

Translate the witness to all cells whose disks do not cross a periodic
seam. There are $F/2-O(\sqrt F)$ distinct selected elastic slots.
Since

$$S=\frac12\sum_{a\in E}
\mathbb E[(u(T_aX)-u(X))^2],$$

nonnegativity and linearity of expectation give

$$\boxed{\liminf_{F\to\infty}\frac SF
\geq s_*(\rho):=576\frac{z_\rho^{60}}{p_0(z_\rho)^{60}}>0.}$$

Overlaps among the translated disks do not affect this inequality.
At density one, $z_\rho=1/\sqrt2$ and

$$s_*(1)=576(3-2\sqrt2)^{60}\approx6.72\times10^{-44}.$$

The limiting disk-orbit probability is approximately
$1.16659\times10^{-46}$. Its exact finite all-sector values at density
one are approximately $1.59925,1.31099,1.24157,1.19820$ times $10^{-46}$
for $F=128,288,512,1152$, respectively.

## A volume-independent cost bound

The preceding [locality proof](triangle-gauge-current-corrector.md) gives
$T\leq12C_1C_2S$, where $T=\langle ELb_3,L ELb_3\rangle$.
To avoid extrapolating a measured finite-size stencil count, use the
following looser, uniform geometric bounds.

The read set of $h_e=E_eu$ lies in boundaries of faces at dual distance
at most four from either face incident to $e$: a fan containing $e$
extends at most two steps, and a current term sharing one of its faces
extends at most two more. A degree-three dual ball of radius four has
at most $1+3+6+12+24=46$ faces. Thus each $h_e$ reads at most
$2\cdot46\cdot3=276$ edges.

By the symmetric distance condition, a single written edge can touch at
most 276 such $h_e$ terms. Each primitive writes at most two edges, giving
$C_1\leq552$. A given primal edge is written by two rooted pair patches
with three rules each, and by four three-face paths with twelve rules
each: $6+48=54$ slots. Therefore $C_2\leq276\cdot54=14904$, and

$$\boxed{T\leq C_*S,\qquad C_*=98\,724\,096.}$$

This uses only the specified triangular-torus incidence and local bank;
it is not a statement about arbitrary dynamically changing hypergraphs.

## Strict improvement of the bulk trial bound

Let $v(\rho)=\lim V/F$, $w(\rho)=\lim W/F$, and $U_1(\rho)$ be the
already derived [one-function bulk bound](triangle-winding-transport.md).
Its $v,w$ are finite and strictly positive at every interior density.
Use $f=(V/W)(b_3-h/C_*)$, with no fitted coefficient. The finite-volume
variational inequality and the lower bound on $S/F$ imply

$$\boxed{\limsup_{F\to\infty}\sigma_F
\leq U_1(\rho)-\epsilon_*(\rho),\qquad
\epsilon_*(\rho)=\frac{2}{810C_*}
\left(\frac{v(\rho)}{w(\rho)}\right)^2s_*(\rho)>0.}$$

At density one, $\epsilon_*\approx1.07\times10^{-57}$. This proves
strictness, not a useful numerical magnitude; it cannot change the
reported digits of the previous bound. The gain is over the specific
charge-drift trial, not the optimal space of all charge-only correctors.
It does not prove that the mobility limit exists, that an individual
trajectory samples the reference, or that this process has a hydrodynamic
diffusion limit. No quantum or molecular behavior follows from this result.

Calculation: [triangle_gauge_bulk_witness.py](../tools/triangle_gauge_bulk_witness.py).
