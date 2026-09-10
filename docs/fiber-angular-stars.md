# Exact angular orbits and the fast-angular effective generator

The angular generator on the actual shared-link mesh factorizes into
independent constrained walks, conditional on its frozen data. This gives
an exact description of its kernel, a volume-independent conditional
spectral-gap bound, and a computable angular-averaging projection.

This advances the [long-time transport question](fiber-bulk-transport.md):
we can now specify precisely which variables survive angular relaxation
and which effective generator must enter the current-corrector equation.
It does not yet solve that global equation or establish a physical force.

## Raw-link factorization, not independent face evolution

Write each actual $D_n$ link permutation in affine form

$$U_e(v)=s_ev+a_e\pmod n,\qquad s_e\in\{+1,-1\}.$$

An angular move preserves every face charge and every link sign $s_e$.
It changes only the shift $a_e$ of a shared edge between one reflection
face and one nonidentity rotation face. Call these edges active. Every
other link remains fixed under angular evolution.

Each active edge belongs to exactly one rotation face. For a rotation
face $z$, let its active edges be $e_1,\ldots,e_d$, where $1\le d\le3$.
Its actual based rotation label is

$$z=c+\sum_{i=1}^d\epsilon_i a_{e_i}\pmod n,\qquad\epsilon_i\in\{+1,-1\}.$$

The coefficients follow from oriented link composition with the fixed
signs; $c$ uses the inactive links. The only restriction is $z\ne0$.
Reflection-face holonomies have arbitrary axes and impose no extra
restriction on these shifts. Identity faces have no active edges.

Therefore the shift variables partition into disjoint rotation-face
stars. Stars may share a reflection neighbor, but do not share an active
link or a nonzero-rotation constraint. Their generators commute and add.
This is a factorization of actual raw-link variables, retaining all
boundary transports, not a product assumption about evolving faces.

After replacing $a_{e_i}$ by $x_i=\epsilon_i a_{e_i}$, a degree-$d$ star is

$$\Omega_{n,d}=\{x\in(\mathbb Z/n\mathbb Z)^d:
c+\sum_i x_i\ne0\}.$$

Each allowed $x_i\to x_i\pm1$ has rate two: the two rooted pair
descriptions of that primal edge give the same raw move. A move to zero
rotation is absent. Thus the local graph is a discrete torus with one
sum hyperplane removed, with

$$|\Omega_{n,d}|=n^{d-1}(n-1).$$

It is connected for every odd $n\ge3$. At fixed nonzero sum, transfers
$x_i\to x_i+1$, $x_j\to x_j-1$ can be ordered to avoid zero at the
intermediate step; at least one ordering works. These transfers connect
the entire fixed-sum layer. Single-coordinate steps connect its nonzero
sum labels along the path $1,\ldots,n-1$.

Consequently a complete angular component is the Cartesian product of
these stars. Its frozen labels are the charge field, all link signs,
and the full values of inactive links. They are framed coordinates, not
all independent physical observables. Gauge transformations map components
to components, and the averaging projection commutes with gauge changes.

The saved 128-face $C_5$ example has 17 degree-one, 13 degree-two, and
13 degree-three stars. Its component contains approximately
$1.407\times10^{53}$ raw states, but its uniform conditional measure is
sampled by drawing each star separately: draw $d-1$ shifts uniformly,
draw its nonzero sum uniformly, then solve for the last shift. This
sampler is an exact conditional projection, not a finite-rate replacement
for angular evolution.

## Fourier reduction and an explicit mixing bound

Put $N=n-1$ and use the sum $z\in\{1,\ldots,N\}$ plus $d-1$ tangent
coordinates. Fourier transform the tangent coordinates with
$k\in(\mathbb Z/n\mathbb Z)^{d-1}$. Define

$$S(k)=1+\sum_{i=1}^{d-1}e^{2\pi i k_i/n},\qquad\rho(k)=|S(k)|.$$

The phase of $S$ can be removed along the open $z$ path. Each block of
the positive generator is the real tridiagonal matrix

$$H_k=2dD_N-2\rho(k)A_N
=2\rho(k)L_N+2[d-\rho(k)]D_N,$$

where $A_N$ is path adjacency, $D_N$ has endpoint entries one and
interior entries two, and $L_N=D_N-A_N$. The zero-frequency block is
$2dL_N$, with first nonzero eigenvalue

$$\lambda_\parallel=4d[1-\cos(\pi/(n-1))].$$

For nonzero $k$, at least $d-1$ pairs of phases differ. Therefore

$$d-\rho(k)\ge\frac{d-1}{d}[1-\cos(2\pi/n)],$$

and $D_N\succeq I$ gives

$$\gamma_{n,d}\ge\min\left\{
4d[1-\cos(\pi/(n-1))],
\frac{2(d-1)}d[1-\cos(2\pi/n)]\right\},\quad d\ge2.$$

Degree one has exactly the longitudinal gap. Combining these bounds
for $d=1,2,3$ yields a simple bound for every angular component:

$$\boxed{B\succeq\frac8{n^2}(I-P),\qquad
\|e^{-\alpha tB}f-Pf\|_2\le e^{-8\alpha t/n^2}\|f-Pf\|_2,}$$

where $P$ averages on the component. A singleton component has $P=I$.
The number of base faces does not enter this conditional $L^2$ bound.
Whole-configuration total-variation mixing can still depend on volume,
and this is not a bound on relaxation of the full reacting model.

The spectra are reduced exactly, not fitted to a diffusion equation:

| Cycle size | $n^2\gamma_{n,1}$ | $n^2\gamma_{n,2}$ | $n^2\gamma_{n,3}$ |
| ---: | ---: | ---: | ---: |
| 3 | 36.000 | 18.000 | 22.823 |
| 5 | 29.289 | 28.086 | 37.078 |
| 81 | 20.233 | 38.975 | 51.966 |
| 243 | 19.902 | 39.314 | 52.419 |

The degree-three $n=243$ star has 14,289,858 states. Its gap comes from
a 242-dimensional tridiagonal block after maximizing the transverse
hopping over the finite Fourier labels. Full small-star graph spectra
agree with this reduction for $n=3,5,7$.

These are raw spectra and include gauge-frame directions. A genuinely
gauge-invariant angular observable for $n\ge5$ is the even path mode

$$f(z)=\cos\!\left(\frac{2\pi(z-1/2)}{n-1}\right).$$

It is invariant under $z\to n-z$ and has eigenvalue
$4d[1-\cos(2\pi/(n-1))]$. Thus gauge-invariant internal relaxation can
also scale as $n^2/\alpha$. The tangent Fourier coordinates are internal
link variables, not newly detected spatial dimensions.

## Exact conditional reaction rates

Fast angular averaging does not make all states with the same charge
field equivalent. The frozen holonomies remain in conditional reaction
probabilities.

There is a useful exact local counting formula. For a uniform star,
$\omega=e^{2\pi i/n}$ and $v\in(\mathbb Z/n\mathbb Z)^d$,

$$\mathbb E[\omega^{v\cdot a}]=
\begin{cases}
1,&v=0,\\
-\omega^{-tc}/(n-1),&v=t\epsilon,\ t\ne0,\\
0,&\text{otherwise}.
\end{cases}$$

Products of these factors evaluate linear holonomy constraints without
discarding the independent frozen data. For prime $n$, the implemented
alternative uses inclusion–exclusion over the touched stars' forbidden
hyperplanes and exact modular Gaussian elimination. Each consistent
affine system of rank $r$ in $D$ variables has $n^{D-r}$ solutions.
Only stars touching the actual based reflection loops enter the count.

If their reflection axes are $a,b,c$, averaged activity is

$$\bar a=\Pr(a=b)+\Pr(b=c)-2\Pr(a=b=c).$$

For the four-face family

$$w_k=(r_a,r_b,r_{k-a+b},t_k),\qquad 1\le k\le n-1,$$

the last pair forms one degree-one star while the first two reflections
stay frozen under $B$. All members have the same charge tuple and boundary
$r_0$. The four kinds of angular component give:

| Frozen anchors $(a,b)$ | Averaged triple activity | Forward reaction rate at $C_5$ |
| --- | ---: | ---: |
| $(0,1)$ | $0$ | $0$ |
| $(1,0)$ | $1/(n-1)$ | $3$ |
| $(1,1)$ | $(n-2)/(n-1)$ | $9$ |
| $(0,0)$ | $1$ | $12$ |

These are realized on actual periodic raw connections with identical
whole-mesh charge fields and one movable star. Including all mesh
reactions, their exact averaged rotation-count drifts at $C_5$ are
$0,-9,+9,0$, respectively. Hence even the infinitely fast-angular process
cannot close on the charge field alone. This is not just a distinction
between marked reaction channels.

## The effective generator and the remaining transport criterion

At fixed finite volume, the fast-angular limit of $A+\alpha B$ on
observables constant on angular components has positive generator

$$\boxed{H=PAP.}$$

Its transitions average old primitive events over the exact product-star
measure. An old event can change signs, charges, and the active-edge
partition, so the component labels evolve. Stars are independent only
during angular evolution with those labels held fixed. Reactions do not
evolve independently in the full model.

The local product description allows exact conditional sampling of the
stars needed by an old event, including stars whose edges will become
frozen or active. This supplies a route to simulating the limiting
generator without enumerating an entire angular component. Replacing
finite-rate evolution by such refreshes would be an approximation and
is not done in the existing bulk runs.

For the winding-current source $b$, $Pb=b$. Write $Q=I-P$ and define

$$C=QAP,\quad D=QAQ,\quad B_Q=QBQ.$$

After removing common stationary constants, block elimination gives

$$P(A+\alpha B)^{-1}P
=\left[H-C^\dagger(D+\alpha B_Q)^{-1}C\right]^{-1}.$$

The star gap ensures $\|B_Q^{-1}\|\le n^2/8$. Let $u_\infty=H^+b$ and
$r=Cu_\infty$. The current-corrector obstruction is now

$$\boxed{r=0\ \Longleftrightarrow\
\text{angular rate leaves the exact DC mobility unchanged}.}$$

For nonzero $r$, the approach to saturation at fixed finite volume is

$$\sigma_F(\infty)-\sigma_F(\alpha)
=\frac{\langle r,B_Q^{-1}r\rangle}{\alpha F}
+O(\alpha^{-2}).$$

The local inverse is controlled, but $u_\infty$ still solves the spatial
effective-generator problem. The conditional gap does not by itself
bound that corrector or justify exchanging the fast-angular and volume
limits. Determining $r$ for the winding source remains the unresolved
low-frequency step.

The new result is the exact angular component structure, its spectrum,
and the conditional rates feeding that effective problem—not a newly
postulated force law. All geometry and stochastic clocks remain supplied;
no quantum amplitude or molecular binding has been derived.

Apparatus: [fiber_angular_stars.py](../tools/fiber_angular_stars.py).
Raw components, conditional samples, spectra, and rational reaction rates:
[fiber-angular-stars.json](../data/fiber-angular-stars.json).

The [refinement-limit derivation](fiber-refinement-limit.md) identifies a
different obstruction: under diffuse raw-shift ensembles, reactions vanish
at fixed time uniformly in angular speed, and the charge process converges
to exclusion diffusion. Fast conditional angular mixing alone cannot
prevent this loss of interaction measure.
