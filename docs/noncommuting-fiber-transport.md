# Noncommuting transport from a derived triangle fiber

Two closed transport circuits, performed in opposite orders, can return every defect
to its original face and preserve every face-curvature class, yet produce **different
gauge connections and different full bundle-graph spectra**. The existing vacancy rule
now realizes this on actual links with a triangle fiber. Its automorphism group is
derived from adjacency, $\operatorname{Aut}(C_3)\cong S_3$; no braid map is installed as
a microscopic update.

The two circuits reach exactly four complete gauge orbits. Their induced permutations
generate $A_4$. A fixed based-loop observer intertwines every measured action edge with
an independently evaluated Hurwitz word, and an exact finite-field calculation explains
the four-state action as $\operatorname{PSL}(2,\mathbb F_3)\cong A_4$.

This is a checked realization of **classical** nonabelian holonomy transport, not a new
discovery of braid groups, finite-group gauge theory, or quantum statistics. The physical
comparison is flux metamorphosis, described in
[de Wild Propitius, section 1.4.2](https://staff.science.uva.nl/f.a.bais/scripties/mdwp1995.pdf).
There are no quantum amplitudes, physical Hamiltonian, self-organized bound states, or
emergent geometry in this experiment. The lattice, fiber, seed, and circuit schedules
are supplied. The microscopic gauge group is $S_3$; $A_4$ is an **action on gauge states**,
not an additional internal symmetry imposed on each fiber.

![Measured primitive transport paths and their complete gauge-state action](images/fiber-transport-order.png)

The left panel uses supplied mesh indices, with periodic seams omitted. The right panel
is a layout of the measured state-transition graph, not a spatial projection of matter.

## 1. One rule, derived groups, actual connections

The generic runner enumerates fiber-graph automorphisms, compiles their products and
inverses, and uses the same vacancy rule for every derived group:

$$T(a,b)=\begin{cases}(b,a),&\text{exactly one of }a,b\text{ is }1,\\
(a,b),&\text{otherwise}.\end{cases}$$

The pair is based using explicit connector transport. The existing boundary-fixed
shared-edge lift realizes $T$ by changing only the common link of two faces. It preserves
the pair's ordered product, leaves all boundary links and spectator face holonomies
unchanged, and is gauge covariant and exactly reversible. The runner does not use the
abstract Hurwitz formulas below to evolve the links.

On a $12\times12$ triangular torus, a two-link noncommuting-reflection seed creates four
defects. Thirty-four vacancy moves park them at faces $(52,62,172,208)$. Circuit $A$ moves
the first around the second; $B$ moves it around the third. Each program has a lead path,
twelve clockwise steps around its target, and the reversed lead: 32 primitive updates.
Integer winding predicates check winding $-1$ around the intended defect and zero around
the other stationary defects. Routing depends on occupancy and mesh geometry, not on the
desired holonomy response.

We replay preparation, both circuit orders, and every edge in the reachable component
through independent C++ and Python link implementations. The entire schedule reverses
to the exact raw initial links. Gauge comparisons separately use a full spanning-forest
signature, a legacy quotient, and explicit vertex-frame witnesses.

## 2. Four states, not merely four labels

Breadth-first closure under $A$ and $B$ stores actual reachable raw connections and
deduplicates them by their **complete** gauge signatures. No canonicalized links are
injected into evolution. A budget cutoff raises an error rather than declaring closure.
In discovery order $0,1,2,3$, the measured permutations are

$$A=(1\;2\;3),\qquad B=(0\;1\;3).$$

Both have order three. Their generated permutation group is exactly the twelve even
permutations on four states; their commutator has order two. Starting at state zero,
$A$ then $B$ reaches state one, whereas $B$ then $A$ reaches state two.

The same preparation paths and circuit schedules run in two controls:

- A square fiber, with derived group $D_4$, reaches two states; $A$ is identity and $B$
  swaps them. The two orders are gauge-equivalent and have identical spectra.
- A triangle fiber with commuting seed reflections reaches one state. Both actions are
  identity, despite having the same defect positions and individual curvature classes
  as the noncommuting triangle case.

This closure concerns these **two controlled programs**. It does not exhaust all
microscopic updates, all torus connections, or all possible braid paths.

## 3. The raw action equals a fixed Hurwitz action

Transport each face loop to the fixed spanning-tree root before comparing holonomies.
The chosen perimeter order is lower-left, lower-right, upper-right, upper-left:

$$\mathbf H=(H_{52},H_{62},H_{208},H_{172}).$$

Every observed representative satisfies $H_{52}H_{62}H_{208}H_{172}=1$. Its simultaneous
conjugacy class distinguishes all states in the controlled component. Neither property
is assumed for arbitrary based loops or arbitrary connections on the torus. In particular,
simply sorting face IDs is **not** this ordered-loop convention.

Use the adjacent Hurwitz action

$$\sigma_i(\ldots,a,b,\ldots)=(\ldots,aba^{-1},a,\ldots),$$
$$\sigma_i^{-1}(\ldots,a,b,\ldots)=(\ldots,b,b^{-1}ab,\ldots).$$

For pair $(i,j)$, first execute $\sigma_{j-1},\ldots,\sigma_{i+1}$; execute
$\sigma_i^{-1}$ twice; then undo the preparation. Indices are zero-based, and the
sequence is written in execution order. This defines the inverse pure word $P_{ij}^{-1}$.
The clockwise raw circuits realize $A=P_{01}^{-1}$ and $B=P_{03}^{-1}$.

On **every** action edge, the predicted tuple equals the tuple read from the evolved raw
links exactly, without an extra fitted root conjugation. The target stored representative
can differ by a gauge frame; its tuple is simultaneously conjugate as required. This
intertwining check holds on all eight triangle edges, four square-control edges, and
two commuting-control edges. It is stronger than identifying an abstract group of the
same order, but is still a finite-component result, not a global isotopy theorem.

## 4. Why the triangle produces a projective line

Label triangle vertices by $\mathbb F_3$. Every derived reflection is
$h_x:y\mapsto-y+x$. Four such reflections have identity product precisely when

$$h_a h_b h_c h_d=1\quad\Longleftrightarrow\quad a-b+c-d=0\pmod3.$$

Conjugation by the affine permutation $y\mapsto sy+t$, $s\in\{1,-1\}$, sends
$x\mapsto sx+2t$. Thus common conjugation removes $a$ and leaves

$$u=b-a,\qquad v=c-a,\qquad d-a=v-u,$$
$$ (u,v)\sim(-u,-v). $$

There are $3^3=27$ neutral raw tuples. The three homogeneous tuples form one commuting
orbit, represented by $(0,0)$. The remaining 24 tuples give four gauge orbits, each of
size six, indexed by the nonzero projective points

$$\mathbb P^1(\mathbb F_3)=\{[0:1],[1:0],[1:1],[1:2]\}.$$

Direct substitution of the inverse pure words yields

$$A:\binom uv\mapsto\begin{pmatrix}1&0\\1&1\end{pmatrix}\binom uv,
\qquad B:\binom uv\mapsto\begin{pmatrix}2&2\\1&0\end{pmatrix}\binom uv\pmod3.$$

The code checks these formulas on **all 27 raw neutral tuples** and checks invariance
under every common gauge frame. Enumerating the matrix closure gives all 24 determinant-one
matrices. The kernel of their action on the four projective points is exactly $\{I,-I\}$;
the twelve resulting permutations are precisely $A_4$. The four measured connection
states map bijectively onto those four projective points, with the measured circuit action
matching the matrices. No field-valued coordinates or matrices enter the microscopic rule:
these are derived coordinates for its observed classical response.

## 5. Local deformations and an exact graph consequence

For each circuit, replace one step on its winding portion by the other five edges around
an adjacent empty primal-vertex star. On every state in each control's component, the
deformed circuit has the same target gauge orbit. Explicit frame witnesses locate the
difference at that one primal vertex. That is 14 independently replayed checks, covering
two chosen deformations across all seven component states. It does not claim that every
possible path deformation was tested. The local justification is the
[flat-star argument](transport-braid-memory.md#3-the-effect-survives-checked-path-deformations).

For the two triangle outcomes, the common-root relative word $H_{52}H_{172}^{-1}$ lies
in the nontrivial three-cycle class after $A$ then $B$, and is identity after $B$ then $A$.
This distinguishes their hidden relation while all individual face classes coincide.

Construct the full **unweighted 432-vertex bundle graph**, including vertical triangle
edges and horizontal permutation edges. Its integer Laplacians $L_{AB},L_{BA}$ have
different exact characteristic polynomials. Their first distinct trace power is twelve:

$$\operatorname{tr}L_{AB}^{12}=547642051586352,\qquad
\operatorname{tr}L_{BA}^{12}=547642051586712.$$

The difference is 360. Exact FLINT characteristic polynomials and Newton identities
give the result; a separate signed-integer matrix-power calculation verifies all moments
through twelve, with the explicit overflow bound $432\cdot16^{12}<2^{63}$. Both controls
have equal characteristic polynomials between orders. This is a consequence for the actual
graph, not just a selected fiber-sector operator or a frame-dependent visualization.
It remains a graph diagnostic, not an assumed physical wave frequency or energy.

## Reproduce and limitations

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j 4
export OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
uv run --with numpy --with scipy --with python-flint python tests/fiber_transport.py
uv run --with numpy --with scipy --with python-flint python tools/fiber_transport.py \
  --output out/fiber-transport.json
diff -u data/fiber-transport-order.json out/fiber-transport.json
uv run --with numpy --with scipy --with matplotlib python tools/plot_fiber_transport.py
```

The runner is deliberately bounded to sides 3–24, fibers with 1–6 vertices and at most
24 derived automorphisms, and 100,000 scheduled attempts per invocation. The four-defect
preparation requires side at least twelve. These are verification-runner limits, not
theoretical restrictions on fiber dynamics. Regression tests also cover singleton, path,
pentagon, and complete-four-vertex fibers, including non-square transport and local frames.

The next physical gap is not a prettier braid picture: it is whether unconstrained local
evolution creates useful, persistent relational structure and couples it to localization
or other graph-native observables. Imposing these circuits as autonomous rules, assigning
complex amplitudes by hand, or calling spectral memory a bound particle would not close it.
