# Complete loop information and a graph-derived two-component mode probe

The [equilibrium class-binding restriction](reachable-gauge-equilibrium.md) concerns local
curvature classes. It does not say those classes determine the gauge connection. We now
have a complete, reconstructible observer for the square-fiber connection on a fixed
labeled graph, an actual conversion witness for the information classes omit, and a
graph-native spectral diagnostic that distinguishes the witness branches.

These additions do **not** install a Hamiltonian, a wave equation, complex amplitudes,
a matter field, or an attractive force in the state-update rules. The two-component
space below is derived from the fiber adjacency and used as an observer. Its selection
does not assert that the simulated dynamics favor it as a physical low-energy sector.

## 1. Derive the central extension from the fiber automorphisms

Enumerate $G=\operatorname{Aut}(C_4)$ from the square adjacency, then find its center and
two noncommuting involutions $r,s$. The center is $\{1,z\}$ and every element has a unique
normal form

$$g=z^c r^a s^b,\qquad a,b,c\in\mathbb F_2.$$

Both implementations verify the multiplication and frame-action formulas against every
entry of the derived permutation tables:

$$
(a,b,c)(d,f,h)=(a+d,b+f,c+h+bd),
$$

$$
(a,b,c)\longmapsto(a,b,c+ub+va)
$$

under conjugation by a frame with quotient coordinates $(u,v)$. Arithmetic is modulo two.
The quotient is $G/\langle z\rangle\cong C_2\times C_2$, but the central term $bd$ cannot
be discarded: it is the noncommuting part of the multiplication. Tests verify its cocycle
identity and exhaust all eight normalized changes of section, none of which remove it.
This is a finite-group central extension, not a declaration of $U(1)$ or a continuum gauge field.

## 2. Count and recover the gauge-invariant information

Consider $m$ loops based at the **same** vertex. Let $r_*$ be the rank over $\mathbb F_2$
of their quotient-coordinate vectors, so $r_*\in\{0,1,2\}$. Simultaneous frame changes
shift the vector of central bits through a rank-$r_*$ subspace. Choosing at most two
independent pivot loops fixes those pivot central bits to zero, leaving

$$2m+(m-r_*)=3m-r_*$$

binary positions in a direct encoding of the tuple modulo simultaneous conjugation.
At fixed quotient vectors there are exactly $2^{m-r_*}$ central possibilities; the
quotient vectors themselves obey rank constraints, so $2^{3m-r_*}$ is **not** an orbit
count for the rank stratum. Nor are the current output arrays bit-packed: the C++ API
stores small group indices as `uint16_t`.

Individual loop classes already reveal the central bit of central holonomies. If $n$
loops are noncentral, precisely $n-r_*$ further binary relations remain hidden. A complete
set of probes uses:

1. The derived-group conjugacy class of every individual loop.
2. For each noncentral quotient type, pick one representative loop. For every other loop
   of that type, measure the central holonomy $g_i g_{\rm rep}^{-1}$.
3. If all three nonzero quotient types occur, add one central word $g_i g_j g_k^{-1}$
   using their representatives.

There are exactly $n-r_*$ relation bits. Reconstruction first chooses two representative
central bits as frame conventions, solves the cross-type word for the third representative
when necessary, and propagates the same-type differences. It recovers a representative
of the **entire** simultaneous-conjugation orbit, not just a collection of marginal labels.

Pairwise probes alone are insufficient. The triples $(1,2,3)$ and $(1,2,6)$ in the
exported group indexing have identical individual classes and identical classes of every
pair product, but opposite values of the three-type central relation. Exact enumeration
and reconstruction cover all raw tuples of zero through four loops. The resulting orbit
counts are $1,5,28,176,1216$, matching the independent Burnside count

$$\#(G^m/G)=\frac{8^m+3\cdot4^m}{4}.$$

Closed holonomy observables and their completeness are established topics; see, for example,
[Lévy, Wilson loops in the light of spin networks](https://arxiv.org/abs/math-ph/0306059).
That paper's compact-group results are not invoked as a proof for this finite specialization.
Our completeness argument is the explicit central-extension reconstruction above. We use
full classes in $\operatorname{Aut}(C_4)$, not symmetric-group cycle type or one chosen trace
that can merge distinct classes. No general novelty or priority claim is made.

## 3. Make the based loops from actual connections in a linear pass

Compile a deterministic spanning forest of the base graph. Within each component, compute
the tree transport $T_v$ from its root to each vertex once. For every non-tree edge $u\to v$,

$$H_{uv}=T_v^{-1}U_{uv}T_u$$

is the corresponding root-based loop holonomy. Normalize the resulting tuple with the
derived central coordinates. There are $E-V+C$ such loops for $C$ connected components,
and the probe uses $O(V+E)$ work after forest compilation; it never constructs every full
walk separately. Setting tree links to identity and chord links to the normalized tuple
reconstructs a gauge-equivalent connection.

Python and C++ independently implement this procedure. Tests compare actual fundamental
and relation walks with raw edge transports, change frames independently at every vertex,
include disconnected graphs and isolated vertices, reconstruct raw representatives, and
compare against the original general-purpose gauge quotient.

The forest and loop basis are coordinate conventions. This is **not** a graph-isomorphism
quotient, an embedding, a dimension estimate, or a time coordinate. Comparisons across
base rewrites need explicit loop transport or recompilation and a basis map; comparing
the coordinate arrays from different graphs directly would be invalid.

## 4. A microscopic witness with identical face classes but different responses

From the existing two-link seed on the 18-face torus, the same fan input $(1,0,2)$ admits
two bank operators giving

$$
(1,0,2)\longmapsto(2,2,6),\qquad
(1,0,2)\longmapsto(2,7,3).
$$

The labels are raw derived-group indices; the identity is 0 and the central involution is 5.
Both outputs have face classes $(2,2,3)$, the exact exterior holonomy is 6 in both cases,
and every exterior link is unchanged. But the same-type two-loop word is identity on
one branch and $z$ on the other. The full connections are not gauge equivalent.

The branches come from registered rules 33 and 34 in the first bank, corresponding to
minimal census rules 364 and 366. They are replayed by the existing mixed-arity C++ engine,
independently reconstructed on Python links, and reversed exactly. The new C++ observer
and the old gauge quotient both distinguish them.

Applying rule 33 again restores the first branch's initial state but leaves the second
branch's face classes different. Thus the hidden relative bit matters for a **specified
microscopic rule**. Averaging uniformly over all 48 rules gives the same next-class
distribution on both branches: 36 identity-class outcomes and two of each of the six
ordered permutations of $(0,1,2)$. This preserves, rather than evades, the previously
proved class-Markov reduction. The branch alternatives are classical choices, not
coherent quantum amplitudes.

## 5. The hidden distinction changes a graph-native mode spectrum

The square adjacency has a two-dimensional zero eigenspace. Rational row reduction gives
a basis matrix $Q$ with $Q^TQ=2I$. Every fiber permutation preserves this space, and its
induced action is derived as

$$R_g=\tfrac12 Q^T P_g Q,\qquad P_gQ=QR_g.$$

The resulting real matrices include

$$R_r=\begin{pmatrix}1&0\\0&-1\end{pmatrix},\qquad
  R_s=\begin{pmatrix}0&1\\1&0\end{pmatrix},\qquad R_z=-I.$$

They obey $R_rR_s=-R_sR_r$ and all 64 group multiplication identities. No common real
line is invariant under both reflection matrices, so this is the nontrivial irreducible
two-component action. These familiar matrix identities do **not** establish spin one-half,
fermionic statistics, a qubit, the Born rule, or an emergent quantum theory.

Construct the actual unweighted total bundle graph: each base vertex carries the square
fiber, and a base edge connects fiber vertex $i$ to $U_{uv}(i)$ in the neighboring fiber.
Restrict its combinatorial Laplacian to the derived two-component space at each vertex.
The implementation checks the exact integer intertwining relation

$$L_{\rm bundle}\mathcal Q=\mathcal Q(L_U+2I),$$

where the vertical shift 2 is the square **Laplacian** eigenvalue of this adjacency-zero
subspace. The horizontal block has

$$\psi^T L_U\psi=\sum_{u\to v}\|\psi_v-R_{U_{uv}}\psi_u\|^2.$$

This identity is verified from graph incidence, alongside symmetry, local-frame covariance,
and the constant-section zero modes of a flat connection. It is a spectral diagnostic of
the existing graph; the quadratic form is **not assigned physical energy** and does not
govern any simulation update.

The two actual conversion branches have equal $\operatorname{tr}(L_U^k)$ for $k=1,2,3$,
but at the fourth power give **43,444** and **43,284**. The full 36-vertex bundle-graph
Laplacians likewise give different fourth moments, **250,968** and **250,808**. These
exact integers prove that the branches have different mode spectra and different full
bundle-graph spectra, despite identical individual face classes and exterior data.
No eigensolver tolerance or fitted orbital shape is involved.

This supplies a concrete bridge from hidden gauge relations to graph mode structure.
It does not yet demonstrate localized modes, persistence under evolution, propagating
matter, binding, or a physical wave equation. Those remain tests to perform, not labels
to assign to the eigenvectors.

## 6. Measured CPU scaling, with the comparison boundaries retained

The checked Apple M1 Max benchmark uses deterministic random connections and a Release
build. Median warm probe timings are:

| Edges | Independent loops | Complete probe | Forest compilation |
| ---: | ---: | ---: | ---: |
| 768 | 513 | 0.002834 ms | 0.128 ms |
| 3,072 | 2,049 | 0.009375 ms | 0.479 ms |
| 49,152 | 32,769 | 0.144208 ms | 7.50 ms |
| 1,080,000 | 720,001 | 3.44362 ms | 199 ms |

The older quotient took median 3.27 and 26.89 ms in the two small cases. This is an
implementation comparison, not a pure algebra speedup: the new path caches topology,
avoids repeatedly scanning all edges to find neighbors, avoids rebuilding full walks,
and emits compact group-index representatives instead of permutation vectors. Probe
times include output allocation, but exclude graph/forest construction and the subsequent
validation comparisons. The largest cases are checked under local frame changes and
representative reconstruction; the legacy quotient was **not** run at that size.

The measured result is CPU **observation** of a fixed topology—not million-edge evolution,
dynamic forest maintenance, a renderer, or GPU performance. The small exact algebra could
support parallel prefix/tree kernels later; no such GPU implementation is claimed here.

## Reproduce

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j 4
python3 tests/d4_loop_observer.py
python3 tests/fiber_mode_probe.py
python3 tools/d4_loop_observer.py --output out/d4-loop-observer.json
diff -u data/d4-loop-observer.json out/d4-loop-observer.json
python3 tools/fiber_mode_probe.py --output out/fiber-mode-probe.json
diff -u data/d4-fiber-mode-probe.json out/fiber-mode-probe.json
python3 tools/run_loop_bench.py --output out/loop-observer-cpu.json
```

Exact observer and mode artifacts are regenerated in Linux/macOS CI. Performance results
are measurements with individual samples retained, not byte-stable outputs or CI thresholds.
