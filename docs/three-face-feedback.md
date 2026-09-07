# Conserved transport with relative-loop feedback

A three-face neighborhood now escapes both earlier limitations: its conserved face densities
can move, and their evolution is not an autonomous function of individual face classes. The
counterexample keeps **every exterior link identical**, not just the boundary's conjugacy class.
It is realized on actual shared links and verified in a larger mesh.

The law is selected from a finite symmetry search; its conservation equations are solved afterward.
No target energy, continuum force, particle template, or molecular geometry enters that search.
The finite square fiber and oriented base triangulations remain supplied model choices.

![Actual fixed-exterior witness before and after the same update](images/three-face-feedback.png)

Colors indicate gauge-invariant face classes. The letter labels distinguish holonomies transported
to one common frame; changing that frame conjugates them together. The pictures are actual selected
fan states from the link simulation. Their coordinates come from the supplied triangulation.

## 1. Exhaustive minimal three-input closures

For $G=\mathrm{Aut}(C_4)$, enumerate proposed transpositions of triples in $G^3$ with the same
ordered product. Close each under simultaneous conjugation and orientation inversion

$$J(a,b,c)=(c^{-1},b^{-1},a^{-1}).$$

Reject overlapping unequal transpositions and any move of $(1,1,1)$. The resulting map is identity
off that closure. C++ breadth-first closure and an independent Python finite-action enumeration
agree on **945 nonidentity minimal closures**. Each map is rechecked for involution, ordered-product
preservation, local-frame covariance, and orientation inversion. This does **not** enumerate all
disjoint unions of closures, all reversible triple maps, or any three-input analogue of a braid law.

For every map, solve the complete rational equations

$$q(a)+q(b)+q(c)=q(a')+q(b')+q(c'),\qquad q(1)=0,$$

over class functions. Of the 945 minimal maps, 561 transport their maximal additive charge vector
without an autonomous charge-tuple description. Among these, 453 still have nonclosure after
conditioning on the exact ordered boundary product. **354** have different outgoing charge tuples
even when the incoming individual face classes and exact boundary product both match.

We select minimum moved-state count among those 354, then lexicographic transpositions to break
ties. This is a reproducible experimental choice, not a uniqueness claim about physical laws.
The selected rule moves eight triples. Its conservation equations give

$$q(z)=2q(r),$$

with independent weights on the other reflection class and the rotation class. Here $r$ is a
vertex-axis reflection, $z$ the central half-turn, and $r'=rz$ the other element in the same
reflection conjugacy class. The complete nontrivial exchanges are

$$
\begin{aligned}
(1,r,z)&\leftrightarrow(r,r,r'), & (1,r',z)&\leftrightarrow(r',r',r),\\
(r,r',r')&\leftrightarrow(z,r',1), & (r',r,r)&\leftrightarrow(z,r,1).
\end{aligned}
$$

The middle holonomy is unchanged. At the level of face populations, two unit-weight reflection
states can convert reversibly to one double-weight central state and one flat state, conditioned
on the neighboring middle reflection. The relative alignments determine whether the move occurs.
The rule's active triples lie in a Klein-four subgroup of $G$; this is not evidence that a fully
nonabelian continuum interaction has emerged.

## 2. Unique boundary-fixed realization on a three-face disk

Use actual oriented triangles

$$H_1=(u,v_0,v_1,u),\quad H_2=(u,v_1,v_2,u),\quad H_3=(u,v_2,v_3,u),$$

with five distinct vertices. Write $S_i=U_{u v_{i-1}}$ and $R_i=U_{v_{i-1}v_i}$. Then

$$H_i=S_{i+1}^{-1}R_iS_i,\qquad
H_3H_2H_1=S_4^{-1}R_3R_2R_1S_1.$$

The rule input is $(H_3,H_2,H_1)$ in that order. All exterior links, including $S_1,S_4$, stay fixed.
For the requested target holonomies, solve only the two internal spokes:

$$S_2'=R_1S_1(H_1')^{-1},\qquad S_3'=R_2S_2'(H_2')^{-1}.$$

Product preservation forces the third target. The solution is unique for those boundary links
and holonomy targets. Independent vertex-frame covariance follows from the endpoint transformation
of each factor. Inverting the triple map restores all holonomies, and uniqueness restores both
raw internal links.

Each written spoke must have exactly two incident faces, both in the fan. Therefore the update
reads every face it affects. The constructor rejects extra incident faces, non-cell loops,
incompatible orientation, and fans that close into a three-vertex link rather than a disk.
Three consecutive faces around a degree-six vertex in the periodic triangulation are valid;
the updates are not isolated degree-three compartments.

The compiled kernel uses the spoke recurrence. The independent connection and Python oracles
instead multiply unchanged exterior-prefix paths and cumulative target holonomies. These are
separate algebraic implementations of the same boundary condition.

## 3. Fixed-exterior feedback witness

In the common root frame, compare

$$x=(r,r,r'),\qquad y=(r,r',r).$$

They have the same ordered product $r'$, the same individual class tuple $([r],[r],[r])$, and
the same initial charge tuple. The selected rule gives

$$F(x)=(1,r,z),\qquad F(y)=y.$$

For the conserved weight satisfying $q(r)=1$, $q(z)=2$, the outgoing tuples are $(0,1,2)$
and $(1,1,1)$. Thus relative loop information controls the movement of a conserved density.

The link witness sets the outer chain transports to identity and the final exterior spoke to
$(r')^{-1}$. Only the two internal spokes encode the differing inputs. Embedded in a torus, the
two connections have identical **global** initial face classes and identical links everywhere
outside the two-spoke write set. A spanning-tree frame reconstruction exhausts every possible
root frame and proves that they are not gauge copies on the fixed labeled base. No quotient by
base-vertex permutations or an unmarked choice of update event is being asserted here.

## 4. Transport control, positivity, and sublattice constraints

The minimal collision law alone does not move an isolated defect. We separately implement the
vacuum-transport control

$$(g,1,1)\leftrightarrow(1,1,g),\qquad g\ne1.$$

Its fourteen moved states are disjoint from the eight collision states. Their union is therefore
their composition and remains an equivariant, orientation-compatible involution. All constraints
and charges are solved again for the combined 22-state-support table; the transport addition is
not silently included in the minimal census.

A basis of global conserved quantities is

$$Q_2=N_{[s]},\qquad Q_3=N_{[\rho]},\qquad Q_r=N_{[r]}+2N_{[z]},$$

where $[s]$ is the other reflection class and $[\rho]$ the rotation class. Their sum supplies a
positive conserved weight $W$: every nonflat face has weight at least one, so

$$N_{\rm nonflat}\le W.$$

This prevents an unbounded zero-weight central-curvature trail from a finite-$W$ seed. It does
not identify $W$ as physical energy, fix a mass scale, or single out one preferred positive
combination of the three charges. Flat connections can also carry nontrivial global holonomy.

There is a further restriction: the middle face is unchanged, and the outer faces of a fan lie
on the same dual sublattice when that graph is bipartite. Exact sublattice-dependent conservation
equations have a six-dimensional solution space after normalizing each sublattice's flat value to zero: the three charges
are separately conserved on each part. Independent replay checks those six totals too. This is
a constraint of the selected dynamics and topology, not evidence of emergent spin or dimensions.

## 5. Controlled evolution and retained negative results

The experiment uses the same initial raw connections and fixed schedule for three variants:
collision-only, transport-only, and combined. A common first fan update is an explicit witness
probe; subsequent updates cycle through seeded, conflict-free layers. The witness's first
conversion is therefore a controlled local test, not claimed to arise from a random encounter.

Two meshes are checked: 288 faces for 256 layers with eight schedules, and 1,152 faces for
512 layers with four independently seeded schedules. Eight initial conditions and three rule
variants give **288 runs**. All **15,309,480 forward updates** are independently replayed from
raw links in Python, checking every target triple, each layer's histogram, conversion/transport
event counts, final raw links, global charges, and sublattice charges. C++ separately reverses
every update and verifies every intermediate histogram and the initial links exactly.

The primitive test covers 3,072 exhaustive-triple/frame/control cases, including an attached
spectator face, plus overlapping torus updates, causal parents, commuting layers, and inversion.

In the combined rule's random-link runs, there are 114–367 conversion events per smaller-mesh
run and 1,532–2,189 per larger-mesh run. Final face classes differ from the matched transport-only
control at 101–213 of 288 faces and 804–829 of 1,152 faces, respectively. These are verified
finite-history responses, not measurements of a continuum scattering amplitude or mixing rate.

The dilute controls matter:

- Flat, single-reflection, single-rotation, and aligned-reflection-pair seeds have no conversion
  events. In their initial frame, subgroup closure explains the absence of the necessary inputs;
  combined and transport-only raw trajectories agree.
- Opposite-reflection pairs have the same initial face classes as aligned pairs and a larger
  generated subgroup, but **also have no conversion events in these runs**. Algebraic permission
  is not evidence of an encounter or a scattering process.
- The active fixed-exterior witness has 1–4 conversion events depending on size/schedule; the
  inactive witness has none. One smaller-mesh active trial even ends with the same face classes
  as its transport control despite its two intermediate conversions.

No stable bound state, calibrated energy, thermal ensemble, continuum limit, quantum amplitude,
or schedule-independent spacetime has been demonstrated. Read/write independence proves only
commutation of disjoint updates, not consistency of arbitrary overlapping schedules.

The subsequent [sparse-encounter audit](three-face-encounters.md) finds transport-generated
reactions under other schedules, proves an exact subgroup face factor, and completely classifies
two small charge sectors. It retains the original null results above and explains why the
selected rule's catalytic requirement strongly suppresses dilute conversion.

## Relation to the field and next discriminating tests

Gauge-invariant reversible automata already exist; for example,
[Arrighi, Di Molfetta, and Eon](https://arxiv.org/abs/1802.07644) give a procedure for gauging an
automaton. Finite-group Hamiltonian gauge theories also have established constructions, including
a group-Laplacian electric operator; see [Mariani, Pradhan, and Ercolessi](https://arxiv.org/abs/2301.12224).
Neither the existence of reversible gauge dynamics nor the interpretation of conserved observables
is claimed as new here. Literature priority for this particular finite construction is unestablished.

The concrete contribution is a reproducible minimal triple census, a boundary-fixed link lift,
and a same-exterior counterexample combining conserved transport with relative-loop feedback.
The next tests should measure which encounters are dynamically reachable, distinguish transient
conversion from persistent localized structure, and test sensitivity to overlapping schedules
and bipartite separation. Geometry-changing engine events and phase-bearing amplitudes remain
separate missing constructions; a force law should not be inserted to bypass them.

## Reproduce

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build -L wgphysics --output-on-failure
python3 tools/triple_rule_search.py --output out/triple-rule-search.json
diff -u data/d4-triple-rule-search.json out/triple-rule-search.json
python3 tests/three_face.py
python3 tools/run_three_face.py --full-replay --output out/three-face.json
diff -u data/d4-three-face.json out/three-face.json
python3 tools/run_three_face.py --side 24 --layers 512 --trials 4 --seed 1290831 --full-replay --output out/three-face-replication.json
diff -u data/d4-three-face-replication.json out/three-face-replication.json
uv run --with matplotlib python tools/plot_three_face.py
```

The datasets retain all minimal closures, charge bases and witnesses, the explicitly composed
tables, schedules, initial/final links, layer histories, and independent replay results.
