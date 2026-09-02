# From gauge-covariant motion to physical curvature transport

This note records both a finite no-go result and the first constructive primitive beyond it. The
distinction matters: changing a connection table is not automatically a physical field evolution.

## Unary cell maps: an exact D4 no-go census

Let `G = Aut(F)` and let one oriented cell have based holonomy `H in G`. A microscopic unary update
`phi : G -> G` is admitted only when it is:

- bijective, so the time step is reversible;
- identity-fixing, so a flat cell cannot create curvature from nothing;
- inversion-equivariant, `phi(H^-1) = phi(H)^-1`, so reversing the loop is consistent;
- conjugation-equivariant, `phi(g H g^-1) = g phi(H) g^-1`, so it commutes with local frame changes.

The implementation validates all four conditions at the application boundary, not merely during
candidate generation. For the square fiber, `Aut(C4)=D4`, exhaustive search finds exactly eight
maps. Some change the based group element, including a quarter turn, but all eight induce the
identity map on the five conjugacy classes. Consequently they produce zero changes in the exact
gauge-quotiented curvature sector across the engine product evolution.

This is a scoped computational no-go result: for this fiber and this unary rule family, apparent
motion is entirely within gauge-equivalent representatives. It is not a theorem about every finite
group or every local interaction. The complete deterministic result is in
`data/d4-cell-dynamics-census.json`.

## Minimal two-cell transport: the Hurwitz move

Physical motion first appears when two adjacent oriented cells interact. With their loops based at
the same vertex and holonomies `(A,B)`, the implemented forward update is

`(A,B) -> (A B A^-1, A)`.

Its exact inverse is

`(C,D) -> (D, D^-1 C D)`.

This is the elementary Hurwitz or braid action. It has three properties needed for a useful
bottom-up primitive:

1. It is reversible.
2. It commutes with simultaneous conjugation at the common basepoint and its link realization
   commutes with arbitrary local fiber-frame changes.
3. It conserves ordered total holonomy because `(A B A^-1) A = A B`.

If `(A,B)=(R,1)` for a nontrivial quarter turn `R`, the result is `(1,R)`: curvature moves from one
cell to its neighbor without inserting a potential, continuum gauge group, or named force.

The current exact realization requires each cell loop to have an exclusive closing edge. It updates
those two connection transports and verifies the target holonomies and conservation law. This
requirement is explicit because inferring faces from a bare graph would be ambiguous.

`OrientedCellComplex` now makes those faces first-class cyclic boundary words. Its product identity
jointly quotients base relabelings, face basepoint rotations and ordering, and local fiber frames,
while retaining face orientation and incidence. Subdividing a shared base edge updates every
incident face boundary. The official-engine product runner now carries this object beside every raw
state and event and computes curvature on the explicit face provenance. The next engine-level step
is to introduce independent gauge events and derive their causal edges from the links and faces they
actually read and write.

## Claim boundary

The Hurwitz move is a mathematically natural, nontrivial transport primitive, not yet a selected law
of nature. A physical candidate still has to survive multiway causal alignment, orientation and
foliation tests, larger fiber/rule censuses, localization, scattering, and continuum-limit studies.
The unary no-go census is retained precisely so a gauge-representative change cannot be mistaken for
such evidence.
