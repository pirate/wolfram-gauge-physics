# Infrageometric gauge foundations

This layer follows the finite combinatorial hierarchy implemented by the Wolfram Institute's
experimental `InfraGaugeTheory` project, rather than naming a continuum gauge group at the
microscopic level. The conceptual dependency is pinned for this audit to repository commit
`edd9bdca46b7838d6b3e940e8ae8cde90b60ef2c` (paclet version 1.0.7). The implementation here is
an independent C++ adaptation for the hypergraph rewrite engine; no Wolfram Language source is
vendored.

## Primitive hierarchy

1. A base graph `B` represents one discrete spatial state.
2. A finite graph `F_x` lies over every base vertex `x`; together these form the total space `E`.
3. Projection `pi: E -> B` says which total-space vertices belong to each fiber.
4. A connection supplies unique horizontal lifts over base edges.
5. Lifting a path composes fiber isomorphisms and defines parallel transport.
6. Lifting a loop gives a fiber automorphism: holonomy, the discrete curvature observable.

The current C++ representation specializes to isomorphic copies of a single microscopic fiber
graph `F`. This is the complete-lift case of `InfraGaugeTheory`'s graph connections. General
non-isomorphic and partial fibers can be added after rewrite behavior is stable.

## Gauge group is derived

The allowed change of local fiber frame is not supplied by a string such as `U(1)`. It is

`G = Aut(F)`,

computed by exhaustive graph-automorphism enumeration. A connection map `U_xy: F_x -> F_y`
changes under local frames `g_x, g_y` as

`U_xy -> g_y U_xy g_x^-1`.

Consequently loop holonomy changes only by conjugation at its basepoint. Its conjugacy class,
represented initially by the permutation cycle signature, is gauge invariant. Different fiber
graphs therefore produce different microscopic symmetry groups. Whether an effective continuous
group appears after coarse-graining is a result to measure, not an input.

## Coupling to base rewrites

For the elementary spatial rewrite

`(x--y) -> (x--w--y)`,

the old transport is factored through the fresh fiber:

`U_xy = U_wy U_xw`.

There is one factorization for every element of `Aut(F)`, but all factorizations form one orbit
under changes of frame at `w`. They must therefore not be counted as different physical multiway
branches. The quantum rewrite is the normalized state over this orbit. Distinct boundary
holonomies remain orthogonal.

## Steps beyond the upstream building blocks

The first implemented extensions are:

- exact derivation and a noncommutativity test for the local automorphism group;
- gauge transformations of entire connections and conjugacy-class observables;
- a rewrite extension that preserves boundary parallel transport;
- proof that fresh-fiber connection factorizations form a single gauge orbit;
- a normalized quantum amplitude over that orbit;
- orthogonality of distinct boundary-holonomy sectors.

An additional exact quotient removes local frame redundancy before general graph
canonicalization. A deterministic spanning forest sets every tree transport to identity. Chord
transports become fundamental-cycle holonomies, and simultaneous conjugation at each component
root is canonicalized exhaustively over `Aut(F)`. This turns the apparent local-frame space from
`|Aut(F)|^|V|` assignments into a compact signature containing one group element per independent
cycle.

For execution, `Aut(F)` is also compiled into dense `uint16` multiplication, inverse, and action
tables. The derivation of a fiber's automorphisms remains an exact CPU-side setup operation; the
repeated evolution operations become small integer lookups suitable for device constant or shared
memory.

This still defines kinematics, curvature sectors, and rewrite amplitudes—not matter or a force
law. The next physical search layer must construct local update operators from total-space graph
structure and causal rewrite data, then test which fiber families have propagating curvature
defects and stable localized sectors.
