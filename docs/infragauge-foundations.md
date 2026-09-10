# Finite-fiber construction and its assumptions

This is a supplied-fiber laboratory, not a derivation of internal fibers
from a bare hypergraph rewrite rule. The base, fiber adjacency, and
connection data are inputs. Computing their symmetries determines what
follows from those inputs; it does not explain why nature would select them.
The missing rewrite-derived construction is the first question in the
[research direction](research-direction.md).

This layer follows the finite combinatorial hierarchy implemented by the Wolfram Institute's
experimental [InfraGaugeTheory](https://github.com/WolframInstitute/InfraGaugeTheory) project,
rather than naming a continuum gauge group at the microscopic level. The initial conceptual
audit used [this source revision](https://github.com/WolframInstitute/InfraGaugeTheory/tree/edd9bdca46b7838d6b3e940e8ae8cde90b60ef2c)
(paclet version 1.0.7). The implementation here is
an independent C++ adaptation for the hypergraph rewrite engine; no Wolfram Language source is
vendored.

## Hierarchy within the supplied-fiber model

1. A base graph `B` is the supplied discrete base, not yet physical space.
2. A finite graph `F_x` lies over every base vertex `x`; together these form the total space `E`.
3. Projection `pi: E -> B` says which total-space vertices belong to each fiber.
4. A connection supplies unique horizontal lifts over base edges.
5. Lifting a path composes fiber isomorphisms and defines parallel transport.
6. Lifting a loop gives a fiber automorphism: holonomy, the discrete curvature observable.

The current C++ representation specializes to isomorphic copies of a single microscopic fiber
graph `F`. This is the complete-lift case of `InfraGaugeTheory`'s graph connections. The research
layer also represents non-isomorphic fibers and partial lifts; the evolution kernels still
specialize to homogeneous fibers.

## Calculate the symmetry group of the chosen fiber

The allowed change of local fiber frame is not supplied by a string such as `U(1)`. It is

`G = Aut(F)`,

computed by exhaustive graph-automorphism enumeration. A connection map `U_xy: F_x -> F_y`
changes under local frames `g_x, g_y` as

`U_xy -> g_y U_xy g_x^-1`.

Consequently loop holonomy changes only by conjugation at its basepoint. Its conjugacy class
inside `Aut(F)` is gauge invariant. Permutation cycle type alone can merge distinct classes
of this smaller group, so the current loop probe canonicalizes over the actual derived group.
Different fiber
graphs therefore produce different microscopic symmetry groups. Whether an effective continuous
group appears after coarse-graining is a result to measure, not an input.

## Coupling to base rewrites

For the elementary spatial rewrite

`(x--y) -> (x--w--y)`,

the old transport is factored through the fresh fiber:

`U_xy = U_wy U_xw`.

There is one factorization for every element of `Aut(F)`, but all factorizations form one orbit
under changes of frame at `w`. They must therefore not be counted as different physical multiway
branches. A separate amplitude construction assigns a normalized state over
this orbit. Distinct boundary holonomies remain orthogonal in that chosen
state space. This is a kinematic isometry, not a derived quantum evolution
law or a derivation of probabilities for distinct rewrite histories.

## Steps beyond the upstream building blocks

The first implemented extensions are:

- exact derivation and a noncommutativity test for the local automorphism group;
- gauge transformations of entire connections and conjugacy-class observables;
- a rewrite extension that preserves boundary parallel transport;
- proof that fresh-fiber connection factorizations form a single gauge orbit;
- a normalized amplitude construction over that orbit;
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

This defines kinematics, loop invariants, and a possible amplitude
representation—not matter or a force law. Its role in the research
program is to provide exact reference constructions against which to
compare candidate fibers and transport obtained from actual rewrites.
Choosing local operators on these fibers is a separate modeling step,
not a completed connection to Wolfram's underlying rewrite proposal.
