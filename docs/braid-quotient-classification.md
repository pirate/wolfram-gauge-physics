# An exact limit of the current square-fiber interaction family

Increasing rule-table complexity inside the current census does **not** let hidden fiber
representatives feed back into one-cell conjugacy-sector dynamics. For all 5,731 strict braid
solutions, the maximal additive class-charge description is an autonomous block-exchange
process. This includes the identity and subgroup-preserving rules; no minimum-support filter
is applied. It is a finite-model classification, not a no-go theorem for Wolfram models or
gauge physics in general.

## What was checked

`tools/classify_braid_quotients.py` independently reconstructs $G=\operatorname{Aut}(C_4)$,
checks all stored strict solutions for flat fixing, bijectivity, involution, ordered-product
conservation, conjugation covariance, orientation reversal, and all 512 braid triples.
For each table it solves the rational additive class-charge equations afresh, then tests
factor closure on every raw pair. All five one-cell conjugacy sectors also close.

The resulting 66 labeled charge-update tables, their charges, partitions, and complete member
rule IDs are recorded in `data/d4-braid-quotient-classification.json`. The original admitted
family and why its enumeration is exhaustive are described in the
[interaction search](equivariant-rule-search.md).

| Additive charge rank, with $q(1)=0$ | Charge symbols | Number of microscopic rules | Distinct labeled charge tables |
|---|---:|---:|---:|
| 3 | 4 | 4,262 | 14 |
| 4 | 5 | 1,469 | 52 |

The three-charge case merges the identity and central half-turn; the four-charge case
distinguishes all conjugacy classes. The counts of distinct labeled tables are the partitions
of four symbols except the single-block partition, and all partitions of five symbols.
Relabeling symmetry is **not** divided out in the 66 count. Microscopic tables can remain
different even when their charge-update tables coincide.

## Why block exchange follows

The charge vectors in this census force preservation of the unordered pair of projected
symbols. Consequently the effective local update either fixes $(a,b)$ or exchanges it to
$(b,a)$. Bijectivity makes the exchange decision symmetric. Write $a\sim b$ when that pair
does not exchange. This relation is reflexive and symmetric.

It must also be transitive under the braid relation. Suppose $a\sim b$ and $b\sim c$ but
$a\not\sim c$. Applying the two braid words to $(a,c,b)$ gives $(a,c,b)$ and $(c,a,b)$,
respectively, a contradiction. Conversely, a partition into non-exchange blocks satisfies
the braid relation: triples belong to one block, two blocks, or three blocks, and the two
words agree in each case. Therefore

$$f(a,b)=\begin{cases}(a,b),&[a]=[b],\\(b,a),&[a]\ne[b],\end{cases}$$

where $[a]$ denotes the block containing $a$. The tests independently enumerate all
$2^{10}=1,024$ symmetric swap/no-swap gates on five symbols. Exactly 52 pass the braid
relation, and these are exactly the partition gates.

This block-exchange structure is related to the unsigned generalized-permutator construction
in the integrable-model literature; see [Anfossi, Dolcini and Montorsi, *Recent results on
integrable electronic models*](https://arxiv.org/abs/cond-mat/0412532). We do not claim the
partition construction itself is new, or infer a quantum electronic Hamiltonian from a
classical permutation gate. The specific result here is the exhaustive map from our constrained
square-fiber gauge tables to these effective processes.

## The stronger sector-closure obstruction does not require braiding

Every one of the 29 minimal symmetry-closed involutions (including the identity) already has
an autonomous **full conjugacy-sector** factor. These local factors are stored with the
classification. Every admitted involution is a disjoint union of compatible minimal closures.
Since their supports are disjoint, this union equals their composition. Factor maps compose.
Therefore all **1,769,472 admitted involutions** have autonomous one-cell sector dynamics,
including those that fail the braid relation.

This extension uses the exhaustive minimal-closure construction, not a second brute-force
scan of 1,769,472 tables in Python. It does not assert that all those non-braid maps have
block-exchange charge dynamics. Dropping the braid condition alone cannot unlock feedback
from relative fiber representatives into the one-cell sectors in this family.

## What this changes about the research direction

The [diffusion experiment](emergent-diffusion.md) is a real, exactly explained statistical
transport phenomenon inside one of these factors. But millions more cells or a more complicated
table from the same family cannot change its autonomous effective law.

This does **not** classify the full gauge-invariant state. Relative orientations, transported
multi-loop products, and other multi-cell observables are not determined by single-loop
conjugacy classes. Nor does it cover changing cell incidence, connector dynamics coupled to
additional degrees of freedom, changing fibers, non-involutive maps, or multiway amplitudes.
Closure survives arbitrary static connector choices because class functions are conjugation
invariant. Merely randomizing those connectors cannot evade the obstruction.

Next experiments should use a discriminating test rather than larger pictures: two states
with the same recorded local sectors but different gauge-invariant multi-loop information;
then determine whether an explicitly justified extension can turn that difference into a
different local observable history. Coupling changing cell/face incidence must carry an
explicit transport map and causal read/write dependencies. No continuum force, target
molecular shape, or prescribed potential should supply the missing dynamics.

## Reproduce

```bash
python3 tools/classify_braid_quotients.py --output out/braid-quotient-classification.json
diff -u data/d4-braid-quotient-classification.json out/braid-quotient-classification.json
python3 tests/effective_dynamics.py
```
