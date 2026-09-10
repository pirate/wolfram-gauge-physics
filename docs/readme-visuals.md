# README figure sources

## Interview-first explanations

Run `node tools/render_primitives_visuals.mjs` to regenerate three
specified examples. They are explanatory diagrams, not measured physical
states or a claim that the open constructions have been completed.

- **Rewrite primitives:** the exact labeled rule in the README's export
  command. It consumes `(0,1),(0,2)` and produces
  `(0,2),(0,3),(1,3),(2,3)`. The shared endpoint pair is consumed and
  reproduced as a new edge occurrence. Layout coordinates have no units.
- **History versus causality:** A subdivides `(0,1)` through node 3;
  B subdivides `(1,2)` through node 4. The generator verifies that either
  order reaches the same labeled edge collection. A later joining event C
  consumes `(3,1),(1,4)` and produces `(3,4)`. It therefore depends on both
  A and B, although A and B do not depend on one another. This illustrates
  two specified rules, not every possible history. Graph-isomorphism
  reduction can also identify the two intermediate states.
- **Fiber projection:** two connected graphs with endpoints 0 and 3,
  one a four-node path and one a four-node cycle. A deliberately supplied
  observer retains only the endpoints and their connectivity. Both map
  to one coarse description; the generator checks their connectivity.
  Their common projection does not imply gauge equivalence or equal
  future behavior. The figure distinguishes a set-theoretic fiber from
  the unresolved derivation of useful internal structure and transport.

All three figures have SVG titles/descriptions and use text as well as
color. Their labeled graphs can be reconstructed without reading pixel
coordinates. The main text pairs them with actual engine exports and an
attributed Wolfram documentation screenshot.

The README distinguishes specified mathematical examples, recorded model
states, and screenshots of external tools. Drawing coordinates are not
inferred physical coordinates.

## Reproduce the four explanatory figures

From the repository root, run:

```bash
node tools/render_readme_visuals.mjs
```

The script uses only Node.js standard libraries and the checked-in
[visual data extract](images/readme-visual-data.json). It does not run or
modify the simulation. SVG titles and descriptions provide text alternatives;
colors are accompanied by labels or numbers.

- **What is a fiber:** an exploded schematic of three base nodes, each
  with a triangle fiber. Dashed lines denote projection to the base, not
  physical edges or an extra spatial dimension.
- **Fiber holonomy:** a specified three-link connection with a triangle
  fiber. Two identity maps and the rotation `v -> v + 1 (mod 3)` compose to
  nonidentity holonomy. The fourth panel returns to the first base vertex.
  This is an explanatory example, not sampled evolution.
- **Transport memory:** the two generator permutations from
  `data/fiber-transport-order.json`, for the noncommuting triangle-fiber
  preparation. The rows evaluate the two circuit orders starting at gauge
  state zero. They do not draw a path through physical space. See the
  [complete calculation](noncommuting-fiber-transport.md).
- **Reaction transfer:** the `two_reaction_transport_witness` record from
  `data/fiber-reaction-carriers.json`. The script reconstructs the documented
  side-four mesh ordering and reads each face holonomy from the three saved
  raw-link vectors. All 96 face charges must match the saved values before
  rendering. Cyan edges are the raw links that differ between consecutive
  states. The five-face patch crosses a periodic seam and is unwrapped for
  display. The bottom arrow denotes net charge balance, not a microscopic
  trajectory. See the [reaction calculation](fiber-reaction-bursts.md).

The data extract keeps these figures reproducible independently of later
changes to the larger experiment files. It is a subset of recorded data,
not a replacement simulation dataset.

## Wolfram website screenshot

[Wolfram's Branchial Graphs and Multiway Causal Graphs](https://www.wolframphysics.org/technical-introduction/the-updating-process-in-our-models/branchial-graphs-and-multiway-causal-graphs/)
is the source of [the official-page screenshot](images/wolfram-official-multiway.png),
captured in a browser on September 9, 2026. The page branding and navigation
are retained. It is a screenshot of Wolfram's visual documentation, not a
live playground or an output from this repository.

The screenshot is reproduced for attributed discussion of the displayed
construction. The original page and artwork belong to their respective
rights holders; this repository's MIT license does not relicense them.

Wolfram Cloud showed a maintenance notice during capture, and the registry
index was unavailable. No screenshot of an operational cloud notebook or
interactive registry is claimed.

## Other figures

The existing rewrite history, multiway graph, reaction, mode, and debugger
images retain their original data and appearance. Their surrounding README
captions explain what to read and link the corresponding mathematical notes.
The user's opening video and notebook screenshot table is unchanged.
