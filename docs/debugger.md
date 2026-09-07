# Simulation debugger

Serve the repository with `python3 -m http.server 8765 --bind 127.0.0.1` and open [the debugger](http://127.0.0.1:8765/viewer/). Load an evolution JSON or use the included real five-step export. The viewer has no external JavaScript dependencies.

![Actual debugger showing a bounded branchial slice at rewrite depth five, with 240 of 348 states displayed](images/simulation-debugger.png)

- **Hypergraph:** ordered hyperedge incidence in a selected vertex neighborhood. Click vertices to move the neighborhood, or diamonds to inspect exact ordered slots. Added and removed hyperedges are compared against the producing event's input state; removed edges are visible when all their vertices remain in the displayed patch.
- **Causal history:** the producing event and its causal ancestors, arranged by rewrite depth. Directed links are prerequisites. Click an event to select its output state. Play lineage follows the selected state's ancestry and then the first exported descendant at each fork; it does not combine mutually exclusive branches.
- **Branchial slice:** raw states at one rewrite depth. Exported branchial endpoints are event IDs, so the adapter maps them to output states and requires both outputs to belong to the slice. Click a state to synchronize the hypergraph, causal history, and projection views.
- **Euclidean projection:** classical metric multidimensional scaling of shortest-path distances on a connected, bounded induced patch of the hypergraph 2-section. Three positive spectral modes supply display coordinates. Rotate with dragging and zoom with the wheel. Edges are omitted to make the recovered distance arrangement inspectable.

The projection reports normalized distance stress, $\sqrt{\sum_{i<j}(\|x_i-x_j\|-d_{ij})^2/\sum_{i<j}d_{ij}^2}$. Low stress means that this patch's graph distances admit a good three-dimensional representation; it does not establish physical space or a continuum limit. Three display coordinates do not determine intrinsic dimension. Intrinsic volume and spectral dimension estimates come from the simulation export independently. Cropping can change shortest paths, and each state is embedded independently, so rotations between frames are not measured motion.

Causal dependencies, spatial locality, and entanglement are separate questions. This debugger shows event dependencies and exported branchial relations. The export has no amplitudes, reduced states, or quantitative entanglement observable. Rewrite depth supplies navigation through an execution history; it is not recovered physical proper time.

The renderer caps views at 240 nodes and 1,200 links; metric projection caps its patch at 96 vertices. Truncation is reported explicitly. Selectors show at most 500 ordinary choices plus the current selection. JSON imports are limited to 32 MiB and the whole export is indexed in memory. These are inspection budgets, not a claim of million-node support. Large simulation debugging still needs an engine-side query/window API, multiscale summaries, and stable landmark embeddings.

Run adapter checks with `node --test viewer/model.test.mjs`.
