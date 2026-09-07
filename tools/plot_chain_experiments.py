#!/usr/bin/env python3
"""Render measured gauge sectors from the chain experiment, without embedding a graph."""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm, ListedColormap

parser = argparse.ArgumentParser()
parser.add_argument("--input", type=Path, default=Path("data/d4-chain-experiments.json"))
parser.add_argument("--output", type=Path, default=Path("docs/images/chain-schedule-control.png"))
args = parser.parse_args()
data = json.loads(args.input.read_text())
fig, axes = plt.subplots(1, 3, figsize=(13, 4.6), constrained_layout=True)
palette = ListedColormap(["#f2f4f6", "#2976ba", "#e88938", "#7654a2", "#c04456", "#208579", "#9a8528", "#484d59"])
norm = BoundaryNorm([i - .5 for i in range(9)], palette.N)
for phase in (0, 1):
    run = next(r for r in data["runs"] if r["cells"] == 128 and r["seed"] == "rotation"
               and r["rule"] == "boundary_shear" and r["schedule_phase"] == phase)
    profiles = [sample["sector_profile"] for sample in run["samples"]]
    axes[phase].imshow(profiles, origin="lower", interpolation="nearest", aspect="auto",
                       extent=(-64.5, 63.5, -2, 34), cmap=palette, norm=norm)
    axes[phase].set(xlim=(-36, 36), ylim=(0, 32), title=f"Boundary shear · phase {phase}",
                    xlabel="Cell offset from seed", ylabel="Update layers (sampled every 4)")
for rule, phase, label, style in (("boundary_shear", 0, "Shear, phase 0", "-"),
                                  ("boundary_shear", 1, "Shear, phase 1", "--"),
                                  ("hurwitz", 0, "Hurwitz, phase 0", ":")):
    run = next(r for r in data["runs"] if r["cells"] == 1024 and r["seed"] == "rotation"
               and r["rule"] == rule and r["schedule_phase"] == phase)
    axes[2].plot([s["layer"] for s in run["samples"]], [s["active_cells"] for s in run["samples"]],
                 style, label=label, linewidth=2)
axes[2].set(title="Occupied support · 1,024 cells", xlabel="Update layers", ylabel="Nonidentity cell holonomies")
axes[2].legend(frameon=False)
axes[2].grid(alpha=.2)
fig.suptitle("One quarter-turn seed: schedule-dependent spreading on a supplied chain", fontsize=14)
fig.supxlabel("Heatmap colors distinguish exact gauge sectors; background is flat. No continuum geometry or physical clock inferred.", fontsize=9)
args.output.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(args.output, dpi=180)
print(args.output)
