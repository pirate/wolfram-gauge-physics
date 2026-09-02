#!/usr/bin/env python3
"""Render README figures from wgphysics JSON exports using Graphviz."""

from __future__ import annotations

import argparse
import json
import subprocess
from collections import defaultdict
from pathlib import Path


PALETTE = [
    "#dbeafe",
    "#dcfce7",
    "#fef3c7",
    "#fce7f3",
    "#ede9fe",
    "#cffafe",
    "#ffedd5",
    "#e2e8f0",
]


def quote(value: str) -> str:
    return '"' + value.replace('"', '\\"') + '"'


def render(dot: str, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        ["dot", "-Tpng", "-Gdpi=180", "-o", str(destination)],
        input=dot,
        text=True,
        check=True,
    )


def rewrite_history(evolution: dict) -> str:
    states = {state["id"]: state for state in evolution["states"]}
    outgoing: dict[int, list[int]] = defaultdict(list)
    for event in evolution["events"]:
        outgoing[event["input"]].append(event["output"])

    lineage = [0]
    while outgoing.get(lineage[-1]):
        lineage.append(min(outgoing[lineage[-1]]))
        if len(lineage) == evolution["steps"] + 1:
            break

    lines = [
        "digraph rewrite_history {",
        '  graph [rankdir=LR, bgcolor="white", pad=0.35, nodesep=0.28, ranksep=1.0,',
        '         fontname="Helvetica", fontsize=18, labelloc=t,',
        '         label="One exact rewrite history exported by wgphysics_evolve"];',
        '  node [shape=circle, fixedsize=true, width=0.30, height=0.30, label="",',
        '        style=filled, fillcolor="#f8fafc", color="#334155", penwidth=1.7];',
        '  edge [color="#64748b", penwidth=1.8];',
        "  compound=true;",
    ]

    for state_id in lineage:
        state = states[state_id]
        lines.extend(
            [
                f"  subgraph cluster_{state_id} {{",
                '    color="#cbd5e1"; penwidth=1.2; style="rounded";',
                f"    label={quote(f'Step {state["step"]} · {state["observables"]["vertices"]} vertices · {state["observables"]["edges"]} edges')};",
                '    fontname="Helvetica"; fontsize=13; fontcolor="#334155";',
            ]
        )
        for vertex in sorted({v for edge in state["edges"] for v in edge}):
            lines.append(f"    s{state_id}_{vertex};")
        for left, right in state["edges"]:
            lines.append(f"    s{state_id}_{left} -> s{state_id}_{right} [dir=none];")
        lines.append("  }")

    for left, right in zip(lineage, lineage[1:]):
        left_vertex = states[left]["edges"][0][0]
        right_vertex = states[right]["edges"][0][0]
        lines.append(
            f"  s{left}_{left_vertex} -> s{right}_{right_vertex} "
            f"[ltail=cluster_{left}, lhead=cluster_{right}, color=\"#2563eb\", "
            'penwidth=2.4, label=" rewrite ", fontname="Helvetica", '
            'fontcolor="#2563eb", fontsize=11];'
        )

    lines.append("}")
    return "\n".join(lines)


def multiway_graph(evolution: dict) -> str:
    states = evolution["states"]
    by_step: dict[int, list[dict]] = defaultdict(list)
    canonical_ids = sorted({state["canonical_id"] for state in states})
    color_for = {
        canonical_id: PALETTE[index % len(PALETTE)]
        for index, canonical_id in enumerate(canonical_ids)
    }

    lines = [
        "digraph multiway {",
        '  graph [rankdir=LR, bgcolor="white", pad=0.35, nodesep=0.28, ranksep=0.95,',
        '         fontname="Helvetica", fontsize=18, labelloc=t,',
        f'         label="Multiway evolution: {evolution["counts"]["raw_states"]} raw states → {evolution["counts"]["canonical_states"]} canonical states"];',
        '  node [shape=ellipse, style=filled, color="#475569", penwidth=1.5,',
        '        fontname="Helvetica", fontsize=10, fontcolor="#0f172a", margin="0.08,0.05"];',
        '  edge [color="#94a3b8", penwidth=1.4, arrowsize=0.7,',
        '        fontname="Helvetica", fontsize=8, fontcolor="#64748b"];',
    ]

    canonical_counts: dict[int, int] = defaultdict(int)
    for state in states:
        by_step[state["step"]].append(state)
        canonical_counts[state["canonical_id"]] += 1
        duplicate = canonical_counts[state["canonical_id"]] > 1
        style = '"filled,dashed"' if duplicate else "filled"
        label = f'raw {state["id"]}\\nclass {state["canonical_id"]}'
        lines.append(
            f"  s{state['id']} [label={quote(label)}, fillcolor={quote(color_for[state['canonical_id']])}, style={style}];"
        )

    for step, step_states in sorted(by_step.items()):
        members = "; ".join(f"s{state['id']}" for state in step_states)
        lines.append(
            f"  {{ rank=same; step_{step} [shape=plaintext, style=\"\", "
            f"fontcolor=\"#475569\", label={quote(f'Step {step}')}]; {members}; }}"
        )

    for event in evolution["events"]:
        lines.append(
            f"  s{event['input']} -> s{event['output']} [label={quote(f'e{event["id"]}')}];"
        )

    lines.extend(
        [
            '  legend [shape=plaintext, style="", fontcolor="#475569", label="Dashed nodes are raw labelings already represented by an earlier canonical class"];',
            "  step_3 -> legend [style=invis];",
            "}",
        ]
    )
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("evolution", type=Path)
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()

    evolution = json.loads(args.evolution.read_text())
    render(rewrite_history(evolution), args.output_dir / "rewrite-history.png")
    render(multiway_graph(evolution), args.output_dir / "multiway-evolution.png")


if __name__ == "__main__":
    main()
