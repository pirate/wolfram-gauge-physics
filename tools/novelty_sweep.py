#!/usr/bin/env python3
"""Explore a bounded local-rule family by intrinsic phenotype novelty.

The selector never receives a target dimension, geometry, or molecule template.
It keeps rules that span the observed space of graph growth, branching,
canonical compression, diffusion, and intrinsic volume-growth behavior.
"""

from __future__ import annotations

import argparse
import itertools
import json
import math
import statistics
import subprocess
import tempfile
from pathlib import Path


LHS = ((0, 1), (0, 2))
INITIAL = "1,2;1,3"
POSSIBLE_RHS_EDGES = tuple(itertools.combinations(range(4), 2))


def connected(edges: tuple[tuple[int, int], ...]) -> bool:
    reached = {0}
    while True:
        expanded = reached | {
            endpoint
            for left, right in edges
            if left in reached or right in reached
            for endpoint in (left, right)
        }
        if expanded == reached:
            return expanded == {0, 1, 2, 3}
        reached = expanded


def candidate_rules() -> list[str]:
    candidates = []
    lhs = ";".join(f"{left},{right}" for left, right in LHS)
    for edge_count in range(3, len(POSSIBLE_RHS_EDGES) + 1):
        for rhs in itertools.combinations(POSSIBLE_RHS_EDGES, edge_count):
            if not connected(rhs) or not any(3 in edge for edge in rhs):
                continue
            encoded_rhs = ";".join(f"{left},{right}" for left, right in rhs)
            candidates.append(f"{lhs}->{encoded_rhs}")
    return candidates


def finite_mean(values: list[float | None]) -> float:
    finite = [value for value in values if value is not None and math.isfinite(value)]
    return statistics.fmean(finite) if finite else 0.0


def phenotype(result: dict) -> tuple[list[float], dict]:
    states_by_class = {}
    for state in result["states"]:
        states_by_class.setdefault(state["canonical_id"], state)
    terminal_step = max(state["step"] for state in states_by_class.values())
    terminal = [state for state in states_by_class.values() if state["step"] == terminal_step]
    initial_vertices = len({vertex for edge in result["states"][0]["edges"] for vertex in edge})

    mean_vertices = statistics.fmean(state["observables"]["vertices"] for state in terminal)
    mean_degree = statistics.fmean(state["observables"]["mean_degree"] for state in terminal)
    volume_dimension = statistics.fmean(
        finite_mean(state["intrinsic_geometry"]["volume_dimension"]) for state in terminal
    )
    spectral_dimension = statistics.fmean(
        finite_mean(state["intrinsic_geometry"]["spectral_dimension"]) for state in terminal
    )
    inhomogeneity = statistics.fmean(
        max(state["intrinsic_geometry"]["ball_volume_cv"], default=0.0) for state in terminal
    )
    counts = result["counts"]
    compression = 1.0 - counts["canonical_states"] / counts["raw_states"]
    branching = counts["events"] / max(1, counts["canonical_states"] - 1)
    causal_density = counts["causal_edges"] / max(1, counts["events"])

    named = {
        "growth_ratio": mean_vertices / max(1, initial_vertices),
        "mean_degree": mean_degree,
        "volume_dimension": volume_dimension,
        "spectral_dimension": spectral_dimension,
        "inhomogeneity": inhomogeneity,
        "canonical_compression": compression,
        "branching": branching,
        "causal_density": causal_density,
    }
    return list(named.values()), named


def standardized(vectors: list[list[float]]) -> list[list[float]]:
    columns = list(zip(*vectors))
    means = [statistics.fmean(column) for column in columns]
    deviations = [statistics.pstdev(column) for column in columns]
    return [
        [
            (value - means[index]) / deviations[index] if deviations[index] > 1e-12 else 0.0
            for index, value in enumerate(vector)
        ]
        for vector in vectors
    ]


def distance(left: list[float], right: list[float]) -> float:
    return math.sqrt(sum((a - b) ** 2 for a, b in zip(left, right)))


def diverse_subset(vectors: list[list[float]], keep: int) -> tuple[list[int], dict[int, float]]:
    if not vectors or keep <= 0:
        return [], {}
    normalized = standardized(vectors)
    pairwise = [
        [distance(normalized[left], normalized[right]) for right in range(len(normalized))]
        for left in range(len(normalized))
    ]
    medoid = min(range(len(normalized)), key=lambda index: sum(pairwise[index]))
    selected = [medoid]
    novelty = {medoid: 0.0}
    while len(selected) < min(keep, len(normalized)):
        candidate, separation = max(
            (
                (index, min(pairwise[index][chosen] for chosen in selected))
                for index in range(len(normalized))
                if index not in selected
            ),
            key=lambda item: (item[1], -item[0]),
        )
        selected.append(candidate)
        novelty[candidate] = separation
    return selected, novelty


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", type=Path, default=Path("build/wgphysics_evolve"))
    parser.add_argument("--output", type=Path, default=Path("out/novelty-sweep.json"))
    parser.add_argument("--steps", type=int, default=3)
    parser.add_argument("--probe-radius", type=int, default=8)
    parser.add_argument("--diffusion-steps", type=int, default=16)
    parser.add_argument("--probe-sources", type=int, default=64)
    parser.add_argument("--maximum-rules", type=int, default=0)
    parser.add_argument("--keep", type=int, default=8)
    parser.add_argument("--timeout-seconds", type=float, default=15.0)
    args = parser.parse_args()
    if args.steps < 0 or args.keep < 0 or args.timeout_seconds <= 0:
        raise SystemExit("steps and keep must be nonnegative; timeout must be positive")

    rules = candidate_rules()
    if args.maximum_rules:
        rules = rules[: args.maximum_rules]
    records = []
    with tempfile.TemporaryDirectory(prefix="wgphysics-novelty-") as temporary:
        temporary_path = Path(temporary)
        for index, rule in enumerate(rules):
            output = temporary_path / f"candidate-{index}.json"
            command = [
                str(args.binary),
                "--rule",
                rule,
                "--init",
                INITIAL,
                "--steps",
                str(args.steps),
                "--probe-radius",
                str(args.probe_radius),
                "--diffusion-steps",
                str(args.diffusion_steps),
                "--probe-sources",
                str(args.probe_sources),
                "--output",
                str(output),
            ]
            try:
                completed = subprocess.run(
                    command, capture_output=True, text=True, timeout=args.timeout_seconds
                )
            except subprocess.TimeoutExpired:
                records.append({"rule": rule, "status": "timeout"})
                continue
            if completed.returncode != 0:
                records.append(
                    {
                        "rule": rule,
                        "status": "failed",
                        "error": completed.stderr.strip().splitlines()[-1]
                        if completed.stderr.strip()
                        else "unknown failure",
                    }
                )
                continue
            result = json.loads(output.read_text())
            vector, named = phenotype(result)
            records.append(
                {
                    "rule": rule,
                    "status": "ok",
                    "counts": result["counts"],
                    "phenotype": named,
                    "_vector": vector,
                }
            )

    successful = [record for record in records if record["status"] == "ok"]
    selected, novelty = diverse_subset(
        [record["_vector"] for record in successful], args.keep
    )
    for index, record in enumerate(successful):
        record["selected_for_diversity"] = index in selected
        record["novelty_distance_at_selection"] = novelty.get(index)
        del record["_vector"]

    document = {
        "schema": 1,
        "semantics": "target-free novelty sweep over intrinsic graph phenotypes",
        "rule_family": {
            "lhs": "0,1;0,2",
            "rhs": "connected simple graphs on variables 0,1,2,3 with fresh variable 3",
            "initial": INITIAL,
        },
        "configuration": {
            "steps": args.steps,
            "probe_radius": args.probe_radius,
            "diffusion_steps": args.diffusion_steps,
            "probe_sources": args.probe_sources,
            "timeout_seconds": args.timeout_seconds,
            "candidate_rules": len(rules),
            "diverse_rules_requested": args.keep,
        },
        "selection": {
            "method": "greedy max-min distance after per-observable standardization",
            "seed": "phenotype medoid",
            "target_geometry": None,
        },
        "candidates": records,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(document, indent=2) + "\n")
    print(
        f"wrote {args.output} with {len(successful)} successful rules and "
        f"{len(selected)} diverse selections"
    )


if __name__ == "__main__":
    main()
