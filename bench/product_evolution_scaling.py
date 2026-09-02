#!/usr/bin/env python3
"""Measure exact raw-to-physical product evolution without changing its semantics."""

from __future__ import annotations

import argparse
import json
import platform
import statistics
import subprocess
import tempfile
import time
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--executable", default="build/wgphysics_product_evolve")
    parser.add_argument("--max-steps", type=int, default=4)
    parser.add_argument("--repetitions", type=int, default=3)
    args = parser.parse_args()
    if args.max_steps < 0 or args.repetitions < 1:
        parser.error("max-steps must be nonnegative and repetitions must be positive")

    executable = str(Path(args.executable).resolve())
    rows: list[dict[str, object]] = []
    with tempfile.TemporaryDirectory(prefix="wgphysics-product-bench-") as temporary:
        for steps in range(args.max_steps + 1):
            durations: list[float] = []
            expected_counts: dict[str, int] | None = None
            for repetition in range(args.repetitions):
                destination = Path(temporary) / f"depth-{steps}-{repetition}.json"
                started = time.perf_counter()
                subprocess.run(
                    [executable, "--steps", str(steps), "--output", str(destination)],
                    check=True,
                    stdout=subprocess.DEVNULL,
                )
                durations.append(time.perf_counter() - started)
                counts = json.loads(destination.read_text())["counts"]
                if expected_counts is None:
                    expected_counts = counts
                elif counts != expected_counts:
                    raise RuntimeError(f"non-deterministic counts at depth {steps}")
            assert expected_counts is not None
            rows.append(
                {
                    "steps": steps,
                    **expected_counts,
                    "minimum_seconds": min(durations),
                    "median_seconds": statistics.median(durations),
                }
            )

    print(
        json.dumps(
            {
                "date": "2026-09-02",
                "platform": platform.platform(),
                "rule": "{x,y}->{x,w},{w,y}",
                "initial_base": "triangle",
                "fiber": "C4",
                "derived_group_order": 8,
                "repetitions": args.repetitions,
                "scope": "engine evolution plus exact CPU product propagation and canonicalization",
                "measurements": rows,
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
