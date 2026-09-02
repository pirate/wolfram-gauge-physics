#!/usr/bin/env python3
"""Metal throughput probe for compact combinatorial gauge operations.

This is not an end-to-end engine benchmark. It validates the device data shape
used by src/infragauge.hpp on Apple Silicon: uint16 group tables, fused local
frame transforms, and ordered path products.
"""

from __future__ import annotations

import argparse
import itertools
import json
import time

import mlx.core as mx
import numpy as np


def compose(after: tuple[int, ...], before: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(after[before[i]] for i in range(len(after)))


def inverse(value: tuple[int, ...]) -> tuple[int, ...]:
    result = [0] * len(value)
    for i, image in enumerate(value):
        result[image] = i
    return tuple(result)


def square_automorphisms() -> list[tuple[int, ...]]:
    edges = {tuple(sorted(edge)) for edge in [(0, 1), (1, 2), (2, 3), (3, 0)]}
    result = []
    for candidate in itertools.permutations(range(4)):
        mapped = {
            tuple(sorted((candidate[u], candidate[v])))
            for u, v in itertools.combinations(range(4), 2)
            if (u, v) in edges
        }
        if mapped == edges:
            result.append(candidate)
    return result


def group_tables() -> tuple[np.ndarray, np.ndarray]:
    group = square_automorphisms()
    index = {value: i for i, value in enumerate(group)}
    multiplication = np.empty((len(group), len(group)), dtype=np.uint16)
    inverses = np.empty(len(group), dtype=np.uint16)
    for left, left_value in enumerate(group):
        inverses[left] = index[inverse(left_value)]
        for right, right_value in enumerate(group):
            multiplication[left, right] = index[compose(left_value, right_value)]
    return multiplication.reshape(-1), inverses


FRAME_SOURCE = r"""
    uint i = thread_position_in_grid.x;
    uint inner = uint(mul[uint(transport[i]) * ORDER + uint(inv[source_frame[i]])]);
    out[i] = mul[uint(target_frame[i]) * ORDER + inner];
"""

PATH_SOURCE = r"""
    uint path = thread_position_in_grid.x;
    uint value = 0;
    for (uint step = 0; step < PATH_LENGTH; ++step) {
        uint edge_value = uint(paths[path * PATH_LENGTH + step]);
        value = uint(mul[edge_value * ORDER + value]);
    }
    out[path] = ushort(value);
"""


def elapsed_ms(operation, iterations: int) -> float:
    for _ in range(5):
        mx.eval(operation())
    start = time.perf_counter()
    for _ in range(iterations):
        mx.eval(operation())
    return (time.perf_counter() - start) * 1000.0 / iterations


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--items", type=int, default=1_000_000)
    parser.add_argument("--path-length", type=int, default=8)
    parser.add_argument("--iterations", type=int, default=40)
    args = parser.parse_args()

    multiplication, inverses = group_tables()
    order = len(inverses)
    rng = np.random.default_rng(42)
    transport_np = rng.integers(0, order, args.items, dtype=np.uint16)
    source_np = rng.integers(0, order, args.items, dtype=np.uint16)
    target_np = rng.integers(0, order, args.items, dtype=np.uint16)
    paths_np = rng.integers(
        0, order, (args.items, args.path_length), dtype=np.uint16
    )

    mul = mx.array(multiplication)
    inv = mx.array(inverses)
    transport = mx.array(transport_np)
    source = mx.array(source_np)
    target = mx.array(target_np)
    paths = mx.array(paths_np.reshape(-1))

    frame_kernel = mx.fast.metal_kernel(
        name="wgphysics_frame_transform_v1",
        input_names=["mul", "inv", "transport", "source_frame", "target_frame"],
        output_names=["out"],
        source=FRAME_SOURCE,
    )
    path_kernel = mx.fast.metal_kernel(
        name="wgphysics_path_transport_v1",
        input_names=["mul", "paths"],
        output_names=["out"],
        source=PATH_SOURCE,
    )

    def frame_operation():
        return frame_kernel(
            inputs=[mul, inv, transport, source, target],
            template=[("ORDER", order)],
            grid=(args.items, 1, 1),
            threadgroup=(256, 1, 1),
            output_shapes=[(args.items,)],
            output_dtypes=[mx.uint16],
        )[0]

    def path_operation():
        return path_kernel(
            inputs=[mul, paths],
            template=[("ORDER", order), ("PATH_LENGTH", args.path_length)],
            grid=(args.items, 1, 1),
            threadgroup=(256, 1, 1),
            output_shapes=[(args.items,)],
            output_dtypes=[mx.uint16],
        )[0]

    frame_result = np.asarray(frame_operation())
    frame_reference = multiplication.reshape(order, order)[
        target_np,
        multiplication.reshape(order, order)[transport_np, inverses[source_np]],
    ]
    if not np.array_equal(frame_result, frame_reference):
        raise RuntimeError("Metal frame transform disagrees with CPU table oracle")

    path_result = np.asarray(path_operation())
    path_reference = np.zeros(args.items, dtype=np.uint16)
    table = multiplication.reshape(order, order)
    for step in range(args.path_length):
        path_reference = table[paths_np[:, step], path_reference]
    if not np.array_equal(path_result, path_reference):
        raise RuntimeError("Metal path transport disagrees with CPU table oracle")

    frame_ms = elapsed_ms(frame_operation, args.iterations)
    path_ms = elapsed_ms(path_operation, args.iterations)
    report = {
        "device": "Apple Metal via MLX",
        "fiber": "cycle graph C4",
        "derived_group_order": order,
        "items": args.items,
        "path_length": args.path_length,
        "iterations": args.iterations,
        "frame_transform_ms": frame_ms,
        "frame_transforms_per_second": args.items / (frame_ms / 1000.0),
        "path_transport_ms": path_ms,
        "group_multiplies_per_second":
            args.items * args.path_length / (path_ms / 1000.0),
        "verification": "exact CPU table equality",
        "scope": "kernel microbenchmark, not end-to-end rewriting",
    }
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
