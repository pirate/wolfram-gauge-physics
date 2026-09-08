#!/usr/bin/env python3
"""Compare Python coordinate-update costs on identical microscopic schedules."""
import argparse
import hashlib
import json
import platform
import random
import statistics
import time
from pathlib import Path

from triangle_charge_current import CurrentProbe
from triangle_lazy_gauge import GaugeConnection


def measure(side, bank, attempts, repeats):
    e = CurrentProbe(side, bank).e
    rng = random.Random(920410+side)
    initial = [rng.randrange(e.geometry.group.n) for _ in e.geometry.edges]
    schedule = [divmod(rng.randrange(13*len(e.factor.pairs)), len(e.factor.pairs)) for _ in range(attempts)]
    timings = {name: [] for name in ('raw_links', 'lazy_tree_gauge', 'eager_tree_gauge')}
    reference = None
    for _ in range(repeats):
        for name in timings:
            raw = initial[:]
            lazy = GaugeConnection(e, initial) if name == 'lazy_tree_gauge' else None
            if name == 'eager_tree_gauge':
                raw = e.forest.representative([e.forest.based_loops(initial)[0]])
            changed = 0
            start = time.perf_counter()
            for rule, support in schedule:
                if lazy is not None:
                    changed += lazy.step(rule*len(e.factor.pairs)+support)['changed']
                else:
                    code = e.factor.oracle.code(raw, rule, support)
                    if code != e.tables[rule][code]:
                        e.factor.oracle.update(raw, rule, support, e.tables[rule])
                        changed += 1
                        if name == 'eager_tree_gauge':
                            raw = e.forest.representative([e.forest.based_loops(raw)[0]])
            timings[name].append(time.perf_counter()-start)
            # Initialization and the common final verification are outside the
            # timed evolution loop. Eager in-loop gauge restoration is included.
            loops = lazy.loops() if lazy is not None else e.forest.based_loops(raw)[0]
            result = (changed, tuple(loops))
            if reference is None:
                reference = result
            if result != reference:
                raise ValueError('benchmark adapters do not evolve the same full connection orbit')
    return {'side': side, 'vertices': e.geometry.size, 'links': len(e.geometry.edges),
            'attempts': attempts, 'events': reference[0], 'repeats': repeats,
            'final_loop_sha256': hashlib.sha256(bytes(reference[1])).hexdigest(),
            'seconds': timings, 'median_seconds': {name: statistics.median(times) for name, times in timings.items()}}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--attempts', type=int, default=20000)
    parser.add_argument('--repeats', type=int, default=5)
    parser.add_argument('--output', type=Path, default=Path('out/triangle-lazy-gauge-benchmark.json'))
    args = parser.parse_args()
    if args.attempts < 1 or args.repeats < 1:
        parser.error('attempts and repeats must be positive')
    bank = json.loads(Path('data/triangle-feedback.json').read_text())['bank']
    records = []
    for side in (3, 6, 12, 24):
        row = measure(side, bank, args.attempts, args.repeats)
        records.append(row)
        print(side, row['median_seconds'], flush=True)
    data = {'python': platform.python_version(), 'platform': platform.platform(),
            'records': records,
            'scope': 'Python adapter timings on identical schedules; not a native C++ or GPU benchmark. Initialization and final snapshots excluded. A gain over eager gauge restoration is not a gain over raw-link evolution.'}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(data, separators=(',', ':'))+'\n')


if __name__ == '__main__':
    main()
