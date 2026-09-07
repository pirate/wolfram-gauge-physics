#!/usr/bin/env python3
"""Record CPU observer timings and their actual validation scope."""
import argparse
import json
import platform
import statistics
import subprocess
from pathlib import Path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, default=Path('out/loop-observer-cpu.json'))
    args = parser.parse_args()
    process = subprocess.run(['build/wgphysics_loop_bench'], capture_output=True, text=True, check=True, timeout=300)
    result = json.loads(process.stdout)
    result['machine'] = {'system': platform.system(), 'architecture': platform.machine()}
    if platform.system() == 'Darwin':
        result['machine']['processor'] = subprocess.check_output(['sysctl', '-n', 'machdep.cpu.brand_string'], text=True).strip()
    for row in result['runs']:
        row['median_probe_ms'] = statistics.median(row['probe_ms'])
        row['median_legacy_ms'] = statistics.median(row['legacy_ms']) if row['legacy_ms'] else None
        row['warm_probe_speedup_over_legacy'] = row['median_legacy_ms']/row['median_probe_ms'] if row['legacy_ms'] else None
    result['scope'] = ('Precompiled-topology CPU gauge observation, not graph construction, rewrite evolution, rendering, '
                       'or GPU throughput. Timed probes include output allocations but exclude validation comparisons. '
                       'Legacy comparisons use different canonical encodings and include repeated topology/path work; '
                       'million-edge cases have frame and reconstruction checks, not an executed legacy comparison.')
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2)+'\n')
    print([(r['edges'], r['median_probe_ms']) for r in result['runs']])


if __name__ == '__main__':
    main()
