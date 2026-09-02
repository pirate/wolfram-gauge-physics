# Contributing

Contributions are welcome, especially small constructions with exact finite tests.

Please preserve these rules:

1. Distinguish gauge-dependent representations from gauge-invariant observables.
2. Never count a local-frame copy as a distinct physical branch.
3. Keep exact CPU behavior as the oracle for accelerated implementations.
4. State approximation error explicitly; do not silently prune branches.
5. Avoid assigning particle or force names to a structure before behavioral tests justify them.
6. Cite the upstream definition or public research result a new construction extends.

Before opening a pull request:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build -L wgphysics --output-on-failure
```
