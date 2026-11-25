# Geometric-probe spin temperature (GPST)

**Known limitation:** `test_uniform_incompatible_force` has a roundoff-scale
combined point estimate whose sign controls rejection in the retained v1.0.0
code. Its status can vary across platforms; each tested Linux path passes
39 of the original 40 tests. `IDENTIFIED` is insufficient for use as a
temperature when the point is at roundoff scale, SE is zero, or the denominator
interval is unresolved. See [raw values, provenance and scope](KNOWN-LIMITATIONS.md)
and [validation results](../../VALIDATION.md). Sources and assertions are unchanged.

The reported observable is T_GP. The reference `tspin` implementation is
v1.0.0; `tspin_v1_1` is the independently tested accelerated implementation.
Source manifests are retained. NumPy and SciPy are required; pytest is needed
for tests. Neither implementation evaluates a model. Sampling qualification
and physical acceptance are separate from the estimator status.

This is an offline tool. The 66.5–510.9 microsecond measurements in the
manuscript refer to an online compute that is not included here. They are not
performance measurements of these offline scripts.

## From a CFI dump

The CFI examples write exactly the columns consumed by this adapter.
Match the `accept` interval and dump interval. Read native LAMMPS hplanck/(2*pi)
in eV ps and boltz in eV/K for your unit system; supply them explicitly.
Do not divide or multiply by magnetic moment magnitude. Unit spin directions
and conservative angular-frequency fm must refer to the same accepted state.
The adapter checks unit norms, IDs, box geometry and interval alignment; the
caller must explicitly confirm the conservative-force contract.

```sh
python tools/gpst/dump_to_moments.py trajectory.dump chain.npz \
  --hbar "$HBAR_EV_PS" --cutoff "$PROBE_CUTOFF_ANGSTROM" \
  --accept-every 10 --force-contract-confirmed
python tools/gpst/infer.py --kb "$KB_EV_K" --dt "$SAVED_INTERVAL_PS" \
  --seed-ids seed-a seed-b seed-c seed-d seed-e seed-f \
  --records a.npz b.npz c.npz d.npz e.npz f.npz --output estimate.json
```

The first command saves sufficient moments per chain. Supply independent chains
to the second command. `--dt` is the interval between saved frames in ps, not
the MD timestep. A returned status other than IDENTIFIED must not be presented
as a certified temperature. No bath target is supplied to the estimator.
The minimal fixture and adapter-negative tests are in tests/test_adapter.py.

## Tests and optional native kernel

```sh
python -m pytest -q -p no:cacheprovider tools/gpst/tests/test_reference.py
cc -std=c99 -O3 -fno-fast-math -ffp-contract=off -shared -fPIC \
  tools/gpst/accelerated/tspin_v1_1/edge_kernel.c -lm \
  -o tools/gpst/accelerated/tspin_v1_1/_edge_kernel.so
python -m pytest -q -p no:cacheprovider tools/gpst/tests/test_accelerated.py
TSPIN_NUMPY_ONLY=1 python -m pytest -q -p no:cacheprovider \
  tools/gpst/tests/test_accelerated.py
```

The shared object is generated locally and is not distributed. Do not compile
with fast-math. A native test must also assert that native._LIB is loaded;
otherwise it may silently test the NumPy fallback. Source provenance is in
PROVENANCE.json.

## License

The GPST tools in this directory, including the reference and accelerated
implementations, the dump adapter, the inference script and the tests, are
licensed under the GNU Affero General Public License, version 3 only
(SPDX: `AGPL-3.0-only`); see [LICENSE](LICENSE).
Copyright (C) 2026 Yuanmai Zheng. The source files are kept identical to their
recorded manifests, so they carry no per-file license header.
