# Validation scope

**Accepted with a documented GPST limitation.** The original
`test_uniform_incompatible_force` status assertion is platform-dependent for
a degenerate roundoff-scale combined point estimate. Each tested Linux path
still has 39 passes and one failure; this is not recorded as 40 passing tests.
The archived baseline replay identifies the sign-sensitive rejection path in
the retained v1.0.0 reference implementation. See
[GPST known limitations](tools/gpst/KNOWN-LIMITATIONS.md) for both original
numeric records, their sources and the user rule. No algorithm, test,
assertion or tolerance has been changed. The manuscript readings do not fall
in this roundoff-scale regime and are unaffected by this specific limitation.

## Installation and compilation

The validated build used SpinGNN development commit
`e227650a2ea113aff16047abfcd59d3dac815066`. The public SpinGNN tag `v0.1.0`
(`0ef17be1714725af563737de36248e6ff3e32395`) contains byte-identical copies of
all 14 baseline files that the installer checks. Installation into those files
is tested; a full build from the public tag has not been repeated.
The installer gates every changed baseline file and relocates the four ML
pair source/header files unchanged. Missing or changed baseline files and a
completed installation are refused before modifying the target. Final payload
hashes and the installed file set are verified separately.

The new CMake build uses GCC/G++ 10.5.0, OpenMPI 4.1.6, C++17, Release,
`-O3 -DNDEBUG -pthread -fopenmp`, one execution thread and no fast-math.
SPIN and SPINGNN are enabled. A compatible Torch build remains a dependency;
linking it is not a GPU model evaluation. Execution uses the explicitly selected
compatible dynamic loader. No global replacement of the C library is required.

A small-target precheck did not cover the full Torch configuration. Earlier
technical attempts failed on library resolution, quota account-name resolution
and the MKL include path. The completed build then exceeded the original
startup smoke timeout; its logs and binary were retained for a separate
completion-only check. Each failed attempt is retained separately. A later
successful build does not erase those failures.

## New CPU results

| Check | Result |
|---|---|
| Full CPU CMake compilation/link and separate true-MPI startup | PASS |
| Historical executable relinks, including the additional package baseline | 9 matching SHA256 values |
| Retained analytic parse/initialization modes | 13 PASS |
| Explicit error-line differences | 15 PASS |
| Style registration | Five retained entries occur once; three removed entries absent |
| Original-answer byte regressions | 185/185 identical, source and destination counts equal and nonzero |
| Call and communication ledgers | Five cases PASS |
| Accepted-output dynamics invariance | 64-step dumps identical with refresh off versus every four steps |
| Pressure conversion | 100 bar / 1602176.5 = 0.00006241509596477042 eV/Angstrom^3 |
| Atom-variable spin mapping | 128 rows, maximum absolute error zero |
| Documented dump adapter and six-record CLI; short single-record negative | PASS |

The byte-regression groups contain 30, 45, 24, 12, 31, 36, 6 and 1 cases,
covering the manuscript integrator versions, frozen precession/hybrid cases,
deterministic and stochastic box fixtures, native baselines and variable-spin
mapping. Historical answers were not regenerated. The short CLI records are
interface demonstrations and are not physically qualified samples.

Historical archives are relinked using their recorded link recipes before
comparing the new package. The old executable hashes must match their frozen
values; the new executable is not required to have an old hash. Original
inputs and original answers are checked by SHA256 before use. Regression
acceptance requires the entire sorted dump to match in bytes. It is not a
floating-point-tolerance comparison. Full stdout is not compared because the
call ledger includes additional diagnostics.

Parser tests of ML pair settings do not load coefficients or a model. New GPU
execution, NequIP moving-lattice execution, ML NPT, multiple MPI ranks, anisotropic
boxes and a full traditional-make build are not certified. NPT admits only its
explicitly validated analytic pair whitelist.

## GPST results and small fixture

| Path | Result |
|---|---|
| Reference original suite | 39 passed, 1 failed |
| Accelerated C original suite | 39 passed, 1 failed |
| Forced NumPy original suite | 39 passed, 1 failed |
| Adapter and fixed synthetic fixture | 5 passed |

All three original-suite failures are `test_uniform_incompatible_force`.
The accelerated C run completed before its collection job timed out while
starting the fallback work. The completed C result is retained; a new
remaining-only job collected fallback and adapter results without rerunning
the completed C suite. Collection-job completion does not make the failing assertion a PASS.
Release acceptance records this known limitation without changing that result.

The retained reference manifest is
`bd4655946a560db099433796fd0b128766392aea8813c5cf8f276a318a080480`;
the accelerated manifest is
`ee986f6d0be8802a1ec443aed317f1b4f2fc4acc5c0d654e8ee7df16765dfd09`.
Accelerated C availability and forced NumPy fallback are asserted explicitly.
The optional shared object is built for testing and is not distributed.
Historical cross-implementation equivalence is separate evidence and does not
replace the currently failing original test. A returned finite estimate or
`IDENTIFIED` alone is not sampling qualification or scientific certification.

The small adapter fixture has four fixed unit directions, explicit synthetic
hbar conversion, a pooled-ratio answer of two and an exact degenerate Fieller
interval at two. Negative cases identify nonunit directions, unsynchronized
frames, wrong fields, missing force-contract confirmation and invalid hbar.
The documented CFI example is exercised through the dump adapter. Six
independent short demonstration records are used for the inference-CLI
interface check. A single eleven-frame record is rejected with the original
`UNCALIBRATED_SINGLE_RECORD_N:11` error. These are interface demonstrations,
not qualified physical sampling.

### Platform-dependent degenerate fixture

The archived baseline replay, `compatibility-risk.json`
(`baseline_primary`), stores mean_N = -6.264279218805137e-19,
mean_D = 1.000000000000002, point = -6.264279218805124e-19, SE = 0.0,
and `MOMENT_MODEL_MISFIT`. Its `COMPATIBILITY-CLARIFICATION.md` establishes
that this observed rejection follows the negative-point branch in
`inference.fieller`, rather than proving a training-time null-conflict path.

The current Linux records (reference job 2386636, C job 2386640 and forced
NumPy/diagnostic job 2386642) return `IDENTIFIED`; the latter's
`INCOMPATIBLE-FORCE-DIAGNOSIS.json` stores combined mean_N =
5.889493282637308e-20, mean_D = 1.000000000000002, point =
5.889493282637296e-20 and SE = 1.0766275688880601e-34. These are combined
CF_GMM values, not the uncombined N0/D0 means.

Both point estimates are at floating-point roundoff scale. Version 1.0.0
has no numerical-zero gate: an otherwise identified negative ratio is marked
`MOMENT_MODEL_MISFIT`, whereas the tiny positive result remains `IDENTIFIED`.
This platform-sensitive behavior predates packaging. No v1.2 zero gate is
introduced. An `IDENTIFIED` result must not be used as a temperature if its
point is at roundoff scale, SE is zero, or its denominator interval remains
unresolved. Independent sampling qualification remains necessary.

This specific limitation concerns roundoff-scale combined estimates. None
of the manuscript's 14 physical-scale temperature readouts (150–1500 K,
nonzero SE) belongs to this regime; those readings are unaffected. This does
not certify arbitrary incompatible forces or replace sampling diagnostics.

## Reading the call ledger

With screen output enabled, every native Verlet run prints one `OBSERVE85`
JSON object, including zero-count runs. CFI prints `S85_NATIVE_CONSTANTS`; box mode prints
`S85N_UNITS`. This fixed-storage observer is retained byte-for-byte from the
manuscript implementation. It increments integer counters around existing pair,
forward/reverse communication and neighbor operations; it adds no such calls.
Formatting and output occur at run end. This is instrumentation overhead, not
an independently measured zero-cost claim.

Read `entered` and `completed` in each of `setup`, `core`, `fresh`, `sib_mid`
and `accept`. Setup is reported separately. CFI has one core evaluation
per step, centered SIB adds a fresh predictor evaluation, and native SIB adds
its midpoint evaluation. Accepted-state refresh is extra output work and does
not enter the dynamics or overwrite the cached predictor field. Forward,
reverse, neighbor-build and requested-virial counters describe those same
existing calls. They must not be folded into a core-only efficiency claim.

## Scope and provenance

Scientific claims and performance qualifications follow the manuscript fact
list. No new physical sampling or efficiency experiment is performed for this
package. Online probe timings do not measure the included offline CLI.
Source lineage, raw process logs, original-answer manifests, resource usage and
failure receipts are retained outside the overlay. The repository contains
neither models nor production trajectories or binary artifacts. See
KNOWN-LIMITATIONS in the package source and the top-level README.
