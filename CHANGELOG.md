# Changelog

## 0.1.0 (2026-11-25)

First public release, accompanying arXiv:2611.32129.

- Rewrite the README for new users: what the package does, installation from
  the public SpinGNN tag `v0.1.0`, a quick start, GPST usage, limitations and
  citation. Add `CITATION.cff`.
- Point the installation baseline to SpinGNN tag `v0.1.0`, whose 14 files
  checked by the installer are byte-identical to the validated baseline.
- Drop references to development-history commits from the documentation and
  the LAMMPS modification notice; update the installer's payload and
  installed-file hashes for the four changed package files and NOTICE.md.
  No source file is changed.

## Licensing (2026-11-23)

- License the offline GPST tools in `tools/gpst/` under AGPL-3.0-only
  (`tools/gpst/LICENSE`). The LAMMPS overlay and the rest of the repository
  remain GPL-2.0-only.

## Documentation update (2026-11-23)

- Rename the method in all documentation to the cached-field integrator (CFI),
  following the manuscript. The fix style `spin/lpsib`, which keeps the earlier
  working name (lagged-predictor SIB), and all source files and class names
  are unchanged.

## Manuscript code candidate

- Unify frozen, moving-lattice and isotropic canonical-volume CFI, retaining
  the centered SIB predictor control and independent built-in OU streams.
- Keep the manuscript call-count observer unchanged, including native SIB and
  Verlet wrappers. Accepted-state output evaluations are counted separately.
- Include the mechanical-force provider for moving CFI. PairSpin resets its
  lattice flag on each initialization. Removing a moving integrator before a
  later run now restores the default flag rather than reusing stale state.
  Allegro does not read this flag; analytic spin pairs may change behavior.
- Retain native SIB-NPT's original cache behavior and the set-spin atom-variable
  mapping correction. Remove experimental fixes and their dedicated name gates;
  the former opt-in cache repair is not part of this release.
- Relocate the unchanged ML pair sources during installation without shipping
  models. Add explicit support boundaries and runnable analytic dump examples.
- Collect offline GPST reference and accelerated sources with their manifests,
  original tests and explicit-unit dump adapters. The incompatible-force status
  assertion is platform-dependent at floating-point roundoff scale, as shown
  by the archived baseline replay and current Linux records. Document the
  user rule without changing algorithms, tests, assertions or tolerances;
  see tools/gpst/KNOWN-LIMITATIONS.md.
