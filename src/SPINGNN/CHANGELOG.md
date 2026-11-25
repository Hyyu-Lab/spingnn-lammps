# Package changelog

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
  regression remains unresolved; see VALIDATION.md before using this candidate.
