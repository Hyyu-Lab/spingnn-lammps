# Using another potential

Install into the exact baseline shown in README, then enable SPIN and SPINGNN.
Keep models outside the overlay. Start with a separately validated atom/spin
configuration and a conservative spin-force contract.

Analytic potentials may use frozen or moving CFI. Start in deterministic
mode and audit energy/force consistency before adding independent baths.
``predictor fresh`` selects the centered SIB control; the default ``predictor lagged`` uses the cached field.
The examples demonstrate syntax and do not certify their timesteps for a new
material. Potential cutoff continuity and stability need separate checks.

For spin/allegro_spin, supply a compatible model path and atom-type mapping.
The CrSBr example is a template: replace both input-data and model placeholders.
For spin/nequip_spin, retain its newton-off requirement. Moving CFI requires
newton on, so that combination is refused. Do not infer support from an
Allegro test. KOKKOS variants remain in the baseline; no new KOKKOS certification
is supplied here. ML NPT is rejected by the current pair-virial guard.

Use internal ``ou T damp seed`` for atomic thermalization in moving CFI.
Do not add fix langevin to that fix. Native spin/sib moving inputs instead use
a separate lattice Langevin thermostat as shown in the retained SIB example.

For GPST, record unit spins and conservative fm at the same accepted state.
Use ``accept`` at the dump interval and feed the dump directly to
``tools/gpst/dump_to_moments.py`` with explicit native hbar and kB.
This reads stored frames without evaluating the potential. Sampling
qualification and independent-chain inference remain the caller's task.
