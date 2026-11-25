# Scope and known limitations

- CFI requires one MPI rank, metal units and Verlet. Moving mode requires
  group all, newton on, atomic spin atoms, no kspace/BEC and an orthogonal box.
- The supported lattice thermostat is the fix's internal OU bath. Overlapping
  integration, external lattice baths and nonconservative STT are refused.
- NPT is isotropic, orthogonal, periodic and single-rank. Testing covers analytic
  Fe only. ML NPT is untested and rejected by the pair-virial guard. NequIP's
  newton-off requirement is incompatible with moving CFI.
- Native spin/sib_npt retains its original cached-field behavior (the historical
  C3 defect). An opt-in cache-refresh variant from earlier development is
  not part of this release.
- Native SIB is the normalized predicted-midpoint implementation. A theorem
  for a different SIB equation does not establish its coupled moving-lattice
  order. CFI's second-order statement is deterministic NVE only, with zero
  damping and no baths; no stochastic weak-second-order claim is made.
- Alpha equal to one is for thermodynamic sampling. Do not tune bath targets to
  force agreement with a temperature readout.
- Conservative spin force is f = hbar * fm only for unit spin directions.
  The simplified Nurdin readout in compute spin is valid only when energy is
  linear in each spin and is not a general ML-potential thermometer.
- The historical dipole/long fmi-versus-fi mechanical-force finding is pending
  separate investigation. It is neither repaired nor certified here.
- ML execution is not revalidated on a GPU in this packaging task. Source
  identity and prior accepted execution are separate from new package tests.
- GPST confidence/status codes are retained unchanged. A finite number alone
  is not a sampling qualification or proof of an identified estimator.
