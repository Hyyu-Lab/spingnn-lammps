# Input demonstrations

The native SIB NVE/NVT/NPT inputs demonstrate its inherited syntax. CFI
examples cover frozen, moving NVE/NVT, centered SIB and canonical-volume NPT.
They use analytic exchange; moving CFI examples add shifted LJ. Numerical
parameters are demonstrations drawn from the validation fixtures, not a
certified timestep recommendation. CrSBr requires user-supplied data and model.
No model or production trajectory is included.

CFI dumps include id, position, unit spin, moment magnitude and conservative
fm in the field order expected by tools/gpst/dump_to_moments.py. The accept
and dump intervals match. Pass native hbar explicitly; do not substitute the
magnetic moment magnitude for unit spin length. See the GPST README.
