# GPST known limitations

## Sign-sensitive rejection at floating-point roundoff scale

The retained reference is v1.0.0, with the corresponding accelerated v1.1
implementation. Each original suite on the tested Linux reference, native C
and forced NumPy paths has 39 passes and one failure:
`test_uniform_incompatible_force` expects `MOMENT_MODEL_MISFIT` but returns
`IDENTIFIED`. The failure remains recorded; it is not changed to a pass.

The fixture's combined point estimate should be zero in exact arithmetic.
The stored estimates differ in the sign of floating-point roundoff:

| Record | Combined mean_N | Combined mean_D | Point | SE | Status |
|---|---|---|---|---|---|
| Archived baseline replay | -6.264279218805137e-19 | 1.000000000000002 | -6.264279218805124e-19 | 0.0 | MOMENT_MODEL_MISFIT |
| Current Linux diagnostic | 5.889493282637308e-20 | 1.000000000000002 | 5.889493282637296e-20 | 1.0766275688880601e-34 | IDENTIFIED |

Sources: the archived baseline replay `compatibility-risk.json`, field
`baseline_primary`, and its `COMPATIBILITY-CLARIFICATION.md`; current package
reference job 2386636, accelerated C job 2386640, and forced NumPy/diagnostic
job 2386642. The second row is copied from job 2386642's
`INCOMPATIBLE-FORCE-DIAGNOSIS.json`, `CF_GMM` fields. It is the combined
estimate, not the raw N0/D0 channel mean. These validation records are
retained outside the distributed overlay.

Version 1.0.0 has no numerical-zero gate. In `inference.fieller`, a finite
identified ratio with a negative point is changed to `MOMENT_MODEL_MISFIT`;
a tiny positive point remains `IDENTIFIED`. Thus rejection of this degenerate
input depends on a rounding sign and can vary across platforms. The archived
clarification establishes the observed negative-point rejection path, rather
than a training-time null-conflict exception. This limitation predates
packaging. Algorithms, tests, assertions and tolerances are unchanged; the
later v1.2 numerical-zero gate is not part of this manuscript version.

## User rule and scope

Do not use an `IDENTIFIED` result as a temperature if any of the following
holds: the combined point estimate is at floating-point roundoff scale,
its SE is zero, or the denominator interval is unresolved. A status label or
a narrow interval alone is insufficient; independent sampling qualification
and the force/unit contract remain necessary. This is a documented use rule,
not a newly implemented threshold or automatic rejection gate.

This specific issue concerns roundoff-scale combined estimates. None of the
manuscript's 14 temperature readouts (physical scale 150–1500 K, nonzero SE)
is in that regime, so those readings are unaffected by this limitation.
This scope statement is not a guarantee for arbitrary incompatible forces,
other degenerate inputs, or unqualified samples.
