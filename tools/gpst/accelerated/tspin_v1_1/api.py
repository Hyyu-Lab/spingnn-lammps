"""Independent-seed inference, with explicit physical units and no target input."""

import numpy as np

from .geometry import FIXED
from .inference import summary

STATUS_CODES = (
    "BAD_INPUT",
    "NONUNIT_SPIN",
    "ID_MISMATCH",
    "FORCE_CONTRACT_MISSING",
    "ZERO_SIGNAL",
    "WEAK_DENOMINATOR",
    "UNBOUNDED_CI",
    "EMPTY_CI",
    "RANK_DEFICIENT",
    "MOMENT_MODEL_MISFIT",
    "JACKKNIFE_FAILED",
    "INSUFFICIENT_SEEDS",
    "INSUFFICIENT_BLOCKS",
    "DIAGNOSTIC_ONLY",
    "NUMERICAL_NOISE_UNCERTIFIED",
    "UNCALIBRATED",
    "IDENTIFIED",
)


def pooled_ratio(n, d, *, kb, method="pooled"):
    """Pool numerator and denominator before division; frame ratios are refused."""
    if method != "pooled":
        raise ValueError("BAD_INPUT: only pooled moment ratios are supported")
    a, b = np.asarray(n, float), np.asarray(d, float)
    if (
        a.shape != b.shape
        or not a.size
        or not np.isfinite(a).all()
        or not np.isfinite(b).all()
    ):
        raise ValueError("BAD_INPUT")
    if not np.isscalar(kb) or not np.isfinite(kb) or kb <= 0:
        raise ValueError("BAD_INPUT")
    if abs(b.sum()) <= 1e-14 * max(float(abs(b).sum()), 1e-300):
        return {"status": "ZERO_SIGNAL", "temperature_K": None}
    if np.max(abs(a)) == 0:
        return {"status": "ZERO_SIGNAL", "temperature_K": None}
    return {"status": "DIAGNOSTIC_ONLY", "temperature_K": float(a.sum() / b.sum() / kb)}


def analyze(records, seed_ids, *, kb, dt=1.0):
    """Whole-seed delete-one inference, with retraining in every replicate.

    Records are made by compute_moments. No sampling reference is read here;
    reference qualification and physical acceptance remain with the caller.
    """
    if (
        not records
        or len(records) != len(seed_ids)
        or len(set(seed_ids)) != len(seed_ids)
        or not np.isscalar(kb)
        or not np.isfinite(kb)
        or kb <= 0
        or not np.isscalar(dt)
        or not np.isfinite(dt)
        or dt <= 0
    ):
        raise ValueError("BAD_INPUT")
    required = (
        "N",
        "D",
        "U_fixed",
        "shell_N",
        "shell_D",
        "shell_w",
        "SS",
        "FS",
        "shell_SS",
        "shell_AN",
        "ns",
        "force_norm",
    )
    schema = None
    for row in records:
        if any(key not in row for key in required):
            raise ValueError("BAD_INPUT")
        count = len(row["N"])
        if (
            count < 4
            or row["N"].shape != (count, len(FIXED))
            or row["D"].shape != row["N"].shape
        ):
            raise ValueError("BAD_INPUT")
        shape = {key: np.asarray(row[key]).shape[1:] for key in required}
        if schema is not None and shape != schema:
            raise ValueError("ID_MISMATCH: basis columns differ")
        schema = shape
        for key in required:
            val = np.asarray(row[key])
            if len(val) != count or not np.isfinite(val).all():
                raise ValueError("BAD_INPUT")
        if "shell_cuts" in row and not np.array_equal(
            row["shell_cuts"], records[0]["shell_cuts"]
        ):
            raise ValueError("ID_MISMATCH: shell identities differ")
    if all(np.max(row["force_norm"]) == 0 for row in records):
        return {
            "n_seeds": len(records),
            "seed_ids": list(seed_ids),
            "channels": {
                name: {
                    "status": "ZERO_SIGNAL",
                    "temperature_K": None,
                    "SE_K": None,
                    "confidence_set_K": [],
                }
                for name in FIXED + ["CF_GMM"]
            },
            "learned_status": "ZERO_SIGNAL",
            "weights": [],
            "analysis_wall_s": 0.0,
        }
    return summary(records, seed_ids, kb=kb, dt=dt)
