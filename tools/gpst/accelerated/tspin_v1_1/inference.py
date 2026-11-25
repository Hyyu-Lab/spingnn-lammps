"""Seed-level Fieller sets and fully re-trained cross-fitted moment inference."""

import time

import numpy as np
from scipy.stats import f as fdist
from scipy.stats import t as student

from .basis import evaluate_basis, learn_basis, upper_bounds
from .basis import learned_names as basis_names
from .geometry import FIXED

# IID zero-null calibration (10,000 fixed-seed simulations, v1.2 freeze).
# The statistic is max_k(range(SE_8,SE_16,SE_32)/max(SE_k)); it is used only
# for single-record replay where block independence is part of the SE gate.
IID_RANGE95_BY_N = {50: 0.458, 512: 0.426, 2000: 0.429}


def covariance_jackknife(x):
    x = np.asarray(x)
    z = x - x.mean(0)
    return (len(x) - 1) / len(x) * (z.T @ z)


def pooled_pair(nlist, dlist):
    count = sum(len(n) for n in nlist)
    return np.stack(
        [sum(n.sum(0) for n in nlist) / count, sum(d.sum(0) for d in dlist) / count]
    )


def fieller(a, b, cov, df, alpha=0.05):
    q = float(student.ppf(1 - alpha / 2, df))
    q2 = q * q
    aa = b * b - q2 * cov[1, 1]
    bb = -2 * a * b + 2 * q2 * cov[0, 1]
    cc = a * a - q2 * cov[0, 0]
    scale = max(abs(aa), abs(bb), abs(cc), 1e-300)
    eps = 1e-13 * scale
    if abs(aa) < eps:
        if abs(bb) < eps:
            confidence_sets = [[None, None]] if cc <= eps else []
        else:
            root = -cc / bb
            confidence_sets = [[None, float(root)]] if bb > 0 else [[float(root), None]]
    else:
        disc = bb * bb - 4 * aa * cc
        if disc < -1e-13 * max(bb * bb, abs(4 * aa * cc), 1e-300):
            confidence_sets = [[None, None]] if aa < 0 else []
        else:
            disc = max(0.0, disc)
            roots = sorted(
                [(-bb - np.sqrt(disc)) / (2 * aa), (-bb + np.sqrt(disc)) / (2 * aa)]
            )
            confidence_sets = (
                [roots] if aa > 0 else [[None, roots[0]], [roots[1], None]]
            )
    zero = abs(b) <= q * np.sqrt(max(0.0, cov[1, 1]))
    finite = len(confidence_sets) == 1 and all(
        v is not None for v in confidence_sets[0]
    )
    if abs(a) < 1e-13 and abs(b) < 1e-13 and np.max(abs(cov)) < 1e-24:
        status = "ZERO_SIGNAL"
    elif not confidence_sets:
        status = "EMPTY_CI"
    elif zero:
        status = "WEAK_DENOMINATOR"
    elif not finite:
        status = "UNBOUNDED_CI"
    else:
        status = "IDENTIFIED"
    point = float(a / b) if abs(b) > 1e-30 else None
    if point is not None and point < 0 and status == "IDENTIFIED":
        status = "MOMENT_MODEL_MISFIT"
    var = (
        None
        if point is None
        else (cov[0, 0] - 2 * point * cov[0, 1] + point * point * cov[1, 1]) / (b * b)
    )
    return {
        "status": status,
        "point": point,
        "SE": float(np.sqrt(max(0, var))) if var is not None else None,
        "confidence_set": confidence_sets,
        "mean_N": float(a),
        "mean_D": float(b),
        "covariance": np.asarray(cov).tolist(),
        "df": int(df),
        "critical": q,
        "denominator_ci": [
            float(b - q * np.sqrt(max(0, cov[1, 1]))),
            float(b + q * np.sqrt(max(0, cov[1, 1]))),
        ],
    }


def contains(result, value):
    return any(
        (lo is None or lo <= value) and (hi is None or value <= hi)
        for lo, hi in result["confidence_set"]
    )


def geyer(x, dt=1.0):
    x = np.asarray(x, float) - np.mean(x)
    n = len(x)
    if n < 4 or np.dot(x, x) < 1e-28:
        return {"tau": None, "status": "ZERO_VARIANCE"}
    size = 1 << (2 * n - 1).bit_length()
    v = np.fft.rfft(x, n=size)
    ac = np.fft.irfft(v * v.conj(), n=size)[:n] / n
    ac /= ac[0]
    pairs = ac[: 2 * (n // 2)].reshape(-1, 2).sum(1)
    good = []
    prev = float("inf")
    for z in pairs:
        if z <= 0:
            break
        prev = min(prev, z)
        good.append(prev)
    tau = max(0.5, -0.5 + sum(good)) * dt
    return {
        "tau": float(tau),
        "status": "RESOLUTION_LIMITED" if tau <= 0.50000001 * dt else "ESTIMATED",
        "truncated": len(good) == len(pairs),
        "span_over_tau": n * dt / tau,
    }


def blocks(x, k):
    return np.array([a.mean(0) for a in np.array_split(x, min(k, len(x)))])


def fit_weight(
    records, *, snr_threshold=5.0, denominator_bounds=None, zero_fraction=1e-10
):
    """Select identifiable channels and fit a training-fold GMM.

    The denominator SNR is computed from independent training-seed means.  The
    selection is deliberately inside this function so every outer delete-one
    replicate reruns it and no evaluation-fold information can enter.
    """
    meansN = np.array([r[0].mean(0) for r in records])
    meansD = np.array([r[1].mean(0) for r in records])
    p = meansD.shape[1]
    dn = meansD.mean(0)
    if len(records) < 2:
        raise ValueError("RANK_DEFICIENT: need independent training seeds")
    sdd = meansD.std(0, ddof=1) / np.sqrt(len(records))
    snr = np.divide(abs(dn), sdd, out=np.zeros(p), where=sdd > 1e-30)
    snr[(sdd <= 1e-30) & (abs(dn) > 1e-20)] = np.inf
    if denominator_bounds is None:
        denominator_bounds = np.ones(p)
    zero_limit = zero_fraction * np.maximum(np.asarray(denominator_bounds), 1e-300)
    numerical_zero = abs(dn) < zero_limit
    selected = np.flatnonzero((snr >= snr_threshold) & ~numerical_zero)
    rejected = np.flatnonzero(~((snr >= snr_threshold) & ~numerical_zero))
    if not len(selected):
        raise ValueError("RANK_DEFICIENT: no channel passed denominator SNR threshold")
    c = np.zeros(p)
    if len(selected) == 1:
        c[selected[0]] = 1.0
        scale = float(c[selected[0]] * dn[selected[0]])
        if abs(scale) < 1e-300:
            raise ValueError("RANK_DEFICIENT: zero training denominator scale")
        c[selected[0]] /= scale
        return c, {
            "rank": 1,
            "selection_threshold": snr_threshold,
            "zero_fraction": zero_fraction,
            "selected": selected.tolist(),
            "rejected": rejected.tolist(),
            "snr": [None if not np.isfinite(z) else float(z) for z in snr],
            "numerical_zero": numerical_zero.tolist(),
            "direct_channel": int(selected[0]),
            "normalized_scale": scale,
        }
    # Pilot from the strongest selected denominator, followed by a shrinkage
    # inverse covariance on block residuals for the selected channels only.
    idx = int(selected[np.argmax(snr[selected])])
    energy_scale = float(meansN[:, idx].mean() / dn[idx])
    rr = np.array([blocks((n - energy_scale * d)[:, selected], 16) for n, d in records])
    flat = rr.reshape(-1, len(selected))
    cov = np.atleast_2d(np.cov(flat, rowvar=False, ddof=1))
    if cov.shape != (len(selected), len(selected)):
        cov = np.diag(np.full(len(selected), float(np.var(flat))))
    scale_vec = np.sqrt(np.maximum(0, np.diag(cov)))
    active = scale_vec > max(scale_vec.max() * 1e-12, 1e-30)
    if not active.any():
        raise ValueError("RANK_DEFICIENT: covariance")
    std = cov[np.ix_(active, active)] / np.outer(scale_vec[active], scale_vec[active])
    v, u = np.linalg.eigh(std)
    keep = np.flatnonzero(v > max(v.max() * 1e-6, 1e-14))
    if not len(keep):
        raise ValueError("RANK_DEFICIENT: no kept covariance direction")
    null = v < max(1e-14, v.max() * 1e-12)
    residual = (meansN.mean(0) - energy_scale * dn)[selected][active] / scale_vec[
        active
    ]
    if null.any() and np.max(abs(u[:, null].T @ residual)) > 1e-10:
        raise ValueError(
            "MOMENT_MODEL_MISFIT: exact covariance null-space contradiction"
        )
    vals = 0.9 * v[keep] + 0.1 * np.mean(v[keep])
    win = (u[:, keep] / vals) @ u[:, keep].T
    w = np.zeros_like(cov)
    w[np.ix_(active, active)] = win / np.outer(scale_vec[active], scale_vec[active])
    c[selected] = w @ dn[selected]
    scale = float(c @ dn)
    if abs(scale) < 1e-300:
        raise ValueError("RANK_DEFICIENT: zero training denominator scale")
    c /= scale
    return c, {
        "rank": len(keep),
        "pilot": energy_scale,
        "pilot_channel": idx,
        "selection_threshold": snr_threshold,
        "zero_fraction": zero_fraction,
        "selected": selected.tolist(),
        "rejected": rejected.tolist(),
        "snr": [None if not np.isfinite(z) else float(z) for z in snr],
        "numerical_zero": numerical_zero.tolist(),
        "scale": scale_vec.tolist(),
        "normalized_scale": scale,
        "null_conflict": False,
    }


def crossfit(records, labels, return_bounds=False):
    nout = [None] * len(records)
    dout = [None] * len(records)
    uout = [None] * len(records)
    info = []
    for fold in (0, 1):
        train = [i for i, l in enumerate(labels) if l != fold]
        test = [i for i, l in enumerate(labels) if l == fold]
        if len(train) < 2 or not test:
            raise ValueError("INSUFFICIENT_SEEDS")
        basis = learn_basis([records[i] for i in train])
        tr = [evaluate_basis(records[i], basis) for i in train]
        ub = np.mean(
            np.concatenate([upper_bounds(records[i], basis) for i in train], axis=0),
            axis=0,
        )
        c, fit_info = fit_weight(tr, denominator_bounds=ub)
        fit_info.update(
            fold=fold,
            Q=basis["q"].tolist(),
            shell_weights=basis["weights"].tolist(),
            A_shell=basis["A"].tolist(),
        )
        info.append(fit_info)
        for i in test:
            n, d = evaluate_basis(records[i], basis)
            nout[i] = np.column_stack([n, n @ c])
            dout[i] = np.column_stack([d, d @ c])
            u = upper_bounds(records[i], basis)
            uout[i] = np.column_stack([u, u @ np.abs(c)])
    return (nout, dout, info, uout) if return_bounds else (nout, dout, info)


def summary(records, seed_ids, *, kb, dt=1.0, reference_names=None):
    started = time.perf_counter()
    m = len(records)
    labels = np.arange(m) % 2
    NAMES = FIXED + basis_names(records)
    LEARNED = basis_names(records)
    cf = m >= 6
    names = NAMES + ["CF_GMM"] if cf else FIXED
    if cf:
        try:
            nlist, dlist, weight_info, ulist = crossfit(
                records, labels, return_bounds=True
            )
        except ValueError as exc:
            names = NAMES + ["CF_GMM"]
            channels = {
                name: {"status": "RANK_DEFICIENT", "reason": str(exc)} for name in names
            }
            return {
                "n_seeds": m,
                "seed_ids": seed_ids,
                "channels": channels,
                "comparisons": {},
                "weights": [],
                "joint_diagnostic": {"status": "RANK_DEFICIENT", "reason": str(exc)},
                "uncertainty": "training-fold channel selection failed",
                "learned_status": "RANK_DEFICIENT",
                "analysis_wall_s": time.perf_counter() - started,
            }
    else:
        nlist = [r["N"] for r in records]
        dlist = [r["D"] for r in records]
        weight_info = []
        ulist = [r.get("U_fixed") for r in records]
    # Historical comparisons are appended and share exactly the selected frames.
    if reference_names:
        names = names + reference_names
        nlist = [np.column_stack([n, r["reference_N"]]) for n, r in zip(nlist, records)]
        dlist = [np.column_stack([d, r["reference_D"]]) for d, r in zip(dlist, records)]
    means = np.array([[n.mean(0), d.mean(0)] for n, d in zip(nlist, dlist)])
    point = pooled_pair(nlist, dlist)
    leaves = []
    leaf_failure = None
    leaf_attempted = 0
    if m > 1:
        for drop in range(m):
            use = [i for i in range(m) if i != drop]
            leaf_attempted += 1
            try:
                if cf:
                    # return_bounds=False has exactly three return values.
                    # Keep this call separate from the full fit so a failed
                    # retraining replicate cannot silently become a different
                    # covariance estimate.
                    nn, dd, _ = crossfit([records[i] for i in use], labels[use])
                    ls = pooled_pair(nn, dd)
                    if reference_names:
                        refs = np.mean(
                            [
                                [
                                    records[i]["reference_N"].mean(0),
                                    records[i]["reference_D"].mean(0),
                                ]
                                for i in use
                            ],
                            axis=0,
                        )
                        ls = np.column_stack([ls, refs])
                else:
                    ls = pooled_pair([nlist[i] for i in use], [dlist[i] for i in use])
                leaves.append(ls)
            except ValueError as exc:
                leaf_failure = str(exc)
                break
    else:
        bn = blocks(nlist[0], 10)
        bd = blocks(dlist[0], 10)
        bm = np.stack([bn, bd], axis=1)
        point = bm.mean(0)
        leaves = [bm[np.arange(len(bm)) != i].mean(0) for i in range(len(bm))]
    df = m - 1 if m > 1 else 9
    channels = {}
    leafarr = np.asarray(leaves) if not leaf_failure else None
    learned_names = set(LEARNED) | {"CF_GMM"}
    for j, name in enumerate(names):
        if leaf_failure and name in learned_names:
            channels[name] = {
                "status": "JACKKNIFE_FAILED",
                "reason": leaf_failure,
                "jackknife_replicates_attempted": leaf_attempted,
                "jackknife_replicates_completed": len(leaves),
                "uncertainty": "whole-seed delete-one with retraining failed",
            }
            continue
        if leaf_failure:
            cov = np.cov(means[:, :, j], rowvar=False, ddof=1) / m
        else:
            cov = covariance_jackknife(leafarr[:, :, j])
        out = fieller(*point[:, j], cov, df)
        out["temperature_K"] = out["point"] / kb if out["point"] is not None else None
        out["SE_K"] = out["SE"] / kb if out["SE"] is not None else None
        out["confidence_set_K"] = [
            [None if a is None else a / kb, None if b is None else b / kb]
            for a, b in out["confidence_set"]
        ]
        out["uncertainty"] = (
            "seed-mean covariance after retraining failure"
            if leaf_failure
            else (
                "whole-seed delete-one with retraining"
                if m > 1
                else "10 contiguous blocks; no independent training"
            )
        )
        if out["point"] is not None:
            ac = [
                geyer(n[:, j] - out["point"] * d[:, j], dt)
                for n, d in zip(nlist, dlist)
            ]
            out["residual_acf"] = ac
            checks = []
            for n, d, aa in zip(nlist, dlist, ac):
                r = n[:, j] - out["point"] * d[:, j]
                ses = [
                    float(np.std(blocks(r, k), ddof=1) / np.sqrt(min(k, len(r))))
                    for k in (8, 16, 32)
                ]
                spread = (max(ses) - min(ses)) / max(ses) if max(ses) > 0 else None
                n_obs = len(r)
                if m == 1:
                    if n_obs not in IID_RANGE95_BY_N:
                        raise ValueError(f"UNCALIBRATED_SINGLE_RECORD_N:{n_obs}")
                    thresh = IID_RANGE95_BY_N[n_obs]
                    calibrated = True
                else:
                    thresh = IID_RANGE95_BY_N.get(n_obs)
                    calibrated = thresh is not None
                passed = bool(
                    calibrated
                    and spread is not None
                    and spread <= thresh
                    and aa["tau"] is not None
                    and n_obs / 16 * dt >= 5 * aa["tau"]
                    and aa.get("span_over_tau", 0) >= 50
                )
                checks.append(
                    {
                        "block_SE": ses,
                        "relative_range": spread,
                        "threshold": thresh,
                        "calibrated": calibrated,
                        "diagnostic_status": "PASS"
                        if passed
                        else ("UNCALIBRATED" if not calibrated else "FAIL"),
                        "passed": passed,
                    }
                )
            out["block_diagnostics"] = checks
            out["sampling_pass_count"] = sum(c["passed"] for c in checks)
            out["sampling_diagnostic"] = (
                "DIAGNOSTIC"
                if m >= 2
                else (
                    "PASS"
                    if all(c["passed"] for c in checks)
                    else "INSUFFICIENT_BLOCKS"
                )
            )
        # The scale-free bound is evaluated on the same full group. A tiny
        # denominator is kept only as a diagnostic and never converted to T.
        if ulist and all(u is not None for u in ulist) and j < ulist[0].shape[1]:
            ub = float(np.mean([u[:, j].mean() for u in ulist]))
            db = float(abs(point[1, j]))
            out["denominator_bound_mean"] = ub
            if ub > 0 and db < 1e-10 * ub:
                out["status"] = "NUMERICAL_NOISE_UNCERTIFIED"
                out["temperature_K"] = None
                out["SE_K"] = None
                out["confidence_set_K"] = []
        channels[name] = out
    comparisons = {}
    if reference_names and leafarr is not None:
        for name in names:
            if name in reference_names:
                continue
            for ref in reference_names:
                j = names.index(name)
                k = names.index(ref)
                o = channels[name]
                b = channels[ref]
                key = name + " vs " + ref
                if o.get("status") != "IDENTIFIED" or b.get("status") != "IDENTIFIED":
                    comparisons[key] = {
                        "status": "INCONCLUSIVE",
                        "reason": "UNIDENTIFIED",
                    }
                    continue
                ratios = leafarr[:, 0, :] / np.where(
                    abs(leafarr[:, 1, :]) > 1e-30, leafarr[:, 1, :], np.nan
                )
                delta = (ratios[:, j] - ratios[:, k]) / kb
                if not np.isfinite(delta).all():
                    comparisons[key] = {
                        "status": "INCONCLUSIVE",
                        "reason": "DELETE_ONE_DENOMINATOR",
                    }
                    continue
                value = o["temperature_K"] - b["temperature_K"]
                se = float(
                    np.sqrt(
                        (len(delta) - 1)
                        / len(delta)
                        * np.sum((delta - delta.mean()) ** 2)
                    )
                )
                half = float(student.ppf(0.975, df) * se)
                bonf = float(student.ppf(1 - 0.05 / 22, df) * se)
                status = "CONSISTENT" if abs(value) <= half else "RESOLVED_DIFFERENCE"
                if m == 1 and (
                    o.get("sampling_diagnostic") != "PASS"
                    or b.get("sampling_diagnostic") != "PASS"
                ):
                    status = "INCONCLUSIVE"
                comparisons[key] = {
                    "status": status,
                    "delta_K": value,
                    "SE_K": se,
                    "ci95_K": [value - half, value + half],
                    "bonferroni11_ci_K": [value - bonf, value + bonf],
                    "precision": "PRECISION_SUFFICIENT"
                    if half <= 0.01 * abs(b["temperature_K"])
                    else "WEAK_PRECISION",
                    "nominal_difference_detected": bool(abs(value) > half),
                }
    # Fixed-channel common-energy_scale diagnostic. Exact null-space discrepancies are
    # checked before the kept covariance subspace is inverted.
    joint = {"status": "INSUFFICIENT_SEEDS"}
    target = channels.get("CF_GMM", channels["pair_l1"]).get("point")
    if m > 4 and target is not None:
        residual = means[:, 0, :6] - target * means[:, 1, :6]
        mu = residual.mean(0)
        cov = np.cov(residual, rowvar=False, ddof=1) / m
        scales = np.sqrt(np.maximum(np.diag(cov), 1e-300))
        sc = cov / np.outer(scales, scales)
        v, u = np.linalg.eigh(sc)
        z = mu / scales
        null = v < v.max() * 1e-12
        conflict = bool(null.any() and np.max(abs(u[:, null].T @ z)) > 1e-10)
        keep = np.flatnonzero(v > v.max() * 1e-6)[-min(4, m - 2) :]
        rank = len(keep)
        if rank:
            q = float(np.sum((u[:, keep].T @ z) ** 2 / v[keep]))
            critical = float(
                rank * (m - 1) / (m - rank) * fdist.ppf(0.95, rank, m - rank)
            )
            joint = {
                "status": "MOMENT_MODEL_MISFIT"
                if conflict or q > critical
                else "NO_MISFIT_DETECTED",
                "Q": q,
                "rank": rank,
                "critical": critical,
                "null_conflict": conflict,
                "note": "Fitted-energy_scale diagnostic, not a calibrated fixed-energy_scale confidence set; primary confidence_set is Fieller.",
            }
    if leaf_failure:
        uncertainty = (
            "whole-seed delete-one with retraining failed for CF/learned "
            "channels; fixed channels use seed-mean covariance"
        )
    elif m > 1:
        uncertainty = "whole-seed delete-one with retraining"
    else:
        uncertainty = "10 contiguous blocks; no independent training"
    return {
        "n_seeds": m,
        "seed_ids": seed_ids,
        "channels": channels,
        "comparisons": comparisons,
        "weights": weight_info,
        "joint_diagnostic": joint,
        "uncertainty": uncertainty,
        "learned_status": "EVALUATED" if cf else "INSUFFICIENT_SEEDS",
        "jackknife_replicates_attempted": leaf_attempted,
        "jackknife_replicates_completed": len(leaves),
        "analysis_wall_s": time.perf_counter() - started,
    }
