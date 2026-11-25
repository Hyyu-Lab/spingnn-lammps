"""Training-only geometry weights; evaluation uses frozen coefficients."""

import numpy as np

from .geometry import ORDERS


def learn_basis(records):
    q = (
        np.mean(
            [np.mean(r["SS"] / r["ns"][:, None, None], axis=0) for r in records], axis=0
        )
        - np.eye(3) / 3
    )
    d = np.mean([r["shell_D"].mean(0) for r in records], axis=0)
    w = np.mean([r["shell_w"].mean(0) for r in records], axis=0)
    factor = np.array([2 * l * (l + 1) for l in ORDERS])[:, None]
    weights = np.divide(
        d / factor, w[None, :], out=np.zeros_like(d), where=w[None, :] > 0
    )
    shell_ss = np.mean([r["shell_SS"].mean(0) for r in records], axis=0)
    shell_w = np.mean([r["shell_w"].mean(0) for r in records], axis=0)
    A = np.divide(
        shell_ss,
        shell_w[:, None, None],
        out=np.zeros_like(shell_ss),
        where=shell_w[:, None, None] > 0,
    )
    return {"q": q, "weights": weights, "A": A}


def evaluate_basis(r, basis):
    q = basis["q"]
    weights = basis["weights"]
    A = basis["A"]
    n_a = np.einsum("tsab,sab->ts", r["shell_AN"], A)
    d_a = 4 * np.einsum("tsab,sab->ts", r["shell_SS"], A)
    n = np.column_stack(
        [
            r["N"],
            np.einsum("tij,ij->t", r["FS"], q),
            np.einsum("tls,ls->tl", r["shell_N"], weights),
            n_a,
        ]
    )
    d = np.column_stack(
        [
            r["D"],
            3 * np.einsum("tij,ij->t", r["SS"], q) - r["ns"] * np.trace(q),
            np.einsum("tls,ls->tl", r["shell_D"], weights),
            d_a,
        ]
    )
    return n, d


def upper_bounds(r, basis):
    """Per-frame deterministic D bounds for the v1.2 numerical-zero gate."""
    q = basis["q"]
    weights = basis["weights"]
    A = basis["A"]
    np.mean(r["U_fixed"], axis=0)
    uq = np.full(len(r["ns"]), 3 * np.linalg.norm(q, 2) * float(np.mean(r["ns"])))
    ushell = np.sum(
        2
        * np.asarray(ORDERS)[:, None]
        * (np.asarray(ORDERS)[:, None] + 1)
        * np.abs(weights)
        * r["shell_w"][:, None, :],
        axis=2,
    )
    ua = 4 * np.linalg.norm(A, 2, axis=(1, 2))[None, :] * r["shell_w"]
    return np.column_stack([r["U_fixed"], uq, ushell, ua])


def learned_names(records):
    count = records[0]["shell_w"].shape[1]
    return (
        ["learn_Q2"]
        + ["learn_shell_l" + str(k) for k in ORDERS]
        + ["learn_A_shell_" + str(k) for k in range(count)]
    )
