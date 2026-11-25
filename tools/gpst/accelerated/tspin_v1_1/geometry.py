"""Distance-only probes and compact learned-basis sufficient moments."""

from itertools import product

import numpy as np
from scipy.spatial import cKDTree

from .native import edge_arrays

ORDERS = (1, 2, 4, 8)
FIXED = ["pair_l" + str(k) for k in ORDERS] + [
    "geometry_Q2",
    "geometry_zonal4",
    "geometry_chirality",
]


def legendre(x, max_l=8):
    x = np.asarray(x, dtype=float)
    if not np.isfinite(x).all() or np.any(abs(x) > 1 + 1e-10):
        raise ValueError("BAD_INPUT")
    x = np.clip(x, -1, 1)
    p, d = [np.ones_like(x), x], [np.zeros_like(x), np.ones_like(x)]
    for k in range(1, max_l):
        p.append(((2 * k + 1) * x * p[-1] - k * p[-2]) / (k + 1))
        d.append(((2 * k + 1) * (p[k] + x * d[-1]) - k * d[-2]) / (k + 1))
    return p, d


def minimum_image(delta, cell, pbc):
    """Exact closest image, including skew boxes; setup/oracle utility."""
    delta = np.asarray(delta, float)
    inv = np.linalg.inv(cell)
    flat = delta.reshape(-1, 3)
    frac = flat @ inv
    base = frac - np.rint(frac) * pbc
    out = base @ cell
    norms = np.linalg.norm(out, axis=1)
    widths = np.linalg.norm(inv, axis=0)
    for k, u in enumerate(base):
        bound = norms[k] * widths + 1e-12
        ranges = [
            range(int(np.ceil(u[j] - bound[j])), int(np.floor(u[j] + bound[j])) + 1)
            if pbc[j]
            else [0]
            for j in range(3)
        ]
        if np.prod([len(r) for r in ranges]) > 100000:
            raise ValueError("BAD_INPUT: box image search too large")
        shifts = np.array(list(product(*ranges)))
        candidates = (u - shifts) @ cell
        out[k] = candidates[np.argmin(np.sum(candidates * candidates, axis=1))]
    return out.reshape(delta.shape)


def nearest_distance(pos, cell, pbc):
    """Nearest distinct magnetic-site spacing; no species or length presets."""
    if len(pos) < 2:
        raise ValueError("BAD_INPUT: no neighbor for relative cutoff")
    inv = np.linalg.inv(cell)
    frac = pos @ inv
    frac -= np.floor(frac) * pbc
    wrapped = frac @ cell
    guess = float(np.min(cKDTree(wrapped).query(wrapped, k=2)[0][:, 1]))
    if guess <= 1e-12:
        raise ValueError("BAD_INPUT: coincident sites")
    widths = np.linalg.norm(inv, axis=0)
    reach = np.ceil(guess * widths).astype(int) + 1
    ranges = [
        range(-int(reach[j]), int(reach[j]) + 1) if pbc[j] else [0] for j in range(3)
    ]
    if np.prod([len(r) for r in ranges]) * len(pos) > 10000000:
        raise ValueError("BAD_INPUT: nearest search too large")
    shifts = np.array(list(product(*ranges))) @ cell
    cloud = (wrapped[None] + shifts[:, None]).reshape(-1, 3)
    dist, _ = cKDTree(cloud).query(wrapped, k=2)
    nearest = float(np.min(dist[:, 1]))
    if nearest <= 1e-12:
        raise ValueError("BAD_INPUT: coincident sites")
    return nearest


def cutoff_radius(pos, box, periodic, cutoff):
    """Resolve a caller cutoff for both new and prepared geometry paths."""
    if isinstance(cutoff, dict):
        if set(cutoff) != {"nearest_multiplier"}:
            raise ValueError("BAD_INPUT")
        mult = cutoff["nearest_multiplier"]
        if not np.isscalar(mult) or not np.isfinite(mult) or mult <= 0:
            raise ValueError("BAD_INPUT")
        rc = float(mult) * nearest_distance(pos, box, periodic)
    else:
        rc = float(cutoff)
    if not np.isfinite(rc) or rc <= 0:
        raise ValueError("BAD_INPUT")
    return rc


def build_graph(positions, cell, pbc, cutoff, *, shell_cuts=None):
    """Expected-linear linked cells; unique periodic images are required."""
    pos, box = np.asarray(positions, float), np.asarray(cell, float)
    periodic = np.asarray(pbc)
    if (
        pos.ndim != 2
        or pos.shape[1] != 3
        or box.shape != (3, 3)
        or periodic.shape != (3,)
        or periodic.dtype.kind != "b"
        or not np.isfinite(pos).all()
        or not np.isfinite(box).all()
        or abs(np.linalg.det(box)) < 1e-12
        or np.linalg.cond(box) > 1e6
    ):
        raise ValueError("BAD_INPUT")
    inv = np.linalg.inv(box)
    rc = cutoff_radius(pos, box, periodic, cutoff)
    height = 1 / np.linalg.norm(inv, axis=0)
    if np.any(periodic & (2 * rc >= height * (1 - 1e-10))):
        raise ValueError("BAD_INPUT: cutoff exceeds unique image range")
    frac = pos @ inv
    frac -= np.floor(frac) * periodic
    width = rc / height
    counts = np.maximum(1, np.floor(1 / width).astype(int))
    binsize = np.where(periodic, 1 / counts, width)
    keys = np.floor(frac / binsize).astype(int)
    keys[:, periodic] %= counts[periodic]
    bins = {}
    for idx, key in enumerate(keys):
        bins.setdefault(tuple(key), []).append(idx)
    left, right = [], []
    for key, vals in bins.items():
        around = set()
        for offset in product((-1, 0, 1), repeat=3):
            target = np.array(key) + offset
            target[periodic] %= counts[periodic]
            around.add(tuple(target))
        for other in around:
            if other not in bins or other < key:
                continue
            a, b = np.array(vals), np.array(bins[other])
            if other == key:
                u, v = np.triu_indices(len(a), 1)
                left.extend(a[u])
                right.extend(a[v])
            else:
                left.extend(np.repeat(a, len(b)))
                right.extend(np.tile(b, len(a)))
    a, b = np.asarray(left, int), np.asarray(right, int)
    i, j = np.minimum(a, b), np.maximum(a, b)
    delta = frac[j] - frac[i]
    # Inside rc, every periodic fractional difference must be in (-1/2,1/2).
    delta -= np.rint(delta) * periodic
    delta = delta @ box
    dist = np.linalg.norm(delta, axis=1)
    if np.any(dist <= 1e-12):
        raise ValueError("BAD_INPUT: coincident sites")
    use = dist < rc
    i, j, delta, dist = i[use], j[use], delta[use], dist[use]
    order = np.lexsort((j, i))
    i, j, delta, dist = i[order], j[order], delta[order], dist[order]
    weight = (1 - (dist / rc) ** 2) ** 3
    direction = delta / dist[:, None]
    if shell_cuts is None:
        lengths = np.sort(dist)
        if len(lengths):
            split = np.flatnonzero(np.diff(lengths) > 0.01 * lengths[0])
            shell_cuts = (lengths[split] + lengths[split + 1]) / 2
        else:
            shell_cuts = np.empty(0)
    shell_cuts = np.asarray(shell_cuts, float)
    if (
        shell_cuts.ndim != 1
        or not np.isfinite(shell_cuts).all()
        or np.any(shell_cuts <= 0)
        or np.any(np.diff(shell_cuts) <= 0)
    ):
        raise ValueError("BAD_INPUT: shell boundaries")
    group = np.searchsorted(shell_cuts, dist)
    group_order = np.argsort(group, kind="stable")
    occupied, starts = np.unique(group[group_order], return_index=True)
    q = np.zeros((len(pos), 3, 3))
    outer = weight[:, None, None] * direction[:, :, None] * direction[:, None, :]
    np.add.at(q, i, outer)
    np.add.at(q, j, outer)
    return {
        "i": i,
        "j": j,
        "w": weight,
        "direction": direction,
        "distance": dist,
        "shell": group,
        "shell_cuts": shell_cuts,
        "Q": q,
        "rc": rc,
        "n_sites": len(pos),
        "group_order": group_order,
        "occupied": occupied,
        "starts": starts,
        "positions": pos.copy(),
        "cell": box.copy(),
        "pbc": periodic.copy(),
    }


def grouped(values, graph, axis):
    """One edge pass, including empty shells; no shell-by-edge nested scan."""
    shape = list(values.shape)
    shape[axis] = len(graph["shell_cuts"]) + 1
    out = np.zeros(shape)
    if len(graph["i"]):
        reduced = np.add.reduceat(
            np.take(values, graph["group_order"], axis=axis), graph["starts"], axis=axis
        )
        slices = [slice(None)] * values.ndim
        slices[axis] = graph["occupied"]
        out[tuple(slices)] = reduced
    return out


def moment_arrays(s, f, graph):
    """All fixed channels and training summaries, O(N + edges*max_order)."""
    s, f = np.asarray(s, float), np.asarray(f, float)
    if (
        s.ndim != 3
        or s.shape != f.shape
        or not np.isfinite(s).all()
        or not np.isfinite(f).all()
    ):
        raise ValueError("BAD_INPUT")
    if np.max(abs(np.sum(s * s, axis=-1) - 1), initial=0) > 1e-10:
        raise ValueError("NONUNIT_SPIN")
    ft = f - s * np.sum(s * f, axis=-1, keepdims=True)
    t, ns, _ = s.shape
    if ns != graph["n_sites"]:
        raise ValueError("ID_MISMATCH")
    i, j, w = graph["i"], graph["j"], graph["w"]
    u, _sh = graph["direction"], graph["shell"]
    groups = len(graph["shell_cuts"]) + 1
    fused = edge_arrays(s, ft, graph)
    if fused is None:
        si, sj = s[:, i], s[:, j]
        fi, fj = ft[:, i], ft[:, j]
        x = np.sum(si * sj, axis=2)
        c = np.sum(sj * fi + si * fj, axis=2)
        ps, ds = legendre(x)
        sn = np.stack([grouped(w * ds[k] * c, graph, 1) for k in ORDERS], axis=1)
        sd = np.stack(
            [2 * k * (k + 1) * grouped(w * ps[k], graph, 1) for k in ORDERS], axis=1
        )
    else:
        raw_n, raw_d, zn, zd, ph, cn, raw_ss, raw_an = fused
        sn = grouped(raw_n, graph, 2)
        factor = np.array([2 * k * (k + 1) for k in ORDERS])
        sd = factor[None, :, None] * grouped(raw_d, graph, 2)
    sw = np.broadcast_to(grouped(w, graph, 0), (t, groups)).copy()
    n = np.zeros((t, len(FIXED)))
    d = np.zeros_like(n)
    n[:, :4] = sn.sum(2)
    d[:, :4] = sd.sum(2)
    q = graph["Q"]
    qs = np.einsum("nij,tnj->tni", q, s)
    n[:, 4] = np.sum(qs * ft, axis=(1, 2))
    d[:, 4] = 3 * np.sum(qs * s, axis=(1, 2)) - np.trace(q, axis1=1, axis2=2).sum()
    if fused is None:
        for directions, forces in ((si, fi), (sj, fj)):
            z = np.einsum("tej,ej->te", directions, u)
            p, dp = legendre(z, 4)
            n[:, 5] += np.sum(w * dp[4] * np.einsum("tej,ej->te", forces, u), axis=1)
            d[:, 5] += 20 * np.sum(w * p[4], axis=1)
        phi = np.sum(w[None] * np.sum(u[None] * np.cross(si, sj), axis=2), axis=1)
        n[:, 6] = np.sum(
            w[None]
            * (
                np.sum(fi * np.cross(sj, u[None]), axis=2)
                + np.sum(fj * np.cross(u[None], si), axis=2)
            ),
            axis=1,
        )
        d[:, 6] = 4 * phi
    else:
        for endpoint in (0, 1):
            n[:, 5] += np.sum(zn[:, endpoint], axis=1)
            d[:, 5] += 20 * np.sum(zd[:, endpoint], axis=1)
        n[:, 6] = np.sum(cn, axis=1)
        d[:, 6] = 4 * np.sum(ph, axis=1)
    ss = np.einsum("tna,tnb->tab", s, s)
    fs = np.einsum("tna,tnb->tab", ft, s)
    if fused is None:
        sss = grouped(
            w[None, :, None, None] * si[:, :, :, None] * sj[:, :, None, :], graph, 1
        )
        san = grouped(
            w[None, :, None, None]
            * (fi[:, :, :, None] * sj[:, :, None, :] + si[:, :, :, None] * fj[:, :, None, :]),
            graph,
            1,
        )
    else:
        sss = grouped(raw_ss, graph, 1)
        san = grouped(raw_an, graph, 1)
    total_w = float(w.sum())
    trq = float(np.trace(q, axis1=1, axis2=2).sum())
    bound = np.column_stack(
        [np.full(t, 2 * k * (k + 1) * total_w) for k in ORDERS]
        + [np.full(t, 4 * trq), np.full(t, 40 * total_w), np.full(t, 4 * total_w)]
    )
    return {
        "N": n,
        "D": d,
        "shell_N": sn,
        "shell_D": sd,
        "shell_w": sw,
        "SS": ss,
        "FS": fs,
        "shell_SS": sss,
        "shell_AN": san,
        "ns": np.full(t, ns),
        "U_fixed": bound,
        "force_norm": np.sum(ft * ft, axis=(1, 2)),
    }
