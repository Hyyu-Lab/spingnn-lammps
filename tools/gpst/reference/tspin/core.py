"""Array contract and first-order moment API; no model evaluation."""

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from .geometry import FIXED, ORDERS, build_graph, cutoff_radius, moment_arrays

CHANNELS = tuple(FIXED)


class FrameError(ValueError):
    """Stable input rejection code."""


@dataclass(frozen=True)
class FrameInput:
    atom_ids: np.ndarray
    magnetic_mask: np.ndarray
    positions: np.ndarray
    cell: np.ndarray
    pbc: np.ndarray
    spins: np.ndarray
    spin_forces: np.ndarray
    frame_index: np.ndarray
    force_contract_confirmed: bool


def _real(value):
    a = np.asarray(value)
    if a.dtype.kind not in "fiu" or not np.isfinite(a).all():
        raise FrameError("BAD_INPUT")
    return a.astype(np.float64, copy=False)


def validate_frame(
    *,
    atom_ids,
    magnetic_mask,
    positions,
    cell,
    pbc,
    spins,
    spin_forces,
    frame_index,
    force_contract_confirmed,
):
    """Confirm eV-valued -dH/ds for unit directions, with aligned accepted states."""
    if force_contract_confirmed is not True:
        raise FrameError("FORCE_CONTRACT_MISSING")
    ids, mask, periodic, steps = map(
        np.asarray, (atom_ids, magnetic_mask, pbc, frame_index)
    )
    for a in (ids, mask, periodic, steps):
        if a.dtype.kind in "fciu" and not np.isfinite(a).all():
            raise FrameError("BAD_INPUT")
    if (
        ids.ndim != 1
        or ids.size == 0
        or ids.dtype.kind not in "iu"
        or np.any(ids <= 0)
        or len(np.unique(ids)) != len(ids)
    ):
        raise FrameError("ID_MISMATCH")
    if mask.shape != ids.shape or mask.dtype.kind != "b" or not mask.any():
        raise FrameError("ID_MISMATCH")
    pos, box, s, f = map(_real, (positions, cell, spins, spin_forces))
    if pos.ndim != 3 or pos.shape[0] == 0 or pos.shape[1:] != (len(ids), 3):
        raise FrameError("ID_MISMATCH")
    count = len(pos)
    if (
        box.shape != (count, 3, 3)
        or s.shape != (count, int(mask.sum()), 3)
        or s.shape != f.shape
        or periodic.shape != (3,)
        or periodic.dtype.kind != "b"
        or steps.shape != (count,)
        or steps.dtype.kind not in "iu"
        or np.any(steps[1:] <= steps[:-1])
    ):
        raise FrameError("ID_MISMATCH")
    if np.max(abs(np.sum(s * s, axis=-1) - 1)) > 1e-10:
        raise FrameError("NONUNIT_SPIN")
    if np.any(abs(np.linalg.det(box)) < 1e-12) or np.any(np.linalg.cond(box) > 1e6):
        raise FrameError("BAD_INPUT")
    atom_order = np.argsort(ids)
    spin_order = np.argsort(ids[mask])
    return FrameInput(
        ids[atom_order],
        mask[atom_order],
        pos[:, atom_order],
        box,
        periodic,
        s[:, spin_order],
        f[:, spin_order],
        steps,
        True,
    )


def load_frame_npz(path: str | Path, *, force_contract_confirmed):
    """Read only the common NPZ arrays. No sidecar file is opened."""
    keys = (
        "atom_ids",
        "magnetic_mask",
        "positions",
        "cell",
        "pbc",
        "spins",
        "spin_forces",
        "frame_index",
    )
    with np.load(path, allow_pickle=False) as data:
        if any(k not in data for k in keys):
            raise FrameError("ID_MISMATCH")
        args = {k: data[k] for k in keys}
    return validate_frame(**args, force_contract_confirmed=force_contract_confirmed)


def from_lammps_fm(fm, hbar_eVps):
    """Caller-confirmed angular-frequency conversion; no default constant."""
    a = _real(fm)
    if not np.isscalar(hbar_eVps) or not np.isfinite(hbar_eVps) or hbar_eVps <= 0:
        raise FrameError("BAD_INPUT")
    return a * hbar_eVps


def compute_moments(frame, cutoff, *, graph=None, chunk_size=32):
    """All fixed and learned-basis sufficient moments, bounded in memory.

    Reference shells stay aligned across a changing box. Prepared geometry is
    reusable only with an exactly matching accepted state geometry.
    """
    if not isinstance(frame, FrameInput) or frame.force_contract_confirmed is not True:
        raise FrameError("FORCE_CONTRACT_MISSING")
    frame = validate_frame(**vars(frame))
    if not isinstance(chunk_size, (int, np.integer)) or chunk_size < 1:
        raise FrameError("BAD_INPUT")
    pos = frame.positions[:, frame.magnetic_mask]
    if graph is None:
        graph = build_graph(pos[0], frame.cell[0], frame.pbc, cutoff)
    elif (
        not np.array_equal(pos[0], graph["positions"])
        or not np.array_equal(frame.cell[0], graph["cell"])
        or not np.array_equal(frame.pbc, graph["pbc"])
    ):
        raise FrameError("ID_MISMATCH")
    elif cutoff_radius(pos[0], frame.cell[0], frame.pbc, cutoff) != graph["rc"]:
        raise FrameError("ID_MISMATCH: prepared graph cutoff differs")
    fixed = np.all(pos == pos[0]) and np.all(frame.cell == frame.cell[0])
    chunks = []
    if fixed:
        for start in range(0, len(pos), chunk_size):
            chunks.append(
                moment_arrays(
                    frame.spins[start : start + chunk_size],
                    frame.spin_forces[start : start + chunk_size],
                    graph,
                )
            )
    else:
        for k in range(len(pos)):
            current = build_graph(
                pos[k], frame.cell[k], frame.pbc, cutoff, shell_cuts=graph["shell_cuts"]
            )
            chunks.append(
                moment_arrays(
                    frame.spins[k : k + 1], frame.spin_forces[k : k + 1], current
                )
            )
    result = {key: np.concatenate([row[key] for row in chunks]) for key in chunks[0]}
    result.update(
        channels=CHANNELS,
        frame_index=frame.frame_index.copy(),
        shell_cuts=graph["shell_cuts"].copy(),
    )
    return result


__all__ = [
    "CHANNELS",
    "ORDERS",
    "FrameError",
    "FrameInput",
    "build_graph",
    "compute_moments",
    "from_lammps_fm",
    "load_frame_npz",
    "validate_frame",
]
