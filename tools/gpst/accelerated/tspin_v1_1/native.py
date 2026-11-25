"""Optional fused edge arithmetic; a missing library uses NumPy."""

import ctypes
import os
from pathlib import Path

import numpy as np

_LIB = None
if os.environ.get("TSPIN_NUMPY_ONLY") != "1":
    try:
        _LIB = ctypes.CDLL(str(Path(__file__).with_name("_edge_kernel.so")))
    except OSError:
        pass
if _LIB is not None:
    _LIB.tspin_edges.argtypes = [ctypes.c_int64] * 3 + [ctypes.c_void_p] * 16
    _LIB.tspin_edges.restype = ctypes.c_int


def edge_arrays(s, ft, graph):
    """Return un-reduced arrays, or None for the NumPy implementation."""
    if _LIB is None:
        return None
    if s.shape[-1] != 3:
        raise ValueError("BAD_INPUT")
    t, ns, _ = s.shape
    edges = len(graph["i"])
    inputs = [
        np.ascontiguousarray(s, dtype=np.float64),
        np.ascontiguousarray(ft, dtype=np.float64),
        np.ascontiguousarray(graph["i"], dtype=np.int64),
        np.ascontiguousarray(graph["j"], dtype=np.int64),
        np.ascontiguousarray(graph["w"], dtype=np.float64),
        np.ascontiguousarray(graph["direction"], dtype=np.float64),
    ]
    if (
        inputs[2].shape != (edges,) or inputs[3].shape != (edges,)
        or inputs[4].shape != (edges,) or inputs[5].shape != (edges, 3)
        or any(np.any(a < 0) or np.any(a >= ns) for a in inputs[2:4])
    ):
        raise ValueError("ID_MISMATCH")
    # Keep the platform's NumPy dot-three summation order for the zonal term.
    z = np.empty((t, 2, edges))
    zf = np.empty_like(z)
    for k, idx in enumerate(inputs[2:4]):
        z[:, k] = np.einsum("tej,ej->te", s[:, idx], graph["direction"])
        zf[:, k] = np.einsum("tej,ej->te", ft[:, idx], graph["direction"])
    inputs += [z, zf]
    outputs = [
        np.empty((t, 4, edges)), np.empty((t, 4, edges)),
        np.empty((2, edges, t)).transpose(2, 0, 1),
        np.empty((2, edges, t)).transpose(2, 0, 1),
        np.empty((t, edges)), np.empty((t, edges)),
        np.empty((t, edges, 3, 3)), np.empty((t, edges, 3, 3)),
    ]
    code = _LIB.tspin_edges(t, ns, edges, *(a.ctypes.data for a in inputs + outputs))
    if code:
        raise ValueError("BAD_INPUT")
    return outputs
