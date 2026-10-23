r"""Shared builders of the S2b gates on the typed realized boundary operators (wave 2, slice S2b).

Not a test module (no ``test_`` prefix). Spec: ``.claude/plans/wave2_boundary_carve.md`` sections 2.6, 3 (laws L-A, L-E,
L-L) and 6.2 (S2b gates 4, 9, 13 to 20).

Why this module imports no production module at the top
=======================================================

At the base tip of slice S2b ``import orpheus.transport`` raises ``ImportError`` (``method.py`` still imports the
retired ``orpheus.geometry.BC``) and ``orpheus.sn.boundary.realized_operators`` does not exist. A collection-time
import would turn every gate into the same collection error and hide which gate fails for which reason (``vv-principles``
anti-pattern 17c). So every production module is read at first use, by the accessors below, and each gate's first red
names its own line: the missing module for the typed operators (gates 14 to 18, 20), the missing table for the realizer
rows (gates 4, 9, 13).

What it builds
==============

* the trace ladder of a quadrature (:func:`trace_of`), the one place the layout incantation is written;
* the certified mirror permutation of a face (:func:`mirror_pairing`), through the production
  ``ordinate_pairing`` (the certification is L4's and is not under test here);
* the dense matrix of a realized face operator by unit vectors (:func:`dense`), and the reversal matrices of a
  slab face (:func:`reversal`), which the reciprocity law L-A (i) and the naturality law L-E read;
* :func:`asymmetric_trace`, the L-A (iv) fixture: a delegating quadrature proxy whose OUTFLOW-side weights are scaled,
  so the trace has ``M_+ != M_-`` (``[M]`` ``scratch/layer_dag/wave2_spec/_d3_trace_axes.py``: GL-8 scaled by 1.1 builds
  with ``M_+ = 0.5563404348779897`` and ``M_- = 0.5057640317072634``; the unperturbed trace has ``M_+ == M_-``).
"""

from __future__ import annotations

import importlib
import types
from typing import Any

import numpy as np

#: The four quadratures of the B3.4c measurement, plus the S16 slab of laws L-A (i) and L-E, lazily built.
QUADRATURES: dict[str, Any] = {
    "gl8": lambda: _quad().gauss_legendre(n_ordinates=8),
    "gl16": lambda: _quad().gauss_legendre(n_ordinates=16),
    "product(2,4)": lambda: _quad().product(n_mu=2, n_phi=4),
    "ls6": lambda: _quad().level_symmetric(sn_order=6),
    "lebedev17": lambda: _quad().lebedev(17),
}


def quadrature_class() -> Any:
    """The ``Quadrature`` class (read at first use)."""
    return _quad()


def _quad() -> Any:
    from orpheus.numerics.quadrature import Quadrature

    return Quadrature


def ro() -> types.ModuleType:
    """``orpheus.sn.boundary.realized_operators`` (the file S2b adds)."""
    return importlib.import_module("orpheus.sn.boundary.realized_operators")


def sn_table() -> Any:
    """``SN_REALIZERS``, the sealed SN realizer table."""
    return importlib.import_module("orpheus.sn.boundary.realizer").SN_REALIZERS


def bc() -> types.ModuleType:
    """``orpheus.transport.boundary_conditions`` (S2a)."""
    return importlib.import_module("orpheus.transport.boundary_conditions")


def trace_of(quadrature: Any, faces: tuple[str, ...] = ("xmin", "xmax")) -> Any:
    """The ``AngularTraceSpace`` of ``quadrature`` over ``faces`` (geometry-blind: the normals come from the names)."""
    from orpheus.numerics.face_layout import FaceLayout
    from orpheus.numerics.spaces.angular_trace_space import AngularTraceSpace

    layout = FaceLayout.from_named_shapes([(f, (int(quadrature.n_ordinates),)) for f in faces])
    return AngularTraceSpace.from_quadrature_and_layout(quadrature, layout)


def method_space(quadrature: Any, face: str, trace: Any) -> Any:
    """The SN method space of ``face`` over ``trace`` (what ``realize`` receives as ``space``)."""
    from orpheus.sn.mesh.method_space import SNMethodSpace

    return SNMethodSpace.for_face(quadrature=quadrature, face=face, trace=trace)


def site(face: str, *, curved: bool = False, dimension: int = 1) -> Any:
    """The ``FaceSite`` of a Cartesian face, or of a radial body's outer face when ``curved``."""
    from orpheus.geometry import CoordSystem

    return bc().FaceSite(face, CoordSystem.CYLINDRICAL if curved else CoordSystem.CARTESIAN, dimension)


def mirror_pairing(quadrature: Any, face: str, dimension: int = 1) -> Any:
    """The certified mirror permutation of ``face``'s own plane (``ordinate_pairing`` of ``mirror_motion``)."""
    from orpheus.geometry.deck import mirror_motion
    from orpheus.geometry.faces import AXIS_NAMES, face_normal
    from orpheus.numerics.quadrature.ordinate_pairing import ordinate_pairing

    axis, _ = face_normal(face)
    return ordinate_pairing(quadrature, mirror_motion(AXIS_NAMES[axis], dimension))


def identity_pairing(quadrature: Any) -> Any:
    """The identity permutation: what the translation wrap's image carries (a translation moves no direction)."""
    from orpheus.geometry.transformation import Permutation

    return Permutation(np.arange(int(quadrature.n_ordinates)))


def scaled_ordinate_permutation(alpha: float, permutation: Any, face: str, partner: str, trace: Any) -> Any:
    """``ScaledOrdinatePermutation(alpha, permutation, face, partner, trace)``: the ONE place the constructor is spelled by
    keyword, so a ratified rename of ``permutation`` (rider (b): a certified ``OrdinatePairing``) is a one-line edit here."""
    return ro().ScaledOrdinatePermutation(alpha=alpha, permutation=permutation, face=face, partner=partner, trace=trace)


def n_outflow(trace: Any, face: str) -> int:
    return int(len(trace.outflow_indices_for_face(face)))


def dense(op: Any, n_outflow_rows: int) -> np.ndarray:
    """The dense ``(n_inflow, n_outflow)`` matrix of a realized face operator, by unit vectors of ``Gamma_+``."""
    return np.stack([np.asarray(op.apply(np.eye(n_outflow_rows)[j])) for j in range(n_outflow_rows)], axis=1)


def positive_outflow(n: int, seed: int = 1234) -> np.ndarray:
    """A positive, non-uniform vector on ``Gamma_+`` (the synthetic outflow of the ratio gate)."""
    return 0.5 + np.random.default_rng(seed).random(n)


def trace_weights(trace: Any, face: str) -> tuple[np.ndarray, np.ndarray]:
    """``(w |mu| on Gamma_+, w |mu| on Gamma_-)``: the ladder's metric on the two half-spaces (read, never re-derived)."""
    return (
        np.asarray(trace.outflow_space(face).inner_product_weights, dtype=float),
        np.asarray(trace.inflow_space(face).inner_product_weights, dtype=float),
    )


def reversal(quadrature: Any, trace: Any, face: str) -> np.ndarray:
    r"""The slab's time reversal :math:`\mu \to -\mu` from ``Gamma_-(face)`` to ``Gamma_+(face)``, as an
    ``(n_out, n_in)`` permutation matrix: entry ``[k, j]`` is 1 iff outflow ordinate ``k`` has the direction cosine
    opposite to inflow ordinate ``j``. Built by matching the cosines (``mean_axis_cosine``), not from the realizer.
    """
    mu = np.asarray(quadrature.mean_axis_cosine(0), dtype=float)
    out_ix = np.asarray(trace.outflow_indices_for_face(face))
    in_ix = np.asarray(trace.inflow_indices_for_face(face))
    r = np.zeros((out_ix.size, in_ix.size))
    for j, g in enumerate(in_ix):
        k = int(np.argmin(np.abs(mu[out_ix] + mu[g])))
        if abs(mu[out_ix[k]] + mu[g]) > 1e-14:
            raise AssertionError(f"no reversed partner for ordinate {int(g)} on {face}")
        r[k, j] = 1.0
    return r


def reciprocity_defect(b: np.ndarray, quadrature: Any, trace: Any, face: str) -> float:
    r""":math:`\max |M_+^{-1} B^{\mathsf T} M_- - P_{\rm rev} B P_{\rm rev}|` (law L-A (i)): zero for a reciprocal
    response, with ``M_\pm`` the half-space metrics and ``P_rev`` the slab reversal."""
    m_out, m_in = trace_weights(trace, face)
    p = reversal(quadrature, trace, face)
    left = b.T * (m_in[None, :] / m_out[:, None])  # the metric ratio first: 1.0 exactly where the paired weights are bit-equal
    right = p @ b @ p
    return float(np.max(np.abs(left - right)))


def asymmetric_trace(scale: float = 1.1) -> tuple[Any, Any]:
    """``(proxy_quadrature, trace)`` with the OUTFLOW side of ``xmax`` scaled by ``scale`` (``M_+ != M_-``), law L-A (iv)."""
    from orpheus.numerics.spaces.angular_trace_space import build_omega_dot_n

    base = QUADRATURES["gl8"]()
    omega_dot_n = build_omega_dot_n(base, ("xmax",))[0]
    weights = np.array(base.weights, dtype=float)
    weights[omega_dot_n > 0] *= scale

    class Proxy:
        def __init__(self) -> None:
            self.weights = weights

        def __getattr__(self, name: str) -> Any:
            return getattr(base, name)

    proxy = Proxy()
    return proxy, trace_of(proxy)
