r"""The codim-1 face measure: the spatial-trace partial-current weight.

L1 primitive (mathematics, knows no neutrons). A layout module holds layouts
only (:mod:`orpheus.numerics.coproduct_layout`); the measure :math:`|\Omega\cdot\hat n|\,w`
lives here, beside no layout.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray


__all__ = ["face_streaming_normal"]


def face_streaming_normal(
    normal_coefficient: NDArray, quadrature_weight: NDArray,
) -> NDArray:
    r"""The spatial-trace partial-current measure :math:`|\Omega\cdot\hat n_f|\,w`.

    The magnitude of the *normal component of the streaming flux through a
    codim-1 spatial boundary face*, weighted by the angular quadrature:
    with :math:`c = \Omega\cdot\hat n_f` it is the partial-current metric
    :math:`G_s = |\Omega\cdot\hat n_f|\,w_n` (Lewis & Miller §3.7), the
    boundary inner product under which the
    :class:`~orpheus.numerics.spaces.angular_trace_space.AngularTraceSpace`
    boundary operators are self-adjoint. It vanishes exactly at
    **grazing** incidence (:math:`\Omega \perp \hat n_f`, so
    :math:`|\Omega\cdot\hat n_f| = 0`) — a tangential ordinate carries no
    partial current across the face.

    L1 primitive (mathematics, knows no neutrons). Its single production
    consumer is the trace-metric build in
    :class:`~orpheus.numerics.spaces.angular_trace_space.AngularTraceSpace`.

    .. note::

       **This is the SPATIAL-trace measure only — it is NOT the ψ½
       starting-direction (angular-pole) block's metric.**  An earlier
       design (``facefield_codim1_design.md`` §3.2) proposed unifying both
       codim-1 metrics through this one kernel, reading the pole's angular
       through-flux coefficient :math:`(1-\mu^2)\big|_{\mu=\pm1} = 0` as
       the ψ½ block's Hilbert inner product.  That unification is
       **refuted** (ERR-067): :math:`(1-\mu^2)|_{\rm pole}` is an
       *operator coefficient* (the angular-redistribution strength that
       makes :math:`\mu=\pm1` a straight characteristic), **not** a state
       metric.  ψ½ is a first-class radial state field whose Hilbert
       metric is the radial cell volume :math:`G_{\rm sd} = V_{\rm cell}`
       (SPD, nonzero) — set by its operator role, not by an integration
       weight.  The two codim-1 metrics do **not** share a kernel; the
       metric is a per-leaf property (spatial trace → this partial
       current; pole → :math:`V_{\rm cell}`), exactly as the bulk's metric
       is its own per-leaf :math:`V\,w`.  See
       :mod:`orpheus.numerics.spaces.radial_characteristic_space`.

    Parameters
    ----------
    normal_coefficient : NDArray
        The ordinate's signed normal streaming coefficient
        :math:`c = \Omega\cdot\hat n_f` on a spatial boundary face. Any
        shape.
    quadrature_weight : NDArray
        The angular quadrature weight :math:`w`, broadcast-compatible with
        ``normal_coefficient``.

    Returns
    -------
    NDArray
        :math:`|c|\,w`, the elementwise normal streaming measure.

    References
    ----------
    * Lewis, E.E. & Miller, W.F. (1993). *Computational Methods of Neutron
      Transport*. §3.7 (the partial-current boundary metric).
    """
    return np.abs(normal_coefficient) * quadrature_weight
