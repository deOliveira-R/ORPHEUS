r"""Prescribed-inflow source for the affine boundary law.

The §16A.1 affine form

.. math::

    \gamma_- \psi \;=\; R\,G\,\gamma_+ \psi \;+\; q

has three terms; this module ships the :math:`q` term -- a
vector-valued quantity on the incoming trace :math:`\Gamma_-`. For
the homogeneous case (vacuum / albedo / reflective / white /
periodic without external injection), :class:`NoSource` is the
sentinel. For a constant scalar inflow, use
:class:`ConstantInflowSource`.

The :class:`InflowSourceSpec` Protocol is :class:`runtime_checkable`
so consumers can ``isinstance(s, InflowSourceSpec)`` to validate the
source contract without subclassing. Custom sources (e.g. a
spatially-varying beam injection) only need to implement
:meth:`evaluate` with the documented signature.

References
----------

* Grand Report v3 §16A.1 (affine boundary form).
* ``.claude/plans/transient-giggling-cake.md`` -- Wave 3 brief.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Protocol, runtime_checkable

import numpy as np

if TYPE_CHECKING:  # pragma: no cover - typing only
    # Typed against the SHARED base, not SN's AngularFaceTraceSpace: diffusion's
    # Γ₋ is one number per face (a ScalarTraceSpace) and will fill this slot
    # when #290 P5 wires its fixed-source arm. A source that reads
    # ``space.directions`` therefore requires an ANGULAR trace, and on a P1
    # trace it fails at the attribute — acceptable while diffusion refuses
    # PrescribedInflow outright; #290 P5 owns the typed refusal.
    from orpheus.numerics.space import FunctionSpace


__all__ = [
    "InflowSourceSpec",
    "NoSource",
    "ConstantInflowSource",
]


@runtime_checkable
class InflowSourceSpec(Protocol):
    r"""Prescribed inflow source :math:`q \in \Gamma_-(f)`.

    :math:`q` is an **element of a specific space** — the inflow half-trace of
    one face — and there are only two ways to specify an element of a discrete
    space: *extensionally* (the numbers, in row order) or *intensionally* (a
    function on the continuum, plus a rule for discretizing it). This Protocol
    is the intensional one: **the implementor is a function, and the space it
    is handed is the discretization.**

    ⛔ **This took a ``shape`` tuple until campaign P6** (see
    ``.claude/plans/affine_boundary_source_channel.md``), with the rationale
    "the source is a property of the boundary, not of any trace space, so it
    takes the shape it must fill rather than a trace object". That is a
    *half-extensional* interface — it demands the numbers while refusing to say
    which row is which — and it is defensible for exactly one implementor,
    :class:`ConstantInflowSource`, because a constant is row-independent.
    Anything angular had to smuggle the row order in through its constructor:
    `[M]` the §4.6 manufactured inflow needed the per-row :math:`\mu` **in trace
    order**, the face's coordinate, and a rank sniff to tell the probe call from
    the delivery call, none of which a shape can carry.

    What the space supplies, and the old signature could not
    ========================================================

    ``space`` is the face's :math:`\Gamma_-(f)` — for SN an
    :class:`~orpheus.numerics.spaces.angular_trace_space.AngularFaceTraceSpace`:

    ==========================  ==================================================
    ``space.directions``        :math:`\Omega_n` per row, **in this space's own
                                row order** — so an angular source is spellable
                                at all, and cannot be silently permuted
    ``space.shape``             the exact array to return, trailing axes and all
                                (`[M]` ``(12, 2, 6)`` on a 2-D face)
    ``space.face``              which face, without parsing a name
    ``space.ordinate_indices``  which global ordinates these rows are
    ==========================  ==================================================

    Contract
    ========

    Return ``np.ndarray`` of ``dtype=float`` with **exactly** ``space.shape``.

    ⭐ **The returned values are the angular flux** :math:`\gamma_-\psi` **on
    those rows, in the same units as** :math:`\psi` — NOT a
    :math:`W`-normalised density. This is the one part of the contract no
    signature can enforce, and it is stated here because getting it wrong is a
    `[M]` **×2652** error on Gauss–Legendre that is numerically
    indistinguishable from a double delivery. If your source is defined through
    a scalar flux (``φ = Σ w ψ``), the :math:`1/W` that converts it belongs to
    **your** definition — and it is the FULL weight sum, not the inflow-
    restricted one (`[M]` 2.0 vs 1.0 on GL-8; using the restricted sum is a
    silent factor of two in the same direction as the double delivery).

    ⭐ **Why ``q \in \Gamma_-`` now holds by construction.** The source is
    handed :math:`\Gamma_-` and returns its shape, so there is nowhere off the
    incoming trace to write. Until P6 that guarantee rested on a presence probe
    (``assert_source_lives_on_incoming_trace``) which a source could **decline**
    — `[M]` returning zeros at the probe shape and ``7.0`` at the delivery shape
    skipped the certification and delivered anyway. The probe was retired with
    this signature; see ERR-047.
    """

    def evaluate(self, space: "FunctionSpace") -> np.ndarray:
        r"""Return :math:`q` on ``space``, with exactly ``space.shape``."""
        ...


@dataclass(frozen=True)
class NoSource:
    r"""Sentinel for the homogeneous case :math:`q = 0`.

    Default for all rank-1 BCs that don't inject external flux
    (vacuum, albedo, reflective, white, periodic). The evaluate
    method returns a fresh zero array of the requested shape;
    callers may either skip the source addition entirely (an
    optimization) or add the zeros (uniform code path).
    """

    def evaluate(self, space: "FunctionSpace") -> np.ndarray:
        return np.zeros(space.shape, dtype=float)


@dataclass(frozen=True)
class ConstantInflowSource:
    r"""Constant scalar inflow :math:`q(\Omega) = \text{value}` on every inflow ordinate.

    Useful for fixed-source problems where a uniform incident flux
    is imposed at one face. The output array has the requested
    shape; downstream consumers are responsible for masking out the
    outflow side (the source LIVES on :math:`\Gamma_-` only -- see
    :class:`~orpheus.geometry.boundary._errors.BoundarySourceNotOnIncomingTraceError`).

    Parameters
    ----------
    value : float
        The scalar inflow level applied to every entry.
    """

    value: float

    def evaluate(self, space: "FunctionSpace") -> np.ndarray:
        # Reads nothing but the shape — and that is the point: a CONSTANT is
        # row-independent, so it is the one source the retired ``shape``
        # signature was ever sufficient for. It is unchanged in behaviour here
        # (``space.shape`` IS the shape it used to be handed) and is the
        # regression control for the migration: if this one moves, the
        # signature change broke the delivery, not the source.
        return np.full(space.shape, float(self.value), dtype=float)
