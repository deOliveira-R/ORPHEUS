r"""The affine part :math:`q` of a wall: what flows in from the exterior (``inflow.py`` of the conditions package).

**What this module is.** The exterior of a face is modelled by the affine map :math:`\gamma_-\psi|_f = R\,\gamma_+\psi|_f
+ q` with :math:`q \ge 0` a function on the inflow half-trace :math:`\Gamma_-`. :mod:`.network` holds :math:`R`; this
module holds :math:`q`, as a weighted sum of profiles: :class:`InflowTerms`.

**The frame is the type, the content is open.** The three tree implementors of the old ``evaluate(space)`` (an MMS
inflow reading ``space.directions``, a deliberately lying profile and a row-order probe) are Python callables over a
space, not closed expressions, so no closed ``Symbolic`` profile union exists. What is decidable, and what the
binding's invariance refusal needs, is the FRAME the profile is written in:

* :class:`FaceFrameInflow`: written in the face's own frame (:math:`\mu = \Omega\cdot n_f` and tangential position);
  :attr:`~InflowTerms.face_local` is True, so every isometry carries it from face to face;
* :class:`LabFrameInflow`: written in laboratory directions; :attr:`~InflowTerms.face_local` is False, so
  ``BoundaryCondition.invariant_under`` holds only for the identity (conservative and deterministic:
  ``function o h == function`` is undecidable for an open callable, and sampling the two faces at a tolerance is a
  check keyed on nothing numerical).

Both wrappers are ``ContentIdentity``-free frozen dataclasses keyed by their sampler: a sampler that is a plain object
has identity equality and no content, so an :class:`InflowTerms` holding one has identity equality too and its digest
raises ``ContentlessError`` (the uniform-only case is ordinary content, and a ``GeometrySpecification`` cannot hold a
sampled term because it states plain data).

**A uniform inflow is a number.** Every uniform contribution folds into :attr:`InflowTerms.uniform` (``c * inflow(v)``
contributes ``c * v``), so ``inflow(1.0) + inflow(2.0) == inflow(3.0)`` without a profile object for a constant.

**Sampling is the method's, not the declaration's.** :class:`InflowSampler` is the open protocol (the old
``InflowSourceSpec.evaluate(space)``, renamed ``sample``: a realization welded onto a declaration is gone). The single
verb that evaluates :math:`q` on a method space is ``RealizerTable.sample_inflow`` in :mod:`.realization`, whose default
body is ``sample_terms``.

**Gotchas.**

* The drop rule is the network's: a zero weight is not stored, so ``0.0 * WallCondition.framed(f)`` holds the empty
  inflow and equals ``WallCondition.vacuum()``. ``-0.0`` is ``0.0``.
* Sampled terms keep first-appearance order (nothing sorts by digest: a digest of an open sampler raises), and two terms
  of EQUAL profile merge their weights.
* ``NoSource`` retires: an inflow of zero is :meth:`InflowTerms.empty`.
* Weights may be negative here: a displacement (:mod:`.displacement`) holds signed inflow. The positive cone is the
  condition's (``q >= 0``).

Theory: :ref:`bc-architecture-of-record`.

Layer: ``transport`` (L2). Imports ``numerics`` (the sampler's space, under ``TYPE_CHECKING``) and ``values`` only.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Protocol, TypeAlias, get_args

import numpy as np

from orpheus.values.content import ContentIdentity
from orpheus.values.scalars import parse_finite_real

if TYPE_CHECKING:
    from orpheus.numerics.space import FunctionSpace

__all__ = ["FaceFrameInflow", "InflowProfile", "InflowSampler", "InflowTerms", "LabFrameInflow"]


class InflowSampler(Protocol):
    """The open protocol of a profile: its values on a method's boundary-source space."""

    def sample(self, space: FunctionSpace) -> np.ndarray:
        r""":math:`q` on ``space`` with exactly ``space.shape``, in the units of :math:`\psi` (NOT W-normalised)."""
        ...


@dataclass(frozen=True)
class FaceFrameInflow:
    """A profile written in the face's own frame; face-local."""

    sampler: InflowSampler


@dataclass(frozen=True)
class LabFrameInflow:
    """A profile written in laboratory directions; NOT face-local (the invariance check of ``bind`` owes it)."""

    sampler: InflowSampler


InflowProfile: TypeAlias = FaceFrameInflow | LabFrameInflow


@dataclass(frozen=True, eq=False)
class InflowTerms(ContentIdentity):
    r""":math:`q` as a weighted sum of profiles: ``uniform`` (all uniform contributions folded) plus ``sampled``.

    ``sampled`` is a tuple of ``(weight, profile)`` in first-appearance order; the parse drops a zero weight, merges
    equal profiles and refuses a profile that is not an :data:`InflowProfile`.
    """

    uniform: float
    sampled: tuple[tuple[float, InflowProfile], ...]

    def __post_init__(self) -> None:
        object.__setattr__(self, "uniform", parse_finite_real(self.uniform, "InflowTerms.uniform"))
        merged: list[list] = []
        for weight, profile in self.sampled:
            if not isinstance(profile, get_args(InflowProfile)):
                raise TypeError(f"InflowTerms: a sampled term holds one of {[p.__name__ for p in get_args(InflowProfile)]}, got {profile!r}")
            weight = parse_finite_real(weight, "InflowTerms sampled weight")
            for entry in merged:
                if entry[1] == profile:
                    entry[0] += weight
                    break
            else:
                merged.append([weight, profile])
        object.__setattr__(self, "sampled", tuple((weight, profile) for weight, profile in merged if weight != 0.0))

    @classmethod
    def empty(cls) -> InflowTerms:
        """No inflow: ``uniform == 0.0`` and nothing sampled."""
        return cls(0.0, ())

    def __add__(self, other: object) -> InflowTerms:
        """Term-by-term sum: the uniform parts add, equal profiles merge, new profiles append."""
        if not isinstance(other, InflowTerms):
            return NotImplemented
        return InflowTerms(self.uniform + other.uniform, self.sampled + other.sampled)

    def __mul__(self, scalar: object) -> InflowTerms:
        """Every weight multiplied by a real scalar (a ``bool`` is not one)."""
        if isinstance(scalar, bool) or not isinstance(scalar, int | float):
            return NotImplemented
        return InflowTerms(scalar * self.uniform, tuple((scalar * weight, profile) for weight, profile in self.sampled))

    @property
    def is_zero(self) -> bool:
        """Nothing flows in."""
        return self.uniform == 0.0 and not self.sampled

    @property
    def face_local(self) -> bool:
        """No term is written in laboratory directions."""
        return not any(isinstance(profile, LabFrameInflow) for _, profile in self.sampled)
