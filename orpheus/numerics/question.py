r"""The question values: what is asked of a system, with no physics in it (#405 P1 step 7).

A question is posed on a family of operators :math:`E(p)` over a parameter
space, at a base point :math:`p_0` (the posing sequence, layer 2,
``.claude/plans/posing_sequence.md``). Three questions are minted:

* :class:`Eigen` asks where, on the line through ``point`` along one
  ``parameter`` direction, the family is singular, and which of those poles
  (``mode``) is wanted. The k-eigenvalue, the classical c-eigenvalue, a boron
  search and a critical extent are each one direction;
* :class:`FixedSource` asks for :math:`E(p_0)^{-1} q`, the response to a
  source :math:`q` that does not depend on the unknown;
* :class:`Response` asks for :math:`E(p_0)^{-\dagger} R`, the importance of a
  detector :math:`R`: the adjoint fixed-source question.

**Physics-free.** A question names a direction of the system only as an
opaque, digestable KEY. What the key means (a set of reaction-grid cells
scaled together, a geometric extent, a nuclide density) is declared by
whoever resolves it: the reference specification first, the system later.
So nothing here names a reaction, and the values live in ``numerics``.

**The point** is a frozen mapping from a parameter key to its offset from
the physical value; the empty mapping is the physical point and the default.
A zero offset is kept, not dropped: ``{key: 0.0}`` names a key the resolver
must still see, so it is a different question (a cache miss, never a wrong
hit).

**No forward/adjoint flag.** The role is the type: a source is held by
:class:`FixedSource` and a detector by :class:`Response`, so a
detector-less adjoint or a forward question given a detector cannot be
spelled. The eigen adjoint :math:`\psi^\dagger` belongs to the eigen ANSWER,
so :class:`Eigen` has no adjoint either. The source and the detector are
mesh-free functions (:mod:`orpheus.numerics.mesh_free_function`); neither
carries its role, and the lift each takes into phase space (the section for
a source rate, the retraction's adjoint for a detector) is picked by the
question that holds it.

**Content identity.** Every value is a
:class:`~orpheus.numerics.content.ContentIdentity`, admitted EAGERLY: an
undigestable key or a non-finite offset is refused by the constructor, never
first by ``hash``.

Not minted here (the posing sequence's later units, #529): the pencil and the
spectral map a resolved parameter derives, the mode law that finds a pole,
the ``Enclosed(region)`` mode, the time direction (α) and the evolution
question.
"""

from __future__ import annotations

from collections.abc import Hashable, Mapping
from dataclasses import dataclass, field
from typing import Any, TypeAlias, get_args

from orpheus.numerics.content import ContentIdentity, ContentlessError, FrozenMapping, content_digest
from orpheus.numerics.mesh_free_function import MeshFreeFunction, parse_mesh_free_function
from orpheus.numerics.scalars import parse_finite_real, parse_member



def _admit_point(point: Mapping[Any, float]) -> FrozenMapping[Any, float]:
    """The point as a frozen mapping of finite real offsets, refusals naming the key."""
    if not isinstance(point, Mapping):
        raise TypeError(f"the point is a mapping from parameter key to offset, got a {type(point).__name__}")
    return FrozenMapping((key, parse_finite_real(offset, f"the offset of {key!r}")) for key, offset in point.items())


def _admit_key(key: Any, where: str) -> None:
    """A parameter key is hashable, so it cannot change after the question is keyed.

    A point key is hashed by the mapping that holds it; the parameter is
    checked here. An unhashable key (a list, an array, a view over a
    caller's ``dict``) is mutable, and its content would not be fixed.
    """
    try:
        hash(key)
    except TypeError:
        raise ContentlessError(
            f"{where}: an unhashable {type(key).__name__} is not a key (it is mutable, so its content is not fixed)"
        ) from None


@dataclass(frozen=True, eq=False)
class Fundamental(ContentIdentity):
    """The first pole met from the removal-dominated end of the direction."""


@dataclass(frozen=True, eq=False)
class Nearest(ContentIdentity):
    """The pole nearest ``tau``, read in the parameter's chart."""

    tau: float

    def __post_init__(self) -> None:
        object.__setattr__(self, "tau", parse_finite_real(self.tau, "Nearest: tau"))


Mode: TypeAlias = Fundamental | Nearest


@dataclass(frozen=True, eq=False)
class Eigen(ContentIdentity):
    """Where, on the line through ``point`` along ``parameter``, is the family singular?

    ``gauge`` declares the functional an eigen flux is scaled by (an
    eigenvector has no scale of its own): an opaque key, as the parameter is,
    which the specification resolves. ``None`` asks for the specification's
    declared default, which its canonical question writes in.
    """

    parameter: Hashable
    point: Mapping[Any, float] = field(default_factory=FrozenMapping)
    mode: Mode = field(default_factory=Fundamental)
    gauge: Hashable | None = None

    def __post_init__(self) -> None:
        _admit_key(self.parameter, "Eigen.parameter")
        if self.gauge is not None:
            _admit_key(self.gauge, "Eigen.gauge")
        object.__setattr__(self, "point", _admit_point(self.point))
        parse_member(self.mode, get_args(Mode), "Eigen", "the mode", "a mode")
        content_digest(self)


@dataclass(frozen=True, eq=False)
class FixedSource(ContentIdentity):
    r""":math:`E(p_0)^{-1} q`, for a source ``q`` that does not depend on the unknown."""

    source: MeshFreeFunction
    point: Mapping[Any, float] = field(default_factory=FrozenMapping)

    def __post_init__(self) -> None:
        parse_mesh_free_function(self.source, "FixedSource", "the source")
        object.__setattr__(self, "point", _admit_point(self.point))
        content_digest(self)


@dataclass(frozen=True, eq=False)
class Response(ContentIdentity):
    r"""The importance of a detector ``R``, :math:`E(p_0)^{-\dagger} R`."""

    detector: MeshFreeFunction
    point: Mapping[Any, float] = field(default_factory=FrozenMapping)

    def __post_init__(self) -> None:
        parse_mesh_free_function(self.detector, "Response", "the detector")
        object.__setattr__(self, "point", _admit_point(self.point))
        content_digest(self)


Question: TypeAlias = Eigen | FixedSource | Response

__all__ = ["Eigen", "FixedSource", "Fundamental", "Mode", "Nearest", "Question", "Response"]
