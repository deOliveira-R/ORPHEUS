r"""The reference's wall value ``WallAmplitudes(specular, diffuse)``.

The reference specification states a wall by the coordinates :math:`(s, d)` of
its linear part in the triangle K of sub-Markov walls (the positive cone cut by
:math:`s + d \le 1`; its vertices are the vacuum, the mirror and the white wall): the fraction :math:`s` of the incoming
current returned specularly and the fraction :math:`d` returned diffusely, as
plain data (the interface between references and production is plain data).
A wall that returns more than it receives is not physical, so both are in
:math:`[0, 1]` with :math:`s + d \le 1`; a value outside that is a parse
refusal (:class:`ValueError`, :class:`TypeError` for a non-real).

The vocabulary is closed: a wall is vacuum :math:`(0, 0)`, specular
:math:`(s, 0)` or diffuse :math:`(0, d)`. A combined wall (both coordinates
positive) is a scope boundary that issue #650 lifts by widening the
vocabulary, so it raises :class:`NotImplementedError`, not :class:`ValueError`:
no better input lifts it. The bounds are checked first, so an unphysical
combined pair such as :math:`(0.6, 0.6)` is a :class:`ValueError`.

The vertices of K are named once, here: ``WallAmplitudes.VACUUM`` :math:`(0, 0)`,
``WallAmplitudes.MIRROR`` :math:`(1, 0)` and ``WallAmplitudes.WHITE`` :math:`(0, 1)`.

**The exterior is an exact cover, admitted at one door.** A reference states the
walls behind a geometry's integer tags as a mapping ``exterior`` from tag to
:class:`WallAmplitudes`, and that mapping must cover ``geometry.tags`` exactly: a tag
without amplitudes and an entry without a tag are each refused, by tag number.
:func:`admit_exterior` is that predicate, written once; the specification
(``GeometrySpecification``) and every reference consumer of the pair (the characteristic reference's ``Walls``, the
singular-eigenfunction spectrum, the :math:`F_N` moment space, the Galerkin basis space and the shared reference body) call it, so a state the specification refuses cannot
be built by a reference, and the frozen copy it returns cannot be changed by the
caller's later mutation of the mapping it passed.

**The production counterpart.** The production wall is ``WallCondition`` of
:mod:`orpheus.transport.boundary_conditions`, the same point of K spelled as a response network, and
``WallCondition.from_amplitudes(specular, diffuse)`` builds it from these two coordinates; this module is plain data so
that a reference never imports the transport layer. The architecture of record is :ref:`bc-architecture-of-record`.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from typing import Any, ClassVar

from orpheus.geometry.structured_geometry import StructuredGeometry
from orpheus.values.content import ContentIdentity, FrozenMapping
from orpheus.values.scalars import parse_integer, parse_real


@dataclass(frozen=True, eq=False)
class WallAmplitudes(ContentIdentity):
    """The linear part of a wall as coordinates in the triangle K: ``specular`` and ``diffuse`` returned fractions."""

    specular: float
    diffuse: float

    VACUUM: ClassVar["WallAmplitudes"]
    """The vertex :math:`(0, 0)`: nothing returned."""
    MIRROR: ClassVar["WallAmplitudes"]
    """The vertex :math:`(1, 0)`: everything returned specularly."""
    WHITE: ClassVar["WallAmplitudes"]
    """The vertex :math:`(0, 1)`: everything returned diffusely."""

    def __post_init__(self) -> None:
        object.__setattr__(self, "specular", parse_real(self.specular, "WallAmplitudes: the specular amplitude"))
        object.__setattr__(self, "diffuse", parse_real(self.diffuse, "WallAmplitudes: the diffuse amplitude"))
        # The simplex: s >= 0, d >= 0, s + d <= 1 (the per-coordinate upper bounds follow from these).
        if not (self.specular >= 0.0 and self.diffuse >= 0.0 and self.specular + self.diffuse <= 1.0):
            raise ValueError(
                f"WallAmplitudes(specular={self.specular}, diffuse={self.diffuse}): a returned amplitude outside "
                "[0, 1], or returning more than it receives, is not a physical wall"
            )
        if self.specular > 0.0 and self.diffuse > 0.0:
            # SCOPE-BOUNDARY[guard] -- machinery: the reference vocabulary widened to the two-chain wall
            # s*Specular + d*Diffuse (plan 2.9, from_amplitudes). ruling: D3, DRAFT 3.1, approved by the
            # user 2026-10-23. revisit: #650 (state prescribed inflow, combined walls and zero flux).
            raise NotImplementedError(
                f"WallAmplitudes(specular={self.specular}, diffuse={self.diffuse}): a wall is specular or diffuse "
                "in the closed vocabulary; the combined wall waits on #650"
            )


WallAmplitudes.VACUUM = WallAmplitudes(0.0, 0.0)
WallAmplitudes.MIRROR = WallAmplitudes(1.0, 0.0)
WallAmplitudes.WHITE = WallAmplitudes(0.0, 1.0)


def admit_exterior(exterior: Any, geometry: StructuredGeometry, *, owner: str) -> FrozenMapping[int, WallAmplitudes]:
    """The exterior as a frozen copy: an exact cover of ``geometry.tags`` by :class:`WallAmplitudes`, each refusal naming its tag.

    Raises :class:`TypeError` for a non-mapping or a value that is not a ``WallAmplitudes``, and
    :class:`ValueError` for a geometry tag the mapping does not state and for an entry no face carries.
    ``owner`` names the class or function that asked, as the prefix of every refusal.
    """
    where = f"{owner}: the exterior"
    if not isinstance(exterior, Mapping):
        raise TypeError(f"{where} is a mapping from tag to WallAmplitudes, got a {type(exterior).__name__}")
    held: dict[int, WallAmplitudes] = {}
    for tag, wall in exterior.items():
        tag = parse_integer(tag, owner, "an exterior tag")
        if not isinstance(wall, WallAmplitudes):
            raise TypeError(f"{where} states tag {tag} as a {type(wall).__name__}; a wall is a WallAmplitudes (specular, diffuse)")
        held[tag] = wall
    tags = geometry.tags
    unstated = sorted(tags - held.keys())
    unused = sorted(held.keys() - tags)
    if unstated:
        raise ValueError(f"{where} states no amplitudes for the tag(s) {unstated} that the geometry's faces carry (tags {sorted(tags)})")
    if unused:
        raise ValueError(f"{where} states the tag(s) {unused} that no face of the geometry carries (tags {sorted(tags)})")
    return FrozenMapping(held)
