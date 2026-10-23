r"""The refusals of the conditions package (``errors.py`` of the conditions package).

**The root.** :class:`BoundaryError` is a ``ValueError`` carrying ``.declaration``: the object, or its type, that was
refused. It replaces the string ``.law`` of the retiring ``orpheus.geometry.boundary`` errors, so a test (and a caller)
compares the object (``err.declaration == MirrorDeck()``) and not a spelling. ``BoundaryError`` is transport's: geometry
has ``DeckError`` and numerics has ``OrdinatePairingError``, because geometry holds no module named for boundary
conditions (the eleventh exchange of the boundary-law ontology).

**The parse refusals** (raised at construction of a condition, so a wall outside the cone is unspellable):

* :class:`SubmarkovViolationError`: the return fraction exceeds one (``WallCondition.specular(1.2)``);
* :class:`BoundaryResponseNotPositiveError`: a response weight or an inflow weight is negative (``white(-0.1)``).

The two names are kept from the retiring vocabulary (retirement-audit B.11: two message vocabularies survive).

**The door refusals:**

* :class:`BoundaryRefused`: ``RealizerTable.realize`` / ``require`` met an :class:`~.admission.Refused` admission; it
  carries the ``method``, the ``face`` and the ``admission`` (so the axis and the issue are data, not prose);
* :class:`UnboundTagError`: a face carries a ``tag`` for which ``bind`` was given no condition;
* :class:`UnusedConditionError`: a condition is declared for a ``tag`` no face carries (the exact cover catches a typo
  that a reachable-subset rule would hide);
* :class:`NonInvariantConditionError`: a condition is not invariant under a non-identity element of the problem's
  colouring group that moves a face of its tag; it carries the ``isometry`` (or ``None`` where the group is not
  computed) and the ``faces``.

The three certification errors that check the SN realizer's OWN index sets stay SN-side (``orpheus/sn/boundary/errors.py``,
slice S2b), rooted here. The three that certify the quadrature (ERR-042, 044, 045) become numerics' errors.

Layer: ``transport`` (L2). Imports nothing from ``orpheus`` (the types in the annotations are under ``TYPE_CHECKING``).
"""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from orpheus.geometry.colouring import FaceIsometry
    from orpheus.transport.boundary_conditions.admission import Refused

__all__ = [
    "BoundaryError",
    "BoundaryRefused",
    "BoundaryResponseNotPositiveError",
    "NonInvariantConditionError",
    "SubmarkovViolationError",
    "UnboundTagError",
    "UnusedConditionError",
]


class BoundaryError(ValueError):
    """The root of every refusal about a boundary declaration; ``.declaration`` is the refused object or its type."""

    def __init__(self, message: str, declaration: object = None) -> None:
        super().__init__(message)
        self.declaration = declaration


class SubmarkovViolationError(BoundaryError):
    """The return fraction of a wall exceeds one: not a physical wall."""


class BoundaryResponseNotPositiveError(BoundaryError):
    """A response or inflow weight is negative: outside the positive cone."""


class BoundaryRefused(BoundaryError):
    """A method's realizer table refused a declaration; the ``admission`` says on which axis and under which issue."""

    def __init__(self, message: str, declaration: object = None, *, method: str, face: str, admission: Refused) -> None:
        super().__init__(message, declaration)
        self.method = method
        self.face = face
        self.admission = admission


class UnboundTagError(BoundaryError):
    """A face carries a tag with no declared condition."""

    def __init__(self, message: str, declaration: object = None, *, tag: int, face: str) -> None:
        super().__init__(message, declaration)
        self.tag = tag
        self.face = face


class UnusedConditionError(BoundaryError):
    """A condition is declared for a tag no face carries."""

    def __init__(self, message: str, declaration: object = None, *, tag: int) -> None:
        super().__init__(message, declaration)
        self.tag = tag


class NonInvariantConditionError(BoundaryError):
    """A condition is not invariant under the colouring group of its carrier (``isometry`` is None where the group is not computed)."""

    def __init__(self, message: str, declaration: object = None, *, tag: int, isometry: FaceIsometry | None, faces: tuple[str, ...]) -> None:
        super().__init__(message, declaration)
        self.tag = tag
        self.isometry = isometry
        self.faces = faces
