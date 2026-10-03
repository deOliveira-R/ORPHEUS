r"""The observables: what is read off an answer, with no mesh and no physics in it (#405 P2 step 3).

An observable names one number to read from the answer of a question, the
same number whichever answer is asked: a production solution on a mesh, a
reference solution in its own representation, or a published table. Every
answer will read it its own way, ``answer.read(observable)`` (the user's
ruling "G1" of 2026-10-02: the answers grow, the observables do not, so the
answer is the receiver; the first ``read`` lands at #405 P2 step 6), and
return a reading typed by whose claim it is
(:data:`orpheus.reference.ReferenceReading` or
:data:`orpheus.numerics.outcome.ProductionReading`).

The closed set (the user's ruling of 2026-10-03):

* :class:`FluxIntegral` ``(weight)``: the linear functional
  :math:`\langle w, \phi\rangle = \sum_g \int w_g(\vec r)\,\phi_g(\vec r)\,dV`
  of the scalar flux, the weight a mesh-free function over position and
  group. It is the primitive: a reaction rate :math:`\langle w, T\phi\rangle`
  is the flux integral of the weight :math:`T^\top w` (for a removal cell, the
  weight times its cross section; for an emission cell, the weight contracted
  with the emission spectrum). The spelling ``Rate(cells, weight)`` will be a
  constructor the SPECIFICATION resolves into a flux integral from its
  materials, as it resolves a question's keys (P2 step 8, with the removal
  cells of #526); it is not a member here, so one functional has one
  canonical spelling. A region-averaged group flux is the
  flux integral of a weight that is the indicator of the region and the group
  divided by the region's volume;
* :class:`Ratio` ``(numerator, denominator)``: the quotient of two LINEAR
  observables (a flux integral or a point value), such as a spectral index or
  a normalised shape. Only a quotient of two linear functionals of the flux is
  independent of the scale an eigen answer's representative was picked at,
  so the operands are :data:`Linear` and a ratio does not nest. A reference
  reads it by dividing its two enclosures, the bound carried outward
  (:meth:`orpheus.numerics.enclosure.Enclosure.__truediv__`);
* :class:`Eigenvalue` ``()``: the eigenvalue of an eigen question's answer, a
  datum of that answer and not a functional of its flux, read in the chart of
  the question's parameter (as :class:`~orpheus.numerics.question.Nearest`'s
  ``tau`` is), so k and 1/k cannot be confused;
* :class:`PointValue` ``(position, group)``: the scalar flux of one group at
  one position on a one-dimensional geometry's axis. Its weight is a Dirac
  delta, which is no function, so it is its own member, and an answer with no
  point evaluation (a production answer of cell averages) refuses it.

Whether an observable FITS a problem (a weight's group count, region count
and the coordinates a ``Symbolic`` weight may depend on; a point's group and
position; an eigenvalue asked of an eigen question) is decided where the
problem is known, once, by the specification (P2 step 6), never here and never
again inside each answer's ``read``.

**Closed.** Each member is ``@final``: a subclass would pass every
``isinstance`` door while a ``match`` over the sum sends it to its parent's
arm.

**Content identity.** Every observable is a
:class:`~orpheus.numerics.content.ContentIdentity` admitted eagerly, so two
spellings of one observable are one value and an observable can key a cache.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TypeAlias, final, get_args

from orpheus.numerics.content import ContentIdentity, content_digest
from orpheus.numerics.mesh_free_function import MeshFreeFunction, parse_mesh_free_function
from orpheus.numerics.scalars import parse_finite_real, parse_index, parse_member

__all__ = ["Eigenvalue", "FluxIntegral", "Linear", "Observable", "PointValue", "Ratio"]


@final
@dataclass(frozen=True, eq=False)
class FluxIntegral(ContentIdentity):
    r"""The flux paired with a weight, :math:`\langle w, \phi\rangle`."""

    weight: MeshFreeFunction

    def __post_init__(self) -> None:
        parse_mesh_free_function(self.weight, "FluxIntegral", "the weight")
        content_digest(self)


@final
@dataclass(frozen=True, eq=False)
class Ratio(ContentIdentity):
    """The quotient of two linear observables."""

    numerator: Linear
    denominator: Linear

    def __post_init__(self) -> None:
        for noun, operand in (("the numerator", self.numerator), ("the denominator", self.denominator)):
            parse_member(operand, get_args(Linear), "Ratio", noun, "a linear observable")
        content_digest(self)


@final
@dataclass(frozen=True, eq=False)
class Eigenvalue(ContentIdentity):
    """The eigenvalue of an eigen question's answer."""


@final
@dataclass(frozen=True, eq=False)
class PointValue(ContentIdentity):
    """The scalar flux of one group at one position on a one-dimensional axis."""

    position: float
    group: int

    def __post_init__(self) -> None:
        object.__setattr__(self, "position", parse_finite_real(self.position, "PointValue: the position"))
        object.__setattr__(self, "group", parse_index(self.group, "PointValue", "the group"))
        content_digest(self)


Linear: TypeAlias = FluxIntegral | PointValue
"""The observables linear in the flux: the operands of a :class:`Ratio`."""

Observable: TypeAlias = FluxIntegral | Ratio | Eigenvalue | PointValue
"""What is read off an answer: a closed set (#405 P2, the user's ruling of 2026-10-03)."""
