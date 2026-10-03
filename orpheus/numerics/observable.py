r"""The observables: what is read off an answer, with no mesh and no physics in it (#405 P2 step 3).

An observable names one number to read from the answer of a question, the
same number whichever answer is asked: a production solution on a mesh, a
reference solution in its own representation, or a published table. Every
answer reads it its own way, ``answer.read(observable)`` (the user's ruling
"G1" of 2026-10-02: the answers grow, the observables do not, so the answer
is the receiver), and returns a reading typed by whose claim it is
(:data:`orpheus.reference.ReferenceReading` or
:data:`orpheus.numerics.outcome.ProductionReading`).

The closed set (the user's ruling of 2026-10-03):

* :class:`FluxIntegral` ``(weight)``: the linear functional
  :math:`\langle w, \phi\rangle = \sum_g \int w_g(\vec r)\,\phi_g(\vec r)\,dV`
  of the scalar flux, the weight a mesh-free function over position and
  group. It is the primitive: a reaction rate :math:`\langle w, T\phi\rangle`
  is the flux integral of the weight :math:`T^\top w` (for a removal cell, the
  weight times its cross section; for an emission cell, the weight contracted
  with the emission spectrum). The spelling ``Rate(cells, weight)`` is a
  constructor the SPECIFICATION resolves into a flux integral from its
  materials, as it resolves a question's keys; it is not a member here, so one
  functional has one canonical spelling. A region-averaged group flux is the
  flux integral of a weight that is the indicator of the region and the group
  divided by the region's volume;
* :class:`Ratio` ``(numerator, denominator)``: the quotient of two
  observables (a spectral index, a normalised shape). A reference reads it by
  dividing its two enclosures, the bound carried outward
  (:meth:`orpheus.numerics.enclosure.Enclosure.__truediv__`);
* :class:`Eigenvalue` ``()``: the eigenvalue of an eigen question's answer, a
  datum of that answer and not a functional of its flux;
* :class:`PointValue` ``(position, group)``: the scalar flux of one group at
  one position on a one-dimensional geometry's axis. Its weight is a Dirac
  delta, which is no function, so it is its own member, and an answer with no
  point evaluation (a production answer of cell averages) refuses it.

Whether a weight FITS a problem (its group count, its region count, the
coordinates a ``Symbolic`` weight may depend on) is decided where the problem
is known: by the specification's admission of a mesh-free datum, reused when
an answer reads, never re-implemented here.

**Content identity.** Every observable is a
:class:`~orpheus.numerics.content.ContentIdentity` admitted eagerly, so two
spellings of one observable are one value and an observable can key a cache.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TypeAlias, get_args

from orpheus.numerics.content import ContentIdentity, content_digest
from orpheus.numerics.mesh_free_function import MeshFreeFunction, parse_mesh_free_function
from orpheus.numerics.scalars import parse_finite_real, parse_integer

__all__ = ["Eigenvalue", "FluxIntegral", "Observable", "PointValue", "Ratio"]


@dataclass(frozen=True, eq=False)
class FluxIntegral(ContentIdentity):
    r"""The flux paired with a weight, :math:`\langle w, \phi\rangle`."""

    weight: MeshFreeFunction

    def __post_init__(self) -> None:
        parse_mesh_free_function(self.weight, "weight")
        content_digest(self)


@dataclass(frozen=True, eq=False)
class Ratio(ContentIdentity):
    """The quotient of two observables."""

    numerator: Observable
    denominator: Observable

    def __post_init__(self) -> None:
        for role, operand in (("numerator", self.numerator), ("denominator", self.denominator)):
            if not isinstance(operand, Observable):
                kinds = ", ".join(kind.__name__ for kind in get_args(Observable))
                raise TypeError(f"Ratio: the {role} is an observable ({kinds}), got a {type(operand).__name__}")
        content_digest(self)


@dataclass(frozen=True, eq=False)
class Eigenvalue(ContentIdentity):
    """The eigenvalue of an eigen question's answer."""


@dataclass(frozen=True, eq=False)
class PointValue(ContentIdentity):
    """The scalar flux of one group at one position on a one-dimensional axis."""

    position: float
    group: int

    def __post_init__(self) -> None:
        object.__setattr__(self, "position", parse_finite_real(self.position, "PointValue: the position"))
        group = parse_integer(self.group, "PointValue", "the group")
        if group < 0:
            raise ValueError(f"PointValue: the group is a non-negative index, got {group}")
        object.__setattr__(self, "group", group)
        content_digest(self)


Observable: TypeAlias = FluxIntegral | Ratio | Eigenvalue | PointValue
"""What is read off an answer: a closed set (#405 P2, the user's ruling of 2026-10-03)."""
