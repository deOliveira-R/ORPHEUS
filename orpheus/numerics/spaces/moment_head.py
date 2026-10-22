r"""The angular HEAD of a moment space — the one door its consumers read it through.

A moment field's space is ``<head> ⊗ cells``: the head is the coefficient
space of the basis the quadrature's frame bound (#429 tracker 2.5, 2026-09-02),
one axis whose layout is a degree layout. Two families ship. The real
spherical harmonics' head is the rectangular ``(L+1, 2L+1)`` table with the
addition-theorem-shifted ``[l + m]`` column and zero padding outside
:math:`|m| \le \ell`; the Legendre basis on :math:`S^2/O(2)_a` has a FLAT
head, ``(L+1,)`` — one coefficient per degree, because the trivial isotypic
component of :math:`O(2)_a` (equivalently of its rotation half, which has the
same orbits) is one-dimensional in every degree (`[M]` 2026-09-02, a rank test
about every axis).

Until #429's fused commit every consumer that indexed ``values[0, 0]`` or
sliced ``values[l, :2l+1]`` read the FIRST family's layout as if it were the
contract. It was the family's: on a flat head ``values[0, 0]`` is group 0's
spatial slice and raises nothing (`[M]` the verification memo's H-15 —
``scalar_flux``, ``isotropic_part``, ``anisotropic_part``, ``l_block``, the
fission ℓ=0 dyad, and the material field's per-degree contraction, which
spelled the m-axis into its einsum). So the LAYOUT is the head's to say —
which index tuple is the isotropic slot, which selects the degree-:math:`\ell`
block, how many leading axes the head owns — and a consumer READS it rather
than assuming it.

Since S7 of the numerics untangling (rulings 72 and 81) the layout is a value
the basis carries (:mod:`orpheus.numerics.coefficient_layout`) and the head
AXIS carries it (:attr:`~orpheus.numerics.axis.Axis.layout`), so the dressed
head of a frame and the continuum head of its basis read one layout. Until S7
it lived on two space classes (``SphericalHarmonicSpace``, ``LegendreSpace``)
behind a runtime-checkable ``MomentHead`` protocol, which four production
``isinstance`` doors checked; they read :func:`moment_head_axis` now.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from orpheus.numerics.axis import Axis
    from orpheus.numerics.coefficient_layout import DegreeLayout
    from orpheus.numerics.space import FunctionSpace

__all__ = ["head_layout", "moment_head_axis", "truncated_head"]


def moment_head_axis(head: FunctionSpace | None, *, owner: str | None = None) -> Axis:
    r"""The head axis of a moment head space: its one axis, carrying a degree layout.

    The one door the moment carriers read their layout through (the four
    production sites that narrowed on ``MomentHead`` until S7). Refuses a
    space that is not a moment head — no axes, several axes, or an axis whose
    layout has no degree structure (a sampled axis, a cell head) — naming it
    and, when given, the ``owner`` that asked.
    """
    from orpheus.numerics.coefficient_layout import DegreeLayout

    axes = None if head is None else head.axes
    if axes is not None and len(axes) == 1 and isinstance(axes[0].layout, DegreeLayout):
        return axes[0]
    asker = f"{owner}: " if owner else ""
    what = "None" if head is None else f"{type(head).__name__} {head.name!r}"
    raise TypeError(
        f"{asker}{what} is not a moment head "
        f"(one axis carrying a degree layout: a truncation order, an "
        f"isotropic slot and degree blocks)."
    )


def head_layout(head: FunctionSpace | None, *, owner: str | None = None) -> DegreeLayout:
    r"""The degree layout of a moment head: :func:`moment_head_axis`'s layout,
    the one test (no second dispatch on the family)."""
    return moment_head_axis(head, owner=owner).layout  # type: ignore[return-value]


def truncated_head(head: FunctionSpace, L_new: int, /) -> FunctionSpace:
    r"""This family's head at the lower order ``L_new`` — re-minted by the head
    axis's ORIGIN (CS4c step 6 item 6.2c-ii, ruling O-3: *re-mint AND re-axis,
    never slice*).

    An axis-built head has one axis, and that axis remembers what minted it:

    * a :class:`~orpheus.numerics.basis.base.TruncatedBasis` (the head is
      the basis's own continuum coefficient space,
      :meth:`FunctionSpace.of_basis <orpheus.numerics.space.FunctionSpace.of_basis>`)
      re-spans the family at ``L_new`` — same functions, same domain — and
      hands back ITS space;
    * a :class:`~orpheus.numerics.frame.GalerkinFrame` (the head is the
      frame's Parseval-dressed :func:`~orpheus.numerics.frame_faces.basis_space`,
      the ONE moment space the tree binds) re-poses itself at ``L_new`` and
      hands back its dressed space at that order — the metric is the
      frame's to install at EVERY order, because the discrete Gram's
      verdict can FLIP with :math:`L` (`[M]` ``folded_product(2,4)`` is
      DENSE at :math:`L = 2` and DIAGONAL at :math:`L = 1`), so a slice of
      the parent's dressing is undefined there.

    A truncated moment field therefore lands on the SAME space its carrier
    mints at ``L_new`` (structurally equal — ruling O-5), never on a
    second-metric twin.
    """
    from orpheus.numerics.axis import FromBasis, FromFrame
    from orpheus.numerics.basis.base import TruncatedBasis
    from orpheus.numerics import frame_faces
    from orpheus.numerics.frame import GalerkinFrame
    from orpheus.numerics.space import FunctionSpace

    axis = moment_head_axis(head)
    L = head_layout(head).L  # the same door, read for the order
    if not 0 <= L_new <= L:
        raise ValueError(f"truncated_head: L_new={L_new} must lie in [0, {L}].")
    match axis.origin:
        case FromFrame(GalerkinFrame() as frame):
            return frame_faces.basis_space(frame.at_order(L_new))
        case FromBasis(TruncatedBasis() as basis):
            return FunctionSpace.of_basis(basis.at_order(L_new))
        case origin:
            raise ValueError(
                f"truncated_head: the head axis of {head.name!r} was minted "
                f"through neither a truncated basis nor a Galerkin frame "
                f"(origin={type(origin).__name__}), so nothing can re-mint "
                f"the family at L_new={L_new}."
            )
