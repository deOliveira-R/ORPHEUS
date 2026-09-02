r"""The angular HEAD of a moment space — the surface every moment carrier reads.

A moment field's space is ``<head> ⊗ cells``: the head is the coefficient
space of the basis the quadrature's frame bound (#429 tracker 2.5, 2026-09-02).
Two families ship. The real spherical harmonics' head is the rectangular
``(L+1, 2L+1)`` table with the addition-theorem-shifted ``[l + m]`` column and
zero padding outside :math:`|m| \le \ell`; the Legendre basis on
:math:`S^2/SO(2)_a` has a FLAT head, ``(L+1,)`` — one coefficient per degree,
because the trivial isotypic component of :math:`SO(2)` is one-dimensional in
every degree (`[M]` 2026-09-02, a rank test about every axis).

Until #429's fused commit every consumer that indexed ``values[0, 0]`` or
sliced ``values[l, :2l+1]`` read the FIRST family's layout as if it were the
contract. It was the family's: on a flat head ``values[0, 0]`` is group 0's
spatial slice and raises nothing (`[M]` the verification memo's H-15 —
``scalar_flux``, ``isotropic_part``, ``anisotropic_part``, ``l_block``, the
fission ℓ=0 dyad, and the material field's per-degree contraction, which
spelled the m-axis into its einsum). So the LAYOUT is the head's to say —
which index tuple is the isotropic slot, which selects the degree-:math:`\ell`
block, how many leading axes the head owns (``len(shape)``), and how it
truncates within its own family — and a consumer READS it rather than
assuming it.

Structural (``Protocol``) and ``runtime_checkable``: the two space classes
satisfy it by carrying the members, and a consumer holding ``space.factors[0]``
narrows with ``isinstance`` — the same idiom as
:class:`~orpheus.numerics.basis.base.TruncatedBasis` on the basis side.
"""
from __future__ import annotations

from typing import TYPE_CHECKING, Protocol, runtime_checkable

if TYPE_CHECKING:
    from orpheus.numerics.space import FunctionSpace

__all__ = ["MomentHead"]


@runtime_checkable
class MomentHead(Protocol):
    r"""The head factor of a moment space: a truncated family that knows its own layout."""

    @property
    def L(self) -> int:
        """The truncation order — the degree the family is cut at."""
        ...

    @property
    def shape(self) -> tuple[int, ...]:
        """The head's own axes (``(L+1, 2L+1)`` rectangular, ``(L+1,)`` flat)."""
        ...

    @property
    def name(self) -> str:
        ...

    @property
    def isotropic_slot(self) -> tuple[int, ...]:
        r"""The index tuple of the :math:`\ell = 0` coefficient within the head's axes."""
        ...

    def degree_block(self, l: int, /) -> tuple[int | slice, ...]:
        r"""The index tuple selecting the degree-:math:`\ell` block within the head's axes."""
        ...

    def truncated(self, L_new: int, /) -> "FunctionSpace":
        """This family's space at the lower order ``L_new``, under this head's own name."""
        ...
