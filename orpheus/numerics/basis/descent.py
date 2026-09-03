r"""``Descent`` — the two realizations of the functions on an orbit space, as ONE object.

The functions on :math:`M/H` have two honest realizations (#429 Part V.5 §C,
user-ruled 2026-08-31): **upstairs**, the :math:`H`-invariant subspace of a
basis on the base — for the real spherical harmonics and :math:`O(2)_x`,
the :math:`m = 0` column, :math:`\{Y_\ell^0\}`; for a coordinate mirror, the
σ-even slots — and **downstairs**, the quotient's OWN classical basis when
one exists — :math:`\{P_\ell(\mu)\}` on :math:`S^2/O(2)_a`. They are
isomorphic (:math:`\mathrm{Funcs}(M)^H \cong \mathrm{Funcs}(M/H)`), and
without a witness they are a Cardinal-Rule-2 twin. This class is the witness
and the DISCRIMINATOR, written where both realizations can see it:

    **downstairs when the quotient has a classical named basis
    (**:math:`S^2/O(2)_a \to \{P_\ell\}`**), upstairs otherwise (**:math:`S^2/\sigma_a`
    **has no classical family — its σ-even harmonics are the only spellable
    realization).**

That sentence is what :attr:`Descent.frame_basis` returns and what
:meth:`Quadrature._harmonic_basis <orpheus.numerics.quadrature.directional.Quadrature._harmonic_basis>`
binds — the basis a frame on THIS orbit space carries is DERIVED from the
entry, never chosen at the call site.

Which slots of a basis on the base descend is the ENTRY's question
(:meth:`Quotient.descending_slots <orpheus.numerics.manifold.Quotient.descending_slots>`,
user-ruled 2026-09-02): a function descends iff it is constant on the fibres
of the entry's own :attr:`~orpheus.numerics.manifold.Quotient.quotient_map`,
sampled at generic base points and their images under the group's generic
elements. This class reads it.

The isomorphism is checkable, and at the BIT tier rather than "to machine
precision": `[M]` 2026-09-02, ``downstairs.evaluate(π(Ω)) == parent.evaluate(Ω)[…, slots]``
is ``array_equal`` (:math:`\max|\Delta| = 0.0`) on 7 of 7 shipped sphere rules
at :math:`L = 4` and on GL8/GL16 at :math:`L \le 2` — because
:class:`~orpheus.numerics.basis.legendre_basis.LegendreBasis` spells the
polynomial exactly as the harmonic table spells its :math:`m = 0` column.

⚠ The upstairs face is slot-ALIGNED only about the harmonics' own polar
axis (x): `[M]` about y and z the invariant subspace is one-dimensional in
every degree but spreads over several slots from :math:`\ell \ge 2` — and is
still aligned at :math:`\ell \le 1`, so a refusal keyed on measured
alignment would be silently inert exactly where every solve mints
(:math:`L = 0`). The refusal is keyed on the AXIS.
"""
from __future__ import annotations

from dataclasses import dataclass
from functools import cached_property

import numpy as np
from numpy.typing import NDArray

from orpheus.numerics.basis.base import Basis
from orpheus.numerics.basis.legendre_basis import LegendreBasis
from orpheus.numerics.basis.spherical_harmonic_basis import (
    MirrorEvenSphericalHarmonicBasis,
    SphericalHarmonicBasis,
)
from orpheus.numerics.manifold import Quotient, SPHERE
from orpheus.numerics.symmetry import AXIS_LETTER

__all__ = ["Descent"]

#: The polar axis of the real spherical harmonics (``cos θ = μ_x`` in
#: ``_evaluate_real_sh``) — the ONE axis about which the :math:`SO(2)`
#: upstairs face is a set of SLOTS rather than a subspace.
_HARMONIC_POLAR_AXIS = "x"


@dataclass(frozen=True)
class Descent:
    r"""The descent of a basis on the base to an orbit-space ENTRY.

    Parameters
    ----------
    entry : Quotient
        The catalogue entry :math:`M/H` (carrying its quotient map and its
        pushforward reference).
    parent : Basis
        A basis on ``entry.base`` — the full real harmonics for every
        shipped entry.
    """

    entry: Quotient
    parent: Basis

    def __post_init__(self) -> None:
        if self.parent.domain != self.entry.base:
            raise ValueError(
                f"Descent: the parent basis eats {self.parent.domain.name} "
                f"but the entry {self.entry.name} is a quotient of "
                f"{self.entry.base.name}; a descent pulls functions on the "
                f"orbit space back to its BASE."
            )

    @classmethod
    def for_entry(cls, entry: Quotient, L: int) -> "Descent":
        r"""The descent of the degree-:math:`L` real harmonics to a sphere entry."""
        if entry.base != SPHERE:
            raise NotImplementedError(
                f"Descent.for_entry: only quotients of S^2 have a shipped "
                f"parent basis (the real spherical harmonics); "
                f"{entry.name} is a quotient of {entry.base.name}."
            )
        return cls(entry=entry, parent=SphericalHarmonicBasis(L=L))

    # ── the entry's probe, read ───────────────────────────────────────

    @cached_property
    def descending_slots(self) -> NDArray:
        r"""Boolean mask over the parent's mode axes: the slots that descend to the entry."""
        return self.entry.descending_slots(self.parent)

    # ── the two realizations ──────────────────────────────────────────

    @cached_property
    def downstairs(self) -> Basis | None:
        r"""The quotient's OWN classical basis, or ``None`` when it has none.

        :math:`S^2/O(2)_a` → :class:`LegendreBasis` about ``a``; a mirror
        quotient has no classical family; the trivial quotient's downstairs
        IS the parent.
        """
        by = self.entry.by
        axis = by.rotation_axis
        if axis is not None:
            L = _truncation_order(self.parent)
            return LegendreBasis(L=L, axis=AXIS_LETTER[axis])
        if self.entry.base == self.entry.realization and by.is_trivial:
            return self.parent
        return None

    @cached_property
    def upstairs(self) -> Basis:
        r"""The parent's invariant SUB-basis, kept in the parent's layout.

        A coordinate mirror: :class:`MirrorEvenSphericalHarmonicBasis`
        (structurally zeroed σ-odd columns). :math:`O(2)_x`: the
        :math:`m = 0` column is the whole face, and the tree ships no masked
        harmonic basis for it — the downstairs realization is the one every
        consumer binds — so the face is answered by :attr:`upstairs_columns`
        (the descended table), not by a class; asking for it as a BASIS is
        refused naming the work.
        """
        by = self.entry.by
        mirror = by.mirror_axis
        if mirror is not None:
            return MirrorEvenSphericalHarmonicBasis(
                L=_truncation_order(self.parent), mirror_axis=mirror,
            )
        if by.is_trivial:
            return self.parent
        raise NotImplementedError(
            f"Descent.upstairs: no masked harmonic basis ships for "
            f"{self.entry.name}; its upstairs face is the descended table "
            f"(`upstairs_columns`), and every consumer binds the downstairs "
            f"realization (`frame_basis`)."
        )

    @property
    def discriminator(self) -> str:
        r"""The sentence that says which realization a frame binds — and why."""
        return (
            "downstairs when the quotient has a classical named basis "
            "(S^2/O(2)_a -> {P_l}), upstairs otherwise (S^2/sigma_a has no "
            "classical family; its sigma-even harmonics are the only "
            "spellable realization)"
        )

    @property
    def frame_basis(self) -> Basis:
        r"""The basis a frame on THIS orbit space binds — the discriminator, applied."""
        down = self.downstairs
        if down is not None:
            return down
        return self.upstairs

    # ── the isomorphism, checkable at the bit ─────────────────────────

    def upstairs_columns(self, points: NDArray) -> NDArray:
        r"""The parent's table at ``points`` restricted to the descending slots — the upstairs face, tabulated.

        Shape ``(N, K)`` with the slots in row-major mask order. For
        :math:`O(2)_x` that is one slot per degree in ascending :math:`\ell`
        — the :math:`m = 0` column. Refused for an axial group about any
        other axis: there the invariant subspace is not slot-aligned from
        :math:`\ell \ge 2`, and a restriction by slots would be WRONG rather
        than incomplete (axis-keyed, so the refusal fires at every
        :math:`L`, including the :math:`L \le 1` where alignment happens to
        hold).
        """
        axis = self.entry.by.rotation_axis
        if axis is not None and AXIS_LETTER[axis] != _HARMONIC_POLAR_AXIS:
            raise NotImplementedError(
                f"Descent.upstairs_columns: the real harmonics' polar axis is "
                f"{_HARMONIC_POLAR_AXIS!r}; about {AXIS_LETTER[axis]!r} the "
                f"{self.entry.by.name}-invariant subspace is not a set of slots (it spreads "
                f"over several from l >= 2). The downstairs face is "
                f"available at every axis; an upstairs face there needs a "
                f"rotated harmonic table, which no consumer has asked for."
            )
        table = self.parent.evaluate(points)
        # Boolean-mask a table of shape (N, *modes): the mask ranges over the
        # mode axes whatever their rank (`[M]` every shipped entry's mask is
        # 2-D — the harmonics' (L+1, 2L+1); a 1-D arm had no witness).
        return table[:, self.descending_slots & self.live_slots]

    @cached_property
    def live_slots(self) -> NDArray:
        r"""The parent's REAL slots — a padded layout's ``|m| > ℓ`` slots are identically zero and
        descend vacuously, so they are excluded from every count and every restriction here."""
        live = getattr(self.parent, "live_slot_mask", None)
        if live is None:
            return np.ones(self.descending_slots.shape, dtype=bool)
        return np.asarray(live, dtype=bool)

    @cached_property
    def descending_real_slots(self) -> NDArray:
        r"""``descending_slots`` restricted to the parent's live slots — the honest denominator."""
        return self.descending_slots & self.live_slots

    def max_abs_difference(self, points: NDArray) -> float:
        r"""``max |downstairs(π(points)) − upstairs(points)|`` — ``0.0`` at the bit on every shipped rule."""
        down = self.downstairs
        if down is None:
            raise NotImplementedError(
                f"Descent.max_abs_difference: {self.entry.name} has no "
                f"downstairs realization to compare against."
            )
        descended = down.evaluate(self.entry.quotient_map(points))
        upstairs = self.upstairs_columns(points)
        if descended.shape != upstairs.shape:
            raise ValueError(
                f"Descent: the two realizations disagree in SIZE — "
                f"downstairs {descended.shape}, upstairs {upstairs.shape}."
            )
        return float(np.max(np.abs(descended - upstairs)))

    def is_isomorphism(self, points: NDArray) -> bool:
        r"""``True`` iff the two realizations agree bit-for-bit at ``points``."""
        return self.max_abs_difference(points) == 0.0


def _truncation_order(basis: Basis) -> int:
    order = getattr(basis, "L", None)
    if not isinstance(order, int):
        raise TypeError(
            f"Descent: the parent basis {type(basis).__name__} carries no "
            f"truncation order L."
        )
    return order
