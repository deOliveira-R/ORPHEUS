r"""Function space for Legendre moment coefficients on an :math:`O(2)_a` orbit space.

The coefficient space of :class:`~orpheus.numerics.basis.legendre_basis.LegendreBasis`
— :math:`\{P_\ell(\mu)\}_{\ell \le L}` on :math:`S^2/O(2)_a` — a **FLAT** head of
shape ``(L+1,)``. The sibling of
:class:`~orpheus.numerics.spaces.spherical_harmonic_space.SphericalHarmonicSpace`
in the :class:`~orpheus.numerics.spaces.moment_head.MomentHead` family, and the
second member of it: the first family has a rectangular head, this one has
one coefficient per degree, and the moment carriers read which through the
protocol rather than assuming either (#429 tracker 3.4, 2026-09-02).

The continuum metric it carries is :math:`4\pi/(2\ell+1)` — the Gram of
:math:`P_\ell` against the **pushforward** of :math:`d\Omega` along the
quotient map, :math:`\pi_* d\Omega = 2\pi\,d\mu` (Archimedes' hat-box;
:attr:`~orpheus.numerics.manifold.Quotient.reference`):
:math:`\int_{-1}^{1} P_\ell^2 \, 2\pi\,d\mu = 4\pi/(2\ell+1)`. It coincides
exactly with the spherical-harmonic
:attr:`~orpheus.numerics.basis.spherical_harmonic_basis.SphericalHarmonicBasis.metric_per_ell`,
which is what makes the descent :math:`\{Y_\ell^0\} \cong \{P_\ell\}` an
ISOMETRY and not merely an isomorphism. ⚠ Not the bare ``LEGENDRE`` mass-2
normalisation :math:`2/(2\ell+1)`: that is a factor :math:`2\pi` away and
would move every operator end's metric (`[M]` the verification memo §3).
"""
from __future__ import annotations

from dataclasses import dataclass, replace

from orpheus.numerics.basis.legendre_basis import LegendreBasis
from orpheus.numerics.space import FunctionSpace

__all__ = ["LegendreSpace"]


@dataclass(frozen=True)
class LegendreSpace(FunctionSpace):
    r"""Function space of Legendre moment coefficients up to degree :math:`L` on :math:`S^2/O(2)_a`.

    Parameters
    ----------
    name : str
        Inherited. Convention: ``"legendre_space(S^2/O2_<axis>)"`` — the
        orbit space's own name, READ off the basis's domain by
        :meth:`from_L` rather than spelled a second time, because two axes
        are two spaces (the tree carries two poles; tracker 2.4).
    shape : tuple[int, ...]
        Inherited. MUST equal ``(L + 1,)``; ``__post_init__`` checks.
    inner_product_weights : NDArray, optional
        The per-degree continuum Gram :math:`4\pi/(2\ell+1)` (module
        docstring); a frame's dressed ``basis_space`` REPLACES it with the
        Parseval inverse exactly as for the spherical-harmonic space.
    L : int, default 0
        Maximum degree retained.
    spent_axis : str, default ``"x"``
        The axis of the spent stabiliser :math:`O(2)_a` (``axis`` is a
        :class:`FunctionSpace` method, hence the longer name).

    Notes
    -----
    Equality and hashing are by ``(name, shape)`` alone, inherited from
    :class:`FunctionSpace` — ``L`` is encoded in ``shape`` and the axis in
    ``name``.
    """

    L: int = 0
    spent_axis: str = "x"

    def __post_init__(self) -> None:
        super().__post_init__()
        if self.spent_axis not in ("x", "y", "z"):
            raise ValueError(
                f"LegendreSpace: spent_axis must be x/y/z, got {self.spent_axis!r}."
            )
        expected = (self.L + 1,)
        if self.shape != expected:
            raise ValueError(
                f"LegendreSpace: shape={self.shape} inconsistent with "
                f"L={self.L}; expected shape={expected} (one coefficient "
                f"per degree — the flat head)."
            )

    def __eq__(self, other: object) -> bool:
        return FunctionSpace.__eq__(self, other)

    def __hash__(self) -> int:
        return FunctionSpace.__hash__(self)

    @classmethod
    def from_L(cls, L: int, axis: str = "x") -> "LegendreSpace":
        r"""The canonical Legendre space for degree :math:`L` about ``axis``.

        The metric is sourced from :class:`LegendreBasis` so the
        :math:`4\pi/(2\ell+1)` formula lives in exactly one place.
        """
        basis = LegendreBasis(L=L, axis=axis)
        return cls(
            name=f"legendre_space({basis.domain.name})",
            shape=(L + 1,),
            inner_product_weights=basis.metric_per_ell,
            L=L,
            spent_axis=axis,
        )

    # ── the MomentHead surface ───────────────────────────────────────

    @property
    def isotropic_slot(self) -> tuple[int, ...]:
        r"""``(0,)`` — the :math:`\ell = 0` coefficient of a flat head."""
        return (0,)

    def degree_block(self, l: int, /) -> tuple[int | slice, ...]:
        r"""``(l,)`` — the degree-:math:`\ell` block of a flat head is one coefficient."""
        if not 0 <= l <= self.L:
            raise ValueError(
                f"LegendreSpace.degree_block: l={l} out of range [0, {self.L}]."
            )
        return (l,)

    def truncated(self, L_new: int, /) -> "LegendreSpace":
        r"""This family's space at the lower order ``L_new``, about the same axis, under this head's own name."""
        if not 0 <= L_new <= self.L:
            raise ValueError(
                f"LegendreSpace.truncated: L_new={L_new} must lie in [0, {self.L}]."
            )
        return replace(type(self).from_L(L_new, self.spent_axis), name=self.name)

    # ── delegated convention (single source in the basis) ──────────

    @property
    def basis(self) -> LegendreBasis:
        """The associated :class:`LegendreBasis`, determined by ``(L, axis)``."""
        return LegendreBasis(L=self.L, axis=self.spent_axis)
