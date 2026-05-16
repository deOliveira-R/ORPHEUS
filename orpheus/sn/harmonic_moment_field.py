r"""Spherical-harmonic moment field of an angular quantity.

Issue #197 PR-TYPED-4 — the typed wrapper for a real-spherical-harmonic
moment field :math:`\phi_\ell^m(\vec r, g)` of an angular quantity
(typically the angular flux :math:`\psi`).  Sits between
:class:`~orpheus.sn.angular_flux.AngularFlux` and the scattering
operator :math:`\Lambda` as the natural data carrier of the
:math:`R \cdot \Lambda \cdot M \cdot \psi` Galerkin pipeline.

Frozen dataclass (``coding-elegance`` Pattern 4 — illegal states
unrepresentable): the underlying ``(L+1, 2L+1, ng, nx, ny)`` ndarray is
wrapped with its :class:`SNMesh` reference and its truncation order
:math:`L` so shape mismatches surface at construction time, not at the
consumer of a downstream einsum.

Reads as the math via dunder arithmetic (``coding-elegance`` Pattern 1)
and named slicing / reduction primitives (Pattern 3):

* ``phi_a + phi_b`` returns a new :class:`HarmonicMomentField`.
* ``alpha * phi`` returns a new :class:`HarmonicMomentField`.
* ``phi.l_block(l)`` returns the per-:math:`\ell` block view shape
  ``(2ℓ+1, ng, nx, ny)`` — retires the explicit
  ``moments[l, :n_m][..., ix, iy]`` slicing.
* ``phi.isotropic_part()`` / ``.anisotropic_part()`` decompose by
  :math:`\ell = 0` vs :math:`\ell \ge 1`.
* ``phi.scalar_flux()`` extracts the :math:`(\ell=0, m=0)` moment as
  a :class:`~orpheus.sn.scalar_flux.ScalarFlux` (under the
  no-prefactor SH convention used by
  :func:`~orpheus.numerics.spherical_harmonics.evaluate_real_sh`,
  :math:`Y_0^0 = 1`, so the isotropic moment IS the scalar flux
  directly — no normalisation factor).
* ``phi.truncate(L_new)`` drops :math:`\ell > L_{\rm new}` blocks.

Units: same as the source angular quantity.  When produced by
:meth:`HarmonicMomentProjection.apply` from an
:class:`AngularFlux`, the moment field inherits flux units
:math:`[1/(\rm cm^2\,s\,sr\,eV)]`.

Why distinct from :class:`AngularFlux` / :class:`ScalarFlux`
============================================================

The moment field lives in **moment space**
(:math:`(L+1) \cdot (2L+1)` coefficients per cell + group); the
angular flux lives in **per-ordinate space** (:math:`N` directions per
cell + group).  Cross-type addition between the two is undefined (the
spaces have different dimensions); the conversion is the
:math:`M / R` Galerkin pair living in
:mod:`orpheus.numerics.projection`.  Keeping the types distinct
enforces this at the dunder-algebra layer — ``moments + psi`` does
not type-check; the legitimate route is ``moments = M.apply(psi)`` or
``psi = R.apply(moments)``.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from .geometry import SNMesh
    from .scalar_flux import ScalarFlux


__all__ = ["HarmonicMomentField"]


@dataclass(frozen=True)
class HarmonicMomentField:
    r"""Real-spherical-harmonic moment field of an angular quantity.

    Stores coefficients :math:`\phi_\ell^m(g, \vec r)` as an ndarray
    of shape ``(L+1, 2L+1, ng, nx, ny)``, with the
    addition-theorem-shifted :math:`m`-index where slot ``l + m``
    holds the :math:`(\ell, m)` entry; entries outside
    :math:`|m| \le \ell` are zero by convention.

    Parameters
    ----------
    values : np.ndarray
        Moment coefficients of shape ``(L+1, 2L+1, ng, nx, ny)``.
    mesh : SNMesh
        The SN phase-space carrier — supplies the canonical
        ``(ng, nx, ny)`` triple for shape validation.
    L : int
        Maximum harmonic order retained.  Determines the leading two
        axes' sizes: ``values.shape[:2] == (L+1, 2L+1)``.

    Notes
    -----
    Same role here as :class:`~orpheus.sn.angular_flux.AngularFlux`
    plays for the per-ordinate field: a typed wrapper that pins the
    layout invariant (Pattern 4) and exposes domain operations
    (Pattern 1).  The cross-type ``__add__`` with
    :class:`AngularFlux` is **deliberately undefined** — moment-space
    and angular-space algebra happens through the
    :math:`M \cdot \psi \to \phi` projection and
    :math:`R \cdot \phi \to q` reconstruction, not via dunder.
    """

    values: np.ndarray
    mesh: "SNMesh"
    L: int

    def __post_init__(self) -> None:
        expected = (
            self.L + 1, 2 * self.L + 1,
            self.mesh.ng, self.mesh.nx, self.mesh.ny,
        )
        if self.values.shape != expected:
            raise ValueError(
                f"HarmonicMomentField expects shape "
                f"(L+1, 2L+1, ng, nx, ny) = {expected}; "
                f"got {self.values.shape}"
            )

    # ── Slicing / decomposition (Pattern 3 — named intermediates) ──────

    def l_block(self, l: int) -> np.ndarray:
        r"""View of one :math:`\ell`-block, shape ``(2ℓ+1, ng, nx, ny)``.

        Returns the slice ``values[l, :2*l+1]`` — the legitimate
        :math:`m`-entries for that :math:`\ell` (the trailing
        zero-padding outside :math:`|m| \le \ell` is excluded).  Use
        this to retire the explicit ``moments[l, :n_m][..., ix, iy]``
        slicing pattern (``coding-elegance`` Pattern 3).
        """
        if not 0 <= l <= self.L:
            raise ValueError(
                f"HarmonicMomentField.l_block: l={l} out of range [0, {self.L}]"
            )
        return self.values[l, : 2 * l + 1]

    def isotropic_part(self) -> "HarmonicMomentField":
        r"""Return the :math:`\ell = 0` (isotropic) projection.

        Same shape as ``self``; all :math:`\ell \ge 1` blocks zeroed.
        Used by the foldable-vs-residual scattering split when the
        consumer wants the :math:`P_0` content alone.
        """
        out = np.zeros_like(self.values)
        out[0, 0] = self.values[0, 0]
        return HarmonicMomentField(out, self.mesh, self.L)

    def anisotropic_part(self) -> "HarmonicMomentField":
        r"""Return the :math:`\ell \ge 1` (anisotropic) projection.

        Same shape as ``self``; the :math:`\ell = 0, m = 0` slot
        zeroed.  Pairs with :meth:`isotropic_part` to partition the
        moment field; ``self.values == isotropic_part().values +
        anisotropic_part().values`` bit-exactly.

        Mirrors the ``skip_l0`` pattern in
        :class:`~orpheus.sn.scattering.LegendreMomentScattering`.
        """
        out = self.values.copy()
        out[0, 0] = 0.0
        return HarmonicMomentField(out, self.mesh, self.L)

    def scalar_flux(self) -> "ScalarFlux":
        r"""Extract the isotropic moment as a :class:`ScalarFlux`.

        Under the no-prefactor SH convention used by
        :func:`~orpheus.numerics.spherical_harmonics.evaluate_real_sh`
        (where :math:`Y_0^0 = 1`), the addition-theorem moment
        :math:`\phi_0^0 = \sum_n w_n Y_0^0 \psi_n = \sum_n w_n
        \psi_n` IS the scalar flux directly — no
        :math:`1/Y_0^0` factor.  This identity is what makes
        ``HarmonicMomentProjection(\psi).scalar_flux()`` agree with
        ``\psi.integrate_angular()`` bit-exactly.

        Returns
        -------
        ScalarFlux
            The :math:`(\ell=0, m=0)` slice ``values[0, 0]``, wrapped
            with the same mesh.
        """
        from .scalar_flux import ScalarFlux
        # values[0, 0] is the per-cell, per-group isotropic moment
        # shape (ng, nx, ny) — the principled ScalarFlux storage.
        return ScalarFlux(self.values[0, 0].copy(), self.mesh)

    def truncate(self, L_new: int) -> "HarmonicMomentField":
        r"""Return a new :class:`HarmonicMomentField` truncated to
        :math:`L_{\rm new} \le L`.

        Drops the :math:`\ell > L_{\rm new}` blocks and the
        corresponding zero-padded :math:`m`-columns; result has
        shape ``(L_new+1, 2*L_new+1, ng, nx, ny)``.

        Parameters
        ----------
        L_new : int
            Target order, must satisfy ``0 <= L_new <= self.L``.
        """
        if L_new > self.L:
            raise ValueError(
                f"HarmonicMomentField.truncate: L_new={L_new} > "
                f"current L={self.L}"
            )
        if L_new < 0:
            raise ValueError(
                f"HarmonicMomentField.truncate: L_new={L_new} < 0"
            )
        new_shape = (
            L_new + 1, 2 * L_new + 1,
            self.mesh.ng, self.mesh.nx, self.mesh.ny,
        )
        new_values = np.zeros(new_shape, dtype=self.values.dtype)
        # Copy the legitimate (l, m) entries; zero-padded slots outside
        # |m|<=l in the source carry zero so the slice copy is safe.
        new_values[: L_new + 1, : 2 * L_new + 1] = (
            self.values[: L_new + 1, : 2 * L_new + 1]
        )
        return HarmonicMomentField(new_values, self.mesh, L_new)

    # ── Dunder algebra (within-type) ──────────────────────────────────

    def __add__(self, other: "HarmonicMomentField") -> "HarmonicMomentField":
        self._validate_partner(other)
        return HarmonicMomentField(
            self.values + other.values, self.mesh, self.L,
        )

    def __sub__(self, other: "HarmonicMomentField") -> "HarmonicMomentField":
        self._validate_partner(other)
        return HarmonicMomentField(
            self.values - other.values, self.mesh, self.L,
        )

    def __mul__(self, scalar: float) -> "HarmonicMomentField":
        return HarmonicMomentField(
            self.values * float(scalar), self.mesh, self.L,
        )

    __rmul__ = __mul__

    def __truediv__(self, scalar: float) -> "HarmonicMomentField":
        return HarmonicMomentField(
            self.values / float(scalar), self.mesh, self.L,
        )

    def __neg__(self) -> "HarmonicMomentField":
        return HarmonicMomentField(-self.values, self.mesh, self.L)

    def _validate_partner(self, other: "HarmonicMomentField") -> None:
        if not isinstance(other, HarmonicMomentField):
            raise TypeError(
                "HarmonicMomentField arithmetic requires a "
                "HarmonicMomentField partner; "
                f"got {type(other).__name__}"
            )
        if other.mesh is not self.mesh:
            raise ValueError(
                "HarmonicMomentField arithmetic across distinct "
                "SNMesh instances is forbidden — the field is "
                "mesh-bound."
            )
        if other.L != self.L:
            raise ValueError(
                "HarmonicMomentField arithmetic requires matching "
                f"L; got self.L={self.L}, other.L={other.L}."
            )

    # ── Diagnostics ───────────────────────────────────────────────────

    def linf(self) -> float:
        r"""Return :math:`\lVert\phi\rVert_\infty` over all entries."""
        return float(np.abs(self.values).max())

    def copy(self) -> "HarmonicMomentField":
        """Return a deep copy carrying an owned ndarray."""
        return HarmonicMomentField(self.values.copy(), self.mesh, self.L)

    # ── Metadata read-throughs ────────────────────────────────────────

    @property
    def ng(self) -> int:
        """Energy group count."""
        return self.mesh.ng

    @property
    def nx(self) -> int:
        """Spatial extent in x."""
        return self.mesh.nx

    @property
    def ny(self) -> int:
        """Spatial extent in y."""
        return self.mesh.ny
