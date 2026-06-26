r"""Real-spherical-harmonic moment field on a tensor-product space.

The L2 typed wrapper for :math:`\phi_\ell^m(\vec r, g)` — a moment
field that sits between :class:`~orpheus.transport.fields.angular_flux.AngularFlux`
and the scattering operator :math:`\Lambda` as the natural data carrier
of the :math:`R \cdot \Lambda \cdot M \cdot \psi` Galerkin pipeline. Its
**flux**-role sibling on the source/sink side is
:class:`~orpheus.transport.source_sinks.harmonic_moment_source_sink.HarmonicMomentSourceSink`
(:math:`\Lambda` maps one to the other — the role-changing edge of the
``(angular ⊗ moment) × (flux ⊗ source)`` carrier grid).

Stores coefficients in an ``(L+1, 2L+1, ng, nx, ny)`` ndarray, with the
addition-theorem-shifted :math:`m`-index where slot ``l + m`` holds the
:math:`(\ell, m)` entry; entries outside :math:`|m| \le \ell` are zero
by convention.

Migration status (Depth B step D-E; renamed in the Frame campaign P4)
====================================================================

This class moved from ``orpheus.sn.harmonic_moment_field`` to
``orpheus.transport.fields`` (Depth B step D-E — a moment field is a
method-agnostic transport concept, not an SN-specific one) and now
inherits from :class:`~orpheus.numerics.field.Field`. In the Frame
campaign **P4** it was renamed ``HarmonicMomentField`` →
``HarmonicMomentFlux`` (and the module ``harmonic_moment_field`` →
``harmonic_moment_flux``): it is a ``(FluxRole, MomentField)`` carrier,
so the rename makes the role token explicit and the whole moment family
uniformly greppable (``grep HarmonicMoment`` / ``grep Flux`` /
``grep SourceSink`` now each find every member, matching
``AngularFlux``/``AngularSourceSink``). The earlier ``orpheus.sn``
re-export shim was retired in the same pass (it had zero importers). The
D-E Field-inheritance migration:

* Drops the hand-coded dunder skeleton (Cardinal Rule 2 — the algebra
  is now inherited via :func:`dataclasses.replace`).
* Adds the ``space: FunctionSpace`` field; the canonical space is a
  :class:`~orpheus.numerics.space.TensorProductSpace` of the form
  :math:`\mathrm{SphericalHarmonicSpace}(L) \otimes
  \mathrm{CellGroupSpace}(ng, nx, ny)` — the **first
  TensorProductSpace consumer in a typed Field** (D-B's L1 primitive
  is now load-bearing).
* Keeps ``mesh: SNMesh`` as an additive field under ``TYPE_CHECKING``
  (same pattern as :class:`~orpheus.transport.fields.scalar_flux.ScalarFlux`).
* Preserves the mesh-identity strict semantic via a
  :meth:`_check_partner` override.
* The ``L`` parameter is encoded in ``space.shape`` (and queryable via
  the composition-tree walk per Issue #207); the redundant ``L`` field
  is kept as a top-level attribute for ergonomic access — equivalent
  to ``self.space.find_factor(SphericalHarmonicSpace).L`` but avoiding
  the traversal at hot-path read sites.
* Introduces :meth:`from_mesh_and_L` for ergonomic 3-arg construction
  (the kw_only constructor requires explicit ``space``; the classmethod
  derives the space from ``mesh`` and ``L``).

Why distinct from :class:`AngularFlux` / :class:`ScalarFlux`
============================================================

The moment field lives in **moment space**
(:math:`(L+1) \cdot (2L+1)` coefficients per cell + group); the angular
flux lives in **per-ordinate space** (:math:`N` directions per cell +
group). Cross-type addition between the two is undefined — Field's
Layer 1 class-identity gate (`coding-elegance` Pattern 4) rejects it
by construction. The legitimate route is
``moments = M.apply(psi)`` (projection) or ``psi = R.apply(moments)``
(reconstruction) via the
:mod:`~orpheus.numerics.projection` Galerkin pair.

Units (B.4 — declared as the ``UNITS`` class constant)
======================================================

:math:`[1/(\mathrm{cm^2 \cdot s})]` — the SAME as
:class:`~orpheus.transport.fields.scalar_flux.ScalarFlux`
(:data:`~orpheus.numerics.units.SCALAR_FLUX_UNITS`), NOT the angular-flux
units. A moment :math:`\phi_\ell^m = \sum_n w_n Y_\ell^m \psi_n` is
**angle-integrated**: the quadrature weights carry ``sr`` (they sum to
:math:`4\pi`), cancelling the ``sr`` of the angular flux, while the
:math:`Y_\ell^m` and addition-theorem :math:`(2\ell+1)` factors are
dimensionless (the :math:`(2\ell+1)` lives on the reconstruction
operator, the :math:`4\pi/(2\ell+1)` metric on the space — neither on the
stored value). The :math:`\ell=0` moment IS the scalar flux exactly
(:meth:`scalar_flux` returns ``values[0, 0]``). The earlier
``1/(cm²·s·sr·eV)`` label was **wrong** — it forgot the angular
integration. eV-free per the binned-energy convention; see
:mod:`orpheus.numerics.units`.

References
----------

* Lewis, E.E. & Miller, W.F. (1993). *Computational Methods of Neutron
  Transport*. ANS. §3.5 — spherical-harmonic moments of the angular
  flux.
* Depth B plan §3.3, §6 step D-E.
* Issue #207 — architectural pattern: composition queries traverse the
  tensor-product tree; ``space.find_factor(SphericalHarmonicSpace).L``
  is the composition-aware way to read the truncation order.
* Issue #197 PR-TYPED-4 — original typed-wrapper introduction (now
  superseded by this Field-inheriting form).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, ClassVar

import numpy as np
from numpy.typing import NDArray

from orpheus.numerics.units import SCALAR_FLUX_UNITS, Unit
from orpheus.transport.fields._bases import MomentField
from orpheus.transport.fields._flux_role import FluxRole

if TYPE_CHECKING:
    from orpheus.sn.geometry import SNMesh
    from orpheus.transport.fields.scalar_flux import ScalarFlux


__all__ = ["HarmonicMomentFlux"]


@dataclass(frozen=True, eq=False, kw_only=True, repr=False)
class HarmonicMomentFlux(FluxRole, MomentField):
    r"""Real-spherical-harmonic moment field :math:`\phi_\ell^m(\vec r, g)`.

    Parameters
    ----------
    values : NDArray
        Moment coefficients of shape ``(L+1, 2L+1, ng, nx, ny)``.
    space : FunctionSpace
        The function space this field lives on. Canonically a
        :class:`TensorProductSpace` of the form
        :math:`\mathrm{SphericalHarmonicSpace}(L) \otimes
        \mathrm{CellGroupSpace}`. Construction via
        :meth:`from_mesh_and_L` is the canonical path; direct kw-only
        construction is for callers that already hold a constructed
        space.
    mesh : SNMesh
        The SN phase-space carrier.
    L : int
        Maximum harmonic order retained. Determines the leading two
        axes' sizes: ``values.shape[:2] == (L+1, 2L+1)``. Encoded in
        ``space.shape`` AND kept as a top-level field for ergonomic
        hot-path read access (avoids a per-read composition-tree
        traversal).

    Notes
    -----
    Algebra is inherited from :class:`~orpheus.numerics.field.Field`
    (dunders ``+``, ``-``, unary ``-``, scalar ``*``, scalar ``/``,
    plus diagnostics ``linf``, ``l2``, ``inner_product``, ``copy``).
    The :meth:`_check_partner` override adds the SN-specific
    mesh-identity check on top of Field's class-and-space gate. The
    ``L`` match is implicit via the space check (different ``L`` values
    produce different ``SphericalHarmonicSpace`` shapes, so different
    ``space`` instances) but checked explicitly in :meth:`_check_partner`
    for a clearer error message at the L-mismatch site.

    Cross-class arithmetic with :class:`AngularFlux` / :class:`ScalarFlux`
    is forbidden by Field's Layer 1 gate (`coding-elegance` Pattern 4).
    The legitimate route is through the discrete spherical-harmonic
    :class:`~orpheus.numerics.frame.GalerkinFrame`'s analysis / reconstruction faces.
    """

    #: Dimensional identity (View-G, B.4): a moment is angle-integrated, so
    #: ``1/(cm²·s)`` — :data:`~orpheus.numerics.units.SCALAR_FLUX_UNITS`,
    #: shared with ``ScalarFlux`` (the ``ℓ=0`` moment IS the scalar flux).
    #: Same units, different class — the gate is class identity. See the
    #: "Units" section above and :mod:`orpheus.numerics.units`.
    UNITS: ClassVar[Unit] = SCALAR_FLUX_UNITS

    # ── Slicing / decomposition (Pattern 3 — named intermediates) ────

    def l_block(self, l: int) -> NDArray:
        r"""View of one :math:`\ell`-block, shape ``(2ℓ+1, ng, nx, ny)``.

        Returns the slice ``values[l, :2*l+1]`` — the legitimate
        :math:`m`-entries for that :math:`\ell` (the trailing
        zero-padding outside :math:`|m| \le \ell` is excluded). Use
        this to retire the explicit ``moments[l, :n_m][..., ix, iy]``
        slicing pattern (``coding-elegance`` Pattern 3).
        """
        if not 0 <= l <= self.L:
            raise ValueError(
                f"HarmonicMomentFlux.l_block: l={l} out of range "
                f"[0, {self.L}]"
            )
        return self.values[l, : 2 * l + 1]

    def isotropic_part(self) -> "HarmonicMomentFlux":
        r"""Return the :math:`\ell = 0` (isotropic) projection.

        Same shape as ``self``; all :math:`\ell \ge 1` blocks zeroed.
        Used by the foldable-vs-residual scattering split when the
        consumer wants the :math:`P_0` content alone.
        """
        out = np.zeros_like(self.values)
        out[0, 0] = self.values[0, 0]
        return HarmonicMomentFlux(
            values=out, space=self.space, mesh=self.mesh, L=self.L,
        )

    def anisotropic_part(self) -> "HarmonicMomentFlux":
        r"""Return the :math:`\ell \ge 1` (anisotropic) projection.

        Same shape as ``self``; the :math:`\ell = 0, m = 0` slot zeroed.
        Pairs with :meth:`isotropic_part` to partition the moment field;
        ``self.values == isotropic_part().values + anisotropic_part().values``
        bit-exactly.

        Mirrors the ``skip_l0`` pattern in
        :class:`~orpheus.transport.operators.scattering.LegendreMomentScattering`.
        """
        out = self.values.copy()
        out[0, 0] = 0.0
        return HarmonicMomentFlux(
            values=out, space=self.space, mesh=self.mesh, L=self.L,
        )

    def scalar_flux(self) -> "ScalarFlux":
        r"""Extract the isotropic moment as a :class:`ScalarFlux`.

        Under the no-prefactor SH convention used by
        :class:`~orpheus.numerics.basis.spherical_harmonic_basis.SphericalHarmonicBasis`
        (where :math:`Y_0^0 = 1`), the addition-theorem moment
        :math:`\phi_0^0 = \sum_n w_n Y_0^0 \psi_n = \sum_n w_n \psi_n`
        IS the scalar flux directly — no :math:`1/Y_0^0` factor. This
        identity is what makes the frame analysis face's
        :math:`\phi_0^0` moment agree with ``\psi.integrate_angular()``
        bit-exactly.

        Returns
        -------
        ScalarFlux
            The :math:`(\ell=0, m=0)` slice ``values[0, 0]``, wrapped
            with the same mesh.
        """
        from orpheus.transport.fields.scalar_flux import ScalarFlux
        # ``values[0, 0]`` is ``(ng, *spatial[, 2^d])`` — the φ̂ spatial-moment
        # axis (if any) rides on the ℓ=0 moment and is propagated as a TYPED
        # factor (#240 D5b-S3; ``spatial_moments`` is this field's stored width).
        return ScalarFlux.from_mesh(
            self.values[0, 0].copy(), self.mesh,
            spatial_moments=self.spatial_moments,
        )

    # ── Truncation ───────────────────────────────────────────────────

    def truncate(self, L_new: int) -> "HarmonicMomentFlux":
        r"""Return a new :class:`HarmonicMomentFlux` truncated to
        :math:`L_{\rm new} \le L`.

        Drops the :math:`\ell > L_{\rm new}` blocks and the
        corresponding zero-padded :math:`m`-columns; result has shape
        ``(L_new+1, 2*L_new+1, ng, nx, ny)``.

        Parameters
        ----------
        L_new : int
            Target order, must satisfy ``0 <= L_new <= self.L``.
        """
        if L_new > self.L:
            raise ValueError(
                f"HarmonicMomentFlux.truncate: L_new={L_new} > "
                f"current L={self.L}"
            )
        if L_new < 0:
            raise ValueError(
                f"HarmonicMomentFlux.truncate: L_new={L_new} < 0"
            )
        new_shape = (
            L_new + 1, 2 * L_new + 1,
            self.mesh.ng, *self.mesh.spatial_shape,
        )
        new_values = np.zeros(new_shape, dtype=self.values.dtype)
        new_values[: L_new + 1, : 2 * L_new + 1] = (
            self.values[: L_new + 1, : 2 * L_new + 1]
        )
        return HarmonicMomentFlux.from_mesh_and_L(
            new_values, self.mesh, L_new,
        )
