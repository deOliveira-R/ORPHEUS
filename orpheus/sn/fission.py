r"""Multigroup fission source operator as a :class:`LinearOperator`.

This module owns the **fission source operator** :math:`F` from the
operator-algebra view of the Boltzmann transport equation,

.. math::

    (L - S - F)\,\psi = q
    \qquad\text{(fixed source)}

.. math::

    (L - S)\,\psi = \tfrac{1}{k}\,F\,\psi
    \qquad\text{(eigenvalue)}

where :math:`L` is the streaming + collision operator, :math:`S` is
the scattering operator (see :mod:`orpheus.sn.scattering`), and
:math:`F` is the **fission emission operator**

.. math::

    (F\,\phi)_g(\vec r) = \chi_g\,
    \sum_{g'} \nu\Sigma_{f,g'}(\vec r)\,\phi_{g'}(\vec r).

The structure is **rank-1 in energy**: the production rate
:math:`\sum_{g'} \nu\Sigma_{f,g'}\phi_{g'}` is a scalar per cell and
the emission spectrum :math:`\chi_g` redistributes it across groups via
an outer product. This rank-1 factorisation is the algebraic origin of
why power iteration converges geometrically in a critical reactor
(dominance ratio :math:`|k_1/k_0|`) — the eigenvalue problem
:math:`A^{-1}F\phi = k\phi` has :math:`F` of rank
:math:`O(N_{\rm cells})` per group, but only one **per cell** in
energy because :math:`\chi` is a rank-1 spectrum.

Per Cardinal Rule 2 (architecture) this lifts the
``SNSolver.compute_fission_source`` math out of the solver and into a
single operator object. The math is **moved verbatim** (Wave D Issue 13
is a bit-identical extraction). The eigenvalue division by :math:`k`
**stays at the solver level** — :meth:`apply` returns
:math:`F\,\phi`, NOT :math:`F\,\phi / k`. The caller (the EigenvalueSolver
Protocol's ``compute_fission_source``, which is :class:`SNSolver`'s
delegator) divides by :math:`k`. Two reasons:

1. The :class:`LinearOperator` Protocol contract is *linear*:
   :meth:`apply` returns :math:`F\,\phi`, not :math:`F\,\phi/k`. The
   :math:`1/k` factor is a scalar multiple in the algebra
   :math:`F^{eff} = (1/k)\,F` — express it via :class:`ScaledOperator`
   if the algebra needs it explicitly.

2. The eigenvalue iteration owns :math:`k`. The fission operator
   should not be re-built every outer iteration just because :math:`k`
   changed; the operator state is :math:`(\chi, \nu\Sigma_f)`, which
   is fixed across the outer iteration.

Capability advertisement
========================

:pydata:`capabilities = frozenset({CAP_APPLY})`. No ``solve``: the
operator is **rank-1 in energy** per cell — its inverse does not exist.
No ``apply_transpose``: the adjoint fission operator is
:math:`F^T\,\phi^* = \nu\Sigma_f \cdot (\chi \cdot \phi^*)`, structurally
distinct (it sums over the **adjoint** energy distribution); the
current ORPHEUS solver does not consume :math:`F^T`. Both can be added
in a future wave when the adjoint transport machinery lands.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from functools import singledispatchmethod
from typing import TYPE_CHECKING

import numpy as np

from orpheus.numerics.operator import (
    CAP_APPLY,
    LinearOperatorMixin,
)

# Runtime imports for :func:`singledispatchmethod.register` — see
# ``scattering.py`` for the same pattern.  These three types form a
# leaf in the SN dependency graph (they do not import fission.py).
from .angular_flux import AngularFlux
from .scalar_flux import ScalarFlux
from .sources import IsotropicSource

if TYPE_CHECKING:
    from .material_xs_field import MaterialXSField


__all__ = ["FissionOperator"]


@dataclass
class FissionOperator(LinearOperatorMixin):
    r"""Fission source operator :math:`F = \chi\,\otimes\,\nu\Sigma_f`.

    Reads :math:`\chi(\vec r)` and :math:`\nu\Sigma_{f,g}(\vec r)`
    through a :class:`MaterialXSField` (Issue #197 PR-TYPED-1) — the
    same per-cell typed views every other operator (L, C, S)
    consumes.  The action is a contraction over groups (the production
    rate) followed by a broadcast across the emission spectrum.

    Use :meth:`from_solver_data` to build instances; pass
    ``mat_xs=sn_mesh.material_xs_field()``.

    Attributes
    ----------
    mat_xs : MaterialXSField
        Macroscopic XS field carrying ``emission_spectrum`` (χ) and
        ``fission_production`` (νΣ_f) per-cell views.
    capabilities : frozenset[str]
        ``{"apply"}`` — the rank-1 structure forbids a useful inverse,
        and the adjoint surface is not yet a consumer.
    """

    mat_xs: "MaterialXSField"

    capabilities: frozenset[str] = field(
        default_factory=lambda: frozenset({CAP_APPLY})
    )

    # ── Read-through properties for backwards compatibility ────────────

    @property
    def chi(self) -> np.ndarray:
        r""":math:`\chi(\vec r)` per-cell view, shape ``(ng, nx, ny)``.

        Read-through onto :attr:`mat_xs.emission_spectrum`.
        """
        return self.mat_xs.emission_spectrum

    @property
    def sig_p(self) -> np.ndarray:
        r""":math:`\nu\Sigma_f(\vec r)` per-cell view, shape ``(ng, nx, ny)``.

        Read-through onto :attr:`mat_xs.fission_production`.
        """
        return self.mat_xs.fission_production

    @classmethod
    def from_solver_data(
        cls, *, mat_xs: "MaterialXSField",
    ) -> "FissionOperator":
        """Construct from a :class:`MaterialXSField`.

        Issue #197 PR-TYPED-1 — the constructor surface collapses the
        ``(chi, sig_p)`` ndarray pair into one :class:`MaterialXSField`
        handle that carries both views consistently with the rest of
        the four-operator algebra.
        """
        return cls(mat_xs=mat_xs)

    @singledispatchmethod
    def apply(self, phi):
        r"""Apply :math:`F\,\phi`: emission rate × spectrum, no :math:`1/k`.

        .. math::

            (F\,\phi)_g(\vec r) = \chi_g(\vec r)\,
              \sum_{g'} \nu\Sigma_{f,g'}(\vec r)\,\phi_{g'}(\vec r)

        Dispatched on input type via
        :func:`functools.singledispatchmethod`:

        * :class:`~orpheus.sn.angular_flux.AngularFlux` →
          :class:`AngularFlux` — fission emission in **per-ordinate
          magnitude** (the trailing :math:`1/W` projection lives at the
          producer boundary; R-1 Step 4 A1).  Consumers are SN
          sweep / GMRES / source-iteration loops operating on per-
          ordinate angular flux.
        * :class:`~orpheus.sn.scalar_flux.ScalarFlux` →
          :class:`~orpheus.sn.sources.IsotropicSource` — fission
          emission in **iso scalar magnitude**.  For consumers in
          scalar-flux equations (eigenvalue outer / depletion /
          diffusion) that do not project to per-ordinate.
          Symmetric with :meth:`ScatteringOperator.apply`'s
          :class:`ScalarFlux` variant.
        * :class:`numpy.ndarray` — legacy bare-ndarray path,
          ``(ng, nx, ny)`` iso scalar fission source.  Preserved for
          packed-vector / FD-matvec consumers (e.g. ``SNStreamingOperator``,
          legacy GMRES on ``L.apply``).

        See ``coding-elegance`` SKILL.md §"Convention crosswalk template"
        and lesson L18 for the Pattern 7 producer-side normalisation
        discipline.  Fission has no :math:`P_\ell` aniso component
        (the emission spectrum is isotropic by construction); the
        :math:`1/W` projection in the :class:`AngularFlux` variant is
        the only producer-side normalisation.

        Notes
        -----
        R-1 Step 3a introduced the :class:`AngularFlux` overload for
        the operator algebra ``(L + C - S - F).apply(ψ)``.  R-1 Step 4
        A1 added the producer-side :math:`1/W` projection.
        """
        raise TypeError(
            f"FissionOperator.apply: unsupported input type "
            f"{type(phi).__name__}; expected AngularFlux, ScalarFlux, "
            f"or numpy.ndarray.  Dispatch table is registered via "
            f"@singledispatchmethod."
        )

    @apply.register
    def _(self, psi: AngularFlux) -> "AngularFlux":
        r"""Typed AngularFlux variant — per-ordinate magnitude output.

        Math: reduce angular → scalar via :math:`\phi = \sum_n w_n
        \psi_n`, compute iso fission source :math:`F\phi`, then
        project to per-ordinate via :math:`/W` and broadcast across
        the :math:`N` ordinates.  The :math:`1/W` is the producer-side
        normalisation (R-1 Step 4 A1).

        .. note:: **Tactical-legacy dimensional convention**.

           Returns :class:`AngularFlux`, but the returned value is
           dimensionally a *source*
           (:class:`~orpheus.sn.sources.PerOrdinateSource` would be the
           honest type).  See
           :meth:`ScatteringOperator.apply` (AngularFlux variant) and
           `#205 <https://github.com/deOliveira-R/ORPHEUS/issues/205>`_
           for the cross-method field/role architecture that will
           split the dimensional roles cleanly.
        """
        phi_scalar = psi.integrate_angular()
        # Compute iso fission source via the ScalarFlux variant.
        fission_iso: IsotropicSource = self.apply(phi_scalar)
        # Project to per-ordinate: ``Q/W`` broadcast across N.
        # Reuses the canonical :class:`PerOrdinateSource.from_isotropic`
        # factory — single source of truth for the iso→per-ord
        # projection (Pattern 2).
        from .sources import PerOrdinateSource
        per_ord = PerOrdinateSource.from_isotropic(
            fission_iso.values, psi.mesh,
        )
        return AngularFlux(per_ord.values, psi.mesh)

    @apply.register
    def _(self, phi: ScalarFlux) -> "IsotropicSource":
        r"""Typed ScalarFlux variant — iso scalar magnitude output.

        Math:
        :math:`Q_g(\vec r) = \chi_g(\vec r)\,
        \sum_{g'} \nu\Sigma_{f,g'}(\vec r)\,\phi_{g'}(\vec r)`.
        Iso scalar magnitude — no :math:`1/W` (scalar consumers do
        not project).

        Returns :class:`IsotropicSource` (R-1 Step 4 A1 — was
        :class:`ScalarFlux` pre-A1; the return type now matches
        :meth:`ScatteringOperator.apply`'s ScalarFlux variant by
        symmetry, and reflects the dimensional truth that the
        fission output is a source quantity, not a flux).
        """
        sig_p = self.sig_p
        chi = self.chi
        # Production rate: per-cell sum over groups, shape (nx, ny).
        fission_rate = np.einsum("gxy,gxy->xy", sig_p, phi.values)
        # Spectrum × rate: chi is (ng, nx, ny), rate is (nx, ny);
        # broadcast rate across the leading group axis.
        out = chi * fission_rate[None, :, :]
        return IsotropicSource(out, phi.mesh)

    @apply.register
    def _(self, phi_arr: np.ndarray) -> np.ndarray:
        r"""Bare-ndarray legacy variant — iso scalar fission source.

        Shape contract: input ``(ng, nx, ny)`` scalar flux, output
        ``(ng, nx, ny)`` iso fission source.  Preserved for packed-
        vector / FD-matvec consumers (``SNStreamingOperator``, legacy
        GMRES on ``L.apply``).  No type wrapping; the bare path
        bypasses the type layer entirely.
        """
        sig_p = self.sig_p
        chi = self.chi
        # Production rate: per-cell sum over groups, shape (nx, ny).
        fission_rate = np.einsum("gxy,gxy->xy", sig_p, phi_arr)
        # Spectrum × rate: chi is (ng, nx, ny), rate is (nx, ny);
        # broadcast rate across the leading group axis.
        return chi * fission_rate[None, :, :]
