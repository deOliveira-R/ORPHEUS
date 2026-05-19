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
from typing import TYPE_CHECKING

import numpy as np

from orpheus.numerics.operator import (
    CAP_APPLY,
    LinearOperatorMixin,
)

if TYPE_CHECKING:
    from .angular_flux import AngularFlux
    from .material_xs_field import MaterialXSField
    from .scalar_flux import ScalarFlux


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

    def apply(
        self, phi: "np.ndarray | ScalarFlux | AngularFlux",
    ) -> "np.ndarray | ScalarFlux | AngularFlux":
        r"""Apply :math:`F\,\phi`: emission rate × spectrum, no :math:`1/k`.

        .. math::

            (F\,\phi)_g(\vec r) = \chi_g(\vec r)\,
              \sum_{g'} \nu\Sigma_{f,g'}(\vec r)\,\phi_{g'}(\vec r)

        Parameters
        ----------
        phi : np.ndarray or ScalarFlux or AngularFlux
            Flux.  Three input types are accepted:

            * Bare ``np.ndarray`` shape ``(ng, nx, ny)`` — returns bare
              ``np.ndarray`` (legacy packed-vector / FD-matvec consumers,
              SNStreamingOperator / GMRES).
            * :class:`~orpheus.sn.scalar_flux.ScalarFlux` — returns
              :class:`ScalarFlux` (operator reads as the math at the
              eigenvalue layer; PR-TYPED-2).
            * :class:`~orpheus.sn.angular_flux.AngularFlux` — returns
              :class:`AngularFlux` (R-1 Step 3 — operator-algebra
              consumer at the within-group layer; lifts via internal
              ``ψ → φ → F·φ → broadcast(N, ·, ·, ·)``).

        Returns
        -------
        np.ndarray or ScalarFlux or AngularFlux
            Fission source *without* the :math:`1/k` eigenvalue
            division.  Type matches the input (typed → typed, raw → raw).

        Notes
        -----
        R-1 Step 3a — the :class:`AngularFlux` overload is added for
        completeness of the operator algebra ``(L + C - S - F).apply(ψ)``.
        In ORPHEUS the within-group operator triple has ``F =
        ZeroOperator``; ``F.apply(AngularFlux)`` is consumed by the
        future Phase A :class:`KEigenvalue` outer iteration, NOT by
        the within-group inner solve.  The result's ``.boundary`` is
        the auto-allocated zero :class:`BoundaryFlux` — fission is a
        volumetric emission with no face-trace contribution.
        """
        from .angular_flux import AngularFlux
        from .scalar_flux import ScalarFlux
        # Typed AngularFlux lift — reduce angular → scalar, apply scalar
        # path, broadcast back across the ordinate axis.  The per-
        # ordinate fission source is the per-cell emission rate (no
        # 1/W normalization — same convention as ScatteringOperator's
        # typed branch returns Q_iso[None, :, :, :] broadcast).
        if isinstance(phi, AngularFlux):
            phi_scalar: "ScalarFlux" = phi.integrate_angular()
            fission_scalar: "ScalarFlux" = self.apply(phi_scalar)
            N = phi.mesh.quad.N
            values = np.broadcast_to(
                fission_scalar.values[None, :, :, :],
                (N, phi.ng, phi.nx, phi.ny),
            ).copy()
            return AngularFlux(values, phi.mesh)
        is_typed = isinstance(phi, ScalarFlux)
        values = phi.values if is_typed else phi
        sig_p = self.sig_p
        chi = self.chi
        # Production rate: per-cell sum over groups, shape (nx, ny).
        fission_rate = np.einsum("gxy,gxy->xy", sig_p, values)
        # Spectrum × rate: chi is (ng, nx, ny), rate is (nx, ny);
        # broadcast rate across the leading group axis.
        out = chi * fission_rate[None, :, :]
        if is_typed:
            return ScalarFlux(out, phi.mesh)
        return out
