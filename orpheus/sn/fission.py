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

import numpy as np

from orpheus.numerics.operator import (
    CAP_APPLY,
    LinearOperatorMixin,
)


__all__ = ["FissionOperator"]


@dataclass
class FissionOperator(LinearOperatorMixin):
    r"""Fission source operator :math:`F = \chi\,\otimes\,\nu\Sigma_f`.

    Holds the per-cell-flattened arrays the SN solver already builds in
    ``__init__`` from the per-material :class:`Mixture` data: the
    emission spectrum :math:`\chi(\vec r)` and the production
    cross-section :math:`\nu\Sigma_{f,g}(\vec r)`. The action is a
    contraction over groups (the production rate) followed by a
    broadcast across the emission spectrum.

    Use :meth:`from_solver_data` to build instances from the same
    precomputed structures :class:`SNSolver` already holds.

    Attributes
    ----------
    chi : np.ndarray
        Emission spectrum shape ``(ng, nx, ny)`` (Issue #196 PR-INDEX-3 —
        principled layout).  Per-cell because the per-material
        :math:`\chi` is broadcast onto the cell grid via
        :class:`assemble_cell_xs`.
    sig_p : np.ndarray
        Production cross-section :math:`\nu\Sigma_f` shape
        ``(ng, nx, ny)``.
    capabilities : frozenset[str]
        ``{"apply"}`` — the rank-1 structure forbids a useful inverse,
        and the adjoint surface is not yet a consumer.

    Notes
    -----
    The constructor receives ``(chi, sig_p)`` already shaped
    ``(ng, nx, ny)`` — the same per-cell tensors :class:`SNSolver`
    holds. No re-broadcasting happens here; the operator simply applies
    the action.
    """

    chi: np.ndarray  # (ng, nx, ny) — Issue #196 PR-INDEX-3
    sig_p: np.ndarray  # (ng, nx, ny)

    capabilities: frozenset[str] = field(
        default_factory=lambda: frozenset({CAP_APPLY})
    )

    @classmethod
    def from_solver_data(
        cls, *, chi: np.ndarray, sig_p: np.ndarray,
    ) -> "FissionOperator":
        """Construct from the precomputed structures held by :class:`SNSolver`.

        ``chi`` and ``sig_p`` are the per-cell arrays
        :class:`SNSolver` builds in ``__init__`` via
        :func:`assemble_cell_xs` — already broadcast from per-material
        :class:`Mixture` data onto the cell grid.
        """
        return cls(chi=np.asarray(chi), sig_p=np.asarray(sig_p))

    def apply(self, phi: np.ndarray) -> np.ndarray:
        r"""Apply :math:`F\,\phi`: emission rate × spectrum, no :math:`1/k`.

        .. math::

            (F\,\phi)_g(\vec r) = \chi_g(\vec r)\,
              \sum_{g'} \nu\Sigma_{f,g'}(\vec r)\,\phi_{g'}(\vec r)

        Parameters
        ----------
        phi : np.ndarray
            Scalar flux shape ``(ng, nx, ny)`` (Issue #196 PR-INDEX-4 —
            principled layout).

        Returns
        -------
        np.ndarray
            Fission source shape ``(ng, nx, ny)`` (principled — Issue
            #196 PR-INDEX-4) *without* the :math:`1/k` eigenvalue
            division. The caller divides by :math:`k` (see module
            docstring for why).

        Notes
        -----
        Issue #196 PR-INDEX-4: ``self.sig_p`` / ``self.chi`` AND ``phi``
        AND the return value all live in the principled ``(ng, nx, ny)``
        layout — energy ``g`` is the leading axis end-to-end.  No
        transposes; ``np.einsum`` names the per-cell production-rate
        contraction (Pattern 3) — ``fission_rate`` has units ``[1/s]``
        per cell.
        """
        # Production rate: per-cell sum over groups, shape (nx, ny).
        # Both sig_p and phi are (ng, nx, ny); einsum names the
        # contraction over the leading group axis explicitly.
        fission_rate = np.einsum("gxy,gxy->xy", self.sig_p, phi)
        # Spectrum × rate: chi is (ng, nx, ny), rate is (nx, ny);
        # broadcast rate across the leading group axis.
        return self.chi * fission_rate[None, :, :]
