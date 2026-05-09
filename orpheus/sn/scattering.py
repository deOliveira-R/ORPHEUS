r"""Multigroup scattering source as a :class:`LinearOperator`.

This module owns the **scattering source operator** :math:`S` from the
operator-algebra view of the Boltzmann transport equation,

.. math::

    (L - S - F)\,\psi = q

where :math:`L` is the streaming + collision operator, :math:`F` is the
fission operator, and :math:`S` here aggregates **all secondary-emission
channels that depend on the in-cell flux**:

* **P0 isotropic in-scatter** :math:`\Sigma_s^0(g'\to g)\,\phi_{g'}`,
* **Pℓ Galerkin reconstruction on real spherical harmonics**
  :math:`Y_\ell^m` for :math:`\ell\ge 1` from the angular flux moments,
* **(n,2n) doubling** :math:`2\,\Sigma_{2n}(g'\to g)\,\phi_{g'}`, which
  ORPHEUS folds into the scattering side because its cell-by-cell
  bookkeeping is identical to in-scatter (vectorise-by-material,
  add-into-:math:`Q`).

The operator advertises only :pydata:`CAP_APPLY`. There is no useful
:meth:`solve`: the discrete :math:`S` is rank :math:`O(N_{\text{cells}}
\cdot N_{\text{groups}})` with no tractable inverse — it is *applied*,
never *inverted*. This is the canonical example the capability-set
design (Wave A Issue 1) was built for: forcing :math:`S` to provide a
``solve`` stub that raises ``NotImplementedError`` is harmful, so we
declare what :math:`S` **can** do and let composers refuse to build
:math:`S^{-1}` at construction time.

Per Cardinal Rule 2 (architecture) this lifts the
``SNSolver._add_scattering_source`` / ``SNSolver._build_aniso_scattering``
/ ``SNSolver._add_n2n_source`` math out of the solver and into a
single operator object. The math is **moved verbatim** (Wave D Issue 13
is a bit-identical extraction); the only change is the home of the code.
:class:`SNSolver` retains thin delegator methods so the EigenvalueSolver
Protocol surface and the ``test_solver_components.py`` probe-tests
keep working unchanged.

Mathematical content
====================

P0 in-scatter
-------------

Per material region :math:`m`, with :math:`\Sigma_{s,0}^m \in
\mathbb{R}^{G\times G}` the row-from / column-to scattering matrix at
order :math:`\ell = 0`,

.. math::

    Q^{(0)}_g(\vec r)
    = \sum_{g'} \Sigma_{s,0}^m(g'\to g)\,\phi_{g'}(\vec r)
    = \bigl(\phi(\vec r)\,\cdot\,\Sigma_{s,0}^m\bigr)_g

with :math:`\phi(\vec r) \in \mathbb{R}^G` and the dot product in
matrix-row order. This is isotropic: it appears identically in every
ordinate via the :math:`1/W = 1/\sum_n w_n` normalisation that the
sweep applies.

Pℓ Galerkin projection on :math:`Y_\ell^m`
------------------------------------------

For anisotropic scattering of order :math:`L \ge 1`, the angular flux
is projected onto the real spherical harmonics
:math:`Y_\ell^m(\hat\Omega)`,

.. math::

    \phi^{\ell m}_g(\vec r)
    = \int_{4\pi}\psi_g(\vec r,\hat\Omega)\,Y_\ell^m(\hat\Omega)\,
      d\hat\Omega
    \;\approx\; \sum_n w_n\,\psi_{n,g}(\vec r)\,Y_\ell^m(\hat\Omega_n)

(the right-hand side is the discrete quadrature; canonical equation
label :eq:`flux-moments` lives in :ref:`theory-discrete-ordinates`).
The :math:`\ell\ge 1` contribution to the per-ordinate scattering
source uses the addition theorem
:math:`\sum_m Y_\ell^m(\hat\Omega)\,Y_\ell^m(\hat\Omega') =
P_\ell(\hat\Omega\cdot\hat\Omega')`,

.. math::

    Q^{(\ell\ge 1)}_{n,g}(\vec r)
    = \sum_{\ell=1}^{L}\,(2\ell+1)\,\sum_m Y_\ell^m(\hat\Omega_n)\,
      \sum_{g'}\Sigma_{s,\ell}^m(g'\to g)\,\phi^{\ell m}_{g'}(\vec r).

The :math:`(2\ell+1)` factor is the discrete-orthogonality
normalisation
:math:`\langle Y_\ell^m | Y_{\ell'}^{m'}\rangle = (4\pi/(2\ell+1))\,
\delta_{\ell\ell'}\delta_{mm'}` working out across both projection
and reconstruction. The Galerkin frame is real spherical harmonics, not
complex — see :class:`orpheus.sn.quadrature.LebedevSphere`.

(n,2n) doubling
---------------

The (n,2n) reaction emits **two** neutrons per absorption, producing a
secondary source :math:`2\,\Sigma_{2n}^m(g'\to g)\,\phi_{g'}`. Because
this depends on the scalar flux only and adds isotropically (in this
discretisation), it folds into the scattering side of the algebra:
:meth:`apply` returns the sum of P0 + (n,2n) + Pℓ contributions.

Capability advertisement
========================

:pydata:`capabilities = frozenset({CAP_APPLY})`. No ``solve`` (rank
prevents efficient inversion); no ``apply_transpose`` (would require an
adjoint Pℓ transposition that the current ORPHEUS solver does not
need; can be added in a future wave).
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from orpheus.numerics.operator import (
    CAP_APPLY,
    LinearOperatorMixin,
)


__all__ = ["ScatteringOperator"]


@dataclass
class ScatteringOperator(LinearOperatorMixin):
    r"""Scattering source operator :math:`S` (P0 + Pℓ + (n,2n)).

    Holds the precomputed per-material Legendre scattering matrices
    :math:`\Sigma_{s,\ell}`, the (n,2n) matrices :math:`\Sigma_{2n}`,
    the precomputed real spherical harmonics :math:`Y_\ell^m` evaluated
    at every quadrature direction, and a per-material cell-index map
    for vectorised assembly.

    Use :meth:`from_solver_data` to build instances from the same
    precomputed structures :class:`SNSolver` already holds (the canonical
    construction path used by Wave D Issue 13's bit-identical
    extraction).

    Attributes
    ----------
    n_ordinates : int
        Number of angular ordinates :math:`N` in the quadrature.
    nx, ny, ng : int
        Spatial grid sizes and group count (the operator's action shape
        is ``(N, nx, ny, ng) -> (N, nx, ny, ng)``).
    scattering_order : int
        Maximum Legendre order :math:`L` retained. ``0`` means P0 only.
    sig_s : dict[int, list[np.ndarray]]
        Per-material list of ``(ng, ng)`` Legendre scattering matrices,
        one per Legendre order ``[0..L]``. ``sig_s[mid][l][g_from, g_to]``
        is the scattering cross-section at order :math:`\ell` from
        ``g_from`` to ``g_to``.
    sig2 : dict[int, np.ndarray]
        Per-material ``(ng, ng)`` (n,2n) cross-section matrices.
    sig_s0 : dict[int, np.ndarray]
        Convenience alias ``sig_s[mid][0]`` for the P0 fast path.
    Y : np.ndarray | None
        Real spherical harmonics evaluated at the quadrature ordinates,
        shape ``(N, L+1, 2L+1)``. ``None`` when ``L == 0``.
    weights : np.ndarray
        Quadrature weights ``(N,)``.
    cells_by_mat : dict[int, tuple[np.ndarray, np.ndarray]]
        Per-material ``(ix, iy)`` arrays for vectorised assembly.
    capabilities : frozenset[str]
        ``{"apply"}`` — :math:`S` has no efficient inverse and the
        current ORPHEUS algebra has no consumer for :math:`S^T`.

    Notes
    -----
    The constructor argument order mirrors how :class:`SNSolver` builds
    the equivalent state in ``__init__`` (Wave D Issue 13's bit-identical
    extraction). See the module docstring for the mathematical content.
    """

    n_ordinates: int
    nx: int
    ny: int
    ng: int
    scattering_order: int
    sig_s: dict[int, list[np.ndarray]]
    sig2: dict[int, np.ndarray]
    sig_s0: dict[int, np.ndarray]
    Y: np.ndarray | None
    weights: np.ndarray
    cells_by_mat: dict[int, tuple[np.ndarray, np.ndarray]]

    capabilities: frozenset[str] = field(
        default_factory=lambda: frozenset({CAP_APPLY})
    )

    @classmethod
    def from_solver_data(
        cls,
        *,
        n_ordinates: int,
        nx: int,
        ny: int,
        ng: int,
        scattering_order: int,
        sig_s: dict[int, list[np.ndarray]],
        sig2: dict[int, np.ndarray],
        sig_s0: dict[int, np.ndarray],
        Y: np.ndarray | None,
        weights: np.ndarray,
        cells_by_mat: dict[int, tuple[np.ndarray, np.ndarray]],
    ) -> "ScatteringOperator":
        """Construct from the precomputed structures held by :class:`SNSolver`.

        This is the canonical construction path: :class:`SNSolver`
        precomputes ``sig_s``, ``sig2``, ``sig_s0``, ``_Y`` (when
        :math:`L > 0`), ``_cells_by_mat`` and the angular-quadrature
        weights during ``__init__``, and then builds the
        :class:`ScatteringOperator` from those structures. Sharing the
        same precomputed handles is what lets the bit-identical
        extraction round-trip exactly.
        """
        return cls(
            n_ordinates=n_ordinates,
            nx=nx,
            ny=ny,
            ng=ng,
            scattering_order=scattering_order,
            sig_s=sig_s,
            sig2=sig2,
            sig_s0=sig_s0,
            Y=Y,
            weights=weights,
            cells_by_mat=cells_by_mat,
        )

    # ── In-place helpers (preserve bit-identity vs SNSolver pre-Wave-D) ─

    def add_iso_source(self, Q: np.ndarray, phi: np.ndarray) -> None:
        """Add P0 in-scatter source to ``Q`` in-place.

        Vectorised by material: ``Q[ix, iy, :] += phi[ix, iy, :] @
        sig_s[mid][0]`` for every cell of every material.

        Parameters
        ----------
        Q : np.ndarray
            Isotropic source array shape ``(nx, ny, ng)``. Modified in
            place.
        phi : np.ndarray
            Scalar flux shape ``(nx, ny, ng)``.
        """
        for mid, (ix, iy) in self.cells_by_mat.items():
            Q[ix, iy, :] += phi[ix, iy, :] @ self.sig_s0[mid]

    def add_n2n_source(self, Q: np.ndarray, phi: np.ndarray) -> None:
        """Add (n,2n) source to ``Q`` in-place.

        Vectorised by material: ``Q[ix, iy, :] += 2 * (phi[ix, iy, :] @
        sig2[mid])`` for every cell of every material.

        Parameters
        ----------
        Q : np.ndarray
            Isotropic source array shape ``(nx, ny, ng)``. Modified in
            place.
        phi : np.ndarray
            Scalar flux shape ``(nx, ny, ng)``.
        """
        for mid, (ix, iy) in self.cells_by_mat.items():
            Q[ix, iy, :] += 2.0 * (phi[ix, iy, :] @ self.sig2[mid])

    def build_aniso_source(
        self, angular_flux: np.ndarray | None,
    ) -> np.ndarray | None:
        r"""Build per-ordinate Pℓ (:math:`\ell \ge 1`) scattering source.

        Implements the Galerkin reconstruction :eq:`pn-scatter` from
        the angular-flux moments :eq:`flux-moments` (declared in
        :ref:`theory-discrete-ordinates`).

        Parameters
        ----------
        angular_flux : np.ndarray or None
            Angular flux shape ``(N, nx, ny, ng)`` from the most recent
            sweep, or ``None`` on the first source iteration before any
            sweep has been run.

        Returns
        -------
        np.ndarray or None
            ``(N, nx, ny, ng)`` per-ordinate :math:`\ell \ge 1`
            contribution, or ``None`` when ``scattering_order == 0`` or
            ``angular_flux is None``.
        """
        if self.scattering_order == 0 or angular_flux is None:
            return None

        N = self.n_ordinates
        nx, ny, ng = self.nx, self.ny, self.ng
        L = self.scattering_order
        Y = self.Y  # (N, L+1, 2L+1)
        w = self.weights

        # Compute Legendre moments: fiL[x, y, g, l, l+m] = Σ_n w_n ψ_n Y_l^m(n)
        fiL = np.zeros((nx, ny, ng, L + 1, 2 * L + 1))
        for l in range(L + 1):
            for m in range(-l, l + 1):
                # (N,) * (N, nx, ny, ng) summed over N → (nx, ny, ng)
                fiL[:, :, :, l, l + m] = np.einsum(
                    'n,nxyg->xyg', w * Y[:, l, l + m], angular_flux,
                )

        # Build anisotropic source: only l >= 1 terms (P0 is in Q_iso)
        Q_aniso = np.zeros((N, nx, ny, ng))
        for mid, (ix, iy) in self.cells_by_mat.items():
            sig_s_l = self.sig_s[mid]
            for l in range(1, L + 1):  # skip l=0 (handled by add_iso_source)
                # Σ_m fiL[..., l, m] * Y_l^m(n) → reconstruct angular moment at ordinate n
                for m in range(-l, l + 1):
                    moment = fiL[ix, iy, :, l, l + m]  # (n_cells, ng)
                    # (n_cells, ng) @ (ng, ng) → (n_cells, ng)
                    scattered = moment @ sig_s_l[l]  # Σ_s^l @ fiL_lm
                    for n in range(N):
                        Q_aniso[n, ix, iy, :] += (2 * l + 1) * Y[n, l, l + m] * scattered

        return Q_aniso

    # ── LinearOperator surface ─────────────────────────────────────────

    def apply(self, psi: np.ndarray) -> np.ndarray:
        r"""Apply the full scattering source operator :math:`S\,\psi`.

        Returns a per-ordinate isotropic-equivalent source plus the
        :math:`\ell \ge 1` per-ordinate contribution, packaged together
        as the :math:`(N, n_x, n_y, n_g)` array the sweep would receive
        from the legacy code path.

        For the P0 + (n,2n) parts, the action is the same scalar source
        broadcast across every ordinate (the sweep's ``1/W`` factor is
        not applied here — the caller, ``transport_sweep``, applies it).
        For :math:`\ell \ge 1`, the action carries genuine per-ordinate
        directional content.

        Parameters
        ----------
        psi : np.ndarray
            Angular flux shape ``(N, nx, ny, ng)``.

        Returns
        -------
        np.ndarray
            Per-ordinate scattering source shape
            ``(N, nx, ny, ng)``: the isotropic P0 + (n,2n) part is
            broadcast across the ``N`` axis; the :math:`\ell \ge 1`
            part is the genuine per-ordinate Galerkin reconstruction.

        Notes
        -----
        ORPHEUS's existing :class:`SNSolver` source iteration consumes
        the isotropic and anisotropic pieces separately (the isotropic
        :math:`Q` is a ``(nx, ny, ng)`` scalar source; the anisotropic
        :math:`Q_{\rm aniso}` is a ``(N, nx, ny, ng)`` per-ordinate
        source). For internal use by :class:`SNSolver` the helpers
        :meth:`add_iso_source`, :meth:`add_n2n_source`, and
        :meth:`build_aniso_source` expose the same math without forcing
        a wasteful broadcast. :meth:`apply` is the LinearOperator
        Protocol surface and is intended for downstream operator-algebra
        consumers (Krylov-on-apply, adjoint sweeps, composition with
        :math:`L` and :math:`F`).
        """
        # Reduce angular flux to scalar flux: φ = Σ_n w_n ψ_n
        # (Equivalent to angular_flux_to_scalar from operator.py.)
        phi = np.einsum('n,nxyg->xyg', self.weights, psi)

        # Isotropic (P0 + (n,2n)) — vectorised by material into a fresh
        # (nx, ny, ng) buffer, then broadcast across the N axis.
        Q_iso = np.zeros((self.nx, self.ny, self.ng))
        self.add_iso_source(Q_iso, phi)
        self.add_n2n_source(Q_iso, phi)

        # Broadcast across ordinates and add the Pℓ (l>=1) contribution.
        Q = np.broadcast_to(Q_iso[None, :, :, :], psi.shape).copy()
        Q_aniso = self.build_aniso_source(psi)
        if Q_aniso is not None:
            Q += Q_aniso
        return Q
