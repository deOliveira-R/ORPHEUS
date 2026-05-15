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
from orpheus.numerics.projection import (
    HarmonicMomentProjection,
    HarmonicMomentReconstruction,
)


__all__ = ["LegendreMomentScattering", "ScatteringOperator"]


@dataclass(frozen=True)
class LegendreMomentScattering(LinearOperatorMixin):
    r"""Per-ℓ block-diagonal scattering on harmonic-moment space.

    The §15.2 / §10 sum-of-tensor-products form

    .. math::

        \Lambda \;=\; \sum_{\ell=1}^{L}\, \mathbf{P}_\ell \otimes \Sigma_{s,\ell}

    where :math:`\mathbf{P}_\ell` selects the :math:`\ell`-th harmonic
    block and :math:`\Sigma_{s,\ell}` is the per-material per-ℓ
    Legendre cross-section matmul on the energy-group axis.

    For an input moment field
    :math:`\phi_\ell^m(\vec r) \in \mathbb{R}^{(L+1) \times (2L+1) \times N_g \times N_x \times N_y}`
    (Issue #196 PR-INDEX-4 — principled ``g`` leading-after-moments),
    the action is

    .. math::

        (\Lambda \phi)_\ell^m(\vec r)\bigg|_{g}
        \;=\; \sum_{g'} \Sigma_{s,\ell}^{m(\vec r)}(g' \to g)\,
              \phi_\ell^m(\vec r)\bigg|_{g'},

    where :math:`m(\vec r)` is the material id at cell :math:`\vec r`
    (the per-material structure is folded into the cell axis via
    ``cells_by_mat``).

    The :math:`\ell = 0` block is **skipped by default** because the
    project's P0 isotropic in-scatter is handled by the separate
    :meth:`ScatteringOperator.add_iso_source` fast path (a reaction-rate
    formulation that does not go through the moment tensor). Setting
    ``skip_l0 = False`` includes the :math:`\ell = 0` block (used by
    the LinearOperator-surface :meth:`ScatteringOperator.apply` when a
    full :math:`R \Lambda M \psi` composition is required).

    Implementation
    --------------

    The two outer Python loops are over **materials** (1-10 typically)
    and **Legendre order ℓ** (0-3 typically) — both small structural
    dimensions, NOT smoking guns. Inside each ``(mid, ℓ)`` pair, ONE
    :func:`numpy.einsum` handles all :math:`m`, cells, and groups in
    a single contraction. The total flop count is identical to the
    legacy hand-rolled version inside
    :meth:`ScatteringOperator.build_aniso_source`; the iteration over
    ordinates that *was* there is no longer needed because moment-space
    scattering does not couple ordinates.

    Capability set: :pydata:`{"apply"}`. There is no efficient
    :meth:`solve` (the moment-space scattering matrix is rank-deficient
    on the :math:`\ell = 0` block by design); :meth:`apply_transpose`
    is not currently consumed by any ORPHEUS solver and is deferred.

    Parameters
    ----------
    sig_s : dict[int, list[np.ndarray]]
        Per-material list of ``(ng, ng)`` Legendre scattering matrices,
        one per Legendre order ``[0..L]``. ``sig_s[mid][l][g_from, g_to]``
        is the scattering cross-section at order :math:`\ell` from
        ``g_from`` to ``g_to``.
    cells_by_mat : dict[int, tuple[np.ndarray, np.ndarray]]
        Per-material ``(ix, iy)`` index arrays for vectorised assembly
        into the cell axis.
    L : int
        Maximum Legendre order :math:`L` retained.
    skip_l0 : bool, default ``True``
        Skip the :math:`\ell = 0` block (handled separately by P0
        in-scatter). Set to ``False`` for the full LinearOperator
        composition :math:`R \Lambda M \psi`.
    """

    sig_s: dict[int, list[np.ndarray]]
    cells_by_mat: dict[int, tuple[np.ndarray, np.ndarray]]
    L: int
    skip_l0: bool = True
    capabilities: frozenset = field(
        default_factory=lambda: frozenset({CAP_APPLY})
    )

    def apply(self, moments: np.ndarray) -> np.ndarray:
        r"""Apply :math:`\Lambda` to a moment field.

        Parameters
        ----------
        moments : np.ndarray
            Moment field of shape ``(L+1, 2L+1, ng, nx, ny)`` (Issue
            #196 PR-INDEX-4 — principled).  The :math:`m`-axis is the
            addition-theorem-shifted index where slot ``l + m`` holds
            the :math:`(\ell, m)` entry; entries outside
            :math:`|m| \le \ell` are conventionally zero.

        Returns
        -------
        np.ndarray
            Same shape as ``moments``.  The :math:`\ell = 0` block is
            zero when ``skip_l0`` is ``True``; otherwise the P0 in-scatter
            contribution is included.
        """
        out = np.zeros_like(moments)
        l_start = 1 if self.skip_l0 else 0
        for mid, (ix, iy) in self.cells_by_mat.items():
            sig_s_mid = self.sig_s[mid]
            for l in range(l_start, self.L + 1):
                # Only the 2ℓ+1 valid m-slots are physical at order ℓ;
                # the remaining slots (2ℓ+1..2L) are zero by convention
                # in the moment tensor produced by HarmonicMomentProjection.
                # Restricting the einsum to the valid range matches the
                # legacy ``for m in range(-l, l+1)`` iteration in
                # `build_aniso_source` and avoids scattering out-of-band
                # garbage when the operator is fed a synthetic random input.
                n_m = 2 * l + 1
                # Indexing pattern: ``moments[l, :n_m][..., ix, iy]``
                # keeps advanced indices CONTIGUOUS and at the trailing
                # position → numpy preserves axis order, output shape
                # ``(n_m, ng, n_cells)``.  If we instead wrote
                # ``moments[l, :n_m, :, ix, iy]`` (advanced indices
                # separated from ``:n_m`` by ``:``), numpy would move
                # the advanced-index axis to the FRONT (output
                # ``(n_cells, n_m, ng)``) — a silent shape rearrangement
                # that breaks the einsum below.  See numpy advanced-
                # indexing rules.
                # moments_view shape: (n_m, ng, n_cells).
                # sig_s_mid[l] shape: (ng, ng), indexed [g_from, g_to].
                # einsum: out_g_to = Σ_g_from moments_g_from · sig_s[g_from, g_to]
                moments_view = moments[l, :n_m][..., ix, iy]
                out_block = np.einsum(
                    "mfc,fg->mgc", moments_view, sig_s_mid[l],
                )                                       # (n_m, ng, n_cells)
                # Same trailing-contiguous indexing on out for the
                # accumulation.  Use the long form ``out += np.einsum``
                # via temporary because numpy fancy-index in-place
                # assignment is well-defined for this pattern.
                out[l, :n_m][..., ix, iy] = (
                    out_block + out[l, :n_m][..., ix, iy]
                )
        return out


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
        is ``(N, ng, nx, ny) -> (N, ng, nx, ny)`` under PR-INDEX-4 —
        principled).
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

        Vectorised by material: per cell ``c`` of material ``mid``,
        ``Q[:, ic, jc] += sig_s0[mid].T @ phi[:, ic, jc]`` where
        ``sig_s0[mid]`` is ``(ng, ng)`` indexed ``[g_from, g_to]``.

        Parameters
        ----------
        Q : np.ndarray
            Isotropic source array shape ``(ng, nx, ny)`` (Issue #196
            PR-INDEX-4 — principled layout). Modified in place.
        phi : np.ndarray
            Scalar flux shape ``(ng, nx, ny)``.

        Notes
        -----
        ``np.einsum("fg,fxy->gxy", sig_s0[mid], phi[:, ix, iy])`` names
        the [g_from → g_to] contraction explicitly — the in-scatter
        source per cell is the source-spectrum-to-sink-spectrum
        contraction of the per-material P0 matrix with the per-cell
        scalar flux (Pattern 3).
        """
        for mid, (ix, iy) in self.cells_by_mat.items():
            # phi[:, ix, iy] shape (ng, n_cells); sig_s0[mid] shape
            # (g_from, g_to). Source-to-sink contraction over g_from.
            Q[:, ix, iy] += np.einsum(
                "fg,fc->gc", self.sig_s0[mid], phi[:, ix, iy],
            )

    def add_n2n_source(self, Q: np.ndarray, phi: np.ndarray) -> None:
        """Add (n,2n) source to ``Q`` in-place.

        Vectorised by material: per cell ``c`` of material ``mid``,
        ``Q[:, ic, jc] += 2 · sig2[mid].T @ phi[:, ic, jc]``.

        Parameters
        ----------
        Q : np.ndarray
            Isotropic source array shape ``(ng, nx, ny)`` (Issue #196
            PR-INDEX-4 — principled). Modified in place.
        phi : np.ndarray
            Scalar flux shape ``(ng, nx, ny)``.
        """
        for mid, (ix, iy) in self.cells_by_mat.items():
            Q[:, ix, iy] += 2.0 * np.einsum(
                "fg,fc->gc", self.sig2[mid], phi[:, ix, iy],
            )

    def build_aniso_source(
        self, angular_flux: np.ndarray | None,
    ) -> np.ndarray | None:
        r"""Build per-ordinate Pℓ (:math:`\ell \ge 1`) scattering source.

        Implements the Galerkin reconstruction :eq:`pn-scatter` from
        the angular-flux moments :eq:`flux-moments` (declared in
        :ref:`theory-discrete-ordinates`) as the **literal operator
        composition**

        .. math::

            Q^{\rm aniso}_n(\vec r) \;=\; (R \, \Lambda \, M \, \psi)_n(\vec r)

        of the three primitives shipped in
        :mod:`orpheus.numerics.projection` and
        :class:`LegendreMomentScattering` (Wave 1 / C1.1):

        * :math:`M` :class:`~orpheus.numerics.projection.HarmonicMomentProjection`:
          :math:`\psi(N, \cdot) \mapsto \phi^{\ell m}(L+1, 2L+1, \cdot)`
          via :math:`\phi_\ell^m = \sum_n w_n Y_\ell^m(\Omega_n) \psi_n`.
        * :math:`\Lambda` :class:`LegendreMomentScattering`: per-ℓ
          block-diagonal cross-section action on moment space (skip ℓ=0,
          handled by P0 in-scatter).
        * :math:`R` :class:`~orpheus.numerics.projection.HarmonicMomentReconstruction`:
          :math:`\phi^{\ell m} \mapsto q_n` via the addition-theorem
          reconstruction :math:`q_n = \sum_\ell (2\ell+1) \sum_m
          Y_\ell^m(\Omega_n) \phi_\ell^m`.

        The :math:`(2\ell+1)` factor in :math:`R` is the addition-theorem
        scaling (cf. Wave 0 ERR-039: :math:`R` is the addition-theorem
        reconstruction, NOT the W-weighted Hilbert adjoint of :math:`M`,
        which differs by exactly this factor).

        Total flop count is identical to the legacy hand-rolled
        ``for n in range(N)`` triple loop; the iteration over ordinates
        is now internal to :func:`numpy.einsum` inside :math:`R` and
        :math:`M`, **not** a Python loop.

        Parameters
        ----------
        angular_flux : np.ndarray or None
            Angular flux shape ``(N, ng, nx, ny)`` (Issue #196
            PR-INDEX-4 — principled) from the most recent sweep, or
            ``None`` on the first source iteration before any sweep
            has been run.

        Returns
        -------
        np.ndarray or None
            ``(N, ng, nx, ny)`` per-ordinate :math:`\ell \ge 1`
            contribution, or ``None`` when ``scattering_order == 0`` or
            ``angular_flux is None``.

        Notes
        -----
        :class:`HarmonicMomentProjection` and
        :class:`HarmonicMomentReconstruction` are layout-agnostic in
        their trailing axes — only the leading ordinate / harmonic axes
        are consumed by the einsums.  Switching ψ from ``(N, nx, ny,
        ng)`` to ``(N, ng, nx, ny)`` therefore transparently produces a
        moment field shape ``(L+1, 2L+1, ng, nx, ny)`` (instead of
        ``(L+1, 2L+1, nx, ny, ng)``), which is exactly what the
        principled :class:`LegendreMomentScattering` consumes after
        PR-INDEX-4.
        """
        if self.scattering_order == 0 or angular_flux is None:
            return None
        L = self.scattering_order
        # Build the §9 "S = R Λ M" pipeline. The constituent primitives
        # are cheap dataclass instantiations; the actual work is in the
        # three np.einsum calls inside their .apply methods.
        M = HarmonicMomentProjection(weights=self.weights, Y=self.Y, L=L)
        Lam = LegendreMomentScattering(
            sig_s=self.sig_s,
            cells_by_mat=self.cells_by_mat,
            L=L,
            skip_l0=True,
        )
        R = HarmonicMomentReconstruction.from_Y(self.Y)
        return R.apply(Lam.apply(M.apply(angular_flux)))

    # ── Foldable / residual split ─────────────────────────────────────
    #
    # The Phase G four-operator unification splits scattering as
    #
    #     S = S_foldable + S_residual
    #
    # where ``S_foldable`` is the **P0 within-group self-scatter**
    # only: per ordinate the source is :math:`\Sigma_{s,0}^{g\to g}
    # \phi_g`, direction-independent, so the cell-balance denominator
    # can absorb it as :math:`\sigma_r = \sigma_t -
    # \Sigma_{s,0}^{g\to g}` (the "removal" cross-section). Every
    # other piece of :math:`S` lives in ``S_residual``:
    #
    # * **Cross-group P0** (off-diagonal of ``sig_s[mid][0]``) — couples
    #   distinct energy groups, cannot collapse into a per-cell scalar.
    # * **All :math:`P_\ell \ge 1`** — direction-dependent through
    #   :math:`Y_\ell^m(\Omega_n)`, fundamentally unfoldable into
    #   :math:`\sigma_r`. This is NOT an SN-specific quirk; it is
    #   inherent to anisotropic scattering's angular structure for
    #   every transport solver family (SN, CP, MoC, diffusion-P_N).
    # * **(n,2n) doubling** — emits two neutrons per absorption;
    #   folding into a "removal" cross-section :math:`\sigma_r` is
    #   conceptually wrong even when ``sig2[mid][g, g]`` is non-zero.
    #
    # Both halves are :class:`ScatteringOperator` siblings (same class,
    # different per-material arrays). No new wrapper type. The split
    # is purely additive: ``S.apply(\psi) =
    # S.foldable_part().apply(\psi) + S.residual_part().apply(\psi)``
    # at FP-non-associativity precision (typically ``rtol=1e-14``).
    #
    # The split is consumed by the within-group operator algebra
    # ``A_wg = L + C - S.foldable_part()`` (Phase G substep 3+4.b)
    # and the within-group accelerator's iteration operator
    # ``scatter_cycle = A_wg.inverse() @ S.residual_part()``
    # (substep 3+4.c). This substep lands the data API only; no
    # solver/sweep/iteration code consumes the new methods yet.

    def foldable_part(self) -> "ScatteringOperator":
        r"""Return the P0 within-group self-scatter sibling of :math:`S`.

        Carries only the diagonal of ``sig_s[mid][0]`` per material —
        the within-group self-scatter cross-section
        :math:`\Sigma_{s,0}^{g\to g}` per energy group. All other
        scattering channels (cross-group P0, all :math:`P_\ell \ge 1`,
        and (n,2n)) live in :meth:`residual_part`.

        The foldability criterion is :math:`\ell = 0` within-group
        only because :math:`Y_\ell^m(\Omega_n)` makes
        :math:`P_\ell \ge 1` self-scatter inherently
        direction-dependent and unfoldable into a per-cell
        :math:`\sigma_r` (the cell-balance removal cross-section).

        Returns
        -------
        ScatteringOperator
            A sibling with ``scattering_order = 0``, ``Y = None``,
            diagonal-only ``sig_s[mid][0]``, zero ``sig2``, and the
            same ``n_ordinates / nx / ny / ng / weights /
            cells_by_mat`` as ``self``. Pure function — calling
            twice returns instances with equal per-material arrays.
        """
        # New per-material arrays — never mutate self.sig_s / self.sig2.
        sig_s_foldable: dict[int, list[np.ndarray]] = {}
        sig_s0_foldable: dict[int, np.ndarray] = {}
        sig2_foldable: dict[int, np.ndarray] = {}
        for mid, mats in self.sig_s.items():
            diag_only = np.diag(np.diag(mats[0]))  # (ng, ng), diagonal only
            sig_s_foldable[mid] = [diag_only]
            sig_s0_foldable[mid] = diag_only
            sig2_foldable[mid] = np.zeros_like(self.sig2[mid])
        return ScatteringOperator(
            n_ordinates=self.n_ordinates,
            nx=self.nx,
            ny=self.ny,
            ng=self.ng,
            scattering_order=0,
            sig_s=sig_s_foldable,
            sig2=sig2_foldable,
            sig_s0=sig_s0_foldable,
            Y=None,
            weights=self.weights,
            cells_by_mat=self.cells_by_mat,
        )

    def residual_part(self) -> "ScatteringOperator":
        r"""Return the non-foldable sibling of :math:`S`.

        Carries everything :meth:`foldable_part` does not: the
        off-diagonal of ``sig_s[mid][0]`` (cross-group P0), every
        :math:`P_\ell \ge 1` block verbatim, and ``sig2[mid]``
        verbatim. These channels cannot collapse into a per-cell
        :math:`\sigma_r`: cross-group P0 couples distinct energy
        groups, :math:`P_\ell \ge 1` is direction-dependent via
        :math:`Y_\ell^m(\Omega_n)`, and (n,2n) doubling emits two
        neutrons per absorption (folding into a "removal"
        cross-section is conceptually wrong).

        Returns
        -------
        ScatteringOperator
            A sibling with ``scattering_order ==
            self.scattering_order``, ``Y is self.Y`` (precomputed
            harmonics are reusable), zero-diagonal P0 matrix, and
            unchanged :math:`P_\ell \ge 1` + (n,2n) data.

        Notes
        -----
        Algebraic contract:
        ``S.apply(\psi) \approx S.foldable_part().apply(\psi) +
        S.residual_part().apply(\psi)`` at ``rtol=1e-14``. The
        residual is FP-non-associativity (sum-of-two-applies vs
        single-apply differs at machine precision).
        """
        sig_s_residual: dict[int, list[np.ndarray]] = {}
        sig_s0_residual: dict[int, np.ndarray] = {}
        for mid, mats in self.sig_s.items():
            p0 = mats[0]
            cross_group = p0 - np.diag(np.diag(p0))
            # Pℓ ≥ 1 blocks carried verbatim — anisotropic scattering is
            # unconditionally residual (Y_ℓ^m direction dependence).
            sig_s_residual[mid] = [cross_group, *mats[1:]]
            sig_s0_residual[mid] = cross_group
        # (n,2n) carried verbatim — same dict; safe to share since
        # apply() never mutates the matrices.
        sig2_residual = self.sig2
        return ScatteringOperator(
            n_ordinates=self.n_ordinates,
            nx=self.nx,
            ny=self.ny,
            ng=self.ng,
            scattering_order=self.scattering_order,
            sig_s=sig_s_residual,
            sig2=sig2_residual,
            sig_s0=sig_s0_residual,
            Y=self.Y,
            weights=self.weights,
            cells_by_mat=self.cells_by_mat,
        )

    def is_foldable_into_sigma_r(self) -> bool:
        r"""Return ``True`` iff this operator is structurally the
        :meth:`foldable_part` of some parent :math:`S`.

        Mechanical predicate consumed by substep 3+4.b.ii's
        :class:`~orpheus.numerics.operator.OperatorSum` fusion hook:
        when the hook sees ``L + C - S_some`` and
        ``S_some.is_foldable_into_sigma_r()`` is True, it knows the
        within-group sum has the shape that admits fusion into a
        removal cross-section :math:`\sigma_r = \sigma_t -
        \Sigma_{s,0}^{g\to g}`, and routes :meth:`solve` through the
        within-group sweep.

        Structural test on the operator's data, not an identity claim
        about its action — every ``ScatteringOperator`` instance whose
        ``scattering_order`` is 0 with diagonal-only ``sig_s[mid][0]``
        and zero ``sig2`` is, by definition, the foldable part of
        itself (its :meth:`residual_part` is the zero operator).

        Returns
        -------
        bool
            ``True`` iff (a) ``scattering_order == 0``,
            (b) every material's ``sig_s[mid][0]`` is diagonal
            (off-diagonal ≈ zero by ``np.allclose``), AND
            (c) every material's ``sig2[mid]`` is the zero matrix
            (``np.allclose`` to 0). ``False`` otherwise.

        Notes
        -----
        Uses ``np.allclose`` (default ``rtol=1e-5, atol=1e-8``) for the
        diagonal / zero checks: floating-point construction
        (e.g. ``np.diag(np.diag(M))``) introduces FP rounding so
        bit-equality is fragile. The tolerance is permissive enough
        to accept any physical input but strict enough to reject
        genuine off-diagonal entries.

        Contract verified by ``TestIsFoldableIntoSigmaR`` in
        :file:`tests/sn/test_scattering_operator.py`:

        * Full ``ScatteringOperator`` (non-zero off-diagonal P0 OR
          Pℓ ≥ 1 OR non-zero sig2) → ``False``.
        * ``S.foldable_part()`` → ``True`` (round-trip).
        * ``S.residual_part()`` → ``False`` (cross-group P0).
        * scattering_order=0 with non-diagonal P0 → ``False``.
        * scattering_order=0 with diagonal P0 but non-zero sig2 →
          ``False``.
        """
        if self.scattering_order != 0:
            return False
        for mid, mats in self.sig_s.items():
            p0 = mats[0]
            # Diagonal-only test: M ≈ diag(diag(M)).
            if not np.allclose(p0, np.diag(np.diag(p0))):
                return False
            # sig2 must be the zero matrix — (n,2n) is unconditionally
            # residual (see ``residual_part`` docstring).
            if not np.allclose(self.sig2[mid], 0.0):
                return False
        return True

    def foldable_sigma(self) -> dict[int, np.ndarray]:
        r"""Return the per-material foldable cross-section :math:`(\sigma_{s,0}^{g\to g})_g`.

        For every material ``mid``, returns the ``(ng,)`` array
        ``np.diag(sig_s[mid][0])`` — the per-group within-group
        self-scatter cross-section :math:`\Sigma_{s,0}^{g\to g}`
        that the cell-balance denominator absorbs into
        :math:`\sigma_r = \sigma_t - \Sigma_{s,0}^{g\to g}`.

        Substep 3+4.b's ``OperatorSum.solve`` fusion hook consumes
        these arrays directly (lazy :math:`\sigma_r` cache); the
        full :meth:`foldable_part` operator is consumed by the
        within-group algebra ``A_wg = L + C - S.foldable_part()``.
        Both are equivalent re-expressions of the same data.

        Returns
        -------
        dict[int, np.ndarray]
            ``{mid: (ng,) array}`` keyed by material id. Each array
            is a fresh copy (mutating the returned dict's values
            does not affect ``self``).
        """
        return {
            mid: np.diag(mats[0]).copy()
            for mid, mats in self.sig_s.items()
        }

    # ── LinearOperator surface ─────────────────────────────────────────

    def apply(self, psi: np.ndarray) -> np.ndarray:
        r"""Apply the full scattering source operator :math:`S\,\psi`.

        Returns a per-ordinate isotropic-equivalent source plus the
        :math:`\ell \ge 1` per-ordinate contribution, packaged together
        as the :math:`(N, n_g, n_x, n_y)` array the operator-algebra
        consumers expect under PR-INDEX-4.

        For the P0 + (n,2n) parts, the action is the same scalar source
        broadcast across every ordinate (the sweep's ``1/W`` factor is
        not applied here — the caller, ``transport_sweep``, applies it).
        For :math:`\ell \ge 1`, the action carries genuine per-ordinate
        directional content.

        Parameters
        ----------
        psi : np.ndarray
            Angular flux shape ``(N, ng, nx, ny)`` (Issue #196 PR-INDEX-4 —
            principled).

        Returns
        -------
        np.ndarray
            Per-ordinate scattering source shape
            ``(N, ng, nx, ny)``: the isotropic P0 + (n,2n) part is
            broadcast across the ``N`` axis; the :math:`\ell \ge 1`
            part is the genuine per-ordinate Galerkin reconstruction.

        Notes
        -----
        ORPHEUS's existing :class:`SNSolver` source iteration consumes
        the isotropic and anisotropic pieces separately (the isotropic
        :math:`Q` is a ``(ng, nx, ny)`` scalar source; the anisotropic
        :math:`Q_{\rm aniso}` is a ``(N, ng, nx, ny)`` per-ordinate
        source).  For internal use by :class:`SNSolver` the helpers
        :meth:`add_iso_source`, :meth:`add_n2n_source`, and
        :meth:`build_aniso_source` expose the same math without forcing
        a wasteful broadcast. :meth:`apply` is the LinearOperator
        Protocol surface and is intended for downstream operator-algebra
        consumers (Krylov-on-apply, adjoint sweeps, composition with
        :math:`L` and :math:`F`).
        """
        # Reduce angular flux to scalar flux: φ_g(r) = Σ_n w_n ψ_n,g(r)
        # — leading-axis contraction over ordinates.  Output is
        # principled ``(ng, nx, ny)``.
        phi = np.einsum('n,ngxy->gxy', self.weights, psi)

        # Isotropic (P0 + (n,2n)) — vectorised by material into a fresh
        # (ng, nx, ny) buffer, then broadcast across the N axis.
        Q_iso = np.zeros((self.ng, self.nx, self.ny))
        self.add_iso_source(Q_iso, phi)
        self.add_n2n_source(Q_iso, phi)

        # Broadcast across ordinates and add the Pℓ (l>=1) contribution.
        # Q_iso (ng, nx, ny) → (1, ng, nx, ny) → broadcast to psi.shape
        # (N, ng, nx, ny).  Same numpy idiom as the legacy
        # (1, nx, ny, ng) broadcast — just with the axis order updated.
        Q = np.broadcast_to(Q_iso[None, :, :, :], psi.shape).copy()
        Q_aniso = self.build_aniso_source(psi)
        if Q_aniso is not None:
            Q += Q_aniso
        return Q
