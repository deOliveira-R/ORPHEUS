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
from functools import singledispatchmethod
from typing import TYPE_CHECKING

import numpy as np

from orpheus.numerics.operator import (
    CAP_APPLY,
    LinearOperatorMixin,
)
from orpheus.numerics.projection import (
    MomentProjection,
    HarmonicMomentReconstruction,
)

# Runtime imports of the flux / source types — required at module load
# time because :func:`singledispatchmethod.register` dispatches on the
# runtime class.  These three modules form a leaf in the SN package
# dependency graph (they do not import scattering.py), so the imports
# are circular-import-safe.
from orpheus.transport.fields.scalar_flux import ScalarFlux
from orpheus.transport.fields.angular_flux import AngularFlux
from orpheus.transport.fields.boundary_flux import BoundaryFlux
from orpheus.transport.sources import ScalarSource, AngularSource
from orpheus.transport.timed_full_field import TimedFullField

if TYPE_CHECKING:
    from orpheus.transport.fields.harmonic_moment_field import HarmonicMomentField
    from .material_xs_field import MaterialXSField
    from orpheus.numerics.quadrature import Quadrature


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
    mat_xs : MaterialXSField
        Macroscopic XS field carrying both the per-material Legendre
        scattering matrices and the cell-to-material map.  Issue #197
        PR-TYPED-1 — the per-material dispatch loop lives ONLY inside
        :meth:`MaterialXSField.apply_legendre_scattering_moments`.
    L : int
        Maximum Legendre order :math:`L` retained.
    skip_l0 : bool, default ``True``
        Skip the :math:`\ell = 0` block (handled separately by P0
        in-scatter). Set to ``False`` for the full LinearOperator
        composition :math:`R \Lambda M \psi`.
    """

    mat_xs: "MaterialXSField"
    L: int
    skip_l0: bool = True
    capabilities: frozenset = field(
        default_factory=lambda: frozenset({CAP_APPLY})
    )

    def apply(
        self, moments: "np.ndarray | HarmonicMomentField",
    ) -> "np.ndarray | HarmonicMomentField":
        r"""Apply :math:`\Lambda` to a moment field.

        Parameters
        ----------
        moments : np.ndarray or HarmonicMomentField
            Moment field of shape ``(L+1, 2L+1, ng, nx, ny)`` (Issue
            #196 PR-INDEX-4 — principled).  The :math:`m`-axis is the
            addition-theorem-shifted index where slot ``l + m`` holds
            the :math:`(\ell, m)` entry; entries outside
            :math:`|m| \le \ell` are conventionally zero.

            Issue #197 PR-TYPED-4 — typed
            :class:`~orpheus.sn.harmonic_moment_field.HarmonicMomentField`
            is accepted; when supplied, the return is a typed field
            with matching ``L`` and ``mesh``.  Bare ndarray is still
            accepted for legacy probe / test callers.

        Returns
        -------
        np.ndarray or HarmonicMomentField
            Same shape as ``moments``.  The :math:`\ell = 0` block is
            zero when ``skip_l0`` is ``True``; otherwise the P0
            in-scatter contribution is included.  Type matches the
            input.

        Notes
        -----
        Per-material dispatch is encapsulated by
        :meth:`MaterialXSField.apply_legendre_scattering_moments` —
        Issue #197 PR-TYPED-1 collapses the ``for mid, (ix, iy) in
        cells_by_mat.items()`` loop into a typed verb.
        """
        from orpheus.transport.fields.harmonic_moment_field import HarmonicMomentField
        if isinstance(moments, HarmonicMomentField):
            out_values = self.mat_xs.apply_legendre_scattering_moments(
                moments.values, L=self.L, skip_l0=self.skip_l0,
            )
            return HarmonicMomentField.from_mesh_and_L(out_values, moments.mesh, moments.L)
        return self.mat_xs.apply_legendre_scattering_moments(
            moments, L=self.L, skip_l0=self.skip_l0,
        )


@dataclass(frozen=True)
class _PerLegendreOrderScattering(LinearOperatorMixin):
    r"""Wave T step T.3 — per-:math:`\ell` summand of the anisotropic
    in-scatter source operator.

    Implements one term of the §15.2 sum-of-products form:

    .. math::

        S_\ell\,\psi \;=\; R_\ell\,\Lambda_\ell\,M_\ell\,\psi

    where

    * :math:`M_\ell` projects the angular flux onto the :math:`\ell`-th
      harmonic block (a single ``(2\ell+1, n_g, n_x, n_y)`` slice of
      the full moment tensor);
    * :math:`\Lambda_\ell` applies the per-material per-:math:`\ell`
      scattering matrix :math:`\Sigma_{s,\ell}^{m(\vec r)}[g'\to g]`
      to that block;
    * :math:`R_\ell` reconstructs the per-ordinate output from the
      :math:`\ell`-th moment block back to ordinate space.

    The full anisotropic in-scatter source is the sum across
    :math:`\ell \in [1, L]`:

    .. math::

        S_{\text{aniso}}\,\psi
        \;=\; \frac{1}{W}\,\sum_{\ell=1}^{L}\,R_\ell\,\Lambda_\ell\,M_\ell\,\psi

    where the producer-side :math:`1/W` normalisation (Pattern 7) is
    applied OUTSIDE the kernel at the
    :meth:`~orpheus.sn.scattering.ScatteringOperator.apply` boundary.

    Wave T design context
    ---------------------

    Plan §6 T.3 originally targeted a
    :class:`~orpheus.numerics.operator.SumOfTensorProductsOperator`
    per Grand Report v3 §15.2.  However, the per-material per-:math:`\ell`
    einsum at :meth:`MaterialXSField.apply_legendre_scattering_moments`
    couples the group axis (matrix multiply on
    :math:`\Sigma_{s,\ell}[g'\to g]`) with the spatial axis (via
    :attr:`cells_by_material` indexing).  That coupling violates the
    :class:`~orpheus.numerics.operator.TensorProductOperator`
    disjoint-axes contract — no factor design satisfies "each factor
    acts on one tensor axis and broadcasts on the rest" without one
    factor doing all the work and the other(s) being degenerate
    identity.

    Per the test-architect spec §6 Q6 (math-honest fallback) the
    kernel is an :class:`~orpheus.numerics.operator.OperatorSum` of
    per-:math:`\ell` summands of this class — NOT a
    :class:`~orpheus.numerics.operator.SumOfTensorProductsOperator`.
    The §15.2 form is preserved at the summation level (one summand
    per Legendre order); the per-summand decomposition into
    :math:`R_\ell \circ \Lambda_\ell \circ M_\ell` is a procedural
    composition, not a tensor product.

    Capability set: :pydata:`{CAP_APPLY}` only.  Per-:math:`\ell`
    scattering is non-invertible (rank-deficient on most ordinate
    inputs) and the adjoint surface is deferred to a future wave.

    Parameters
    ----------
    mat_xs : MaterialXSField
        The macroscopic XS field.  Read-through for
        :attr:`emission_spectrum`, :attr:`fission_production`,
        :meth:`sig_s_legendre`, :attr:`cells_by_material`.
    ell : int
        The single Legendre order this summand handles.  Must satisfy
        ``1 <= ell <= L`` (the :math:`\ell = 0` block is in the P0 +
        ``(n,2n)`` fast path, NOT the SOTP-style kernel).
    L : int
        The total scattering order — needed to size the moment tensor
        ``(L+1, 2L+1, ng, nx, ny)`` that :math:`M` projects onto.
    weights : np.ndarray
        ``(N,)`` ordinate weights — needed for :math:`M`.
    Y : np.ndarray
        ``(N, L+1, 2L+1)`` real spherical harmonics — needed for both
        :math:`M` and :math:`R`.

    Notes
    -----
    The apply operates on the bare ``(N, ng, nx, ny)`` per-ordinate
    angular flux ndarray and returns a same-shape per-ordinate
    contribution.  The :math:`1/W` projection is OUTSIDE the kernel
    (Pattern 7 producer-side normalisation at the apply boundary).
    """

    mat_xs: "MaterialXSField"
    ell: int
    L: int
    weights: np.ndarray
    Y: np.ndarray
    capabilities: frozenset = field(
        default_factory=lambda: frozenset({CAP_APPLY})
    )

    def apply(self, psi_values: np.ndarray) -> np.ndarray:
        r"""Apply :math:`R_\ell \Lambda_\ell M_\ell` to per-ordinate flux.

        Parameters
        ----------
        psi_values : np.ndarray
            ``(N, ng, nx, ny)`` per-ordinate angular flux.

        Returns
        -------
        np.ndarray
            ``(N, ng, nx, ny)`` per-ordinate in-scatter contribution
            from the single Legendre order ``self.ell``.  The
            :math:`1/W` projection is NOT applied here — that's the
            kernel-consumer's responsibility at the apply boundary.
        """
        # M: project ψ onto moment tensor (L+1, 2L+1, ng, nx, ny).  The
        # full moment tensor is built; only the ell-th block is used by
        # this summand.
        M_op = MomentProjection(weights=self.weights, Y=self.Y, L=self.L)
        moments_full = M_op.apply(psi_values)

        # Λ_ℓ: per-material per-ℓ scattering on ONLY the ell-th block.
        # Build a moment tensor with all blocks zero except ℓ=self.ell.
        out_moments = np.zeros_like(moments_full)
        n_m = 2 * self.ell + 1
        for mid, (ix, iy) in self.mat_xs.cells_by_material.items():
            sig_s_mid = self.mat_xs.sig_s_legendre(mid)
            moments_view = moments_full[self.ell, :n_m][..., ix, iy]
            # Trailing-contiguous einsum (same pattern as the full
            # MaterialXSField.apply_legendre_scattering_moments
            # body).  Idempotent additivity isn't needed here — the
            # output starts at zero and only the ell-th block is touched.
            out_moments[self.ell, :n_m][..., ix, iy] = np.einsum(
                "mfc,fg->mgc",
                moments_view,
                sig_s_mid[self.ell],
            )

        # R: reconstruct ordinate space from the moment tensor.  Only
        # the ell-th block is non-zero, so the output captures ONLY the
        # per-ell in-scatter contribution.
        R_op = HarmonicMomentReconstruction.from_Y(self.Y)
        return R_op.apply(out_moments)


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
    mat_xs : MaterialXSField
        Macroscopic XS field — the single source of truth for both
        per-material scattering / (n,2n) data AND the cell-to-material
        topology (Issue #197 PR-TYPED-1).  Every per-material loop
        that previously lived on this class now routes through one of
        the typed verbs (``mat_xs.apply_p0_in_scatter``, ``apply_n2n``,
        etc.) which encapsulates the dispatch.
    quadrature : AngularQuadrature
        The angular quadrature.  Carries ``N``, ``weights``, and
        :meth:`spherical_harmonics` — the previously-leaked
        ``n_ordinates`` / ``weights`` / ``Y`` constructor parameters
        all collapse into ``self.quadrature``.
    scattering_order : int
        Maximum Legendre order :math:`L` retained. ``0`` means P0 only.
    capabilities : frozenset[str]
        ``{"apply"}`` — :math:`S` has no efficient inverse and the
        current ORPHEUS algebra has no consumer for :math:`S^T`.
    """

    mat_xs: "MaterialXSField"
    quadrature: "AngularQuadrature"
    scattering_order: int

    capabilities: frozenset[str] = field(
        default_factory=lambda: frozenset({CAP_APPLY})
    )

    # Lazy cache for the precomputed spherical harmonics — only
    # populated when ``scattering_order > 0`` (avoids paying the
    # cost on P0-only problems).
    _Y_cached: np.ndarray | None = field(
        default=None, init=False, repr=False,
    )

    # ── Convenience read-throughs (the legacy n_ordinates / nx / ny
    # / ng / weights / Y / cells_by_mat / sig_s / sig2 / sig_s0
    # attributes the test suite + internal call sites consumed).
    # Issue #197 PR-TYPED-1: these become read-throughs onto
    # mat_xs + quadrature so consumers don't break.  Marked TRANSIENT
    # in source comments where appropriate; full retirement deferred
    # to a follow-up PR once consumers all read mat_xs.* directly.

    @property
    def n_ordinates(self) -> int:
        """Number of angular ordinates :math:`N`."""
        return self.quadrature.N

    @property
    def nx(self) -> int:
        """Spatial extent in x — read-through from :attr:`mat_xs`."""
        return self.mat_xs.nx

    @property
    def ny(self) -> int:
        """Spatial extent in y."""
        return self.mat_xs.ny

    @property
    def ng(self) -> int:
        """Energy group count."""
        return self.mat_xs.ng

    @property
    def weights(self) -> np.ndarray:
        """Quadrature weights ``(N,)``."""
        return self.quadrature.weights

    @property
    def Y(self) -> np.ndarray | None:
        r"""Real spherical harmonics :math:`Y_\ell^m(\Omega_n)`,
        shape ``(N, L+1, 2L+1)``, or ``None`` when ``L == 0``."""
        if self.scattering_order == 0:
            return None
        if self._Y_cached is None:
            self._Y_cached = self.quadrature.spherical_harmonics(
                self.scattering_order,
            )
        return self._Y_cached

    @property
    def sig_s(self) -> dict[int, list[np.ndarray]]:
        """TRANSIENT — per-material dense Legendre scattering dict.

        Routes through :meth:`MaterialXSField.sig_s_legendre`.  Kept
        as a shim for the four ``_build_rhs_*`` helpers in
        :mod:`orpheus.sn.solver` that still consume the per-material
        dict directly.  PR-TYPED-2 will rewire those helpers to
        consume ``mat_xs`` directly and retire this shim.
        """
        return {
            mid: self.mat_xs.sig_s_legendre(mid)
            for mid in self.mat_xs.materials
        }

    @property
    def sig2(self) -> dict[int, np.ndarray]:
        """TRANSIENT — per-material dense (n,2n) dict.  See :attr:`sig_s`."""
        return {
            mid: self.mat_xs.n2n_matrix(mid)
            for mid in self.mat_xs.materials
        }

    @property
    def sig_s0(self) -> dict[int, np.ndarray]:
        """TRANSIENT — per-material P0 scattering matrix dict.
        See :attr:`sig_s`."""
        return {
            mid: self.mat_xs.sig_s_legendre(mid)[0]
            for mid in self.mat_xs.materials
        }

    @property
    def cells_by_mat(self) -> dict[int, tuple[np.ndarray, np.ndarray]]:
        """TRANSIENT — per-material cell-index dict.  See :attr:`sig_s`."""
        return self.mat_xs.cells_by_material

    # ── Wave T step T.3 — per-ℓ kernel (anisotropic in-scatter) ────────

    @property
    def kernel_summands(
        self,
    ) -> tuple["_PerLegendreOrderScattering", ...]:
        r"""Wave T step T.3 — per-:math:`\ell` summands of the
        anisotropic in-scatter kernel.

        Returns one :class:`_PerLegendreOrderScattering` instance per
        Legendre order :math:`\ell \in [1, L]`.  The summation
        ``sum(s.apply(psi.values) for s in kernel_summands)`` reproduces
        :meth:`build_aniso_source` (pre :math:`/W`) at FP-non-associativity
        precision.

        Empty when ``self.scattering_order == 0`` (no anisotropic
        in-scatter; the P0 + (n,2n) fast path handles all of S).  The
        :math:`\ell = 0` block is NOT included in this list per the
        test-architect spec §6 Q3 Option (β): P0 stays in the
        :meth:`add_iso_source` + :meth:`add_n2n_source` fast path,
        architecturally distinct from the moment-tensor path.
        """
        if self.scattering_order == 0:
            return ()
        return tuple(
            _PerLegendreOrderScattering(
                mat_xs=self.mat_xs,
                ell=ell,
                L=self.scattering_order,
                weights=self.weights,
                Y=self.Y,  # type: ignore[arg-type]
            )
            for ell in range(1, self.scattering_order + 1)
        )

    @property
    def kernel(self) -> "LinearOperatorMixin":
        r"""Wave T step T.3 — the anisotropic-in-scatter kernel.

        Returns the :class:`~orpheus.numerics.operator.OperatorSum`
        tree of the per-:math:`\ell` summands at
        :attr:`kernel_summands`, OR a
        :class:`~orpheus.numerics.operator.ZeroOperator` when
        ``scattering_order == 0``.

        Math-honest fallback (test-architect spec §6 Q6).  The §15.2
        target form was a
        :class:`~orpheus.numerics.operator.SumOfTensorProductsOperator`
        but the per-material per-:math:`\ell` einsum couples group +
        spatial axes — the disjoint-axes contract fails.  Per the
        spec's Q6 escalation path, the kernel is an
        :class:`OperatorSum` of custom per-:math:`\ell` summands
        (NOT a SOTP), preserving the §15.2 form at the summation
        level (one summand per Legendre order) while staying
        algebraically honest about the per-summand composition
        (:math:`R_\ell \circ \Lambda_\ell \circ M_\ell`, not a tensor
        product).

        The :math:`1/W` Pattern 7 producer-side normalisation lives
        OUTSIDE this kernel at the
        :meth:`~ScatteringOperator.apply` boundary — not inside any
        summand.

        Substep T.3b ships this property without rewiring any
        :meth:`apply` arm; the call-site migration lands in T.3c
        (`build_aniso_source` rewire) and T.3d (per-arm migration).

        Built fresh on each access to honor the read-through semantics
        for :attr:`mat_xs` (depletion / thermal feedback may update
        cross-sections in-place).
        """
        from functools import reduce
        from operator import add as op_add
        from orpheus.numerics.operator import ZeroOperator

        summands = self.kernel_summands
        if not summands:
            return ZeroOperator()
        return reduce(op_add, summands)

    @classmethod
    def from_solver_data(
        cls,
        *,
        mat_xs: "MaterialXSField",
        quadrature: "AngularQuadrature",
        scattering_order: int,
    ) -> "ScatteringOperator":
        """Construct from a :class:`MaterialXSField` + quadrature.

        Issue #197 PR-TYPED-1 — the constructor surface collapses the
        eight separate per-material handles (``sig_s``, ``sig2``,
        ``sig_s0``, ``cells_by_mat``, ``Y``, ``weights``, ``n_ordinates``,
        ``nx``, ``ny``, ``ng``) into two typed objects.  The
        :class:`MaterialXSField` carries everything per-material plus
        the spatial topology; the :class:`AngularQuadrature` carries
        ``N`` / ``weights`` / harmonics.
        """
        return cls(
            mat_xs=mat_xs,
            quadrature=quadrature,
            scattering_order=scattering_order,
        )

    # ── In-place helpers (preserve bit-identity vs SNSolver pre-Wave-D) ─

    def add_iso_source(
        self,
        Q: "np.ndarray | ScalarSource",
        phi: "np.ndarray | ScalarFlux",
    ) -> "np.ndarray | ScalarSource | None":
        r"""Add P0 in-scatter source to :math:`Q`.

        Vectorised by material: per cell ``c`` of material ``mid``,
        ``Q[:, ic, jc] += sig_s0[mid].T @ phi[:, ic, jc]`` where
        ``sig_s0[mid]`` is ``(ng, ng)`` indexed ``[g_from, g_to]``.

        Issue #197 PR-TYPED-3 introduces the typed-action overload:

        * **Raw-in (legacy)** — ``Q: np.ndarray`` is mutated in place
          and ``None`` is returned (the Wave A–D / PR-TYPED-1 contract).
        * **Typed-in (return-new)** — ``Q: ScalarSource`` is treated
          as an immutable input; a NEW :class:`ScalarSource` is
          returned carrying ``Q.values + in_scatter`` (Pattern 4 —
          frozen typed inputs stay immutable; Pattern 1 — the caller
          spells the algebra as ``Q = scattering.add_iso_source(Q,
          phi)``).

        Parameters
        ----------
        Q : np.ndarray or ScalarSource
            Isotropic source carrier.  Raw ``(ng, nx, ny)`` ndarray is
            mutated in place; typed :class:`ScalarSource` returns a
            new instance.
        phi : np.ndarray or ScalarFlux
            Scalar flux.  Either form is accepted; the underlying
            values are unwrapped before the per-material dispatch.

        Returns
        -------
        np.ndarray or ScalarSource or None
            Raw-in returns ``None`` (legacy in-place); typed-in returns
            a fresh :class:`ScalarSource`.

        Notes
        -----
        Issue #197 PR-TYPED-1 — the per-material dispatch lives ONLY
        inside :meth:`MaterialXSField.apply_p0_in_scatter`.
        """
        from orpheus.transport.fields.scalar_flux import ScalarFlux
        from orpheus.transport.sources import ScalarSource
        phi_values = phi.values if isinstance(phi, ScalarFlux) else phi
        if isinstance(Q, ScalarSource):
            Q_values = Q.values.copy()
            self.mat_xs.apply_p0_in_scatter(Q_values, phi_values)
            return ScalarSource.from_mesh(Q_values, Q.mesh)
        self.mat_xs.apply_p0_in_scatter(Q, phi_values)
        return None

    def add_n2n_source(
        self,
        Q: "np.ndarray | ScalarSource",
        phi: "np.ndarray | ScalarFlux",
    ) -> "np.ndarray | ScalarSource | None":
        r"""Add (n,2n) source to :math:`Q`.

        Vectorised by material: per cell ``c`` of material ``mid``,
        ``Q[:, ic, jc] += 2 · sig2[mid].T @ phi[:, ic, jc]``.

        Issue #197 PR-TYPED-3 introduces the same typed-action overload
        as :meth:`add_iso_source` — raw-in mutates in place and returns
        ``None``; typed-in returns a fresh :class:`ScalarSource`.

        Parameters
        ----------
        Q : np.ndarray or ScalarSource
            Isotropic source carrier.
        phi : np.ndarray or ScalarFlux
            Scalar flux.

        Returns
        -------
        np.ndarray or ScalarSource or None
            Raw-in returns ``None`` (legacy in-place); typed-in returns
            a fresh :class:`ScalarSource`.

        Notes
        -----
        Issue #197 PR-TYPED-1 — per-material dispatch lives ONLY
        inside :meth:`MaterialXSField.apply_n2n`.
        """
        from orpheus.transport.fields.scalar_flux import ScalarFlux
        from orpheus.transport.sources import ScalarSource
        phi_values = phi.values if isinstance(phi, ScalarFlux) else phi
        if isinstance(Q, ScalarSource):
            Q_values = Q.values.copy()
            self.mat_xs.apply_n2n(Q_values, phi_values)
            return ScalarSource.from_mesh(Q_values, Q.mesh)
        self.mat_xs.apply_n2n(Q, phi_values)
        return None

    def build_aniso_source(
        self,
        angular_flux: "np.ndarray | AngularFlux | None",
    ) -> "np.ndarray | AngularSource | None":
        r"""Build per-ordinate Pℓ (:math:`\ell \ge 1`) scattering source.

        Implements the Galerkin reconstruction :eq:`pn-scatter` from
        the angular-flux moments :eq:`flux-moments` (declared in
        :ref:`theory-discrete-ordinates`) as the **literal operator
        composition**

        .. math::

            Q^{\rm aniso}_n(\vec r) \;=\;
            \frac{1}{W}\,(R \, \Lambda \, M \, \psi)_n(\vec r)

        of the three primitives shipped in
        :mod:`orpheus.numerics.projection` and
        :class:`LegendreMomentScattering` (Wave 1 / C1.1):

        * :math:`M` :class:`~orpheus.numerics.projection.MomentProjection`:
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
        scaling sourced from
        :attr:`~orpheus.numerics.spaces.SphericalHarmonicSpace.addition_theorem_factor`
        — the single canonical home of the SH convention's
        :math:`(2\ell+1)` literal (ERR-039 endpoint).  :math:`R` is the
        addition-theorem reconstruction, NOT the W-weighted Hilbert
        adjoint of :math:`M` (the two differ by the
        :class:`SphericalHarmonicSpace` Gram :math:`g_C`); the
        relationship is pinned by
        ``tests/numerics/test_spherical_harmonic_space.py``.

        The trailing :math:`1/W = 1/\sum_n w_n` is the **producer-side
        per-ordinate projection** introduced in R-1 Step 4 A1: the Pℓ
        aniso source enters the per-ordinate transport equation
        :math:`(\Omega\cdot\nabla + \sigma_t)\psi_n = Q/W +
        q_{\rm aniso,n}` already in per-ordinate magnitude, so the
        sweep does NOT need to apply ``/W`` again
        (cf. :func:`~orpheus.sn.sweep.transport_sweep` whose
        ``aniso_source`` parameter is documented as per-ordinate
        magnitude post-A1).  See ``coding-elegance`` SKILL.md
        §"Convention crosswalk template" / lesson L18.

        Total flop count is identical to the legacy hand-rolled
        ``for n in range(N)`` triple loop; the iteration over ordinates
        is now internal to :func:`numpy.einsum` inside :math:`R` and
        :math:`M`, **not** a Python loop.

        Parameters
        ----------
        angular_flux : np.ndarray or AngularFlux or None
            Angular flux shape ``(N, ng, nx, ny)`` (Issue #196
            PR-INDEX-4 — principled) from the most recent sweep, or
            ``None`` on the first source iteration before any sweep
            has been run.  Issue #197 PR-TYPED-3 — typed
            :class:`AngularFlux` is accepted; the result is then a
            typed :class:`AngularSource` (preserving the type chain
            through the scattering composition).

        Returns
        -------
        np.ndarray or AngularSource or None
            ``(N, ng, nx, ny)`` per-ordinate :math:`\ell \ge 1`
            contribution in **per-ordinate magnitude** (the trailing
            ``/W`` is applied here, R-1 Step 4 A1).  Returns ``None``
            when ``scattering_order == 0`` or ``angular_flux is None``.
            Type matches the input.

        Notes
        -----
        :class:`MomentProjection` and
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
        # D-I.2 (2026-05-29): only the typed :class:`AngularFlux` is
        # accepted.  The bare-ndarray fallthrough retired alongside the
        # singledispatch ``np.ndarray`` arm on :meth:`apply`; the typed
        # leaf carries the mesh, which is required to construct the
        # output :class:`AngularSource`.
        mesh = angular_flux.mesh
        # Wave T step T.3c — delegate the §9 "S = R Λ M" inner numerics
        # to the typed :attr:`kernel` (an :class:`OperatorSum` of per-ℓ
        # :class:`_PerLegendreOrderScattering` summands per the §15.2
        # form, Q6 honest-fallback resolution).  The kernel's apply
        # produces the pre-:math:`/W` per-ordinate output; the
        # :math:`/sum_w` Pattern 7 producer-side normalisation lives
        # HERE at the apply boundary, OUTSIDE the kernel (per spec
        # §6 Q5).
        sum_w = float(self.weights.sum())
        out_values = self.kernel.apply(angular_flux.values) / sum_w
        return AngularSource.from_mesh(out_values, mesh)

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

        Returns
        -------
        ScatteringOperator
            A sibling with ``scattering_order = 0``, diagonal-only
            ``sig_s[mid][0]``, zero ``sig2``, and a derived
            :class:`MaterialXSField` carrying the foldable overrides.

        Notes
        -----
        Issue #197 PR-TYPED-1 — the per-material loop that built the
        ``sig_s_foldable`` / ``sig2_foldable`` dicts lives ONLY inside
        :meth:`MaterialXSField.foldable_sig_s`.
        """
        sig_s_foldable = self.mat_xs.foldable_sig_s()
        sig2_foldable = {
            mid: np.zeros_like(self.mat_xs.n2n_matrix(mid))
            for mid in self.mat_xs.materials
        }
        derived_mat_xs = self.mat_xs.with_overridden_sig_s_and_n2n(
            sig_s_dense=sig_s_foldable,
            n2n_dense=sig2_foldable,
        )
        return ScatteringOperator(
            mat_xs=derived_mat_xs,
            quadrature=self.quadrature,
            scattering_order=0,
        )

    def residual_part(self) -> "ScatteringOperator":
        r"""Return the non-foldable sibling of :math:`S`.

        Carries everything :meth:`foldable_part` does not: the
        off-diagonal of ``sig_s[mid][0]`` (cross-group P0), every
        :math:`P_\ell \ge 1` block verbatim, and ``sig2[mid]``
        verbatim.

        Returns
        -------
        ScatteringOperator
            A sibling with ``scattering_order ==
            self.scattering_order``, zero-diagonal P0 matrix, and
            unchanged :math:`P_\ell \ge 1` + (n,2n) data.

        Notes
        -----
        Algebraic contract:
        ``S.apply(\psi) \approx S.foldable_part().apply(\psi) +
        S.residual_part().apply(\psi)`` at ``rtol=1e-14``.

        Issue #197 PR-TYPED-1 — the per-material loop that built the
        ``sig_s_residual`` dict lives ONLY inside
        :meth:`MaterialXSField.residual_sig_s`.
        """
        sig_s_residual = self.mat_xs.residual_sig_s()
        # (n,2n) carried verbatim — pull from mat_xs.
        sig2_residual = {
            mid: self.mat_xs.n2n_matrix(mid)
            for mid in self.mat_xs.materials
        }
        derived_mat_xs = self.mat_xs.with_overridden_sig_s_and_n2n(
            sig_s_dense=sig_s_residual,
            n2n_dense=sig2_residual,
        )
        return ScatteringOperator(
            mat_xs=derived_mat_xs,
            quadrature=self.quadrature,
            scattering_order=self.scattering_order,
        )

    def is_foldable_into_sigma_r(self) -> bool:
        r"""Return ``True`` iff this operator is structurally the
        :meth:`foldable_part` of some parent :math:`S`.

        Returns
        -------
        bool
            ``True`` iff ``scattering_order == 0`` AND every material's
            P0 is diagonal AND every material's (n,2n) is zero.

        Notes
        -----
        Issue #197 PR-TYPED-1 — the per-material allclose checks live
        ONLY inside :meth:`MaterialXSField.is_p0_diagonal_with_zero_n2n`.
        """
        if self.scattering_order != 0:
            return False
        return self.mat_xs.is_p0_diagonal_with_zero_n2n()

    def foldable_sigma(self) -> dict[int, np.ndarray]:
        r"""Return the per-material foldable cross-section :math:`(\sigma_{s,0}^{g\to g})_g`.

        Returns
        -------
        dict[int, np.ndarray]
            ``{mid: (ng,) array}`` keyed by material id. Each array
            is a fresh copy.

        Notes
        -----
        Issue #197 PR-TYPED-1 — delegates to
        :meth:`MaterialXSField.foldable_sigma`.
        """
        return self.mat_xs.foldable_sigma()

    # ── LinearOperator surface ─────────────────────────────────────────

    @singledispatchmethod
    def apply(self, psi):
        r"""Apply the full scattering source operator :math:`S\,\psi`.

        Dispatched on input type via
        :func:`functools.singledispatchmethod`:

        * :class:`~orpheus.sn.angular_flux.AngularFlux` (legacy) →
          :class:`AngularFlux` — full :math:`P_\ell` Galerkin
          reconstruction in **per-ordinate magnitude** (the trailing
          :math:`1/W` projection lives at this producer boundary;
          R-1 Step 4 A1).  Consumers are SN sweep / GMRES / source-
          iteration loops operating on legacy per-ordinate angular
          flux.  Retires in D-H.1c alongside the legacy ``AngularFlux``.
        * :class:`~orpheus.transport.timed_full_field.TimedFullField`
          → :class:`TimedFullField` — composite bulk + boundary
          variant for the D-H.1+ migration path.  Bulk follows the
          same full :math:`P_\ell` Galerkin path; boundary is the
          implicit-zero :class:`BoundaryFlux` (Option β3 —
          scattering is volumetric; Wave O Issue #208 will encode the
          bulk-only nature in the type via :class:`BulkOperator`).
        * :class:`~orpheus.sn.scalar_flux.ScalarFlux` →
          :class:`~orpheus.sn.sources.ScalarSource` — :math:`P_0`
          in-scatter + :math:`(n,2n)` doubling only, in **iso scalar
          magnitude**.  No :math:`P_\ell` (scalar flux has no angular
          info); no :math:`1/W` (scalar consumers — diffusion, CP,
          kinetics — do not project to per-ordinate).
        * :class:`numpy.ndarray` — legacy bare-ndarray path,
          ``(N, ng, nx, ny)`` per-ordinate magnitude (consistent with
          the typed :class:`AngularFlux` variant post-A1).  Preserved
          for FD-matvec / probe-test call sites that bypass the type
          layer.

        The single :class:`ScatteringOperator` instance now serves
        every consumer kind via type-directed dispatch — Pattern 1
        (read-as-the-math via dunder) + Pattern 7 (producer-side
        normalisation, the :math:`1/W` lives at the apply boundary,
        decided by consumer type).  See ``coding-elegance`` SKILL.md
        §"Convention crosswalk template" and lesson L18.

        Parameters
        ----------
        psi : AngularFlux or ScalarFlux or np.ndarray
            Flux argument.  Dispatch on runtime type.

        Returns
        -------
        AngularFlux or ScalarSource or np.ndarray
            Scattering source in the convention appropriate to the
            consumer:

            * AngularFlux input → AngularFlux output, per-ordinate
              magnitude, full P_ℓ.
            * ScalarFlux input → ScalarSource output, iso scalar
              magnitude, P_0 + (n,2n) only.
            * np.ndarray input → np.ndarray output ``(N, ng, nx, ny)``,
              per-ordinate magnitude.

        Notes
        -----
        Internal helpers :meth:`add_iso_source`, :meth:`add_n2n_source`,
        and :meth:`build_aniso_source` remain available for callers
        that need the iso / aniso pieces separately.
        :meth:`build_aniso_source` returns **per-ordinate magnitude**
        post-A1 (the :math:`1/W` is applied there); the per-ordinate
        :class:`AngularFlux` variant of :meth:`apply` combines the iso
        piece (in iso magnitude, projected by :math:`1/W` at the
        boundary) with the aniso piece (already per-ordinate) into a
        single :class:`AngularSource`.
        """
        raise TypeError(
            f"ScatteringOperator.apply: unsupported input type "
            f"{type(psi).__name__}; expected AngularFlux, ScalarFlux, "
            f"or numpy.ndarray.  Dispatch table is registered via "
            f"@singledispatchmethod."
        )

    @apply.register
    def _(self, psi: TimedFullField) -> "TimedFullField":
        r"""Composite :class:`TimedFullField` variant — bulk-only scattering source.

        Math: identical to the :class:`AngularFlux` branch above —
        reduce ``psi.bulk`` angular → scalar, build iso :math:`P_0 +
        (n,2n)` source on the typed
        :class:`~orpheus.transport.sources.ScalarSource` accumulator,
        build the per-ordinate :math:`P_\ell\ge 1` Galerkin
        contribution via :meth:`build_aniso_source` (now accepting
        the L2 :class:`AngularFlux` on the composite's bulk), then
        combine via the producer-side :math:`1/W` projection
        (Pattern 7).  The output bulk is a pure-Field
        :class:`AngularFlux`; the output boundary is the
        **implicit-zero** :class:`BoundaryFlux` — scattering is
        volumetric, no face-trace contribution.

        Per Option β3 (`#208
        <https://github.com/deOliveira-R/ORPHEUS/issues/208>`_) the
        implicit-zero boundary is a transitional shim: Wave O will
        introduce :class:`BulkOperator` / :class:`FullOperator`
        Protocols so that scattering's volumetric nature is encoded
        in the *type*, not in a zero-valued boundary member.  Until
        then the composite return enables ``L.apply(state) -
        S.apply(state) - F.apply(state)`` to compose under
        :meth:`TimedFullField.__sub__` once all four operators expose
        the composite branch (D-H.1c).
        """
        bulk = psi.bulk
        mesh = bulk.mesh
        # Reduce bulk → scalar via the L2 typed reduction (the L2
        # ``AngularFlux.integrate_angular`` returns the same
        # :class:`ScalarFlux` type the legacy branch consumes).
        phi = bulk.integrate_angular()
        # iso = P0 in-scatter + (n,2n) doubling on the typed accumulator.
        iso: ScalarSource = mesh.zeros_scalar_source()
        iso = self.add_iso_source(iso, phi)
        iso = self.add_n2n_source(iso, phi)
        # Pℓ (ℓ≥1) — per-ordinate after R-1 Step 4 A1 producer-side /W.
        # ``build_aniso_source`` now accepts the L2 ``AngularFlux`` on
        # the composite's bulk (the type-check widening above).
        aniso_or_none = self.build_aniso_source(bulk)
        if aniso_or_none is None:
            aniso = mesh.zeros_angular_source()
        else:
            aniso = aniso_or_none
        # Producer-side projection at the apply boundary (Pattern 7).
        sum_w = float(mesh.quad.weights.sum())
        combined: AngularSource = (iso / sum_w) + aniso
        return TimedFullField(
            bulk=AngularFlux.from_mesh(combined.values, mesh),
            boundary=BoundaryFlux.zeros_for_sn_mesh(mesh),
            _history=(),
            history_depth=psi.history_depth,
        )

    @apply.register
    def _(self, phi: ScalarFlux) -> "ScalarSource":
        r"""Typed ScalarFlux variant — iso scalar magnitude output (P0 + n2n only).

        Math:
        :math:`Q(\vec r, g) = \Sigma_{s,0}(g'\to g)\,\phi_{g'}(\vec r)
        + 2\,\Sigma_{2n}(g'\to g)\,\phi_{g'}(\vec r)`.
        No :math:`P_\ell` (scalar flux lacks the angular info needed
        for the Galerkin reconstruction); no :math:`1/W` (consumers in
        scalar-flux equations — diffusion, CP, kinetics outer — do not
        project to per-ordinate).
        """
        mesh = phi.mesh
        iso: ScalarSource = mesh.zeros_scalar_source()
        iso = self.add_iso_source(iso, phi)
        iso = self.add_n2n_source(iso, phi)
        return iso

    @apply.register
    def _(self, psi: AngularFlux) -> "AngularSource":
        r"""Typed :class:`AngularFlux` variant — per-ordinate magnitude output.

        Math: identical to the :class:`TimedFullField` arm above on the
        bulk axis — reduce ``psi`` angular → scalar, build iso :math:`P_0
        + (n,2n)` source on the typed
        :class:`~orpheus.transport.sources.ScalarSource` accumulator,
        build the per-ordinate :math:`P_\ell\ge 1` Galerkin contribution
        via :meth:`build_aniso_source`, then combine via the producer-side
        :math:`1/W` projection (Pattern 7).

        D-I.2 (2026-05-29) added this typed leaf alongside the existing
        :class:`TimedFullField` composite arm; the legacy bare-ndarray
        :class:`numpy.ndarray` arm retired in the same commit.  Direct
        :class:`AngularFlux` consumers (the within-group iteration when
        scattering is computed on the bulk-only carrier) now have a
        first-class dispatch arm instead of routing through a transient
        :class:`TimedFullField` wrap.
        """
        mesh = psi.mesh
        phi = psi.integrate_angular()
        iso: ScalarSource = mesh.zeros_scalar_source()
        iso = self.add_iso_source(iso, phi)
        iso = self.add_n2n_source(iso, phi)
        aniso_or_none = self.build_aniso_source(psi)
        if aniso_or_none is None:
            aniso = mesh.zeros_angular_source()
        else:
            aniso = aniso_or_none
        sum_w = float(mesh.quad.weights.sum())
        combined: AngularSource = (iso / sum_w) + aniso
        return combined
