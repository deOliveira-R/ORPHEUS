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
from typing import TYPE_CHECKING

import numpy as np

from orpheus.numerics.operator import (
    CAP_APPLY,
    LinearOperatorMixin,
)
from orpheus.numerics.projection import (
    HarmonicMomentProjection,
    HarmonicMomentReconstruction,
)

if TYPE_CHECKING:
    from .angular_flux import AngularFlux
    from .material_xs_field import MaterialXSField
    from .quadrature import AngularQuadrature
    from .scalar_flux import ScalarFlux
    from .sources import IsotropicSource, PerOrdinateSource


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

        Notes
        -----
        Per-material dispatch is encapsulated by
        :meth:`MaterialXSField.apply_legendre_scattering_moments` —
        Issue #197 PR-TYPED-1 collapses the ``for mid, (ix, iy) in
        cells_by_mat.items()`` loop into a typed verb.
        """
        return self.mat_xs.apply_legendre_scattering_moments(
            moments, L=self.L, skip_l0=self.skip_l0,
        )


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
        Q: "np.ndarray | IsotropicSource",
        phi: "np.ndarray | ScalarFlux",
    ) -> "np.ndarray | IsotropicSource | None":
        r"""Add P0 in-scatter source to :math:`Q`.

        Vectorised by material: per cell ``c`` of material ``mid``,
        ``Q[:, ic, jc] += sig_s0[mid].T @ phi[:, ic, jc]`` where
        ``sig_s0[mid]`` is ``(ng, ng)`` indexed ``[g_from, g_to]``.

        Issue #197 PR-TYPED-3 introduces the typed-action overload:

        * **Raw-in (legacy)** — ``Q: np.ndarray`` is mutated in place
          and ``None`` is returned (the Wave A–D / PR-TYPED-1 contract).
        * **Typed-in (return-new)** — ``Q: IsotropicSource`` is treated
          as an immutable input; a NEW :class:`IsotropicSource` is
          returned carrying ``Q.values + in_scatter`` (Pattern 4 —
          frozen typed inputs stay immutable; Pattern 1 — the caller
          spells the algebra as ``Q = scattering.add_iso_source(Q,
          phi)``).

        Parameters
        ----------
        Q : np.ndarray or IsotropicSource
            Isotropic source carrier.  Raw ``(ng, nx, ny)`` ndarray is
            mutated in place; typed :class:`IsotropicSource` returns a
            new instance.
        phi : np.ndarray or ScalarFlux
            Scalar flux.  Either form is accepted; the underlying
            values are unwrapped before the per-material dispatch.

        Returns
        -------
        np.ndarray or IsotropicSource or None
            Raw-in returns ``None`` (legacy in-place); typed-in returns
            a fresh :class:`IsotropicSource`.

        Notes
        -----
        Issue #197 PR-TYPED-1 — the per-material dispatch lives ONLY
        inside :meth:`MaterialXSField.apply_p0_in_scatter`.
        """
        from .scalar_flux import ScalarFlux
        from .sources import IsotropicSource
        phi_values = phi.values if isinstance(phi, ScalarFlux) else phi
        if isinstance(Q, IsotropicSource):
            Q_values = Q.values.copy()
            self.mat_xs.apply_p0_in_scatter(Q_values, phi_values)
            return IsotropicSource(Q_values, Q.mesh)
        self.mat_xs.apply_p0_in_scatter(Q, phi_values)
        return None

    def add_n2n_source(
        self,
        Q: "np.ndarray | IsotropicSource",
        phi: "np.ndarray | ScalarFlux",
    ) -> "np.ndarray | IsotropicSource | None":
        r"""Add (n,2n) source to :math:`Q`.

        Vectorised by material: per cell ``c`` of material ``mid``,
        ``Q[:, ic, jc] += 2 · sig2[mid].T @ phi[:, ic, jc]``.

        Issue #197 PR-TYPED-3 introduces the same typed-action overload
        as :meth:`add_iso_source` — raw-in mutates in place and returns
        ``None``; typed-in returns a fresh :class:`IsotropicSource`.

        Parameters
        ----------
        Q : np.ndarray or IsotropicSource
            Isotropic source carrier.
        phi : np.ndarray or ScalarFlux
            Scalar flux.

        Returns
        -------
        np.ndarray or IsotropicSource or None
            Raw-in returns ``None`` (legacy in-place); typed-in returns
            a fresh :class:`IsotropicSource`.

        Notes
        -----
        Issue #197 PR-TYPED-1 — per-material dispatch lives ONLY
        inside :meth:`MaterialXSField.apply_n2n`.
        """
        from .scalar_flux import ScalarFlux
        from .sources import IsotropicSource
        phi_values = phi.values if isinstance(phi, ScalarFlux) else phi
        if isinstance(Q, IsotropicSource):
            Q_values = Q.values.copy()
            self.mat_xs.apply_n2n(Q_values, phi_values)
            return IsotropicSource(Q_values, Q.mesh)
        self.mat_xs.apply_n2n(Q, phi_values)
        return None

    def build_aniso_source(
        self,
        angular_flux: "np.ndarray | AngularFlux | None",
    ) -> "np.ndarray | PerOrdinateSource | None":
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
        angular_flux : np.ndarray or AngularFlux or None
            Angular flux shape ``(N, ng, nx, ny)`` (Issue #196
            PR-INDEX-4 — principled) from the most recent sweep, or
            ``None`` on the first source iteration before any sweep
            has been run.  Issue #197 PR-TYPED-3 — typed
            :class:`AngularFlux` is accepted; the result is then a
            typed :class:`PerOrdinateSource` (preserving the type chain
            through the scattering composition).

        Returns
        -------
        np.ndarray or PerOrdinateSource or None
            ``(N, ng, nx, ny)`` per-ordinate :math:`\ell \ge 1`
            contribution, or ``None`` when ``scattering_order == 0`` or
            ``angular_flux is None``.  Type matches the input.

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
        from .angular_flux import AngularFlux
        from .sources import PerOrdinateSource
        if self.scattering_order == 0 or angular_flux is None:
            return None
        is_typed = isinstance(angular_flux, AngularFlux)
        values = angular_flux.values if is_typed else angular_flux
        L = self.scattering_order
        # Build the §9 "S = R Λ M" pipeline. The constituent primitives
        # are cheap dataclass instantiations; the actual work is in the
        # three np.einsum calls inside their .apply methods.
        Y = self.Y  # cached on first access
        M = HarmonicMomentProjection(weights=self.weights, Y=Y, L=L)
        Lam = LegendreMomentScattering(
            mat_xs=self.mat_xs,
            L=L,
            skip_l0=True,
        )
        R = HarmonicMomentReconstruction.from_Y(Y)
        result = R.apply(Lam.apply(M.apply(values)))
        if is_typed:
            return PerOrdinateSource(result, angular_flux.mesh)
        return result

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

    def apply(
        self, psi: "np.ndarray | AngularFlux",
    ) -> "np.ndarray | AngularFlux":
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
        psi : np.ndarray or AngularFlux
            Angular flux.  When typed as
            :class:`~orpheus.sn.angular_flux.AngularFlux` (Issue #197
            PR-TYPED-2), the result is also an :class:`AngularFlux` —
            the operator reads as the math.  When passed as bare
            ``np.ndarray`` shape ``(N, ng, nx, ny)``, the result is a
            bare ``np.ndarray`` — preserved for legacy
            packed-vector / FD-matvec call sites.

        Returns
        -------
        np.ndarray or AngularFlux
            Per-ordinate scattering source shape
            ``(N, ng, nx, ny)``.  Type matches the input.

        Notes
        -----
        ORPHEUS's existing :class:`SNSolver` source iteration consumes
        the isotropic and anisotropic pieces separately (the isotropic
        :math:`Q` is a ``(ng, nx, ny)`` scalar source; the anisotropic
        :math:`Q_{\rm aniso}` is a ``(N, ng, nx, ny)`` per-ordinate
        source).  For internal use by :class:`SNSolver` the helpers
        :meth:`add_iso_source`, :meth:`add_n2n_source`, and
        :meth:`build_aniso_source` expose the same math without forcing
        a wasteful broadcast.
        """
        from .angular_flux import AngularFlux
        from .scalar_flux import ScalarFlux
        from .sources import IsotropicSource, PerOrdinateSource
        is_typed = isinstance(psi, AngularFlux)
        values = psi.values if is_typed else psi
        if is_typed:
            mesh = psi.mesh
            # Reduce angular flux to scalar flux via the typed reduction
            # (Pattern 3 — named intermediate).
            phi: "ScalarFlux | np.ndarray" = psi.integrate_angular()
            # Build the isotropic source as a typed IsotropicSource
            # accumulator (Pattern 1 — read-as-the-math).  P0 + (n,2n)
            # both return-new under the typed overload.
            iso: IsotropicSource = mesh.zeros_isotropic_source()
            iso = self.add_iso_source(iso, phi)
            iso = self.add_n2n_source(iso, phi)
            # Pℓ (l>=1) contribution as a typed PerOrdinateSource.
            aniso_or_none = self.build_aniso_source(psi)
            if aniso_or_none is None:
                aniso = mesh.zeros_per_ordinate_source()
            else:
                aniso = aniso_or_none  # type: PerOrdinateSource
            # Cross-type dunder algebra: iso + aniso broadcasts the
            # isotropic source across the N axis and adds the
            # anisotropic source — the load-bearing PR-TYPED-3
            # pattern (replaces ``np.broadcast_to(...).copy() + Q``).
            combined: PerOrdinateSource = iso + aniso
            return AngularFlux(combined.values, mesh)
        # Bare-ndarray path preserves the legacy in-place contract for
        # packed-vector / FD-matvec consumers.  The implicit-broadcast
        # add (``Q_iso[None, :, :, :] + Q_aniso``) replaces the explicit
        # ``np.broadcast_to(...).copy()`` followed by ``+=`` — numpy
        # produces a fresh ``(N, ng, nx, ny)`` allocation from the
        # broadcast add in one step.  When ``Q_aniso is None`` (P0-only),
        # add against an implicit zero (``Q_iso[None, :, :, :] +
        # np.zeros((1,1,1,1))``) — still elementwise add, no broadcast
        # view, owns its data, ``(N, ng, nx, ny)`` shape.
        phi_arr = np.einsum('n,ngxy->gxy', self.weights, values)
        Q_iso = np.zeros((self.ng, self.nx, self.ny))
        self.add_iso_source(Q_iso, phi_arr)
        self.add_n2n_source(Q_iso, phi_arr)
        Q_aniso = self.build_aniso_source(values)
        if Q_aniso is None:
            # Explicit allocation + broadcast assignment.  Pattern 3 —
            # the result owns its memory; downstream callers may mutate.
            N = values.shape[0]
            Q = np.empty((N, self.ng, self.nx, self.ny), dtype=Q_iso.dtype)
            Q[...] = Q_iso[None, :, :, :]
            return Q
        return Q_iso[None, :, :, :] + Q_aniso
