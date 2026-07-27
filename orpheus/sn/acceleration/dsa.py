r"""Consistent Diffusion Synthetic Acceleration — the low-order half (#2).

The production realization of the slab-DD **consistent** low-order
system and the correction operator that wires it into the within-group
iteration. "Consistent" is a theorem, not a slogan: the low-order
operator here is the one *derived* from the assembled DD transport
system by Larsen's four-step moment reduction — reduce-the-discrete,
never discretize-the-reduced — proven symbolically in
:mod:`orpheus.derivations.discrete.sn.dsa` (the algebra of record) and
realized numerically by its ``build_consistent_dd_system``. This module
is pinned entry-for-entry against that reference builder
(``tests/sn/acceleration/test_dsa_low_order.py``); it transcribes no
formula the derivation has not proven.

Why the edge-centered system and not the diffusion module's operator
====================================================================

The R4 ruling (2026-07-26, evidence in
``.claude/plans/dsa_d2_characterization.md``): the landed cell-centered
diffusion loss (:mod:`orpheus.diffusion`) is the right *standalone*
discretization (RT0/harmonic-mean, M-matrix) but the WRONG correction
operator — measured :math:`\rho` up to 54.7 (divergent) for
:math:`\sigma_t h \ge 2`, the Adams & Larsen (3.43)/(3.44) inconsistent
class. The derived edge-centered system is unconditionally stable
(measured :math:`\rho \le 0.181` over :math:`\sigma_t h \in [0.1, 30]`,
:math:`c \le 0.99`, vacuum and reflective; S2 one-iteration exactness at
machine zero). The two coexist by defining law — standalone accuracy vs
consistency-by-derivation — not as twin paths.

The system lives on the SN side (not in :mod:`orpheus.diffusion`)
because its coefficients are properties of the SN discretization: the
Marshak row carries the DISCRETE half-range moment :math:`\gamma_N =
\sum_{\mu_m>0}\mu_m\omega_m` (≈ 0.2606 for S4 GL — measurably not the
continuum ¼), the closure carries the quadrature moment :math:`W_2`, and
the weighted-diamond generalization enters through the scheme's
:math:`\rho = L_1[\alpha]` (zero for DD, arm 1's scope).

The algebra
===========

Per group :math:`g`, the correction unknowns are edge values
:math:`f_0` on the :math:`K+1` mesh edges. The interior row at the
shared edge between cells *lo* and *hi* is Larsen (27) with the DD
(:math:`\alpha=0`, :math:`\rho=0`) coefficients (23a–f):

.. math::

    -\frac{D_{hi}}{h_{hi}}(f_{0,R} - f_{0,S})
    + \frac{D_{lo}}{h_{lo}}(f_{0,S} - f_{0,L})
    + \tfrac14\bigl[\hat\sigma_{R,hi}h_{hi}(f_{0,R}+f_{0,S})
                   + \hat\sigma_{R,lo}h_{lo}(f_{0,S}+f_{0,L})\bigr]
    = \tfrac12(g_{0,hi} + g_{0,lo}) - (g_{1,hi} - g_{1,lo})

with :math:`\hat\sigma_R = \sigma_t - \sigma_{s0}^{g\to g}`,
:math:`D = 1/[3(\sigma_t - \sigma_{s1}^{g\to g})]`, and the sources
:math:`g_0 = \hat\sigma_S h\, d_0`, :math:`g_1 = a\, d_1` built from the
**raw** sweep displacement moments :math:`d_n = \phi_n^{l+1/2} -
\phi_n^{l}` (the :math:`\sigma`-weighting lives in :math:`G`, once).
Boundary rows close one-sidedly: Marshak
:math:`\gamma_N f_0 + (W_2^+/W_2) f_1 = 0` under vacuum, :math:`f_1 = 0`
under reflection. The accelerated cell updates are Larsen (28): (28a)
collapses for DD to the edge average
:math:`f_{0i} = \tfrac12(f_{0,i-1/2} + f_{0,i+1/2})`, and — when the
sweep retains :math:`\ell \ge 1` (the P1-DSA arm, R5 ruling) — (28b)
gives the moment-1 update :math:`f_{1i} = -(D_i/h_i)\,\Delta f_0 +
a_i d_{1,i}`, injected through Larsen's (33) synthesis
:math:`\Psi_m = f_0 + (\mu_m/W_2) f_1` (measured payoff: the
anisotropy ladder's 86-iteration worst rung returns to 15, the flat
Adams-Larsen band).

Note the structural consequences of the edge home: the conductances
:math:`D_i/h_i` are **per cell** — an edge unknown straddles cells,
which are material-homogeneous by construction, so no harmonic-mean
interpolation exists anywhere in this operator (the cell-centered
standalone solver needs one because its unknowns straddle faces). The
removal is the consistent :math:`\tfrac14(1,2,1)` mass — NOT lumped —
which costs the M-matrix sign pattern at thick cells and buys
unconditional stability.

Two convention boundaries, each guarded numerically:

* the ORPHEUS raw slab quadrature (:math:`\sum w = 2`) maps once to
  Larsen's normalization :math:`\omega = w/2` (asserted, with
  :math:`W_2 = 1/3` — any rule integrating :math:`\mu^2` exactly);
* the within-group data row: :math:`\sigma_{s0}^{g\to g}` via
  :meth:`~orpheus.transport.mesh.material_xs_field.MaterialXSField.foldable_sigma`
  and :math:`\sigma_{s1}^{g\to g}` via
  :meth:`~orpheus.transport.mesh.material_xs_field.MaterialXSField.residual_sig_s`
  — the foldable accessors' first production consumer, reached through
  the LOW-ORDER BUILD only, never the sweep (the #215 trap: folding
  :math:`\sigma_r` into the sweep changes the fixed point for
  anisotropic flux; the low-order operator is correction→0-safe by
  construction, so the fold is legitimate HERE and only here).

The D definition deliberately differs from the standalone diffusion
module's: (23c) uses the within-group P1 **self**-scatter, the outflow
transport approximation uses the **total** P1 out-scatter row — they
coincide only for isotropic data (D2 report, structural diff 4).

References
----------
* E. W. Larsen, "Unconditionally Stable Diffusion-Synthetic
  Acceleration Methods for the Slab Geometry Discrete Ordinates
  Equations. Part I: Theory", *Nucl. Sci. Eng.* **82**, 47–63 (1982) —
  eqs. (23a–f), (27), (28), (38), (39).
* M. L. Adams, E. W. Larsen, "Fast iterative methods for discrete
  ordinates particle transport calculations", *Prog. Nucl. Energy*
  **40**, 3–159 (2002) — §III stability taxonomy, (3.65) rate bound.
* ``.claude/plans/dsa_literature_memo.md`` §6 (build synthesis),
  ``.claude/plans/dsa_d2_characterization.md`` (the R4 evidence).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING

import numpy as np
from scipy.linalg import lu_factor, lu_solve

from orpheus.numerics.operator import LinearOperator

if TYPE_CHECKING:
    from orpheus.sn.mesh.augmented_mesh import SNMesh
    from orpheus.transport.full_field import FullField

__all__ = ["DSALowOrderSystem", "DSACorrection"]


#: Boundary kinds the arm-1 low-order system realizes. Albedo/white walls
#: need the Marshak-albedo generalization of the (38) row — a documented
#: seam, not a silent approximation.
_SUPPORTED_BC = frozenset({"vacuum", "reflective"})


@dataclass(frozen=True)
class DSALowOrderSystem:
    r"""Per-group consistent low-order edge systems for the DD slab sweep.

    Holds, for every group :math:`g`, the LU-factored edge operator
    :math:`A_g` (``(K+1, K+1)``, Larsen (27) interior rows + one-sided
    Marshak/reflecting boundary rows) and the residual-source map
    :math:`G_g` (``(K+1, 2K)``) taking the per-cell displacement moments
    ``[d0; d1]`` to the row sources. Built by :meth:`from_sn_mesh`; the
    correction solve is :meth:`solve_correction`.

    Frozen and immutable — the system is a function of the mesh's
    geometry, data, quadrature, and boundary declarations alone, all of
    which are construction-time constants of the ``SNMesh``.
    """

    #: ``(ng, K+1, K+1)`` assembled edge operators (kept for gates /
    #: ``as_matrix``-style introspection; the solve uses the LU factors).
    a_low: np.ndarray
    #: ``(ng, K+1, 2K)`` displacement→source maps.
    g_map: np.ndarray
    #: Per-group LU factorizations of ``a_low`` (scipy ``lu_factor``).
    _lu: tuple = field(repr=False)
    #: ``(ng, K)`` per-cell :math:`D_i/h_i` — retained for the (28b)
    #: moment-1 update (:meth:`moment1_update`).
    _dh: np.ndarray = field(repr=False)
    #: ``(ng, K)`` per-cell :math:`a = \sigma_{s1}/(\sigma_t-\sigma_{s1})`
    #: — the (23f) weight, shared by the :math:`g_1` sources and (28b).
    _a_coef: np.ndarray = field(repr=False)

    # ── Construction ─────────────────────────────────────────────────

    @classmethod
    def from_sn_mesh(
        cls, sn_mesh: "SNMesh", *, scattering_order: int = 0,
    ) -> "DSALowOrderSystem":
        r"""Build the per-group systems from the SN phase space.

        Admission is arm 1's proven scope, refused loudly outside it:
        1-D Cartesian slab (curvilinear is #282-blocked — no stability
        theory exists; 2-D needs the Alcouffe corner moments), the
        diamond-difference scheme (a weighted-diamond member changes the
        (23) coefficients through :math:`\rho = L_1[\alpha] \ne 0` — the
        R5 seam), and vacuum/reflective walls.

        ``scattering_order`` is the SOLVER's retained Legendre order —
        consistency is with the discrete system BEING ITERATED, so the
        (23c) :math:`\sigma_{s1}^{g\to g}` enters the low-order D and
        the :math:`g_1` weight only when the sweep itself retains
        :math:`\ell \ge 1`; a P0-truncated sweep's consistent partner
        carries the bare-P0 :math:`D = 1/(3\sigma_t)` even when the
        MIXTURE has P1 rows (the data is not the operator).
        """
        from orpheus.geometry.mesh import Mesh1D

        if not (
            sn_mesh.is_cartesian
            and sn_mesh.ndim == 1
            and isinstance(sn_mesh.mesh, Mesh1D)
        ):
            raise NotImplementedError(
                "DSALowOrderSystem: consistent DSA is realized for the "
                "1-D Cartesian slab only (arm 1). Curvilinear is blocked "
                "on #282 (no WDD stability theory exists — Adams & "
                "Larsen p. 79); 2-D Cartesian needs the Alcouffe 9-point "
                "corner moments (follow-up issue at Phase-3 close)."
            )
        scheme_key = getattr(sn_mesh.scheme, "key", None)
        if scheme_key != "diamond_difference":
            raise NotImplementedError(
                f"DSALowOrderSystem: the low-order coefficients are "
                f"proven for the diamond member (α = 0, ρ = 0); scheme "
                f"{scheme_key!r} requires the weighted-diamond (23a–f) "
                f"instantiation at its own ρ (the R5 seam — "
                f"orpheus/derivations/discrete/sn/dsa.py carries the "
                f"general-α theorems)."
            )
        bc_kinds = (
            str(sn_mesh.bc["xmin"].kind), str(sn_mesh.bc["xmax"].kind),
        )
        unsupported = set(bc_kinds) - _SUPPORTED_BC
        if unsupported:
            raise NotImplementedError(
                f"DSALowOrderSystem: boundary kind(s) {sorted(unsupported)!r} "
                f"have no proven low-order row (arm 1 realizes vacuum "
                f"(Marshak) and reflective; albedo walls need the "
                f"Marshak-albedo generalization of Larsen (38))."
            )

        h = np.diff(np.asarray(sn_mesh.mesh.edges, dtype=float))
        mat_xs = sn_mesh.material_xs_field()
        sigma_t = np.asarray(
            mat_xs.total_cross_section_field.values, dtype=float
        )  # (ng, K)
        ng = sigma_t.shape[0]

        mat_ids = np.asarray(sn_mesh.mat_map, dtype=int).ravel()
        fold = mat_xs.foldable_sigma()  # {mid: (ng,)} — σ_s0^{g→g}
        sigma_s0 = np.stack([fold[int(m)] for m in mat_ids], axis=1)
        if scattering_order >= 1:
            residual = mat_xs.residual_sig_s()  # {mid: [cross_P0, Σ_s1, …]}
            s1_diag = {
                mid: (
                    np.asarray(np.diag(mats[1]), dtype=float)
                    if len(mats) > 1
                    else np.zeros(ng)
                )
                for mid, mats in residual.items()
            }
            sigma_s1 = np.stack([s1_diag[int(m)] for m in mat_ids], axis=1)
        else:
            # P0-truncated sweep ⟹ the consistent low-order carries NO
            # σ_s1 (bare-P0 D) — the mixture's P1 rows are data the
            # iterated operator does not contain.
            sigma_s1 = np.zeros_like(sigma_s0)

        mu = np.asarray(sn_mesh.quad.mu_x, dtype=float)
        w = np.asarray(sn_mesh.quad.weights, dtype=float)
        return cls._build(h, sigma_t, sigma_s0, sigma_s1, mu, w, bc_kinds)

    @classmethod
    def _build(
        cls,
        h: np.ndarray,
        sigma_t: np.ndarray,
        sigma_s0: np.ndarray,
        sigma_s1: np.ndarray,
        mu: np.ndarray,
        w: np.ndarray,
        bc: tuple[str, str],
    ) -> "DSALowOrderSystem":
        r"""Assemble the proven rows, vectorized over groups.

        The formula body mirrors the derivation-side reference builder
        (``orpheus.derivations.discrete.sn.dsa.build_consistent_dd_system``)
        row for row — the production tie pins the two entry-for-entry,
        so a drift in either is a red gate, not a silent fork.
        """
        n_cells = h.shape[0]
        ng = sigma_t.shape[0]

        # The ONE normalization map: raw ORPHEUS slab GL (Σw = 2) →
        # Larsen's Σω = 1, guarded numerically at the boundary.
        if not np.isclose(w.sum(), 2.0, rtol=0, atol=1e-12):
            raise ValueError(
                f"DSALowOrderSystem: expected the raw ORPHEUS slab "
                f"quadrature with Σw = 2, got Σw = {w.sum()!r}"
            )
        omega = w / 2.0
        w2 = float(omega @ mu**2)
        if not np.isclose(w2, 1.0 / 3.0, rtol=0, atol=1e-12):
            raise ValueError(
                f"DSALowOrderSystem: the printed convention requires "
                f"W2 = 1/3 exactly (any rule integrating μ² exactly); "
                f"got W2 = {w2!r}"
            )
        up = mu > 0
        gamma_n = float(omega[up] @ mu[up])  # the DISCRETE γ_N
        c1 = float(omega[up] @ mu[up] ** 2) / w2  # W2⁺/W2 (1/2 for GL)

        transport = sigma_t - sigma_s1
        if not (transport > 0.0).all():
            raise ValueError(
                "DSALowOrderSystem: σ_t − σ_s1^{g→g} must be positive "
                "in every (group, cell) — the (23c) D is undefined "
                "otherwise."
            )

        # The proven DD coefficients (23a–f at α = 0) as (ng, K) rows.
        dh = (1.0 / (3.0 * transport)) / h  # D_i / h_i
        quarter = 0.25 * (sigma_t - sigma_s0) * h  # ¼ σ̂_R h
        half_gs = 0.5 * sigma_s0 * h  # ½ σ̂_S h  (g0 weight)
        a_coef = sigma_s1 / transport  # a       (g1 weight)

        a_low = np.zeros((ng, n_cells + 1, n_cells + 1))
        g_map = np.zeros((ng, n_cells + 1, 2 * n_cells))

        # Interior rows (27): edge j between cells lo = j−1 and hi = j.
        j = np.arange(1, n_cells)
        lo, hi = j - 1, j
        a_low[:, j, j - 1] = -dh[:, lo] + quarter[:, lo]
        a_low[:, j, j] = (
            dh[:, lo] + dh[:, hi] + quarter[:, lo] + quarter[:, hi]
        )
        a_low[:, j, j + 1] = -dh[:, hi] + quarter[:, hi]
        g_map[:, j, hi] = half_gs[:, hi]
        g_map[:, j, lo] = half_gs[:, lo]
        g_map[:, j, n_cells + hi] = -a_coef[:, hi]
        g_map[:, j, n_cells + lo] = +a_coef[:, lo]

        # Boundary rows — the verified one-sided f1 forms:
        # f1 = −(D/h)(f0R − f0L) ± ¼σ̂_R h (f0R + f0L) + g1 ∓ ½g0.
        for row, cell_i, side, kind in (
            (0, 0, "left", bc[0]),
            (n_cells, n_cells - 1, "right", bc[1]),
        ):
            sgn = +1.0 if side == "left" else -1.0
            f0_cols = {
                cell_i: +dh[:, cell_i] + sgn * quarter[:, cell_i],
                cell_i + 1: -dh[:, cell_i] + sgn * quarter[:, cell_i],
            }
            g_cols = {
                n_cells + cell_i: +a_coef[:, cell_i],
                cell_i: -sgn * half_gs[:, cell_i],
            }
            if kind == "vacuum":
                # Marshak (38): γ_N f0 + (W2⁺/W2)·(±f1) = 0, outward-
                # consistent orientation.
                orient = +1.0 if side == "left" else -1.0
                a_low[:, row, row] += gamma_n
                for col, val in f0_cols.items():
                    a_low[:, row, col] += orient * c1 * val
                for col, val in g_cols.items():
                    g_map[:, row, col] -= orient * c1 * val
            else:  # "reflective" — admission guarantees the pair
                # (39): f1 at the wall edge = 0.
                for col, val in f0_cols.items():
                    a_low[:, row, col] += val
                for col, val in g_cols.items():
                    g_map[:, row, col] -= val

        lu = tuple(lu_factor(a_low[g]) for g in range(ng))
        return cls(
            a_low=a_low, g_map=g_map, _lu=lu, _dh=dh, _a_coef=a_coef,
        )

    # ── The correction solve ─────────────────────────────────────────

    @property
    def n_groups(self) -> int:
        return self.a_low.shape[0]

    @property
    def n_cells(self) -> int:
        return self.a_low.shape[1] - 1

    def solve_correction(
        self, d0: np.ndarray, d1: np.ndarray | None = None
    ) -> np.ndarray:
        r"""Solve the correction systems for the EDGE values :math:`f_0`.

        Parameters
        ----------
        d0 : (ng, K)
            The raw scalar-moment sweep displacements
            :math:`d_0 = \phi_0^{l+1/2} - \phi_0^{l}` per (group, cell).
        d1 : (ng, K), optional
            The raw moment-1 (current-like) displacements
            :math:`d_1 = \phi_1^{l+1/2} - \phi_1^{l}`. Default zero —
            the P0-DSA arm. The columns carry the (23f) weight
            :math:`a = \sigma_{s1}/(\sigma_t-\sigma_{s1})` (zero for
            isotropic data); :class:`DSACorrection` feeds them whenever
            the sweep retains :math:`\ell \ge 1` (the P1-DSA arm — the
            same consistency-with-the-iterated-operator rule that
            gates the :math:`\sigma_{s1}` data row).

        Returns
        -------
        (ng, K+1) ndarray
            The edge corrections :math:`f_{0,i+1/2}` — the system's
            primary unknowns. The wall columns (``[:, 0]`` / ``[:, -1]``)
            ARE the boundary-flux corrections (the trace arm of
            :class:`DSACorrection` — the lagged reflective gain must
            read the CORRECTED outflow); the cell-average update derives
            via :meth:`cell_update`.
        """
        d0 = np.asarray(d0, dtype=float)
        k = self.n_cells
        if d0.shape != (self.n_groups, k):
            raise ValueError(
                f"solve_correction: d0 must be (ng, K) = "
                f"{(self.n_groups, k)}; got {d0.shape}"
            )
        if d1 is None:
            d1 = np.zeros_like(d0)
        d_stack = np.concatenate([d0, np.asarray(d1, dtype=float)], axis=1)
        rhs = np.einsum("gij,gj->gi", self.g_map, d_stack)
        return np.stack(
            [lu_solve(self._lu[g], rhs[g]) for g in range(self.n_groups)]
        )

    def cell_update(self, f0_edges: np.ndarray) -> np.ndarray:
        r"""The accelerated cell-average update from the edge solution.

        Larsen (28a) at the DD member's :math:`\rho = 0`:
        :math:`\Delta\phi_i = \tfrac12(f_{0,i-1/2} + f_{0,i+1/2})` —
        the :math:`\rho`-proportional terms vanish for diamond
        (``derive_update_relations`` proves the general form).
        """
        return 0.5 * (f0_edges[:, :-1] + f0_edges[:, 1:])

    def moment1_update(
        self, f0_edges: np.ndarray, d1: np.ndarray
    ) -> np.ndarray:
        r"""The accelerated moment-1 cell update — Larsen (28b) at
        :math:`\rho = 0`:

        .. math::

            \Delta\phi_{1,i} = -\frac{D_i}{h_i}
                \bigl(f_{0,i+1/2} - f_{0,i-1/2}\bigr) + a_i\, d_{1,i},

        with the same raw-moment homogeneity as the rest of the system
        (the map is degree-1 homogeneous, so the raw :math:`\Sigma w`
        normalization flows through untouched). The :math:`\rho/2` term
        of the proven general form (``derive_update_relations``)
        vanishes for diamond. The P1-DSA arm's second half: the edge
        solve already carries the :math:`g_1 = a\,d_1` sources; this
        realizes the correction the low-order predicts for the
        current-like moment the :math:`\ell \ge 1` sweep iterates.
        """
        return (
            -self._dh * (f0_edges[:, 1:] - f0_edges[:, :-1])
            + self._a_coef * np.asarray(d1, dtype=float)
        )


class DSACorrection(LinearOperator["FullField", "FullField"]):
    r"""The DSA correction operator on the within-group iterate composite.

    Maps the iterate **displacement** :math:`\Delta\psi = \psi^{l+1/2} -
    \psi^{l}` (a full-angular ``FullField``/``TimedFullField``
    composite) to the additive angular correction:

    .. math::

        \Delta\psi \;\mapsto\; P\,\bigl[(28)\bigr]\,
        A_{\rm low}^{-1}\,G\,R\,\Delta\psi,

    where :math:`R` is the moment restriction, :math:`(28)` the proven
    cell updates, and :math:`P` the moment synthesis — with the arm
    count decided by the SAME consistency rule as the data row
    (``scattering_order``, consistency with the ITERATED operator):

    * **P0 arm** (:math:`\ell = 0` sweep): :math:`R` is the canonical
      moment-0 reduction (:meth:`AngularFlux.integrate_angular
      <orpheus.transport.fields.angular_flux.AngularFlux.integrate_angular>`
      — the frame's ℓ=0 analysis face; pinned 0-ULP in the D8 gate),
      the update is (28a), and :math:`P` the normalized isotropic
      injection :math:`\Delta\phi_0/\sum_n w_n` (the frame's ℓ=0
      reconstruction; moment-0 of :math:`P\,x` is :math:`x` exactly).
    * **P1 arm** (:math:`\ell \ge 1` sweep): the restriction is the
      moment PAIR — the ℓ=1 row is :math:`\sum_n w_n\mu_n(\cdot)` (the
      frame's ℓ=1 analysis row: the SH table's slab component is
      :math:`\mu` bit-exactly) — the solve feeds :math:`[d_0; d_1]`,
      the updates are (28a) AND (28b) (:meth:`DSALowOrderSystem
      .moment1_update`), and :math:`P` is Larsen's (33) ℓ≤1 synthesis
      :math:`\Psi_m = f_0 + (\mu_m/W_2) f_1` under the one raw
      normalization map — :math:`R \circ P = I` on the moment pair by
      the :math:`W_2` quadrature exactness. Without this arm a
      P0-only correction leaves the :math:`\ell = 1` gain iterating
      plainly (measured: the anisotropy ladder climbs 14 → 86
      iterations as :math:`\sigma_{s1}/\sigma_{s0} \to 0.9`).

    No new angular reduction is minted — the anti-mint verdict of the
    3-P0 frame analysis: the frame's ℓ≤1 faces already exist as
    ``integrate_angular``, the :math:`w\mu` row, and the (33)
    synthesis weights; a third spelling would be the Smell-16 twin.

    Both acceleration postures consume this ONE operator:

    * **SI+DSA** — :class:`~orpheus.numerics.iteration.SourceIteration`
      applies it to the sweep displacement each iteration (the
      ``corrector`` parameter);
    * **Krylov-DSA** — the GMRES left preconditioner is
      ``sweep + correction-of-sweep`` (the swept vector IS the
      displacement from zero).

    The returned composite corrects BOTH blocks:

    * **bulk** — the isotropic injection of the (28a) cell update;
    * **trace** — the isotropic injection of the WALL-EDGE solutions
      :math:`f_{0,1/2}, f_{0,K+1/2}`. This arm is load-bearing for
      reflective walls: the production splitting lags the reflection as
      the ``B`` gain reading the ITERATE's outgoing trace, while
      Larsen's reflecting low-order row (:math:`f_1 = 0`) models an
      error equation whose reflected inflow tracks the corrected state
      — leaving the trace uncorrected feeds the next reflection from
      the UN-corrected iterate, and the mismatch measurably diverges
      (:math:`\rho > 1` observed at :math:`c = 0.9`; with the trace arm
      the measured rate matches the D2 scan's stable
      :math:`\rho \approx 0.15`). On vacuum faces the trace arm is
      inert (no boundary gain reads it). An isotropic trace correction
      carries zero net current — consistent with the reflecting row by
      construction.

    At the converged fixed point the displacement vanishes, so the
    correction vanishes — the correction→0 safety property (D6): a bug
    here degrades the RATE, never the answer.
    """

    def __init__(
        self,
        low_order: DSALowOrderSystem,
        sn_mesh: "SNMesh",
        *,
        scattering_order: int = 0,
    ):
        self._low_order = low_order
        self._mesh = sn_mesh
        self._scattering_order = int(scattering_order)
        w = np.asarray(sn_mesh.quad.weights, dtype=float)
        # The ℓ=1 analysis/synthesis coefficients come off the FRAME's
        # own ℓ=1 table row (its slab component IS μ bit-exactly — the
        # D8-family pin), so "the frame's ℓ=1 row" is a CALLED single
        # source, not a claim: one spelling shared with the scattering
        # kernel's moment machinery, never a re-derived w·μ twin.
        mu_row = np.asarray(
            sn_mesh.quad.angular_frame(1).table, dtype=float
        )[:, 1, 1]
        self._sum_w = float(w.sum())
        self._w_mu = w * mu_row
        self._mu = mu_row
        #: 1/W₂ of Larsen's (33) synthesis Ψ = f₀ + (μ/W₂)f₁ — computed
        #: from the quadrature (= 3 exactly under the _build W₂ guard),
        #: never a transcribed constant.
        self._inv_w2 = float(1.0 / ((w / 2.0) @ mu_row**2))

    @classmethod
    def from_sn_mesh(
        cls, sn_mesh: "SNMesh", *, scattering_order: int = 0,
    ) -> "DSACorrection":
        r"""Build the correction operator for an admitted SN phase space
        (admission — geometry, scheme, walls — and the
        ``scattering_order`` consistency rule live on
        :meth:`DSALowOrderSystem.from_sn_mesh`; the same order decides
        whether the ℓ=1 arm of the correction is live — consistency is
        with the iterated operator, in the restriction exactly as in
        the data row)."""
        return cls(
            DSALowOrderSystem.from_sn_mesh(
                sn_mesh, scattering_order=scattering_order,
            ),
            sn_mesh,
            scattering_order=scattering_order,
        )

    @property
    def low_order(self) -> DSALowOrderSystem:
        r"""The per-group low-order systems (gate surface)."""
        return self._low_order

    def apply(self, displacement: "FullField") -> "FullField":
        r"""The correction of one iterate displacement — a DISPLACEMENT.

        The correction is a tangent vector, never a state: the returned
        composite carries :class:`~orpheus.transport.displacements
        .angular_displacement.AngularDisplacement` interior (zero
        boundary displacement), so the update step is the torsor action
        ``ψ ⊕ Δψ_corr`` — ``flux + flux`` stays unspellable (the #208
        affine algebra). The composite's concrete type is preserved
        (``_recombine``; a ``TimedFullField``'s history is reborn empty,
        correct for an addend).

        Two admitted interiors, one per posture:

        * :class:`~orpheus.transport.displacements.angular_displacement
          .AngularDisplacement` — the SI sweep displacement
          ``ψ_{n+1/2} ⊖ ψ_n``; the restriction is its
          ``integrate_angular`` (the tangent map of the canonical
          reduction).
        * :class:`~orpheus.transport.fields.angular_flux.AngularFlux` —
          the Krylov posture's swept vector: GMRES vectors are Krylov
          directions whose role is erased at the scipy ``from_flat``
          boundary (the template rebuilds them flux-typed), and the
          swept vector IS the displacement from zero.

        A moment-windowed carrier (``HarmonicMomentFlux``) refuses
        loudly — 2-D-only, outside the arm-1 admission."""
        from orpheus.transport.displacements.angular_displacement import (
            AngularDisplacement,
        )
        from orpheus.transport.fields.angular_flux import AngularFlux

        interior = displacement.interior
        if not isinstance(interior, (AngularDisplacement, AngularFlux)):
            raise TypeError(
                f"DSACorrection: the input's interior must be a "
                f"full-angular AngularDisplacement (SI posture) or "
                f"AngularFlux (Krylov swept vector); got "
                f"{type(interior).__name__} (a moment-windowed iterate "
                f"is 2-D-only, outside the arm-1 admission)."
            )
        from orpheus.transport.displacements.angular_boundary_displacement import (
            AngularBoundaryDisplacement,
        )

        d0 = interior.integrate_angular().values  # R — the ONE reduction
        if self._scattering_order >= 1:
            # The P1-DSA arm (live iff the sweep retains ℓ ≥ 1 — the
            # SAME consistency rule as the low-order data row): the
            # moment-1 restriction is the frame's ℓ=1 analysis row
            # (w·μ — the SH table's slab component is μ bit-exactly),
            # and the injection is Larsen's (33) ℓ≤1 synthesis
            # Ψ = f₀ + (μ/W₂)f₁ under the one raw-normalization map
            # (moment-0 of the μ-arm vanishes and moment-1 recovers
            # d₁ exactly by the W₂ quadrature guard — R∘P = I on the
            # moment pair).
            d1 = np.einsum("n,ng...->g...", self._w_mu, interior.values)
            f0_edges = self._low_order.solve_correction(d0, d1)
            delta_phi0 = self._low_order.cell_update(f0_edges)
            delta_phi1 = self._low_order.moment1_update(f0_edges, d1)
            angular_values = (
                delta_phi0[None]
                + self._inv_w2 * self._mu[:, None, None] * delta_phi1[None]
            ) / self._sum_w
        else:
            f0_edges = self._low_order.solve_correction(d0)
            delta_phi0 = self._low_order.cell_update(f0_edges)
            angular_values = np.broadcast_to(
                delta_phi0[None] / self._sum_w, interior.values.shape
            ).copy()  # P — normalized isotropic injection
        # The trace arm: the wall-edge solutions injected isotropically
        # per face (see the class docstring — load-bearing under the
        # lagged reflective gain; inert on vacuum faces). The arm stays
        # ℓ=0 by THEOREM even when the ℓ=1 interior arm is live: the
        # reflecting row (39) forces the wall-edge f₁ to zero, and a
        # vacuum wall's trace is read by nothing — so an ℓ=1 trace
        # component is identically zero where it would matter. Minted as
        # the displacement role regardless of the input's role, so the
        # torsor add ``ψ ⊕ Δψ`` is well-formed on both postures (a
        # flux-typed boundary would trip the affine flux+flux gate).
        n_ord = interior.values.shape[0]
        ng = f0_edges.shape[0]
        trace = AngularBoundaryDisplacement.from_face_arrays(
            self._mesh,
            {
                "xmin": np.broadcast_to(
                    f0_edges[None, :, 0] / self._sum_w, (n_ord, ng)
                ).copy(),
                "xmax": np.broadcast_to(
                    f0_edges[None, :, -1] / self._sum_w, (n_ord, ng)
                ).copy(),
            },
        )
        return displacement._recombine(
            interior=AngularDisplacement.from_mesh(angular_values, self._mesh),
            boundary=trace,
        )
