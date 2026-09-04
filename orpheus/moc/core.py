"""MOC transport solver satisfying the EigenvalueSolver protocol.

Solves the multi-group neutron transport equation via the Method of
Characteristics on a 2-D pin-cell geometry with reflective boundary
conditions.  The flat-source approximation is used within each FSR,
and the angular discretisation uses a product quadrature (azimuthal x
Tabuchi-Yamamoto polar).

Reference: Boyd et al. (2014) "The OpenMOC Method of Characteristics
Neutral Particle Transport Code", Ann. Nucl. Energy 68, 43-52.
"""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np

from orpheus.transport.kernels import N2NKernel

from orpheus.data.macro_xs.mixture import Mixture
from orpheus.numerics.convergence import (
    IterationBudget,
    IterationRecord,
    StoppingCriterion,
)

from .geometry import MOCMesh


def _relative_increment(new: np.ndarray, old: np.ndarray) -> float:
    r"""``‖new − old‖_F / ‖new‖_F`` — MoC's convergence metric, one home.

    Relative Frobenius (no ``axis=``, so the flat ℓ² over every index
    jointly), with the CURRENT iterate in the denominator and a ``1e-30``
    floor.  Both of MoC's levels judge themselves with it: the outer on the
    scalar flux (:meth:`MOCSolver.measure_stopping_criteria`) and the inner
    on the scalar flux AND the boundary angular flux
    (:meth:`MOCSolver.solve_fixed_source`).  It lives here so "the same
    metric" is a fact about the code rather than a claim about two
    expressions that happen to match today.
    """
    return float(
        np.linalg.norm(new - old) / max(float(np.linalg.norm(new)), 1e-30)
    )


class MOCSolver:
    """2-D Method of Characteristics eigenvalue solver.

    Satisfies the :class:`~numerics.eigenvalue.EigenvalueSolver` protocol.

    Flux representation: ``(n_regions, ng)`` — scalar flux per flat-source
    region and energy group.
    """

    def __init__(
        self,
        moc_mesh: MOCMesh,
        materials: dict[int, Mixture],
        keff_tol: float = 1e-6,
        flux_tol: float = 1e-5,
        inner_tol: float = 1e-8,
        max_inner_sweeps: int = 200,
    ) -> None:
        self.geom = moc_mesh
        self.keff_tol = keff_tol
        self.flux_tol = flux_tol
        self.inner_tol = inner_tol
        self.max_inner_sweeps = max_inner_sweeps
        # #340 N4: one record per outer's sweep loop, in order. APPEND-ONLY
        # across solves (``power_iteration`` slices from its own high-water
        # mark), so never reset this.
        self._inner_records: list[IterationRecord] = []

        nr = moc_mesh.n_regions
        _any_mat = next(iter(materials.values()))
        self.ng = _any_mat.ng

        # Per-region cross sections
        self.sig_t = np.empty((nr, self.ng))
        self.sig_a = np.empty((nr, self.ng))
        self.sig_p = np.empty((nr, self.ng))
        self.chi = np.empty((nr, self.ng))
        self.sig_s0: list[np.ndarray] = []
        self.sig2: list[np.ndarray] = []

        for k in range(nr):
            mat_id = int(moc_mesh.region_mat_ids[k])
            mix = materials[mat_id]
            self.sig_t[k, :] = mix.SigT
            self.sig_a[k, :] = mix.absorption_xs
            self.sig_p[k, :] = mix.SigP
            self.chi[k, :] = mix.chi
            self.sig_s0.append(mix.SigS[0].toarray())
            self.sig2.append(mix.Sig2[0].toarray())  # P0 block: MoC's Q_g is isotropic

        # Persistent boundary angular fluxes (survive between outer iters)
        n_tracks = len(moc_mesh.tracks)
        n_polar = moc_mesh.quad.n_polar
        ng = self.ng
        self._fwd_bflux = np.zeros((n_tracks, n_polar, ng))
        self._bwd_bflux = np.zeros((n_tracks, n_polar, ng))

    # ── EigenvalueSolver protocol ────────────────────────────────────

    def initial_flux_distribution(self) -> np.ndarray:
        return np.ones((self.geom.n_regions, self.ng))

    def compute_fission_source(
        self,
        flux_distribution: np.ndarray,
        keff: float,
    ) -> np.ndarray:
        """Fission source: chi * (SigP . phi) / keff."""
        nr = self.geom.n_regions
        fission_source = np.empty((nr, self.ng))
        for i in range(nr):
            production = self.sig_p[i, :] @ flux_distribution[i, :]
            fission_source[i, :] = self.chi[i, :] * production / keff
        return fission_source

    def solve_fixed_source(
        self,
        fission_source: np.ndarray,
        flux_distribution: np.ndarray,
    ) -> np.ndarray:
        r"""MOC transport sweeps, iterated until the inner is CERTIFIED.

        Each inner sweep:

        1. Build isotropic source Q from fission + scattering + (n,2n)
        2. Sweep all tracks (forward + backward), accumulate delta_phi
        3. Update scalar flux from delta_phi (Boyd Eq. 45)
        4. Measure both increments; stop when BOTH clear ``inner_tol``

        ⛔ **Was a fixed ``n_inner_sweeps`` count with no tolerance, no
        residual and no break until 2026-08-11 (#340 N4).** The
        :class:`~orpheus.numerics.convergence.IterationRecord` primitive
        refuses to let a level that ran and measured nothing claim
        convergence, which is what exposed the schedule as wrong in BOTH
        directions at once: `[M]` at a converged outer the loop needs **1**
        sweep and took **15**, while from a cold boundary flux it needs
        **80–110** and took 15.  A magic constant with measured evidence on
        each side.

        **Why TWO readings.** The loop exists to converge the boundary
        angular flux (the cyclic-track closure), and the scalar flux cannot
        see that.  `[M]` 2026-08-11, ``moc_cyl1D_1eg_1rg`` cold: ``‖Δφ‖`` is
        machine zero from sweep 4 and clears ``1e-5`` at sweep **2**, while
        ``‖Δψ_b‖`` clears only at sweep **18** — φ is a volume moment, so a
        1-group problem (no group coupling to iterate) collapses it in one
        pass while the geometric feedback around reflective tracks is
        untouched.  A break on ``‖Δφ‖`` alone would return with this method's
        own convergence claim four orders unfinished.  Both are declared on
        the record, so *which mode binds* is readable rather than assumed.
        """
        geom = self.geom
        quad = geom.quad
        nr = geom.n_regions
        ng = self.ng

        phi = flux_distribution.copy()

        dphi_trajectory: list[float] = []
        dpsi_trajectory: list[float] = []
        sweeps_run = 0
        # Bound for every path: ``max_inner_sweeps <= 0`` leaves the range
        # empty and the body never runs, so the record below would read an
        # unbound name.  Empty trajectories then say "never measured", which
        # is what actually happened.
        criteria: tuple[StoppingCriterion, ...] = ()

        for sweeps_run in range(1, self.max_inner_sweeps + 1):
            # The boundary angular flux is what this loop EXISTS to converge
            # (the cyclic-track closure), and it is carried in solver state,
            # so snapshot it before the sweep mutates it.
            psi_boundary_old = np.concatenate(
                (self._fwd_bflux.ravel(), self._bwd_bflux.ravel())
            )

            # 1. Build source: Q[i,g] = (1/4pi) * [fission + scatter + n2n]
            Q = np.empty((nr, ng))
            for i in range(nr):
                scatter = self.sig_s0[i].T @ phi[i, :]
                n2n = N2NKernel.multiplicity * self.sig2[i].T @ phi[i, :]
                Q[i, :] = (fission_source[i, :] + scatter + n2n) / (4.0 * np.pi)

            # Q / sig_t (asymptotic angular flux per region)
            q_over_sigt = np.zeros((nr, ng))
            mask = self.sig_t > 0
            q_over_sigt[mask] = Q[mask] / self.sig_t[mask]

            # 2. Sweep all tracks, accumulate delta_phi
            delta_phi = np.zeros((nr, ng))

            for a_idx in range(quad.n_azi):
                ts = geom.effective_spacing(a_idx)
                omega_a = quad.omega_azi[a_idx]

                for t_idx in geom.tracks_per_azi[a_idx]:
                    track = geom.tracks[t_idx]

                    for p_idx in range(quad.n_polar):
                        sin_p = quad.sin_polar[p_idx]
                        omega_p = quad.omega_polar[p_idx]
                        weight = 4.0 * np.pi * omega_a * omega_p * ts * sin_p

                        # --- Forward sweep ---
                        psi = self._fwd_bflux[t_idx, p_idx, :].copy()
                        for seg in track.segments:
                            rid = seg.region_id
                            for g in range(ng):
                                if self.sig_t[rid, g] <= 0:
                                    continue
                                tau = self.sig_t[rid, g] * seg.length / sin_p
                                if tau < 1e-10:
                                    one_minus_exp = tau * (1.0 - 0.5 * tau)
                                else:
                                    one_minus_exp = 1.0 - np.exp(-tau)
                                dpsi = (psi[g] - q_over_sigt[rid, g]) * one_minus_exp
                                psi[g] -= dpsi
                                delta_phi[rid, g] += weight * dpsi

                        # Pass outgoing to linked track
                        fwd_t = track.fwd_link
                        if fwd_t >= 0:
                            if track.fwd_link_fwd:
                                self._fwd_bflux[fwd_t, p_idx, :] = psi
                            else:
                                self._bwd_bflux[fwd_t, p_idx, :] = psi

                        # --- Backward sweep ---
                        psi = self._bwd_bflux[t_idx, p_idx, :].copy()
                        for seg in reversed(track.segments):
                            rid = seg.region_id
                            for g in range(ng):
                                if self.sig_t[rid, g] <= 0:
                                    continue
                                tau = self.sig_t[rid, g] * seg.length / sin_p
                                if tau < 1e-10:
                                    one_minus_exp = tau * (1.0 - 0.5 * tau)
                                else:
                                    one_minus_exp = 1.0 - np.exp(-tau)
                                dpsi = (psi[g] - q_over_sigt[rid, g]) * one_minus_exp
                                psi[g] -= dpsi
                                delta_phi[rid, g] += weight * dpsi

                        bwd_t = track.bwd_link
                        if bwd_t >= 0:
                            if track.bwd_link_fwd:
                                self._fwd_bflux[bwd_t, p_idx, :] = psi
                            else:
                                self._bwd_bflux[bwd_t, p_idx, :] = psi

            # 3. Update scalar flux: Boyd Eq. 45
            #    phi_i = (4pi/sig_t) * [Q_i + delta_phi_i / A_i]
            #    where 4pi is already absorbed into the weight.
            #    So: phi_i = (4pi*Q_i + delta_phi_i/A_i) / sig_t_i
            phi_new = np.zeros((nr, ng))
            for i in range(nr):
                for g in range(ng):
                    if self.sig_t[i, g] > 0:
                        phi_new[i, g] = (
                            4.0 * np.pi * Q[i, g]
                            + delta_phi[i, g] / geom.region_areas[i]
                        ) / self.sig_t[i, g]

            # 4. #340 N4 — the two readings this loop is driven by, in the
            #    SAME metric the outer uses (relative Frobenius, current
            #    iterate in the denominator).  Both, not just the flux:
            #    `[M]` 2026-08-11 on ``moc_cyl1D_1eg_1rg`` from a cold start,
            #    ||dphi|| is MACHINE ZERO from sweep 4 while ||dpsi_b|| is
            #    still 3.5e-02 and does not clear 1e-5 until sweep 18 — a
            #    1-group problem has no group coupling to iterate, so the
            #    scalar flux (a VOLUME MOMENT) is blind to the cyclic-track
            #    closure's slow mode.  Breaking on ||dphi|| alone would stop
            #    with this loop's own stated job four orders unfinished.
            psi_boundary_new = np.concatenate(
                (self._fwd_bflux.ravel(), self._bwd_bflux.ravel())
            )
            dphi_trajectory.append(_relative_increment(phi_new, phi))
            dpsi_trajectory.append(
                _relative_increment(psi_boundary_new, psi_boundary_old)
            )
            phi = phi_new

            # The loop stops on the very objects the record then reports, so
            # "cleared" has ONE definition (StoppingCriterion.cleared) instead
            # of a hand-spelled comparison here that could drift from it.
            criteria = (
                StoppingCriterion(
                    name="dphi", trajectory=tuple(dphi_trajectory),
                    tolerance=self.inner_tol,
                ),
                StoppingCriterion(
                    name="dpsi_boundary", trajectory=tuple(dpsi_trajectory),
                    tolerance=self.inner_tol,
                ),
            )
            if all(criterion.cleared for criterion in criteria):
                break

        self._inner_records.append(IterationRecord(
            label="inner(transport sweeps)",
            criteria=criteria,
            # One unit of ``max_inner_sweeps`` is one transport sweep, so the
            # identity conversion is the honest statement.
            budget=IterationBudget(self.max_inner_sweeps, "max_inner_sweeps"),
            iterations_run=sweeps_run,
        ))

        # Normalise flux (prevent overflow between outer iterations)
        total_prod = 0.0
        for i in range(nr):
            sig2_out = np.array(self.sig2[i].sum(axis=1)).ravel()
            total_prod += (
                (self.sig_p[i, :] + N2NKernel.multiplicity * sig2_out)
                @ phi[i, :] * geom.region_areas[i]
            )
        if total_prod > 0:
            phi *= 1.0 / total_prod
            # ⭐ The boundary angular flux is the OTHER HALF of this solution,
            # and it is carried across outers, so the normalisation has to act
            # on the whole vector.  Rescaling phi alone splits one solution's
            # scale in two: the next outer then inherits a boundary condition
            # inconsistent with its own flux and the inner has to walk the gap
            # back, every time.  (The fixed-source problem is linear, so a
            # global rescale is a symmetry of the pair (phi, psi_b) — acting
            # on one factor only is not.)
            #
            # `[M]` 2026-08-11, total sweeps to a certified solve:
            # 1eg 64 -> 32, 2eg 351 -> 177, 4eg 259 -> 128, with keff error
            # unchanged-or-better (1eg reaches EXACTLY 0.0).  Found because
            # #340 N4 made per-outer sweep counts readable at all.
            self._fwd_bflux *= 1.0 / total_prod
            self._bwd_bflux *= 1.0 / total_prod

        return phi

    def compute_keff(self, flux_distribution: np.ndarray) -> float:
        r"""k of the posed eigenproblem: fission production / net removal.

        .. math::

            k \;=\; \frac{\sum_i A_i\, \nu\Sigma_{f,i}\cdot\phi_i}
                    {\sum_i A_i\, (\Sigma_{a,i} - 2\Sigma_{2,\mathrm{out},i})
                     \cdot\phi_i}

        R7 (#259, 2026-07-03): the MoC inner solve poses the eigenproblem
        with ONLY fission scaled by :math:`1/k` — the (n,2n) emission is
        a plain gain in the sweep source — so the estimator is
        fission-only production over NET removal (``absorption_xs``
        counts the (n,2n) collision once; the emission ``2Σ₂`` is a gain
        that reduces net removal). The previous spelling
        ``(νΣf + 2Σ₂)/Σa`` equals the posed problem's eigenvalue only
        when ``Σ₂ = 0`` or ``k = 1``. Leakage is structurally zero
        (reflective track-linked pin cell). On ``Σ₂ = 0`` mixtures the
        float arithmetic is unchanged (adding/subtracting exact zeros).
        """
        nr = self.geom.n_regions
        p_rate = 0.0
        removal_rate = 0.0
        for i in range(nr):
            A_i = self.geom.region_areas[i]
            phi_i = flux_distribution[i, :]
            sig2_out = np.array(self.sig2[i].sum(axis=1)).ravel()
            p_rate += self.sig_p[i, :] @ phi_i * A_i
            removal_rate += (
                self.sig_a[i, :] - N2NKernel.multiplicity * sig2_out
            ) @ phi_i * A_i
        return p_rate / removal_rate if removal_rate > 0 else 1.0

    def measure_stopping_criteria(
        self,
        keff: float,
        keff_old: float,
        flux_distribution: np.ndarray,
        flux_old: np.ndarray,
    ) -> tuple[StoppingCriterion, ...]:
        r"""``|Δk|`` against ``keff_tol``, relative ``‖Δφ‖₂`` against ``flux_tol``.

        ⛔ Was ``converged(...) -> bool`` until 2026-08-09 (#340 N2b); its
        leading ``if iteration <= 2`` is now
        :data:`~orpheus.numerics.eigenvalue.MINIMUM_OUTER_ITERATIONS`, stated
        once for every solver instead of five times.
        """
        dk = abs(keff - keff_old)
        dphi = _relative_increment(flux_distribution, flux_old)
        return (
            StoppingCriterion.reading("dk", float(dk), self.keff_tol),
            StoppingCriterion.reading("dphi", dphi, self.flux_tol),
        )

    @property
    def inner_records(self) -> "Sequence[IterationRecord]":
        """Every inner sweep loop this instance has run, in order (#340 N4).

        Satisfies the optional
        :class:`~orpheus.numerics.eigenvalue.RecordingSolver` member, so
        :func:`~orpheus.numerics.eigenvalue.power_iteration` splices these in
        as the outer record's children — one per outer.
        """
        return self._inner_records
