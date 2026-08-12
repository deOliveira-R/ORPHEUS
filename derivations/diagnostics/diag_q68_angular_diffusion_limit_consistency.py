r"""Diagnostic: is the shipped curvilinear angular weighting τ the
DIFFUSION-LIMIT-CONSISTENT one, and is τ ≡ ½ measurably not?

Created by numerics-investigator on 2026-08-12 for issue #319 (the flux-dip
discriminator); full evidence in ``scratch/q68_flux_dip_discriminator.md``.

Reference (PRIMARY, structurally independent of ORPHEUS): Morel, J. E. &
Montry, G. R. (1984), *Analysis and Elimination of the Discrete-Ordinates
Flux Dip*, Transport Theory and Statistical Physics 13(5):615–633.  Local
scan ``scratch/literature/Morel-Montry(1984)…pdf``; OCR sidecar
``scratch/literature_ocr/Morel-Montry(1984)…md``.

Three legs, in increasing cost:

``TestBetaIdentity``
    SOLVE-FREE.  M&M Eq. (6a) evaluated with the cell-edge cosines IMPLIED by
    τ must vanish identically for the shipped weighting — that is their
    Eqs. (17)–(19), the whole point of the scheme — and must NOT vanish for
    the plain angular diamond τ ≡ ½ (the negative leg, vv-principles #11/#19:
    a positive reading alone cannot show the gate is loaded on τ).

``TestEffectiveStartingCosineAgainstMorelMontryFigure4``
    L1 against a PUBLISHED FIGURE.  On M&M's own fixture (pure-scattering
    10-mfp sphere, Gauss-S2, 100 zones) their Fig. 4 reports the
    angular-diamond effective starting cosine reaching ≈ −1.35 at the origin,
    and their Fig. 6 reports the weighted-diamond curve at an "essentially
    constant value of −1".

``TestAngularDefectVanishesUnderSpatialRefinement``
    THE DISCRIMINATOR.  M&M's summary (p. 16): the S_N equations become
    equivalent to the diffusion equation in the diffusive limit **as the
    spatial mesh is refined**.  So the shipped τ's starting-cosine defect must
    CONVERGE TO ZERO with h, while τ ≡ ½'s must converge to a finite non-zero
    limit.  This is the axis that separates the two schemes — NOT optical
    thickness, which separates them not at all (measured: the defect is
    constant to four figures over Σ_t·R = 5…50 at fixed cells-per-mfp).

If these catch a real bug, promote to ``tests/sn/`` — the natural homes are
``tests/sn/sweep/`` (the closure seam) for the first leg and
``tests/sn/verification/analytical/`` for the other two.
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.derivations.common.xs_library import make_mixture
from orpheus.geometry.mesh import BC, CoordSystem, Mesh1D
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn import solve_sn_fixed_source
from orpheus.sn.sweep.pole_angular_closure import morel_montry_tau_per_level


# ─────────────────────────────────────────────────────────────────────────
# the instruments
# ─────────────────────────────────────────────────────────────────────────
def morel_montry_beta(mu: np.ndarray, w: np.ndarray, tau: np.ndarray) -> float:
    r"""Morel–Montry Eq. (6a) with the cell-edge cosines IMPLIED by ``tau``.

    .. math::

        \tilde\mu_{m+1/2} = \frac{\mu_m - (1-\tau_m)\,\tilde\mu_{m-1/2}}{\tau_m},
        \qquad \tilde\mu_{1/2} = -1

        \beta = 3 \sum_m \mu_m\bigl(\alpha_{m+1/2}\tilde\mu_{m+1/2}
                                  - \alpha_{m-1/2}\tilde\mu_{m-1/2}\bigr)

    with the standard weight recursion :math:`\alpha_{m+1/2} =
    \alpha_{m-1/2} - \mu_m W_m`, :math:`\alpha_{1/2} = 0`.  ``beta`` is the
    :math:`2\beta/r` corruption of the S\ :sub:`N` diffusion coefficient
    :math:`D = 1/(3(\sigma_t + 2\beta/r))` (their Eq. 7a) — the flux dip.

    ``mu`` ascending; weights renormalised to :math:`\sum W = 1` (M&M's
    convention, in which the diffusion condition reads
    :math:`\sum \mu^2 W = 1/3`).

    ⚠ A ``beta`` built instead from the STANDARD weight-partition edges
    :math:`\mu_{m+1/2} = \mu_{m-1/2} + 2W_m` is τ-BLIND by construction —
    substituting them into Eq. (6a) IS M&M's proof that β = 0.  The τ-loaded
    quantity is the one above.
    """
    mu = np.asarray(mu, dtype=float)
    tau = np.asarray(tau, dtype=float)
    weights = np.asarray(w, dtype=float)
    normalized = weights / weights.sum()
    n = mu.size
    alpha = np.zeros(n + 1)
    for m in range(n):
        alpha[m + 1] = alpha[m] - mu[m] * normalized[m]
    implied = np.zeros(n + 1)
    implied[0] = -1.0
    for m in range(n):
        implied[m + 1] = (mu[m] - (1.0 - tau[m]) * implied[m]) / tau[m]
    return float(
        3.0 * np.sum(mu * (alpha[1:] * implied[1:] - alpha[:-1] * implied[:-1]))
    )


def effective_starting_cosine(solution, quad) -> np.ndarray:
    r"""Morel–Montry's effective starting cosine, Eq. (10), level-locally.

    The diffusion limit makes the angular flux affine in the level's own
    angular coordinate, :math:`\psi_{\ell,m} = a_\ell + k_\ell c_m` with
    :math:`c_m = \mu_m/\sin\theta_\ell`, and puts the starting direction at
    :math:`c = -1`.  So

    .. math::

        c_s(r) = \frac{\psi_{1/2}(r) - a(r)}{k(r)} \;\longrightarrow\; -1 .

    On the sphere this is M&M Eq. (10) verbatim.  The level-LOCAL moments are
    mandatory on a cylinder: its on-axis flux is azimuth-independent but
    genuinely polar-angle dependent (M&M appendix), so the all-level form
    reads a large artefact there.
    """
    psi = np.asarray(solution.angular_flux.interior.values, dtype=float).sum(axis=1)
    weights = np.asarray(quad.weights, dtype=float)
    mu_radial = np.asarray(quad.mu_x, dtype=float)
    interior = solution.radial_characteristic.interior
    level = sorted(solution.mesh.pole_angular_closure._carrying_levels)[0]
    idx = (np.arange(psi.shape[0]) if quad.level_structure is None
           else np.asarray(quad.level_indices[level]))
    sin_theta = (1.0 if quad.level_structure is None
                 else float(np.sqrt(1.0 - float(quad.level_mu[level]) ** 2)))
    cosine = mu_radial[idx] / sin_theta
    w_level = weights[idx]
    mean = np.einsum("n,nx->x", w_level, psi[idx]) / w_level.sum()
    slope = (np.einsum("n,nx->x", w_level * cosine, psi[idx])
             / float(np.sum(w_level * cosine**2)))
    psi_start = np.asarray(interior.cells(level, -1), dtype=float).sum(axis=0)
    return (psi_start - mean) / slope


# ─────────────────────────────────────────────────────────────────────────
# the fixture: Morel–Montry's own test problem
# ─────────────────────────────────────────────────────────────────────────
def _pure_scatterer(sigma_t: float):
    total = np.array([sigma_t])
    scatter = np.array([[sigma_t]])
    mixture = make_mixture(
        sig_t=total, sig_c=total - scatter.sum(axis=1), sig_f=np.zeros(1),
        nu=np.zeros(1), chi=np.zeros(1), sig_s=scatter,
        sig_s1=np.zeros((1, 1)),
    )
    mixture.assert_balanced()
    return mixture


def solve_uniform_source_ball(*, sigma_t, radius, nx, quad, tau_transform=None):
    """Uniform isotropic source in a homogeneous pure-scattering ball.

    ``tau_transform`` monkeypatches the angular weighting IN PROCESS (the
    only supported swap); it is always restored.
    """
    from orpheus.sn.sweep import pole_angular_closure as pac

    shipped = pac.morel_montry_tau_per_level
    if tau_transform is not None:
        def producer(q, coord):
            return tuple(tau_transform(np.asarray(t, dtype=float))
                         for t in shipped(q, coord))
        pac.morel_montry_tau_per_level = producer
    try:
        edges = np.linspace(0.0, radius, nx + 1)
        mesh = Mesh1D(edges=edges, mat_ids=np.zeros(nx, dtype=int),
                      coord=CoordSystem.SPHERICAL, bc_right=BC("vacuum"))
        source = np.full((quad.N, 1, nx), 1.0 / float(quad.weights.sum()))
        solution = solve_sn_fixed_source(
            {0: _pure_scatterer(sigma_t)}, mesh, quad, source,
            max_inner=400_000, inner_tol=1e-12,
        )
    finally:
        pac.morel_montry_tau_per_level = shipped
    converged = solution.converged
    assert bool(converged() if callable(converged) else converged), (
        "the fixed-source solve did not converge — every number below would "
        "be a statement about the stopping test, not about the scheme"
    )
    return solution


_DIAMOND = staticmethod(lambda t: np.full(t.size, 0.5))


# ─────────────────────────────────────────────────────────────────────────
class TestBetaIdentity:
    r"""Solve-free: β ≡ 0 for the shipped τ, β ≠ 0 for τ ≡ ½."""

    @pytest.mark.parametrize("n", [2, 4, 8, 16, 32])
    def test_shipped_tau_annihilates_beta(self, n: int) -> None:
        quad = Quadrature.gauss_legendre(n)
        order = np.argsort(quad.mu_x)
        tau = np.asarray(
            morel_montry_tau_per_level(quad, CoordSystem.SPHERICAL)[0], float
        )
        beta = morel_montry_beta(quad.mu_x[order], quad.weights[order], tau)
        assert abs(beta) < 1e-13, (
            f"GL-S{n}: the shipped tau gives beta = {beta:.6e}, but "
            "Morel-Montry Eqs. (17)-(19) prove it is identically zero — a "
            "non-zero value means the shipped tau is NOT the barycentric "
            "coordinate of the ordinate between the standard weight-partition "
            "edges, i.e. the angular scheme is no longer the first-order "
            "diffusion-limit-consistent one."
        )

    @pytest.mark.parametrize(
        "n, beta_min",
        [(2, 1e-1), (4, 1e-3), (8, 1e-5), (16, 1e-7)],
    )
    def test_plain_diamond_does_not(self, n: int, beta_min: float) -> None:
        """NEGATIVE leg — without it the gate above is not shown to be
        LOADED on tau (vv-principles #19)."""
        quad = Quadrature.gauss_legendre(n)
        order = np.argsort(quad.mu_x)
        beta = morel_montry_beta(
            quad.mu_x[order], quad.weights[order], np.full(n, 0.5)
        )
        assert abs(beta) > beta_min, (
            f"GL-S{n}: plain angular diamond gives beta = {beta:.6e}, which "
            f"should exceed {beta_min:.0e}. beta falls ~4 orders per doubling "
            "of N, so a too-small value here means the instrument stopped "
            "reading tau, not that the diamond became consistent."
        )


class TestEffectiveStartingCosineAgainstMorelMontryFigure4:
    r"""L1 against the published figures of Morel & Montry (1984).

    Fixture is theirs: pure scattering (c = 1), 10-mfp sphere, Gauss-S2,
    100 spatial zones, spatial diamond, vacuum surface.
    """

    @pytest.mark.slow
    def test_weighted_diamond_reaches_minus_one_and_plain_diamond_does_not(
        self,
    ) -> None:
        quad = Quadrature.gauss_legendre(2)
        kw = dict(sigma_t=10.0, radius=1.0, nx=100, quad=quad)
        shipped = effective_starting_cosine(
            solve_uniform_source_ball(**kw), quad)
        diamond = effective_starting_cosine(
            solve_uniform_source_ball(**kw, tau_transform=_DIAMOND), quad)

        # M&M Fig. 6: the weighted-diamond curve sits at "an essentially
        # constant value of -1" in the interior.
        interior = slice(1, 60)          # 0.15 .. 5.95 mfp from the origin
        assert np.max(np.abs(shipped[interior] + 1.0)) < 2e-2, (
            "shipped tau: max |c_s + 1| over the interior = "
            f"{np.max(np.abs(shipped[interior] + 1.0)):.4e} (measured 1.1e-3 "
            "on 2026-08-12); Morel-Montry Fig. 6 requires ~ -1."
        )
        # M&M Fig. 4: the angular-diamond curve reaches ~ -1.35 at the origin.
        assert -1.45 < diamond[0] < -1.25, (
            f"plain diamond: c_s(r_0) = {diamond[0]:.4f}; Morel-Montry Fig. 4 "
            "reports 'a minimum value at the origin of approximately -1.35' "
            "(measured -1.3157 on 2026-08-12)."
        )
        assert abs(diamond[0] + 1.0) > 25.0 * abs(shipped[0] + 1.0), (
            "the two schemes' origin defects are not separated: "
            f"shipped {shipped[0] + 1.0:+.4e} vs diamond "
            f"{diamond[0] + 1.0:+.4e} (measured 1.16e-2 vs 3.16e-1, 27x)"
        )


class TestAngularDefectVanishesUnderSpatialRefinement:
    r"""THE discriminator: refine h, not the optical thickness.

    Morel-Montry p. 16 (SUMMARY): the S_N equations become equivalent to the
    diffusion equation in the diffusive limit **as the spatial mesh is
    refined**.  ⟹ the shipped tau's defect must go to zero with h (measured:
    first order, rate -0.95); tau = 1/2's must saturate (measured: it GROWS,
    rate +0.24, to a limit of ~0.39).
    """

    @pytest.mark.slow
    def test_shipped_converges_and_diamond_saturates(self) -> None:
        quad = Quadrature.gauss_legendre(2)
        meshes = (40, 80, 160)          # 4, 8, 16 cells/mfp at sigma_t*R = 10
        shipped, diamond = [], []
        for nx in meshes:
            kw = dict(sigma_t=10.0, radius=1.0, nx=nx, quad=quad)
            shipped.append(abs(
                effective_starting_cosine(
                    solve_uniform_source_ball(**kw), quad)[0] + 1.0))
            diamond.append(abs(
                effective_starting_cosine(
                    solve_uniform_source_ball(**kw, tau_transform=_DIAMOND),
                    quad)[0] + 1.0))

        h = np.log(1.0 / np.asarray(meshes, dtype=float))
        rate_shipped = float(np.polyfit(h, np.log(shipped), 1)[0])
        rate_diamond = float(np.polyfit(h, np.log(diamond), 1)[0])

        assert rate_shipped > 0.85, (
            "the shipped tau's starting-cosine defect must vanish at ~first "
            f"order in h; measured order {rate_shipped:.3f} over "
            f"{meshes} (values {['%.3e' % v for v in shipped]}; 2026-08-12 "
            "read 0.94)"
        )
        assert rate_diamond < 0.10, (
            "plain angular diamond's defect must NOT converge — it saturates "
            f"at a finite limit; measured order {rate_diamond:.3f} over "
            f"{meshes} (values {['%.3e' % v for v in diamond]}; 2026-08-12 "
            "read -0.24, i.e. it grows)"
        )
        assert diamond[-1] / shipped[-1] > 20.0, (
            "at the finest mesh the two schemes must be separated by more "
            f"than 20x; measured {diamond[-1] / shipped[-1]:.1f}x "
            f"({diamond[-1]:.3e} vs {shipped[-1]:.3e}; 2026-08-12 read 46x)"
        )

    @pytest.mark.slow
    def test_optical_thickness_is_NOT_the_discriminating_axis(self) -> None:
        """The refutation of #319's stated discriminator, pinned.

        At CONSTANT cells-per-mfp the angular-consistency defect is
        thickness-independent to four figures for BOTH schemes, so their
        decay RATES with optical thickness are identical (zero) and cannot
        separate them.  If this ever stops holding, the framing in
        ``scratch/q68_flux_dip_discriminator.md`` needs revisiting.
        """
        quad = Quadrature.gauss_legendre(2)
        defects = {"shipped": [], "diamond": []}
        thicknesses = (5.0, 10.0, 20.0)
        for tau_r in thicknesses:
            nx = int(round(4 * tau_r))          # 4 cells/mfp, held fixed
            kw = dict(sigma_t=tau_r, radius=1.0, nx=nx, quad=quad)
            defects["shipped"].append(
                effective_starting_cosine(
                    solve_uniform_source_ball(**kw), quad)[0] + 1.0)
            defects["diamond"].append(
                effective_starting_cosine(
                    solve_uniform_source_ball(**kw, tau_transform=_DIAMOND),
                    quad)[0] + 1.0)
        for name, values in defects.items():
            spread = np.ptp(values) / np.abs(np.mean(values))
            assert spread < 5e-3, (
                f"{name}: the starting-cosine defect varied by {spread:.2e} "
                f"relative over sigma_t*R = {thicknesses} at fixed "
                f"cells-per-mfp ({['%+.5e' % v for v in values]}). It was "
                "constant to four figures on 2026-08-12; a thickness "
                "dependence here would reopen #319's decay-rate framing."
            )
