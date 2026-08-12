r"""L1: the shipped curvilinear angular weighting recovers the diffusion limit
under SPATIAL refinement, and plain angular diamond does not — against Morel &
Montry's own published figures.

Promoted 2026-08-12 from
``derivations/diagnostics/diag_q68_angular_diffusion_limit_consistency.py``
(retired in the same commit — no twin). Evidence
``scratch/q68_flux_dip_discriminator.md`` (251 solves, all ``converged`` AND
``fully_converged``, zero warnings, ``max|balance_defect| = 0``); campaign
``.claude/plans/issue_235_angular_accuracy_campaign.md``.

Reference (PRIMARY, structurally independent of ORPHEUS): Morel, J. E. &
Montry, G. R. (1984), *Analysis and Elimination of the Discrete-Ordinates
Flux Dip*, Transport Theory and Statistical Physics **13**(5):615–633. Local
scan ``scratch/literature/Morel-Montry(1984)…pdf``. Fig. 4 reports the
angular-diamond effective starting cosine reaching ≈ −1.35 at the origin;
Fig. 6 reports the weighted-diamond curve at "an essentially constant value
of −1"; p. 16 states that the S_N equations approach the diffusion equation
**as the spatial mesh is refined**.

⭐⭐ **THE AXIS IS h, NOT OPTICAL THICKNESS — and the second row here pins
that as a refutation.** GitHub #319 proposed the dip's *decay rate with
optical thickness* as the discriminator. `[M]` it is not one: at constant
cells-per-mfp the defect is thickness-independent to four figures for the
shipped τ, for τ ≡ ½, **and for a garbage control**. Worse, #319's literal
design self-destructs — at fixed ``c`` the absorption optical depth grows
with thickness, the interior becomes a source plateau, and all schemes agree
to three figures, so a naive run reports "equal ⟹ REFUTED" for a reason
having nothing to do with angular differencing.
:func:`test_optical_thickness_is_NOT_the_discriminating_axis` exists so that
framing cannot be silently re-adopted.

**What these rows CANNOT do.** They are a *constraint* — "does the scheme
preserve the diffusion limit?" — never a ranker. `[M]` on a fixture far from
the diffusion limit (the shipped aniso-cyl MMS is ``c = 0.5``, ``Σ_a·R =
2.5``) the plain diamond can be the *more accurate* scheme. Do not cite
these in an accuracy argument.
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.derivations.common.xs_library import make_mixture
from orpheus.geometry import BC, CoordSystem, Mesh1D
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn import solve_sn_fixed_source

pytestmark = pytest.mark.l1


# ─────────────────────────────────────────────────────────────────────────
# the instrument: Morel–Montry's effective starting cosine, their Eq. (10)
# ─────────────────────────────────────────────────────────────────────────
def effective_starting_cosine(solution, quad) -> np.ndarray:
    r"""M&M Eq. (10), evaluated LEVEL-LOCALLY.

    The diffusion limit makes the angular flux affine in the level's own
    angular coordinate, :math:`\psi_{\ell,m} = a_\ell + k_\ell c_m` with
    :math:`c_m = \mu_m/\sin\theta_\ell`, and puts the starting direction at
    :math:`c = -1`. So

    .. math::

        c_s(r) = \frac{\psi_{1/2}(r) - a(r)}{k(r)} \;\longrightarrow\; -1 .

    On the sphere this is Eq. (10) verbatim. ⚠ The level-LOCAL moments are
    **mandatory** on a cylinder: its on-axis flux is azimuth-independent but
    genuinely polar-angle dependent (M&M's appendix), so the all-level form
    reads a large artefact there (`[M]` ``+2.76``).
    """
    psi = np.asarray(
        solution.angular_flux.interior.values, dtype=float
    ).sum(axis=1)
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


def solve_uniform_source_ball(*, sigma_t, radius, nx, quad,
                              tau_transform=None):
    """Uniform isotropic source in a homogeneous pure-scattering ball.

    ``tau_transform`` monkeypatches the angular weighting IN PROCESS (the
    only supported swap — never a git-level revert); it is always restored.

    ⚠ ``max_inner`` is deliberately enormous. Source iteration's spectral
    radius here is ``c → 1``, so a "reasonable" budget silently returns a
    mid-descent iterate and every number below becomes a statement about
    the stopping test. The explicit ``converged`` assertion is the backstop.
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


def _diamond(t):
    return np.full(t.size, 0.5)


# ─────────────────────────────────────────────────────────────────────────
class TestEffectiveStartingCosineAgainstMorelMontryFigure4:
    r"""L1 against the published figures. Fixture is theirs: pure scattering
    (c = 1), 10-mfp sphere, Gauss-S2, 100 zones, vacuum surface."""

    @pytest.mark.slow
    def test_weighted_diamond_reaches_minus_one_and_plain_diamond_does_not(
        self,
    ) -> None:
        quad = Quadrature.gauss_legendre(2)
        kw = dict(sigma_t=10.0, radius=1.0, nx=100, quad=quad)
        shipped = effective_starting_cosine(
            solve_uniform_source_ball(**kw), quad)
        diamond = effective_starting_cosine(
            solve_uniform_source_ball(**kw, tau_transform=_diamond), quad)

        interior = slice(1, 60)          # 0.15 .. 5.95 mfp from the origin
        assert np.max(np.abs(shipped[interior] + 1.0)) < 2e-2, (
            "shipped tau: max |c_s + 1| over the interior = "
            f"{np.max(np.abs(shipped[interior] + 1.0)):.4e} (measured 1.1e-3 "
            "on 2026-08-12); Morel-Montry Fig. 6 requires ~ -1."
        )
        assert -1.45 < diamond[0] < -1.25, (
            f"plain diamond: c_s(r_0) = {diamond[0]:.4f}; Morel-Montry Fig. 4 "
            "reports 'a minimum value at the origin of approximately -1.35' "
            "(measured -1.3157 on 2026-08-12, i.e. 2.5 % from the figure)."
        )
        assert abs(diamond[0] + 1.0) > 25.0 * abs(shipped[0] + 1.0), (
            "the two schemes' origin defects are not separated: "
            f"shipped {shipped[0] + 1.0:+.4e} vs diamond "
            f"{diamond[0] + 1.0:+.4e} (measured 1.16e-2 vs 3.16e-1, 27x)"
        )


class TestAngularDefectVanishesUnderSpatialRefinement:
    r"""THE discriminator: refine h, not the optical thickness."""

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
                    solve_uniform_source_ball(**kw, tau_transform=_diamond),
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
        separate them. If this ever stops holding, the framing in
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
                    solve_uniform_source_ball(**kw, tau_transform=_diamond),
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
