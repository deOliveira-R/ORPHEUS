r"""The shipped curvilinear angular weight τ IS the diffusion-limit-consistent
one — proved solve-free, on the sphere, by Morel & Montry's own β identity.

Promoted 2026-08-12 from
``derivations/diagnostics/diag_q68_angular_diffusion_limit_consistency.py``
(retired in the same commit — no twin). Evidence:
``scratch/q68_flux_dip_discriminator.md``; campaign
``.claude/plans/issue_235_angular_accuracy_campaign.md``.

Reference (PRIMARY, structurally independent of ORPHEUS): Morel, J. E. &
Montry, G. R. (1984), *Analysis and Elimination of the Discrete-Ordinates
Flux Dip*, Transport Theory and Statistical Physics **13**(5):615–633.
Local scan ``scratch/literature/Morel-Montry(1984)…pdf``. Their
Eqs. (17)–(19) prove β ≡ 0 for the weighted-diamond τ; that proof is the
whole point of the scheme.

⭐ **Why this gate is not in the instrument graveyard with the other five.**
The campaign has killed six plausible τ instruments, all of which turned out
to measure something other than what they were credited with — including
BMC's β, which is **τ-blind**. This β is the *same equation* on a
**different edge set**: the cell edges IMPLIED by τ rather than the standard
weight partition. Substituting the standard edges IS M&M's β = 0 proof, so a
β built that way is blind **by construction**. See
:func:`~orpheus.derivations.discrete.sn.angular_differencing.morel_montry_beta`,
which carries the distinction, and which this module is the gate for.

⛔ **Sphere only.** `[M]` on a σ_y-folded cylinder arc this β is identically
zero for BOTH the shipped τ and τ ≡ ½, at every ``n_phi`` — refuted there,
for the same antisymmetry reason that kills BMC's. The cylinder's
diffusion-limit evidence is the solve-based
``tests/sn/verification/analytical/test_angular_diffusion_limit_consistency.py``.

**What this CANNOT do.** It is a *constraint*, not a ranker: it answers "does
the scheme preserve the diffusion limit?" (pass/fail), never "which scheme is
more accurate". `[M]` the two questions genuinely diverge — on a fixture far
from the diffusion limit the plain diamond can be the more accurate scheme.
Do not cite this gate in an accuracy argument.
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.derivations.discrete.sn.angular_differencing import (
    morel_montry_beta,
)
from orpheus.geometry import CoordSystem
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.angular.closure import morel_montry_tau_per_level
from orpheus.sn.angular.redistribution import alpha_dome

pytestmark = pytest.mark.foundation


def _sphere_ascending(n: int):
    """Gauss–Legendre S_n with μ and w sorted by ascending μ (M&M's order)."""
    quad = Quadrature.gauss_legendre(n)
    order = np.argsort(quad.mu_x)
    tau = np.asarray(
        morel_montry_tau_per_level(quad, CoordSystem.SPHERICAL)[0], dtype=float
    )
    mu, w = quad.mu_x[order], quad.weights[order]
    # ``morel_montry_beta`` used to normalise w and build alpha itself; L0 may
    # not import orpheus.sn, so the caller does it now (P4.2). Same arithmetic.
    alpha = alpha_dome(mu, np.asarray(w, dtype=float) / np.sum(w))
    return mu, w, tau, alpha


class TestBetaIdentity:
    r"""β ≡ 0 for the shipped τ, β ≠ 0 for τ ≡ ½. No solve."""

    @pytest.mark.parametrize("n", [2, 4, 8, 16, 32])
    def test_shipped_tau_annihilates_beta(self, n: int) -> None:
        mu, w, tau, alpha = _sphere_ascending(n)
        beta = morel_montry_beta(mu, tau, alpha=alpha)
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
        """The NEGATIVE leg. Without it the row above is not shown to be
        LOADED on τ (``vv-principles`` #19: a positive reading alone cannot
        discriminate a loaded gate from a blind one)."""
        quad = Quadrature.gauss_legendre(n)
        order = np.argsort(quad.mu_x)
        mu, w = quad.mu_x[order], quad.weights[order]
        beta = morel_montry_beta(
            mu, np.full(n, 0.5),
            alpha=alpha_dome(mu, np.asarray(w, dtype=float) / np.sum(w)),
        )
        assert abs(beta) > beta_min, (
            f"GL-S{n}: plain angular diamond gives beta = {beta:.6e}, which "
            f"should exceed {beta_min:.0e}. beta falls ~4 orders per doubling "
            "of N, so a too-small value here means the instrument stopped "
            "reading tau, not that the diamond became consistent."
        )

    @pytest.mark.parametrize(
        "n, beta_min",
        [(2, 1e-1), (4, 1e-4), (8, 1e-6), (16, 1e-8), (32, 1e-10)],
    )
    def test_beta_catches_the_tau_reflection(
        self, n: int, beta_min: float,
    ) -> None:
        r"""⭐ ``τ → 1−τ`` is the Mode-12 stabiliser hole this gate exists for.

        The reflection preserves the ``[0,1]`` membership guard, preserves
        the fold box, and preserves ``Π(1−τ)/τ = 1`` **exactly** (it inverts
        every factor and the product telescopes to 1 either way), so the
        neutral-stability gate in ``test_psi_half_positivity.py`` is GREEN
        for it. β is not.

        `[M]` 2026-08-12, ``|β(1−τ)|`` = ``2.679e-1 / 1.162e-3 / 1.673e-5 /
        2.818e-7 / 4.698e-9`` at N = 2/4/8/16/32 — it falls ~4 orders per
        doubling, exactly like the diamond leg, so the floors are per-N.
        (A fixed floor was tried first and was a **false red at N=16**;
        that is why they are tabulated rather than guessed.) The shipped τ
        reads ``≤ 8.3e-17`` at every N, so the separation is ≥10 orders
        throughout. Curiosity worth its line: `[M]` ``β(1−τ) ≈
        β(diamond)/2`` at every N ≥ 4, to 3 figures.
        """
        mu, w, tau, alpha = _sphere_ascending(n)
        reflected = 1.0 - tau
        ratios = (1.0 - reflected) / reflected
        assert abs(float(np.prod(ratios)) - 1.0) < 1e-12, (
            "precondition: the reflection must preserve the end-to-end "
            "product, or this row is not testing the hole it claims to"
        )
        beta = morel_montry_beta(mu, reflected, alpha=alpha)
        assert abs(beta) > beta_min, (
            f"GL-S{n}: beta = {beta:.6e} did not catch tau -> 1-tau above the "
            f"{beta_min:.0e} floor, and every other shipped tau gate passes "
            "that reflection. This gate's whole justification is that it "
            "sees it."
        )
