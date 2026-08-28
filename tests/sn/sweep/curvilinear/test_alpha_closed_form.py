r"""The cylindrical :math:`\alpha`-recursion has a CLOSED FORM: it is the
azimuthal cosine :math:`\xi` sampled at the half-angle boundaries.

Promoted 2026-08-01 from ``derivations/diagnostics/diag_326_alpha_closed_form.py``
(numerics-investigator, GitHub issue **#326** — "per-level ordinate ordering is
under-determined").

WHAT THIS FILE PROVES
---------------------
:eq:`alpha-cylindrical` is stated in the corpus only as a *recursion*,
``alpha_{p,m+1/2} = alpha_{p,m-1/2} - w_m eta_m``.  It also has a closed form,
and that closed form is the literature's own **definition** of alpha (Hebert
2009 Eqs. 3.398/3.399 define ``alpha_{p,q+-1/2} = W_p * eta_tangential`` at REAL
half-angle boundaries in the azimuthal angle ``omega``).  Two structurally
independent grounds agree on it: the literature derivation, and the
Dirichlet-kernel identity measured here to 3e-16.

For the ``product`` rule the per-level weight is constant
(``w = w_gl(mu_p) * d_omega``, ``rules_product.py``), so

    alpha_k = -w_gl * d_omega * sin(theta) * sum_{j<k} cos(omega_j)

and the Dirichlet kernel closes that sum EXACTLY:

    alpha_k = -w_gl * kappa * [ xi(omega_{k-1/2}) - xi(omega_{-1/2}) ],
    kappa   = d_omega / (2 sin(d_omega/2)) = 1 + d_omega**2/24 + ...

The rows that carry ``verifies("alpha-cylindrical")`` drive the PRODUCTION
producer :func:`~orpheus.sn.mesh.reduced_operator.cylindrical_streaming`
(``reduced_operator.py``, the ``alpha[m+1] = alpha[m] - w[m]*eta[m]`` loop) —
they are not a local re-implementation of it.  The ordering is supplied as an
INPUT via :func:`~tests.sn._test_helpers.product_level_ordering`.

WHAT THIS FILE IS BLIND TO
--------------------------
* It says NOTHING about which per-level ordering is correct in production.  The
  closed form holds under the AZIMUTHAL ordering and fails by O(1) under the
  production ascending-``eta`` ordering — but the azimuthal ordering cannot be
  used (it breaks the alpha-dome positivity every downstream consumer assumes,
  and the cylindrical solve NaNs).  That adjudication lives in
  ``test_azimuthal_mirror_symmetry.py``.
* The half-range rows are a DERIVATION of the constructive target (Hebert
  §3.9.3's ``0 <= omega <= pi`` level).  Since Q5.5/Q5.6 production DOES run
  that configuration (``Quadrature.folded_product`` — and `[M]` 2026-08-08
  the production folded α closes exactly, α[0] = α[M] = 0 machine-zero on
  every level, with each level MONOTONE in ω so the 2-to-1 η collision that
  breaks the closed form under the full-circle production ordering cannot
  arise).  The closed-FORM identity on the production folded producer — the
  arc's own staggered Dirichlet-kernel normalization — is the Q5.6.4 "T3 α
  gate" deliverable; until it lands, these rows stay ``foundation`` and the
  full-circle rows above remain the record of the ordering pathology the
  fold dissolved.

THE ``kappa`` TRAP
------------------
``alpha == -xi`` is FALSE at production quadrature orders: ``kappa`` is 2.6 %
off unity at ``n_phi = 8``.  Every closed-form assertion below carries
``kappa`` explicitly; a bare ``alpha == -xi`` gate would be asserting something
untrue.
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.geometry import CoordSystem
from orpheus.sn.mesh.reduced_operator import cylindrical_streaming
from orpheus.numerics.quadrature import Quadrature

from tests.sn._test_helpers import (
    curvilinear_homogeneous_mesh,
    product_level_ordering,
)


# ── production alpha, with the level ordering as a controllable input ───

def _production_alpha_per_level(n_mu: int, n_phi: int, *, ordering: str | None):
    """Return ``(alpha_per_level, quad)`` from the PRODUCTION factory.

    ``ordering=None`` means production's own ordering (no patch).  Anything
    else is one of :data:`~tests.sn._test_helpers.PRODUCT_LEVEL_ORDERINGS`.
    """
    mesh = curvilinear_homogeneous_mesh(
        4, 2.0, coord=CoordSystem.CYLINDRICAL
    )
    if ordering is None:
        quad = Quadrature.product(n_mu=n_mu, n_phi=n_phi)
    else:
        with product_level_ordering(ordering):
            quad = Quadrature.product(n_mu=n_mu, n_phi=n_phi)
    reduced = cylindrical_streaming(mesh, quad)
    return list(reduced.angular.alpha_per_level), quad


def _level_geometry(quad: Quadrature, p: int):
    """``(eta, xi, w_gl, sin_theta)`` for level ``p``, in level order."""
    level = np.asarray(quad.level_indices[p])
    sin_theta = float(np.sqrt(1.0 - quad.mu_z[level[0]] ** 2))
    n_phi = level.size
    w_gl = float(quad.weights[level[0]]) / (2.0 * np.pi / n_phi)
    return quad.eta[level], quad.xi[level], w_gl, sin_theta


def _kappa(n_phi: int) -> float:
    """The Dirichlet-kernel prefactor of the discrete cumulative sum."""
    d_omega = 2.0 * np.pi / n_phi
    return float(d_omega / (2.0 * np.sin(d_omega / 2.0)))


def _alpha_closed_form(n_phi: int, w_gl: float, sin_theta: float) -> np.ndarray:
    r"""``alpha_k = -w_gl * kappa * [xi(omega_{k-1/2}) - xi(omega_{-1/2})]``.

    Half-angle boundary ``k`` sits at ``omega_{k-1/2} = (k - 1/2) d_omega`` —
    the midpoint between node ``k-1`` and node ``k``.
    """
    d_omega = 2.0 * np.pi / n_phi
    k = np.arange(n_phi + 1)
    xi_half = sin_theta * np.sin((k - 0.5) * d_omega)
    xi_start = sin_theta * np.sin(-0.5 * d_omega)
    return -w_gl * _kappa(n_phi) * (xi_half - xi_start)


# ── A. the closed form, on the PRODUCTION recursion ─────────────────────

@pytest.mark.l1
@pytest.mark.parametrize("n_phi", [4, 8, 16, 32])
@pytest.mark.parametrize("n_mu", [2, 4])
@pytest.mark.verifies("alpha-cylindrical")
def test_production_alpha_equals_xi_at_the_half_angle_boundaries(n_mu, n_phi):
    r"""``alpha == -w_gl * kappa * (xi - xi_start)`` EXACTLY, every level.

    The structurally-independent ground is Hebert 2009 Eq. 3.399 (alpha IS the
    tangential cosine at a real half-angle boundary); the numerical identity is
    the Dirichlet kernel.  Driven through the production
    :func:`cylindrical_streaming` recursion with the AZIMUTHAL ordering, under
    which the level's half-angle boundaries really are contiguous in ``omega``.
    """
    alphas, quad = _production_alpha_per_level(
        n_mu, n_phi, ordering="azimuthal"
    )
    for p, alpha in enumerate(alphas):
        _eta, _xi, w_gl, sin_theta = _level_geometry(quad, p)
        want = _alpha_closed_form(n_phi, w_gl, sin_theta)
        err = float(np.max(np.abs(alpha - want)))
        scale = max(float(np.max(np.abs(want))), 1e-300)
        assert err / scale < 1e-13, (
            f"n_mu={n_mu} n_phi={n_phi} level={p}: the production alpha "
            f"recursion departs from the closed form by {err:.3e} "
            f"(rel {err / scale:.3e})\n got ={alpha}\n want={want}"
        )


@pytest.mark.foundation
def test_kappa_prefactor_is_not_unity_and_converges_to_one():
    """The closed form carries a genuine O(d_omega^2) prefactor.

    ``alpha == -xi`` is exact only up to ``kappa``; at the production
    ``n_phi = 8`` the correction is 2.6 %, so a gate asserting a bare
    ``alpha == -xi`` would be asserting something FALSE.  This row exists so
    that trap cannot be re-introduced silently.
    """
    got = {n_phi: _kappa(n_phi) for n_phi in (4, 8, 16, 32, 64)}
    # kappa = 1 + d_omega**2 / 24 + O(d_omega**4).
    assert got[8] > 1.02, f"kappa(8) = {got[8]}"
    assert 1.0003 < got[64] < 1.0005, f"kappa(64) = {got[64]}"
    ratios = [
        (got[a] - 1.0) / (got[b] - 1.0)
        for a, b in ((4, 8), (8, 16), (16, 32), (32, 64))
    ]
    assert all(3.7 < r < 4.3 for r in ratios), (
        f"kappa-1 should fall as O(d_omega^2) (ratio 4 per doubling); "
        f"got {ratios} from {got}"
    )


# ── B. what the PRODUCTION ordering gives instead ───────────────────────

@pytest.mark.foundation
@pytest.mark.parametrize("n_phi", [8, 16, 32])
def test_alpha_is_NOT_the_omega_closed_form_under_production_ordering(n_phi):
    """Production sorts by ``eta``, which is 2-to-1 in ``omega``.

    The closed form is then wrong by O(1), and the discrepancy does NOT shrink
    with ``n_phi`` — it is structural (the level's "half-angle boundaries" are
    not contiguous in ``omega`` once the sort interleaves the two
    half-circles), not a quadrature error.  Characterisation only: no
    ``verifies`` edge, because a negative result verifies no equation.
    """
    alphas, quad = _production_alpha_per_level(4, n_phi, ordering=None)
    rels = []
    for p, alpha in enumerate(alphas):
        _eta, _xi, w_gl, sin_theta = _level_geometry(quad, p)
        want = _alpha_closed_form(n_phi, w_gl, sin_theta)
        scale = max(float(np.max(np.abs(want))), 1e-300)
        rels.append(float(np.max(np.abs(alpha - want))) / scale)
    assert min(rels) > 0.5, (
        f"n_phi={n_phi}: the eta-sorted alpha unexpectedly MATCHES the omega "
        f"closed form (rel {rels}) — the ordering premise changed"
    )


@pytest.mark.l1
@pytest.mark.parametrize("n_phi", [4, 8, 16, 32])
@pytest.mark.parametrize("n_mu", [2, 4, 8])
@pytest.mark.verifies("alpha-cylindrical", "alpha-recursion")
def test_production_alpha_is_a_non_negative_closing_dome(n_mu, n_phi):
    r"""Ascending-``eta`` ordering <=> ``alpha >= 0`` with ``alpha[0] == alpha[M] == 0``.

    This is the real content of the sort, and the property
    ``curvilinear_one_group.rst`` asserts for :eq:`alpha-cylindrical`: the
    prefix sums of an ascending-sorted zero-sum sequence are non-positive, so
    ``alpha = -w * cumsum(eta)`` is a non-negative dome closing at both ends.
    It is an if-and-only-if — the azimuthal companion below goes negative.
    """
    alphas, _quad = _production_alpha_per_level(n_mu, n_phi, ordering=None)
    for p, alpha in enumerate(alphas):
        assert alpha.min() > -1e-14, f"level {p}: alpha dome negative {alpha}"
        assert abs(alpha[0]) < 1e-14 and abs(alpha[-1]) < 1e-13, (
            f"level {p}: dome does not close, alpha={alpha}"
        )


@pytest.mark.foundation
@pytest.mark.parametrize("n_phi", [8, 16, 32])
def test_azimuthal_ordering_breaks_the_alpha_dome_positivity(n_phi):
    """The omega ordering makes ``alpha`` change sign — so it is NOT a drop-in.

    Downstream, ``morel_montry_tau_per_level`` builds ``eta_edge`` as
    midpoints of CONSECUTIVE ordinates bracketed by ``-sin_theta`` and
    ``+sin_theta``, and ``cylindrical_streaming`` publishes
    ``mu_start_per_level = -sin_theta``; both presume a monotone ascending
    ``eta`` march.  The azimuthal ordering violates that premise, which is why
    the closed-form criterion cannot simply be adopted.
    """
    alphas, _quad = _production_alpha_per_level(
        4, n_phi, ordering="azimuthal"
    )
    for p, alpha in enumerate(alphas):
        assert alpha.min() < -1e-3, (
            f"level {p}: expected the omega ordering to break dome "
            f"positivity, got alpha={alpha}"
        )


# ── C. on the HALF range the two criteria COINCIDE ──────────────────────
#
# psi is even in xi (the 1-D cylindrical mirror symmetry — see
# ``test_azimuthal_mirror_symmetry.py``), so only ``omega in [0, pi]`` is
# independent.  There ``eta = sin(theta) cos(omega)`` is strictly MONOTONE:
# ascending-eta and descending-omega are the SAME ordering, there are no ties,
# and alpha is BOTH a non-negative dome AND the closed form with NO additive
# constant (the march starts at omega = pi, where xi = 0).  The whole #326
# degeneracy is an artefact of carrying the redundant half of the circle.
#
# These rows are ``foundation``: they pin the DERIVATION of the constructive
# target (Hebert §3.9.3), not a production path — production does not run a
# half-range level, so they carry no ``verifies`` edge.

def _half_range_level(n_mu: int, n_phi: int, p: int):
    """Fold the production product level onto ``omega in [0, pi]``.

    Nodes ``omega_j = j d_omega``, ``j = 0..n_phi/2`` (endpoints INCLUDED, as
    the full rule has them); weight ``2w`` for the interior partners that were
    folded together, ``w`` at the two self-paired endpoints.  Total level
    weight is preserved exactly.
    """
    mu_gl, w_gl = np.polynomial.legendre.leggauss(n_mu)
    sin_theta = float(np.sqrt(1.0 - mu_gl[p] ** 2))
    d_omega = 2.0 * np.pi / n_phi
    j = np.arange(n_phi // 2 + 1)
    omega = j * d_omega
    eta = sin_theta * np.cos(omega)
    w = np.where((j == 0) | (j == n_phi // 2), 1.0, 2.0) * (w_gl[p] * d_omega)
    order = np.argsort(eta)              # == descending omega; MONOTONE
    return eta[order], w[order], omega[order], float(w_gl[p]), sin_theta


def _alpha_recursion(eta: np.ndarray, w: np.ndarray) -> np.ndarray:
    """:eq:`alpha-cylindrical`, applied to a hand-built level."""
    alpha = np.zeros(eta.size + 1)
    for m in range(eta.size):
        alpha[m + 1] = alpha[m] - w[m] * eta[m]
    return alpha


@pytest.mark.foundation
@pytest.mark.parametrize("n_phi", [8, 16, 32, 64])
def test_half_range_march_has_no_ties_and_preserves_total_weight(n_phi):
    """The degeneracy #326 adjudicates does not exist on the half range."""
    for p in range(4):
        eta, w, _omega, w_gl, _sin_theta = _half_range_level(4, n_phi, p)
        assert np.min(np.abs(np.diff(eta))) > 1e-6, (
            f"half-range level {p} still has a near-tie: eta={eta}"
        )
        np.testing.assert_allclose(w.sum(), w_gl * 2.0 * np.pi, rtol=1e-14)


@pytest.mark.foundation
@pytest.mark.parametrize("n_phi", [8, 16, 32, 64])
def test_half_range_alpha_is_a_dome_AND_equals_xi(n_phi):
    """Both criteria hold simultaneously, and the additive constant is zero.

    Measured 2026-08-01: ``alpha / (w_gl * xi)`` is CONSTANT across every
    interior half-angle boundary and equals ``2 * kappa_half`` -> 2.0523
    (n_phi=8), 2.0129 (16), 2.0032 (32) — converging to the folded weight
    factor 2.
    """
    for p in range(4):
        eta, w, omega, w_gl, sin_theta = _half_range_level(4, n_phi, p)
        alpha = _alpha_recursion(eta, w)
        assert alpha.min() > -1e-14, f"not a dome: {alpha}"
        assert abs(alpha[0]) < 1e-14 and abs(alpha[-1]) < 1e-13, (
            f"dome does not close: {alpha}"
        )
        omega_half = 0.5 * (omega[:-1] + omega[1:])   # half-angle boundaries
        xi_half = sin_theta * np.sin(omega_half)
        ratio = alpha[1:-1] / (w_gl * xi_half)
        assert np.ptp(ratio) < 1e-12, (
            f"alpha is not proportional to xi at the half-angle boundaries: "
            f"ratio={ratio}"
        )
        assert abs(float(ratio[0]) - 2.0) < 0.06, (
            f"alpha/(w_gl*xi) should approach the folded weight factor 2, "
            f"got {ratio[0]} at n_phi={n_phi}"
        )
