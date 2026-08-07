r"""The 1-D cylindrical SN solution must be EVEN in the azimuthal cosine
:math:`\xi` — and the Morel--Montry :math:`\eta`-march is not.

Promoted 2026-08-01 from
``derivations/diagnostics/diag_326_azimuthal_mirror_symmetry.py``
(numerics-investigator, GitHub issue **#326**; the defect itself is issue
**#229**'s azimuthal floor, seen without a reference solution).

WHAT THIS FILE PROVES
---------------------
1-D cylindrical geometry is invariant under the mirror through the plane
spanned by ``z_hat`` and ``r_hat``.  That mirror maps the direction cosines

    (eta, xi, mu_z)  ->  (eta, -xi, mu_z)

and leaves the mesh, the cross sections, and any ``xi``-even source invariant.
Hence the exact angular flux satisfies

    psi(r, eta, xi, mu_z) == psi(r, eta, -xi, mu_z).

The ``product`` quadrature is closed under that mirror
(the :math:`\sigma_y` ordinate pairing is an involution pairing ``phi`` with
``2*pi - phi``) and the two partners carry EQUAL weights — so the SEMI-discrete
(discrete-ordinate, continuous-space) problem inherits the symmetry EXACTLY.  A
symmetry-respecting closure must therefore reproduce it to solver tolerance.

This is a **structurally independent** criterion in the ``vv-principles`` sense:
a closed-form statement about the operator's own symmetry group.  No reference
solver, no manufactured source, no second implementation.  It is also the ONE
criterion available that the curvilinear MMS cannot see — both production
cylindrical MMS ansatzes are functions of ``eta`` and ``xi**2`` only, i.e. they
live entirely INSIDE the symmetric sector (proved in
``tests/sn/verification/mms/test_mms_ordering_blindness.py``).

WHAT THIS FILE IS BLIND TO
--------------------------
* It does NOT rank the candidate per-level orderings.  The defect MAGNITUDE is
  ordering-INVARIANT (measured identical to 1e-16 across three tie-breaks);
  only the LABEL moves.  That is the finding, not a limitation: there is no
  correct ordering to choose.
* It says nothing about the SCALAR flux, which is very nearly invariant here —
  the two mirror partners' errors largely cancel in ``sum_n w_n psi_n``.  The
  claim is per-ordinate.

MEASURED 2026-08-01 (branch ``refactor/operator-strategy-layers``)
------------------------------------------------------------------
``max_n |psi_n - psi_{mirror_y(n)}| / max|psi|`` on a cylindrical fixed-source
problem with an ISOTROPIC (hence ``xi``-even) source::

    product(n_mu=4, n_phi=8),  nx=20, homogeneous    ->  1.19e-1  (30 % local)
    product(n_mu=4, n_phi=8),  nx=20, heterogeneous  ->  5.14e-2
    level_symmetric(4)                               ->  3.08e-1

MECHANISM.  ``morel_montry_tau_raw_per_level`` builds ``eta_edge`` as midpoints
of CONSECUTIVE ordinates in the level order.  Two mirror partners share ``eta``
bit-exactly, so their shared edge collapses ONTO the node: the lower partner
gets ``tau_raw = 1`` and the upper ``tau_raw = 0``, which the structural
``[1/2, 1]`` clamp turns into ``tau = 1`` and ``tau = 1/2``.  Two ordinates the
geometry says are identical receive DIFFERENT angular weights.

WHY NO RE-ORDERING FIXES IT.  Any ordering ascending in ``eta`` puts the
partners adjacent, hence forces the ``{1, 1/2}`` split; any ordering that is
NOT ascending in ``eta`` breaks the alpha-dome and the ``eta_edge``
monotonicity the rest of the closure assumes (the azimuthal ordering NaNs — see
``test_alpha_closed_form.py``).  This is a CLOSURE defect, not an ordering
choice.  The constructive exit is the half-range level (Hebert §3.9.3), under
which only one member of each pair exists and the symmetry holds by
construction.

DELIBERATE RED GATES — and how the incidental-failure hole is closed
--------------------------------------------------------------------
The three ``xfail(strict=True)`` rows below document the LIVE defect; the
XPASS-failure is what forces their deletion when the closure is repaired.  Each
body is exactly ONE assertion on a float computed by a FIXTURE, so nothing but
the documented inequality can fail *in the call phase* — the ``vv-principles``
Mode-8 class-4 (misattributed strict-xfail) discipline.

**MEASURED 2026-08-01, and it corrects a plausible assumption:** moving the
solve into a fixture is NOT by itself sufficient.  pytest's ``xfail`` swallows
**fixture setup errors too** — an in-process plugin replacing
``_solve_isotropic_source`` with a raising stub left all three rows reporting
``3 xfailed`` in 0.32 s (the solve never ran) instead of erroring.  So the
incidental-failure hole is closed by a reddenable SIBLING instead:
:func:`test_mirror_defect_fixtures_are_well_formed` consumes the same three
fixtures WITHOUT an xfail, so any construction/solve breakage reds there,
loudly, while the xfail rows keep asserting only their documented inequality.

Both directions are control-verified:

* break the solve  -> the un-xfailed sibling ERRORs (the breakage is loud);
* simulate the fix -> all three rows flip to ``XPASS(strict)`` -> FAILED,
  proving the gates measure exactly what the repair will change.

No row carries ``verifies(...)``: the equation-level property is currently
BROKEN, and a ``verifies`` edge on a red gate is a phantom coverage claim.
Add it when the xfails flip.
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.derivations.continuous.mms.sn import _make_1g_mixture
from orpheus.geometry import BC, CoordSystem, Mesh1D
from orpheus.geometry.boundary import SelfPairedDeck
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn import solve_sn_fixed_source


def _mirror_pairing(quad: Quadrature, axis: str) -> np.ndarray:
    """The ordinate pairing :math:`\\sigma_{axis}` induces — read from
    production's own source (``ordinate_permutation``, the same derivation
    the coupled-pole seed and the BC realizer consume). The control legs
    below verify its ξ-mirror facts INDEPENDENTLY (η/μ_z held, ξ negated,
    weights equal), so the defect fixtures do not merely trust it."""
    pi = quad.ordinate_permutation(
        SelfPairedDeck.mirror(axis=axis, dimension=3).motion
    )
    if pi is None:
        raise AssertionError(f"no {axis}-mirror ordinate pairing on this rule")
    return pi.indices

_XFAIL_326 = pytest.mark.xfail(
    strict=True,
    reason=(
        "ISSUE #326/#229 — the cylindrical Morel-Montry eta-march breaks the "
        "geometry's xi -> -xi mirror symmetry at O(1e-1). Delete this marker "
        "when the azimuthal closure is repaired."
    ),
)


def _solve_isotropic_source(nx: int, quad: Quadrature, *, het: bool):
    """Cylindrical fixed source, ISOTROPIC (``xi``-even) source, ``xi``-even BCs."""
    mat_ids = (
        np.array([0] * (nx // 2) + [1] * (nx - nx // 2), dtype=int)
        if het else np.zeros(nx, dtype=int)
    )
    mesh = Mesh1D(
        edges=np.linspace(0.0, 2.0, nx + 1),
        mat_ids=mat_ids,
        coord=CoordSystem.CYLINDRICAL,
        bc_left=BC("reflective"),          # r = 0 symmetry axis
        bc_right=BC("vacuum"),
    )
    materials = {0: _make_1g_mixture(1.0, 0.5)}
    if het:
        materials[1] = _make_1g_mixture(4.0, 3.6)
    Q = np.full((quad.N, 1, nx), 1.0 / float(quad.weights.sum()))
    result = solve_sn_fixed_source(
        materials, mesh, quad, Q, max_inner=2000, inner_tol=1e-13,
    )
    return np.asarray(result.angular_flux.interior.values, dtype=float)


def _mirror_defect(quad: Quadrature, psi: np.ndarray) -> float:
    partner = _mirror_pairing(quad, "y")
    return float(np.max(np.abs(psi - psi[partner]))) / max(
        float(np.max(np.abs(psi))), 1e-300
    )


# ── fixtures: the solve lives OUTSIDE the xfail'd call phase ────────────

@pytest.fixture(scope="module")
def defect_homogeneous() -> float:
    quad = Quadrature.product(n_mu=4, n_phi=8)
    return _mirror_defect(quad, _solve_isotropic_source(20, quad, het=False))


@pytest.fixture(scope="module")
def defect_heterogeneous() -> float:
    quad = Quadrature.product(n_mu=4, n_phi=8)
    return _mirror_defect(quad, _solve_isotropic_source(20, quad, het=True))


@pytest.fixture(scope="module")
def defect_level_symmetric() -> float:
    quad = Quadrature.level_symmetric(sn_order=4)
    return _mirror_defect(quad, _solve_isotropic_source(20, quad, het=False))


# ── the invariant's PRECONDITION (control legs — MUST stay green) ───────

@pytest.mark.foundation
@pytest.mark.parametrize("n_phi", [4, 8, 16])
def test_y_mirror_pairing_is_the_xi_mirror_involution(n_phi):
    """The σ_y pairing really pairs ``(eta, xi)`` with ``(eta, -xi)``.

    Control leg.  Without it the defect measurement below would be comparing
    the wrong ordinates and could manufacture a false positive.
    """
    quad = Quadrature.product(n_mu=4, n_phi=n_phi)
    partner = _mirror_pairing(quad, "y")
    np.testing.assert_array_equal(partner[partner], np.arange(quad.N))
    np.testing.assert_allclose(quad.eta[partner], quad.eta, atol=1e-15)
    np.testing.assert_allclose(quad.mu_z[partner], quad.mu_z, atol=1e-15)
    np.testing.assert_allclose(quad.xi[partner], -quad.xi, atol=1e-15)
    np.testing.assert_array_equal(quad.weights[partner], quad.weights)


@pytest.mark.foundation
def test_slab_quadrature_has_a_trivial_xi_mirror_control():
    """Control leg: on a 1-D slab GL rule the ``xi`` mirror is the identity.

    Pins that the defect below is specific to the cylindrical AZIMUTHAL
    structure, not an artefact of how the pairing is derived.
    """
    quad = Quadrature.gauss_legendre(8)
    np.testing.assert_array_equal(
        _mirror_pairing(quad, "y"), np.arange(quad.N)
    )


# ── the reddenable sibling that closes the xfail's incidental-error hole ─

@pytest.mark.foundation
def test_mirror_defect_fixtures_are_well_formed(
    defect_homogeneous, defect_heterogeneous, defect_level_symmetric
):
    """Every solve the xfail rows below depend on actually ran and converged.

    This row is deliberately NOT xfailed.  Measured: pytest's ``xfail``
    swallows fixture SETUP errors as well as call-phase failures, so without
    this sibling a broken solve would report three quiet ``xfailed`` rows and
    look like committed coverage of the #326 defect while asserting nothing
    about it.  Here a broken solve reds loudly instead.

    The band is deliberately wide — this is an activation check, not the
    defect measurement.  ``> 1e-13`` also pins that the defect is REACHABLE:
    if it ever collapses to solver noise the xfail rows have stopped measuring
    anything and must be re-adjudicated (they would XPASS anyway).
    """
    for name, value in (
        ("homogeneous", defect_homogeneous),
        ("heterogeneous", defect_heterogeneous),
        ("level_symmetric", defect_level_symmetric),
    ):
        assert np.isfinite(value), f"{name}: mirror defect is not finite ({value})"
        assert 1e-13 < value <= 2.0, (
            f"{name}: mirror defect {value:.6e} outside the plausible band "
            f"(1e-13, 2.0] — the solve did not produce a usable angular flux, "
            f"or the defect vanished and the xfail rows below are stale"
        )


# ── the defect (deliberate red gates) ───────────────────────────────────

@_XFAIL_326
@pytest.mark.l1
def test_cylindrical_solution_is_even_in_xi_homogeneous(defect_homogeneous):
    """The discrete solution must inherit the geometry's ``xi``-mirror symmetry."""
    assert defect_homogeneous < 1e-11, (
        f"1-D cylindrical angular flux is NOT even in xi: relative mirror "
        f"defect {defect_homogeneous:.6e} (isotropic source, homogeneous)."
    )


@_XFAIL_326
@pytest.mark.l1
def test_cylindrical_solution_is_even_in_xi_heterogeneous(defect_heterogeneous):
    """Same invariant with a two-region radial material map."""
    assert defect_heterogeneous < 1e-11, (
        f"1-D cylindrical angular flux is NOT even in xi: relative mirror "
        f"defect {defect_heterogeneous:.6e} (isotropic source, heterogeneous)."
    )


@_XFAIL_326
@pytest.mark.l1
def test_level_symmetric_cylinder_is_even_in_xi(defect_level_symmetric):
    """Same invariant on ``level_symmetric``, whose ``eta`` degeneracy is 4-to-1.

    Measured WORSE than the product rule (3.08e-1 vs 1.19e-1) — the degeneracy
    MULTIPLICITY, not the rule family, sets the severity.
    """
    assert defect_level_symmetric < 1e-11, (
        f"level_symmetric(4) cylindrical flux is NOT even in xi: "
        f"{defect_level_symmetric:.6e}"
    )


# ── characterisation (GREEN — these pin what the defect IS) ─────────────

@pytest.mark.l1
@pytest.mark.slow
def test_mirror_defect_is_an_angular_not_a_spatial_defect():
    """It does NOT vanish under mesh refinement — so it is the closure.

    Measured 2026-08-01, ``product(4, 8)``, homogeneous:
    ``nx 10/20/40/80 -> 1.127e-1, 1.191e-1, 1.252e-1, 1.287e-1`` (it GROWS
    slightly).  A spatial artefact would fall as O(h) or better.
    """
    quad = Quadrature.product(n_mu=4, n_phi=8)
    defect = {
        nx: _mirror_defect(quad, _solve_isotropic_source(nx, quad, het=False))
        for nx in (10, 20, 40, 80)
    }
    print(f"mirror defect vs nx: {defect}")
    ratios = [defect[a] / defect[b] for a, b in ((10, 20), (20, 40), (40, 80))]
    assert all(r < 2.0 for r in ratios), (
        f"mirror defect fell at least O(h) under refinement {defect} "
        f"(ratios {ratios}) — it would then be a spatial artefact, not the "
        f"angular-closure defect"
    )


@pytest.mark.l1
@pytest.mark.slow
def test_mirror_defect_scales_with_the_AZIMUTHAL_order_only():
    """The defect is the #229 azimuthal-thread floor wearing a new face.

    It falls with ``n_phi`` and is FLAT in ``n_mu`` — exactly the axis
    ``test_cyl_aniso_floor_scales_with_quadrature`` already identified for the
    anisotropic-MMS floor.  Measured 2026-08-01 (nx=20):
    ``n_phi 4/8/16/32/64 -> 1.341e-1, 1.191e-1, 9.665e-2, 6.119e-2, 3.373e-2``;
    ``n_mu  2/4/8/16     -> 1.251e-1, 1.191e-1, 1.205e-1, 1.208e-1``.
    """
    by_phi, by_mu = {}, {}
    for n_phi in (8, 64):
        quad = Quadrature.product(n_mu=4, n_phi=n_phi)
        by_phi[n_phi] = _mirror_defect(
            quad, _solve_isotropic_source(20, quad, het=False)
        )
    for n_mu in (2, 16):
        quad = Quadrature.product(n_mu=n_mu, n_phi=8)
        by_mu[n_mu] = _mirror_defect(
            quad, _solve_isotropic_source(20, quad, het=False)
        )
    print(f"defect by n_phi={by_phi}   by n_mu={by_mu}")
    assert by_phi[64] < by_phi[8] / 3.0, (
        f"defect did not fall with the AZIMUTHAL order: {by_phi}"
    )
    assert 0.8 < by_mu[16] / by_mu[2] < 1.25, (
        f"defect should be flat in the POLAR order: {by_mu}"
    )


# ── the same unenforced symmetry, a second time: the pole map ───────────

@pytest.mark.foundation
def test_cylinder_pole_map_and_axis_crossing_differ_by_exactly_the_xi_mirror():
    r"""The r=0 continuity map is ``(eta,xi) -> (-eta,+xi)``; the straight-line
    axis crossing is ``(eta,xi) -> (-eta,-xi)``.  They agree IFF psi is xi-even.

    ``orpheus/sn/loss_representation/__init__.py`` seeds each outward
    ordinate's r=0 inflow at the x-mirror partner
    (``pole_outflow[mirror[n]]``, the pairing ``_ensure_pole_mirror``
    derives) — the map ``omega -> pi - omega``.  A ray that actually crosses
    the axis keeps its LAB direction while the local ``(r_hat, phi_hat)``
    frame rotates by ``pi``, so BOTH components flip: ``omega -> omega + pi``.

    The two maps compose to exactly the ``sigma_y`` pairing — the ``xi``
    mirror.  So the shipped pole map is correct up to the ``xi``-symmetry the
    solution is ASSUMED to have and (per the xfailed rows above) measurably
    does not.  FLAGGED as a structural observation, NOT a confirmed bug: it
    needs its own adjudication.
    """
    n_phi = 8
    quad = Quadrature.product(n_mu=2, n_phi=n_phi)
    mirror_eta = _mirror_pairing(quad, "x")
    mirror_xi = _mirror_pairing(quad, "y")
    rotate_pi = np.array([                        # omega -> omega + pi
        (n // n_phi) * n_phi + ((n % n_phi) + n_phi // 2) % n_phi
        for n in range(quad.N)
    ])
    np.testing.assert_allclose(quad.eta[mirror_eta], -quad.eta, atol=1e-15)
    np.testing.assert_allclose(quad.xi[mirror_eta], +quad.xi, atol=1e-15)
    np.testing.assert_allclose(quad.eta[rotate_pi], -quad.eta, atol=1e-15)
    np.testing.assert_allclose(quad.xi[rotate_pi], -quad.xi, atol=1e-15)
    np.testing.assert_array_equal(mirror_eta, mirror_xi[rotate_pi])


@pytest.mark.foundation
def test_tie_break_permutation_does_not_commute_with_the_pole_map():
    """WHY the tie-break leaks into the answer even on ``xi``-EVEN data.

    A tie-break swap ``sigma`` acts inside ONE ``eta`` class.  The pole seed
    couples ordinate ``n`` to its x-mirror partner, which lives in the
    ``-eta`` class — and nothing forces the tie-break to have made the same
    choice there.  So ``sigma`` and the pole map do NOT commute, and the pole
    seed feeds an ordinate from the wrong ``xi`` branch.

    This is the mechanism behind issue #326's own table: an ISOTROPIC
    (``xi``-even) source still moves 2.6e-2 in psi / 6.6e-4 in phi on a
    heterogeneous cylinder, even though the tie-break permutes only ordinates
    carrying identical data.  Without the pole coupling it would be an exact
    relabeling and phi would be bit-invariant.
    """
    quad = Quadrature.product(n_mu=2, n_phi=8)
    mirror_eta = _mirror_pairing(quad, "x")
    mirror_xi = _mirror_pairing(quad, "y")

    # A tie-break swap acting on ONE eta class only: swap the (+eta, +-xi)
    # partners and leave the (-eta, +-xi) partners alone.
    sigma = np.arange(quad.N)
    level = np.asarray(quad.level_indices[0])
    eta_level = quad.eta[level]
    for m in range(eta_level.size - 1):
        if abs(eta_level[m] - eta_level[m + 1]) < 1e-15 and eta_level[m] > 0.1:
            a, b = int(level[m]), int(level[m + 1])
            sigma[a], sigma[b] = b, a
            break
    assert not np.array_equal(sigma, np.arange(quad.N)), (
        "no one-eta-class swap could be constructed — the level no longer has "
        "an exact positive-eta tie, so this row is vacuous"
    )
    print(f"sigma swaps {np.flatnonzero(sigma != np.arange(quad.N))}")

    assert not np.array_equal(mirror_eta[sigma], sigma[mirror_eta]), (
        "a one-eta-class tie-break swap unexpectedly COMMUTES with the pole "
        "map — the leak mechanism would then be something else"
    )
    # ... whereas it always commutes with the xi mirror (both act within a
    # xi-pair), which is why phi IS invariant once the pole map is xi-even.
    np.testing.assert_array_equal(mirror_xi[sigma], sigma[mirror_xi])
