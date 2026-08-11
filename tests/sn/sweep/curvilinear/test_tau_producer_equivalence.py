r"""Issue #236 Phase 2 — τ producer-equivalence gate (Leg 1).

The :class:`~orpheus.sn.sweep.pole_angular_closure.MorelMontryAngularSweep`
closure now PRODUCES the Morel–Montry angular weight τ from the quadrature
``(μ, w, levels)`` it already binds (an angular-scheme property), instead of
reading it back from the streaming-GEOMETRY factory (``reduced_operator.py``).

Step C (the geometry-τ retirement) deleted the geometry-side τ producer, so
the former ``*_equals_geometry_factory_0ulp`` legs (closure τ vs
``spherical_streaming(...).tau_mm`` / ``cylindrical_streaming(...).tau_mm_per_level``)
are now vacuous and have been removed.  The independence floor is carried
entirely by the surviving legs below, and **the two arms now use two
DIFFERENT references** — a distinction that must not be flattened again:

⛔ ``angular_differencing.morel_montry_weights`` **DELEGATES to the
production producer** (Q5.6.4) — deliberately, so a "reference" can never
drift into a second definition of the angular cell, which is exactly how
its cylinder arm had gone wrong.  It is therefore **no longer an
independent reference for τ at all**, and both arms below now use
hand-authored references instead:

* **SPHERE** → a hand-written cumulative-weight expression, inline.  BMC
  2010 Eq. 12/42: ``mu_edge[n+1] = mu_edge[n] + w[n]`` from ``-1``, then
  P2.  Authored here, independent of the producer.
* **CYLINDER** → the **analytic closed form**
  :math:`\tau_m = \tfrac12 + \tfrac12\cot\omega_m\tan(\Delta\omega/4)`,
  hand-derived from the arc geometry.  STRUCTURAL independence — it shares
  no code path with the producer, a strictly stronger footing than the
  procedural twin it replaces.

This gate proves the closure-produced τ:

* equals a structurally- (cylinder) or procedurally- (sphere) INDEPENDENT
  reference;
* on the cylinder, DIFFERS from the retired chord convention — a negative
  control (`vv-principles` #19), without which the row would pass equally
  against the partition the Q5.6.4 carve replaced and so could not be
  evidence about the partition choice;
* on Cartesian, is the neutral τ = 1.0 produced WITHOUT a geometry branch
  (the closure TYPE is the dispatch — Pattern 4).

Mode-8 discipline: every numeric assertion is ``np.testing`` /
``np.array_equal`` / ``pytest.fail`` so it FIRES under the canonical ``-O``
invocation.  The two guard-the-guard predicates use bare ``assert`` (pytest
rewrites those in test modules so they fire regardless of ``-O``); the VALUE
checks never rely on a bare assert.
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.geometry import CoordSystem, Mesh1D
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.mesh.augmented_mesh import SNMesh
from orpheus.sn.sweep.pole_angular_closure import (
    IdentityAngularClosure,
    morel_montry_tau_per_level,
)
from tests.sn._test_helpers import placeholder_materials


# ═══════════════════════════════════════════════════════════════════════
# SPHERE — closure producer == structurally-independent reference (0-ULP)
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.foundation
@pytest.mark.parametrize("N", [8, 16])
def test_sphere_tau_matches_independent_reference(N):
    r"""Closure-produced τ == a HAND-WRITTEN cumulative-weight reference.

    ⛔ **Reference replaced at Q5.6.4 (2026-08-11).** This row compared
    against ``contamination.morel_montry_weights``, described as "a
    DIFFERENT code path to the SAME unclamped BMC weight".  That module's
    successor **delegates to production**, so the comparison would now be
    τ against itself through a wrapper — green forever, and unable to
    detect the drift its name advertises (`coding-standards`: a rewire can
    demote a gate's claim class without touching one line of the body).

    The reference is therefore authored HERE: BMC 2010 Eq. 12 —
    :math:`\mu_{1/2} = -1`, :math:`\mu_{m+1/2} = \mu_{m-1/2} + w_m` — then
    P2 (Eq. 42).

    ⚠ **Not bit-exact, and the reason is the independence itself.**
    ``np.cumsum`` sums pairwise while the producer accumulates
    sequentially, so the two edge ladders differ by FP association — which
    is precisely what makes this a second computation rather than a copy.
    `[M]` the gap GROWS with N, because the edges near :math:`\mu = 0` are
    built by cancellation:

    ===== ===================  ==========
    N     max|ref − prod|      ULP
    ===== ===================  ==========
    4     ``0.000e+00``            0
    8     ``1.776e-15``           16
    16    ``6.550e-15``           59
    32    ``6.106e-15``           55
    64    ``2.247e-13``         2024
    ===== ===================  ==========

    Asserted at ``atol=1e-13``, comfortable for the N ∈ {8, 16} rows here.
    ⛔ **A new row at N ≥ 64 must widen it** — do not read the current
    bound as a claim that holds at every order (`vv-principles` #16: never
    assert tighter than the producer achieves).

    ⚠ It is worth stating what this leg does NOT establish: the sphere's
    cumulative-weight partition is the convention the literature confirms
    verbatim, so this pins the *implementation*, not the *choice*.  The
    choice is argued in
    :func:`~orpheus.sn.sweep.pole_angular_closure.angular_cell_edges_per_level`.
    """
    quad = Quadrature.gauss_legendre(N)
    # Hand-written reference: cumulative-weight edges from −1, then P2.
    w = np.asarray(quad.weights)
    mu_edge = np.concatenate(([-1.0], -1.0 + np.cumsum(w)))
    tau_ref = (np.asarray(quad.mu_x) - mu_edge[:-1]) / np.diff(mu_edge)
    (tau_close,) = morel_montry_tau_per_level(quad, CoordSystem.SPHERICAL)

    np.testing.assert_allclose(
        tau_close, tau_ref, rtol=0.0, atol=1e-13,
        err_msg=(
            f"sphere closure τ != the hand-written cumulative-weight "
            f"reference (N={N}); both are unclamped BMC-2010 Eq.12/42 "
            f"weights and must agree to FP-association noise"
        ),
    )


# ═══════════════════════════════════════════════════════════════════════
# CYLINDER — closure producer == clamp(independent reference) + neg control
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.foundation
@pytest.mark.parametrize("n_phi", [8, 16])
def test_cyl_tau_equals_the_ANALYTIC_closed_form_not_the_chord_convention(
    n_phi,
):
    r"""Cylinder: closure τ == the ANALYTIC P2 closed form on the arc.

    ⛔ **Re-posed at Q5.6.4 (2026-08-11).** This row was
    ``test_cyl_tau_clamp_is_the_only_difference_from_reference``, and its
    whole thesis — *"closure τ == clip(reference τ_raw, ½, 1), and must
    DIFFER from the raw reference where the clamp bites"* — became
    vacuous when the ``[½, 1]`` absorber retired: there is no clamp left
    to be the only difference.  It is **kept, not deleted**, because it
    is the cylinder's only vv-L11 producer gate (two independently
    produced sides).

    ⭐ **What changed about the INDEPENDENCE, and it is an upgrade.** The
    old reference was ``contamination.morel_montry_weights`` — a second
    *procedural* implementation of the same edge recursion, i.e. only
    procedural independence (`vv-principles` L11), and one that silently
    became WRONG when the partition moved.  The reference is now
    the **analytic closed form** obtained by putting the ω-midpoint
    partition through P2 (BMC Eq. 43):

    > :math:`\tau_m = \tfrac12 + \tfrac12\cot\omega_m\,
    > \tan(\Delta\omega/4)`

    derived by hand from the arc geometry, sharing no code path with the
    producer.  That is *structural* independence.

    ⚠ **NEGATIVE control (`vv-principles` #19 — the positive leg alone
    cannot show a gate is loaded on the structure it is credited with).**
    The retired **chord** (η-midpoint) convention is computed inline
    below and the closed form must DIFFER from it — otherwise this row
    would pass just as happily against the partition the carve replaced,
    and would be certifying nothing about the partition choice.

    ⚠ ``angular_differencing.morel_montry_weights`` is NOT used here: it
    delegates to production, so comparing against it would be comparing τ
    with itself through a wrapper (see the module-level note).
    """
    quad = Quadrature.folded_product(n_mu=4, n_phi=n_phi)
    tau_close = morel_montry_tau_per_level(quad, CoordSystem.CYLINDRICAL)

    for p, level_idx in enumerate(quad.level_indices):
        eta = quad.mu_x[level_idx]
        xi = quad.mu_y[level_idx]
        sin_theta = float(np.sqrt(1.0 - quad.mu_z[level_idx[0]] ** 2))
        M = len(level_idx)
        omega = np.arctan2(xi, eta)
        d_omega = np.pi / M

        # (a) the ANALYTIC reference — hand-derived, no shared code path.
        closed_form = 0.5 + 0.5 / np.tan(omega) * np.tan(d_omega / 4.0)
        np.testing.assert_allclose(
            tau_close[p], closed_form, rtol=0.0, atol=1e-13,
            err_msg=(
                f"cylinder level {p} (n_phi={n_phi}): production τ departs "
                f"from the analytic P2 closed form "
                f"½ + ½·cot(ω)·tan(Δω/4)\n"
                f"production ={tau_close[p]}\nclosed form={closed_form}"
            ),
        )

        # (b) NEGATIVE control: the RETIRED chord (η-midpoint) partition.
        chord_edge = np.empty(M + 1)
        chord_edge[0] = -sin_theta
        chord_edge[M] = sin_theta
        chord_edge[1:M] = 0.5 * (eta[:-1] + eta[1:])
        chord_tau = (eta - chord_edge[:-1]) / np.diff(chord_edge)
        gap = float(np.max(np.abs(closed_form - chord_tau)))
        assert gap > 1e-3, (
            f"cylinder level {p} (n_phi={n_phi}): the analytic reference is "
            f"indistinguishable from the RETIRED chord convention "
            f"(max gap {gap:.3e}) — this row cannot then be evidence about "
            f"the partition choice, which is the thing the carve changed"
        )



# ═══════════════════════════════════════════════════════════════════════
# CARTESIAN — IdentityAngularClosure supplies neutral τ = 1.0 (no branch)
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.foundation
def test_identity_closure_tau_is_neutral_one():
    """Cartesian: IdentityAngularClosure supplies τ = 1.0 (no redistribution).

    POSITIVE control: the neutral element is produced WITHOUT a geometry
    branch — the closure TYPE IS the dispatch (Pattern 4).  τ = 1.0 makes
    ``c_out = α_out/1`` and ``c_in = α_in`` with α ≡ 0 (slab), so the M-M
    contribution vanishes identically.
    """
    nx = 8
    edges = np.linspace(0.0, 5.0, nx + 1)
    mesh = Mesh1D(
        edges=edges,
        mat_ids=np.zeros(nx, dtype=int),
        coord=CoordSystem.CARTESIAN,
    )
    quad = Quadrature.gauss_legendre(8)
    sn_mesh = SNMesh(mesh, quad, placeholder_materials())
    closure = IdentityAngularClosure(sn_mesh)

    (tau,) = closure._tau_per_level
    np.testing.assert_array_equal(
        tau,
        np.ones_like(tau),
        err_msg=(
            "IdentityAngularClosure τ is not the neutral 1.0 — Cartesian must "
            "carry no angular redistribution"
        ),
    )
    # No geometry dispatch: morel_montry_tau_per_level raises on Cartesian;
    # Identity supplies the neutral τ by its TYPE, not by a coord branch.
    with pytest.raises(ValueError):
        morel_montry_tau_per_level(quad, CoordSystem.CARTESIAN)
