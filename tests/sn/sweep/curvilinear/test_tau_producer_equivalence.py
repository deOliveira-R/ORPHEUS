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

* **SPHERE** → ``contamination.morel_montry_weights``.  Still valid: that
  module's spherical arm is the cumulative-weight recursion (BMC 2010
  Eq. 12/42), the same convention production uses, reached by a different
  code path.  This is PROCEDURAL independence (`vv-principles` L11).
* **CYLINDER** → the **analytic closed form**
  :math:`\tau_m = \tfrac12 + \tfrac12\cot\omega_m\tan(\Delta\omega/4)`.
  ⛔ ``contamination.morel_montry_weights`` is NOT usable here: `[M]` its
  cylindrical arm still builds the **RETIRED η-midpoint (chord)** edges
  (``contamination.py:64-66``), so as of Q5.6.4 it disagrees with
  production by construction.  Migrating it is tracked with the analysis
  module; until then it serves as this gate's NEGATIVE control rather
  than its reference.  The closed form is STRUCTURAL independence — it
  shares no code path with the producer at all, which is a strictly
  stronger footing than the procedural twin it replaced.

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

from orpheus.derivations.discrete.sn.contamination import morel_montry_weights
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
    """Closure-produced τ == ``contamination.morel_montry_weights`` (sphere).

    Structurally-independent leg (vv-principles L11): contamination.py is a
    DIFFERENT code path to the SAME unclamped BMC weight.  Sphere is
    UNCLAMPED on both sides ⇒ exact (0-ULP) equality.
    """
    quad = Quadrature.gauss_legendre(N)
    tau_ref = morel_montry_weights(quad, "spherical")
    (tau_close,) = morel_montry_tau_per_level(quad, CoordSystem.SPHERICAL)

    np.testing.assert_array_equal(
        tau_close,
        tau_ref,
        err_msg=(
            f"sphere closure τ != independent reference τ (N={N}); both are "
            f"unclamped BMC-2010-Eq.43 weights and must agree exactly"
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
    procedural independence (`vv-principles` L11).  The reference is now
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

    ⚠ ``contamination.morel_montry_weights`` still carries the RETIRED
    η-midpoint convention, so it is deliberately NOT used as the
    reference here; it appears only as the negative control's sibling
    (see the module-level note).
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
