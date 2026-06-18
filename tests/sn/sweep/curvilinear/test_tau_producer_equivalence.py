r"""Issue #236 Phase 2 — τ producer-equivalence gate (Leg 1).

The :class:`~orpheus.sn.spatial.pole_angular_closure.MorelMontryAngularSweep`
closure now PRODUCES the Morel–Montry angular weight τ from the quadrature
``(μ, w, levels)`` it already binds (an angular-scheme property), instead of
reading it back from the streaming-GEOMETRY factory (``reduced_operator.py``).

Step C (the geometry-τ retirement) deleted the geometry-side τ producer, so
the former ``*_equals_geometry_factory_0ulp`` legs (closure τ vs
``spherical_streaming(...).tau_mm`` / ``cylindrical_streaming(...).tau_mm_per_level``)
are now vacuous and have been removed.  The structural-independence floor is
carried entirely by the surviving legs below — the closure τ producer is
pinned against the STRUCTURALLY-INDEPENDENT ``contamination.morel_montry_weights``
(a DIFFERENT code path to the same BMC-2010-Eq.43 weight; vv-principles L11),
which never depended on the geometry producer.

This gate proves the closure-produced τ:

* equals a STRUCTURALLY-INDEPENDENT reference
  (``contamination.morel_montry_weights`` — a different code path to the SAME
  BMC-2010-Eq.43 weight; vv-principles L11);
* on the cylinder, equals ``clamp(reference τ_raw)`` and NOT the raw reference
  — the clamp is a real producer-side transform, pinned by a NEGATIVE control
  (vv anti-pattern #11);
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
from orpheus.sn.geometry import SNMesh
from orpheus.sn.spatial.pole_angular_closure import (
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
def test_cyl_tau_clamp_is_the_only_difference_from_reference(n_phi):
    """Cylinder: closure τ == ``clamp(contamination τ_raw)``, NOT the raw τ.

    NEGATIVE-control companion (vv anti-pattern #11): the independent
    reference is the UNCLAMPED τ_raw; the closure τ must equal
    ``clip(τ_raw, ½, 1)`` everywhere, and must DIFFER from the raw reference
    exactly where the clamp bites (τ_raw < ½, at the most-inward ordinate
    where ``eta == eta_edge == -sinθ`` so τ_raw = 0).  This pins the clamp
    as a real transform, not an accident, AND proves the producer is not the
    naked unclamped reference.
    """
    quad = Quadrature.product(n_mu=4, n_phi=n_phi)
    tau_raw_ref = morel_montry_weights(quad, "cylindrical")  # list[(M,)], raw
    tau_clamped_ref = [np.clip(t, 0.5, 1.0) for t in tau_raw_ref]
    tau_close = morel_montry_tau_per_level(quad, CoordSystem.CYLINDRICAL)

    # Guard-the-guard: the clamp MUST bite on at least one level, else the
    # negative control is vacuous (it would compare τ against itself).
    assert any(t.min() < 0.5 for t in tau_raw_ref), (
        "fixture too weak — clamp never bites; pick a quadrature where the "
        f"most-inward ordinate has τ_raw < ½ (n_phi={n_phi}); "
        f"mins={[float(t.min()) for t in tau_raw_ref]}"
    )

    clamp_bites_somewhere = False
    for p, (tc, t_raw, t_clamped) in enumerate(
        zip(tau_close, tau_raw_ref, tau_clamped_ref)
    ):
        # The closure τ equals the CLAMPED reference, level by level.
        assert np.array_equal(tc, t_clamped), (
            f"cylinder level {p}: closure τ != clamp(reference τ_raw) "
            f"(n_phi={n_phi})\nclosure={tc}\nclamp(ref)={t_clamped}"
        )
        # Where the clamp bit, the closure τ must DIFFER from the raw
        # reference — proving the closure is NOT the naked unclamped τ.
        bit = t_raw < 0.5
        if np.any(bit):
            clamp_bites_somewhere = True
            assert not np.array_equal(tc[bit], t_raw[bit]), (
                f"cylinder level {p}: closure τ == raw reference where the "
                f"clamp should bite (n_phi={n_phi}) — the clamp is missing "
                f"from the producer\nclosure={tc[bit]}\nraw={t_raw[bit]}"
            )

    assert clamp_bites_somewhere, (
        "no level exercised the clamp-difference branch; negative control "
        f"vacuous (n_phi={n_phi})"
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
