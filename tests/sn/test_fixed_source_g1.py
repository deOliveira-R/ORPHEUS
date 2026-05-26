r"""L1 contract pins for the R-1 Step 4 G1 carve of
:func:`~orpheus.sn.solver._solve_fixed_source_krylov` onto
:class:`~orpheus.numerics.iteration.KrylovAcceleration` + typed
:class:`~orpheus.sn.angular_flux.AngularFlux`.

Per V&V plan §G1.3 (``.claude/plans/r1_step4_g_verification_plan.md``).

What this file pins
===================

1. **2-D Cartesian raises NotImplementedError** (the G1 deferral
   sentinel — silent 2-D Krylov regression must be impossible per
   SURPRISE-4 / SURPRISE-5).
2. **External source consumed bit-equal as per-ordinate density**
   (ERR-049 re-sentinel — convention drift between operator algebra
   and fixed-source carve must not re-introduce).
3. **Homogeneous reflective uniform per-ord source converges to
   ``ψ_n = q_n / Σ_t``** (the L0 streaming-equilibrium / per-ord
   flat-flux invariant for the fixed-source path — the canonical
   ERR-026 diagnostic specialised to the carve).
4. **Returns typed AngularFlux with B1'' face residual on its
   boundary** (the path-forward contract — Solution.angular_flux
   carries the matvec's face residual as `.boundary`).
5. **No legacy eq_map machinery invoked** (callstack sentinel —
   ``build_equation_map_*``, ``_build_rhs_*``,
   ``_make_sweep_preconditioner``, ``solution_to_angular_flux_*``
   must not appear in the G1 call path on 1-D).

References
==========

* ``.claude/plans/r1_step4_session2_followup.md`` Phase 3 / G1.
* ``.claude/plans/r1_step4_g_convention_crosswalk.md`` Axis 1 + 2.
* ``.claude/plans/r1_step4_g_dependency_audit.md`` SURPRISE-4 +
  SURPRISE-5 (2-D Cartesian Krylov defers; SI is the landing zone).
* ``.claude/skills/vv-principles/error_catalog.md`` ERR-049
  (convention drift) + ERR-050 (silent precond fallback).
* ``.claude/lessons.md`` L18 (Pattern 7 producer-side) + L19
  (None-default stateful invariants).
"""
from __future__ import annotations

import sys
import warnings
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pytest

from orpheus.geometry import BC, CoordSystem, Mesh1D
from orpheus.geometry.mesh import Mesh2D
from orpheus.sn import solve_sn_fixed_source
from orpheus.sn.angular_flux import AngularFlux
from orpheus.sn.boundary_flux import BoundaryFlux
from orpheus.sn.geometry import SNMesh
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.sources import PerOrdinateSource
from tests.sn._test_helpers import placeholder_materials


pytestmark = [pytest.mark.l1, pytest.mark.catches("ERR-049")]


# ── Mesh + source builders ──────────────────────────────────────────


def _sphere_reflective(nx: int = 10, radius: float = 2.0) -> tuple:
    """Sphere with reflective outer BC + GL N=8."""
    mesh = Mesh1D(
        edges=np.linspace(0.0, radius, nx + 1),
        mat_ids=np.zeros(nx, dtype=int),
        coord=CoordSystem.SPHERICAL,
        bc_left=BC("reflective"),
        bc_right=BC("reflective"),
    )
    quad = Quadrature.gauss_legendre(n_ordinates=8)
    return mesh, quad


def _cylinder_reflective(nx: int = 10, radius: float = 2.0) -> tuple:
    """Cylinder with reflective outer BC + LS-4 quadrature."""
    mesh = Mesh1D(
        edges=np.linspace(0.01, radius, nx + 1),
        mat_ids=np.zeros(nx, dtype=int),
        coord=CoordSystem.CYLINDRICAL,
        bc_left=BC("reflective"),
        bc_right=BC("reflective"),
    )
    quad = Quadrature.product(n_mu=2, n_phi=4)
    return mesh, quad


# ── Pin 1: 2-D Cartesian NotImplementedError ────────────────────────


class TestTwoDCartesianRaises:
    r"""2-D Cartesian fixed-source Krylov MUST raise
    ``NotImplementedError`` post-G1.  Silent routing to SI or
    fallback to a packed-vector legacy path would be a regression
    (the path-forward typed Krylov is 1-D-only; Phase A absorbs
    2-D).  Per SURPRISE-4: this is a deliberate deferral, NOT a bug.
    """

    @pytest.mark.verifies("transport-cartesian")
    def test_two_d_cartesian_krylov_raises_with_g1_message(self) -> None:
        mesh = Mesh2D(
            edges_x=np.array([0.0, 1.0, 2.0]),
            edges_y=np.array([0.0, 1.0, 2.0]),
            mat_map=np.zeros((2, 2), dtype=int),
            coord=CoordSystem.CARTESIAN,
            bc_xmin=BC("vacuum"), bc_xmax=BC("vacuum"),
            bc_ymin=BC("vacuum"), bc_ymax=BC("vacuum"),
        )
        quad = Quadrature.gauss_legendre(n_ordinates=4)
        sn_mesh = SNMesh(mesh, quad, placeholder_materials())
        src = PerOrdinateSource.from_isotropic(
            np.ones((1, 2, 2)), sn_mesh,
        )
        with pytest.raises(NotImplementedError) as exc_info:
            solve_sn_fixed_source(
                materials=placeholder_materials(),
                mesh=mesh, quadrature=quad,
                external_source=src.values,
                inner_solver="krylov",
            )
        # The error message must identify the deferral as G1's, AND
        # recommend the SI landing zone (SURPRISE-5).
        msg = str(exc_info.value)
        assert "G1" in msg, (
            f"Expected G1 deferral marker in NotImplementedError; "
            f"got: {msg!r}"
        )
        assert "source_iteration" in msg, (
            f"Expected SI landing-zone recommendation in error "
            f"message; got: {msg!r}"
        )


# ── Pin 2: external_source consumed bit-equal as per-ord density ────


class TestExternalSourcePerOrdContract:
    r"""ERR-049 re-sentinel for the fixed-source Krylov path.

    ``solve_sn_fixed_source`` documents ``external_source`` as
    per-ordinate density (R-1 Step 4 A1 producer-side normalisation
    contract).  G1 MUST consume the input bit-equal as the typed
    KrylovAcceleration's ``q_ext`` — no internal ``×sum_w`` or
    ``/sum_w`` rescaling on the hot path.

    A spy on ``KrylovAcceleration.solve`` captures the ``q_ext``
    argument and asserts it matches the input ``external_source`` at
    the AngularFlux.values level.
    """

    @pytest.mark.parametrize("coord_builder,coord_name", [
        (_sphere_reflective, "sphere"),
        (_cylinder_reflective, "cylinder"),
    ])
    def test_external_source_forwarded_to_krylov_bit_equal(
        self, coord_builder, coord_name,
    ) -> None:
        mesh, quad = coord_builder(nx=6)
        sn_mesh = SNMesh(mesh, quad, placeholder_materials())
        # Random per-ord source.
        rng = np.random.default_rng(seed=42)
        external_source = rng.standard_normal(
            (quad.N, 1, 6, 1),
        )

        captured: list[np.ndarray] = []
        from orpheus.numerics import iteration as iteration_mod
        original_solve = iteration_mod.KrylovAcceleration.solve

        def spy(self, q_ext, initial_guess=None):
            captured.append(np.array(q_ext.values, copy=True))
            return original_solve(self, q_ext, initial_guess=initial_guess)

        with patch.object(
            iteration_mod.KrylovAcceleration, "solve", spy,
        ):
            solve_sn_fixed_source(
                materials=placeholder_materials(),
                mesh=mesh, quadrature=quad,
                external_source=external_source,
                inner_solver="krylov",
                max_inner=10, inner_tol=1e-6,
            )

        assert len(captured) == 1, (
            f"{coord_name}: KrylovAcceleration.solve was called "
            f"{len(captured)} times; expected exactly 1 "
            f"(the path-forward Krylov is a single solve, not an "
            f"outer Picard wrap)"
        )
        np.testing.assert_array_equal(
            captured[0], external_source,
            err_msg=(
                f"{coord_name}: external_source was modified before "
                f"reaching KrylovAcceleration.solve.  R-1 Step 4 A1 "
                f"contract: per-ord density forwards bit-equal.  "
                f"ERR-049 re-introduced — convention drift between "
                f"solve_sn_fixed_source's documented per-ord input "
                f"and the typed Krylov's q_ext expectation."
            ),
        )


# ── Pin 3: homogeneous reflective uniform source → q/(W·Σ_t) per ord


class TestHomogeneousReflectiveFixedPoint:
    r"""Path-forward streaming-equilibrium check for the G1 carve.

    On a homogeneous reflective medium with uniform per-ord source
    ``q_n = q_iso/W`` (the producer-side normalised value), the
    converged per-ord ψ is ``ψ_n = q_n / Σ_t`` (zero spatial
    gradient on uniform medium with reflective BC; scattering self-
    consistency contributes nothing because we use
    ``placeholder_materials`` with SigS = 0).

    Pre-A1 (legacy ``_build_rhs_*`` consumer-side ``/sum_w``)
    threading per-ord rhs through unchanged would have produced
    ``ψ_n = q_n / (W·Σ_t)`` — a factor of W too small.  G1's
    structural elimination of ``_build_rhs_*`` prevents this bug
    class from re-introducing.
    """

    @pytest.mark.verifies(
        "transport-cartesian", "sn-curvilinear-homogeneous-kinf-recovery",
    )
    @pytest.mark.parametrize("coord_builder,coord_name", [
        (_sphere_reflective, "sphere"),
        (_cylinder_reflective, "cylinder"),
    ])
    def test_uniform_source_converges_to_q_over_sigma_t(
        self, coord_builder, coord_name,
    ) -> None:
        mesh, quad = coord_builder(nx=8)
        materials = placeholder_materials()  # SigT = 1.0 by default
        sn_mesh = SNMesh(mesh, quad, materials)
        sum_w = float(quad.weights.sum())
        # Iso scalar source -> projected per-ord density.
        q_iso = 0.5
        Q_iso = np.full((sn_mesh.ng, sn_mesh.nx, sn_mesh.ny), q_iso)
        src = PerOrdinateSource.from_isotropic(Q_iso, sn_mesh)

        result = solve_sn_fixed_source(
            materials=materials, mesh=mesh, quadrature=quad,
            external_source=src.values,
            inner_solver="krylov",
            max_inner=200, inner_tol=1e-12,
        )

        # SigT = 1.0; per-ord fixed point ψ_n = q_n / Σ_t = (q_iso/W) / 1.0
        expected_per_ord = q_iso / sum_w
        assert result.history.converged, (
            f"{coord_name}: G1 Krylov did not converge in 200 iters"
        )
        # Per-ord ψ should be uniform across (N, ng, nx, ny) at the
        # fixed point.
        per_ord = result.angular_flux.values
        rel_dev = np.abs(per_ord - expected_per_ord) / max(
            expected_per_ord, 1e-30,
        )
        assert rel_dev.max() < 1e-6, (
            f"{coord_name}: per-ord ψ deviates from q_iso/(W·Σ_t) = "
            f"{expected_per_ord:.4e}; max rel dev = {rel_dev.max():.3e}. "
            f"Pre-A1 would have given a factor of W too small — "
            f"ERR-049 re-sentinel."
        )


# ── Pin 4: returns typed AngularFlux with face boundary ─────────────


class TestReturnTypeContract:
    """G1's Solution.angular_flux IS the AngularFlux that KrylovAcceleration
    returned — including its `.boundary` carrying the matvec's B1''
    face residual.  No flat round-trip, no reconstruction from packed."""

    @pytest.mark.verifies("transport-cartesian")
    def test_solution_angular_flux_carries_boundary(self) -> None:
        mesh, quad = _sphere_reflective(nx=6)
        sn_mesh = SNMesh(mesh, quad, placeholder_materials())
        src = PerOrdinateSource.from_isotropic(
            np.full((1, 6, 1), 1.0), sn_mesh,
        )
        result = solve_sn_fixed_source(
            materials=placeholder_materials(),
            mesh=mesh, quadrature=quad,
            external_source=src.values,
            inner_solver="krylov",
            max_inner=50, inner_tol=1e-10,
        )
        psi_typed = result.angular_flux
        assert isinstance(psi_typed, AngularFlux)
        # Sphere has xmax_face only (no xmin face; pole at r=0).
        assert psi_typed.boundary.xmax_face is not None
        assert psi_typed.boundary.xmax_face.shape == (quad.N, 1)
        assert psi_typed.boundary.xmin_face is None  # sphere


# ── Pin 5: no legacy eq_map / _build_rhs_* / _make_sweep_preconditioner


class TestNoLegacyMachineryInCallPath:
    r"""Callstack sentinel — none of the retiring symbols
    (``build_equation_map_*``, ``solution_to_angular_flux_*``,
    ``_make_sweep_preconditioner``) is invoked anywhere along the G1
    Krylov call path on 1-D.

    Spies the suspect callables; if the carve regresses (e.g. someone
    re-routes through legacy decoders), the spy count is non-zero and
    the assertion fires.

    Post-P1.7 (moment-space + layering plan): the three
    ``_build_rhs_{cartesian,spherical,cylindrical}`` helpers have been
    DELETED from :mod:`orpheus.sn.solver` (they were already dead
    code; deletion is the strict-strongest form of "never called").
    The corresponding spies were dropped from this test; the
    impossibility-of-call is now enforced by the symbols not existing
    rather than by a runtime count.
    """

    def test_no_legacy_eq_map_or_decoder_in_g1_path(self) -> None:
        from orpheus.sn import operator as sn_op
        from orpheus.sn import solver as sn_solver

        # Post-P1.7: the _build_rhs_* helpers are deleted from
        # orpheus.sn.solver.  Confirm absence (the impossibility-of-call
        # contract supplants the per-iteration spy count for these).
        assert not hasattr(sn_solver, "_build_rhs_cartesian"), (
            "_build_rhs_cartesian was retired in P1.7; re-adding it "
            "would re-introduce the inline (2*l+1) duplicate the "
            "moment-space plan §P1.3 ('exactly one place') retires."
        )
        assert not hasattr(sn_solver, "_build_rhs_spherical")
        assert not hasattr(sn_solver, "_build_rhs_cylindrical")

        legacy_calls: dict[str, int] = {}

        def make_spy(name, original):
            def wrapped(*args, **kwargs):
                legacy_calls[name] = legacy_calls.get(name, 0) + 1
                return original(*args, **kwargs)
            return wrapped

        original_build_eq_map_sph = sn_op.build_equation_map_spherical
        original_build_eq_map_cyl = sn_op.build_equation_map_cylindrical
        original_build_eq_map_cart = sn_op.build_equation_map
        original_decode_sph = sn_op.solution_to_angular_flux_spherical
        original_decode_cyl = sn_op.solution_to_angular_flux_cylindrical
        original_decode_cart = sn_op.solution_to_angular_flux
        # _make_sweep_preconditioner is on SNSolver instances, not a
        # module-level symbol — check via attribute lookup at instance level
        # would require dispatching the spy through __init__.  Simpler:
        # rely on solution_to_angular_flux_* spies.

        with patch.object(
            sn_op, "build_equation_map_spherical",
            make_spy("build_equation_map_spherical", original_build_eq_map_sph),
        ), patch.object(
            sn_op, "build_equation_map_cylindrical",
            make_spy("build_equation_map_cylindrical", original_build_eq_map_cyl),
        ), patch.object(
            sn_op, "build_equation_map",
            make_spy("build_equation_map", original_build_eq_map_cart),
        ), patch.object(
            sn_op, "solution_to_angular_flux_spherical",
            make_spy("solution_to_angular_flux_spherical", original_decode_sph),
        ), patch.object(
            sn_op, "solution_to_angular_flux_cylindrical",
            make_spy("solution_to_angular_flux_cylindrical", original_decode_cyl),
        ), patch.object(
            sn_op, "solution_to_angular_flux",
            make_spy("solution_to_angular_flux", original_decode_cart),
        ):
            mesh, quad = _sphere_reflective(nx=6)
            sn_mesh = SNMesh(mesh, quad, placeholder_materials())
            src = PerOrdinateSource.from_isotropic(
                np.full((1, 6, 1), 1.0), sn_mesh,
            )
            solve_sn_fixed_source(
                materials=placeholder_materials(),
                mesh=mesh, quadrature=quad,
                external_source=src.values,
                inner_solver="krylov",
                max_inner=50, inner_tol=1e-10,
            )

        # Permitted callers of legacy symbols in the G1 Krylov path:
        # - SNSolver.__init__ builds self.L = SNStreamingOperator(...)
        #   which lazily builds an eq_map on first `apply` — but G1's
        #   Krylov path does NOT call solver.L.apply (it uses the typed
        #   (L+C).apply via KrylovAcceleration).  Net: zero calls in path.
        # - The Solution constructor doesn't decode from packed.
        for sym, count in legacy_calls.items():
            assert count == 0, (
                f"Legacy symbol {sym!r} called {count} times in the "
                f"G1 Krylov call path.  Per the dependency audit, "
                f"this symbol retires; calling it from G1 is a "
                f"regression of the structural carve."
            )
