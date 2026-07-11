r"""L1 contract pins for the R-1 Step 4 G1 carve of
:func:`~orpheus.sn.solver._solve_fixed_source_krylov` onto
:class:`~orpheus.numerics.iteration.KrylovAcceleration` + typed
:class:`~orpheus.sn.angular_flux.AngularFlux`.

Per V&V plan §G1.3 (``.claude/plans/r1_step4_g_verification_plan.md``).

What this file pins
===================

1. *(RETIRED — Phase 1 "gap" un-gate, 2026-06-04.)* The former 2-D
   Cartesian ``NotImplementedError`` pin is INVERTED into a positive
   correctness pin in
   ``tests/sn/solve/test_fixed_source_2d_equivalence.py`` (NEW-1: the
   closed-form Q/Σ_t leg + the SI≡Krylov twin).
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
5. **No legacy eq_map machinery exists** (D-J 2026-05-30: the
   ``build_equation_map_*`` / ``solution_to_angular_flux_*`` /
   ``pack_with_traces`` family was deleted from
   :mod:`orpheus.sn.operators.streaming`; ``_build_rhs_*`` and
   ``_make_sweep_preconditioner`` were already retired pre-D-J).
   The impossibility-of-call contract is enforced by the symbols
   not existing — runtime spy is unnecessary.

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
from orpheus.sn import solve_sn_fixed_source
from orpheus.sn.mesh.augmented_mesh import SNMesh
from orpheus.numerics.quadrature import Quadrature
from orpheus.transport.source_sinks import AngularSourceSink
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


# ── Pin 1 RETIRED (Phase 1 "gap" un-gate, 2026-06-04) ───────────────
#
# The former ``TestTwoDCartesianRaises`` pinned that 2-D Cartesian
# fixed-source Krylov RAISES ``NotImplementedError``.  Phase 1 deleted
# that guard — the path is now the geometry-agnostic structural twin of
# the live 2-D eigenvalue Krylov inner (:meth:`SNSolver._solve_krylov`)
# and the 2-D fixed-source SI path (same ``build_within_group_system`` +
# ``_within_group_krylov``, differing only in ``q_ext``).  The pin is
# INVERTED into a positive correctness pin —
# ``tests/sn/solve/test_fixed_source_2d_equivalence.py`` (NEW-1: the
# closed-form Q/Σ_t leg + the SI≡Krylov twin).


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
            (quad.N, sn_mesh.ng, *sn_mesh.spatial_shape),
        )

        captured: list[np.ndarray] = []
        from orpheus.numerics import iteration as iteration_mod
        original_solve = iteration_mod.KrylovAcceleration.solve

        def spy(self, q_ext, initial_guess=None):
            # B.2d: on a carrying mesh (the sphere param) the driver input is
            # the coupled pair — the per-ordinate source lives on System A's
            # interior member; elsewhere (cylinder — seedless, #229) q_ext is
            # the fused TimedFullField. Value UNMOVED either way (F1: an
            # accessor re-point, the source still folds in bit-equal).
            from orpheus.numerics.coupled_system import CoupledField

            member = (
                q_ext.systems[0] if isinstance(q_ext, CoupledField) else q_ext
            )
            captured.append(np.array(member.interior.values, copy=True))
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
        # Sentinel for the ``solve`` capability node: the curvilinear
        # Q/Σ_t fixed-source equilibrium (Signature 4 — the single most
        # powerful curvilinear diagnostic). cylinder only, to keep the
        # sentinel set minimal. See .claude/plans/sn_sentinel_harness.md.
        pytest.param(
            _cylinder_reflective, "cylinder",
            marks=pytest.mark.sentinel,
        ),
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
        Q_iso = np.full((sn_mesh.ng, *sn_mesh.spatial_shape), q_iso)
        src = AngularSourceSink.from_isotropic(Q_iso, sn_mesh)

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
        # fixed point. D-H.1b: angular_flux is now a TimedFullField;
        # the bulk per-ordinate values live at .interior.values.
        per_ord = result.angular_flux.interior.values
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
        src = AngularSourceSink.from_isotropic(
            np.full((sn_mesh.ng, *sn_mesh.spatial_shape), 1.0), sn_mesh,
        )
        result = solve_sn_fixed_source(
            materials=placeholder_materials(),
            mesh=mesh, quadrature=quad,
            external_source=src.values,
            inner_solver="krylov",
            max_inner=50, inner_tol=1e-10,
        )
        # D-H.1b: angular_flux is now a TimedFullField composite. The
        # boundary is the L2 AngularBoundaryFlux with flat-layout storage; the
        # legacy ``.xmin_face`` / ``.xmax_face`` accessors are replaced
        # by ``.face_view(name)``.
        from orpheus.transport.timed_full_field import TimedFullField
        state = result.angular_flux
        assert isinstance(state, TimedFullField)
        # Sphere has xmax face only (no xmin face; pole at r=0).
        assert "xmax" in state.boundary.layout.faces
        assert state.boundary.face_view("xmax").shape == (quad.N, 1)
        assert "xmin" not in state.boundary.layout.faces  # sphere has no inner face


# ── Pin 5: no legacy eq_map machinery — retired with D-J 2026-05-30
#
# The legacy ``build_equation_map_*`` / ``solution_to_angular_flux_*`` /
# ``pack_with_traces`` codec family was deleted from
# :mod:`orpheus.sn.operators.streaming` in D-J (alongside the bare-ndarray
# packed-vector contract).  ``_build_rhs_{cartesian,spherical,cylindrical}``
# were already deleted in P1.7.  ``_make_sweep_preconditioner`` retired
# pre-D-J.  The runtime-spy ``TestNoLegacyMachineryInCallPath`` test
# (which patched these symbols and asserted zero invocations) retires
# WITH the symbols: the impossibility-of-call contract is now enforced
# structurally (the symbols don't exist), not behaviourally (the
# patched spy counts zero calls).
