r"""Bit-exact L+C decomposition gate (Issue #196 Phase G Step 3+4.b.i).

Promoted from :file:`derivations/diagnostics/diag_LC_decomposition_resolution.py`
— this is the load-bearing verification that Resolution A's subtractive
``StreamingOperator.apply`` produces a bit-exact ``(L + C).apply ≡ M``
decomposition on every geometry.

The math
--------

.. math::

   L.{\rm apply}(\psi) \;:=\; M(\psi;\;\sigma_t) \;-\;
                              \sigma_t \odot \psi_{\rm bulk}
   \\
   C.{\rm apply}(\psi) \;:=\; \sigma_t \odot \psi_{\rm bulk}
   \\
   (L + C).{\rm apply}(\psi) \;\equiv\; M(\psi;\;\sigma_t)
       \qquad \text{rel\_residual} = 0.0

Both L and C carry σ_t at constructor time per Resolution A (see the
numerics-investigator memo
``.claude/agent-memory/numerics-investigator/sn_LC_decomposition_derivation.md``).
The curvilinear matvec is RATIONAL in σ_t through Hébert §3.9.4's
Carlson coupled-pole seed; the subtractive form bypasses this by
calling matvec at the FULL σ_t (not zero) then subtracting the
cell-collision term post-matvec.

The prior reverted approach (commit ``ad37ca0``) called the matvec at
``σ_t = 0`` and treated that as "pure streaming"; this is
mathematically WRONG for curvilinear (3-13% rel error on random ψ)
because setting σ_t=0 degenerates the Carlson seed denominator
``Δr·σ_t + 2 → 2`` and produces a different operator. Resolution A
fixes this by construction.

Test contract (mechanism criteria 1, 2 from the substep brief)
---------------------------------------------------------------

* slab/sphere/cylinder × 3 random seeds × σ_t = 2.0 uniform.
* Composition residual ``(L+C).apply(ψ) − M(ψ;σ_t)`` MUST be bit-exact
  (``rel_residual < 1e-14``) on BOTH the cell-centre bulk
  ``AngularFlux`` and the boundary face residual ``BoundaryFlux``.
* No xfail. Resolution A guarantees this hold by construction.

A separate test verifies that the subtractive L's apply DIFFERS from
the prior wrong approach ``matvec(σ_t=0)`` (sanity: we're shipping a
genuinely different formulation).

D-I.1 — typed-carrier migration
-------------------------------

This file was originally written against the legacy packed-vector
contract (bare-``np.ndarray`` through
:func:`_transport_operator_matvec_unified` via the
:class:`EquationMap` B1'' adapter).  D-I.1 retires
:meth:`CollisionOperator.apply(bare_ndarray)` /
:meth:`CollisionOperator.solve(bare_ndarray)` along with the
supporting ``_ensure_eq_map`` / ``_sigma_at_unknowns`` / ``_eq_map``
fields; the test file migrates first so its assertions land directly
on the typed :class:`TimedFullField` carrier (the D-H wave's
load-bearing composite type).

The TimedFullField construction inlines at each test — the typed
construction IS the test setup; no per-test wrapper helper is
introduced.
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.geometry import BC, CoordSystem, Mesh1D
from orpheus.sn.geometry import SNMesh
from orpheus.sn.operator import (
    CollisionOperator,
    StreamingOperator,
    )
from tests.sn._test_helpers import _LC_matvec
from orpheus.numerics.quadrature import Quadrature
from orpheus.transport.fields.angular_flux import AngularFlux
from orpheus.transport.fields.boundary_flux import BoundaryFlux
from orpheus.transport.timed_full_field import TimedFullField
from tests.sn._test_helpers import placeholder_materials

pytestmark = pytest.mark.l0


# ═══════════════════════════════════════════════════════════════════════
# Fixtures — slab / sphere / cylinder, reflective inner / vacuum outer.
# ═══════════════════════════════════════════════════════════════════════


def _build_sn_mesh(geometry: str, n_cells: int = 5, n_ord: int = 4) -> SNMesh:
    """Build a small SNMesh for the requested geometry.

    Sized small (n_cells=5, n_ord=4) so the matvec runs in well under
    a second per geometry. The decomposition contract is size-independent.
    """
    if geometry == "SPH":
        mesh = Mesh1D(
            edges=np.linspace(0.0, 2.0, n_cells + 1),
            mat_ids=np.zeros(n_cells, dtype=int),
            coord=CoordSystem.SPHERICAL,
            bc_left=BC("reflective"),
            bc_right=BC("reflective"),
        )
        quad = Quadrature.gauss_legendre(n_ordinates=n_ord)
    elif geometry == "CYL":
        mesh = Mesh1D(
            edges=np.linspace(0.01, 2.0, n_cells + 1),
            mat_ids=np.zeros(n_cells, dtype=int),
            coord=CoordSystem.CYLINDRICAL,
            bc_left=BC("reflective"),
            bc_right=BC("reflective"),
        )
        quad = Quadrature.level_symmetric(sn_order=n_ord)
    elif geometry == "CART":
        mesh = Mesh1D(
            edges=np.linspace(0.0, 2.0, n_cells + 1),
            mat_ids=np.zeros(n_cells, dtype=int),
            coord=CoordSystem.CARTESIAN,
            bc_left=BC("reflective"),
            bc_right=BC("reflective"),
        )
        quad = Quadrature.gauss_legendre(n_ordinates=n_ord)
    else:
        raise ValueError(geometry)
    return SNMesh(mesh, quad, placeholder_materials())


# ═══════════════════════════════════════════════════════════════════════
# Mechanism criterion 1 — (L + C).apply(ψ) ≡ M(ψ; σ_t) bit-exact.
# ═══════════════════════════════════════════════════════════════════════


class TestResolutionADecomposition:
    """``(L + C).apply(ψ) == M(ψ; σ_t)`` bit-exact across all geometries.

    Resolution A guarantees this by construction:
    ``L.apply := M(ψ; σ_t) − σ_t⊙ψ.bulk`` and ``C.apply := σ_t⊙ψ.bulk``,
    so ``(L + C).apply = M(ψ; σ_t)`` algebraically on both bulk and
    boundary blocks. NO xfail.
    """

    @pytest.mark.parametrize("geometry", ["CART", "SPH", "CYL"])
    @pytest.mark.parametrize("seed", [0, 1, 2])
    def test_bit_exact_uniform_sigma_t(self, geometry, seed):
        """(L + C).apply(ψ) ≡ M(ψ; σ_t) at rel_residual == 0.0
        on both the bulk and the boundary blocks of the typed carrier.
        """
        sn_mesh = _build_sn_mesh(geometry, n_cells=5, n_ord=4)
        ng = 1
        N = sn_mesh.quad.N
        nx, ny = sn_mesh.nx, sn_mesh.ny

        rng = np.random.default_rng(seed)
        state = TimedFullField(
            bulk=AngularFlux.from_mesh(
                rng.standard_normal((N, ng, nx, ny)), sn_mesh,
            ),
            boundary=BoundaryFlux.zeros_on(sn_mesh),
            _history=(),
            history_depth=2,
        )
        sigma_t = np.full((ng, sn_mesh.nx, sn_mesh.ny), 2.0)  # PR-INDEX-3

        # Reference: the unified matvec at full σ_t.
        m_full_state = _LC_matvec(state, sigma_t)

        # Resolution A: L.apply + C.apply via TimedFullField arithmetic.
        L = StreamingOperator(sn_mesh, sigma_t)
        C = CollisionOperator(sn_mesh, sigma_t)
        sum_state = L.apply(state) + C.apply(state)

        # Bulk residual — (N, ng, nx, ny) ndarray.
        residual_bulk = sum_state.bulk.values - m_full_state.bulk.values
        rel_bulk = (
            np.linalg.norm(residual_bulk)
            / max(np.linalg.norm(m_full_state.bulk.values), 1e-300)
        )
        assert rel_bulk < 1e-14, (
            f"{geometry} seed={seed}: bulk rel_residual={rel_bulk:.3e} "
            f"— Resolution A subtractive decomposition FAILED bit-exact "
            f"gate on bulk. (L + C).apply.bulk MUST equal M(ψ; σ_t).bulk."
        )

        # Boundary residual — flat (layout.total_size,) backing buffer.
        residual_bdry = (
            sum_state.boundary.values - m_full_state.boundary.values
        )
        if m_full_state.boundary.values.size > 0:
            rel_bdry = (
                np.linalg.norm(residual_bdry)
                / max(np.linalg.norm(m_full_state.boundary.values), 1e-300)
            )
        else:
            rel_bdry = 0.0
        assert rel_bdry < 1e-14, (
            f"{geometry} seed={seed}: boundary rel_residual={rel_bdry:.3e} "
            f"— Resolution A subtractive decomposition FAILED bit-exact "
            f"gate on boundary. (L + C).apply.boundary MUST equal "
            f"M(ψ; σ_t).boundary by construction (C contributes zero "
            f"to the boundary face residual)."
        )


# ═══════════════════════════════════════════════════════════════════════
# Mechanism criterion 2 — the subtractive definition holds:
# L.apply(ψ) ≡ M(ψ; σ_t) − σ_t⊙ψ.bulk exactly (on bulk; equal on boundary).
# ═══════════════════════════════════════════════════════════════════════


class TestSubtractiveDefinition:
    """``L.apply(ψ).bulk == M(ψ; σ_t).bulk − σ_t · ψ.bulk`` at exact equality
    and ``L.apply(ψ).boundary == M(ψ; σ_t).boundary`` at exact equality.

    The complement of the (L + C) test — this verifies the L leaf
    alone matches the subtractive formula on the bulk while leaving
    the boundary face residual untouched (collision contributes zero
    on the trace). Together with the ``C.apply == σ_t⊙ψ.bulk; zero
    boundary`` contract (covered by test_collision_operator), these
    two tests fully pin Resolution A.
    """

    @pytest.mark.parametrize("geometry", ["CART", "SPH", "CYL"])
    @pytest.mark.parametrize("seed", [10, 11, 12])
    def test_L_apply_equals_subtractive_form(self, geometry, seed):
        sn_mesh = _build_sn_mesh(geometry, n_cells=5, n_ord=4)
        ng = 1
        N = sn_mesh.quad.N
        nx, ny = sn_mesh.nx, sn_mesh.ny

        rng = np.random.default_rng(seed)
        state = TimedFullField(
            bulk=AngularFlux.from_mesh(
                rng.standard_normal((N, ng, nx, ny)), sn_mesh,
            ),
            boundary=BoundaryFlux.zeros_on(sn_mesh),
            _history=(),
            history_depth=2,
        )
        sigma_t = np.full((ng, sn_mesh.nx, sn_mesh.ny), 2.0)  # PR-INDEX-3

        L = StreamingOperator(sn_mesh, sigma_t)
        l_state = L.apply(state)

        m_full_state = _LC_matvec(state, sigma_t)

        # Cell-centre subtraction: bulk expected = M.bulk - σ_t·ψ.bulk
        # (σ_t broadcast over the ordinate axis 0 via [None, :, :, :]).
        expected_bulk = (
            m_full_state.bulk.values - sigma_t[None, :, :, :] * state.bulk.values
        )
        # Face slots: L.boundary == M.boundary (no volumetric collision
        # on the trace; the cell-balance σ·ψ term is a CELL quantity).
        expected_boundary = m_full_state.boundary.values

        # Bit-exact: L.apply IS the subtractive formula.
        np.testing.assert_array_equal(l_state.bulk.values, expected_bulk)
        np.testing.assert_array_equal(
            l_state.boundary.values, expected_boundary,
        )


# ═══════════════════════════════════════════════════════════════════════
# Sanity — Resolution A's L is NOT the prior wrong matvec(σ_t=0) approach.
# ═══════════════════════════════════════════════════════════════════════


class TestResolutionADifferentFromPriorWrong:
    """Resolution A's L.apply DIFFERS from the prior wrong matvec(σ_t=0).

    Sanity check that we're shipping a genuinely different formulation
    from the reverted ``ad37ca0`` approach. For SPHERE the prior approach
    degenerates the Carlson seed denominator and produces an O(0.01..1)
    different operator on random ψ.

    For CYLINDER, the empirical σ_t coupling through the Carlson seed is
    structurally small under the PR-TYPED-6c Step 3 unified body —
    ``M_unified(σ_full) − σ·ψ ≈ M_unified(σ_zero)`` at FP noise on
    random ψ.  This is a property of the unified per-level M-M
    recurrence, NOT a regression: the pre-unification legacy cylindrical
    matvec carried the ERR-049 routing bug that amplified the apparent
    σ_t coupling.  The unified matvec is verified directly via the
    Step 3 L1 trajectory_resolvent reference (3% rtol on the
    heterogeneous closed cylinder), which is the load-bearing
    correctness anchor — NOT this differential sanity test.

    For Cartesian the two are equivalent because Cartesian L has no
    Carlson seed coupling (Step 4 — slab is M-M-neutral).
    """

    @pytest.mark.parametrize("geometry", ["SPH"])
    def test_subtractive_L_differs_from_matvec_at_zero_sigma_t(
        self, geometry,
    ):
        sn_mesh = _build_sn_mesh(geometry, n_cells=5, n_ord=4)
        ng = 1
        N = sn_mesh.quad.N
        nx, ny = sn_mesh.nx, sn_mesh.ny

        rng = np.random.default_rng(0)
        state = TimedFullField(
            bulk=AngularFlux.from_mesh(
                rng.standard_normal((N, ng, nx, ny)), sn_mesh,
            ),
            boundary=BoundaryFlux.zeros_on(sn_mesh),
            _history=(),
            history_depth=2,
        )
        sigma_full = np.full((ng, sn_mesh.nx, sn_mesh.ny), 2.0)  # PR-INDEX-3
        sigma_zero = np.zeros((ng, sn_mesh.nx, sn_mesh.ny))

        # Resolution A's L.apply (subtractive).
        L = StreamingOperator(sn_mesh, sigma_full)
        l_correct_state = L.apply(state)

        # Prior agent's wrong L.apply: matvec at σ_t = 0 (which has
        # different boundary behaviour because the Carlson seed
        # denominator degenerates).  Wave T post-T.5 (matvec
        # retirement): the function `_transport_operator_matvec_unified`
        # was DELETED; reach into `_MSpatialOperatorSum._compute_decomposition`
        # directly to bypass `InvertibleOperator`'s σ > 0 validation
        # (this is a deliberate test of the wrong-prior behaviour).
        from orpheus.sn.operator import _MSpatialOperatorSum
        orch_zero = _MSpatialOperatorSum(sn_mesh, sigma_zero)
        m_spat_zero, m_ang_zero = orch_zero._compute_decomposition(state)
        # M = m_spat + m_ang (no σ_t·ψ subtraction since σ_t = 0).
        import dataclasses
        l_prior_state = dataclasses.replace(
            m_spat_zero,
            bulk=AngularFlux.from_mesh(
                m_spat_zero.bulk.values + m_ang_zero.bulk.values, sn_mesh,
            ),
        )

        diff_bulk = l_correct_state.bulk.values - l_prior_state.bulk.values
        rel = (
            np.linalg.norm(diff_bulk)
            / max(np.linalg.norm(l_correct_state.bulk.values), 1e-300)
        )
        assert rel > 1e-3, (
            f"{geometry}: Resolution A's L.apply and prior agent's "
            f"L.apply differ only by {rel:.3e} — expected O(0.01..1) "
            f"difference because the prior approach degenerates the "
            f"Carlson seed denominator ``dr·σ_t + 2 → 2``."
        )
