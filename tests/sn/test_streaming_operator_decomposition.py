r"""Bit-exact L+C decomposition gate (Issue #196 Phase G Step 3+4.b.i).

Promoted from :file:`derivations/diagnostics/diag_LC_decomposition_resolution.py`
— this is the load-bearing verification that Resolution A's subtractive
``StreamingOperator.apply`` produces a bit-exact ``(L + C).apply ≡ M``
decomposition on every geometry.

The math
--------

.. math::

   L.{\rm apply}(\psi) \;:=\; M(\psi;\;\sigma_t) \;-\;
                              \sigma_t \odot \psi_{\rm packed}
   \\
   C.{\rm apply}(\psi) \;:=\; \sigma_t \odot \psi_{\rm packed}
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
  (``rel_residual < 1e-14``).
* No xfail. Resolution A guarantees this hold by construction.

A separate test verifies that the subtractive L's apply DIFFERS from
the prior wrong approach ``matvec(σ_t=0)`` (sanity: we're shipping a
genuinely different formulation).

D-H.2-C1 note
-------------
All tests in this file exercise the bare-``np.ndarray`` flat-vector
path through :func:`transport_operator_matvec_unified` (the legacy
matvec kernel).  ``_call_matvec`` internally adapts via the legacy
:class:`orpheus.sn.angular_flux.AngularFlux` +
:class:`orpheus.sn.boundary_flux.BoundaryFlux` pair as a transient
bridge to the unified matvec.  The composite migration deferred
this file to C4 because the matvec kernel itself needs the L2-native
rewrite before the test fixtures can express the input as a
:class:`TimedFullField`.  Tests stay legacy until C4.
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.geometry import BC, CoordSystem, Mesh1D
from orpheus.sn.geometry import SNMesh
from orpheus.sn.operator import (
    CollisionOperator,
    StreamingOperator,
    build_equation_map,
    build_equation_map_with_traces,
    pack_with_traces,
    solution_to_angular_flux_with_traces,
    transport_operator_matvec_unified,
)
from orpheus.numerics.quadrature import Quadrature
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


def _eq_map_for(sn_mesh: SNMesh, ng: int):
    """Geometry-dispatched EquationMap factory matching
    :meth:`StreamingOperator._ensure_eq_map` after PR-TYPED-6.5 Phase 3b.

    1-D paths use the B1'' face-aware layout (cell-centres + face
    slots); 2-D Cartesian retains the legacy FD layout.
    """
    nx, ny = sn_mesh.nx, sn_mesh.ny
    quad = sn_mesh.quad
    curv = getattr(sn_mesh, "curvature", None)
    if curv is None and ny > 1:
        return build_equation_map(nx, ny, quad, ng)
    has_inner_bc = (curv is None)
    return build_equation_map_with_traces(
        nx, quad, ng, has_inner_bc=has_inner_bc,
    )


def _call_matvec(sn_mesh: SNMesh, psi_vec: np.ndarray,
                 sigma_t_arr: np.ndarray, eq_map, ng: int) -> np.ndarray:
    """Geometry-dispatched matvec at the supplied σ_t — the reference value.

    Mirrors :meth:`StreamingOperator.apply`'s internal dispatch.  Issue #197
    PR-TYPED-6c Step 5 — 1-D slab / sphere / cylinder all route through
    :func:`transport_operator_matvec_unified`.

    D-H.2-C4e.6 (2026-05-29) — the 2-D Cartesian dead branch (predicate
    ``curv is None and ny > 1``) retired with the legacy
    ``transport_operator_matvec``: this helper is only called from 1-D
    Cartesian (``CART``) tests where ``ny == 1``, so the 2-D fallthrough
    never fired in practice.  The 1-D B1'' face-aware L2-native path
    below is the only path now.

    PR-TYPED-6.5 Phase 3b — 1-D paths consume the B1'' face-aware
    packed layout via :func:`solution_to_angular_flux_with_traces` and
    :func:`pack_with_traces`.
    """
    nx, ny = sn_mesh.nx, sn_mesh.ny
    quad = sn_mesh.quad

    # 1-D B1'' face-aware path.  D-H.2-C4c — L2-native
    # ``transport_operator_matvec_unified`` consumes ``TimedFullField``;
    # decode packed face slots into the L2 ``face_view`` writable arrays.
    from orpheus.transport.fields.angular_flux import (
        AngularFlux,
    )
    from orpheus.transport.fields.boundary_flux import (
        BoundaryFlux,
    )
    from orpheus.transport.timed_full_field import TimedFullField
    psi_cell, psi_face_outer, psi_face_inner = (
        solution_to_angular_flux_with_traces(
            psi_vec, eq_map, nx, ng, N=quad.N,
        )
    )
    boundary_in = BoundaryFlux.zeros_for_sn_mesh(sn_mesh)
    if eq_map.n_face_outer > 0:
        boundary_in.face_view("xmax")[eq_map.face_outer_ordinate, :] = (
            psi_face_outer
        )
    if eq_map.n_face_inner > 0 and "xmin" in boundary_in.layout.faces:
        boundary_in.face_view("xmin")[eq_map.face_inner_ordinate, :] = (
            psi_face_inner
        )
    composite_in = TimedFullField(
        bulk=AngularFlux.from_mesh(psi_cell, sn_mesh),
        boundary=boundary_in,
        _history=(),
        history_depth=2,
    )
    result = transport_operator_matvec_unified(composite_in, sigma_t_arr)
    m_cell = result.bulk.values
    m_face_outer = (
        result.boundary.face_view("xmax")[eq_map.face_outer_ordinate, :]
        if eq_map.n_face_outer > 0 else None
    )
    m_face_inner = (
        result.boundary.face_view("xmin")[eq_map.face_inner_ordinate, :]
        if eq_map.n_face_inner > 0
        and "xmin" in result.boundary.layout.faces
        else None
    )
    return pack_with_traces(m_cell, m_face_outer, m_face_inner, eq_map)


# ═══════════════════════════════════════════════════════════════════════
# Mechanism criterion 1 — (L + C).apply(ψ) ≡ M(ψ; σ_t) bit-exact.
# ═══════════════════════════════════════════════════════════════════════


class TestResolutionADecomposition:
    """``(L + C).apply(ψ) == M(ψ; σ_t)`` bit-exact across all geometries.

    Resolution A guarantees this by construction:
    ``L.apply := M(ψ; σ_t) − σ_t⊙ψ`` and ``C.apply := σ_t⊙ψ``, so
    ``(L + C).apply = M(ψ; σ_t)`` algebraically. NO xfail.
    """

    @pytest.mark.parametrize("geometry", ["CART", "SPH", "CYL"])
    @pytest.mark.parametrize("seed", [0, 1, 2])
    def test_bit_exact_uniform_sigma_t(self, geometry, seed):
        """(L + C).apply(ψ) ≡ M(ψ; σ_t) at rel_residual == 0.0."""
        sn_mesh = _build_sn_mesh(geometry, n_cells=5, n_ord=4)
        ng = 1
        eq_map = _eq_map_for(sn_mesh, ng=ng)

        rng = np.random.default_rng(seed)
        psi_vec = rng.standard_normal(eq_map.n_unknowns).astype(np.float64)
        sigma_t = np.full((ng, sn_mesh.nx, sn_mesh.ny), 2.0)  # PR-INDEX-3

        # Reference: the matvec at full σ_t.
        m_full = _call_matvec(sn_mesh, psi_vec, sigma_t, eq_map, ng)

        # Resolution A: L.apply + C.apply.
        L = StreamingOperator(sn_mesh, sigma_t)
        C = CollisionOperator(sn_mesh, sigma_t)
        sum_apply = L.apply(psi_vec) + C.apply(psi_vec)

        residual = sum_apply - m_full
        rel_residual = (
            np.linalg.norm(residual)
            / max(np.linalg.norm(m_full), 1e-300)
        )
        assert rel_residual < 1e-14, (
            f"{geometry} seed={seed}: rel_residual={rel_residual:.3e} "
            f"— Resolution A subtractive decomposition FAILED bit-exact "
            f"gate. (L + C).apply MUST equal M(ψ; σ_t) by construction."
        )


# ═══════════════════════════════════════════════════════════════════════
# Mechanism criterion 2 — the subtractive definition holds:
# L.apply(ψ) ≡ M(ψ; σ_t) − σ_t⊙ψ exactly.
# ═══════════════════════════════════════════════════════════════════════


class TestSubtractiveDefinition:
    """``L.apply(ψ) == M(ψ; σ_t) − σ_t_packed * ψ`` at exact equality.

    The complement of the (L + C) test — this verifies the L leaf
    alone matches the subtractive formula. Together with the
    ``C.apply == σ_t⊙ψ`` contract (covered by test_collision_operator),
    these two tests fully pin Resolution A.
    """

    @pytest.mark.parametrize("geometry", ["CART", "SPH", "CYL"])
    @pytest.mark.parametrize("seed", [10, 11, 12])
    def test_L_apply_equals_subtractive_form(self, geometry, seed):
        sn_mesh = _build_sn_mesh(geometry, n_cells=5, n_ord=4)
        ng = 1
        eq_map = _eq_map_for(sn_mesh, ng=ng)

        rng = np.random.default_rng(seed)
        psi_vec = rng.standard_normal(eq_map.n_unknowns).astype(np.float64)
        sigma_t = np.full((ng, sn_mesh.nx, sn_mesh.ny), 2.0)  # PR-INDEX-3

        L = StreamingOperator(sn_mesh, sigma_t)
        l_apply = L.apply(psi_vec)

        m_full = _call_matvec(sn_mesh, psi_vec, sigma_t, eq_map, ng)
        # PR-INDEX-3: σ_t shape (ng, nx, ny); advanced index gives (ng, n_eq).
        # PR-TYPED-6.5 Phase 3b: B1'' face slots carry no volumetric
        # collision; build the σ_t subtraction vector for the cell-
        # centre block only and zero-pad the face block.
        sigma_packed_cell = sigma_t[
            :, eq_map.ix, eq_map.iy
        ].ravel(order='F')
        n_cell_scalars = eq_map.n_eq * ng
        expected = m_full.copy()
        expected[:n_cell_scalars] -= sigma_packed_cell * psi_vec[:n_cell_scalars]
        # face slots: expected[n_cell_scalars:] = m_full[face] − 0

        # Bit-exact: L.apply IS the subtractive formula.
        np.testing.assert_array_equal(l_apply, expected)


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
        eq_map = _eq_map_for(sn_mesh, ng=ng)

        rng = np.random.default_rng(0)
        psi_vec = rng.standard_normal(eq_map.n_unknowns).astype(np.float64)
        sigma_full = np.full((ng, sn_mesh.nx, sn_mesh.ny), 2.0)  # PR-INDEX-3
        sigma_zero = np.zeros((ng, sn_mesh.nx, sn_mesh.ny))

        # Resolution A's L.apply (subtractive).
        L = StreamingOperator(sn_mesh, sigma_full)
        l_correct = L.apply(psi_vec)

        # Prior agent's wrong L.apply: matvec at σ_t = 0.
        l_prior = _call_matvec(sn_mesh, psi_vec, sigma_zero, eq_map, ng)

        diff = l_correct - l_prior
        rel = (
            np.linalg.norm(diff) / max(np.linalg.norm(l_correct), 1e-300)
        )
        assert rel > 1e-3, (
            f"{geometry}: Resolution A's L.apply and prior agent's "
            f"L.apply differ only by {rel:.3e} — expected O(0.01..1) "
            f"difference because the prior approach degenerates the "
            f"Carlson seed denominator ``dr·σ_t + 2 → 2``."
        )
