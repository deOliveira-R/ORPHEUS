r"""Pure-L + C decomposition gate (Issue #196 Phase G Step 3+4.b.i; #257 S8b).

Promoted from :file:`derivations/diagnostics/diag_LC_decomposition_resolution.py`
— this is the load-bearing verification that the σ-free pure-``L``
``StreamingOperator.apply`` composes with ``C = M[σ_t]`` to recover the
full within-group loss ``(L + C).apply ≡ M`` on every geometry.

The math (#257 S8b — pure-L)
----------------------------

.. math::

   L.{\rm apply}(\psi) \;:=\; \Omega\cdot\nabla\psi
       \;=\; \text{streaming\_action}(\psi)
   \\
   C.{\rm apply}(\psi) \;:=\; \sigma_t \odot \psi_{\rm bulk}
   \\
   (L + C).{\rm apply}(\psi) \;\equiv\; M(\psi;\;\sigma_t)
       \qquad \text{rel\_residual} = 0.0

``L`` is now pure σ-free streaming (#257 S8b — the ``(L+C)−C`` fold is
retired): it reads NO σ and calls the loss representation's named
``streaming_action`` leaf directly.  ``C = M[σ_t]`` is the separate
shared multiplier.  The within-group WDD matvec is AFFINE in σ in the
forward direction — ``streaming_action(ψ) = loss_action(0, ψ)`` and
``M(ψ; σ_t) = streaming_action(ψ) + σ_t⊙ψ`` — so the composition still
recovers ``M`` bit-exactly (the composite ``InvertibleOperator.apply``
calls ``loss_action(σ_t)`` directly, unchanged from pre-S8b).

The affinity is genuine since ERR-058 (#195) made the curvilinear
Carlson coupled-pole seed σ-independent (``AngularEdgeExtrapolation``):
``loss_action(0, ψ)`` IS pure streaming (the seed denominator
``Δr·σ_t + 2`` no longer degenerates — its σ-contribution is exactly
the collision diagonal it injects, which cancels into σ⊙ψ).  Through
the legacy ``CarlsonInwardSweep`` era the σ=0 matvec WAS a different
operator (the reverted ``ad37ca0`` approach); ERR-058 removed that
coupling, which is what licenses the pure-L carve.

Test contract (mechanism criteria 1, 2 from the substep brief)
---------------------------------------------------------------

* slab/sphere/cylinder × 3 random seeds × σ_t = 2.0 uniform.
* Composition residual ``(L+C).apply(ψ) − M(ψ;σ_t)`` MUST be bit-exact
  (``rel_residual < 1e-14``) on BOTH the cell-centre bulk
  ``AngularFlux`` and the boundary face residual ``AngularBoundaryFlux``.
* No xfail.  ``(L+C).apply`` is the composite's ``loss_action(σ_t)``,
  byte-identical to ``M``.

A separate test verifies the pure-``L`` apply matches the affine
subtractive form ``M − σ_t⊙ψ`` to FP-non-associativity ULP (the named
``streaming_action`` leaf re-associates the reduction tree vs the old
``(L+C)−C`` fold).

D-I.1 — typed-carrier migration
-------------------------------

This file was originally written against the legacy packed-vector
contract (bare-``np.ndarray`` through
:func:`_transport_operator_matvec_unified` via the
:class:`EquationMap` B1'' adapter).  D-I.1 retires the collision
multiplier's legacy bare-``np.ndarray`` ``apply`` / ``solve`` contract
(the multiplier is now a
:class:`~orpheus.transport.operators.multiplication_operator.MultiplicationOperator`,
``C = M[σ_t]``) along with the
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
from orpheus.sn.mesh.augmented_mesh import SNMesh
from orpheus.sn.operators.streaming import (
    StreamingOperator,
    )
from orpheus.transport.operators.multiplication_operator import MultiplicationOperator
from tests.sn._test_helpers import _LC_matvec
from orpheus.numerics.quadrature import Quadrature
from orpheus.transport.fields.angular_flux import AngularFlux
from orpheus.transport.fields.angular_boundary_flux import AngularBoundaryFlux
from orpheus.transport.timed_full_field import TimedFullField
from tests.sn._test_helpers import placeholder_materials, radial_characteristic_edge_seed

pytestmark = pytest.mark.l0


def _l_apply(L, state, seed_leg, sn_mesh):
    """``L.apply`` through the B.2d explicit ψ½ legs (scratch rows out)."""
    if seed_leg is None:
        return L.apply(state)
    from orpheus.transport.radial_characteristic_field import (
        RadialCharacteristicField,
    )

    return L.apply(
        state,
        radial_characteristic_flux=seed_leg,
        radial_characteristic_source=(
            RadialCharacteristicField.source_zeros_on(sn_mesh)
        ),
    )


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
    ``L.apply := M(ψ; σ_t) − σ_t⊙ψ.interior`` and ``C.apply := σ_t⊙ψ.interior``,
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

        rng = np.random.default_rng(seed)
        bulk_arr = rng.standard_normal((N, ng, *sn_mesh.spatial_shape))
        state = TimedFullField(
            interior=AngularFlux.from_mesh(bulk_arr, sn_mesh),
            boundary=AngularBoundaryFlux.zeros_on(sn_mesh),
            _history=(),
            history_depth=2,
        )
        # #282 route (a) → B.2d: the CONSISTENT edge-extrapolated ψ½ seed on
        # a carrying mesh (SPH) rides the walk's EXPLICIT flux leg in EVERY
        # apply below; None on non-carrying (CART/CYL).  The System-A
        # decomposition identity (L+C ≡ L + C on bulk ⊕ trace) holds for ANY
        # seed since both σ-paths consume the SAME leg (the seed's bulk feed
        # is σ-independent); the ray's own σ term is A_BB's, System B.
        seed_leg = radial_characteristic_edge_seed(bulk_arr, sn_mesh)
        sigma_t = np.full((ng, *sn_mesh.spatial_shape), 2.0)

        # Reference: the unified matvec at full σ_t.
        m_full_state = _LC_matvec(
            state, sigma_t, radial_characteristic_flux=seed_leg,
        )

        # Pure-L + C via TimedFullField arithmetic (#257 S8b): L reads no σ.
        L = StreamingOperator(sn_mesh)
        C = MultiplicationOperator.from_mesh(sigma_t, sn_mesh)
        sum_state = _l_apply(L, state, seed_leg, sn_mesh) + C.apply(state)

        # Bulk residual — (N, ng, nx, ny) ndarray.
        residual_bulk = sum_state.interior.values - m_full_state.interior.values
        rel_bulk = (
            np.linalg.norm(residual_bulk)
            / max(np.linalg.norm(m_full_state.interior.values), 1e-300)
        )
        assert rel_bulk < 1e-14, (
            f"{geometry} seed={seed}: bulk rel_residual={rel_bulk:.3e} "
            f"— Resolution A subtractive decomposition FAILED bit-exact "
            f"gate on bulk. (L + C).apply.interior MUST equal M(ψ; σ_t).interior."
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
# Mechanism criterion 2 — the pure-L apply matches the affine subtractive
# form: L.apply(ψ).interior ≈ M(ψ; σ_t).interior − σ_t⊙ψ.interior to FP ULP (the named
# streaming_action leaf re-associates the reduction tree); boundary STRICT.
# ═══════════════════════════════════════════════════════════════════════


class TestSubtractiveDefinition:
    """``L.apply(ψ).interior ≈ M(ψ; σ_t).interior − σ_t · ψ.interior`` to FP ULP
    and ``L.apply(ψ).boundary == M(ψ; σ_t).boundary`` at STRICT exact equality.

    #257 S8b: pure-``L`` ``apply`` is ``streaming_action(ψ) =
    loss_action(0, ψ)`` — no longer the *defining* equation ``M − σ⊙ψ``
    (that bit-exact subtractive form was the retired ``(L+C)−C`` fold).
    Because the WDD matvec is affine in σ, ``streaming_action(ψ)``
    numerically equals ``M − σ⊙ψ`` to a few ULP (the σ-free walk
    re-associates the FP reduction tree relative to subtracting σ⊙ψ off
    the full σ-bearing matvec).  The BOUNDARY block stays STRICT 0-ULP:
    ``C`` contributes zero on the trace and ``loss_action(0)`` produces
    the identical outflow defect, so the face residual must not move
    (a moved boundary would be a real bug).  Together with the
    ``C.apply == σ_t⊙ψ.interior; zero boundary`` contract (covered by
    test_collision_operator), these two tests pin the pure-L + C split.

    The drift bound is ``REDUCTION_DEPTH × ULP`` (a single-step matvec FP
    re-association — the vv-principles direct-kind convention), measured
    ≤ ~64 ULP across geometries.
    """

    # The σ-free streaming walk re-associates the FP reduction tree vs the
    # σ-bearing matvec minus σ⊙ψ; the drift is bounded by the matvec
    # reduction depth.  64 ULP is comfortably above the measured ≤16 ULP
    # at n_cells=5 and below any real-bug magnitude.
    _BULK_NULP = 256

    @pytest.mark.parametrize("geometry", ["CART", "SPH", "CYL"])
    @pytest.mark.parametrize("seed", [10, 11, 12])
    def test_L_apply_equals_subtractive_form(self, geometry, seed):
        sn_mesh = _build_sn_mesh(geometry, n_cells=5, n_ord=4)
        ng = 1
        N = sn_mesh.quad.N

        rng = np.random.default_rng(seed)
        bulk_arr = rng.standard_normal((N, ng, *sn_mesh.spatial_shape))
        state = TimedFullField(
            interior=AngularFlux.from_mesh(bulk_arr, sn_mesh),
            boundary=AngularBoundaryFlux.zeros_on(sn_mesh),
            _history=(),
            history_depth=2,
        )
        # #282 route (a) → B.2d: the CONSISTENT edge-extrapolated ψ½ seed on
        # a carrying mesh (SPH) rides the walk's EXPLICIT flux leg in EVERY
        # apply below; None on non-carrying (CART/CYL).  The System-A
        # decomposition identity (L+C ≡ L + C on bulk ⊕ trace) holds for ANY
        # seed since both σ-paths consume the SAME leg (the seed's bulk feed
        # is σ-independent); the ray's own σ term is A_BB's, System B.
        seed_leg = radial_characteristic_edge_seed(bulk_arr, sn_mesh)
        sigma_t = np.full((ng, *sn_mesh.spatial_shape), 2.0)

        # Pure-L apply (σ-free, #257 S8b): L takes only the mesh.
        L = StreamingOperator(sn_mesh)
        l_state = _l_apply(L, state, seed_leg, sn_mesh)

        m_full_state = _LC_matvec(
            state, sigma_t, radial_characteristic_flux=seed_leg,
        )

        # Affine subtractive form: bulk expected = M.interior - σ_t·ψ.interior
        # (σ_t broadcast over the ordinate axis 0 via [None, :, :, :]).
        expected_bulk = (
            m_full_state.interior.values - sigma_t[None] * state.interior.values
        )
        # Face slots: L.boundary == M.boundary (no volumetric collision
        # on the trace; the cell-balance σ·ψ term is a CELL quantity).
        expected_boundary = m_full_state.boundary.values

        # Pure-L apply == the affine subtractive form to FP ULP on bulk.
        np.testing.assert_array_almost_equal_nulp(
            l_state.interior.values, expected_bulk, nulp=self._BULK_NULP,
        )
        # Boundary STRICT 0-ULP — the face residual must not move.
        np.testing.assert_array_equal(
            l_state.boundary.values, expected_boundary,
        )


# ═══════════════════════════════════════════════════════════════════════
# Defining identity — pure-L apply IS loss_action(σ=0) (the σ-free leaf).
# ═══════════════════════════════════════════════════════════════════════


class TestPureLIsLossActionAtZeroSigma:
    """Pure-``L`` ``apply`` IS ``loss_action(σ=0)`` byte-for-byte (#257 S8b).

    The defining identity of the σ-free carve: ``StreamingOperator.apply``
    routes to the representation's ``streaming_action``, which is
    ``loss_action(0, ψ)`` (Pattern 2 — the streaming discretization is
    single-sourced; there is no twin σ-free walk).  This test pins that
    identity directly: ``L.apply(ψ) == loss_action(0, ψ)`` bit-exactly,
    on bulk AND boundary, for every geometry.

    HISTORY: through the legacy ``CarlsonInwardSweep`` era this class
    pinned the opposite — that the (then-subtractive) ``L.apply`` DIFFERED
    from ``matvec(σ=0)`` on sphere, because the Carlson seed's denominator
    ``dr·σ_t + 2`` degenerated at σ=0.  ERR-058 (#195, 2026-06-12) flipped
    the default seed to the σ-INDEPENDENT ``AngularEdgeExtrapolation``, so
    ``loss_action(0)`` became genuine pure streaming (measured 1.3e-16
    agreement with the subtractive form).  That σ-freedom is precisely
    what #257 S8b builds on: pure ``L`` now IS ``loss_action(0)`` by
    construction.  If a future seed strategy reintroduces σ-dependence
    into the closure, the affinity breaks and this gate (plus the
    decomposition gate above) flips loudly.
    """

    @pytest.mark.parametrize("geometry", ["CART", "SPH", "CYL"])
    def test_pure_L_is_loss_action_at_zero_sigma(self, geometry):
        sn_mesh = _build_sn_mesh(geometry, n_cells=5, n_ord=4)
        ng = 1
        N = sn_mesh.quad.N

        rng = np.random.default_rng(0)
        bulk_arr = rng.standard_normal((N, ng, *sn_mesh.spatial_shape))
        state = TimedFullField(
            interior=AngularFlux.from_mesh(bulk_arr, sn_mesh),
            boundary=AngularBoundaryFlux.zeros_on(sn_mesh),
            _history=(),
            history_depth=2,
        )
        # #282 route (a) → B.2d: the CONSISTENT edge-extrapolated ψ½ seed on
        # a carrying mesh (SPH) rides the walk's EXPLICIT flux leg in EVERY
        # apply below; None on non-carrying (CART/CYL).  The System-A
        # decomposition identity (L+C ≡ L + C on bulk ⊕ trace) holds for ANY
        # seed since both σ-paths consume the SAME leg (the seed's bulk feed
        # is σ-independent); the ray's own σ term is A_BB's, System B.
        seed_leg = radial_characteristic_edge_seed(bulk_arr, sn_mesh)
        sigma_zero = np.zeros((ng, *sn_mesh.spatial_shape))

        # Pure-L apply (σ-free, #257 S8b): L takes only the mesh.
        L = StreamingOperator(sn_mesh)
        l_state = _l_apply(L, state, seed_leg, sn_mesh)

        # The single-sourced σ-free walk: loss_action at σ = 0 directly
        # (same explicit legs — the SAME call under streaming_action).
        from orpheus.sn.loss_representation import default_for
        if seed_leg is None:
            l_action_zero = default_for(sn_mesh).loss_action(sigma_zero, state)
        else:
            from orpheus.transport.radial_characteristic_field import (
                RadialCharacteristicField,
            )

            l_action_zero = default_for(sn_mesh).loss_action(
                sigma_zero, state,
                radial_characteristic_flux=seed_leg,
                radial_characteristic_source=(
                    RadialCharacteristicField.source_zeros_on(sn_mesh)
                ),
            )

        # Byte-exact: pure-L apply IS loss_action(0) (same call under
        # streaming_action) on BOTH bulk and boundary.
        np.testing.assert_array_equal(
            l_state.interior.values, l_action_zero.interior.values,
        )
        np.testing.assert_array_equal(
            l_state.boundary.values, l_action_zero.boundary.values,
        )
