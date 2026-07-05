r"""L0 foundation: operator ``apply(TimedFullField) → TimedFullField`` (D-H.2).

The four SN operators in the algebra :math:`(L + C - S - F)` each
expose a typed apply surface so the algebra reads as the math:

    (L + C - S - F).apply(state)  →  TimedFullField

Per-operator boundary semantics pinned here:

* :class:`~orpheus.sn.operators.streaming.StreamingOperator` (``L``) — the ONLY
  operator that emits a non-zero face residual into
  ``result.boundary``.  The B1'' face equation
  :math:`(WDD\text{-propagated face}) − \psi_{\rm face}` lives in
  the appropriate face slots of ``result.boundary``.
* :class:`~orpheus.sn.operators.streaming.CollisionOperator` (``C``) — volumetric;
  ``result.boundary`` is the implicit-zero L2 AngularBoundaryFlux.
* :class:`~orpheus.transport.operators.scattering.ScatteringOperator` (``S``) —
  volumetric secondary-emission; ``result.boundary`` is zero.
* :class:`~orpheus.sn.fission.FissionOperator` (``F``) — volumetric
  rank-1 emission; ``result.boundary`` is zero.

D-H.2-C1 (2026-05-28) — migrated from the legacy
:class:`orpheus.sn.angular_flux.AngularFlux` input contract to the
composite :class:`~orpheus.transport.timed_full_field.TimedFullField`
carrier.  The legacy class retires in D-H.2-C5; the typed-vs-packed
``to_flat_with_traces`` parity tests retired with the migration
(both branches share the same compute kernel — exercising the
composite branch alone is sufficient).  The composite-vs-legacy
parity tests also retired (tautological once both sides are
composite).
"""
from __future__ import annotations

from dataclasses import replace
from typing import assert_type, cast

import numpy as np
import pytest

from orpheus.geometry import BC, Mesh1D, Region, RegionMesh, StructuredGeometry
from orpheus.transport.operators.fission import FissionOperator
from orpheus.sn.mesh.augmented_mesh import SNMesh
from orpheus.sn.operators.streaming import StreamingOperator
from orpheus.transport.operators.multiplication_operator import MultiplicationOperator
from orpheus.numerics.quadrature import Quadrature
from orpheus.transport.fields.angular_flux import AngularFlux
from orpheus.transport.fields.harmonic_moment_flux import HarmonicMomentFlux
from orpheus.transport.source_sinks import AngularSourceSink, ScalarSourceSink
from orpheus.transport.fields.scalar_flux import ScalarFlux
from orpheus.transport.full_field import FullField
from orpheus.transport.timed_full_field import TimedFullField
from orpheus.transport.operators.scattering import ScatteringOperator
from orpheus.numerics.green_operator import GreenOperator
from orpheus.numerics.iteration import SupportsSeededApply
from orpheus.numerics.matrix_inverse_operator import MatrixInverseOperator
from orpheus.numerics.operator import (
    DiagonalOperator,
    IdentityOperator,
    InverseOperator,
    LinearOperator,
    OperatorProduct,
    OperatorSum,
    PermutationOperator,
    ScaledOperator,
    SupportsAdjoint,
    SupportsInverse,
    TensorProductOperator,
    adjointable,
    invertible,
)
from orpheus.sn.operators.sweep_operator import SweepOperator
from orpheus.sn.operators.windowing import BulkAnalysisOperator, WindowedSweep
from tests.sn._test_helpers import placeholder_materials
from orpheus.transport.fields.angular_boundary_flux import AngularBoundaryFlux
from orpheus.transport.fields.starting_direction_flux import StartingDirectionFlux


pytestmark = [pytest.mark.foundation]


# ───────────────────────────────────────────────────────────────────────
# Geometry fixtures (slab + sphere — the two B1'' face geometries)
# ───────────────────────────────────────────────────────────────────────


def _slab_mesh(nx: int = 5, n_ord: int = 4, ng: int = 2) -> SNMesh:
    geom = StructuredGeometry(
        geometry="SLB",
        regions=(Region(mat_id=0, outer_thickness_cm=2.0),),
        bcs=(BC.reflective, BC.reflective),
    )
    mesh = Mesh1D.from_geometry(geom, region_meshes=(RegionMesh(n_cells=nx),))
    quad = Quadrature.gauss_legendre(n_ordinates=n_ord)
    return SNMesh(mesh, quad, placeholder_materials(ng=ng))


def _sphere_mesh(nx: int = 5, n_ord: int = 4, ng: int = 2) -> SNMesh:
    geom = StructuredGeometry(
        geometry="SPH",
        regions=(Region(mat_id=0, outer_thickness_cm=2.0),),
        bcs=(BC.reflective,),
    )
    mesh = Mesh1D.from_geometry(geom, region_meshes=(RegionMesh(n_cells=nx),))
    quad = Quadrature.gauss_legendre(n_ordinates=n_ord)
    return SNMesh(mesh, quad, placeholder_materials(ng=ng))


GEOMETRIES = [("slab", _slab_mesh), ("sphere", _sphere_mesh)]


def _random_state(sn: SNMesh, seed: int) -> TimedFullField:
    """Build a :class:`TimedFullField` with non-trivial bulk + boundary.

    Bulk values + boundary face values are independently sampled from
    a standard-normal RNG seeded with ``seed`` / ``seed+1000``; the
    composite's ``history_depth`` stays at the default (2).
    """
    rng_bulk = np.random.default_rng(seed)
    rng_bnd = np.random.default_rng(seed + 1000)
    state = TimedFullField.zeros(
        bulk=AngularFlux, boundary=AngularBoundaryFlux, mesh=sn,
        starting_direction=StartingDirectionFlux,
    )
    bulk_values = rng_bulk.standard_normal(state.bulk.values.shape)
    boundary_values = rng_bnd.standard_normal(state.boundary.values.shape)
    state = replace(state, bulk=replace(state.bulk, values=bulk_values))
    return replace(
        state, boundary=replace(state.boundary, values=boundary_values),
    )


# ───────────────────────────────────────────────────────────────────────
# F — FissionOperator.apply(TimedFullField) lift
# ───────────────────────────────────────────────────────────────────────


@pytest.mark.parametrize("name,builder", GEOMETRIES)
def test_F_apply_timed_full_field_returns_composite(name, builder) -> None:
    """F.apply(TimedFullField) returns TimedFullField with zero boundary."""
    sn = builder()
    state = _random_state(sn, seed=1)
    F = FissionOperator.from_solver_data(mat_xs=sn.material_xs_field())

    Fpsi = F.apply(state)
    assert isinstance(Fpsi, FullField)  # #257 S8a: timeless codomain (base arrow)
    assert isinstance(Fpsi.bulk, AngularSourceSink)
    assert Fpsi.bulk.values.shape == state.bulk.values.shape
    assert Fpsi.bulk.mesh is sn
    # F is volumetric; result's boundary is implicit-zero.
    np.testing.assert_array_equal(Fpsi.boundary.values, 0.0)


@pytest.mark.parametrize("name,builder", GEOMETRIES)
def test_F_typed_lift_equivalent_to_scalar(name, builder) -> None:
    """F.apply(TimedFullField) broadcasts F.apply(integrate_angular(bulk))."""
    sn = builder()
    state = _random_state(sn, seed=2)
    F = FissionOperator.from_solver_data(mat_xs=sn.material_xs_field())

    # Composite path: F(state) — TimedFullField
    Fpsi_typed = F.apply(state)
    # Scalar path: F(integrate_angular(state.bulk)) — ScalarFlux, then broadcast.
    # ``state.bulk`` is statically the broad ``BulkField`` (``_random_state``
    # builds it as the concrete ``AngularFlux``), and ``integrate_angular`` is
    # itself under-typed (returns ``object``); the casts name both runtime
    # truths so the typed ``F.apply`` overload resolves.
    phi = cast(
        ScalarFlux, cast(AngularFlux, state.bulk).integrate_angular(),
    )
    F_phi = F.apply(phi)
    expected_values = np.broadcast_to(
        F_phi.values[None], Fpsi_typed.bulk.values.shape,
    )
    np.testing.assert_allclose(
        Fpsi_typed.bulk.values, expected_values, rtol=1e-14,
    )


# ───────────────────────────────────────────────────────────────────────
# S — ScatteringOperator.apply(TimedFullField)
# ───────────────────────────────────────────────────────────────────────


@pytest.mark.parametrize("name,builder", GEOMETRIES)
def test_S_apply_timed_full_field_zero_boundary(name, builder) -> None:
    """S.apply(TimedFullField) returns TimedFullField; boundary stays zero."""
    sn = builder()
    state = _random_state(sn, seed=3)
    S = ScatteringOperator.from_solver_data(
        mat_xs=sn.material_xs_field(),
        quadrature=sn.quad,
        scattering_order=0,
    )

    Spsi = S.apply(state)
    assert isinstance(Spsi, FullField)  # #257 S8a: timeless codomain (base arrow)
    assert isinstance(Spsi.bulk, AngularSourceSink)
    assert Spsi.bulk.values.shape == state.bulk.values.shape
    # S is volumetric; result's boundary is implicit-zero.
    np.testing.assert_array_equal(Spsi.boundary.values, 0.0)


# ───────────────────────────────────────────────────────────────────────
# C — CollisionOperator.apply(TimedFullField)
# ───────────────────────────────────────────────────────────────────────


@pytest.mark.parametrize("name,builder", GEOMETRIES)
def test_C_apply_timed_full_field_zero_boundary(name, builder) -> None:
    """C.apply(TimedFullField) returns TimedFullField; boundary stays zero."""
    sn = builder()
    state = _random_state(sn, seed=5)
    sigma = np.ones((sn.ng, *sn.spatial_shape)) * 0.7
    C = MultiplicationOperator.from_mesh(sigma, sn)

    Cpsi = C.apply(state)
    assert isinstance(Cpsi, FullField)  # #257 S8a: timeless codomain (base arrow)
    assert isinstance(Cpsi.bulk, AngularSourceSink)
    assert Cpsi.bulk.values.shape == state.bulk.values.shape
    np.testing.assert_array_equal(Cpsi.boundary.values, 0.0)


@pytest.mark.parametrize("name,builder", GEOMETRIES)
def test_C_diagonal_action(name, builder) -> None:
    """C.apply(TimedFullField).bulk.values == σ ⊙ bulk.values element-wise."""
    sn = builder()
    state = _random_state(sn, seed=6)
    rng = np.random.default_rng(7)
    sigma = 0.3 + 0.5 * rng.random((sn.ng, *sn.spatial_shape))
    C = MultiplicationOperator.from_mesh(sigma, sn)

    Cpsi = C.apply(state)
    expected = sigma[None] * state.bulk.values
    np.testing.assert_array_equal(Cpsi.bulk.values, expected)


# ───────────────────────────────────────────────────────────────────────
# L — StreamingOperator.apply(TimedFullField) — load-bearing
# ───────────────────────────────────────────────────────────────────────


@pytest.mark.parametrize("name,builder", GEOMETRIES)
def test_L_apply_timed_full_field_returns_composite(name, builder) -> None:
    """L.apply(TimedFullField) returns TimedFullField; .boundary may be non-zero."""
    sn = builder()
    state = _random_state(sn, seed=10)
    rng = np.random.default_rng(11)
    sigma_t = 0.4 + 0.4 * rng.random((sn.ng, *sn.spatial_shape))
    L = StreamingOperator(sn)

    Lpsi = L.apply(state)
    assert isinstance(Lpsi, FullField)  # #257 S8a: timeless codomain (base arrow)
    assert isinstance(Lpsi.bulk, AngularSourceSink)
    assert Lpsi.bulk.values.shape == state.bulk.values.shape
    # L emits the face residual into .boundary; for these random
    # inputs at least one face slot is non-zero.
    assert not np.allclose(Lpsi.boundary.values, 0.0)


# ───────────────────────────────────────────────────────────────────────
# Compose (L + C - S - F).apply(state) under TimedFullField
# ───────────────────────────────────────────────────────────────────────


@pytest.mark.parametrize("name,builder", GEOMETRIES)
def test_full_algebra_returns_timed_full_field(name, builder) -> None:
    """``(L + C - S - F).apply(state)`` returns a TimedFullField end-to-end.

    Load-bearing for the D-H.1b operator-algebra type-uniformity claim:
    every leaf operator (L, C, S, F) accepts :class:`TimedFullField`
    input and returns :class:`TimedFullField` output, so the composed
    :class:`~orpheus.numerics.operator.OperatorSum` tree's inherited
    ``apply`` (``a.apply(x) + b.apply(x)``) type-checks at every node.

    ``L + C`` dispatches to :class:`InvertibleOperator` (subclass of
    OperatorSum); subsequent ``-S - F`` chain through the generic
    LinearOperator's ``__sub__`` → ``OperatorSum(self,
    ScaledOperator(other, -1))``.  The ``-1`` scaling propagates via
    :meth:`TimedFullField.__mul__` (which propagates to bulk +
    boundary members).  All four leaves carry their TimedFullField
    branches as of D-H.1b.3..6.
    """
    sn = builder()
    state = _random_state(sn, seed=24)
    sigma_t = np.full((sn.ng, *sn.spatial_shape), 0.7)
    L = StreamingOperator(sn)
    C = MultiplicationOperator.from_mesh(sigma_t * 0.5, sn)
    S = ScatteringOperator.from_solver_data(
        mat_xs=sn.material_xs_field(),
        quadrature=sn.quad,
        scattering_order=0,
    )
    F = FissionOperator.from_solver_data(mat_xs=sn.material_xs_field())

    A = L + C - S - F  # full operator-algebra composition
    out = A.apply(state)

    # #257 S8a — every matvec leaf is a base arrow ``FullField -> FullField``,
    # so the composed algebra output is the TIMELESS FullField (history-free);
    # the comonad lives on the iteration driver, not the operator.
    assert isinstance(out, FullField)
    assert not isinstance(out, TimedFullField)
    assert out.bulk.mesh is sn


@pytest.mark.parametrize("name,builder", GEOMETRIES)
def test_full_algebra_linearity(name, builder) -> None:
    """``A.apply(α·d₁ + β·d₂) == α·A·d₁ + β·A·d₂`` on flux displacements.

    Linearity of the composed operator on the composite carrier — the
    ``+`` / scalar-``*`` dunders propagate through bulk + boundary, and the
    operator algebra commutes with them (linear-operator-on-vector-space
    invariant). #208: the linear combination is formed on flux DISPLACEMENTS
    (the difference vector space V), since ``α·ψ + β·ψ'`` with ``α+β ≠ 1`` is
    illegal on affine flux STATES (no origin) but legal on displacements.
    """
    sn = builder()
    state1 = _random_state(sn, seed=27)
    state2 = _random_state(sn, seed=28)
    sigma_t = np.full((sn.ng, *sn.spatial_shape), 0.7)
    L = StreamingOperator(sn)
    C = MultiplicationOperator.from_mesh(sigma_t * 0.5, sn)
    S = ScatteringOperator.from_solver_data(
        mat_xs=sn.material_xs_field(),
        quadrature=sn.quad,
        scattering_order=0,
    )
    F = FissionOperator.from_solver_data(mat_xs=sn.material_xs_field())
    A = L + C - S - F
    # #208: a general α·ψ₁ + β·ψ₂ (α+β≠1) is illegal on affine flux STATES;
    # verify linearity via the affine-supported ops — homogeneity
    # A(c·ψ)=c·A(ψ) AND affine additivity in torsor form A(ψ₁ + λ(ψ₂⊖ψ₁)) =
    # (1−λ)A(ψ₁) + λA(ψ₂). Together they imply full linearity; A.apply stays on
    # flux states (S, a summand of A, guards its flux-state domain).
    c, lam = 2.5, 0.7

    hom_lhs = A.apply(c * state1)
    hom_rhs = c * A.apply(state1)
    np.testing.assert_allclose(
        hom_lhs.bulk.values, hom_rhs.bulk.values, rtol=1e-12, atol=1e-13)
    np.testing.assert_allclose(
        hom_lhs.boundary.values, hom_rhs.boundary.values, rtol=1e-12, atol=1e-13)

    lhs = A.apply(state1 + lam * (state2 - state1))   # (1−λ)ψ₁ + λψ₂, a flux
    rhs = (1.0 - lam) * A.apply(state1) + lam * A.apply(state2)
    # Bulk linearity check.
    np.testing.assert_allclose(
        lhs.bulk.values, rhs.bulk.values, rtol=1e-12, atol=1e-13,
    )
    # Boundary linearity check.
    np.testing.assert_allclose(
        lhs.boundary.values, rhs.boundary.values,
        rtol=1e-12, atol=1e-13,
    )


# ───────────────────────────────────────────────────────────────────────
# C6 — #257 S8c: the @singledispatchmethod fibration is honestly typed.
#
# ``FissionOperator`` / ``ScatteringOperator`` are NOT endomorphisms ``V -> V``
# (the ``LinearOperator`` nominal contract): their ``apply`` maps each
# input *carrier* to a DISTINCT output carrier.  S8c made that honest with an
# ``@overload`` surface over the runtime ``@singledispatchmethod`` (the public
# ``apply`` is an alias of the private ``_apply_impl`` dispatcher).  Two halves:
#   * STATIC  — ``_c6_static_typing_pins`` below (pyright-checked, never run):
#     reddens if the overload surface regresses to the dispatcher's untyped
#     ``Any``/``NoReturn`` fallback.
#   * RUNTIME — ``test_c6_apply_dispatch_parity`` (Mode-11): exercises the
#     aliased PUBLIC ``apply`` and pins each arm's output TYPE.
# ───────────────────────────────────────────────────────────────────────


def _c6_static_typing_pins(
    F: FissionOperator,
    S: ScatteringOperator,
    state: TimedFullField,
    phi: ScalarFlux,
    psi: AngularFlux,
    moments: HarmonicMomentFlux,
    arr: np.ndarray,
) -> None:
    """Static typing pins for the ``apply`` fibration (pyright-only, never run).

    ``assert_type`` is the STATIC half of gate C6 — each call forces pyright to
    confirm the per-carrier ``@overload`` resolves to the right output type.
    No ``test_`` prefix → pytest never collects this; only the type checker
    reads it.
    """
    assert_type(F.apply(state), FullField)
    assert_type(F.apply(phi), ScalarSourceSink)
    assert_type(F.apply(arr), np.ndarray)
    assert_type(S.apply(state), FullField)
    assert_type(S.apply(phi), ScalarSourceSink)
    assert_type(S.apply(psi), AngularSourceSink)
    assert_type(S.apply(moments), AngularSourceSink)


def _windowed_product_static_typing_pins(
    product: "WindowedSweep",
    p: "BulkAnalysisOperator",
    rhs: TimedFullField,
) -> None:
    """Static pins for the windowed composition ``P @ A.inverse()``
    (#226 step 2, spec §13.4 factor 3; pyright-only, never run).

    The composition's codomain carrier is the timed composite — the
    ``FullField → moment-bulk`` refinement is NOT statically expressible
    today (``FullField`` is not generic over its bulk type), so the
    moment-bulk half of the codomain fact is pinned at RUNTIME in the
    deforestation gate (``test_2d_windowed_product_equals_post_projection``
    asserts ``isinstance(out.bulk, HarmonicMomentFlux)``); minting a
    bulk-generic composite is a separate typing carve, not smuggled here.
    """
    assert_type(product.apply(rhs), TimedFullField)
    assert_type(p.apply(rhs), TimedFullField)


def test_c6_apply_dispatch_parity() -> None:
    """#257 S8c: the aliased public ``apply`` dispatches each carrier to its arm.

    Runtime half of gate C6.  Mode-11: calls the PUBLIC ``apply`` (the
    ``apply = _apply_impl`` alias), exercising the runtime dispatcher — not
    ``_apply_impl`` directly.  Mode-8: ``pytest.fail`` (a function call), not a
    bare ``assert`` (which ``-O`` strips).  Value-level bit-identity of each arm
    lives in ``test_fission_operator.py`` / ``test_scattering_operator.py``;
    here we pin the carrier → output-TYPE dispatch contract.
    """
    sn = _slab_mesh()
    state = _random_state(sn, seed=57)
    F = FissionOperator.from_solver_data(mat_xs=sn.material_xs_field())
    S = ScatteringOperator.from_solver_data(
        mat_xs=sn.material_xs_field(), quadrature=sn.quad, scattering_order=0,
    )
    # Carriers built directly (independent of the under-typed integrate_angular).
    psi = cast(AngularFlux, state.bulk)
    phi = ScalarFlux.from_mesh(np.ones((sn.ng, *sn.spatial_shape)), sn)

    cases = [
        ("F(TimedFullField)", F.apply(state), FullField),
        ("F(ScalarFlux)", F.apply(phi), ScalarSourceSink),
        ("F(ndarray)", F.apply(phi.values), np.ndarray),
        ("S(TimedFullField)", S.apply(state), FullField),
        ("S(ScalarFlux)", S.apply(phi), ScalarSourceSink),
        ("S(AngularFlux)", S.apply(psi), AngularSourceSink),
    ]
    for label, out, expected in cases:
        if not isinstance(out, expected):
            pytest.fail(
                f"{label}: dispatch returned {type(out).__name__}, "
                f"expected {expected.__name__}"
            )


# ───────────────────────────────────────────────────────────────────────
# The inverse family's SEEDED-APPLY static contract (#285 → STRUCTURAL;
# taxonomy §12 step 4, verification spec §19.2).  The abstract
# ``InverseWrapMixin.apply(x, /, *, initial_guess=None)`` is what makes a
# sibling UNABLE to forget the keyword: pyright rejects a kwarg-less
# override (reportIncompatibleMethodOverride — M-MIXIN-KWARG teeth), and
# these pins make the family's conformance a checked fact rather than a
# docstring claim.  Pyright-only, never run (no ``test_`` prefix).
# ───────────────────────────────────────────────────────────────────────


def _inverse_family_seeded_apply_static_pins(
    sweep: "SweepOperator",
    exact: "InverseOperator",
    green: "GreenOperator",
    matrix: "MatrixInverseOperator",
    rhs: TimedFullField,
    seed: TimedFullField,
    arr: np.ndarray,
) -> None:
    """Every wrap-delegate sibling accepts the canonical seeded-apply
    signature AND structurally satisfies the driver's
    :class:`~orpheus.numerics.iteration.SupportsSeededApply` contract —
    so ``SourceIteration`` can consume any of them with no per-type
    probes (the step-3 design, now mixin-enforced; the 4th sibling
    ``MatrixInverseOperator`` joined at taxonomy §12 step 5 —
    M-MINV-KWARG: dropping its kwarg is a pyright
    reportIncompatibleMethodOverride vs the mixin abstract)."""
    assert_type(sweep.apply(rhs, initial_guess=seed), TimedFullField)
    _ = exact.apply(arr, initial_guess=arr)
    _ = green.apply(arr, initial_guess=arr)
    _ = matrix.apply(arr, initial_guess=arr)
    a: SupportsSeededApply = sweep
    b: SupportsSeededApply = exact
    c: SupportsSeededApply = green
    d: SupportsSeededApply = matrix
    del a, b, c, d


# ───────────────────────────────────────────────────────────────────────
# The TypeGuard-bridge NARROWING contract (carve P4, taxonomy §12 step 6,
# verification spec §39.3).  Design C: the runtime predicate check IS the
# static permission — ``invertible(op)`` / ``adjointable(op)`` narrow a
# ``LinearOperator``-typed operand to ``SupportsInverse`` /
# ``SupportsAdjoint`` (which EXTEND ``LinearOperator``, so the guarded
# branch keeps the whole algebra), and the member access type-checks with
# NO cast.  Deleting a production guard un-narrows its call site — CLI
# pyright REDs (M-GUARD-DELETE-PYRIGHT); widening a bridge's return
# annotation to plain ``bool`` reddens THESE pins (M-BRIDGE-ANNOT).
# Pyright-only, never run (no ``test_`` prefix).  The structural-absent
# converse — ``ZeroOperator().inverse()`` does NOT type-check — lives in
# ``_static_zero_inverse_unspellable.py`` (a static pin cannot share a
# module with the expression it forbids).
# ───────────────────────────────────────────────────────────────────────


def _typeguard_bridge_narrowing_static_pins(
    endo: LinearOperator[TimedFullField, TimedFullField],
    two_space: LinearOperator[ScalarFlux, AngularFlux],
    b_endo: TimedFullField,
    b_cod: AngularFlux,
) -> None:
    """The Design-C checked bridge, pinned as a static fact — no ``cast`` anywhere.

    Three legs per axis: (1) the guarded member access RESOLVES (the
    narrowing exists); (2) the narrowed signature carries the honest
    carrier SWAP (``inverse(): [D, C] → [C, D]``; ``apply_transpose:
    Codomain → Domain``); (3) the narrowed operand still IS a
    ``LinearOperator`` — the Protocol-promotion payoff: ``TypeGuard``
    REPLACES the declared type, so a bare-Protocol target would have
    LOST the algebra inside the guard (spec §44.E).
    """
    if invertible(endo):
        narrowed_inv: SupportsInverse[TimedFullField, TimedFullField] = endo
        assert_type(endo.inverse(), LinearOperator[TimedFullField, TimedFullField])
        assert_type(endo.apply(b_endo), TimedFullField)  # leg 3: algebra retained
        del narrowed_inv
    if adjointable(endo):
        narrowed_adj: SupportsAdjoint[TimedFullField, TimedFullField] = endo
        assert_type(endo.apply_transpose(b_endo), TimedFullField)
        assert_type(endo.apply(b_endo), TimedFullField)  # leg 3: algebra retained
        del narrowed_adj
    # The two-space operand pins leg 2 — the carrier swap — where domain
    # and codomain are DISTINCT (an endomorphism cannot see a swap bug).
    if invertible(two_space):
        assert_type(two_space.inverse(), LinearOperator[AngularFlux, ScalarFlux])
    if adjointable(two_space):
        assert_type(two_space.apply_transpose(b_cod), ScalarFlux)


# ───────────────────────────────────────────────────────────────────────
# The COMPOSITION-ALGEBRA RETURN-TYPE laws (carve P5, taxonomy §12
# step-6 tail; verification spec §9 "P5" + §2.2).  Every composition
# surface of the algebra declares the PRECISE structure it returns —
# ``+`` an :class:`OperatorSum`, ``*`` a :class:`ScaledOperator`, ``@``
# an :class:`OperatorProduct` with the intermediate carrier solved
# away, ``.H`` the carrier SWAP — and every ``.inverse()`` declares
# which of the TWO inverse kinds it returns (the canonical statement
# lives on :meth:`OperatorProduct.inverse`): ALGEBRA-CLOSED structures
# invert into themselves as first-class forwards (scaled → scaled on
# swapped carriers, identity/permutation/tensor → themselves), while
# solve-backed structures return a WRAP-DELEGATE family member
# (leaves and products → the generic :class:`InverseOperator`; a
# generic sum → the deliberately-LOOSE family face
# ``LinearOperator[Codomain, Domain]``, because its subclass overrides
# return sibling kinds — Sweep/Green — not a hierarchy).  These pins
# make each law a pyright-checked fact; un-parametrizing a composed
# return annotation (e.g. ``ScaledOperator.inverse() ->
# "ScaledOperator"``, dropping the swap — M-SCALED-BARE) reddens the
# matching pin via reportAssertTypeFailure.  Pyright-only, never run
# (no ``test_`` prefix).  NOTE the honest boundary: pins parametrize
# only surfaces the algebra PROMISES today — the bare-generic leaves
# (Diagonal/Permutation/TensorProduct) and the wrap-family's
# ``apply`` stay unparametrized until the leaf-parametrization
# campaign, and ``A.H.inverse()`` has NO pin because
# ``_AdjointOperator`` deliberately has no ``inverse()`` yet (the
# swap law ``(A.H)⁻¹ = (A⁻¹).H`` is #280's deliverable — extend this
# bank with its pin when #280 lands).
# ───────────────────────────────────────────────────────────────────────


def _composition_algebra_return_type_static_pins(
    a: LinearOperator[ScalarFlux, AngularFlux],
    a2: LinearOperator[ScalarFlux, AngularFlux],
    b: LinearOperator[TimedFullField, ScalarFlux],
    endo: LinearOperator[TimedFullField, TimedFullField],
    scaled: ScaledOperator[ScalarFlux, AngularFlux],
    product: OperatorProduct[TimedFullField, AngularFlux],
    a_sum: OperatorSum[ScalarFlux, AngularFlux],
    identity: IdentityOperator[ScalarFlux],
    perm: PermutationOperator,
    tensor: TensorProductOperator,
    diag: DiagonalOperator,
    green: GreenOperator,
    matrix_inv: MatrixInverseOperator,
) -> None:
    """The composition-algebra return-type laws, pinned as static facts.

    Two-space operands (``ScalarFlux → AngularFlux``) are used wherever a
    swap could hide, so a parameter-order bug is VISIBLE; the
    endomorphic/bare operands (``endo``, ``identity``, ``perm``,
    ``tensor``, …) carry no swap to see.
    """
    # ── the composition constructors carry the honest parametrization ──
    assert_type(a + a2, OperatorSum[ScalarFlux, AngularFlux])
    assert_type(a - a2, OperatorSum[ScalarFlux, AngularFlux])
    assert_type(2.0 * a, ScaledOperator[ScalarFlux, AngularFlux])
    assert_type(-a, ScaledOperator[ScalarFlux, AngularFlux])
    assert_type(a / 2.0, ScaledOperator[ScalarFlux, AngularFlux])
    # ``a @ b``: the intermediate carrier (ScalarFlux) is solved away —
    # the product is honestly ``TimedFullField → AngularFlux``.
    assert_type(a @ b, OperatorProduct[TimedFullField, AngularFlux])
    # ``.H`` / ``.adjoint()``: the carrier SWAP, at the base's honest face.
    assert_type(a.H, LinearOperator[AngularFlux, ScalarFlux])
    assert_type(a.adjoint(), LinearOperator[AngularFlux, ScalarFlux])
    # ``A**n`` is endomorphism-only and stays on the one carrier.
    assert_type(endo**3, LinearOperator[TimedFullField, TimedFullField])

    # ── the inverse return-type laws: which KIND, on swapped carriers ──
    # ALGEBRA-CLOSED: the inverse is a first-class forward of the same
    # structure — the scaled inverse carries the swap in its parameters.
    assert_type(scaled.inverse(), ScaledOperator[AngularFlux, ScalarFlux])
    assert_type(identity.inverse(), IdentityOperator[ScalarFlux])
    assert_type(perm.inverse(), PermutationOperator)
    assert_type(tensor.inverse(), TensorProductOperator)
    # WRAP-DELEGATE: solve-backed structures return a family member.
    assert_type(product.inverse(), InverseOperator)
    assert_type(diag.inverse(), InverseOperator)
    # A generic sum's inverse is DRIVER-realized; the annotation is the
    # deliberately-loose family face (its subclass overrides return
    # sibling kinds — SweepOperator, GreenOperator — not a hierarchy),
    # so the pin asserts exactly that loose face, on swapped carriers.
    assert_type(a_sum.inverse(), LinearOperator[AngularFlux, ScalarFlux])

    # ── the involution is typed per sibling through the mixin's
    # ``_ForwardT`` — ``(A⁻¹)⁻¹`` returns the wrapped FORWARD's type
    # (object identity, taxonomy §13 I2), not a re-wrapped inverse. ──
    assert_type(green.inverse(), OperatorSum)
    assert_type(matrix_inv.inverse(), LinearOperator)
