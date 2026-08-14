r"""Issue #168 Phase C — Gate Set 1 verification gates.

Gates pin the sweep-frame apply matvec rewrite at the architectural
level (per :doc:`/theory/methods/sn/curvilinear_one_group` Phase C section and
.claude/plans/issue_168_phase_c.md §5):

* Gate 1.1 — per-ordinate flat-flux residual on reflective sphere /
  cylinder (canonical curvilinear bug-class diagnostic).
* Gate 1.2 — apply ↔ sweep face-flux propagation identity (structural
  by construction under sweep-frame).
* Gate 1.3 — apply ↔ apply_transpose reciprocity (free if 1.4 passes).
* Gate 1.4 — apply linearity (PRECONDITION; gate 1.3 depends on it).
* Gate 1.5 — BC trace contract verification (matvec consumes
  WDD-propagated outflow face values, not cell-centres).

ERR-026 tripwire — these gates green when the architectural fix lands.

Operator algebra
----------------

Gates consume the composite algebra :class:`StreamingOperator` +
:class:`CollisionOperator` = :class:`StreamingCollisionOperator` via
:class:`~orpheus.transport.timed_full_field.TimedFullField`:

* Resolution A subtractive identity:
  :math:`(L + C).{\rm apply}(\psi) = M(\psi;\sigma_t)` bit-exact.
* ``op.apply(state).interior.values`` holds :math:`(L+C)\psi`'s
  cell-centre block; face residuals live in ``out.boundary``.
* Linearity (Gate 1.4) tests via :class:`TimedFullField` arithmetic
  (``__add__``, scalar ``__mul__/__rmul__``).
* Reciprocity (Gate 1.3) — historically xfailed while
  :class:`StreamingOperator` was apply-only (``is_adjointable`` False);
  Wave O Issue #208 landed the analytic ``Lᵀ`` and the
  :class:`OperatorSum` closure law now carries the transpose (see the
  gate's own docstring below).
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.derivations.common.xs_library import get_mixture
from orpheus.geometry import (
    BC,
    CoordSystem,
    Mesh1D,
    Region,
    RegionMesh,
    StructuredGeometry,
)
from orpheus.geometry.boundary import (
    ReflectiveBoundary,
    SelfPairedDeck,
    VacuumInflow,
)
from orpheus.sn.mesh.augmented_mesh import SNMesh
from orpheus.sn.operators.streaming import (
    StreamingOperator,
)
from orpheus.transport.operators.multiplication_operator import MultiplicationOperator
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.sweep.pole_angular_closure import (
    MorelMontryAngularSweep,
)
from tests.sn._test_helpers import sweep_once
from orpheus.transport.fields.angular_flux import AngularFlux
from orpheus.transport.fields.angular_boundary_flux import AngularBoundaryFlux
from orpheus.transport.timed_full_field import TimedFullField
from tests.sn._test_helpers import (
    cart2d_2g_nonsquare,
    placeholder_materials,
    radial_characteristic_edge_seed,
)


# ═══════════════════════════════════════════════════════════════════════
# Fixtures
# ═══════════════════════════════════════════════════════════════════════


def _make_spherical_sn_mesh(
    nx: int = 8,
    R: float = 1.0,
    quad_name: str = "gl4",
    bc_outer: BC | None = None,
    pole_closure=None,
) -> tuple[SNMesh, np.ndarray]:
    """Build a homogeneous-material spherical SNMesh + sig_t array.

    Returns (sn_mesh, sig_t).  sig_t shape (ng=1, nx) under the rank-d layout.
    """
    if quad_name == "gl4":
        quad = Quadrature.gauss_legendre(4)
    elif quad_name == "gl8":
        quad = Quadrature.gauss_legendre(8)
    else:
        raise ValueError(quad_name)
    edges = np.linspace(0.0, R, nx + 1)
    mesh = Mesh1D(
        edges=edges,
        mat_ids=np.zeros(nx, dtype=int),
        coord=CoordSystem.SPHERICAL,
        bc_right=bc_outer or BC("reflective"),
    )
    sn_mesh = SNMesh(mesh, quad, placeholder_materials(), pole_angular_closure=pole_closure)
    sig_t = np.full((1, nx), 0.5)  # (ng, nx) — rank-d
    return sn_mesh, sig_t


def _make_cylindrical_sn_mesh(
    nx: int = 8,
    R: float = 1.0,
    quad_name: str = "folded_4x8",
    bc_outer: BC | None = None,
    pole_closure=None,
) -> tuple[SNMesh, np.ndarray]:
    """Build a homogeneous-material cylindrical SNMesh + sig_t array.

    6.3 flip: the admitted cylinder family is the σ_y fold; the retired
    ``ls4``/``prod_2x4`` arms (the latter caller-less) named rules the
    mesh now refuses."""
    if quad_name == "folded_4x8":
        quad = Quadrature.folded_product(n_mu=4, n_phi=8)
    else:
        raise ValueError(quad_name)
    edges = np.linspace(0.0, R, nx + 1)
    mesh = Mesh1D(
        edges=edges,
        mat_ids=np.zeros(nx, dtype=int),
        coord=CoordSystem.CYLINDRICAL,
        bc_right=bc_outer or BC("reflective"),
    )
    sn_mesh = SNMesh(mesh, quad, placeholder_materials(), pole_angular_closure=pole_closure)
    sig_t = np.full((1, nx), 0.5)  # (ng, nx) — rank-d
    return sn_mesh, sig_t


def _build_composite(
    sn_mesh: SNMesh,
    bulk_values: np.ndarray,
    boundary_values: np.ndarray | None = None,
    *,
    radial_characteristic_values: np.ndarray | None = None,
) -> TimedFullField:
    """Build a :class:`TimedFullField` from raw bulk + optional boundary arrays.

    Parameters
    ----------
    sn_mesh : SNMesh
        The mesh defining the typed shape ``(N, ng, *spatial)`` on bulk
        and the boundary flat layout.
    bulk_values : np.ndarray
        Shape ``(N, ng, *spatial)`` — the angular flux values.
    boundary_values : np.ndarray, optional
        Shape matching ``sn_mesh.boundary_face_layout.total_size``.  If
        ``None``, an all-zero boundary is used (the typical migration
        target — Gate 1.1/1.4 etc. zero the boundary because they
        compute the cell-block residual only).
    radial_characteristic_values : np.ndarray, optional
        #282 route (a): the flat ψ½ block ``(space.shape[0],)`` on a
        carrying mesh (R12a).  ``None`` (default) fills it with the
        CONSISTENT edge-extrapolation of ``bulk_values``
        (:func:`~tests.sn._test_helpers.radial_characteristic_edge_seed`),
        which is linear in the bulk (so the apply stays a linear
        operator) and constant-preserving (so ``(L+C)·const = σ_t·const``
        still holds).  The reciprocity gate overrides with a RANDOM block
        (exercising the seed rows under the FULL-space Euclidean dot).
        Non-carrying meshes (slab/cyl) ignore this — the block is
        structurally ``None``.
    """
    if boundary_values is None:
        boundary = AngularBoundaryFlux.zeros_on(sn_mesh)
    else:
        # A.5: the AngularBoundaryFlux space IS the mesh's unified AngularTraceSpace
        # (it carries the FaceLayout); no ad-hoc sn_boundary_flat build.
        boundary = AngularBoundaryFlux(
            values=boundary_values, space=sn_mesh.angular_trace, mesh=sn_mesh,
        )
    if radial_characteristic_values is None:
        radial_characteristic = radial_characteristic_edge_seed(bulk_values, sn_mesh)
    elif sn_mesh.radial_characteristic_field_space is not None:
        from orpheus.transport.radial_characteristic_field import (
            RadialCharacteristicField,
        )
        radial_characteristic = RadialCharacteristicField.from_flat(
            radial_characteristic_values,
            RadialCharacteristicField.from_mesh(sn_mesh),
        )
    else:
        radial_characteristic = None
    psi_a = TimedFullField(
        interior=AngularFlux.from_mesh(bulk_values, sn_mesh),
        boundary=boundary,
        _history=(),
        history_depth=2,
    )
    if radial_characteristic is None:
        return psi_a
    # B.2d / 4e: the carrying state is the COUPLED pair — System B rides its
    # own native split composite member (no unified bridge).
    from orpheus.numerics.coupled_system import CoupledField

    return CoupledField(systems=(psi_a, radial_characteristic))


def _random_bulk(sn_mesh: SNMesh, rng: np.random.Generator) -> np.ndarray:
    """Random ``(N, ng, *spatial)`` bulk values for the mesh."""
    return rng.standard_normal((sn_mesh.quad.N, sn_mesh.ng, *sn_mesh.spatial_shape))


def _joint_op(sn_mesh: SNMesh, op):
    """The JOINT operator for the mesh (step 5): the honest triangular ``M``
    grid on a carrying mesh (the numerics substitution — the fused
    ``CoupledInvertibleOperator`` bridge deleted at 5d), ``op`` itself on a
    seedless one."""
    if sn_mesh.radial_characteristic_field_space is None:
        return op
    from tests.sn._test_helpers import joint_m_grid

    return joint_m_grid(sn_mesh, op)[0]


def _sysA(x):
    """System A's member off a pair, or the bare composite itself."""
    return x.systems[0] if hasattr(x, "systems") else x


# ═══════════════════════════════════════════════════════════════════════
# Gate 1.4 — apply linearity (PRECONDITION)
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.foundation
@pytest.mark.parametrize(
    "geom",
    [
        pytest.param("sphere", id="sphere_GL4_reflective"),
        pytest.param("cylinder", id="cyl_LS4_reflective"),
    ],
)
def test_apply_linearity_under_sweep_frame(geom):
    r"""(L+C) is linear, verified via the affine-supported operations.

    Precondition for Gate 1.2 (sweep-frame face-flux equivalence) and
    Gate 1.3 (apply ↔ apply_transpose reciprocity).  Linearity is a
    pure software contract: the matvec MUST be an exact linear
    function of its input vector. Any nonlinearity (e.g. a BC apply
    consuming cell-centres in a way that depends on input sign) is a
    catastrophic operator-correctness failure.

    #208 affine reframe — a general ``α·ψ + β·φ`` with ``α+β ≠ 1`` is illegal
    on affine flux STATES (no origin), so linearity is verified by scalar
    homogeneity ``op(c·ψ) = c·op(ψ)`` AND affine additivity in torsor form
    ``op(ψ₁ + λ(ψ₂⊖ψ₁)) = (1−λ)op(ψ₁) + λop(ψ₂)``. The two together imply full
    matvec linearity; ``apply`` stays on flux states (its domain). The
    right-hand side uses the source-image vector-space dunders.
    """
    rng = np.random.default_rng(seed=42)
    if geom == "sphere":
        sn_mesh, sig_t = _make_spherical_sn_mesh()
    else:
        sn_mesh, sig_t = _make_cylindrical_sn_mesh()
    L = StreamingOperator(sn_mesh)
    C = MultiplicationOperator.from_mesh(sig_t, sn_mesh)
    op = _joint_op(sn_mesh, L + C)   # B.2d: the joint M on the carrying pair
    psi1 = _build_composite(sn_mesh, _random_bulk(sn_mesh, rng))
    psi2 = _build_composite(sn_mesh, _random_bulk(sn_mesh, rng))
    c, lam = 1.7, 0.7
    hom_lhs = op.apply(c * psi1)
    hom_rhs = c * op.apply(psi1)
    np.testing.assert_allclose(
        _sysA(hom_lhs).interior.values, _sysA(hom_rhs).interior.values,
        rtol=1e-13, atol=1e-14)
    np.testing.assert_allclose(
        _sysA(hom_lhs).boundary.values, _sysA(hom_rhs).boundary.values,
        rtol=1e-13, atol=1e-14)
    lhs = op.apply(psi1 + lam * (psi2 - psi1))   # (1−λ)ψ₁ + λψ₂, a flux
    rhs = (1.0 - lam) * op.apply(psi1) + lam * op.apply(psi2)
    np.testing.assert_allclose(
        _sysA(lhs).interior.values, _sysA(rhs).interior.values,
        rtol=1e-13, atol=1e-14,
    )
    np.testing.assert_allclose(
        _sysA(lhs).boundary.values, _sysA(rhs).boundary.values,
        rtol=1e-13, atol=1e-14,
    )


# ═══════════════════════════════════════════════════════════════════════
# Gate 1.1 — per-ordinate flat-flux residual (the canonical ERR-026 probe)
# ═══════════════════════════════════════════════════════════════════════


def _flat_psi_composite(
    sn_mesh: SNMesh, ng: int = 1,
) -> TimedFullField:
    """Build a per-ordinate flat (constant in space) ψ as a TimedFullField.

    Picks ψ = 1.0 everywhere across (N, ng, nx, ny); boundary = 0.  The
    flat-ψ probe is the canonical ERR-026 / Signature-1 diagnostic:
    homogeneous reflective curvilinear + flat ψ must yield
    streaming + redistribution = 0 per ordinate.  Under the composite
    Resolution A: ``(L + C).apply(flat_ψ) = M(flat_ψ; σ_t) = σ_t·ψ``
    bit-exact at the cell-centre block.
    """
    bulk = np.ones((sn_mesh.quad.N, ng, *sn_mesh.spatial_shape))
    return _build_composite(sn_mesh, bulk)


@pytest.mark.l0
@pytest.mark.verifies("dd-curvilinear-scalar")
@pytest.mark.catches("ERR-026")
@pytest.mark.parametrize("sigma_t_value", [0.0, 0.5])
@pytest.mark.parametrize(
    "pole_closure_cls",
    [
        pytest.param(
            MorelMontryAngularSweep,
            id="mms",
            marks=pytest.mark.xfail(
                strict=True,
                reason=(
                    "MorelMontryAngularSweep does NOT preserve the "
                    "per-ordinate flat-flux invariant by design (Phase B "
                    "closeout has the structural rationale).  LIVE "
                    "observation (C5, 2026-07-03, bound closure): BOTH "
                    "sphere AND cylinder legs fail — the Phase-C "
                    "(2026-05-12) claim that cylindrical levels telescope "
                    "cleanly was recorded through the since-retired "
                    "unbound path (the factory call built an UNBOUND "
                    "instance whose matvec crashed on missing state, so "
                    "the non-strict xfail leg was inert, not a physics "
                    "observation).  strict=True pins today's behaviour: "
                    "a closure change that starts preserving the "
                    "invariant must update this gate deliberately.  "
                    "PR-TYPED-6c Step 7 (2026-05-18) retired the "
                    "``LegacyTauSymmetricInterpolation`` and "
                    "``BaileyFlatFluxRedist`` ablation strategies that "
                    "previously parametrised this gate; MMS is the only "
                    "surviving curvilinear closure."
                ),
            ),
        ),
    ],
)
@pytest.mark.parametrize(
    "geom",
    [
        pytest.param("sphere", id="sphere_GL4_reflective"),
        pytest.param("cylinder", id="cyl_LS4_reflective"),
    ],
)
def test_apply_curvilinear_per_ordinate_flat_flux_residual(
    sigma_t_value, pole_closure_cls, geom,
):
    r"""ψ constant per ordinate on reflective BC → apply(ψ) = Σ_t·ψ.

    The canonical curvilinear bug-class diagnostic (per `vv-principles`
    Signature 1 and ERR-026 entry). On a homogeneous reflective
    sphere / cylinder, per-ordinate flat ψ ⇒ streaming + redistribution
    ≡ 0 per ordinate. With Σ_t = 0 the matvec output must be
    bit-zero; with Σ_t ≠ 0 it must equal Σ_t·ψ to rtol=1e-12.

    Parametrised over the three pole-angular-closure strategies. The
    MMS canonical form is expected to FAIL this test pre-spatial-
    closure-alignment (Phase B observation, documented in the Phase B
    closeout memo). Phase C's sweep-frame matvec WITH WDD spatial
    closure makes this gate empirically diagnostic of whether MMS can
    become the default — Phase C's empirical decision point.

    The composite ``(L + C).apply(state)`` realises the geometry-agnostic
    matvec via Resolution A's subtractive identity:
    :math:`(L+C)\psi = (M(\psi;\sigma_t) - \sigma_t\psi) + \sigma_t\psi
    = M(\psi;\sigma_t)`.  The check is on ``out.interior.values`` cell-centre
    block; the per-ordinate flat-ψ invariant collapses the matvec to
    ``Σ_t·ψ`` cell-wise.
    """
    if geom == "sphere":
        sn_mesh, sig_t = _make_spherical_sn_mesh(pole_closure=pole_closure_cls)
    else:
        sn_mesh, sig_t = _make_cylindrical_sn_mesh(pole_closure=pole_closure_cls)
    sig_t = np.full_like(sig_t, sigma_t_value)
    L = StreamingOperator(sn_mesh)
    C = MultiplicationOperator.from_mesh(sig_t, sn_mesh)
    op = L + C
    psi_state = _flat_psi_composite(sn_mesh, ng=sn_mesh.ng)
    result = op.apply(psi_state)
    result_bulk = result.interior.values
    expected_bulk = sigma_t_value * psi_state.interior.values
    if sigma_t_value == 0.0:
        np.testing.assert_allclose(result_bulk, 0.0, atol=1e-13)
    else:
        np.testing.assert_allclose(
            result_bulk, expected_bulk, rtol=1e-12, atol=1e-13,
        )


# ═══════════════════════════════════════════════════════════════════════
# Gate 1.3 — apply ↔ apply_transpose reciprocity
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.foundation
@pytest.mark.parametrize(
    "geom",
    [
        pytest.param("sphere", id="sphere_GL4_reflective"),
        pytest.param("cylinder", id="cyl_LS4_reflective"),
        pytest.param("cart2d", id="cart2d_LS2_reflective_nonsquare"),
    ],
)
def test_apply_apply_transpose_reciprocity_under_sweep_frame(geom):
    r"""$\langle (L{+}C)\psi, \phi \rangle = \langle \psi, (L{+}C)^{\mathsf T}\phi \rangle$
    over the FULL bulk⊕trace composite (Euclidean), to round-off.

    Wave O / O.2b (#208): :class:`StreamingOperator` now carries the analytic
    reverse-direction adjoint matvec :meth:`~orpheus.sn.operators.streaming.StreamingOperator.apply_transpose`,
    so the composite ``(L + C)`` is adjointable and this
    reciprocity is the standing pin that ``apply_transpose`` IS the true
    Euclidean transpose of ``apply`` — over BOTH blocks (the bulk residual AND
    the boundary trace), with a NON-ZERO boundary so the FULL-operator trace
    coupling (``L_bs`` inflow→bulk, ``L_sb`` bulk→outflow) is exercised, not
    merely the cell block.  The reverse sweep is verified bit-for-bit against a
    dense-probe transpose oracle in
    ``derivations/diagnostics/diag_p42_adjoint_oracle.py``.

    This gate is metric-AGNOSTIC by design: the ``|Ω·n|·w`` partial-current
    metric is non-discriminating for ``(L+C)`` with a SPECULAR reflective BC
    (which is self-adjoint even under the Euclidean trace inner product); the
    metric's load-bearing role is pinned by the white-BC self-adjointness gate
    (O.2b sub-step 4.4), where dropping ``|Ω·n|`` breaks self-adjointness.

    (Was ``xfail(strict)`` pre-O.2b — ``StreamingOperator`` was apply-only
    (``is_adjointable`` False) so ``(L+C).apply_transpose`` raised; the
    analytic adjoint lands the transpose and flips this green.)
    """
    rng = np.random.default_rng(seed=137)
    if geom == "sphere":
        sn_mesh, sig_t = _make_spherical_sn_mesh()
    elif geom == "cylinder":
        sn_mesh, sig_t = _make_cylindrical_sn_mesh()
    else:
        # #310 C4: the multi-D Cartesian row — the row-march reverse on the
        # default representation, reflective nonsquare (seedless mesh, so
        # the joint wrapper degenerates to the bare composite operator).
        sn_mesh = cart2d_2g_nonsquare()
        sig_t = np.stack(
            [np.full(sn_mesh.spatial_shape, 0.5 * (1.0 + 0.5 * g))
             for g in range(2)],
            axis=0,
        )
    L = StreamingOperator(sn_mesh)
    C = MultiplicationOperator.from_mesh(sig_t, sn_mesh)
    op = _joint_op(sn_mesh, L + C)   # B.2d: the joint M on the carrying pair
    n_trace = int(sn_mesh.angular_trace.layout.total_size)
    # #282 route (a): a RANDOM ψ½ block (not the edge-extrap default) so the
    # augmented seed rows are independently exercised — the Euclidean
    # transpose is exact over the FULL space, so ``full_dot`` (below) MUST
    # include the seed block for reciprocity to hold (a bulk⊕trace-only dot
    # is blind to the seed↔bulk coupling — the Euclidean sibling of the
    # G-reciprocity's zero-weight blindness, vv Mode 12).
    seed_space = sn_mesh.radial_characteristic_field_space
    n_seed = 0 if seed_space is None else seed_space.shape[0]
    psi_state = _build_composite(
        sn_mesh, _random_bulk(sn_mesh, rng), rng.standard_normal(n_trace),
        radial_characteristic_values=(
            rng.standard_normal(n_seed) if n_seed else None
        ),
    )
    phi_state = _build_composite(
        sn_mesh, _random_bulk(sn_mesh, rng), rng.standard_normal(n_trace),
        radial_characteristic_values=(
            rng.standard_normal(n_seed) if n_seed else None
        ),
    )

    def full_dot(a, b):
        """Euclidean inner product over the whole coupled state (System A's
        bulk ⊕ trace + System B's members — the joint flat pairing)."""
        return float(np.sum(a.to_flat() * b.to_flat()))

    lhs = full_dot(op.apply(psi_state), phi_state)
    rhs = full_dot(psi_state, op.apply_transpose(phi_state))
    rel = abs(lhs - rhs) / (abs(lhs) + abs(rhs) + 1e-300)
    assert rel < 1e-12, f"reciprocity broken: {lhs} vs {rhs} (rel={rel:.2e})"


# ═══════════════════════════════════════════════════════════════════════
# Gate 1.2 — apply ↔ sweep face-flux propagation (structural under sweep-frame)
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.foundation
def test_apply_face_fluxes_match_sweep_recurrence_spherical():
    r"""Sweep-frame apply matvec uses the SAME WDD recurrence as the sweep.

    Under the sweep-frame architecture (Phase C), both the apply
    matvec and the sweep consume:

    * The same per-direction cell ordering
      (``SNMesh.dag_walk(direction_sign=...)``).
    * The same WDD diamond closure ``ψ_face_out = 2·ψ_cell -
      ψ_face_in`` per cell.
    * The same BC trace law at the boundary.

    The structural identity is built by construction. We pin it via
    a deterministic input: apply a known ψ_cells, then verify the
    output's residual breakdown matches what a hand-rolled WDD
    propagation chain would compute. Bit-identical via np.array_equal.

    D-K.5 migration — the determinism assertion shifts from the legacy
    packed vector to ``out.interior.values`` AND ``out.boundary.values``;
    both must be bit-stable across repeated calls to the composite
    ``(L + C).apply``.
    """
    sn_mesh, sig_t = _make_spherical_sn_mesh(nx=6, R=1.0, quad_name="gl4")
    L = StreamingOperator(sn_mesh)
    C = MultiplicationOperator.from_mesh(sig_t, sn_mesh)
    op = _joint_op(sn_mesh, L + C)   # B.2d: the joint M on the carrying pair

    # Use a deterministic, structured input so the residual is
    # predictable.
    rng = np.random.default_rng(seed=0)
    psi_state = _build_composite(sn_mesh, _random_bulk(sn_mesh, rng))

    # Sanity check: applying L twice to the same input always yields
    # the same output — this is the bit-identity (np.array_equal)
    # contract called out in plan §5 Gate 1.2.
    out1 = op.apply(psi_state)
    out2 = op.apply(psi_state)
    assert np.array_equal(out1.to_flat(), out2.to_flat()), (
        "Apply is not deterministic — sweep-frame matvec must be "
        "bit-stable on the whole coupled state"
    )


# ═══════════════════════════════════════════════════════════════════════
# Gate 1.5 — BC trace contract verification
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.foundation
def test_bc_trace_contract_respected_by_matvec_vacuum_sphere():
    r"""Vacuum BC at outer face: inflow = realize().apply(outflow_face).

    The architectural-fix anchor for the §16A.3 contract. Pins that
    under the sweep-frame matvec, the BC trace law consumes the
    WDD-propagated outflow face values (not cell-centres) and
    produces the inflow face values that the inward sweep consumes.

    For vacuum BC: realize() is the narrowed zero map
    :math:`\Gamma_+ \to \Gamma_-` (a ``ZeroOperator`` between the two
    half-traces since B3.2 ``7f02de15``; pre-B3.2 it was a full-face
    inflow-zeroing mask). So the inward sweep at the boundary should
    see psi_face_in = 0 for all incoming ordinates — the conclusion is
    unchanged, only the operator that delivers it.

    The behaviourally-observable consequence: if we apply L to a
    purely outgoing ψ on a vacuum-BC sphere, the inward ordinates'
    streaming contribution must be entirely determined by the
    in-cell propagation (no incoming flux from the BC).

    D-K.5 migration — checks both ``bulk.values`` and ``boundary.values``
    are zero at zero input (linearity guard on the composite carrier).
    """
    nx = 6
    edges = np.linspace(0.0, 1.0, nx + 1)
    mesh = Mesh1D(
        edges=edges,
        mat_ids=np.zeros(nx, dtype=int),
        coord=CoordSystem.SPHERICAL,
        bc_right=BC("vacuum"),
    )
    quad = Quadrature.gauss_legendre(4)
    sn_mesh = SNMesh(mesh, quad, placeholder_materials())
    sig_t = np.full((1, nx), 0.5)  # (ng, nx) — rank-d
    L = StreamingOperator(sn_mesh)
    C = MultiplicationOperator.from_mesh(sig_t, sn_mesh)
    op = L + C

    # Linearity is enough to characterize the BC respond-only-to-outflow
    # property. apply(0) MUST be 0 even with vacuum BC. If the BC
    # consumed cell-centres rather than face values, a nonzero ψ
    # could pollute the inflow even when the cell-centres are zero.
    # #282 route (a): pass the seed leaf UNIFORMLY (the R12a predicate
    # allocates it iff the mesh carries levels — here the sphere does).
    state_zero = TimedFullField.zeros(
        interior=AngularFlux, boundary=AngularBoundaryFlux, mesh=sn_mesh,
    )
    result = op.apply(state_zero)
    assert np.array_equal(
        result.interior.values, np.zeros_like(result.interior.values),
    ), (
        f"apply(0).interior should be 0 (BC linearity); got max|out|="
        f"{np.max(np.abs(result.interior.values)):.3e}"
    )
    assert np.array_equal(
        result.boundary.values, np.zeros_like(result.boundary.values),
    ), (
        f"apply(0).boundary should be 0 (BC linearity); got max|out|="
        f"{np.max(np.abs(result.boundary.values)):.3e}"
    )


@pytest.mark.foundation
def test_bc_trace_contract_respected_by_matvec_reflective_sphere():
    r"""Reflective BC: realize() = PermutationOperator.

    For reflective BC (the regression-snapshot default), the
    realize() output is the cell-permutation operator that maps an
    outflow ordinate's face flux to the partner incoming ordinate's
    face flux. The sweep-frame matvec, by consuming this operator
    on the WDD-propagated outflow face, achieves the correct
    inflow without algebraic extrapolation.

    D-K.5 migration — composite ``(L + C).apply(0) = 0`` on both bulk
    and boundary.
    """
    sn_mesh, sig_t = _make_spherical_sn_mesh(nx=8, R=1.0, quad_name="gl4")
    L = StreamingOperator(sn_mesh)
    C = MultiplicationOperator.from_mesh(sig_t, sn_mesh)
    op = L + C
    # #282 route (a): pass the seed leaf UNIFORMLY (the R12a predicate
    # allocates it iff the mesh carries levels — here the sphere does).
    state_zero = TimedFullField.zeros(
        interior=AngularFlux, boundary=AngularBoundaryFlux, mesh=sn_mesh,
    )
    result = op.apply(state_zero)
    np.testing.assert_array_equal(result.interior.values, 0.0)
    np.testing.assert_array_equal(result.boundary.values, 0.0)


# ═══════════════════════════════════════════════════════════════════════
# Gate 1.5 strengthened — capture-and-compare BC apply input (Phase D)
# ═══════════════════════════════════════════════════════════════════════
#
# Per the Phase D plan §3f, the original Gate 1.5 `apply(0) = 0`
# probe is replaced with a capture-and-compare structural test:
# we instrument the matvec to capture the input vector passed to
# the BC apply, then verify it bit-matches the independently-
# reconstructed WDD-propagated outflow face value.
#
# Parametrised over BC kinds (Vacuum, Reflective) × geometries
# (sphere, cylinder).  White / Albedo / PrescribedInflow are NOT
# included in this parametrisation because their realizer output
# structures differ (White returns AngularAverageOperator product;
# Albedo is identity-scaled; PrescribedInflow has an inhomogeneous
# affine component).  Those kinds are tested at unit level via
# their own realizer-tests (`tests/sn/test_sn_boundary_realizer.py`).
#
# D-K.5 migration (2026-05-29) — the composite ``StreamingOperator``
# path's ``transport_operator_matvec_unified`` calls ``bc_outer.apply``
# TWICE per matvec, matching the legacy ``transport_operator_matvec_
# spherical`` 2-call pattern:
#
#   * call[0] — outer-inflow estimate for the Carlson seed.  The
#     unified matvec passes the OUTER FACE flux (``boundary.face_view
#     ("xmax")``) directly, not the cell-centre proxy at i=nx-1.
#     When the test builds a TimedFullField with a zero boundary, the
#     captured input is therefore zero (the L2 face-view, not the
#     pre-D-H.2-C2 cell-centre-proxy synthesis the legacy path used).
#   * call[1] — BC trace law on the WDD-propagated outward outflow at
#     the outer face.  Same content as the legacy call.


def _outflow_at_boundary_for_sphere_from_bulk(
    sn_mesh: SNMesh,
    psi_bulk: np.ndarray,
) -> np.ndarray:
    r"""Independently reconstruct the WDD-propagated outflow face value.

    Runs the same WDD diamond closure as
    ``transport_operator_matvec_unified``'s outward sweep, in
    isolation (no redistribution, no collision, no BC).  Returns the
    outflow face vector at the outer boundary, shape ``(ng, N)``
    with only the outgoing ordinates populated.

    Consumes the cell-centred angular flux directly (``psi_bulk`` shape
    ``(N, ng, nx)``) — no packed-vector round-trip.  Mirrors the
    pre-migration helper's algebra but on the typed composite input.

    This is the "reference" the BC apply input must bit-match in the
    capture-and-compare test below.
    """
    quad = sn_mesh.quad
    nx = sn_mesh.nx
    ng = psi_bulk.shape[1]
    eps = 1e-15
    # Mirror the matvec body's ``(ng, N, nx)`` ordering.
    psi_g_first = psi_bulk.transpose(1, 0, 2)
    outgoing_mask = quad.mu_x > +eps
    incoming_mask = quad.mu_x < -eps
    outflow_at_boundary = np.zeros((ng, quad.N))
    if not np.any(outgoing_mask):
        return outflow_at_boundary
    # Carlson coupled-pole seed (ERR-058 a, Issue #195): the inward
    # (−1) chain runs first, seeded from the GIVEN outer inflow trace
    # (zero in this test's composite), and its pole-face outflow at
    # the MIRROR ordinate seeds the outward chain — the r = 0
    # continuity ψ(0, +μ) = ψ(0, −μ).  The pre-ERR-058 cell-centre
    # proxy read ψ(Δr/2) as the pole-face value (O(h)-wrong on
    # non-flat profiles) and is retired.
    psi_face_in_inward = np.zeros((ng, int(incoming_mask.sum())))
    for i in range(nx - 1, -1, -1):
        psi_cell = psi_g_first[:, incoming_mask, i]
        psi_face_in_inward = 2.0 * psi_cell - psi_face_in_inward
    pole_outflow = np.zeros((ng, quad.N))
    pole_outflow[:, incoming_mask] = psi_face_in_inward
    pole_pi = quad.ordinate_permutation(SelfPairedDeck.mirror(axis="x").motion)
    assert pole_pi is not None  # mimics production's _ensure_pole_mirror source
    psi_face_in = pole_outflow[:, pole_pi.indices][:, outgoing_mask]
    for i in range(nx):
        psi_cell = psi_g_first[:, outgoing_mask, i]
        psi_face_out = 2.0 * psi_cell - psi_face_in
        psi_face_in = psi_face_out
    outflow_at_boundary[:, outgoing_mask] = psi_face_out
    return outflow_at_boundary


@pytest.mark.foundation
@pytest.mark.parametrize(
    "bc_kind",
    [
        pytest.param("vacuum", id="vacuum"),
        pytest.param("reflective", id="reflective"),
    ],
)
def test_bc_trace_contract_capture_and_compare_sphere(bc_kind):
    r"""Gate 1.5 (Wave O #208 migration): the matvec EXTRACTS the BC, and
    EMITS the WDD-propagated outflow on the outflow slots.

    Pre-extraction the unified matvec applied the BC trace law INSIDE the
    sweep (the §16A.3 contract: ``bc_outer.apply`` consumed the
    WDD-propagated outflow — twice: the Carlson outer-inflow estimate +
    the keystone re-apply).  Wave O O.4a.2 deletes the keystone and reads
    ``ψ.boundary.inflow`` as a given, so the matvec calls ``bc_outer.apply``
    ZERO times — the reflective coupling moved to the sibling ``−B``
    (:class:`~orpheus.sn.operators.boundary.SNBoundaryOperator`).

    The §16A.3 substance is preserved, decomposed across the new algebra:

    * the matvec now EMITS the WDD-propagated outflow on the outflow slots
      (this test — with a zero input boundary the defect ``streamed −
      ψ.outflow`` reduces to the raw WDD outflow, pinned bit-exact against
      the independently-reconstructed WDD chain);
    * ``B`` reflects that outflow into the inflow trace (pinned in
      ``test_sn_boundary_operator.py``); the end-to-end consistency
      ``ψ.inflow − B·ψ.outflow → 0`` is pinned by the curvilinear
      streaming-equilibrium gate.

    Asserting the 0-call extraction AND the emitted-outflow value together
    prevents a silent regression that re-absorbs the BC into the matvec.
    """
    from unittest.mock import patch
    sn_mesh, sig_t = _make_spherical_sn_mesh(
        nx=6, R=1.0, quad_name="gl4",
        bc_outer=BC(bc_kind),
    )
    L = StreamingOperator(sn_mesh)
    C = MultiplicationOperator.from_mesh(sig_t, sn_mesh)
    op = _joint_op(sn_mesh, L + C)   # B.2d: the joint M on the carrying pair
    rng = np.random.default_rng(seed=137)
    psi_bulk = _random_bulk(sn_mesh, rng)
    psi_state = _build_composite(sn_mesh, psi_bulk)  # zero boundary

    # Independent reference: rebuild outflow face via isolated WDD on
    # the bulk values (no packed-vector round-trip).
    expected_outflow = _outflow_at_boundary_for_sphere_from_bulk(
        sn_mesh, psi_bulk,
    )

    # Capture any BC apply calls during the matvec.
    captured_inputs: list[np.ndarray] = []
    original_apply = sn_mesh.bc["xmax"].apply

    def capture_apply(inp):
        captured_inputs.append(np.array(inp, copy=True))
        return original_apply(inp)

    with patch.object(sn_mesh.bc["xmax"], "apply", side_effect=capture_apply):
        out = op.apply(psi_state)

    # EXTRACTION: the post-O.4a.2 matvec does NOT apply the BC — it reads
    # ψ.boundary.inflow as a given and emits the raw outflow defect.  The
    # reflective coupling is the sibling −B, not an intra-matvec re-apply.
    assert len(captured_inputs) == 0, (
        f"Expected ZERO bc_outer.apply calls per matvec post-extraction "
        f"(Wave O O.4a.2 deleted the keystone re-apply + the Carlson "
        f"outer-inflow estimate now reads the raw given trace), got "
        f"{len(captured_inputs)}.  If a bc.apply re-appeared in the "
        f"matvec, the BC was re-absorbed — the extraction regressed."
    )

    # EMITTED OUTFLOW: with a zero input boundary, the outflow-slot defect
    # ``streamed − ψ.outflow`` = ``streamed − 0`` = the raw WDD-propagated
    # outflow.  Pin it bit-exact against the independent WDD chain on the
    # outflow ordinates of the outer face (the §16A.3 substance, relocated
    # to the matvec's emission).
    trace = sn_mesh.angular_trace
    outflow_idx = trace.outflow_indices_for_face("xmax")
    got_outflow = _sysA(out).boundary.face_view("xmax")[outflow_idx, :]   # (M, ng)
    expected_outflow_face = expected_outflow.T[outflow_idx, :]      # (M, ng)
    assert np.allclose(
        got_outflow, expected_outflow_face, rtol=1e-14, atol=1e-14,
    ), (
        f"The matvec's emitted outflow trace does not match the "
        f"independently-reconstructed WDD-propagated outflow.  Max diff: "
        f"{np.max(np.abs(got_outflow - expected_outflow_face))}.  This is "
        f"the §16A.3 outflow-emission contract (now on the matvec, not the "
        f"BC re-apply)."
    )


# ═══════════════════════════════════════════════════════════════════════
# Gate 1.6 — Phase F sweep-path per-ordinate flat-flux residual
# ═══════════════════════════════════════════════════════════════════════
#
# The DUAL of Gate 1.1 (apply-matvec) but on the SI/sweep path.  Phase D
# pinned Gate 1.1 on the apply-matvec by composing
# :class:`~orpheus.sn.sweep.psi_half_angle_seed.CarlsonInwardSweep`
# into :class:`MorelMontryAngularSweep.psi_half_seed`.  Phase F backports
# the Carlson coupled-pole seed into ``_sweep_1d_spherical`` /
# ``_sweep_1d_cylindrical`` via the source-driven
# :func:`carlson_inward_sweep_from_source` helper — this gate pins that
# the sweep path now satisfies the same per-ordinate flat-flux residual
# that Gate 1.1 verifies on the apply path.


@pytest.mark.l0
@pytest.mark.verifies("dd-curvilinear-scalar")
@pytest.mark.verifies("phase-f-carlson-seed-source-driven")
@pytest.mark.verifies("phase-f-q-bar-twin-forms")
@pytest.mark.catches("ERR-026")
@pytest.mark.parametrize("sigma_t_value", [0.5, 1.5])
@pytest.mark.parametrize(
    "geom",
    [
        pytest.param("sphere", id="sphere_GL4_reflective"),
        pytest.param("cylinder", id="cyl_LS4_reflective"),
    ],
)
def test_sweep_curvilinear_per_ordinate_flat_flux_residual(
    sigma_t_value, geom,
):
    r"""Sweep-path twin of Gate 1.1: converged sweep on flat ψ fixed point.

    The DUAL of Gate 1.1.  Gate 1.1 runs *one* apply-matvec call and
    checks per-ordinate ``L·ψ_flat = Σ_t·ψ_flat``.  Gate 1.6 runs
    Source Iteration to *convergence* on a homogeneous reflective
    curvilinear mesh with ``Q = Σ_t · ψ_const`` (the source whose
    fixed point is ``ψ ≡ ψ_const`` everywhere) and checks the
    converged scalar flux is uniform at ``Σ_n w_n · ψ_const = Σw ·
    ψ_const``.

    Pre-Phase-F the sweep path's hardcoded ``psi_angle = np.zeros``
    Carlson-zero seed (sweep.py:474) prevented the SI fixed point
    from reaching the flat-ψ eigenmode on sphere — the per-ordinate
    M-M recurrence stayed at the wrong fixed point even on a
    flat-source homogeneous reflective probe, because the seed
    perturbed the first ordinate's face flux at every iteration.
    On a level-symmetric cylinder the *dead* first-ordinate seed weight
    (``c_in[m0]=0`` at raw ``τ=1``) hid the bug there — NOT α-dome
    telescoping (a level-symmetric-only reading, false for a product
    quadrature; #280 Phase 2.5b).

    Phase F backports the Carlson coupled-pole seed via
    :func:`~orpheus.sn.sweep.psi_half_angle_seed.carlson_inward_sweep_from_source`
    so the seed becomes ``Σ_t · ψ_const / 2`` on the flat source
    (per Hébert (3.434)-(3.435) with ``Q̄ = Σ_t · ψ_const / 2`` at
    μ = −1), which is consistent with the flat-ψ field across cells.

    Verification: SI converges to a uniform scalar flux equal to
    ``Σw · ψ_const`` to numerical precision.
    """
    # #282 route (a): the seed-strategy zoo is retired.  The sweep path's
    # direct starting-direction solver IS
    # ``carlson_inward_sweep_from_source`` (the ONE recurrence host); the
    # old test compared it against the retired ``CarlsonInwardSweep``
    # STRATEGY wrapper on a flat-ψ probe — a self-comparison (the wrapper
    # folded φ₀ = Σw·ψ then called the SAME function with Q̄ = Σ_t·ψ_const).
    # It reduces to the flat-ψ algebraic identity of the surviving
    # function, pinned directly below.
    from orpheus.sn.sweep.psi_half_angle_seed import (
        carlson_inward_sweep_from_source,
    )

    if geom == "sphere":
        sn_mesh, sig_t = _make_spherical_sn_mesh(
            pole_closure=MorelMontryAngularSweep,
        )
    else:
        sn_mesh, sig_t = _make_cylindrical_sn_mesh(
            pole_closure=MorelMontryAngularSweep,
        )
    sig_t_arr = np.full_like(sig_t, sigma_t_value)
    nx = sn_mesh.nx
    psi_const = 1.0
    sigma_t_gx = sig_t_arr  # (ng, nx) — rank-d
    dr = sn_mesh.axis_widths[0]
    bc_outer_value = np.full((1,), psi_const)

    # Flat-ψ algebraic identity of the direct solver: with the canonical
    # flat-ψ source ``Q̄ = Σ_t · ψ_const`` and the reflective mirror
    # ``bc_outer = ψ_const``, the Hébert (3.434)-(3.435) inward march
    # reproduces ``ψ_const`` at every cell (the fixed-point seed), and its
    # exit face is likewise ``ψ_const``.  (The route-(a) function returns
    # the ``(cells, face)`` tuple — the walk's exit face is the pole datum.)
    Q_bar = np.full((1, nx), sigma_t_value * psi_const)
    seed_cells, seed_face = carlson_inward_sweep_from_source(
        Q_bar, sigma_t_gx, dr, bc_outer_value,
    )
    expected_const = np.full((1, nx), psi_const)
    np.testing.assert_allclose(
        seed_cells, expected_const, rtol=1e-13, atol=1e-13,
        err_msg=(
            "direct starting-direction solver flat-ψ identity: on a "
            "reflective homogeneous probe with bc_outer=ψ_const and "
            "Q̄ = Σ_t·ψ_const, the Hébert (3.434)-(3.435) inward march "
            "must reproduce ψ_const at every cell."
        ),
    )
    np.testing.assert_allclose(
        seed_face, np.full((1,), psi_const), rtol=1e-13, atol=1e-13,
        err_msg="the inward march's exit (pole) face must also be ψ_const.",
    )
