r"""C5 — the #257 S8a apply(FullField) timeless round-trip gate.

#257 S8a re-points every SN operator matvec LEAF (``L`` streaming, ``C``
collision, ``S`` scattering, ``F`` fission, ``B`` boundary) from emitting the
history-bearing
:class:`~orpheus.transport.timed_full_field.TimedFullField` to the **timeless**
:class:`~orpheus.transport.full_field.FullField`.  This is the cofree-comonad
finding (``TimedFullField = Cofree(FullField, d)``): an operator output is a
base arrow ``FullField -> FullField`` and must NOT carry a history tail — only
the iteration DRIVER (``SourceIteration`` / ``KrylovAcceleration``) sees the
comonad.

The gate has three legs (gate spec
``.claude/agent-memory/test-architect/issue_257_s8_streaming_pure_L_verification.md``
catcher C5):

* **C5a — timeless codomain (the defining property).**  ``L.apply(x)`` /
  ``C.apply(x)`` / ``F.apply(x)`` / ``B.apply(x)`` / ``S.apply(x)`` return a
  TIMELESS ``FullField`` (``type(out) is FullField``, NOT the timed subclass),
  history-free, regardless of the input iterate's ``history_depth``.  Per
  geometry (slab / sphere / cylinder).  **Mode-11**: the matvec leaf has ZERO
  graph callers (it is reached only via ``OperatorSum`` / the driver), so C5a
  calls ``op.apply(x)`` DIRECTLY to reach the re-typed leaf body — a solve-only
  test routes through the sweep (``loss_representation``) and never touches the
  matvec emit path.
* **C5b — driver re-attach, byte-identical converged state.**  The SI / Krylov
  driver adds the timeless gain outputs to the TIMED rhs (the ``TimedFullField``
  ``__add__`` recombine re-attaches the timed type), so the converged solution
  is unchanged.  Verified against the STRUCTURALLY-INDEPENDENT closed-form
  ``k_inf = νΣ_f / Σ_a`` (NOT an old-vs-new ULP proximity, NOT a hardcoded
  baseline that L28 would force a later mutation of), per geometry.
* **C5c — the ``advance`` history type-guard still fires.**  The driver-level
  comonad verb (``timed_full_field.py``) still rejects a mismatched-type new
  frame.

``-O``-firing (Mode 8): every check is a ``pytest.fail`` / ``np.testing.*`` /
``pytest.raises`` — NO bare ``assert`` on the load-bearing predicate (the
canonical ORPHEUS invocation is ``-O``, which strips bare ``assert``).
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.derivations.continuous.analytical.homogeneous import (
    derive_1g_continuous,
    derive_2g_continuous,
)
from orpheus.geometry import (
    BC,
    Mesh1D,
    Region,
    RegionMesh,
    StructuredGeometry,
)
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.operators.boundary import SNBoundaryOperator
from orpheus.sn.mesh.augmented_mesh import SNMesh
from orpheus.sn.operators.streaming import StreamingOperator
from orpheus.transport.operators.multiplication_operator import MultiplicationOperator
from orpheus.sn.solver import SNSolver, solve_sn
from orpheus.transport.fields.angular_flux import AngularFlux
from orpheus.transport.fields.angular_boundary_flux import AngularBoundaryFlux
from orpheus.transport.full_field import FullField
from orpheus.transport.timed_full_field import TimedFullField

pytestmark = pytest.mark.foundation


# ─── geometry fixtures (homogeneous fissile, reflective → k_eff = k_inf) ──

_GEOMETRY_TAG = {"slab": "SLB", "sphere": "SPH", "cylinder": "CYL"}
_CASE_BUILDERS = {"1eg": derive_1g_continuous, "2eg": derive_2g_continuous}


def _get_continuous_case(ng_key: str):
    """Name the analytical homogeneous derivation explicitly (V&V hygiene)."""
    return _CASE_BUILDERS[ng_key]()


def _homogeneous_mesh(coord: str, mat_id: int) -> Mesh1D:
    refl = BC("reflective")
    bcs = (refl, refl) if coord == "slab" else (refl,)
    geom = StructuredGeometry(
        geometry=_GEOMETRY_TAG[coord],
        regions=(Region(mat_id=mat_id, outer_thickness_cm=2.0),),
        bcs=bcs,
    )
    return Mesh1D.from_geometry(geom, region_meshes=(RegionMesh(n_cells=6),))


def _quadrature_for(coord: str):
    if coord == "cylinder":
        return Quadrature.folded_product(n_mu=2, n_phi=4)
    return Quadrature.gauss_legendre(n_ordinates=8)


def _solver_for(coord: str, ng_key: str) -> tuple[SNSolver, object]:
    """Build an :class:`SNSolver` on a homogeneous fissile medium + its case."""
    from orpheus.sn.solver import _as_sn_mesh

    case = _get_continuous_case(ng_key)
    mat_id = next(iter(case.problem.materials.keys()))
    mesh = _homogeneous_mesh(coord, mat_id)
    quad = _quadrature_for(coord)
    sn_mesh = _as_sn_mesh(mesh, quad, case.problem.materials)
    solver = SNSolver(sn_mesh)
    return solver, case


def _timed_random_state(sn_mesh: SNMesh, *, history_depth: int, seed: int) -> TimedFullField:
    """A timed iterate with random bulk — the comonad the driver carries."""
    state = TimedFullField.zeros(
        interior=AngularFlux, boundary=AngularBoundaryFlux, space=sn_mesh.full_field_space,
        history_depth=history_depth,
    )
    from dataclasses import replace

    rng = np.random.default_rng(seed)
    state = replace(
        state,
        interior=replace(state.interior, values=rng.standard_normal(state.interior.values.shape)),
    )
    return state


# ═══════════════════════════════════════════════════════════════════════
# C5a — timeless codomain: every matvec leaf is a base arrow.
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.parametrize("coord", ["slab", "sphere", "cylinder"])
def test_c5a_matvec_leaves_emit_timeless_full_field(coord: str) -> None:
    """L / C / S / F / B ``.apply`` outputs are a TIMELESS ``FullField``.

    Mode-11: calls each leaf's ``.apply`` DIRECTLY (the matvec leaf has zero
    graph callers — a solve-only test routes around it through the sweep).
    The input is a history-bearing ``TimedFullField`` iterate (depth=2); the
    output MUST be history-free (``type(out) is FullField``), proving the leaf
    drops the comonad tail.
    """
    solver, _case = _solver_for(coord, "2eg")
    sn_mesh = solver.sn_mesh

    L = StreamingOperator.pose(sn_mesh)
    C = MultiplicationOperator.from_mesh(solver.mat_xs.total_cross_section_field, sn_mesh)
    S = solver.scattering_op
    F = solver.fission_op
    B = SNBoundaryOperator(sn_mesh)

    state = _timed_random_state(sn_mesh, history_depth=2, seed=11)

    leaves = {"L": L.apply, "C": C.apply, "S": S.apply, "F": F.apply, "B": B.apply}
    for name, apply in leaves.items():
        out = apply(state)
        # base-arrow codomain: a timeless FullField, NOT the timed subclass.
        if type(out) is not FullField or isinstance(out, TimedFullField):
            pytest.fail(
                f"[{coord}] {name}.apply must emit a TIMELESS FullField "
                f"(base arrow, history-free); got {type(out).__name__}"
            )
        # No history tail leaked onto the output.
        if hasattr(out, "history_length") or hasattr(out, "history_depth"):
            pytest.fail(
                f"[{coord}] {name}.apply output carries a history attribute — "
                f"the comonad must live on the driver, not the operator leaf."
            )


@pytest.mark.parametrize("coord", ["slab", "sphere", "cylinder"])
def test_c5a_apply_transpose_emits_timeless_full_field(coord: str) -> None:
    """L / C / B ``.apply_transpose`` outputs are a TIMELESS ``FullField``."""
    solver, _case = _solver_for(coord, "2eg")
    sn_mesh = solver.sn_mesh

    L = StreamingOperator.pose(sn_mesh)
    C = MultiplicationOperator.from_mesh(solver.mat_xs.total_cross_section_field, sn_mesh)
    B = SNBoundaryOperator(sn_mesh)

    state = _timed_random_state(sn_mesh, history_depth=3, seed=23)

    for name, op in {"L": L, "C": C, "B": B}.items():
        out = op.apply_transpose(state)
        if type(out) is not FullField or isinstance(out, TimedFullField):
            pytest.fail(
                f"[{coord}] {name}.apply_transpose must emit a TIMELESS "
                f"FullField; got {type(out).__name__}"
            )


def test_c5a_independent_of_input_history_depth() -> None:
    """The timeless codomain holds for EVERY input ``history_depth`` (0..4)."""
    solver, _case = _solver_for("slab", "2eg")
    sn_mesh = solver.sn_mesh
    L = StreamingOperator.pose(sn_mesh)
    for depth in (0, 1, 2, 4):
        state = TimedFullField.zeros(
            interior=AngularFlux, boundary=AngularBoundaryFlux, space=sn_mesh.full_field_space,
            history_depth=depth,
        )
        out = L.apply(state)
        if type(out) is not FullField or isinstance(out, TimedFullField):
            pytest.fail(
                f"depth={depth}: L.apply output must be a timeless FullField, "
                f"got {type(out).__name__}"
            )


# ═══════════════════════════════════════════════════════════════════════
# C5b — driver re-attach: byte-identical converged state vs closed-form.
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.parametrize("inner_solver", ["source_iteration", "krylov"])
@pytest.mark.parametrize("coord", ["slab", "sphere", "cylinder"])
def test_c5b_driver_reattach_recovers_kinf(coord: str, inner_solver: str) -> None:
    """The SI / Krylov driver re-attaches the timeless gain outputs into the
    timed iterate, so the converged eigenvalue is unchanged.

    The matvec leaves now emit a timeless ``FullField`` source; the driver
    builds ``rhs = q_ext + sum_i g_i.apply(psi)`` where ``q_ext`` is the timed
    composite, so the ``TimedFullField.__add__`` recombine re-attaches the
    timed type and the iterate stays a ``TimedFullField`` (the comonad lives on
    the driver).  Verified against the STRUCTURALLY-INDEPENDENT closed-form
    ``k_inf = nu*Sigma_f / Sigma_a`` (the eigenvalue ground — MMS does not
    prove eigenvalues), NOT an old-vs-new ULP proximity.

    History — the retired ``sphere-krylov`` exclusion.  Until 2026-08-10 this
    row imperatively xfailed ``sphere × krylov`` as "ill-conditioned at high
    scattering ratio (issue #200 — same exclusion as test_kinf_homogeneous)".
    It was NOT the same exclusion: that module excluded only
    ``sphere-4eg-krylov``, while this row is fixed at **2eg** and its condition
    named no group count at all, so it excluded a combination its own cited
    authority never did — ``test_kinf_homogeneous[sphere-2eg-krylov]`` has
    passed throughout.  Retired as healed: measured with the imperative xfail
    neutralised, the row passes with ``rel = 1.184e-14`` against the
    closed-form reference (gate: ``rtol=1e-10``) and emits no
    :class:`~orpheus.numerics.convergence.ConvergenceWarning`.  The underlying
    sphere Krylov stall was cured by the GMRES ``restart``-sizing lineage
    (ERR-053, then #282 / #280 route (a)), not by #200, which is still open;
    the "History" section of
    ``tests/sn/verification/analytical/test_kinf_homogeneous.py`` carries the
    full account.  Any future exclusion here uses
    ``@pytest.mark.xfail(strict=True)`` so it retires itself — the imperative
    form can never report ``XPASS`` (issue #340, R5).
    """
    case = _get_continuous_case("2eg")
    mat_id = next(iter(case.problem.materials.keys()))
    mesh = _homogeneous_mesh(coord, mat_id)
    quad = _quadrature_for(coord)
    result = solve_sn(
        case.problem.materials, mesh, quad,
        inner_solver=inner_solver,
        max_outer=1000, keff_tol=1e-14, flux_tol=1e-12,
        max_inner=1000, inner_tol=1e-12,
    )
    keff, k_ref = result.keff, case.k_eff
    if keff is None or k_ref is None:
        pytest.fail(
            f"[{coord}/{inner_solver}] missing keff (solve={keff!r}, "
            f"reference={k_ref!r})"
        )
    np.testing.assert_allclose(
        float(keff), float(k_ref), rtol=1e-10,
        err_msg=(
            f"[{coord}/{inner_solver}] S8a driver re-attach moved the "
            f"converged k_inf: got {keff!r}, closed-form {k_ref!r}"
        ),
    )


def test_c5b_si_driver_iterate_stays_timed() -> None:
    """The SI driver re-attaches the timeless gain outputs into a TIMED iterate.

    The operator leaves are timeless ``FullField`` base arrows; the DRIVER
    carries the timed comonad.  Build the within-group ``SourceIteration``
    directly (``L+C`` resolvent, ``S`` + ``B`` gains) on a timed ``q_ext``: the
    step ``rhs = q_ext + S.apply(psi) + B.apply(psi)`` adds timeless gain
    outputs to the TIMED ``q_ext`` (the ``TimedFullField.__add__`` recombine
    re-attaches the timed type), and ``L.solve(rhs)`` returns the timed iterate.
    The returned iterate MUST be a ``TimedFullField`` — proving the re-attach
    landed on the driver, not the leaf.
    """
    from orpheus.numerics.iteration import SourceIteration
    from orpheus.sn.coupled_system import build_within_group_system
    from orpheus.transport.source_sinks import (
        AngularSourceSink,
        AngularBoundarySourceSink,
    )

    solver, _case = _solver_for("slab", "2eg")
    sn_mesh = solver.sn_mesh
    system = build_within_group_system(
        sn_mesh, solver.mat_xs, scattering_op=solver.scattering_op,
    )
    LC, (S, N2N, B) = system.implicit_operator, system.explicit_gains  # seedless slab record shape (§14.1)

    # A timed external source (the driver's comonad-carrying rhs).
    q_ext = TimedFullField(
        interior=AngularSourceSink.from_isotropic(
            np.ones((solver.ng, *sn_mesh.spatial_shape)), sn_mesh,
        ),
        boundary=AngularBoundarySourceSink.zeros(sn_mesh.angular_trace),
        _history=(),
        history_depth=2,
    )
    # Production always warm-starts the SI with a FLUX composite iterate (so
    # the gains dispatch on a flux bulk, not the source-role q_ext bulk); mirror
    # that here (SNSolver._solve_source_iteration passes a flux initial_guess).
    flux_seed = TimedFullField.zeros(
        interior=AngularFlux, boundary=AngularBoundaryFlux, space=sn_mesh.full_field_space, history_depth=2,
    )
    si = SourceIteration(LC.inverse(), S, B, max_iter=50, tol=1e-10)
    psi, _residuals = si.solve(q_ext, initial_guess=flux_seed)
    if not isinstance(psi, TimedFullField):
        pytest.fail(
            "the SI driver iterate must be a TimedFullField (the comonad lives "
            f"on the driver, re-attached from timeless leaf outputs); got "
            f"{type(psi).__name__}"
        )


# ═══════════════════════════════════════════════════════════════════════
# C5c — the driver-level ``advance`` history type-guard still fires.
# ═══════════════════════════════════════════════════════════════════════


def test_c5c_advance_type_guard_still_fires() -> None:
    """``TimedFullField.advance`` still rejects a mismatched-type new frame.

    The comonad verb (``advance`` — the only place the driver would re-stamp
    history) keeps its type-guard (``timed_full_field.py`` ~:306-317): a
    ``new_bulk`` whose type differs from the current bulk raises ``TypeError``.
    """
    solver, _case = _solver_for("slab", "2eg")
    sn_mesh = solver.sn_mesh
    state = TimedFullField.zeros(
        interior=AngularFlux, boundary=AngularBoundaryFlux, space=sn_mesh.full_field_space, history_depth=2,
    )
    # A bare ndarray is not an AngularFlux — the advance guard must reject it.
    with pytest.raises(TypeError, match="new_bulk type"):
        state.advance(np.zeros(3), state.boundary)  # type: ignore[arg-type]
