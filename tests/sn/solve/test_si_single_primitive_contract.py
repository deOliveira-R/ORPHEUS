r"""Phase 1 (R1) single-primitive contract — structural SSoT spy.

This file pins the *mechanism* of the Phase 1 R1 fold (2026-06-04): the
fixed-source within-group solve :func:`~orpheus.sn.solver._solve_fixed_source_si`
and the eigenvalue within-group inner
:meth:`~orpheus.sn.solver.SNSolver._solve_source_iteration` are ONE primitive
(:class:`~orpheus.numerics.iteration.SourceIteration`) on ONE loss
decomposition ``(L+C, S, [A_BA,] B)`` — the invertible resolvent ``L+C`` plus
the scattering ``S``, the curvilinear bulk→ray coupling ``A_BA`` (seed-carrying
meshes only — campaign step 4c, THE LIFT), and the boundary ``B`` coupling gains
— differing only in ``q_ext`` (external vs fission source) and the returned
contract.

Why a *structural* test and not a numerical one
================================================

A pure-numerical agreement test ("does fixed-source SI give the same ``φ`` as
the retired hand-rolled loop?") would PASS even if the two implementations
stayed separate parallel shapes — it verifies the green, NOT the architecture
(lesson L13).  Coding-elegance Pattern 2 (single source of truth) is an
architectural invariant; the only way to pin it is to prove the SAME primitive
is constructed.  This spy captures the ``(L, *gains)`` operands handed to
``SourceIteration.__init__`` from BOTH paths and asserts they are the same
structural decomposition (resolvent + scattering gain + boundary gain).

If a future edit re-introduces a hand-rolled ``for n_inner`` loop on the
fixed-source path, that path constructs ZERO ``SourceIteration`` instances and
this test fails — the twin-path regression is caught by construction.

Foundation-tagged: this is a software SSoT invariant, NOT an equation
``:label:`` claim (no V&V level, no ``verifies``).
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.derivations.common.xs_library import get_mixture
from orpheus.geometry import BC, CoordSystem, Mesh1D
from orpheus.sn import solve_sn_fixed_source
from orpheus.sn.mesh.augmented_mesh import SNMesh
from orpheus.sn.solver import SNSolver
from orpheus.sn.operators.streaming import InvertibleOperator
from orpheus.sn.operators.sweep_operator import SweepOperator
from orpheus.transport.operators.scattering import ScatteringOperator
from orpheus.sn.operators.radial_characteristic import RadialCharacteristicEmission
from orpheus.numerics.operator import BlockRole, SystemRole
from orpheus.numerics import iteration as _iteration
from orpheus.numerics.quadrature import Quadrature
from orpheus.transport.source_sinks import AngularSourceSink

pytestmark = pytest.mark.foundation


# ── Heterogeneous 2-group mesh builders (non-trivial operands) ────────
#
# ≥2G + 2-region so the captured triple's operands are non-degenerate
# (lesson L2 — a 1G / homogeneous case would not exercise the operand
# structure).  Fuel "A" (fission) + moderator "B".


def _slab_2g_2region(nx: int = 6) -> tuple:
    """1-D slab, 2 regions (fuel | moderator), 2G, vacuum BCs."""
    mat_ids = np.array([1] * (nx // 2) + [0] * (nx - nx // 2), dtype=int)
    mesh = Mesh1D(
        edges=np.linspace(0.0, 2.0, nx + 1),
        mat_ids=mat_ids,
        coord=CoordSystem.CARTESIAN,
        bc_left=BC("vacuum"),
        bc_right=BC("vacuum"),
    )
    quad = Quadrature.gauss_legendre(n_ordinates=4)
    materials = {1: get_mixture("A", "2g"), 0: get_mixture("B", "2g")}
    return mesh, quad, materials


def _sphere_2g_2region(nx: int = 6) -> tuple:
    """1-D sphere, 2 regions, 2G, reflective centre + vacuum outer.

    The curvilinear sibling: ``inner_solver="source_iteration"`` is passed
    explicitly so the fixed-source entry routes to ``_solve_fixed_source_si``
    (overriding the curvilinear ``"krylov"`` default) — i.e. it exercises the
    SAME SI primitive on a curvilinear mesh.
    """
    mat_ids = np.array([1] * (nx // 2) + [0] * (nx - nx // 2), dtype=int)
    mesh = Mesh1D(
        edges=np.linspace(0.0, 2.0, nx + 1),
        mat_ids=mat_ids,
        coord=CoordSystem.SPHERICAL,
        bc_left=BC("reflective"),
        bc_right=BC("vacuum"),
    )
    quad = Quadrature.gauss_legendre(n_ordinates=8)
    materials = {1: get_mixture("A", "2g"), 0: get_mixture("B", "2g")}
    return mesh, quad, materials


_BUILDERS = {"slab": _slab_2g_2region, "sphere": _sphere_2g_2region}


@pytest.mark.parametrize("case", list(_BUILDERS))
def test_fixed_source_si_and_eigenvalue_inner_share_one_primitive(
    case: str, monkeypatch,
) -> None:
    """Both within-group SI paths build the SAME ``SourceIteration`` decomposition.

    Captures the ``(L, *gains)`` operands at ``SourceIteration.__init__`` from
    the eigenvalue inner and the migrated fixed-source SI, and asserts the
    structural identity of the decomposition (mechanism, not numerics).
    """
    mesh, quad, materials = _BUILDERS[case]()
    sn_mesh = SNMesh(mesh, quad, materials)

    # Spy: wrap SourceIteration.__init__ to record the (resolvent, *gains)
    # operands, then delegate to the real __init__ so the solve still runs.
    captured: list[tuple] = []
    real_init = _iteration.SourceIteration.__init__

    def _spy_init(self, L, *gains, **kwargs):
        captured.append((L, gains))
        real_init(self, L, *gains, **kwargs)

    monkeypatch.setattr(_iteration.SourceIteration, "__init__", _spy_init)

    # (a) Eigenvalue inner — the canonical SourceIteration consumer.
    solver = SNSolver(sn_mesh, max_inner=4)
    ng = solver.ng
    fission_source = np.ones((ng, *sn_mesh.spatial_shape))
    flux = np.ones((ng, *sn_mesh.spatial_shape))
    solver._solve_source_iteration(fission_source, flux)
    assert len(captured) == 1, (
        "the eigenvalue inner must build exactly one SourceIteration; "
        f"captured {len(captured)}."
    )
    L_eig, gains_eig = captured[-1]

    # (b) Migrated fixed-source SI — MUST route through the SAME primitive.
    external = AngularSourceSink.from_isotropic(
        np.ones((ng, *sn_mesh.spatial_shape)), sn_mesh,
    ).values
    solve_sn_fixed_source(
        materials=materials, mesh=mesh, quadrature=quad,
        external_source=external,
        inner_solver="source_iteration",
        max_inner=4,
    )
    assert len(captured) == 2, (
        "the fixed-source SI path must build a SourceIteration — a "
        "hand-rolled `for n_inner` loop would capture ZERO (the Pattern 2 "
        "single-source-of-truth regression this test exists to catch)."
    )
    L_fs, gains_fs = captured[-1]

    # ── Structural identity of the decomposition (no numerical tolerance) ──
    # (1) The step operator is the INVERSE of (L + C) — a SweepOperator
    # whose ``inner`` is the InvertibleOperator forward, same concrete
    # types both paths (#226 taxonomy step 3: the solver builds the
    # inverse, SourceIteration applies it — the forward identity moved
    # one level in, onto ``.inner``).
    assert isinstance(L_eig, SweepOperator)
    assert isinstance(L_fs, SweepOperator)
    assert isinstance(L_eig.inner, InvertibleOperator)
    assert isinstance(L_fs.inner, InvertibleOperator)
    assert type(L_eig.inner) is type(L_fs.inner)

    # (2) The coupling gains — SAME structural tuple both paths (Wave O O.2a:
    # the honest (L+C, S, [A_BA,] B); the transitional S+B fold is RETIRED, B a
    # SEPARATE gain). On a seed-carrying (curvilinear) mesh the gains ALSO carry
    # A_BA — the scattering bulk→ray coupling as its OWN lagged gain (campaign
    # step 4c, THE LIFT; user-ruled "own gain slot", the #208 pattern that
    # separated B from S extended to the ψ½ ray). Both SI paths share it, so the
    # single-primitive contract now pins the THREE-gain sphere shape too.
    carrying = sn_mesh.radial_characteristic_space is not None
    expected_n = 3 if carrying else 2
    assert len(gains_eig) == expected_n
    assert len(gains_fs) == expected_n

    # (2a) gain 0 = S: the BULK scattering coupling.
    S_eig, S_fs = gains_eig[0], gains_fs[0]
    assert isinstance(S_eig, ScatteringOperator)
    assert isinstance(S_fs, ScatteringOperator)

    # (2b) On a carrying mesh, gain 1 = A_BA: the bulk→ray emission
    # (RadialCharacteristicEmission = Fold∘K_iso∘integrate; SystemRole.COUPLED),
    # riding the B.2b transient FullField-gain adapter (the block itself speaks
    # System B's composite codomain — the wraps-predicate reads the wrapped
    # ``_emission``). Absent on a seedless mesh (no ψ½ ray). Same primitive
    # both paths.
    if carrying:
        A_BA_eig, A_BA_fs = gains_eig[1], gains_fs[1]
        assert isinstance(
            getattr(A_BA_eig, "_emission", None), RadialCharacteristicEmission)
        assert isinstance(
            getattr(A_BA_fs, "_emission", None), RadialCharacteristicEmission)
        assert A_BA_eig.system_role is SystemRole.COUPLED

    # (2c) last gain = B: the augmented BOUNDARY coupling. On a seed-carrying
    # (sphere) mesh it is the direct sum B_a + B_b (System A trace ⊕ System B ray
    # corner — an OperatorSum, RULING P1; the un-weld of the old welded
    # SNBoundaryOperator). The |Ω·n|·w trace adjoint metric lives on B_a. The
    # single-primitive contract: SI and eigenvalue share the SAME structural
    # boundary primitive — same type both paths, carrying the BOUNDARY block-role
    # (never folded into the bulk S).
    B_eig, B_fs = gains_eig[-1], gains_fs[-1]
    assert type(B_eig) is type(B_fs)
    assert B_eig.block_role is BlockRole.BOUNDARY
    assert B_fs.block_role is BlockRole.BOUNDARY
