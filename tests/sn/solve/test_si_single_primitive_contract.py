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
from orpheus.sn.coupled_system import (
    CoupledInvertibleOperator,
    CoupledSweepOperator,
)
from orpheus.numerics.coupled_system import CoupledOperator
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
    carrying = sn_mesh.radial_characteristic_field_space is not None
    # (1) The step operator is the INVERSE of M — since B.2d, on a carrying
    # mesh the coupled joint sweep (CoupledSweepOperator over the
    # CoupledInvertibleOperator bridge; the fused (L+C) rides ``.fused``);
    # seedless, the plain SweepOperator over the InvertibleOperator. Same
    # concrete types both paths (#226 taxonomy step 3: the solver builds the
    # inverse, SourceIteration applies it).
    if carrying:
        assert isinstance(L_eig, CoupledSweepOperator)
        assert isinstance(L_fs, CoupledSweepOperator)
        assert isinstance(L_eig.inner, CoupledInvertibleOperator)
        assert isinstance(L_eig.inner.fused, InvertibleOperator)
        assert isinstance(L_fs.inner.fused, InvertibleOperator)
    else:
        assert isinstance(L_eig, SweepOperator)
        assert isinstance(L_fs, SweepOperator)
        assert isinstance(L_eig.inner, InvertibleOperator)
        assert isinstance(L_fs.inner, InvertibleOperator)
    assert type(L_eig.inner) is type(L_fs.inner)

    # (2) The coupling gains — SAME structural shape both paths. Since B.2d
    # the record's splitting A = M − N: on a carrying mesh N is ONE
    # CoupledOperator gain grid ``[[S+B_a, ∅],[Emission, B_b]]`` (the block-
    # native successor of the step-4c "own gain slot" — the Emission is a
    # grid BLOCK now, not a tuple entry); seedless, the (S, B_a) pair.
    expected_n = 1 if carrying else 2
    assert len(gains_eig) == expected_n
    assert len(gains_fs) == expected_n

    if carrying:
        # (2a/2b carrying) The ONE gain is the N grid; its (B,A) block is the
        # bulk→ray Emission (Fold∘K_iso∘integrate) and the grid spans systems
        # (SystemRole.COUPLED). Same structural grid both paths.
        N_eig, N_fs = gains_eig[0], gains_fs[0]
        assert isinstance(N_eig, CoupledOperator)
        assert isinstance(N_fs, CoupledOperator)
        assert N_eig.blocks[0][1] is None  # Seeding lives in M, never in N
        assert isinstance(
            N_eig.blocks[1][0], RadialCharacteristicEmission)
        assert isinstance(
            N_fs.blocks[1][0], RadialCharacteristicEmission)
        assert N_eig.system_role is SystemRole.COUPLED
    else:
        # (2a seedless) gain 0 = S: the BULK scattering coupling.
        S_eig, S_fs = gains_eig[0], gains_fs[0]
        assert isinstance(S_eig, ScatteringOperator)
        assert isinstance(S_fs, ScatteringOperator)

    # (2c) the BOUNDARY coupling. RULING P1's per-system direct sum is
    # STRUCTURAL since B.2d: on a carrying mesh B_a rides inside N's (A,A)
    # sum and B_b sits at the (B,B) slot (the ray corner, BOUNDARY-role);
    # seedless, the last gain is B_a itself. Same structural boundary
    # primitive both paths (never folded into the bulk S).
    if carrying:
        B_eig = gains_eig[0].blocks[1][1]
        B_fs = gains_fs[0].blocks[1][1]
    else:
        B_eig, B_fs = gains_eig[-1], gains_fs[-1]
    assert type(B_eig) is type(B_fs)
    assert B_eig.block_role is BlockRole.BOUNDARY
    assert B_fs.block_role is BlockRole.BOUNDARY
