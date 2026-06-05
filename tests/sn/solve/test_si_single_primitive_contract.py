r"""Phase 1 (R1) single-primitive contract — structural SSoT spy.

This file pins the *mechanism* of the Phase 1 R1 fold (2026-06-04): the
fixed-source within-group solve :func:`~orpheus.sn.solver._solve_fixed_source_si`
and the eigenvalue within-group inner
:meth:`~orpheus.sn.solver.SNSolver._solve_source_iteration` are ONE primitive
(:class:`~orpheus.numerics.iteration.SourceIteration`) on ONE operator triple
``(L+C, S+B, F=0)`` — differing only in ``q_ext`` (external vs fission source)
and the returned contract.

Why a *structural* test and not a numerical one
================================================

A pure-numerical agreement test ("does fixed-source SI give the same ``φ`` as
the retired hand-rolled loop?") would PASS even if the two implementations
stayed separate parallel shapes — it verifies the green, NOT the architecture
(lesson L13).  Coding-elegance Pattern 2 (single source of truth) is an
architectural invariant; the only way to pin it is to prove the SAME primitive
is constructed.  This spy captures the ``(L, S, F)`` operands handed to
``SourceIteration.__init__`` from BOTH paths and asserts they are the same
structural triple.

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
from orpheus.sn.geometry import SNMesh
from orpheus.sn.solver import SNSolver, _zero_within_group_fission
from orpheus.sn.operator import InvertibleOperator
from orpheus.numerics import iteration as _iteration
from orpheus.numerics.operator import OperatorSum, ZeroOperator
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
    """Both within-group SI paths build the SAME ``SourceIteration`` triple.

    Captures the ``(L, S, F)`` operands at ``SourceIteration.__init__`` from
    the eigenvalue inner and the migrated fixed-source SI, and asserts the
    structural identity of the triple (mechanism, not numerics).
    """
    mesh, quad, materials = _BUILDERS[case]()
    sn_mesh = SNMesh(mesh, quad, materials)
    nx, ny = sn_mesh.nx, sn_mesh.ny

    # Spy: wrap SourceIteration.__init__ to record the operand triple, then
    # delegate to the real __init__ so the solve still runs.
    captured: list[tuple] = []
    real_init = _iteration.SourceIteration.__init__

    def _spy_init(self, L, S, F, **kwargs):
        captured.append((L, S, F))
        real_init(self, L, S, F, **kwargs)

    monkeypatch.setattr(_iteration.SourceIteration, "__init__", _spy_init)

    # (a) Eigenvalue inner — the canonical SourceIteration consumer.
    solver = SNSolver(sn_mesh, max_inner=4)
    ng = solver.ng
    fission_source = np.ones((ng, nx, ny))
    flux = np.ones((ng, nx, ny))
    solver._solve_source_iteration(fission_source, flux)
    assert len(captured) == 1, (
        "the eigenvalue inner must build exactly one SourceIteration; "
        f"captured {len(captured)}."
    )
    L_eig, S_eig, F_eig = captured[-1]

    # (b) Migrated fixed-source SI — MUST route through the SAME primitive.
    external = AngularSourceSink.from_isotropic(
        np.ones((ng, nx, ny)), sn_mesh,
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
    L_fs, S_fs, F_fs = captured[-1]

    # ── Structural identity of the triple (no numerical tolerance) ────
    # (1) L = (L + C): an InvertibleOperator, same concrete type both paths.
    assert isinstance(L_eig, InvertibleOperator)
    assert isinstance(L_fs, InvertibleOperator)
    assert type(L_eig) is type(L_fs)

    # (2) S = (S + B): the scattering-with-boundary fold — a *plain*
    # OperatorSum (NOT the InvertibleOperator subclass that L is).
    assert type(S_eig) is OperatorSum
    assert type(S_fs) is OperatorSum

    # (3) F = 0_wg: a ZeroOperator carrying the SAME within-group-fission
    # codomain zero (module-level singleton — identical identity both paths).
    assert isinstance(F_eig, ZeroOperator)
    assert isinstance(F_fs, ZeroOperator)
    assert F_eig._codomain_zero is _zero_within_group_fission
    assert F_fs._codomain_zero is _zero_within_group_fission
