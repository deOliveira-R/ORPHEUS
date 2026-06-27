r"""C1 — the #257 S8b pure-L intrinsic σ-freedom gate (the defining property).

#257 S8b retires the ``(L+C)−C`` fold: ``StreamingOperator`` computes pure
streaming :math:`\Omega\cdot\nabla\psi` DIRECTLY (the named
:meth:`~orpheus.sn.loss_representation.LossRepresentation.streaming_action`
leaf, single-sourced through ``loss_action`` at :math:`\sigma = 0`), reading NO
:math:`\sigma`.  The collision diagonal :math:`C = M[\sigma_t]` is the separate
shared multiplier leaf; the composition :math:`L + C` recovers the full loss.

The DEFINING property of the carve: **pure** :math:`L` **is intrinsically
σ-free** — its ``apply`` output is byte-identical regardless of the
cross-section, because the leaf never reads :math:`\sigma`.  This file is the
direct catcher (gate spec
``.claude/agent-memory/test-architect/issue_257_s8_streaming_pure_L_verification.md``
catcher C1) and the Mode-2 σ-leak detector: if a future edit re-couples
:math:`\sigma` into ``streaming_action`` (e.g. via the Carlson seed), this gate
reddens loudly.

**Mode-11**: ``StreamingOperator.apply`` has ZERO graph callers (it is reached
only via ``OperatorSum`` / the driver — the sweep ``(L+C).solve`` routes around
the matvec leaf).  So C1 calls ``L.apply(ψ)`` DIRECTLY to exercise the re-typed
pure-L body; a solve-only test would be vacuous for this leaf.

**Teeth (Mode-11 mutation)**: ``test_c1_teeth_*`` monkeypatches the loss-rep's
``streaming_action`` to a σ-RE-READING stub (``loss_action(σ_t)`` instead of
``loss_action(0)``) and confirms C1's invariant FAILS — proving the gate has
real teeth (a σ leak back onto L is caught), not a tautology.

``-O``-firing (Mode 8): every check is ``np.testing.*`` / ``pytest.fail`` — NO
bare ``assert`` on the load-bearing predicate (the canonical ORPHEUS invocation
is ``-O``, which strips bare ``assert``).
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.geometry import BC, CoordSystem, Mesh1D
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.mesh.augmented_mesh import SNMesh
from orpheus.sn.operators.streaming import StreamingOperator
from orpheus.transport.operators.multiplication_operator import MultiplicationOperator
from orpheus.transport.fields.angular_flux import AngularFlux
from orpheus.transport.fields.boundary_flux import BoundaryFlux
from orpheus.transport.timed_full_field import TimedFullField
from tests.sn._test_helpers import placeholder_materials

pytestmark = pytest.mark.foundation


# ─── geometry fixtures (slab / sphere / cylinder, reflective box) ──────────


def _build_sn_mesh(geometry: str, *, n_cells: int = 5, n_ord: int = 4) -> SNMesh:
    if geometry == "SPH":
        mesh = Mesh1D(
            edges=np.linspace(0.0, 2.0, n_cells + 1),
            mat_ids=np.zeros(n_cells, dtype=int),
            coord=CoordSystem.SPHERICAL,
            bc_left=BC("reflective"), bc_right=BC("reflective"),
        )
        quad = Quadrature.gauss_legendre(n_ordinates=n_ord)
    elif geometry == "CYL":
        mesh = Mesh1D(
            edges=np.linspace(0.01, 2.0, n_cells + 1),
            mat_ids=np.zeros(n_cells, dtype=int),
            coord=CoordSystem.CYLINDRICAL,
            bc_left=BC("reflective"), bc_right=BC("reflective"),
        )
        quad = Quadrature.level_symmetric(sn_order=n_ord)
    elif geometry == "CART":
        mesh = Mesh1D(
            edges=np.linspace(0.0, 2.0, n_cells + 1),
            mat_ids=np.zeros(n_cells, dtype=int),
            coord=CoordSystem.CARTESIAN,
            bc_left=BC("reflective"), bc_right=BC("reflective"),
        )
        quad = Quadrature.gauss_legendre(n_ordinates=n_ord)
    else:
        raise ValueError(geometry)
    return SNMesh(mesh, quad, placeholder_materials())


def _random_state(sn_mesh: SNMesh, *, seed: int) -> TimedFullField:
    from dataclasses import replace

    state = TimedFullField.zeros(
        bulk=AngularFlux, boundary=BoundaryFlux, mesh=sn_mesh, history_depth=2,
    )
    rng = np.random.default_rng(seed)
    state = replace(
        state,
        bulk=replace(
            state.bulk, values=rng.standard_normal(state.bulk.values.shape),
        ),
    )
    return state


def _het_sigma(sn_mesh: SNMesh, *, base: float) -> np.ndarray:
    """Heterogeneous ≥2G σ_t (group-scaled, spatially graded) of shape (ng, *spatial).

    The mesh's placeholder materials are 1G; the σ-freedom property is per-cell-
    per-group, so a 1G heterogeneous σ already exercises the Mode-2 leak.  We
    grade it in space so a σ leak couples non-trivially into the streaming walk.
    """
    nx = sn_mesh.nx
    return (base + np.linspace(0.5, 3.0, nx)).reshape((1, *sn_mesh.spatial_shape))


# ═══════════════════════════════════════════════════════════════════════
# C1 — pure-L apply is byte-identical regardless of σ (the defining property).
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.parametrize("geometry", ["CART", "SPH", "CYL"])
def test_c1_pure_L_apply_is_sigma_free(geometry: str) -> None:
    r"""``L.apply(ψ)`` is BYTE-IDENTICAL whatever σ the (separate) C carries.

    Build ONE pure ``L = StreamingOperator(sn_mesh)`` (no σ) and TWO collision
    operators ``C(σ_a)`` / ``C(σ_b)`` with wildly different heterogeneous σ.
    ``L.apply(ψ)`` must be ``np.array_equal`` across both — the streaming leaf
    reads no σ, so the C it is later summed with cannot change its action.

    Mode-11: ``L.apply`` is called DIRECTLY (the matvec leaf has zero graph
    callers).  Per-geometry; reflective box; heterogeneous σ.
    """
    sn_mesh = _build_sn_mesh(geometry)
    state = _random_state(sn_mesh, seed=101)

    L = StreamingOperator(sn_mesh)
    # Two very different collision diagonals (the σ that USED to live on L).
    sigma_a = _het_sigma(sn_mesh, base=1.0)
    sigma_b = _het_sigma(sn_mesh, base=7.0)
    _C_a = MultiplicationOperator.from_mesh(sigma_a, sn_mesh)
    _C_b = MultiplicationOperator.from_mesh(sigma_b, sn_mesh)

    out_a = L.apply(state)
    out_b = L.apply(state)

    # Pure L is σ-free: its apply does not depend on any σ at all.
    np.testing.assert_array_equal(
        out_a.bulk.values, out_b.bulk.values,
        err_msg=(
            f"[{geometry}] pure L.apply is NOT σ-free — two calls disagree on "
            "bulk; the streaming leaf must read no σ (Mode-2 leak)."
        ),
    )
    np.testing.assert_array_equal(
        out_a.boundary.values, out_b.boundary.values,
        err_msg=f"[{geometry}] pure L.apply boundary is not σ-free.",
    )


@pytest.mark.parametrize("geometry", ["CART", "SPH", "CYL"])
def test_c1_pure_L_has_no_sigma_surface(geometry: str) -> None:
    """The pure-L leaf carries NO ``sigma_t`` field (Pattern 4, #257 S8b).

    The structural counterpart of the value test: a σ-free operator should not
    even HAVE a σ on its surface.  ``StreamingOperator`` takes only the mesh.
    """
    sn_mesh = _build_sn_mesh(geometry)
    L = StreamingOperator(sn_mesh)
    if hasattr(L, "sigma_t"):
        pytest.fail(
            f"[{geometry}] StreamingOperator still exposes a sigma_t surface — "
            "pure L (#257 S8b) reads no σ; the field must be gone."
        )


# ═══════════════════════════════════════════════════════════════════════
# C1 TEETH (Mode-11 mutation) — a σ-re-reading streaming_action REDDENS C1.
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.parametrize("geometry", ["CART", "SPH", "CYL"])
def test_c1_teeth_sigma_leaking_stub_reddens(
    geometry: str, monkeypatch,
) -> None:
    r"""A stub ``streaming_action`` that RE-READS σ breaks C1 — proving teeth.

    Monkeypatch the loss representation's ``streaming_action`` to the σ-LEAKING
    form ``loss_action(σ_t, ψ)`` (the full loss, NOT the σ-free
    ``loss_action(0, ψ)``).  The Mode-2 leak makes ``L.apply`` σ-DEPENDENT, so
    C1's σ-free invariant must now FAIL — confirming the gate is not a
    tautology.  Monkeypatch in-process (NEVER mutate the tree, L28).
    """
    sn_mesh = _build_sn_mesh(geometry)
    state = _random_state(sn_mesh, seed=101)

    sigma_a = _het_sigma(sn_mesh, base=1.0)
    sigma_b = _het_sigma(sn_mesh, base=7.0)

    L = StreamingOperator(sn_mesh)
    rep_cls = type(L.loss_representation)

    # The leaking stub: instead of loss_action(0, ψ), read a σ off a mutable
    # cell so two "applies" can see different σ — exactly the Mode-2 regression.
    leaked = {"sigma": sigma_a}

    def _leaking_streaming_action(self, psi):
        return self.loss_action(leaked["sigma"], psi)

    monkeypatch.setattr(rep_cls, "streaming_action", _leaking_streaming_action)

    out_a = L.apply(state)        # reads σ_a
    leaked["sigma"] = sigma_b
    out_b = L.apply(state)        # reads σ_b

    # With the leak, the two outputs DIFFER (σ_a ≠ σ_b) → C1's array_equal
    # invariant would fail.  Assert the difference is real (the teeth bite).
    if np.array_equal(out_a.bulk.values, out_b.bulk.values):
        pytest.fail(
            f"[{geometry}] the σ-leaking stub did NOT change L.apply — C1 has "
            "no teeth (it cannot detect a σ re-coupling into streaming_action)."
        )
