r"""G-adjoint reciprocity for the SN composite ``(L + C - B)`` (Wave O / O.2b R5).

The Hilbert adjoint ``A.H`` is the **G-adjoint** :math:`A^\dagger = G^{-1}
A^{\mathsf T} G`, NOT the plain Euclidean transpose. For the SN within-group
loss operator :math:`A = L + C - B` acting on the composite carrier
:math:`V_{\rm bulk} \oplus V_{\rm trace}`
(:class:`~orpheus.numerics.spaces.full_field_space.FullFieldSpace`), the metric
:math:`G` is block-diagonal:

* **bulk** :math:`G_{\rm bulk} = V_{\rm cell}\,w_n` (phase-space measure
  :math:`\mathrm{d}V\,\mathrm{d}\Omega`),
* **trace** :math:`G_{\rm trace} = |\Omega\cdot\hat n_f|\,w_n`
  (partial-current surface measure).

These tests pin the *defining* adjoint property — reciprocity
:math:`\langle A\psi,\varphi\rangle_G = \langle\psi, A^\dagger\varphi\rangle_G`
— and an L11 negative control proving the :math:`|\Omega\cdot\hat n|` weighting
is load-bearing.

Anti-R1 (metric-blind false green)
==================================

The reciprocity inner products are evaluated with an **independent**
``g_inner`` built directly from ``sn.trace.omega_dot_n`` + ``sn.quad.weights``
+ ``sn.volumes`` — NOT from the production metric under test
(``op.codomain.inner_product``). If ``op.H`` used a *wrong internal* metric
:math:`G'`, then :math:`\langle A\psi,\varphi\rangle_{G'} \ne \langle\psi,
G'^{-1}A^{\mathsf T}G'\varphi\rangle_{G_{\rm true}}` — the residual is O(1).
Evaluating with the true (independent) metric is what makes the test able to
see a wrong adjoint (test-architect verification plan §2, risk R1). A separate
cross-check pins ``g_inner == op.codomain.inner_product`` so a wrong *metric
population* is also caught.

These are algebraic identities (exact by construction), so the suite is
``@pytest.mark.foundation`` — no ``verifies()`` label (no theory equation; the
claim is a software/algebra invariant).

References
----------

* ``.claude/agent-memory/test-architect/phase4_o2b_g_adjoint_verification_plan.md``
  — §2 (reciprocity), §3 (L11 wrong-metric control), §4 (crosswalk).
* ``.claude/plans/glimmering-launching-lantern.md`` — Phase 4 / O.2b R5.
* Dense-probe oracle (structurally independent ``G^{-1}A^{\mathsf T}G`` fold):
  ``derivations/diagnostics/diag_p42_adjoint_oracle.py::validate_composite_adjoint``.
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.geometry import BC, CoordSystem, Mesh1D
from orpheus.sn.boundary_operator import SNBoundaryOperator
from orpheus.sn.geometry import SNMesh
from orpheus.sn.operator import CollisionOperator, StreamingOperator
from orpheus.numerics.quadrature import Quadrature
from orpheus.transport.fields.angular_flux import AngularFlux
from orpheus.transport.fields.boundary_flux import BoundaryFlux
from orpheus.transport.timed_full_field import TimedFullField
from tests.sn._test_helpers import placeholder_materials


# ═══════════════════════════════════════════════════════════════════════
# Mesh builders (homogeneous, reflective — B ≠ 0, the A_ss block live)
# ═══════════════════════════════════════════════════════════════════════


def _make_slab(nx: int = 4, R: float = 1.0, ng: int = 1, sigma: float = 0.5):
    quad = Quadrature.gauss_legendre(4)
    mesh = Mesh1D(
        edges=np.linspace(0.0, R, nx + 1),
        mat_ids=np.zeros(nx, dtype=int),
        coord=CoordSystem.CARTESIAN,
        bc_left=BC("reflective"),
        bc_right=BC("reflective"),
    )
    sn = SNMesh(mesh, quad, placeholder_materials(ng=ng))
    # Per-group-varying σ_t so the 2g row is non-degenerate in the group axis
    # (exercises the G_bulk broadcast over ng); rank-d (ng, *spatial).
    sig_t = np.stack(
        [np.full(sn.spatial_shape, sigma * (1.0 + 0.5 * g)) for g in range(ng)], axis=0
    )
    return sn, sig_t


def _make_sphere(nx: int = 4, R: float = 1.0, ng: int = 1, sigma: float = 0.5):
    quad = Quadrature.gauss_legendre(4)
    mesh = Mesh1D(
        edges=np.linspace(0.0, R, nx + 1),
        mat_ids=np.zeros(nx, dtype=int),
        coord=CoordSystem.SPHERICAL,
        bc_right=BC("reflective"),
    )
    sn = SNMesh(mesh, quad, placeholder_materials(ng=ng))
    sig_t = np.stack(
        [np.full(sn.spatial_shape, sigma * (1.0 + 0.5 * g)) for g in range(ng)], axis=0
    )
    return sn, sig_t


def _make_cyl(nx: int = 4, R: float = 1.0, ng: int = 1, sigma: float = 0.5):
    quad = Quadrature.level_symmetric(4)
    mesh = Mesh1D(
        edges=np.linspace(0.0, R, nx + 1),
        mat_ids=np.zeros(nx, dtype=int),
        coord=CoordSystem.CYLINDRICAL,
        bc_right=BC("reflective"),
    )
    sn = SNMesh(mesh, quad, placeholder_materials(ng=ng))
    sig_t = np.full((ng, nx), sigma)
    return sn, sig_t


_BUILDERS = {
    "slab": lambda: _make_slab(ng=1),
    "sphere": lambda: _make_sphere(ng=1),
    "cyl": lambda: _make_cyl(ng=1),
    "slab_2g": lambda: _make_slab(ng=2),
    "sphere_2g": lambda: _make_sphere(ng=2),
}


# ═══════════════════════════════════════════════════════════════════════
# Independent G inner product + random composite (structurally independent
# of the production FullFieldSpace metric — anti-R1)
# ═══════════════════════════════════════════════════════════════════════


def _bulk_measure(sn: SNMesh) -> np.ndarray:
    r"""G_bulk = V_cell · w_n on (N, 1, *spatial) — built from raw mesh data."""
    w_n = np.asarray(sn.quad.weights, dtype=float)
    V = np.asarray(sn.volumes, dtype=float)  # (*spatial,)
    # (N, 1) ordinate+group axes ⊗ (*spatial,) volume axes — rank-generic.
    w_b = w_n.reshape((w_n.shape[0], 1) + (1,) * V.ndim)
    return w_b * V[None, None]


def _trace_cosine_weight(sn: SNMesh, face_idx: int, *, with_cosine: bool) -> np.ndarray:
    r"""Per-ordinate trace weight for a face: ``|Ω·n|·w_n`` (true) or ``w_n`` (wrong)."""
    w_n = np.asarray(sn.quad.weights, dtype=float)
    if with_cosine:
        return np.abs(sn.trace.omega_dot_n[face_idx]) * w_n
    return w_n  # the L11 wrong metric: drops |Ω·n|


def _g_inner(a: TimedFullField, b: TimedFullField, sn: SNMesh, *,
             with_cosine: bool = True) -> float:
    r"""``⟨a,b⟩_G = Σ_bulk a·b·(V·w_n) + Σ_trace a·b·(|Ω·n|·w_n)``.

    Built directly from ``omega_dot_n`` / ``quad.weights`` / ``volumes`` — the
    structurally-independent reference inner product. ``with_cosine=False``
    drops the ``|Ω·n|`` factor (the L11 wrong-metric control).
    """
    bulk = float(np.sum(_bulk_measure(sn) * a.bulk.values * b.bulk.values))
    trace = 0.0
    for f_idx, face in enumerate(sn.trace.layout.faces):
        af = a.boundary.face_view(face)
        bf = b.boundary.face_view(face)
        w_face = _trace_cosine_weight(sn, f_idx, with_cosine=with_cosine)
        w_b = w_face.reshape((w_face.shape[0],) + (1,) * (af.ndim - 1))
        trace += float(np.sum(af * bf * w_b))
    return bulk + trace


def _random_composite(sn: SNMesh, rng: np.random.Generator) -> TimedFullField:
    r"""Random NON-FLAT composite (bulk + boundary both random per ordinate)."""
    N, ng = sn.quad.N, sn.ng
    bulk = AngularFlux.from_mesh(rng.standard_normal((N, ng, *sn.spatial_shape)), sn)
    boundary = BoundaryFlux(
        values=rng.standard_normal(int(sn.trace.layout.total_size)),
        space=sn.trace,
        mesh=sn,
    )
    return TimedFullField(bulk=bulk, boundary=boundary, _history=(), history_depth=2)


def _loss_operator(sn: SNMesh, sig_t: np.ndarray):
    r"""The within-group loss ``A = L + C - B`` (the boundary sibling ``-B`` live)."""
    L = StreamingOperator(sn, sig_t)
    C = CollisionOperator(sn, sig_t)
    B = SNBoundaryOperator(sn)
    return L + C - B


# ═══════════════════════════════════════════════════════════════════════
# §2 — G-adjoint reciprocity (the defining adjoint property)
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.foundation
@pytest.mark.parametrize("case", list(_BUILDERS))
def test_g_adjoint_reciprocity_full_block(case):
    r"""``⟨Aψ,φ⟩_G = ⟨ψ, A.H φ⟩_G`` for ``A = L + C - B`` (G = bulk V·w_n ⊕ trace |Ω·n|·w_n).

    The defining property of the G-adjoint, evaluated with the INDEPENDENT
    ``_g_inner`` (anti-R1). Random non-flat ψ, φ on both blocks so the
    partial-current trace metric (which varies across ordinates as
    :math:`|\mu_n|`) is not degenerate.
    """
    sn, sig_t = _BUILDERS[case]()
    A = _loss_operator(sn, sig_t)
    assert "apply_transpose" in A.capabilities  # the adjoint must be reachable

    rng = np.random.default_rng(2026)
    psi = _random_composite(sn, rng)
    phi = _random_composite(sn, rng)

    lhs = _g_inner(A.apply(psi), phi, sn)
    rhs = _g_inner(psi, A.H.apply(phi), sn)
    rel = abs(lhs - rhs) / (abs(lhs) + abs(rhs) + 1e-300)
    assert rel < 1e-12, (
        f"[{case}] G-reciprocity broken: ⟨Aψ,φ⟩_G={lhs:.6e} vs "
        f"⟨ψ,A.Hφ⟩_G={rhs:.6e} (rel={rel:.2e})"
    )


@pytest.mark.foundation
@pytest.mark.parametrize("case", list(_BUILDERS))
def test_full_field_space_metric_matches_independent_reference(case):
    r"""The production composite metric == the independent ``_g_inner``.

    Cross-check that ``op.codomain.inner_product`` (the
    :class:`~orpheus.numerics.spaces.full_field_space.FullFieldSpace` metric
    built from ``sn.full_field_space``) reproduces the reference built from raw
    ``omega_dot_n`` / ``volumes`` — catches a wrong metric *population* (a
    dropped factor, a double-count, a mis-broadcast over the group axis).
    """
    sn, sig_t = _BUILDERS[case]()
    space = sn.full_field_space
    rng = np.random.default_rng(7)
    a = _random_composite(sn, rng)
    b = _random_composite(sn, rng)
    ref = _g_inner(a, b, sn)
    prod = space.inner_product(a, b)
    rel = abs(ref - prod) / (abs(ref) + abs(prod) + 1e-300)
    assert rel < 1e-13, (
        f"[{case}] FullFieldSpace.inner_product {prod:.6e} != independent "
        f"reference {ref:.6e} (rel={rel:.2e}) — wrong metric population"
    )


# ═══════════════════════════════════════════════════════════════════════
# §3 — L11 negative control: the |Ω·n| weighting MUST be load-bearing
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.foundation
@pytest.mark.parametrize(
    "case",
    # slab + sphere break strongly (dense oracle: 6.4e-2 / 8.3e-1); cyl breaks
    # only marginally (1.9e-3) because its single reflective face carries less
    # of the reciprocity residual — so the decisive control lives on slab/sphere.
    ["slab", "sphere"],
)
def test_wrong_trace_metric_breaks_reciprocity(case):
    r"""Dropping ``|Ω·n|`` from the trace metric MUST break reciprocity (L11).

    ``A.H`` is built for the TRUE metric ``G``. Evaluating reciprocity under a
    WRONG metric ``G'`` (trace = ``w_n`` only, no ``|Ω·n|``) gives an O(1)
    residual — proving the ``|Ω·n|`` cosine weighting is load-bearing and that
    §2 is not a vacuous (metric-blind) green. Paired positive(§2)+negative(§3)
    per ERR-051 / vv-principles L11.
    """
    sn, sig_t = _BUILDERS[case]()

    # Precondition (anti-dud): the partial-current weight |Ω·n| must genuinely
    # VARY across ordinates on a reflective face, else dropping it is a no-op
    # and the control cannot fire (G-4 latent-dud lesson).
    spreads = [
        float(np.ptp(np.abs(sn.trace.omega_dot_n[f])))
        for f in range(sn.trace.omega_dot_n.shape[0])
    ]
    assert max(spreads) > 0.1, (
        f"[{case}] |Ω·n| spread {max(spreads):.3f} too small — the wrong-metric "
        f"control would be a dud (quadrature degeneracy)."
    )

    A = _loss_operator(sn, sig_t)
    rng = np.random.default_rng(2026)
    psi = _random_composite(sn, rng)
    phi = _random_composite(sn, rng)

    AH_phi = A.H.apply(phi)  # built for the TRUE metric G
    lhs_wrong = _g_inner(A.apply(psi), phi, sn, with_cosine=False)
    rhs_wrong = _g_inner(psi, AH_phi, sn, with_cosine=False)
    rel_wrong = abs(lhs_wrong - rhs_wrong) / (
        abs(lhs_wrong) + abs(rhs_wrong) + 1e-300
    )
    assert rel_wrong > 1e-3, (
        f"[{case}] wrong metric (drop |Ω·n|) did NOT break reciprocity "
        f"(rel={rel_wrong:.2e}); the |Ω·n|·w_n weighting is NOT load-bearing — "
        f"§2 proves nothing."
    )
