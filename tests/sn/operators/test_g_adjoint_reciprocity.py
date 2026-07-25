r"""G-adjoint reciprocity for the SN composites ``(L + C - B)`` and the full
within-group loss ``(L + C - S - B)`` (Wave O / O.2b R5; full-loss rows Phase
2.5 S0.2, #280).

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
``g_inner`` built directly from ``sn.angular_trace.omega_dot_n`` + ``sn.quad.weights``
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

``-O``-safe (vv Mode 8): every gate is a ``pytest.fail`` function call — a bare
``assert`` is stripped to a NO-OP under the canonical ``python -O``. Migrated off
bare ``assert`` in S6.3 (issue #222) so the curvilinear (sphere/cyl) angular
second-triangular-factor reciprocity — the gate S6.3's walk-move leans on —
ACTUALLY FIRES under ``-O`` instead of being a false green.

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

from scipy.sparse import csr_matrix

from orpheus.derivations.common.xs_library import make_mixture
from orpheus.geometry import BC, CoordSystem, Mesh1D
from orpheus.sn.operators.boundary import SNBoundaryOperator
from orpheus.sn.mesh.augmented_mesh import SNMesh
from orpheus.sn.solver import SNSolver
from orpheus.sn.operators.streaming import StreamingOperator
from orpheus.transport.operators.multiplication_operator import MultiplicationOperator
from orpheus.numerics.quadrature import Quadrature
from orpheus.transport.fields.angular_flux import AngularFlux
from orpheus.transport.fields.angular_boundary_flux import AngularBoundaryFlux
from orpheus.transport.timed_full_field import TimedFullField
from tests.sn._test_helpers import (
    cart2d_2g_nonsquare,
    g_bulk_measure,
    g_inner,
    g_trace_cosine_weight,
    placeholder_materials,
)


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


def _make_cyl_product(nx: int = 4, R: float = 1.0, ng: int = 1, sigma: float = 0.5):
    r"""Cylinder on the equispaced PRODUCT rule — the DEGENERATE-class row.

    ``Quadrature.product(n_mu=2, n_phi=4)`` samples φ ∈ {0, π/2, π, 3π/2},
    so the φ = π/2, 3π/2 ordinates carry :math:`|\mu_x| \approx 6\cdot
    10^{-17}` — genuinely degenerate pure-azimuthal ordinates whose cell
    balance is volumetric (the matvec's degenerate branch, no face march).
    The :func:`_make_cyl` ``level_symmetric`` rule has NO such ordinates,
    so every pre-2.5a reciprocity row was structurally BLIND to the
    degenerate rows of the transpose (which silently dropped them —
    the #280 2.5a completion; vv Mode 7: this builder ACTIVATES the term
    the level-symmetric rows null by quadrature choice).
    """
    quad = Quadrature.product(n_mu=2, n_phi=4)
    mesh = Mesh1D(
        edges=np.linspace(0.0, R, nx + 1),
        mat_ids=np.zeros(nx, dtype=int),
        coord=CoordSystem.CYLINDRICAL,
        bc_right=BC("reflective"),
    )
    sn = SNMesh(mesh, quad, placeholder_materials(ng=ng))
    sig_t = np.stack(
        [np.full(sn.spatial_shape, sigma * (1.0 + 0.5 * g)) for g in range(ng)], axis=0
    )
    return sn, sig_t


def _make_ld_slab(ng: int = 2, sigma: float = 0.5):
    r"""LD slab, NON-UNIFORM h + group- AND space-varying σ_t (#310 C2).

    The spec §4.4 reciprocity config: heterogeneous ≥2G on non-uniform
    widths, so ``M = diag(h, θh)`` varies cell-to-cell (the mass-order
    discriminant is live) and the random composite's slope moments are
    genuinely exercised against a non-flat coefficient field.
    """
    from orpheus.transport.spatial.linear_discontinuous import (
        LinearDiscontinuous,
    )

    quad = Quadrature.gauss_legendre(4)
    edges = np.array([0.0, 0.17, 0.45, 0.62, 1.0])       # non-uniform h
    mesh = Mesh1D(
        edges=edges,
        mat_ids=np.zeros(4, dtype=int),
        coord=CoordSystem.CARTESIAN,
        bc_left=BC("reflective"),
        bc_right=BC("reflective"),
    )
    sn = SNMesh(mesh, quad, placeholder_materials(ng=ng),
                scheme=LinearDiscontinuous())
    nx = sn.spatial_shape[0]
    space_factor = 1.0 + 0.3 * np.arange(nx) / nx        # spatially het
    sig_t = np.stack(
        [sigma * (1.0 + 0.5 * g) * space_factor for g in range(ng)], axis=0,
    )
    return sn, sig_t


def _make_cart2d(ng: int = 2, sigma: float = 0.5):
    r"""2-D Cartesian DD, reflective nonsquare (#310 C4 — the multi-D row).

    The G-metric reciprocity on cart2d is NOT slab-degenerate (the L19
    correction): the trace weight ``|Ω·n|·w_n`` varies per ordinate AND per
    face normal, so a Euclidean-trace ``.H`` reds here at O(1) exactly as
    on the curvilinear rows.  σ_t per-group-varying (anti-#3).
    """
    sn = cart2d_2g_nonsquare()
    sig_t = np.stack(
        [np.full(sn.spatial_shape, sigma * (1.0 + 0.5 * g)) for g in range(ng)],
        axis=0,
    )
    return sn, sig_t


def _make_ld_2d(ng: int = 2, sigma: float = 0.5):
    r"""LD 2-D Cartesian, reflective NONSQUARE + NON-UNIFORM h, het σ (#310 C5).

    The spec §6.1 reciprocity config — both #310 moment-metric faces live
    at once, each on its own committed convention:

    * the BULK metric carries the d=2 moment-mass Kronecker
      ``V·w_n ⊗ [1, θ, θ, θ²]`` (ruling 3's d-generic kron, first
      instantiated at d=2 here — :func:`g_bulk_measure`'s independent
      spelling vs the production ``moment_mass_diagonal``, pinned equal by
      the metric cross-check row);
    * the TRACE is the moment-resolved ``[avg, transverse-slope]`` face
      (#251), whose partial-current metric ``|Ω·n|·w_n`` broadcasts
      UNIFORMLY over the moment axis — the purely-ANGULAR Wave-O
      convention: the trace metric carries no spatial measure at all (no
      face area, hence no spatial moment mass; the θ-mass is a bulk
      phase-space-measure concept).  ``g_inner``'s trace term and the
      production :func:`_build_trace_metric_weights` spell the same
      broadcast, so the cross-check row pins the convention too.

    Non-uniform on BOTH axes (``M = diag`` varies cell-to-cell), reflective
    (``B`` live on the moment-resolved trace), σ_t group- AND space-varying.
    """
    from orpheus.geometry import Mesh2D
    from orpheus.transport.spatial.linear_discontinuous import (
        LinearDiscontinuous,
    )

    geom = Mesh2D(
        edges_x=np.array([0.0, 0.17, 0.45, 0.62, 1.0]),
        edges_y=np.array([0.0, 0.33, 0.8, 1.4]),
        mat_map=np.zeros((4, 3), dtype=int),
        bc_xmin=BC("reflective"), bc_xmax=BC("reflective"),
        bc_ymin=BC("reflective"), bc_ymax=BC("reflective"),
    )
    sn = SNMesh(
        geom, Quadrature.level_symmetric(2), placeholder_materials(ng=ng),
        scheme=LinearDiscontinuous(),
    )
    nx, ny = sn.spatial_shape
    space_factor = (
        1.0 + 0.3 * (np.arange(nx)[:, None] + 2.0 * np.arange(ny)[None, :])
        / (nx + ny)
    )
    sig_t = np.stack(
        [sigma * (1.0 + 0.5 * g) * space_factor for g in range(ng)], axis=0,
    )
    return sn, sig_t


_BUILDERS = {
    "slab": lambda: _make_slab(ng=1),
    "sphere": lambda: _make_sphere(ng=1),
    "cyl": lambda: _make_cyl(ng=1),
    "slab_2g": lambda: _make_slab(ng=2),
    "sphere_2g": lambda: _make_sphere(ng=2),
    "cyl_product_2g": lambda: _make_cyl_product(ng=2),
    "ld_slab_2g": lambda: _make_ld_slab(ng=2),
    "cart2d_2g": lambda: _make_cart2d(ng=2),
    "ld_2d_2g": lambda: _make_ld_2d(ng=2),
}


# ═══════════════════════════════════════════════════════════════════════
# Independent G inner product + random composite (structurally independent
# of the production FullFieldSpace metric — anti-R1).  The oracle trio
# g_bulk_measure / g_trace_cosine_weight / g_inner is SHARED with the
# tests/sn/solve adjoint batteries — single-sourced in
# tests/sn/_test_helpers.py (#276 A4 sweep).
# ═══════════════════════════════════════════════════════════════════════


def _random_composite(
    sn: SNMesh, rng: np.random.Generator,
) -> TimedFullField:
    r"""Random NON-FLAT 2-block composite (bulk + boundary both random).

    B.2d: System A's composite carries no ψ½ block — the reciprocity gates
    below pin the (A,A) block's G-adjoint on the 2-block metric (the walk's
    no-leg surfaces zero-substitute the seed and discard the pullback, i.e.
    the honest ray-decoupled block); the COUPLED reciprocity including
    System B's ``G_sd = V_cell`` gauge is the grid ``.H`` gate
    (``test_psi_half_coupling::TestCoupledBuilder``).
    """
    N, ng = sn.quad.N, sn.ng
    per_axis = sn.scheme.spatial_basis_per_axis
    if per_axis > 1:
        # A multi-moment closure (LD): the bulk carries the trailing 2^d
        # spatial-moment axis with RANDOM (non-zero) slope moments — the
        # §4.4 anisotropy audit's requirement (an all-flat suite is blind
        # to a dropped/mis-signed slope row).
        tail_size = per_axis ** sn.ndim
        bulk = AngularFlux.from_mesh(
            rng.standard_normal((N, ng, *sn.spatial_shape, tail_size)), sn,
            spatial_moments=per_axis,
        )
    else:
        bulk = AngularFlux.from_mesh(
            rng.standard_normal((N, ng, *sn.spatial_shape)), sn,
        )
    boundary = AngularBoundaryFlux(
        values=rng.standard_normal(int(sn.angular_trace.layout.total_size)),
        space=sn.angular_trace,
        mesh=sn,
    )
    return TimedFullField(
        interior=bulk, boundary=boundary,
        _history=(), history_depth=2,
    )


def _loss_operator(sn: SNMesh, sig_t: np.ndarray):
    r"""The within-group loss ``A = L + C - B`` (the boundary sibling ``-B`` live)."""
    L = StreamingOperator(sn)
    C = MultiplicationOperator.from_mesh(sig_t, sn)
    B = SNBoundaryOperator(sn)
    return L + C - B


# ═══════════════════════════════════════════════════════════════════════
# §2 — G-adjoint reciprocity (the defining adjoint property)
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.foundation
@pytest.mark.catches("ERR-066")
@pytest.mark.parametrize("case", list(_BUILDERS))
def test_g_adjoint_reciprocity_full_block(case):
    r"""``⟨Aψ,φ⟩_G = ⟨ψ, A.H φ⟩_G`` for ``A = L + C - B`` (G = bulk V·w_n ⊕ trace |Ω·n|·w_n).

    The defining property of the G-adjoint, evaluated with the INDEPENDENT
    ``g_inner`` (anti-R1). Random non-flat ψ, φ on both blocks so the
    partial-current trace metric (which varies across ordinates as
    :math:`|\mu_n|`) is not degenerate.

    The ``cyl_product_2g`` row is the ERR-066 catcher: it activates the
    degenerate pure-azimuthal class every ``level_symmetric`` row nulls by
    rule choice, so a transpose that drops the volumetric branch reds here
    at O(1).
    """
    sn, sig_t = _BUILDERS[case]()
    A = _loss_operator(sn, sig_t)
    if not A.is_adjointable:  # the adjoint must be reachable (carve P4 rewire)
        pytest.fail(f"[{case}] adjoint unreachable: A.is_adjointable is False")

    rng = np.random.default_rng(2026)
    psi = _random_composite(sn, rng)
    phi = _random_composite(sn, rng)

    lhs = g_inner(A.apply(psi), phi, sn)
    rhs = g_inner(psi, A.H.apply(phi), sn)
    rel = abs(lhs - rhs) / (abs(lhs) + abs(rhs) + 1e-300)
    if not rel < 1e-12:
        pytest.fail(
            f"[{case}] G-reciprocity broken: ⟨Aψ,φ⟩_G={lhs:.6e} vs "
            f"⟨ψ,A.Hφ⟩_G={rhs:.6e} (rel={rel:.2e})"
        )


@pytest.mark.foundation
@pytest.mark.parametrize("case", list(_BUILDERS))
def test_full_field_space_metric_matches_independent_reference(case):
    r"""The production composite metric == the independent ``g_inner``.

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
    ref = g_inner(a, b, sn)
    prod = space.inner_product(a, b)
    rel = abs(ref - prod) / (abs(ref) + abs(prod) + 1e-300)
    if not rel < 1e-13:
        pytest.fail(
            f"[{case}] FullFieldSpace.inner_product {prod:.6e} != independent "
            f"reference {ref:.6e} (rel={rel:.2e}) — wrong metric population"
        )


@pytest.mark.foundation
def test_ld_moment_mass_metric_is_load_bearing(monkeypatch):
    r"""#310 C2 ruling 3, the stabiliser proof: a slope-row transpose error
    is VISIBLE to the θ-carrying metric and EXACTLY INVISIBLE to a
    slope-ghost metric.

    Plant a slope-row sign flip in LD's registered batch VJP (the ψ̂
    cotangent negated — the M-R1c slope-row error class) and measure the
    reciprocity defect ``|⟨Aψ,φ⟩_m − ⟨ψ, A.H φ⟩_m|`` under two bulk
    measures:

    * the θ-metric (``V·w_n ⊗ diag(1, θ)``): the mutation moves the defect
      O(1) — the committed reciprocity row has teeth on the slope rows;
    * a slope-GHOST metric (``V·w_n ⊗ diag(1, 0)`` — the L18 ghost-G
      family): the mutated defect equals the unmutated defect EXACTLY —
      the error class sits in the ghost metric's stabiliser (Mode 12), so
      a metric that drops the moment mass is structurally BLIND here.

    The asymmetry IS the proof the moment mass is load-bearing in the
    reciprocity gate (spec §4.2); the metric-FREE value catchers are the
    Euclidean oracles (the SymPy/dense ``AᵀM⁻¹`` gates + G2).
    """
    from orpheus.transport.spatial.linear_discontinuous import (
        LinearDiscontinuous,
    )

    sn, sig_t = _make_ld_slab(ng=2)
    A = _loss_operator(sn, sig_t)
    rng = np.random.default_rng(3104)
    psi = _random_composite(sn, rng)
    phi = _random_composite(sn, rng)

    def _ghost_inner(a, b) -> float:
        # bulk V·w_n ⊗ diag(1, 0) + the TRUE trace metric.
        w_n = np.asarray(sn.quad.weights, dtype=float)
        V = np.asarray(sn.volumes, dtype=float)
        base = (w_n.reshape((w_n.shape[0], 1, 1)) * V[None, None])[..., None]
        ghost = base * np.array([1.0, 0.0])
        bulk = float(np.sum(ghost * a.interior.values * b.interior.values))
        trace = 0.0
        for f_idx, face in enumerate(sn.angular_trace.layout.faces):
            w_face = g_trace_cosine_weight(sn, f_idx, with_cosine=True)
            af, bf = a.boundary.face_view(face), b.boundary.face_view(face)
            trace += float(np.sum(af * bf * w_face[:, None]))
        return bulk + trace

    def _defect(inner) -> float:
        return abs(inner(A.apply(psi), phi) - inner(psi, A.H.apply(phi)))

    theta_clean = _defect(lambda x, y: g_inner(x, y, sn))
    ghost_clean = _defect(_ghost_inner)

    orig = LinearDiscontinuous.residual_kernel_batch_transpose

    def slope_flipped(self, **kw):
        psi_bar_cot, psi_in_cots = orig(self, **kw)
        psi_bar_cot = psi_bar_cot.copy()
        psi_bar_cot[..., 1] = -psi_bar_cot[..., 1]   # the planted slope flip
        return psi_bar_cot, psi_in_cots

    monkeypatch.setattr(
        LinearDiscontinuous, "residual_kernel_batch_transpose", slope_flipped,
    )
    theta_mut = _defect(lambda x, y: g_inner(x, y, sn))
    ghost_mut = _defect(_ghost_inner)

    scale = abs(float(g_inner(psi, psi, sn)))
    if not theta_mut > 1e-6 * scale:
        pytest.fail(
            f"θ-metric reciprocity did not red on the slope flip "
            f"(defect {theta_mut:.3e}, clean {theta_clean:.3e}) — the "
            "moment-mass metric row has no slope teeth"
        )
    if not np.isclose(ghost_mut, ghost_clean, rtol=0.0, atol=1e-12 * scale):
        pytest.fail(
            f"slope-ghost metric SAW the slope flip (mutated {ghost_mut:.3e} "
            f"vs clean {ghost_clean:.3e}) — the stabiliser asymmetry proof "
            "is broken"
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
        float(np.ptp(np.abs(sn.angular_trace.omega_dot_n[f])))
        for f in range(sn.angular_trace.omega_dot_n.shape[0])
    ]
    if not max(spreads) > 0.1:
        pytest.fail(
            f"[{case}] |Ω·n| spread {max(spreads):.3f} too small — the wrong-metric "
            f"control would be a dud (quadrature degeneracy)."
        )

    A = _loss_operator(sn, sig_t)
    rng = np.random.default_rng(2026)
    psi = _random_composite(sn, rng)
    phi = _random_composite(sn, rng)

    AH_phi = A.H.apply(phi)  # built for the TRUE metric G
    lhs_wrong = g_inner(A.apply(psi), phi, sn, with_cosine=False)
    rhs_wrong = g_inner(psi, AH_phi, sn, with_cosine=False)
    rel_wrong = abs(lhs_wrong - rhs_wrong) / (
        abs(lhs_wrong) + abs(rhs_wrong) + 1e-300
    )
    if not rel_wrong > 1e-3:
        pytest.fail(
            f"[{case}] wrong metric (drop |Ω·n|) did NOT break reciprocity "
            f"(rel={rel_wrong:.2e}); the |Ω·n|·w_n weighting is NOT load-bearing — "
            f"§2 proves nothing."
        )


# ═══════════════════════════════════════════════════════════════════════
# §4 — G3: the FULL within-group loss ``(L + C - S - B)`` (Phase 2.5 S0.2)
# ═══════════════════════════════════════════════════════════════════════
#
# S† landed at 15185e5 (#276 A2b / #118) and F† with it, so composite ``.H``
# reachability extends to the full within-group loss. These rows harden the
# adjoint-matvec surface the #280 apply-loop unification (2.5a) rebuilds —
# the pre-carve composite canary (gate-spec §13, S0 row).
#
# Config discipline (vv §0.6): the group-transfer transpose lives in S, so
# the mixtures carry ASYMMETRIC SigS (strong downscatter, weak upscatter —
# symmetric SigS ⟹ S† = S and the gate is blind, the ERR-002 family). The
# 2G rows additionally carry asymmetric P1 blocks + Sig2 ≠ 0 (the n2n arm
# rides S's kernel, scattering_order=1); the 4G rows are P0-only (the pure
# group-axis trap at width 4). Heterogeneous 2-material meshes; slab +
# sphere (the curvilinear trace metric is the non-degenerate one).
#
# Mode-12 honesty: a POSING-sign mutation (+S for −S) is INVISIBLE to
# reciprocity BY CONSTRUCTION — ⟨Aψ,φ⟩_G = ⟨ψ,A.Hφ⟩_G holds for whatever A
# was posed, so this gate is NEVER credited against posing errors (those
# are pinned by the solver-level k anchors). The S-specific tooth here is
# the transpose-drop (S wired where Sᵀ belongs inside .H) — O(1) under
# asymmetric SigS.

_P0_2G_A = np.array([[0.38, 0.10], [0.05, 0.90]])
_P1_2G_A = np.array([[0.02, 0.01], [0.00, 0.04]])
_P0_2G_B = np.array([[0.55, 0.03], [0.12, 0.40]])
_P1_2G_B = np.array([[0.06, 0.02], [0.01, 0.03]])

# 4G [g_from, g_to]: strong downscatter above the diagonal, weak upscatter
# below — strongly non-normal group transfer at width 4.
_P0_4G_A = np.array(
    [
        [0.30, 0.15, 0.05, 0.01],
        [0.01, 0.35, 0.20, 0.04],
        [0.00, 0.02, 0.40, 0.25],
        [0.00, 0.00, 0.01, 0.50],
    ]
)
_P0_4G_B = np.array(
    [
        [0.45, 0.08, 0.02, 0.00],
        [0.03, 0.28, 0.30, 0.02],
        [0.01, 0.01, 0.35, 0.18],
        [0.00, 0.02, 0.03, 0.42],
    ]
)


def _mix_2g(p0: np.ndarray, p1: np.ndarray, sig2: np.ndarray):
    m = make_mixture(
        sig_t=np.array([0.5, 1.0]), sig_c=np.array([0.01, 0.02]),
        sig_f=np.array([0.0, 0.0]), nu=np.array([0.0, 0.0]),
        chi=np.zeros(2), sig_s=p0,
    )
    m.SigS = [csr_matrix(p0), csr_matrix(p1)]
    m.Sig2 = csr_matrix(sig2)
    return m


def _mix_4g(p0: np.ndarray):
    ng = 4
    return make_mixture(
        sig_t=np.array([0.6, 0.9, 1.1, 1.4]), sig_c=np.full(ng, 0.02),
        sig_f=np.zeros(ng), nu=np.zeros(ng), chi=np.zeros(ng), sig_s=p0,
    )


def _full_loss_case(coord: CoordSystem, ng: int):
    r"""Het 2-material mesh + the full loss ``A = (L + C) - S - A_BA - B``.

    ``S`` is the PRODUCTION scattering operator off a real :class:`SNSolver`
    (never a re-implementation); σ_t comes off the same material field the
    solver consumes (one source).

    **Campaign step 4c (THE LIFT) — the A_BA reciprocity leaf.** On a
    seed-carrying mesh (the sphere, R12a) the within-group loss carries the
    bulk→ray coupling ``A_BA = RadialCharacteristicEmission`` (``S`` is now pure
    bulk). Its ``.H`` sums the ``w·K_isoᵀ(Reconstructionᵀ χ_seed)`` bulk pullback
    the S-adjoint used to carry inline — WITHOUT this leaf the composite ``.H`` is
    mathematically WRONG on a carrying mesh (a silent-correctness landmine). A
    **present-zero** seed hides the loss (``Reconstructionᵀ(0) = 0``), but
    ``_random_composite`` seeds a NONZERO ψ½ block, so the sphere rows are a real
    composite-level catcher of an ``A_BA.apply_transpose`` corruption
    (``test_tooth_a_ba_transpose_drop_reds``). Seedless meshes (slab) have no ray
    coupling, so the leaf is absent there (``RadialCharacteristicEmission`` rejects
    a seedless mesh).

    NOTE (R2): reciprocity ``⟨Aψ,φ⟩_G = ⟨ψ, A.Hφ⟩_G`` holds for WHATEVER A is
    posed, so it cannot catch a *deleted* A_BA leaf (A.apply and A.H drop it
    consistently) — the DELETED-leaf catcher is L4-S (the driver-routing sentinel
    in ``test_psi_half_coupling``); this leaf makes the gate test the CORRECT
    complete loss and catches an A_BA transpose that is present-but-wrong.
    """
    quad = Quadrature.gauss_legendre(4)
    edges = np.linspace(0.0, 1.0, 5)
    mat_ids = np.array([0, 1, 1, 0])
    if coord is CoordSystem.CARTESIAN:
        mesh = Mesh1D(
            edges=edges, mat_ids=mat_ids, coord=coord,
            bc_left=BC("reflective"), bc_right=BC("reflective"),
        )
    else:
        mesh = Mesh1D(
            edges=edges, mat_ids=mat_ids, coord=coord,
            bc_right=BC("reflective"),
        )
    if ng == 2:
        mixtures = {
            0: _mix_2g(_P0_2G_A, _P1_2G_A, np.array([[0.0, 0.03], [0.01, 0.0]])),
            1: _mix_2g(_P0_2G_B, _P1_2G_B, np.array([[0.0, 0.02], [0.02, 0.0]])),
        }
        order = 1
    else:
        mixtures = {0: _mix_4g(_P0_4G_A), 1: _mix_4g(_P0_4G_B)}
        order = 0
    sn = SNMesh(mesh, quad, mixtures)
    S = SNSolver(sn, scattering_order=order).scattering_op
    sig_t = np.asarray(
        sn.material_xs_field().total_cross_section_field.values, dtype=float
    )
    # B.2d: this is System A's (A,A) block ``L + C − S − B_a`` on the honest
    # 2-block composite — the coupling blocks (A_AB/A_BA/B_b) live on the
    # coupled grid, whose full block ``.H`` reciprocity is G-c2.5
    # (``test_psi_half_coupling::TestCoupledBuilder``); a FullField-summable
    # fused spelling is unrepresentable since the eviction.
    A = (
        StreamingOperator(sn)
        + MultiplicationOperator.from_mesh(sig_t, sn)
    ) - S - SNBoundaryOperator(sn)
    return sn, A, S


def _full_loss_case_cart2d():
    r"""cart2d-DD full loss ``A = (L+C) − S − B`` (#310 C4 — the multi-D row).

    The 2-D sibling of :func:`_full_loss_case`: the SAME het 2-material
    P0/P1 mixtures on a reflective NON-SQUARE box (nx=4 ≠ ny=5), so S
    carries live asymmetric group transfer and the trace metric is the
    4-face ``|Ω·n|·w_n``.  ``scattering_order=0`` — the S† content here is
    the group-transfer transpose composing through ``OperatorSum.H`` on
    the multi-D walk; the anisotropic-order composition rows stay with the
    1-D full-loss cases.
    """
    from orpheus.geometry import Mesh2D

    mesh = Mesh2D(
        edges_x=np.linspace(0.0, 2.0, 5),
        edges_y=np.linspace(0.0, 3.0, 6),
        mat_map=np.array([[0, 1, 1, 0, 0], [1, 0, 0, 1, 0],
                          [0, 0, 1, 0, 1], [1, 0, 0, 0, 1]]),
        bc_xmin=BC("reflective"), bc_xmax=BC("reflective"),
        bc_ymin=BC("reflective"), bc_ymax=BC("reflective"),
    )
    mixtures = {
        0: _mix_2g(_P0_2G_A, _P1_2G_A, np.array([[0.0, 0.03], [0.01, 0.0]])),
        1: _mix_2g(_P0_2G_B, _P1_2G_B, np.array([[0.0, 0.02], [0.02, 0.0]])),
    }
    sn = SNMesh(mesh, Quadrature.level_symmetric(4), mixtures)
    S = SNSolver(sn, scattering_order=0).scattering_op
    sig_t = np.asarray(
        sn.material_xs_field().total_cross_section_field.values, dtype=float
    )
    A = (
        StreamingOperator(sn)
        + MultiplicationOperator.from_mesh(sig_t, sn)
    ) - S - SNBoundaryOperator(sn)
    return sn, A, S


_FULL_LOSS_BUILDERS = {
    "slab_2g": lambda: _full_loss_case(CoordSystem.CARTESIAN, 2),
    "sphere_2g": lambda: _full_loss_case(CoordSystem.SPHERICAL, 2),
    "slab_4g": lambda: _full_loss_case(CoordSystem.CARTESIAN, 4),
    "sphere_4g": lambda: _full_loss_case(CoordSystem.SPHERICAL, 4),
    "cart2d_2g": _full_loss_case_cart2d,
}


def _one_hot_group_composite(
    sn: SNMesh, rng: np.random.Generator, g: int,
) -> TimedFullField:
    r"""A composite supported ONLY in group ``g`` (bulk random, trace zero).

    The S group-transfer transpose lives in the BULK group axis (S is
    trace-silent), so the one-hot bulk isolates the g-restricted residual;
    the trace blocks are group-diagonal and fully covered by the random-φ
    row above.
    """
    N, ng = sn.quad.N, sn.ng
    bulk_vals = np.zeros((N, ng, *sn.spatial_shape))
    bulk_vals[:, g] = rng.standard_normal((N, *sn.spatial_shape))
    boundary = AngularBoundaryFlux(
        values=np.zeros(int(sn.angular_trace.layout.total_size)),
        space=sn.angular_trace,
        mesh=sn,
    )
    return TimedFullField(
        interior=AngularFlux.from_mesh(bulk_vals, sn),
        boundary=boundary,
        _history=(), history_depth=2,
    )


@pytest.mark.foundation
@pytest.mark.parametrize("case", list(_FULL_LOSS_BUILDERS))
def test_full_loss_g_adjoint_reciprocity(case):
    r"""``⟨Aψ,φ⟩_G = ⟨ψ, A.H φ⟩_G`` for the FULL loss ``A = (L+C) - S - B``.

    The composite claim (beyond the leaf-level Euclidean S† gates in
    ``test_scattering_adjoint``): S† COMPOSES into the full-loss G-adjoint
    through ``OperatorSum.H`` WITH the block metric.
    """
    sn, A, _ = _FULL_LOSS_BUILDERS[case]()
    if not A.is_adjointable:
        pytest.fail(
            f"[{case}] full-loss adjoint unreachable — S†/F† landed (#276), "
            f"so the a∧b chain must propagate is_adjointable."
        )
    rng = np.random.default_rng(20260704)
    psi = _random_composite(sn, rng)
    phi = _random_composite(sn, rng)
    lhs = g_inner(A.apply(psi), phi, sn)
    rhs = g_inner(psi, A.H.apply(phi), sn)
    rel = abs(lhs - rhs) / (abs(lhs) + abs(rhs) + 1e-300)
    if not rel < 1e-12:
        pytest.fail(
            f"[{case}] full-loss G-reciprocity broken: ⟨Aψ,φ⟩_G={lhs:.6e} vs "
            f"⟨ψ,A.Hφ⟩_G={rhs:.6e} (rel={rel:.2e})"
        )


@pytest.mark.foundation
@pytest.mark.parametrize("case", list(_FULL_LOSS_BUILDERS))
def test_full_loss_reciprocity_per_group_one_hot(case):
    r"""Per-group rows (vv L27): ``φ_g`` one-hot in each group ``g``.

    A g'→g group-transfer transpose swap shows in the g-restricted residual
    — the weight-summed scalar of the random-φ row could partially mask it
    (anti-pattern #8's group-axis sibling).
    """
    sn, A, _ = _FULL_LOSS_BUILDERS[case]()
    rng = np.random.default_rng(97)
    psi = _random_composite(sn, rng)
    for g in range(sn.ng):
        phi_g = _one_hot_group_composite(sn, rng, g)
        lhs = g_inner(A.apply(psi), phi_g, sn)
        rhs = g_inner(psi, A.H.apply(phi_g), sn)
        rel = abs(lhs - rhs) / (abs(lhs) + abs(rhs) + 1e-300)
        if not rel < 1e-12:
            pytest.fail(
                f"[{case}] per-group reciprocity broken in group {g}: "
                f"lhs={lhs:.6e} rhs={rhs:.6e} (rel={rel:.2e}) — a "
                f"group-transfer transpose defect restricted to g={g}."
            )


@pytest.mark.foundation
def test_tooth_s_transpose_drop_reds(monkeypatch):
    r"""TOOTH: wiring S (not Sᵀ) into the adjoint breaks full-loss
    reciprocity O(1) under asymmetric SigS — the gate SEES the
    group-transfer transpose. In-process monkeypatch; ``-O``-safe."""
    sn, A, S = _FULL_LOSS_BUILDERS["slab_2g"]()
    monkeypatch.setattr(
        type(S), "apply_transpose", lambda self, chi: self.apply(chi),
    )
    rng = np.random.default_rng(11)
    psi = _random_composite(sn, rng)
    phi = _random_composite(sn, rng)
    lhs = g_inner(A.apply(psi), phi, sn)
    rhs = g_inner(psi, A.H.apply(phi), sn)
    rel = abs(lhs - rhs) / (abs(lhs) + abs(rhs) + 1e-300)
    if not rel > 1e-6:
        pytest.fail(
            f"transpose-drop mutation NOT caught (rel={rel:.2e}) — the "
            f"full-loss gate is blind to the S group-transfer transpose; "
            f"check the mixture asymmetry has not degraded."
        )


# ``test_tooth_a_ba_transpose_drop_reds`` RETIRED at B.2d d2 with its
# mechanism: the FullField-summable fused spelling (the FusedRayEmissionGain
# shim leaf inside ``_full_loss_case``) is unrepresentable on the 2-block
# composite. The A_BA-transpose corruption catcher is the BLOCK-level
# reciprocity tooth ``test_psi_half_coupling::TestCoupledLift.
# test_L1_adj_pullback_catcher_has_teeth`` (the direct nonzero-χ Euclidean
# fwd↔adj cross-check on the operator itself).


