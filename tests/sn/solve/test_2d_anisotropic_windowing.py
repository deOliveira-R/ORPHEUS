r"""Phase 5a — angular-windowing correctness guards (the moment-path trap).

The Phase 5a carve reduces the 2-D source-iteration iterate from the full
per-ordinate angular field ``(N, ng, nx, ny)`` to the reduced moment state the
scattering operator :math:`S` consumes, performing the angular→moment reduction
:math:`\phi_\ell^m = \sum_n w_n Y_\ell^m \psi_n` per anti-diagonal instead of
materializing the whole field every sweep.  It is a REPRESENTATION-ONLY carve:
the converged solution MUST be bit-identical (only the live working-set shrinks).

These tests exist to make the central trap UNABLE TO PASS:

* **The trap** — "reduce angular→scalar" naively keeps only :math:`\phi_0`
  (the scalar flux).  For :math:`P_\ell` scattering of order :math:`\ge 1`, the
  iterate needs ALL moments :math:`\phi_\ell^m` (:math:`\ell \le L`); a windowing
  that drops :math:`\ell\ge 1` silently breaks anisotropic 2-D while passing
  every isotropic / 1-group / homogeneous test (vv-principles H1/H2/H3).

The bit-identity reference for the converged scalar flux is the frozen snapshot
``2d_2g_p1_aniso_dd_8x4_het_si`` (regenerated at the pre-carve commit; pinned by
``tests/sn/regression/test_dd_regression.py``).  THIS file carries the live
guards the snapshot cannot:

* (a) the P1 ≠ P0 degeneracy precondition (the moment path carries signal) +
  SI ≡ Krylov cross-check (Krylov stays full-angular → genuine independent
  reference, not self-agreement);
* (b) the FULL angular-flux reconstruction byte-identity (the snapshots store
  scalar flux ONLY, so the reconstructed ``(N, ng, nx, ny)`` field is otherwise
  uncovered — a windowing can produce the right :math:`\phi` while corrupting
  :math:`\psi`);
* (c) the reflective-trace-non-zero precondition (windowing is interior-only;
  the boundary trace must survive — guards the latent-dud where windowing
  zeroes the boundary and the interior alone still matches).

See ``.claude/agent-memory/test-architect/``
``phase5a_angular_windowing_verification_plan.md`` for the full plan.

NOTE on the verification posture: these tests PIN the current (pre-carve)
behaviour.  Authored against HEAD, they pass today; they FAIL if the carve
perturbs the converged value, the moment content, the full reconstruction, or
the boundary trace.  The bit-identity tolerance is exact (``np.array_equal``)
unless the implementer chooses a per-anti-diagonal accumulation that changes the
FP reduction order — in which case the plan's §6 principled-equivalence
relaxation (``assert_array_almost_equal_nulp``) applies and these asserts are
narrowed WITH a documented justification, never silently loosened.
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.derivations.common.xs_library import get_mixture
from orpheus.geometry import BC
from orpheus.geometry.mesh import Mesh2D
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn import solve_sn_fixed_source
from orpheus.sn.mesh.augmented_mesh import SNMesh
from orpheus.sn.coupled_system import build_within_group_system
from orpheus.sn.solver import (
    SNSolver,
    _select_si_splitting,
)
from orpheus.transport.fields.angular_flux import AngularFlux
from orpheus.transport.fields.angular_boundary_flux import AngularBoundaryFlux
from orpheus.transport.timed_full_field import TimedFullField

# L1 — equation-level (anisotropic scattering source + 2-D streaming closure).
pytestmark = pytest.mark.l1


# ── The shared config (mirror of the snapshot ``2d_2g_p1_aniso_dd_8x4_het_si``) ──
#
# 2-D Cartesian, 2G, fuel|moderator split across x (genuinely NON-FLAT flux),
# vacuum-x / reflective-y, mixture B (μ̄=0.6, strongly anisotropic, P1 data),
# uniform isotropic source, source-iteration inner.  Single source of truth so
# the live tests and the snapshot cannot drift in configuration.


def _build_config():
    """Return ``(materials, mesh, quadrature, q_ext, kwargs)`` for the case."""
    fuel = get_mixture("A", "2g")     # fissile, but eigen unused (fixed-source)
    mod = get_mixture("B", "2g")      # μ̄=0.6 strongly anisotropic, P1 data
    nx, ny = 8, 4
    mat = np.zeros((nx, ny), dtype=int)
    mat[:4, :] = 2                    # fuel (id 2) | moderator (id 0) split in x
    mesh = Mesh2D(
        edges_x=np.linspace(0.0, 2.0, nx + 1),
        edges_y=np.linspace(0.0, 1.0, ny + 1),
        mat_map=mat,
        bc_xmin=BC("vacuum"), bc_xmax=BC("vacuum"),
        bc_ymin=BC("reflective"), bc_ymax=BC("reflective"),
    )
    quad = Quadrature.level_symmetric(sn_order=4)
    sum_w = float(quad.weights.sum())
    # Per-ordinate density (producer-side /W projection, R-1 Step 4 A1).
    q_ext = np.full((quad.N, 2, nx, ny), 1.0 / sum_w)
    kwargs = dict(
        scattering_order=1, inner_solver="source_iteration",
        max_inner=3000, inner_tol=1e-12,
    )
    return {2: fuel, 0: mod}, mesh, quad, q_ext, kwargs


def _solve(scattering_order: int | None = None, inner_solver: str | None = None):
    """Solve the shared config, optionally overriding order / inner solver."""
    materials, mesh, quad, q_ext, kwargs = _build_config()
    if scattering_order is not None:
        kwargs["scattering_order"] = scattering_order
    if inner_solver is not None:
        kwargs["inner_solver"] = inner_solver
    return solve_sn_fixed_source(materials, mesh, quad, q_ext, **kwargs)


# ─────────────────────────────────────────────────────────────────────────
# (a) Anisotropic scattering — the moment path carries signal, AND SI ≡ Krylov.
# ─────────────────────────────────────────────────────────────────────────


@pytest.mark.verifies("pn-scatter", "transport-cartesian", "multigroup")
def test_2d_p1_aniso_moment_path_carries_signal_and_si_krylov_agree():
    r"""P1 ≠ P0 (moments :math:`\ell\ge1` are load-bearing) AND SI ≡ Krylov.

    The DEGENERACY GUARD for the angular-windowing carve.  Two assertions:

    1. **P1 ≠ P0 precondition.**  If the converged scalar flux were
       insensitive to the :math:`\ell\ge1` moments on this config, the
       snapshot would pin nothing about the moment path — the "windowing
       drops :math:`\phi_\ell`" trap would slip through a degenerate gate.
       Here P1 differs from P0 by ~18 % → the moments genuinely matter.

    2. **SI ≡ Krylov cross-check.**  Krylov stays FULL-ANGULAR (the carve
       windows only the SI iterate), so SI-windowed ≡ Krylov-full is a
       genuine independent reference, NOT twin self-agreement.  Post-carve,
       a windowing that drops :math:`\phi_\ell` makes the SI flux disagree
       with the unwindowed Krylov flux on this anisotropic config.
    """
    sol_p0 = _solve(scattering_order=0)
    sol_p1 = _solve(scattering_order=1)
    phi_p0 = np.asarray(sol_p0.scalar_flux.values, dtype=np.float64)
    phi_p1 = np.asarray(sol_p1.scalar_flux.values, dtype=np.float64)

    # (1) Precondition: the moment path carries signal.
    rel_p1_vs_p0 = np.max(np.abs(phi_p1 - phi_p0)) / np.max(np.abs(phi_p0))
    assert rel_p1_vs_p0 > 1e-2, (
        f"P1 vs P0 rel diff {rel_p1_vs_p0:.2e} too small — the ℓ≥1 moments "
        "carry no signal on this config, so the moment-path gate is "
        "DEGENERATE (cannot catch a φ_ℓ-dropping windowing). Reconfigure."
    )

    # (1b) Precondition: the flux is genuinely NON-FLAT (redistribution active).
    prof = phi_p0[0].mean(axis=1)  # group-0 x-profile
    assert prof.max() / prof.min() > 1.2, (
        f"flux too flat (max/min={prof.max() / prof.min():.3f}); the "
        "2-D streaming redistribution is not exercised — gate degenerate."
    )

    # (2) Cross-check: the windowed SI flux must match the full-angular Krylov.
    sol_kry = _solve(scattering_order=1, inner_solver="krylov")
    phi_kry = np.asarray(sol_kry.scalar_flux.values, dtype=np.float64)
    np.testing.assert_allclose(
        phi_p1, phi_kry, rtol=1e-6, atol=1e-8,
        err_msg=(
            "2-D P1 anisotropic SI flux disagrees with the FULL-ANGULAR "
            "Krylov flux beyond iterative tolerance — a windowing that "
            "dropped the ℓ≥1 moments would land here."
        ),
    )


# ─────────────────────────────────────────────────────────────────────────
# (b) Final full reconstruction — Solution.angular_flux byte-identity.
# ─────────────────────────────────────────────────────────────────────────


@pytest.mark.verifies("pn-scatter", "transport-cartesian")
def test_2d_windowed_si_full_angular_flux_self_consistent():
    r"""The reconstructed FULL angular flux reproduces the scalar flux exactly.

    The snapshots store the scalar flux ONLY, so the reconstructed
    ``(N, ng, nx, ny)`` angular field is otherwise UNVERIFIED.  A windowing
    that discards interior :math:`\psi` between anti-diagonals could converge
    the right :math:`\phi` while producing a wrong reconstructed :math:`\psi`
    (the inverse of vv-principles H3 — :math:`\phi` right, :math:`\psi` wrong).

    This pins the reconstruction by the self-consistency identity
    :math:`\phi = \sum_n w_n \psi_n` (``integrate_angular``): the full
    angular flux returned by ``Solution.angular_flux`` MUST reduce to the
    returned scalar flux to machine precision.  Post-carve, if the final
    reconstruction is not a full sweep, this identity breaks.

    NOTE: this is a NECESSARY structural check, not the bit-identity gate.
    The bit-identity of the full field is the step-3 pin in the plan ladder
    (capture pre-carve ``angular_flux.interior.values``, assert ``np.array_equal``
    post-carve).  That pin lives in the carve session (it needs the pre-carve
    array on disk); this self-consistency test is the always-on guard that
    needs no frozen reference.
    """
    sol = _solve(scattering_order=1)
    psi = np.asarray(sol.angular_flux.interior.values, dtype=np.float64)  # (N,ng,nx,ny)
    phi = np.asarray(sol.scalar_flux.values, dtype=np.float64)         # (ng,nx,ny)
    weights = sol.mesh.quad.weights

    assert psi.shape[0] == weights.shape[0], (
        f"angular flux N={psi.shape[0]} disagrees with quadrature "
        f"N={weights.shape[0]} — the full reconstruction was truncated."
    )
    phi_from_psi = np.einsum("n,ngij->gij", weights, psi)
    np.testing.assert_allclose(
        phi_from_psi, phi, rtol=1e-10, atol=1e-12,
        err_msg=(
            "Σ_n w_n ψ_n (from the reconstructed full angular flux) does NOT "
            "equal the returned scalar flux — the windowed SI corrupted the "
            "final angular-flux reconstruction (φ right, ψ wrong)."
        ),
    )


# ─────────────────────────────────────────────────────────────────────────
# (c) Reflective BC — windowing is interior-only; the boundary trace survives.
# ─────────────────────────────────────────────────────────────────────────


@pytest.mark.verifies("transport-cartesian")
def test_2d_windowed_si_reflective_trace_is_nonzero():
    r"""The reflective-y boundary trace is non-zero post-solve.

    Windowing is INTERIOR-bulk-only; the reflective coupling reads the full
    angular flux at the boundary FACES via the typed
    :math:`\iota_*`/:math:`\iota^*` trace exchange, which the carve must NOT
    touch.  This guards the latent-dud where a windowing accidentally zeroes
    the boundary trace and the test still passes because the interior alone
    matches the reference (the boundary contribution is small here).

    The precondition: on a reflective-y config the converged boundary flux
    MUST carry non-zero outflow on the reflective faces.  A zeroed trace
    means the reflective coupling was dropped — a correctness bug the
    interior-only scalar-flux snapshot would miss.
    """
    sol = _solve(scattering_order=1)
    boundary = sol.angular_flux.boundary
    # Sum |trace| over all faces; on a reflective-y problem with a non-flat
    # interior this is strictly positive.
    total_trace = 0.0
    for face in boundary.layout.faces:
        total_trace += float(np.sum(np.abs(boundary.face_view(face))))
    assert total_trace > 1e-6, (
        f"converged boundary trace ‖trace‖₁={total_trace:.2e} is ~0 — the "
        "reflective coupling was dropped (windowing must be interior-only)."
    )


# ─────────────────────────────────────────────────────────────────────────
# (d) Phase 5c — the in-sweep moment accumulation ≡ the fuller-view oracle.
# ─────────────────────────────────────────────────────────────────────────


def _windowed_product_and_oracle_operands(
    scattering_order: int | None = None, *, inner_schedule: str = "jacobi",
):
    """Shared operands for the deforestation gates + the mutation twin:
    ``(quad, base, product, rhs, oracle_fn)`` on the production config.

    ``inner_schedule`` selects the base forward through the production
    dispatch (:func:`_select_si_splitting`): ``"jacobi"`` → the plain
    ``(L+C)``; ``"gauss_seidel"`` → the reified splitting
    ``M = (L+C) − B_lower`` (#226 §17 W2) — the windowed×G-S corner the
    step-3 pin exercises.
    """
    materials, mesh, quad, _q_ext, kwargs = _build_config()
    L = (
        kwargs["scattering_order"]  # P1 — ℓ≥1 is load-bearing
        if scattering_order is None
        else scattering_order
    )
    solver = SNSolver(SNMesh(mesh, quad, materials), scattering_order=L)

    # The within-group forward + the scattering operator — the SAME
    # operators (and the SAME schedule dispatch) the windowed SI driver
    # consumes (single source of truth).
    system = build_within_group_system(
        solver.sn_mesh, solver.mat_xs, scattering_op=solver.scattering_op,
    )
    LC, (S, N2N, B) = system.implicit_operator, system.explicit_gains  # seedless 2-D record shape (§14.1: S, N2N, B_a)
    sn_mesh = solver.sn_mesh
    base, _boundary_gain = _select_si_splitting(LC, B, sn_mesh, inner_schedule)

    # THE production windowed object (#226 steps 2–3, §17 W1): the typed
    # composition ``P @ A.inverse()`` — the scattering operator's MINTED
    # flux-analysis face on the bulk ⊕ identity on the trace, composed with
    # the sweep inverse, FUSED through the substrate's moment emit.
    # ``S.flux_analysis`` is the SAME minted object the driver's factory
    # (``_maybe_window``) injects, so the SUT exercises the production
    # projection, not a test-local one.
    from orpheus.sn.operators.windowing import BulkAnalysisOperator

    product = BulkAnalysisOperator(S.flux_analysis, sn_mesh.full_field_space) @ base.inverse()

    # A representative per-ordinate source (seeded random ⇒ strong, deterministic
    # ℓ≥1 content in the swept ψ; the projection-order equivalence is
    # input-independent — it is a property of the reduction tree, not of ψ).
    rng = np.random.default_rng(50301)
    rhs = TimedFullField.zeros(
        interior=AngularFlux, boundary=AngularBoundaryFlux, space=sn_mesh.full_field_space,
    )
    rhs.interior.values[...] = rng.uniform(0.0, 1.0, size=rhs.interior.values.shape)

    def oracle_fn():
        # ORACLE (the fuller view, DEFORESTED): flat post-projection (the
        # frame's analysis face) of the full-angular sweep — structurally the
        # inherited ``OperatorProduct.apply`` body, spelled independently.
        return np.asarray(
            S.frame.analysis.apply(base.solve(rhs).interior.values),
            dtype=np.float64,
        )

    return quad, base, product, rhs, oracle_fn


@pytest.mark.verifies(
    "pn-scatter", "harmonic-moment-projection", "transport-cartesian",
)
def test_2d_windowed_product_equals_post_projection():
    r"""W1 deforestation gate (#226 step 2, spec §13.4): the FUSED windowed
    product ``(P @ A.inverse()).apply(rhs)`` (per-anti-diagonal in-sweep
    moment accumulation) ≡ the flat post-sweep ``frame.analysis`` projection
    of the SAME full-angular sweep — the FULLER-VIEW VERIFICATION ORACLE.

    Renamed/migrated from ``test_2d_windowed_moments_in_sweep_equal_post_
    projection`` when the ``solve_moments`` cross-reach retired into the
    typed product (§17 W1): the SUT is now the production
    :class:`~orpheus.sn.operators.windowing.WindowedSweep` object itself.
    Fusion (never materialize the ``(N, ng, nx, ny)`` intermediate) is an
    EVALUATION STRATEGY of the one composition; the deforested oracle is
    structurally the inherited ``OperatorProduct.apply`` body.  Both share
    the SAME cell kernel and the SAME ``frame`` (table + weights), so only
    the ordinate-sum reduction ORDER differs ⇒ ULP-level drift
    (principled-equivalence, ``vv-principles``).

    This is a REPRESENTATION-equivalence pin, NOT a correctness claim — it is
    procedurally (NOT structurally) independent from the oracle (shared kernel).
    Its structural-independence anchor is
    :func:`test_2d_p1_aniso_moment_path_carries_signal_and_si_krylov_agree`
    (the windowed-SI moment ℓ=0 ≡ the FULL-ANGULAR Krylov scalar) + the
    closed-form ``k_inf`` eigenvalue gate; the coisometry factor of the
    honest contract is the LANDED
    ``test_frame.py::test_pi_R_is_4pi_identity_through_the_frame``
    (``analysis ∘ reconstruction = 4π·I`` — cited, NOT re-minted; asserting
    ``= I`` would be ERR-051).  Its value-add: it pins the ℓ≥1
    (anisotropic) moment block that the ℓ=0 scalar cross-check is BLIND to —
    the load-bearing catcher for an ℓ/m index drift (vv-principles Mode 5) or a
    wrong octant-``Y`` slice (Mode 2) in the per-level
    ``"nlm,ngd,n->lmgd"`` accumulation, which blow past 32 ULP.
    M-DEFOREST teeth: :func:`test_mutation_dropped_higher_moments_redden`.

    leg-2 (ℓ≥1 non-vacuity) guards against a degenerate gate: the ℓ≥1 block
    being pinned MUST carry signal (else the comparison is vacuous there) —
    the Mode-7 witness that the anisotropic P1 config is load-bearing.

    Metric (vv-principles §"Bit-identity vs principled-equivalence", criterion
    3 "drift dimensionally explainable"): the max drift RELATIVE TO THE
    MOMENT-TENSOR SCALE, bounded by ``4·N·eps`` (``N`` = reduction depth = the
    ordinate count; 4× headroom for the cross-octant ``+=`` nesting).  NOT
    element-wise ``assert_array_almost_equal_nulp``: the moment tensor spans ~3
    orders of magnitude (the ℓ=0 scalar dominates; ℓ≥1 are ~10⁻³×), so
    element-wise ULP inflates a machine-eps ABSOLUTE diff on a small-magnitude
    ℓ≥1 element into hundreds of ULP even when the reorder is pure FP noise —
    the scale-relative drift (de-risk measured ``2.7e-16``, ~78× under the
    ``4·N·eps ≈ 2.1e-14`` bound; re-measured ``1.8e-16`` on the product SUT)
    is the principled-equivalence quantity.
    """
    quad, _base, product, rhs, oracle_fn = (
        _windowed_product_and_oracle_operands()
    )

    # SUT: the FUSED typed product (substrate per-anti-diagonal moment emit).
    out = product.apply(rhs)
    # The runtime half of the codomain fact (§13.4 factor 3): the composite's
    # bulk IS the moment flux (the static ``FullField → moment-bulk``
    # refinement is not expressible until the composite is bulk-generic —
    # see ``_windowed_product_static_typing_pins``).
    from orpheus.transport.fields.harmonic_moment_flux import (
        HarmonicMomentFlux,
    )

    if not isinstance(out.interior, HarmonicMomentFlux):
        raise AssertionError(
            f"windowed product bulk is {type(out.interior).__name__}, "
            f"expected HarmonicMomentFlux — the codomain fact broke"
        )
    sut = np.asarray(out.interior.values, dtype=np.float64)
    oracle = oracle_fn()

    # leg-2: ℓ≥1 carries signal (anti-degeneracy) — raise fires under -O.
    l0 = float(np.abs(oracle[0]).max())
    l_ge1 = float(np.abs(oracle[1:]).max())
    if not (l_ge1 > 1e-3 * l0):
        raise AssertionError(
            f"ℓ≥1 moment block is vacuous (max|ℓ≥1|={l_ge1:.2e} ≤ "
            f"1e-3·max|ℓ0|={1e-3 * l0:.2e}); the nulp pin on the higher "
            "moments would be degenerate — reconfigure the source/config."
        )

    # leg-1: representation equivalence — max drift RELATIVE to the moment-tensor
    # scale ≤ 4·N·eps (the principled-equivalence metric; element-wise nulp would
    # inflate on the small-magnitude ℓ≥1 block — see the docstring).
    scale = float(np.abs(oracle).max())
    rel_drift = float(np.abs(sut - oracle).max() / scale)
    bound = 4 * quad.N * float(np.finfo(np.float64).eps)
    if not (rel_drift <= bound):
        raise AssertionError(
            f"in-sweep moment accumulation drifts {rel_drift:.2e} (relative to "
            f"scale {scale:.3e}) > 4·N·eps = {bound:.2e} — NOT a pure FP reorder; "
            "a wrong octant-Y slice / missing weight / ℓ-m index drift "
            "(vv-principles Mode 2/3/5) in the per-level "
            '"nlm,ngd,n->lmgd" accumulation.'
        )


@pytest.mark.verifies(
    "pn-scatter", "harmonic-moment-projection", "transport-cartesian",
)
def test_2d_windowed_product_over_gauss_seidel_M_equals_post_projection():
    r"""The step-3 corner pin: the windowed product OVER THE REIFIED G-S
    SPLITTING — ``P @ M.inverse()`` with ``M = (L+C) − B_lower`` — FUSED
    (the SCHEDULED walk's per-anti-diagonal moment emit) ≡ the flat
    post-projection of M's own scheduled solve, within the same
    scale-relative ``4·N·eps``.

    #226 step 2 landed the two compositions separately (W1 pinned
    ``P @ (L+C).inverse()``; W2 pinned M's round-trip, full-angular);
    PRODUCTION composes them — the 2-D Cartesian default
    ``inner_schedule`` IS ``gauss_seidel``, so the driver's actual step
    operator is ``P @ M.inverse()`` — and until this pin the windowed×G-S
    corner was covered only at the ℓ=0/integration level (the end-to-end
    SI value gates).  This is the ℓ≥1 fused-vs-deforested pin ON the
    scheduled walk: the fused leg routes ``moment_frame`` AND
    ``(schedule, reflect)`` through the ONE ``_solve_timed_full_field``
    surface SIMULTANEOUSLY — the interaction no per-feature gate sees.

    leg-3 guards non-degeneracy: the config's reflective-y faces must
    yield a genuinely populated ``B_lower`` (else the split collapses,
    ``M ≡ L+C``, and this pin silently duplicates the Jacobi one).
    """
    from orpheus.sn.operators.scheduled_invertible import (
        ScheduledInvertibleOperator,
    )

    quad, M, product, rhs, oracle_fn = _windowed_product_and_oracle_operands(
        inner_schedule="gauss_seidel",
    )

    # leg-3: the production dispatch selected the reified M, and its split
    # is non-degenerate on this reflective-y config.
    if not isinstance(M, ScheduledInvertibleOperator):
        raise AssertionError(
            f"gauss_seidel dispatch returned {type(M).__name__}, not the "
            f"reified ScheduledInvertibleOperator — the G-S corner is not "
            f"being exercised"
        )
    n_lower = sum(int(rows.size) for rows in M.lower.rows.values())
    if n_lower == 0:
        raise AssertionError(
            "B_lower has no populated rows on the reflective-y config — "
            "the split degenerated and this pin collapsed to the Jacobi one"
        )

    out = product.apply(rhs)
    sut = np.asarray(out.interior.values, dtype=np.float64)
    oracle = oracle_fn()

    # leg-2: ℓ≥1 carries signal (anti-degeneracy).
    l0 = float(np.abs(oracle[0]).max())
    l_ge1 = float(np.abs(oracle[1:]).max())
    if not (l_ge1 > 1e-3 * l0):
        raise AssertionError(
            f"ℓ≥1 moment block is vacuous over M (max|ℓ≥1|={l_ge1:.2e} ≤ "
            f"1e-3·max|ℓ0|={1e-3 * l0:.2e}) — reconfigure the source/config."
        )

    # leg-1: fused ≡ deforested over the SCHEDULED walk (same metric and
    # bound as the Jacobi gate — only the reduction order may differ).
    scale = float(np.abs(oracle).max())
    rel_drift = float(np.abs(sut - oracle).max() / scale)
    bound = 4 * quad.N * float(np.finfo(np.float64).eps)
    if not (rel_drift <= bound):
        raise AssertionError(
            f"windowed×G-S: in-sweep moment accumulation drifts "
            f"{rel_drift:.2e} (relative to scale {scale:.3e}) > 4·N·eps = "
            f"{bound:.2e} — the (moment_frame × schedule/reflect) "
            "interaction through _solve_timed_full_field is NOT a pure FP "
            "reorder."
        )


def test_mutation_dropped_higher_moments_redden(monkeypatch):
    r"""M-DEFOREST teeth (spec §13.4/§13.5): a fused product that drops the
    ℓ≥1 accumulation (emits ℓ=0 only) REDs the deforestation gate on the P1
    config — and the SAME mutation stays green on a P0 config (the Mode-7
    witness that the anisotropic config is load-bearing: an isotropic gate
    is structurally BLIND to a dropped ℓ≥1 block).
    """
    from orpheus.sn.operators.windowing import WindowedSweep

    real = WindowedSweep.apply

    def truncated(self, rhs):
        out = real(self, rhs)
        out.interior.values[1:] = 0.0  # the trap: "reduce angular→scalar"
        return out

    monkeypatch.setattr(WindowedSweep, "apply", truncated)

    # P1: the gate's comparison MUST red.
    quad, _base, product, rhs, oracle_fn = (
        _windowed_product_and_oracle_operands()
    )
    sut = np.asarray(product.apply(rhs).interior.values, dtype=np.float64)
    oracle = oracle_fn()
    bound = 4 * quad.N * float(np.finfo(np.float64).eps)
    rel_drift = float(np.abs(sut - oracle).max() / np.abs(oracle).max())
    if rel_drift <= bound:
        pytest.fail(
            f"M-DEFOREST not caught on P1: drift {rel_drift:.2e} ≤ bound "
            f"{bound:.2e} — the deforestation gate has no ℓ≥1 teeth"
        )

    # P0 control: the same mutation is a no-op (no ℓ≥1 block exists) — the
    # gate stays green, witnessing that P1 anisotropy is what carries teeth.
    quad0, _base0, product0, rhs0, oracle_fn0 = (
        _windowed_product_and_oracle_operands(scattering_order=0)
    )
    sut0 = np.asarray(product0.apply(rhs0).interior.values, dtype=np.float64)
    oracle0 = oracle_fn0()
    rel_drift0 = float(np.abs(sut0 - oracle0).max() / np.abs(oracle0).max())
    if rel_drift0 > 4 * quad0.N * float(np.finfo(np.float64).eps):
        pytest.fail(
            f"P0 control unexpectedly red ({rel_drift0:.2e}) — the Mode-7 "
            "witness is broken (the mutation should be a no-op at L=0)"
        )
