r"""Wave O (#208) O.4b Phase E3 — the 2-D BC-extraction verification gates.

The 2-D Cartesian BC extraction (Phase E1a/E1b/E2, HEAD 02f8f0a) is COMPLETE:
the representation ``loss_action`` matvec (loss_representation.py; since S6.3
the 2-D walk lives on the loss representation, off the operator —
``ScanMarch`` default since S6.9, ``MovingFrontierWindow`` peer, BOTH bare)
and the scheduled sweep are **bare** — they read
``psi.boundary.inflow`` as the GIVEN incoming edge and emit the active
boundary-block residual (OUTFLOW ordinate slots: ``streamed − psi.outflow``;
INFLOW ordinate slots: identity ``psi.inflow``).  The reflective coupling is
the sibling ``−B`` (:class:`~orpheus.sn.boundary_operator.SNBoundaryOperator`),
folded into the 2-D Krylov composed matvec.  The 2-D eigenvalue path is Krylov
(``solve_sn(..., inner_solver="krylov")``); 2-D source-iteration AND 2-D
fixed-source Krylov are both ``NotImplementedError`` (deferred to a later phase).

This file is the E3 verification pass — the gates that PIN the now-complete
extraction.  It implements Gates **V-2D, Q-2D, R-2D, R-2D-NEG** from the
verification spec (``.claude/agent-memory/test-architect/issue_208_o4b_2d_verification_plan.md``
§2, §3, §6).  Gate K (the k_inf catcher) ALREADY exists and is green at
``tests/sn/sweep/cartesian_2d/test_2d_l2_matvec_correctness.py::test_solve_sn_2d_krylov_homogeneous_reflective_recovers_kinf``;
Gates M-2D / D-2D / S-2D are not E3-blocking and live elsewhere.

The pillar / claim-layer / failure-mode for each gate is stated in its
docstring (vv-principles discipline).

TYPE-AGNOSTIC by mandate
========================

The matvec currently emits its boundary residual typed as ``BoundaryFlux``
(the honest retype to ``BoundaryResidual`` is a SEPARATE holistic carve,
B.5.2, deferred).  These gates read the residual ONLY via
``.boundary.face_view(face)`` / ``.boundary.values`` and
``trace.inflow_indices_for_face`` / ``outflow_indices_for_face`` — they do NOT
``isinstance``-assert the boundary type, so they survive the future retype.

Baseline mechanism (Gate V-2D)
==============================

Reuses the O.4a.2 ``--capture-baseline`` flag (``pytest_addoption`` in
``tests/conftest.py``).  With the flag present the V-2D test WRITES the bare
2-D vacuum bulk snapshot under ``tests/sn/_data/bc_extraction_2d_baseline/``
and skips the assert; absent (the default, incl. CI) it READS the committed
snapshot and asserts byte-identity.  Capture command::

    PYTHONPATH=$PWD .venv/bin/python -m pytest \
        tests/sn/operators/test_bc_extraction_2d.py \
        -k vacuum_bulk_bit_identical --capture-baseline

The captured ``.npy`` files are committed alongside this test as the permanent
regression reference (the bare path's vacuum bulk MUST never move a bit).
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.derivations.common.xs_library import get_mixture
from orpheus.geometry import BC, Mesh2D
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.boundary_operator import SNBoundaryOperator
from orpheus.sn.geometry import SNMesh
from orpheus.sn.operator import CollisionOperator, StreamingOperator
from orpheus.transport.fields.angular_flux import AngularFlux
from orpheus.transport.fields.boundary_flux import BoundaryFlux
from orpheus.transport.source_sinks import BoundarySourceSink
from orpheus.transport.timed_full_field import TimedFullField
from tests.sn._test_helpers import SN_TESTS_ROOT

# ─────────────────────────────────────────────────────────────────────
# Baseline snapshot store + the --capture-baseline flag (Gate V-2D).
# ─────────────────────────────────────────────────────────────────────

_BASELINE_DIR = SN_TESTS_ROOT / "_data" / "bc_extraction_2d_baseline"


def _capturing(request) -> bool:
    return bool(request.config.getoption("--capture-baseline", default=False))


# ─────────────────────────────────────────────────────────────────────
# Mesh builders (shared with the E0 diagnostic + Gate K's
# test_2d_l2_matvec_correctness.py).
# ─────────────────────────────────────────────────────────────────────


def _homogeneous_reflective_2d(nx: int = 4, ny: int = 4) -> SNMesh:
    r"""All-reflective 2-D Cartesian mesh; mixture A 2g homogeneous.

    The canonical k_inf / reflective fixture: all-reflective + homogeneous ⇒
    k = νΣ_f/Σ_a = k_inf (1.875 for mixture A 2g) for ALL meshes, and the
    spatially-flat fundamental mode is annihilated by ``(L − B)``.
    """
    geom = Mesh2D(
        edges_x=np.linspace(0.0, 2.0, nx + 1),
        edges_y=np.linspace(0.0, 2.0, ny + 1),
        mat_map=np.zeros((nx, ny), dtype=int),
        bc_xmin=BC("reflective"), bc_xmax=BC("reflective"),
        bc_ymin=BC("reflective"), bc_ymax=BC("reflective"),
    )
    quad = Quadrature.level_symmetric(sn_order=4)
    return SNMesh(geom, quad, {0: get_mixture("A", "2g")})


def _vacuum_2d(nx: int = 4, ny: int = 4) -> SNMesh:
    r"""Vacuum-on-all-faces 2-D Cartesian mesh; mixture A 2g homogeneous."""
    geom = Mesh2D(
        edges_x=np.linspace(0.0, 1.0, nx + 1),
        edges_y=np.linspace(0.0, 1.0, ny + 1),
        mat_map=np.zeros((nx, ny), dtype=int),
        bc_xmin=BC("vacuum"), bc_xmax=BC("vacuum"),
        bc_ymin=BC("vacuum"), bc_ymax=BC("vacuum"),
    )
    quad = Quadrature.level_symmetric(sn_order=4)
    return SNMesh(geom, quad, {0: get_mixture("A", "2g")})


def _sigt_2g(sn_mesh: SNMesh) -> np.ndarray:
    """Per-cell per-group total cross-section field ``(ng, nx, ny)`` for mix A 2g."""
    mix = get_mixture("A", "2g")
    ng = mix.SigT.size
    nx, ny = sn_mesh.spatial_shape
    return np.broadcast_to(mix.SigT[:, None, None], (ng, nx, ny)).copy()


# ═════════════════════════════════════════════════════════════════════
# GATE V-2D — vacuum bit-identity (foundation).
# ═════════════════════════════════════════════════════════════════════
#
# Pillar:  foundation (committed bit-exact snapshot of the bare 2-D matvec).
# Layer:   software-invariant (the vacuum bulk stencil is touched by NOTHING
#          downstream — it must stay byte-frozen forever).
# Catches: ANY future change to the 2-D matvec that perturbs the vacuum bulk
#          path — the carve's own keystone risk (the abandoned active-trace
#          formulation perturbed the bulk and converged to a non-uniform mode
#          ~10% off k_inf; this pin is the permanent tripwire for that class).


@pytest.mark.foundation
class TestVacuum2DBitIdentity:
    """Vacuum 2-D ``L.apply`` bulk is byte-identical to the committed snapshot.

    This gate PINS the bare 2-D vacuum bulk matvec as a permanent strict
    (``np.array_equal``) regression reference.  The committed baseline is the
    **post-#240-D5a frozen value**: E0-T1 proved the bare 2-D vacuum bulk was
    bit-identical to the pre-carve production path, and #240-D5a then
    re-baselined it by ~1 ULP when the 2-D Cartesian matvec moved onto the
    scheme coefficient model (the ``÷D → ×inverse_denom`` FP re-association — a
    sanctioned principled re-baseline per ``vv-principles`` §"Bit-identity vs
    principled-equivalence", NOT a bug).  Henceforth the bulk must not move a
    single bit vs that frozen post-D5a snapshot.  Parametrized over ≥3 fixed
    seeds — a random ψ̄ stresses the diamond-difference bulk stencil (a flat ψ̄
    would NULL the redistribution, vv §H2).
    """

    @pytest.mark.parametrize("seed", [0, 1, 2])
    def test_vacuum_bulk_bit_identical(self, seed, request):
        """[foundation] Bare 2-D vacuum pure-L bulk matvec == committed snapshot, bit-exact.

        Build a vacuum 2-D mesh, seed ``psi.bulk`` with fixed-seed random
        values and a ZERO inflow trace, apply the bare ``StreamingOperator``,
        and assert ``L.apply(psi).bulk.values`` is ``np.array_equal`` to the
        committed baseline.  The baseline was **re-captured at #257 S8b**:
        ``L.apply`` is now pure σ-free ``streaming_action(ψ) = loss_action(0,
        ψ)`` (the ``(L+C)−C`` fold is retired), which re-associates the FP
        reduction tree vs the old subtractive form.  The new frozen value is
        the genuine pure streaming — verified ``== (L+C).apply.bulk − σ_t⊙ψ``
        to ≤64 ULP (the affine relation) with the BYTE-IDENTICAL ``(L+C)``
        composite as the structural ground (NOT old-vs-new proximity).  The
        bit-exact ``array_equal`` gate is RESTORED against this re-baselined
        reference: the pure-L bulk must not move a single bit vs the frozen
        value for any future refactor.
        """
        sn_mesh = _vacuum_2d(nx=4, ny=4)
        sig_t = _sigt_2g(sn_mesh)
        L = StreamingOperator(sn_mesh)

        rng = np.random.default_rng(20260603 + seed)
        state = TimedFullField.zeros(
            bulk=AngularFlux, boundary=BoundaryFlux, mesh=sn_mesh,
        )
        state.bulk.values[...] = rng.standard_normal(state.bulk.values.shape)
        # vacuum: no incoming inflow trace (boundary stays zero).

        out = L.apply(state)

        # B.5.2: the 2-D matvec output boundary is the source/sink role leaf
        # (Aψ, NOT a residual — a residual only arises from from_balance(Aψ, b)).
        # Mirrors the bulk's AngularSourceSink; completes the boundary role grid.
        assert isinstance(out.boundary, BoundarySourceSink)

        key = f"vacuum_bulk_2d_seed{seed}"
        path = _BASELINE_DIR / f"{key}.npy"
        if _capturing(request):
            # Gate the capture on the STRUCTURALLY-INDEPENDENT ground so a
            # re-capture can never launder a streaming bug into the snapshot:
            # pure-L must equal the byte-identical (L+C) composite minus the
            # collision diagonal (the affine relation (L+C)ψ = streaming(ψ) +
            # σ_t⊙ψ; #257 S8b).  A laundered sign/factor error in
            # streaming_action is O(1)-relative and trips this assert long
            # before the bytes are frozen — the re-baseline is NOT
            # self-referential ("freeze whatever pure-L emits") by construction.
            composite = (L + CollisionOperator(sn_mesh, sig_t)).apply(state)
            structural_ground = (
                composite.bulk.values - sig_t[None] * state.bulk.values
            )
            np.testing.assert_allclose(
                out.bulk.values, structural_ground, rtol=1e-11, atol=0.0,
                err_msg=(
                    "capture-time structural guard: pure-L streaming_action does "
                    "NOT match (L+C) − σ_t⊙ψ — refusing to freeze a possibly-"
                    "laundered snapshot."
                ),
            )
            _BASELINE_DIR.mkdir(parents=True, exist_ok=True)
            np.save(path, out.bulk.values)
            pytest.skip(f"captured baseline {key}")
        assert path.exists(), (
            f"missing baseline snapshot {path}; run with --capture-baseline "
            f"to write the bare 2-D vacuum bulk reference, then `git add` it"
        )
        expected = np.load(path)
        np.testing.assert_array_equal(
            out.bulk.values, expected,
            err_msg=(
                f"2-D vacuum seed={seed}: the bare 2-D matvec bulk moved a "
                f"bit vs the committed snapshot — the keystone 2-D path was "
                f"perturbed (the baseline is the post-#240-D5a frozen value; it "
                f"must stay frozen vs that strict reference)."
            ),
        )


# ═════════════════════════════════════════════════════════════════════
# GATE Q-2D — 2-D L0 streaming-equilibrium (per-ordinate flat-flux residual).
# ═════════════════════════════════════════════════════════════════════
#
# Pillar:  closed-form (φ = Q/Σ_t per group; per-ordinate ψ_n = φ/Σw).
# Layer:   flux-shape (the per-ordinate streaming-equilibrium identity).
# Catches: a missing-factor in the extracted boundary coupling (#3); the
#          active-trace non-uniform-mode bug (the ~10% signature) shows as a
#          per-ordinate residual ≠ 0 on flat ψ at the boundary ring.  This is
#          the FORWARD-ONLY matvec-level sentinel that catches the miss WITHOUT
#          a full eigenvalue solve.


@pytest.mark.l0
@pytest.mark.verifies("streaming-equilibrium")
class TestStreamingEquilibrium2D:
    """Uniform Q, uniform Σ_t, fully-reflective 2-D ⇒ flat ψ_n = φ/W is exact.

    Pure-attenuation configuration (no scattering operator in the residual):
    the flat scalar flux ``φ_g = Q_g/Σ_{t,g}`` and the per-ordinate angular
    flux ``ψ_n,g = φ_g/W`` (W = Σw = 4π for level-symmetric).  A consistent
    reflective trace (``ψ.boundary`` uniform = ψ̄ value) makes
    ``L·ψ_flat = 0`` per ordinate, so ``(L_full + C − B)·ψ_flat = C·ψ_flat =
    Σ_t·ψ_flat = Q_n`` per ordinate.

    The assertion is PER-ORDINATE, not particle-balance (vv anti-pattern #8 /
    §H3: telescoping global balance holds by construction even when
    per-ordinate balance is wrong — the canonical ERR-006 hide).  ≥2G per
    Cardinal Rule.
    """

    def _build_flat_state(
        self, sn_mesh: SNMesh, phi: np.ndarray, W: float,
    ) -> TimedFullField:
        """Flat ψ_n,g = φ_g/W everywhere, with a consistent uniform trace."""
        N = sn_mesh.quad.N
        ng = phi.size
        nx, ny = sn_mesh.spatial_shape
        state = TimedFullField.zeros(
            bulk=AngularFlux, boundary=BoundaryFlux, mesh=sn_mesh,
        )
        state.bulk.values[...] = (phi / W)[None, :, None, None] * np.ones(
            (N, ng, nx, ny)
        )
        # Consistent reflective trace: every boundary face carries ψ_n,g = φ_g/W
        # (the reflection of a spatially-flat field is itself).
        for face in sn_mesh.trace.face_names:
            fv = state.boundary.face_view(face)
            for g in range(ng):
                fv[:, g, :] = phi[g] / W
        return state

    def test_per_ordinate_flat_flux_residual(self):
        """[L0] (L_full + C − B)·ψ_flat − Q has per-ordinate residual < 1e-12.

        Builds the flat infinite-medium fixed point and feeds it to the
        COMPOSED ``(L + C − B)`` action via the production operators
        (``StreamingOperator`` + ``CollisionOperator`` + ``SNBoundaryOperator``).
        Asserts the residual is < 1e-12 for EVERY ordinate (not the weighted
        sum).  Pure-attenuation ⇒ no scattering operator in the residual; the
        closed form ``φ = Q/Σ_t`` is exact.
        """
        sn_mesh = _homogeneous_reflective_2d(nx=4, ny=4)
        quad = sn_mesh.quad
        N, ng = quad.N, 2
        nx, ny = sn_mesh.spatial_shape
        W = float(quad.weights.sum())

        # Pure attenuation: treat Σ_t as the full collision (Σ_s = 0 in the
        # residual — no S operator).  Per-group asymmetric Σ_t exercises the
        # group axis (Cardinal Rule ≥2G).
        sig_t = np.zeros((ng, nx, ny))
        sig_t[0] = 0.5
        sig_t[1] = 1.0
        # Asymmetric per-group scalar source so the two groups differ.
        Q_scalar = np.array([2.0, 3.0])
        phi = Q_scalar / np.array([0.5, 1.0])  # φ_g = Q_g/Σ_{t,g}

        state = self._build_flat_state(sn_mesh, phi, W)

        L = StreamingOperator(sn_mesh)
        C = CollisionOperator(sn_mesh, sig_t)
        B = SNBoundaryOperator(sn_mesh)

        l_bulk = L.apply(state).bulk.values
        c_bulk = C.apply(state).bulk.values
        b_bulk = B.apply(state).bulk.values  # B has zero bulk action

        # Per-ordinate source Q_n,g = Q_g/W (the producer-side /W projection).
        Qn = (Q_scalar / W)[None, :, None, None] * np.ones((N, ng, nx, ny))
        resid = l_bulk + c_bulk - b_bulk - Qn  # (N, ng, nx, ny)

        # PER-ORDINATE residual (vv anti-pattern #8): max over (g, x, y) per n.
        per_ordinate = np.abs(resid).reshape(N, -1).max(axis=1)
        worst_n = int(per_ordinate.argmax())
        assert per_ordinate.max() < 1e-12, (
            f"per-ordinate streaming-equilibrium residual {per_ordinate.max():.3e} "
            f"≥ 1e-12 at ordinate {worst_n} — the extracted 2-D boundary "
            f"coupling breaks the flat-flux fixed point per ordinate (the "
            f"active-trace ~10% non-uniform-mode signature; missing-factor "
            f"failure mode #3).  Per-ordinate max residuals: "
            f"min={per_ordinate.min():.3e}, max={per_ordinate.max():.3e}."
        )

    def test_no_boundary_ring_spike(self):
        """[L0] Flat-flux residual has no boundary-ring spike.

        The active-trace miss converged to a non-uniform mode — its signature
        on a FLAT input is a residual SPIKE at the boundary ring where the
        inconsistent trace coupling fails to cancel.  Assert
        ``max(residual over boundary cells) ≤ 2× median(interior-cell
        residual)`` (the 2-D analogue of the curvilinear pole-spike detector).
        """
        sn_mesh = _homogeneous_reflective_2d(nx=4, ny=4)
        quad = sn_mesh.quad
        N, ng = quad.N, 2
        nx, ny = sn_mesh.spatial_shape
        W = float(quad.weights.sum())

        sig_t = np.zeros((ng, nx, ny))
        sig_t[0] = 0.5
        sig_t[1] = 1.0
        Q_scalar = np.array([2.0, 3.0])
        phi = Q_scalar / np.array([0.5, 1.0])

        state = self._build_flat_state(sn_mesh, phi, W)
        L = StreamingOperator(sn_mesh)
        C = CollisionOperator(sn_mesh, sig_t)
        B = SNBoundaryOperator(sn_mesh)

        Qn = (Q_scalar / W)[None, :, None, None] * np.ones((N, ng, nx, ny))
        resid = (
            L.apply(state).bulk.values
            + C.apply(state).bulk.values
            - B.apply(state).bulk.values
            - Qn
        )
        # Per-cell residual magnitude (max over ordinates + groups).
        per_cell = np.abs(resid).max(axis=(0, 1))  # (nx, ny)

        boundary_mask = np.zeros((nx, ny), bool)
        boundary_mask[0, :] = boundary_mask[-1, :] = True
        boundary_mask[:, 0] = boundary_mask[:, -1] = True
        interior = per_cell[~boundary_mask]
        # Guard a denominator floor: the residual is ~machine-eps everywhere,
        # so use an absolute floor to avoid a divide-by-tiny artefact.
        interior_median = float(np.median(interior))
        boundary_max = float(per_cell[boundary_mask].max())
        floor = 1e-14
        assert boundary_max <= 2.0 * interior_median + floor, (
            f"boundary-ring residual {boundary_max:.3e} spikes above 2× "
            f"interior median {interior_median:.3e} — fingerprint of an "
            f"inconsistent boundary trace coupling (the active-trace "
            f"non-uniform-mode signature on a flat input)."
        )


# ═════════════════════════════════════════════════════════════════════
# GATE R-2D — the new 2-D boundary residual drives to zero at convergence.
# ═════════════════════════════════════════════════════════════════════
#
# Pillar:  closed-form (the residual ``inflow − (R·G·outflow + q_inflow)`` is
#          identically zero at the converged ψ of the consistent system).
# Layer:   trace-equation consistency.
# Catches: Mode-8 passive-boundary blind spot — a boundary residual that is
#          absent / mis-indexed / zero-by-construction (the bug O.4b fixes).
#          PAIRED with the L11 negative control below.


@pytest.mark.l1
class TestBoundaryResidual2DDrivesToZero:
    """At a converged 2-D reflective Krylov solve, the boundary balance holds.

    The boundary-balance defect ``ψ.boundary.inflow − (R·G·ψ.boundary.outflow
    + q_inflow)`` must be < 10× solver_tol at convergence (q_inflow = 0 for
    reflective).  ``R·G·ψ.outflow`` is computed via the canonical
    :class:`~orpheus.sn.boundary_operator.SNBoundaryOperator` (single source of
    truth) — ``B.apply(ψ)`` emits the reflective coupling on the inflow slots
    (the ``A_ss`` block ``V_outflow → V_inflow``).
    """

    def test_reflective_boundary_balance_at_convergence(self):
        """[L1] Reflective 2-D converged ψ satisfies ψ.inflow = B·ψ.outflow.

        Converge a homogeneous reflective 2-D Krylov eigenvalue solve to
        keff/flux_tol 1e-10, take the converged angular flux (typed
        ``TimedFullField``), and assert the per-face inflow-balance defect is
        < 10× solver_tol.  Reads the inflow slots type-agnostically via
        ``trace.inflow_indices_for_face`` + ``face_view`` (survives the future
        BoundaryResidual retype).
        """
        from orpheus.sn.solver import solve_sn

        sn_mesh = _homogeneous_reflective_2d(nx=4, ny=4)
        solver_tol = 1e-10
        res = solve_sn(
            materials={0: get_mixture("A", "2g")},
            mesh=sn_mesh.mesh,
            quadrature=sn_mesh.quad,
            inner_solver="krylov",
            keff_tol=solver_tol, flux_tol=solver_tol,
            max_outer=200, max_inner=200,
        )
        # Sanity: the eigenvalue must be k_inf (Gate K's anchor) — if the solve
        # converged to a non-uniform mode the boundary balance is meaningless.
        assert abs(res.keff - 1.875) < 1e-6, (
            f"precondition: 2-D reflective keff = {res.keff:.10f} ≠ k_inf "
            f"1.875 — the converged mode is non-uniform; the boundary-balance "
            f"check below would be checking the wrong fixed point."
        )

        psi = res.angular_flux  # TimedFullField (bulk + boundary)
        mesh = res.mesh
        trace = mesh.trace
        # R·G·ψ.outflow via the canonical operator (inflow slots carry it).
        B = SNBoundaryOperator(mesh)
        coupled = B.apply(psi)  # boundary-only; inflow slots = R·G·ψ.outflow

        worst_face, worst_defect = None, 0.0
        for face in trace.face_names:
            in_idx = trace.inflow_indices_for_face(face)
            if not in_idx.size:
                continue
            given_inflow = psi.boundary.face_view(face)[in_idx]
            rg_outflow = coupled.boundary.face_view(face)[in_idx]
            defect = float(np.abs(given_inflow - rg_outflow).max())
            if defect > worst_defect:
                worst_face, worst_defect = face, defect

        assert worst_defect < 10.0 * solver_tol, (
            f"2-D reflective boundary-balance defect {worst_defect:.3e} "
            f"(face {worst_face}) ≥ 10× solver_tol {10 * solver_tol:.3e} — "
            f"the converged inflow does NOT equal B·ψ.outflow; the boundary "
            f"residual is not being driven to zero by the outer loop."
        )


# ═════════════════════════════════════════════════════════════════════
# GATE R-2D-NEG — the MANDATORY L11 negative control (Mode-8 catcher).
# ═════════════════════════════════════════════════════════════════════
#
# Pillar:  closed-form (the bare matvec's boundary residual is COMPUTED from
#          ψ.boundary — perturbing the inflow MUST move the residual).
# Layer:   trace-equation consistency (the assertion's load-bearingness).
# Catches: Mode-8 passive-boundary blind spot — proves Gate R-2D is not
#          self-satisfying (a zero-by-construction boundary slot would pass the
#          positive assertion vacuously, the ERR-051 self-referential trap).


@pytest.mark.l1
class TestBoundaryResidual2DResponds:
    """The bare 2-D matvec's boundary residual RESPONDS to a ψ.boundary perturbation.

    E0-T3 is the prototype; this promotes it to a permanent
    production-operator test using the real ``StreamingOperator``.  Without
    this negative control, Gate R-2D is self-satisfying: a boundary residual
    hardcoded to zero (the pre-extraction state) would pass any bulk-only test
    vacuously (ERR-051 lesson — the negative control is the test that makes
    Mode-8 observable).

    Paired positive control: perturbing a NON-inflow slot (the outflow) leaves
    the inflow-identity slot unchanged (the inflow-identity slot is exactly
    ``ψ.inflow``, structurally independent of the outflow).
    """

    def _perturbable_state(self, sn_mesh: SNMesh, seed: int) -> TimedFullField:
        rng = np.random.default_rng(seed)
        state = TimedFullField.zeros(
            bulk=AngularFlux, boundary=BoundaryFlux, mesh=sn_mesh,
        )
        state.bulk.values[...] = rng.standard_normal(state.bulk.values.shape)
        state.boundary.values[...] = rng.standard_normal(
            state.boundary.values.shape
        )
        return state

    def _copy_state(self, src: TimedFullField, sn_mesh: SNMesh) -> TimedFullField:
        dst = TimedFullField.zeros(
            bulk=AngularFlux, boundary=BoundaryFlux, mesh=sn_mesh,
        )
        dst.bulk.values[...] = src.bulk.values
        dst.boundary.values[...] = src.boundary.values
        return dst

    def test_inflow_perturbation_moves_residual(self):
        """[L1] Perturbing ψ.boundary.inflow by δ moves the boundary residual.

        Perturb ONLY the xmin inflow ordinate slots by δ ≠ 0; assert:
          (1) the inflow-identity slot moves by EXACTLY δ (the matvec emits
              the identity ``ψ.inflow`` on the inflow slots — so the residual
              is COMPUTED from ψ.boundary, not structurally zero);
          (2) the boundary-ring bulk responds (> 1e-9) — the inflow feeds the
              wavefront, so the residual is genuinely active.
        """
        sn_mesh = _homogeneous_reflective_2d(nx=4, ny=4)
        sig_t = _sigt_2g(sn_mesh)
        L = StreamingOperator(sn_mesh)
        trace = sn_mesh.trace

        base = self._perturbable_state(sn_mesh, seed=7)
        pert = self._copy_state(base, sn_mesh)
        delta = 0.37
        in_idx = trace.inflow_indices_for_face("xmin")
        pert.boundary.face_view("xmin")[in_idx] += delta

        r_base = L.apply(base)
        r_pert = L.apply(pert)

        # (1) inflow-identity slot moves by exactly δ.
        dident = (
            r_pert.boundary.face_view("xmin")[in_idx]
            - r_base.boundary.face_view("xmin")[in_idx]
        )
        np.testing.assert_allclose(
            dident, delta, atol=1e-13,
            err_msg=(
                "the bare 2-D matvec inflow-identity slot did NOT move by δ "
                "under a ψ.boundary.inflow perturbation — the boundary "
                "residual is NOT computed from ψ.boundary (Mode-8: passive / "
                "structurally-zero boundary block)."
            ),
        )

        # (2) the boundary-ring bulk responds (the inflow feeds the wavefront).
        dbulk = float(np.abs(r_pert.bulk.values - r_base.bulk.values).max())
        assert dbulk > 1e-9, (
            f"perturbing ψ.boundary.inflow left the bulk UNCHANGED "
            f"(max|Δbulk| = {dbulk:.3e} ≤ 1e-9) — the bare matvec is NOT "
            f"reading the inflow trace into the wavefront (Mode-8: passive "
            f"boundary)."
        )

    def test_outflow_perturbation_leaves_inflow_identity_unchanged(self):
        """[L1] Positive control: perturbing the OUTFLOW slot does not move the
        inflow-identity slot.

        The inflow-identity slot emits exactly ``ψ.inflow`` — it is structurally
        independent of the outflow.  Perturbing ONLY the xmin outflow slots must
        leave the inflow-identity slot bit-stable.  This pins that the response
        in :meth:`test_inflow_perturbation_moves_residual` is SPECIFIC to the
        inflow (not a global "everything moves" artefact).
        """
        sn_mesh = _homogeneous_reflective_2d(nx=4, ny=4)
        sig_t = _sigt_2g(sn_mesh)
        L = StreamingOperator(sn_mesh)
        trace = sn_mesh.trace

        base = self._perturbable_state(sn_mesh, seed=7)
        pert = self._copy_state(base, sn_mesh)
        delta = 0.37
        in_idx = trace.inflow_indices_for_face("xmin")
        out_idx = trace.outflow_indices_for_face("xmin")
        pert.boundary.face_view("xmin")[out_idx] += delta

        r_base = L.apply(base)
        r_pert = L.apply(pert)

        dident = (
            r_pert.boundary.face_view("xmin")[in_idx]
            - r_base.boundary.face_view("xmin")[in_idx]
        )
        np.testing.assert_allclose(
            dident, 0.0, atol=1e-13,
            err_msg=(
                "perturbing the OUTFLOW slot moved the inflow-identity slot — "
                "the inflow-identity emission is NOT ``ψ.inflow`` alone (it "
                "leaked an outflow dependence, breaking the block structure)."
            ),
        )
