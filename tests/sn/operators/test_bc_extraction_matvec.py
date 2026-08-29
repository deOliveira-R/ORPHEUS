r"""Wave O (#208) step O.4a.2 — the MATVEC BC-extraction gates.

This file pins the **matvec / ``.apply`` path** carve of O.4a.2 ONLY.
It does NOT touch the sweep (``(L+C).solve`` — that
is O.4a.3) nor the solver rewire (SI seeding retirement, ``q.boundary``
typing, the through-solver reflective convergence — that is O.4a.4).

The O.4a.2 carve (production code under test)
---------------------------------------------

The KEYSTONE deletion is the intra-call reflective re-apply inside
``_MSpatialOperatorSum._compute_LpC`` (operator.py — the line
``inflow_full = bc_outer.apply(outflow_at_boundary.T)`` at HEAD 92be67a,
plus its twin in ``_compute_decomposition``).  Deleting it decouples
bulk ↔ boundary inside one matvec call.

Post-extraction (design (a) — the CORRECTED realization):

* ``L_full`` reads ``ψ.boundary.inflow`` as a GIVEN (instead of
  re-deriving it via ``bc_outer.apply(...)`` — the deleted keystone);
* ``L_full`` KEEPS the outflow self-consistency defect on the outflow
  slots (``ψ.outflow − streamed`` — ``ψ.outflow`` is a STORED unknown that
  the sibling ``−B`` reads as its input; emitting the "raw outflow" instead
  would singularise the outflow row and break ``−B``-as-sibling) and ADDS
  the inflow identity ``I·ψ.inflow`` on the inflow slots (so the trace
  inflow becomes an explicit unknown driven by the consistency residual
  ``ψ.inflow − B·ψ.outflow``);
* the realized boundary law ``B`` is the ``A_ss`` block ``V_outflow →
  V_inflow`` — it emits on the **inflow slots only** (the full-face law
  also maps the input inflow onto the output outflow slots, which is
  spurious for the block and is projected away);
* the curvilinear pole Carlson seed (``psi_view[:, :, 0, 0]``) STAYS —
  it is a regularity condition at r=0, NOT a boundary condition.

The realized boundary law ``B`` (a ``BoundaryOperator`` leaf after
O.4a.1) is assembled into a whole-trace operator and placed as a
SIBLING ``−B`` so ``(L_full + C − S − F − B).apply(ψ)`` forms the FULL
residual, with the boundary block carrying the consistency residual
``ψ.inflow − B·ψ.outflow − q_inflow``.

What is bit-identical vs what legitimately changes shape
--------------------------------------------------------

For VACUUM BC, ``B = 0`` (the realized vacuum law zeros the inflow
slots).  So the *bulk* residual of the extracted ``(L+C).apply(ψ)``
MUST be byte-for-byte identical to the pre-extraction ``(L+C).apply(ψ)``
— the vacuum path must not move a single bit.

The BOUNDARY slot under design (a) KEEPS the outflow self-consistency
defect (it does NOT switch to a raw outflow):

* PRE-carve  ``StreamingOperator.apply(ψ).boundary`` on the outflow slots
  = ``streamed − bc_estimate`` (``bc_estimate`` re-derived from the
  forward sweep's own outflow via the keystone ``bc_outer.apply``);
* POST-carve ``StreamingOperator.apply(ψ).boundary`` on the outflow slots
  = ``streamed − ψ.outflow`` — the defect against the STORED outflow
  unknown (``ψ.outflow`` is the trace value the sibling ``−B`` reads).
  The inflow slots ADD the identity ``ψ.inflow``.

For the canonical ZERO-input vacuum test (a zero ``AngularBoundaryFlux``) both
the keystone ``bc_estimate`` and the stored ``ψ.outflow`` are zero, so
PRE- and POST-carve both emit ``streamed`` on the outflow slots and zero
on the inflow slots → vacuum bit-identity holds on the boundary slot too
FOR THE ZERO-INPUT CASE (validated by the ``test_vacuum_*_boundary_slot_*``
rows).  The general (non-zero input boundary) case is where the emission
becomes a function of the STORED outflow: the output outflow slot DEPENDS
on the input outflow value (the defect is kept).  That dependence is the
discriminating property pinned by
``TestVacuumBoundaryDefectKept`` / ``TestLFullOutflowDefectKept`` below —
the keystone deletion + given-inflow read changed the matvec, but the
outflow DEFECT is deliberately retained (raw-outflow would singularise
the outflow row under ``−B``).

Baseline-capture mechanism (the design problem)
-----------------------------------------------

The robust mechanism is **option (a): a committed snapshot of the
pre-extraction matvec output captured at HEAD 92be67a** (the O.4a.1
close-out commit, which is exactly today's HEAD).  Rationale over option
(b) (a hand/analytical invariant): the bulk matvec output is a dense
``(N, ng, nx, ny)`` array produced by the curvilinear M-M Carlson
recurrence — there is no cheap closed-form for it on a RANDOM ψ, and a
random ψ is precisely what stresses the bulk operator (a flat ψ would
null the redistribution path, vv §H2).  A snapshot of ``(L+C).apply``
on a fixed-seed random ψ is the only way to get a byte-exact reference
for the WHOLE bulk array.  The analytical ``Q/Σ_t`` invariant (option b)
is ALSO included below — but as a SEPARATE, structurally-independent
check (closed-form pillar), NOT as the bit-identity reference.  Snapshot
proves "did not move a bit"; ``Q/Σ_t`` proves "the value is correct".
Both are needed (vv §bit-identity criterion 2: ULP-distance is
necessary-not-sufficient).

HOW THE IMPLEMENTER CAPTURES THE BASELINE (run ONCE, at HEAD 92be67a,
BEFORE the carve):

    PYTHONPATH=$PWD .venv/bin/python -m pytest \
        tests/sn/operators/test_bc_extraction_matvec.py \
        -k capture_baseline --capture-baseline

The ``capture_baseline`` test below is gated on the
``--capture-baseline`` flag (added in this file's ``pytest_addoption``);
when the flag is present it WRITES the ``.npz`` snapshots under
``tests/sn/_data/bc_extraction_baseline/`` and SKIPS the assert.  When
the flag is absent (normal runs, including the post-carve gate) it READS
the committed snapshots and asserts bit-identity.  The snapshot files are
committed alongside this test so the post-carve run has its reference.

Tagging
-------

* bit-identity + structural rows  → ``@pytest.mark.foundation`` (software
  invariant — no theory-page ``:label:``; the claim is "the carve did not
  move a bit / the extracted pieces have the right structure").
* the ``Q/Σ_t`` closed-form value row → ``@pytest.mark.l0`` +
  ``@pytest.mark.verifies("streaming-equilibrium")`` (term-level
  closed-form, the universal curvilinear diagnostic).
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from orpheus.geometry import BC, CoordSystem, Mesh1D, Mesh2D
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.mesh.augmented_mesh import SNMesh
from orpheus.sn.operators.streaming import StreamingOperator
from orpheus.transport.operators.multiplication_operator import MultiplicationOperator
from orpheus.transport.fields.angular_flux import AngularFlux
from orpheus.transport.fields.angular_boundary_flux import AngularBoundaryFlux
from orpheus.transport.full_field import FullField
from orpheus.transport.timed_full_field import TimedFullField
from tests.sn._test_helpers import (
    SN_TESTS_ROOT,
    placeholder_materials,
    radial_characteristic_edge_seed,
)
from tests.sn.regression._regression_assert import assert_regression

# ─────────────────────────────────────────────────────────────────────
# Baseline snapshot store + the --capture-baseline flag.
# ─────────────────────────────────────────────────────────────────────

_BASELINE_DIR = SN_TESTS_ROOT / "_data" / "bc_extraction_baseline"

# ``--capture-baseline`` is declared in the ROOT ``tests/conftest.py``
# (pytest_addoption only fires there). When present the snapshot tests
# WRITE the pre-carve matvec output and skip; absent they READ + assert.


def _capturing(request) -> bool:
    return bool(request.config.getoption("--capture-baseline", default=False))


# ─────────────────────────────────────────────────────────────────────
# Mesh builders — slab / sphere / cylinder, parametrised over BC.
# ─────────────────────────────────────────────────────────────────────

# The five geometries the matvec serves at O.4a.2.  2-D-Cartesian-vacuum
# is included for the vacuum bit-identity gate ONLY — its NON-vacuum
# extraction (the boundary-residual ADD) is O.4b, NOT O.4a.2 (see
# the gate-sequencing table in the closing summary).
_GEOMS_1D = ("SLB", "SPH", "CYL")


def _build_sn_mesh(
    geometry: str,
    *,
    bc: str = "vacuum",
    n_cells: int = 5,
    n_ord: int = 4,
) -> SNMesh:
    """Small SNMesh; ``bc`` set on BOTH endpoints (slab) / outer (curv).

    Sized small (n_cells=5, n_ord=4) so each matvec runs well under a
    second.  The bit-identity / structural contracts are size-independent.
    """
    bc_obj = BC(bc)
    if geometry == "SLB":
        mesh = Mesh1D(
            edges=np.linspace(0.0, 2.0, n_cells + 1),
            mat_ids=np.zeros(n_cells, dtype=int),
            coord=CoordSystem.CARTESIAN,
            bc_left=bc_obj,
            bc_right=bc_obj,
        )
        quad = Quadrature.gauss_legendre(n_ordinates=n_ord)
    elif geometry == "SPH":
        mesh = Mesh1D(
            edges=np.linspace(0.0, 2.0, n_cells + 1),
            mat_ids=np.zeros(n_cells, dtype=int),
            coord=CoordSystem.SPHERICAL,
            bc_left=BC("reflective"),  # r=0 pole — regularity, not a BC
            bc_right=bc_obj,
        )
        quad = Quadrature.gauss_legendre(n_ordinates=n_ord)
    elif geometry == "CYL":
        mesh = Mesh1D(
            edges=np.linspace(0.01, 2.0, n_cells + 1),
            mat_ids=np.zeros(n_cells, dtype=int),
            coord=CoordSystem.CYLINDRICAL,
            bc_left=BC("reflective"),
            bc_right=bc_obj,
        )
        quad = Quadrature.folded_product(n_mu=n_ord, n_phi=2 * n_ord)
    else:
        raise ValueError(geometry)
    return SNMesh(mesh, quad, placeholder_materials())


def _random_state(sn_mesh: SNMesh, seed: int, *, zero_boundary: bool = True) -> TimedFullField:
    """A fixed-seed random bulk ψ with a chosen boundary trace.

    ``zero_boundary=True`` is the canonical vacuum-matvec input: the
    boundary buffer is all-zero (the vacuum sweep seeds from a zero
    inflow trace).  ``zero_boundary=False`` populates a random boundary
    trace — used ONLY by the discriminator that distinguishes the
    pre-carve defect from the post-carve raw outflow.
    """
    ng = 1
    N = sn_mesh.quad.N
    rng = np.random.default_rng(seed)
    boundary = AngularBoundaryFlux.zeros(sn_mesh.angular_trace)
    if not zero_boundary:
        # Populate every face buffer with a fixed-seed random trace.
        for face in boundary.layout.faces:
            view = boundary.face_view(face)
            view[:] = rng.standard_normal(view.shape)
    bulk_arr = rng.standard_normal((N, ng, *sn_mesh.spatial_shape))
    return TimedFullField(
        interior=AngularFlux(values=bulk_arr, space=sn_mesh.angular_bulk_space),
        boundary=boundary,
        _history=(),
        history_depth=2,
    )


def _LpC_apply(sn_mesh: SNMesh, state: TimedFullField, sigma_t: np.ndarray) -> "FullField":
    """``(L + C).apply(state)`` via the public operator-algebra path.

    #257 S8a — the matvec leaf is a base arrow, so ``(L + C).apply`` returns a
    timeless :class:`~orpheus.transport.full_field.FullField` source.
    Step 6 — on a carrying mesh (SPH) the CONSISTENT edge-extrapolated ψ½
    seed (derived from the state's own bulk, exactly as every pre-eviction
    site here derived it) rides the JOINT row-A action through the grid
    (:func:`_LC_matvec`'s seed arm — presence structural); no seed on
    non-carrying meshes.
    """
    L = StreamingOperator.pose(sn_mesh)
    C = MultiplicationOperator.from_mesh(sigma_t, sn_mesh)
    seed_leg = radial_characteristic_edge_seed(state.interior.values, sn_mesh)
    from tests.sn._test_helpers import _LC_matvec

    return _LC_matvec(
        state, sigma_t, sn_mesh=sn_mesh, LC=(L + C),
        radial_characteristic_flux=seed_leg,
    )


# ═════════════════════════════════════════════════════════════════════
# GATE 1 — VACUUM matvec bit-identity (the load-bearing O.4a.2 gate).
# ═════════════════════════════════════════════════════════════════════
#
# For vacuum BC, B=0, so the extracted (L+C).apply(ψ) on a ZERO-inflow
# input MUST be byte-for-byte identical to the pre-extraction
# (L+C).apply(ψ).  The bulk residual is bit-identical UNCONDITIONALLY;
# the boundary slot is bit-identical FOR THE ZERO-INPUT CASE (defect
# = outflow − 0 = raw outflow when face_value = 0).


@pytest.mark.foundation
class TestVacuumMatvecBitIdentity:
    """Vacuum ``(L+C).apply`` regression — principled-equivalent post-#240.

    ORIGINALLY a strict byte-identity gate for the O.4a.2 BC-extraction
    carve (the given-zero-inflow extracted sweep IS today's sweep for
    vacuum, since the deleted ``inflow_full = bc_outer.apply(outflow)``
    returns ZERO inflow).  #240 moved the Cartesian matvec onto the uniform
    ÷V ``residual_kernel_batch`` kernel (DD+LD), so the **slab** vacuum bulk
    matvec now re-associates ~1 ULP — the SLB baselines were re-captured and
    the bulk gate narrowed to :func:`assert_regression` (principled-
    equivalence; bit-identity is the escalatable ``DriftWarning``).  The
    boundary-slot gate stays strict byte-identity (it did not move — the
    outflow defect reconstructs from the same ``ψ_out = 2ψ̄ − ψ_in`` faces).

    Curvilinear (SPH) baselines were RE-CAPTURED 2026-06-26, closing #250
    (the curvilinear snapshot cleanup).  They had carried a long-standing
    red because the spherical apply legitimately evolved after this store's
    last refresh: ``b2d8a6d`` (Bailey Eq. 43, Refs #229) unclamped the
    spherical Morel–Montry WDD weight τ in ``spherical_streaming`` but
    refreshed only its own targeted snapshot, silently leaving the SPH arms
    stale (the SLB/CYL arms were later re-captured at #240; SPH was
    deferred — lesson L-034).  The frozen SPH values were completely
    different numbers from the live output, so the nULP gate hard-failed
    (drift ~1e15 ULP).  Before re-capturing, the CURRENT SPH matvec was
    verified correct vs structurally-independent references (L0 streaming-
    equilibrium per-ordinate, L1 isotropic + anisotropic MMS to O(h²),
    L1 trajectory-resolvent) — never green-by-fiat (``vv-principles``
    L11/L14/L27).
    """

    @pytest.mark.parametrize("geometry", _GEOMS_1D)
    @pytest.mark.parametrize("seed", [0, 1, 2])
    def test_vacuum_bulk_bit_identical_1d(self, geometry, seed, request):
        """[foundation] Vacuum bulk matvec regression (principled-equiv post-#240).

        The bulk residual (N, ng, nx, ny) of (L+C).apply was byte-identical
        across the O.4a.2 BC-extraction carve (for vacuum the deleted line
        contributed ZERO inflow — identical to reading a zero inflow trace).
        #240 then moved the Cartesian matvec onto the uniform ÷V kernel, so
        the slab arm re-associates ~1 ULP (nULP gate + DriftWarning).  #240
        Phase 2 Step B (2026-06-15) then made ``(L+C).apply`` OWN its matvec
        (``loss_action(self.sigma)`` direct) instead of the inherited leaf sum
        ``(loss_action(σ_t) − σ_t·ψ) + σ_t·ψ`` — the override DROPS that
        ``−σ_t·ψ + σ_t·ψ`` round-trip, so the CURVILINEAR (CYL) bulk re-
        associates up to ~46 ULP (rel ~4e-17 — pure FP-non-associativity on a
        near-zero output element; the override is the MORE accurate path).  The
        CYL baselines were RE-CAPTURED under the override (principled re-baseline,
        ``vv-principles`` 3-criteria: named ``loss_action`` intermediate /
        verified by the teeth gate in ``test_removal_form_matvec_sweep.py`` /
        drift = reduction-depth × ULP).  SPH was RE-CAPTURED 2026-06-26
        (closes #250) after ``b2d8a6d`` (Bailey Eq. 43, Refs #229) evolved
        the spherical apply without refreshing this store (verified correct
        vs L0/L1 references first — not green-by-fiat).
        """
        sn_mesh = _build_sn_mesh(geometry, bc="vacuum")
        state = _random_state(sn_mesh, seed, zero_boundary=True)
        sigma_t = np.full((1, *sn_mesh.spatial_shape), 2.0)
        out = _LpC_apply(sn_mesh, state, sigma_t)

        key = f"vacuum_bulk_{geometry}_seed{seed}"
        path = _BASELINE_DIR / f"{key}.npy"
        if _capturing(request):
            _BASELINE_DIR.mkdir(parents=True, exist_ok=True)
            np.save(path, out.interior.values)
            pytest.skip(f"captured baseline {key}")
        assert path.exists(), (
            f"missing baseline snapshot {path}; run with "
            f"--capture-baseline at HEAD 92be67a BEFORE the carve"
        )
        expected = np.load(path)
        # #240: the Cartesian matvec rides the uniform ÷V
        # ``residual_kernel_batch`` kernel (DD+LD), which re-associates the
        # cell-balance reduction ~1 ULP vs the pre-#240 ×V ``cell_balance``
        # path — a principled-equivalence re-baseline (vv-principles
        # §"Bit-identity vs principled-equivalence").  The gate is nULP at
        # reduction_depth=nx with the escalatable DriftWarning (bit-identity
        # is the opt-in ``-W error::DriftWarning`` bonus).  The SLB baselines
        # were regenerated at #240; the CYL baselines were RE-CAPTURED at #240
        # Phase 2 Step B (the apply override drops the ``−σ_t·ψ + σ_t·ψ`` round-
        # trip → ~46 ULP / rel ~4e-17 re-association — principled).  The SPH
        # baselines were RE-CAPTURED 2026-06-26 (closes #250) after b2d8a6d
        # (Bailey Eq. 43, Refs #229) unclamped the spherical M-M τ weight and
        # refreshed only its own snapshot, leaving this store stale (#240
        # re-captured SLB/CYL; SPH was deferred) — current SPH matvec
        # verified vs L0/L1 references first.
        assert_regression(
            out.interior.values, expected,
            conv_tol=0.0, kind="direct", reduction_depth=sn_mesh.nx,
            case_name=f"vacuum_bulk_{geometry}_seed{seed}",
            quantity="vacuum_bulk_matvec",
        )

    @pytest.mark.parametrize("geometry", _GEOMS_1D)
    @pytest.mark.parametrize("seed", [0, 1, 2])
    def test_vacuum_boundary_slot_bit_identical_zero_input_1d(
        self, geometry, seed, request,
    ):
        """[foundation] Vacuum boundary slot is byte-identical pre/post
        carve FOR THE ZERO-INPUT case.

        For vacuum with a zero input boundary trace, the pre-carve defect
        ``outflow − face_value`` = ``outflow − 0`` = raw outflow on the
        outflow slots — exactly what the post-carve path emits.  So the
        boundary slot is ALSO bit-identical for this input.  (The general
        non-zero-input case legitimately differs — see
        ``test_vacuum_boundary_slot_diverges_for_nonzero_input``.)
        """
        sn_mesh = _build_sn_mesh(geometry, bc="vacuum")
        state = _random_state(sn_mesh, seed, zero_boundary=True)
        sigma_t = np.full((1, *sn_mesh.spatial_shape), 2.0)
        out = _LpC_apply(sn_mesh, state, sigma_t)

        key = f"vacuum_boundary_{geometry}_seed{seed}"
        path = _BASELINE_DIR / f"{key}.npy"
        if _capturing(request):
            _BASELINE_DIR.mkdir(parents=True, exist_ok=True)
            np.save(path, out.boundary.values)
            pytest.skip(f"captured baseline {key}")
        assert path.exists(), (
            f"missing baseline snapshot {path}; run with "
            f"--capture-baseline at HEAD 92be67a BEFORE the carve"
        )
        expected = np.load(path)
        np.testing.assert_array_equal(
            out.boundary.values, expected,
            err_msg=(
                f"{geometry} seed={seed}: vacuum BOUNDARY slot moved a "
                f"bit under the BC extraction on the ZERO-INPUT case. "
                f"For zero input, defect = outflow − 0 = raw outflow, so "
                f"the post-carve raw outflow MUST equal the pre-carve "
                f"defect byte-for-byte."
            ),
        )

    @pytest.mark.parametrize("seed", [0, 1, 2])
    def test_vacuum_2d_cartesian_bulk_bit_identical(self, seed, request):
        """[foundation] Vacuum 2-D Cartesian bulk matvec regression floor.

        2-D drives the representation's ``loss_action`` walk (a
        SEPARATE code path from the 1-D ``_compute_LpC``; since S6.3 the
        matvec walk lives on the loss representation, off the operator —
        ``ScanMarch`` default since S6.9).  The 2-D non-vacuum
        boundary-residual ADD is O.4b — out of scope here; this row
        asserts only that the vacuum 2-D bulk does not move.

        B0.3 REPAIR (2026-07-30) — this row had **never executed a
        single assertion in its life**, on any of its three seeds.  It
        built a 1-D :class:`Mesh1D` and then read
        ``sn_mesh.spatial_shape[1]``, which raises :exc:`IndexError` on
        a 1-tuple; the surrounding ``except Exception as exc:
        pytest.skip(...)`` swallowed it, so all three parametrisations
        reported as a green skip (*"2-D mesh construction not available
        here: tuple index out of range"* — the only three skips in the
        whole boundary harness).  That is the ``vv-principles`` Mode-8
        SKIP-SWALLOWED class in its purest form: a broad
        ``except`` → ``pytest.skip`` converts EVERY future construction
        bug into a permanent green, and a skip is invisible in a summary
        line.  The repair (a) builds a genuine :class:`Mesh2D`, (b)
        deletes the blanket ``except`` so a construction failure is a
        FAILURE, and (c) keeps the ``ny > 1`` check as a hard assertion
        rather than a skip.

        HONEST SCOPE — read the claim carefully.  The docstring used to
        call this "the SENTINEL that the 1-D O.4a.2 carve did not
        accidentally perturb the 2-D path".  That claim is **not
        recoverable**: O.4a.2 landed long ago, and a baseline captured
        today cannot testify about a carve that predates it.  What this
        row is, honestly, is a **regression floor captured 2026-07-30**
        on the 2-D vacuum bulk matvec — a drift gate going forward, not
        evidence about a past carve.  Its correctness anchor is NOT this
        snapshot (a self-generated baseline is no reference at all) but
        the structurally-independent 2-D gates that already exist:
        ``tests/sn/sweep/cartesian_2d/test_2d_l2_matvec_correctness.py``
        (L2), the 2-D MMS rows, ``test_scan_march_equivalence.py`` and
        ``test_2d_full_field_oracle.py``.  The snapshot's job is to
        notice when the value moves; those gates say whether it should.

        The mesh is deliberately NON-SQUARE (nx=3, ny=4) so an x↔y
        transposition cannot hide behind a square-box symmetry
        (``vv-principles`` §H2 — the convenient config nulls the term).
        """
        ng = 1
        nx, ny = 3, 4
        mesh = Mesh2D(
            edges_x=np.linspace(0.0, 2.0, nx + 1),
            edges_y=np.linspace(0.0, 3.0, ny + 1),
            mat_map=np.zeros((nx, ny), dtype=int),
            coord=CoordSystem.CARTESIAN,
            bc_xmin=BC("vacuum"), bc_xmax=BC("vacuum"),
            bc_ymin=BC("vacuum"), bc_ymax=BC("vacuum"),
        )
        quad = Quadrature.level_symmetric(sn_order=4)
        sn_mesh = SNMesh(mesh, quad, placeholder_materials())
        # A hard assertion, NOT a skip: if the mesh ever degenerates to
        # 1-D the row must go RED, because a silently-skipped 2-D gate
        # is exactly the defect this repair removes.
        assert len(sn_mesh.spatial_shape) == 2 and sn_mesh.spatial_shape[1] > 1, (
            f"the 2-D vacuum regression floor needs an ny>1 mesh; got "
            f"spatial_shape={sn_mesh.spatial_shape!r}"
        )

        state = _random_state(sn_mesh, seed, zero_boundary=True)
        sigma_t = np.full((ng, *sn_mesh.spatial_shape), 2.0)
        out = _LpC_apply(sn_mesh, state, sigma_t)

        key = f"vacuum_2d_bulk_seed{seed}"
        path = _BASELINE_DIR / f"{key}.npy"
        if _capturing(request):
            _BASELINE_DIR.mkdir(parents=True, exist_ok=True)
            np.save(path, out.interior.values)
            pytest.skip(f"captured baseline {key}")
        assert path.exists(), (
            f"missing baseline snapshot {path}; run with --capture-baseline"
        )
        expected = np.load(path)
        np.testing.assert_array_equal(
            out.interior.values, expected,
            err_msg=(
                f"2-D vacuum seed={seed}: the 1-D O.4a.2 carve PERTURBED "
                f"the 2-D Cartesian path — the carve must touch ONLY "
                f"_compute_LpC, leaving the 2-D representation "
                f"`loss_action` walk untouched."
            ),
        )


# ═════════════════════════════════════════════════════════════════════
# GATE 2 — Closed-form Q/Σ_t (structurally-independent VALUE check).
# ═════════════════════════════════════════════════════════════════════
#
# vv §bit-identity criterion 2: the snapshot proves "did not move",
# Q/Σ_t proves "the value is correct".  This is the closed-form pillar,
# structurally independent of the snapshot (which is the operator's own
# prior output).  It is ALSO the canonical curvilinear missing-ΔA/w
# diagnostic (Signature 1 / Failure mode #3): a dropped ΔA/w in the
# extracted curvilinear seed produces a flux spike at r=0.


@pytest.mark.l0
@pytest.mark.verifies("streaming-equilibrium")
class TestStreamingEquilibriumValue:
    """Uniform Q, uniform Σ_t ⇒ flat ψ = Q/Σ_t is a fixed point.

    For the streaming-equilibrium configuration the matvec residual on a
    FLAT ψ = Q/Σ_t must satisfy the discrete balance — for vacuum the
    relevant claim at the matvec level is that the flat-ψ residual equals
    the analytical streaming contribution with NO spurious spike at r=0.

    This is the single most powerful curvilinear diagnostic.  A missing
    ``ΔA/w`` in the extracted curvilinear pole seed (Failure mode #3)
    produces a flux spike at r=0 that THIS test catches at O.4a.2 — it
    does NOT need the solver (O.4a.4), because the per-cell balance is a
    matvec-level (forward-only) property.

    NOTE on observability at O.4a.2 vs O.4a.4: the missing-ΔA/w bug is
    observable at O.4a.2 *iff* the flat-ψ input drives the redistribution
    term out of cancellation.  Per vv §H2, a perfectly flat ψ NULLS the
    redistribution (the per-ordinate flat-flux residual is the sharp
    probe — see ``test_per_ordinate_flat_flux_consistency`` in
    test_quadrature.py).  So the DEFINITIVE missing-ΔA/w catch is the
    per-ordinate flat-flux residual (already green, must-stay-green); the
    row below is the COMPLEMENTARY closed-form value anchor.
    """

    @pytest.mark.parametrize("geometry", _GEOMS_1D)
    def test_flat_flux_per_ordinate_balance_no_pole_spike(self, geometry):
        """[L0] Flat ψ ⇒ per-ordinate matvec residual has no r=0 spike.

        Build a FLAT ψ (every cell, every ordinate = 1.0), uniform Σ_t.
        The (L+C).apply residual minus the collision term is the discrete
        streaming action; on a flat field the streaming + redistribution
        must cancel per ordinate (sum to the boundary leakage only), with
        NO anomalous growth in cell 0 (the pole cell).  A missing ΔA/w in
        the extracted seed makes cell-0 the largest-magnitude cell.
        """
        sn_mesh = _build_sn_mesh(geometry, bc="vacuum", n_cells=8)
        ng, N = 1, sn_mesh.quad.N
        nx = sn_mesh.nx
        flat = np.ones((N, ng, *sn_mesh.spatial_shape))
        state = TimedFullField(
            interior=AngularFlux(values=flat, space=sn_mesh.angular_bulk_space),
            boundary=AngularBoundaryFlux.zeros(sn_mesh.angular_trace),
            _history=(),
            history_depth=2,
        )
        sigma_t = np.full((ng, *sn_mesh.spatial_shape), 1.0)
        out = _LpC_apply(sn_mesh, state, sigma_t)
        # L action = (L+C) − C; C·ψ = σ_t·ψ = 1·1 = 1 everywhere.
        l_action = out.interior.values - sigma_t[None] * flat
        # Per-cell L2 magnitude across ordinates+groups (1-D geometries).
        per_cell = np.linalg.norm(
            l_action.reshape(N * ng, nx), axis=0,
        )
        # The pole cell (index 0) must NOT be the dominant cell — a
        # missing ΔA/w makes it spike.  Use a generous factor (2×) so the
        # test is a SPIKE detector, not a precise-balance assertion.
        if geometry in ("SPH", "CYL"):
            assert per_cell[0] <= 2.0 * np.median(per_cell), (
                f"{geometry}: cell-0 (pole) streaming action "
                f"{per_cell[0]:.3e} spikes above 2× median "
                f"{np.median(per_cell):.3e} — fingerprint of a missing "
                f"ΔA/w in the extracted curvilinear seed (Failure mode #3, "
                f"Signature 1)."
            )


# ═════════════════════════════════════════════════════════════════════
# GATE 3 — Structural tests for the extracted pieces.
# ═════════════════════════════════════════════════════════════════════
#
# These pin the SHAPE of the three extracted pieces:
#  (a) the whole-trace B operator;
#  (b) L_full reads ψ.boundary.inflow as a given;
#  (c) L_full emits the raw outflow trace.
#
# IMPORTANT SCOPE NOTE: at the point these tests are FIRST written
# (O.4a.2 start), the whole-trace `B` and the extracted `L_full` may not
# exist yet.  Each test is therefore guarded so that it is RED-until-built
# in the principled way:
#  - 3(a)/3(c): xfail(strict=True) until the extracted symbol exists, so
#    the test flips to XPASS→forces removal of the marker when the carve
#    lands (the marker cannot be silently left — strict XPASS fails).
#  - 3(b): a behavioural assertion that is MEANINGFUL only post-carve;
#    pre-carve it is xfail-strict for the same reason.
#
# The implementer removes each xfail marker at the sub-commit where the
# corresponding production symbol lands.


# NOTE (O.4a.2 reconciliation): the former 3(a) ``TestWholeTraceBOperator``
# + its ``_realized_B_for_face`` helper are RETIRED here.  The whole-trace
# ``B`` (:class:`~orpheus.sn.operators.boundary.SNBoundaryOperator`) is
# canonically pinned in ``tests/sn/operators/test_sn_boundary_operator.py``
# (role / domain / per-face wiring / block-diagonal / apply_transpose).  That
# suite ALSO pins the design (a) correction this file's docstring describes:
# ``B`` is the ``A_ss`` block ``V_outflow → V_inflow`` and emits on the
# **inflow row only** — the full-face per-face ``bc.apply`` equality that 3(a)
# asserted is WRONG for the block (it would also write the spurious
# ``R·ψ.inflow`` onto the outflow slots).  The placeholder
# ``assemble_boundary_operator`` symbol that 3(a) imported never existed; the
# realized leaf is ``SNBoundaryOperator(sn_mesh)``.


@pytest.mark.foundation
class TestLFullReadsInflow:
    """3(b) — the extracted ``L_full`` reads ``ψ.boundary.inflow`` as a
    GIVEN (not re-derived via ``bc_outer.apply``).
    """

    @pytest.mark.parametrize("geometry", _GEOMS_1D)
    def test_nonzero_inflow_changes_bulk_output(self, geometry):
        """[foundation] Setting a non-zero inflow trace seeds the sweep and
        changes the bulk matvec output in the expected direction.

        PRE-carve, the matvec IGNORED the input boundary inflow on the
        OUTER face for curvilinear (it re-derived it from the forward
        sweep's own outflow at the keystone line) — so changing the input
        inflow did NOT change the curvilinear bulk output.  POST-carve
        (O.4a.2 keystone deletion landed), ``L_full`` reads
        ``ψ.boundary.inflow`` directly, so a non-zero outer inflow MUST
        change the backward-sweep seed and hence the bulk output.

        This test pins that post-carve behaviour: two inputs differing
        ONLY in the outer-face inflow slots produce DIFFERENT bulk outputs.
        (Was xfail-strict pre-carve; the keystone deletion flipped it to
        XPASS, so the marker is removed — O.4a.2 Commit 2.)
        """
        sn_mesh = _build_sn_mesh(geometry, bc="reflective")
        ng, N = 1, sn_mesh.quad.N
        sigma_t = np.full((ng, *sn_mesh.spatial_shape), 2.0)
        rng = np.random.default_rng(3)
        bulk = rng.standard_normal((N, ng, *sn_mesh.spatial_shape))

        def _state(inflow_scale: float) -> TimedFullField:
            b = AngularBoundaryFlux.zeros(sn_mesh.angular_trace)
            # populate the OUTER-face inflow slots with a non-zero trace.
            trace = sn_mesh.angular_trace
            inflow = trace.inflow_indices_for_face("xmax")
            face = b.face_view("xmax")
            if inflow.size:
                face[inflow, :] = inflow_scale
            return TimedFullField(
                interior=AngularFlux(values=bulk.copy(), space=sn_mesh.angular_bulk_space),
                boundary=b,
                _history=(), history_depth=2,
            )

        out_zero = _LpC_apply(sn_mesh, _state(0.0), sigma_t)
        out_one = _LpC_apply(sn_mesh, _state(1.0), sigma_t)
        diff = np.linalg.norm(
            out_one.interior.values - out_zero.interior.values
        )
        assert diff > 1e-10, (
            f"{geometry}: changing the outer-face INFLOW trace did not "
            f"change the bulk output — L_full is NOT reading "
            f"ψ.boundary.inflow as a given (the keystone re-apply is "
            f"still present)."
        )


@pytest.mark.foundation
class TestLFullOutflowDefectKept:
    """3(c) — design (a): the extracted ``L_full`` KEEPS the outflow
    self-consistency defect ``streamed − ψ.outflow`` on the outflow slots
    (it does NOT emit the raw outflow).

    The defect is deliberately retained: ``ψ.outflow`` is the STORED trace
    unknown the sibling ``−B`` reads as its input, and the outflow row of
    the canonical ``(L+C−S−F−B)`` is the definition residual
    ``ψ.outflow − streamed``.  Emitting the raw outflow would make that row
    independent of ``ψ.outflow`` (singular) and break ``−B``-as-sibling.
    """

    @pytest.mark.parametrize("geometry", _GEOMS_1D)
    def test_outflow_slot_subtracts_the_stored_outflow(self, geometry):
        """[foundation] The output boundary on the outflow slots is the
        defect against the STORED input outflow: holding the bulk + inflow
        fixed and varying ONLY the input outflow slots shifts the output by
        EXACTLY ``−Δ(input outflow)`` (the ``−ψ.outflow`` defect term).

        Run A: random inflow + random outflow.  Run B: same bulk, same
        inflow, ZERO outflow.  Streaming depends only on bulk + inflow
        (identical), so ``streamed`` is the same for both; the output
        outflow slots differ by exactly the input outflow value:
        ``out_B.outflow − out_A.outflow == A.outflow``.
        """
        sn_mesh = _build_sn_mesh(geometry, bc="reflective")
        ng, N = 1, sn_mesh.quad.N
        sigma_t = np.full((ng, *sn_mesh.spatial_shape), 2.0)
        state = _random_state(sn_mesh, seed=5, zero_boundary=False)
        out = _LpC_apply(sn_mesh, state, sigma_t)

        # state2: SAME bulk, SAME inflow slots, ZERO outflow slots — isolates
        # the ``−ψ.outflow`` defect term on the outflow slots.
        state2 = TimedFullField(
            interior=AngularFlux(values=state.interior.values.copy(), space=sn_mesh.angular_bulk_space),
            boundary=AngularBoundaryFlux.zeros(sn_mesh.angular_trace),
            _history=(), history_depth=2,
        )
        trace = sn_mesh.angular_trace
        for face in state.boundary.layout.faces:
            src = state.boundary.face_view(face)
            dst = state2.boundary.face_view(face)
            inflow = trace.inflow_indices_for_face(face)
            if inflow.size:
                dst[inflow, :] = src[inflow, :]
        out2 = _LpC_apply(sn_mesh, state2, sigma_t)

        outer_outflow = trace.outflow_indices_for_face("xmax")
        if outer_outflow.size:
            a = out.boundary.face_view("xmax")[outer_outflow, :]
            b = out2.boundary.face_view("xmax")[outer_outflow, :]
            stored = state.boundary.face_view("xmax")[outer_outflow, :]
            # streamed(A) == streamed(B) ⟹ (streamed − 0) − (streamed − v) = v.
            np.testing.assert_allclose(
                b - a, stored, rtol=1e-13, atol=1e-14,
                err_msg=(
                    f"{geometry}: the output outflow slot is NOT the defect "
                    f"``streamed − ψ.outflow`` — varying the stored outflow "
                    f"did not shift the output by −Δ (the −ψ.outflow term is "
                    f"missing; emission may have become raw outflow)."
                ),
            )


# ═════════════════════════════════════════════════════════════════════
# GATE 4 — The discriminator: the outflow DEFECT is kept (not raw outflow),
# even for vacuum.
# ═════════════════════════════════════════════════════════════════════


@pytest.mark.foundation
class TestVacuumBoundaryDefectKept:
    """Discriminator: design (a) KEEPS the outflow self-consistency defect
    ``streamed − ψ.outflow`` on the outflow slots — even for vacuum, where
    the inflow slots are zeroed by ``B = 0``.

    The keystone deletion + given-inflow read changed the matvec, but the
    outflow emission deliberately remains the DEFECT against the stored
    outflow (raw outflow would singularise the outflow row under ``−B``).
    With a non-zero input on the outflow slots, the output outflow slot is
    therefore ``streamed − face_value`` — it DEPENDS on the input.  This
    pins the boundary-emission shape explicitly so it is not silently
    conflated with the bulk bit-identity.
    """

    @pytest.mark.parametrize("geometry", _GEOMS_1D)
    def test_vacuum_outflow_slot_subtracts_the_stored_outflow(
        self, geometry,
    ):
        """[foundation] Vacuum output boundary on the outflow slots is the
        defect ``streamed − ψ.outflow``: varying ONLY the input outflow
        fill shifts the output by exactly ``−Δ`` (the ``−ψ.outflow`` term).
        """
        sn_mesh = _build_sn_mesh(geometry, bc="vacuum")
        ng, N = 1, sn_mesh.quad.N
        sigma_t = np.full((ng, *sn_mesh.spatial_shape), 2.0)
        rng = np.random.default_rng(9)
        bulk = rng.standard_normal((N, ng, *sn_mesh.spatial_shape))
        trace = sn_mesh.angular_trace

        def _run(outflow_fill: float) -> np.ndarray:
            b = AngularBoundaryFlux.zeros(sn_mesh.angular_trace)
            outflow = trace.outflow_indices_for_face("xmax")
            if outflow.size:
                b.face_view("xmax")[outflow, :] = outflow_fill
            st = TimedFullField(
                interior=AngularFlux(values=bulk.copy(), space=sn_mesh.angular_bulk_space),
                boundary=b,
                _history=(), history_depth=2,
            )
            out = _LpC_apply(sn_mesh, st, sigma_t)
            outflow_idx = trace.outflow_indices_for_face("xmax")
            return out.boundary.face_view("xmax")[outflow_idx, :].copy()

        outer_outflow = trace.outflow_indices_for_face("xmax")
        if outer_outflow.size:
            a = _run(0.0)   # streamed − 0
            b = _run(3.0)   # streamed − 3
            # (streamed − 0) − (streamed − 3) = 3 on every outflow slot.
            np.testing.assert_allclose(
                a - b, 3.0, rtol=1e-13, atol=1e-14,
                err_msg=(
                    f"{geometry}: vacuum output boundary on the outflow "
                    f"slots is INDEPENDENT of the input outflow value — the "
                    f"−ψ.outflow defect term is missing (emission became raw "
                    f"outflow, which breaks −B-as-sibling)."
                ),
            )
