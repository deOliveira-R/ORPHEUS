r"""#240 Phase 2 Step B — the removal-form ``InvertibleOperator(L(σ_t), C(σ_r))``
matvec ≡ sweep capability gate (σ_r ≠ σ_t).

SPEC STUB (test-architect, 2026-06-15); IMPLEMENTED (#240 Phase 2 Step B,
2026-06-15). The ``InvertibleOperator.apply`` / ``apply_transpose`` overrides
now realise the FULL within-group loss ``M(σ_C)ψ`` directly via
``loss_representation.loss_action(self.sigma, psi)``. The structural teeth gates
``test_invertible_apply_is_M_of_C_sigma_bit_identical`` and
``..._apply_transpose_...`` were ``xfail(strict=True)`` until the override
landed; the marker is REMOVED (they pass plainly). The value gates pass under
BOTH the pre-fix leak AND the override (see "THE PREMISE CORRECTION" below).

---

THE PREMISE CORRECTION (read before touching this file)
=======================================================

The driving brief framed the leak as a *numerical* defect: "matvec ≠ sweep⁻¹
for σ_r ≠ σ_t — a silent L21 violation". **This is FALSE as a numerical
claim, and the verification design pivots on that fact.**

Established by probe (``/tmp/probe_affine_robust.py``, this session) across
slab / sphere / cylinder / 2-D-Cartesian, vacuum + reflective, ≥2G het
non-flat ψ:

    leaky_apply(σ_t, σ_r) ≡ M(σ_r)ψ   to ≤ 2 ULP (rel ≤ 2.5e-16)

Mechanism: in the **forward (apply) direction** the WDD cell-balance matvec is
``M(σ)ψ = (denom·ψ_cell − numer_upstream)/V`` with ``denom = streaming_term +
σ·V``. So ``M(σ)ψ = streaming_action(ψ) + σ·ψ_cell`` — **affine in σ**.
Therefore the inherited ``OperatorSum.apply`` leak

    L.apply(ψ) + C.apply(ψ)
      = (M(σ_t)ψ − σ_t·ψ) + σ_r·ψ
      = streaming_action(ψ) + σ_r·ψ
      = M(σ_r)ψ                          (the override's value, to ULP)

is VALUE-CORRECT for the removal form. The non-affinity of ``M`` in σ (σ in the
``1/denom`` of the cell-average) lives in the SWEEP/SOLVE direction, NOT the
apply direction — so the round-trip ``op.apply(op.solve(q)) ≈ q`` PASSES under
the pre-fix leak too (probe RT1/RT2: ≤ 1.3e-15). **The proposed round-trip gate
has NO teeth against the leak** (it cannot distinguish leaky from override — a
vv Mode-9 degeneracy IN THE GATE DESIGN).

What the carve actually fixes is therefore a **principled-equivalence /
elegance** defect, not a correctness bug:

* ``coding-elegance`` Pattern 1+2: the composite ``(L+C)`` OWNS a single matvec
  ``loss_action(σ_C)``; realising it as ``L.apply + C.apply`` re-derives the
  streaming action through σ_t and then cancels it, a redundant round-trip and a
  twin-path realisation of one operator's action.
* The leak would become a genuine *latent* trap only if a future refactor made
  ``L.apply`` NON-affine in σ (e.g. an L-leaf that inverts a denominator). The
  override removes that latent coupling by construction (the composite never
  asks the leaf for a σ-bearing action it must then undo).

So the verification problem is the #240-Phase-2-Step-A shape (a value-preserving
re-association): prove (a) the override is the SAME value to nULP [value gate],
(b) the preserved value is CORRECT against a structurally-independent reference
[k_inf / Q-over-Σ_t ground], and (c) the override is STRUCTURALLY the composite's
own ``loss_action(σ_C)`` — NOT the leaf sum [the teeth, a structural gate, not a
numerical discriminator].

vv-principles compliance
========================

* **NOT a new ERR-NNN.** The leak is not a caught numerical bug (no wrong
  value ever shipped); it is a principled-equivalence/elegance carve. Per the
  "Log every caught bug" directive's own scope, ERR entries are for L0-caught
  *bugs*. This is a ``foundation`` software-contract gate, tagged
  ``verifies("loss-rep-resolution-a")`` (the Resolution-A / composite-owns-its-
  matvec convention it pins). NO ``catches(...)``.
* **vv Mode 9 (degeneracy):** the value gates run on the REMOVAL form
  (σ_r ≠ σ_t) — the regime the production-σ tests (σ_C==σ_t) cannot reach —
  AND on ≥2G het non-flat ψ with a genuine-``mu_y`` 2-D cubature. But note the
  Mode-9 degeneracy here is SUBTLER than usual: even on the non-degenerate
  removal regime the leaky and override values COINCIDE (affine-in-σ), so the
  value gate is a same-value cross-check, and the TEETH live in the structural
  gate (c).
* **Cardinal Rule (1G degenerate):** every value row is ≥2G. A 1G removal row
  is included ONLY as the affine-in-σ characterisation (no eigenvalue claim).
* **vv Mode 8 (-O strips assert):** ``np.testing`` / ``pytest.fail`` only;
  fires under the canonical ``python -O``.

Run (HOST → ``.venv/bin/python``; canonical ``-O``)::

    .venv/bin/python -O -m pytest \
      tests/sn/operators/test_removal_form_matvec_sweep.py \
      -p no:cacheprovider -q

Route-arounds for the wider operators sweep (the 7 pre-existing reds —
#195/#209 SPH, #214 mu_y; NEVER all tests/sn, #212 hang)::

    -k "not (vacuum_bulk_bit_identical_1d and SPH) \
        and not (sphere_1g_apply_bit_identical or sphere_2g_apply_bit_identical) \
        and not test_2d_mesh_resolution and not two_d_cartesian_loss_action"
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.geometry import (
    BC, CoordSystem, Mesh1D, Mesh2D, Region, RegionMesh, StructuredGeometry,
)
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.geometry import SNMesh
from orpheus.sn.boundary_operator import SNBoundaryOperator
from orpheus.sn.operator import (
    CollisionOperator,
    InvertibleOperator,
    StreamingOperator,
)
from orpheus.transport.fields.angular_flux import AngularFlux
from orpheus.transport.fields.boundary_flux import BoundaryFlux
from orpheus.transport.timed_full_field import TimedFullField
from tests.sn._test_helpers import placeholder_materials

pytestmark = [
    pytest.mark.foundation,
    # The convention this file pins: the composite (L+C) owns its matvec as
    # the FULL loss M(σ_C); the operator's only glue is Resolution-A. The
    # removal form (σ_r ≠ σ_t) is the regime that makes "which σ does the
    # matvec realise?" observable as a STRUCTURAL fact.
    pytest.mark.verifies("loss-rep-resolution-a"),
]

# #240 Phase 2 Step B (2026-06-15) landed the InvertibleOperator.apply /
# apply_transpose overrides (loss_action(self.sigma)); the strict-xfail-
# until-override marker the teeth gates carried is REMOVED — they pass plainly.


# ── Geometry builders ───────────────────────────────────────────────────


def _slab(nx: int = 6, n_ord: int = 4, ng: int = 2, bc: str = "vacuum") -> SNMesh:
    geom = StructuredGeometry(
        geometry="SLB",
        regions=(Region(mat_id=0, outer_thickness_cm=2.0),),
        bcs=(BC(bc), BC(bc)),
    )
    mesh = Mesh1D.from_geometry(geom, region_meshes=(RegionMesh(n_cells=nx),))
    quad = Quadrature.gauss_legendre(n_ordinates=n_ord)
    return SNMesh(mesh, quad, placeholder_materials(ng=ng))


def _sphere(nx: int = 6, ng: int = 2, bc: str = "vacuum") -> SNMesh:
    geom = StructuredGeometry(
        geometry="SPH",
        regions=(Region(mat_id=0, outer_thickness_cm=2.0),),
        bcs=(BC(bc),),
    )
    mesh = Mesh1D.from_geometry(geom, region_meshes=(RegionMesh(n_cells=nx),))
    quad = Quadrature.level_symmetric(sn_order=4)
    return SNMesh(mesh, quad, placeholder_materials(ng=ng))


def _cyl(nx: int = 6, ng: int = 2, bc: str = "vacuum") -> SNMesh:
    geom = StructuredGeometry(
        geometry="CYL",
        regions=(Region(mat_id=0, outer_thickness_cm=2.0),),
        bcs=(BC(bc),),
    )
    mesh = Mesh1D.from_geometry(geom, region_meshes=(RegionMesh(n_cells=nx),))
    quad = Quadrature.level_symmetric(sn_order=4)
    return SNMesh(mesh, quad, placeholder_materials(ng=ng))


def _cart2d(nx: int = 4, ny: int = 5, ng: int = 2, bc: str = "reflective") -> SNMesh:
    """NON-SQUARE 2-D Cartesian, ``level_symmetric`` (genuine mu_y — avoids the
    #214 ``mu_y==0`` GL rank-mismatch). Non-square is the x↔y-swap catcher."""
    mesh = Mesh2D(
        edges_x=np.linspace(0.0, 2.0, nx + 1),
        edges_y=np.linspace(0.0, 2.0, ny + 1),
        mat_map=np.zeros((nx, ny), dtype=int),
        coord=CoordSystem.CARTESIAN,
        bc_xmin=BC(bc), bc_xmax=BC(bc), bc_ymin=BC(bc), bc_ymax=BC(bc),
    )
    return SNMesh(mesh, Quadrature.level_symmetric(sn_order=4), placeholder_materials(ng=ng))


# Removal-form value gates: slab+sphere+cyl+2D, all ≥2G, vacuum (single solve,
# no reflective iteration — the streaming-anisotropic Mode-9 regime).
_REMOVAL_CASES = {
    "slab_2g": lambda: _slab(ng=2, bc="vacuum"),
    "sphere_2g": lambda: _sphere(ng=2, bc="vacuum"),
    "cyl_2g": lambda: _cyl(ng=2, bc="vacuum"),
    "cart2d_2g": lambda: _cart2d(ng=2, bc="reflective"),
}


def _removal_sigmas(sn: SNMesh, *, seed: int) -> tuple[np.ndarray, np.ndarray]:
    r"""Heterogeneous σ_t and a TRUE removal σ_r = σ_t − Σ_s0 > 0.

    σ_r is built as a numpy array directly (independent of the mesh's
    ``materials``, which carry no scattering) — the operator algebra
    ``InvertibleOperator(StreamingOperator(σ_t), CollisionOperator(σ_r))``
    consumes the arrays, so the removal form is representable on EVERY
    geometry (slab, sphere, cylinder, 2-D) without a scattering mixture.
    ``__init__`` requires ``σ_r > 0`` (operator.py:1244) — guaranteed by the
    bounded self-scatter draw.
    """
    rng = np.random.default_rng([seed, 240])
    sig_t = rng.uniform(1.0, 3.0, size=(sn.ng, *sn.spatial_shape))
    # Self-scatter Σ_s0 strictly below σ_t so σ_r = σ_t − Σ_s0 > 0 everywhere.
    sig_s0 = rng.uniform(0.2, 0.8, size=(sn.ng, *sn.spatial_shape))
    sig_r = sig_t - sig_s0
    # Fixture self-checks via pytest.fail (NOT bare assert) so a malformed σ_r
    # fixture fails LOUDLY under the canonical -O (vv Mode 8), not silently.
    if not np.all(sig_r > 0.0):
        pytest.fail("removal σ_r must stay positive (test-fixture bug)")
    # Guard the regime is genuinely NON-degenerate: σ_r ≠ σ_t materially.
    if not float(np.min(np.abs(sig_t - sig_r))) > 0.1:
        pytest.fail("σ_r ≈ σ_t — degenerate fixture (must break σ_C==σ_t)")
    return sig_t, sig_r


def _random_state(sn: SNMesh, *, seed: int) -> TimedFullField:
    rng = np.random.default_rng([seed, 7])
    state = TimedFullField.zeros(bulk=AngularFlux, boundary=BoundaryFlux, mesh=sn)
    state.bulk.values[...] = rng.standard_normal(state.bulk.values.shape)
    for face in state.boundary.layout.faces:
        fv = state.boundary.face_view(face)
        fv[...] = rng.standard_normal(fv.shape)
    return state


# ═══════════════════════════════════════════════════════════════════════
# (a) THE STRUCTURAL TEETH — the override realises the composite's OWN
#     loss_action(σ_C), NOT the leaf sum L.apply + C.apply.
#     These are the gates that DISTINGUISH the override from the leak. They were
#     xfail(strict=True) until the override landed; #240 Phase 2 Step B
#     (2026-06-15) wired InvertibleOperator.apply = loss_action(self.sigma), so
#     the marker is REMOVED and these now register as plain pass.
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.parametrize("case", list(_REMOVAL_CASES))
def test_invertible_apply_is_M_of_C_sigma_bit_identical(case):
    r"""[L0 structural] ``(L+C).apply(ψ) == M(σ_C)ψ`` BIT-IDENTICALLY (the teeth).

    THE discriminator. The override MUST realise the composite's matvec as the
    FULL loss ``M(σ_C)ψ`` directly. The structurally-independent reference for
    ``M(σ_C)ψ`` (σ_C = σ_r, the diagonal's σ) is a SEPARATE
    ``StreamingOperator`` whose OWN σ_t IS σ_r — its ``loss_action`` is
    unambiguously ``M(σ_r)ψ`` (no leak / removal ambiguity; the leaf reads its
    own σ). Post-#240-Step-B the reference is built with the migrated
    ``loss_action(σ_r, psi)`` signature (σ passed explicitly) — only
    ``op.apply`` drives the discriminator.

    * Under the OVERRIDE, ``op.apply`` IS ``loss_action(σ_r)`` on the SAME
      walk → ``array_equal`` (0 ULP) → xpass.
    * Under the LEAK, ``op.apply`` = ``L.apply(σ_t) + C.apply(σ_r)`` = the leaf
      sum, value-EQUAL to ``M(σ_r)`` (affine-in-σ) but a DIFFERENT reduction
      tree → ≤ 2 ULP off ``array_equal`` → AssertionError → strict-xfail
      registers (proves the gate has teeth, on slab/sphere where the leaf
      diamond walk and the direct walk re-associate; cyl/2D too — probe shows
      2/192, 26/768 ULP-level differences).

    This is the ONE gate that distinguishes the carve. Everything else is
    consistency / value-ground (passes under BOTH states).
    """
    sn = _REMOVAL_CASES[case]()
    sig_t, sig_r = _removal_sigmas(sn, seed=sum(map(ord, case)))
    L = StreamingOperator(sn)
    C_r = CollisionOperator(sn, sig_r)
    op = InvertibleOperator(L, C_r)
    psi = _random_state(sn, seed=sum(map(ord, case)))

    composite_matvec = op.apply(psi).bulk.values
    # M(σ_r)ψ — independent leaf whose OWN σ_t IS σ_r. Post-#240-Step-B
    # signature: loss_action takes the σ array directly (here σ_r), not the
    # operator. The reference's σ IS σ_r, so it is unambiguously M(σ_r)ψ.
    L_ref = StreamingOperator(sn)
    M_sigma_r = L_ref.loss_representation.loss_action(sig_r, psi).bulk.values

    np.testing.assert_array_equal(
        composite_matvec, M_sigma_r,
        err_msg=(
            f"[{case}] (L+C).apply.bulk != M(σ_r)ψ BIT-IDENTICAL — the "
            "composite matvec is NOT its own loss_action(σ_r); it is still the "
            "inherited leaf sum L.apply(σ_t) + C.apply(σ_r) (the #240 Step B "
            "leak). Value-EQUAL (affine-in-σ) but ULP-distinct: the override "
            "must own the matvec directly via loss_action(self.sigma)."
        ),
    )


@pytest.mark.parametrize("case", ["slab_2g", "sphere_2g", "cyl_2g"])  # NOT cart2d
def test_invertible_apply_transpose_is_M_transpose_of_C_sigma_bit_identical(case):
    r"""[L0 structural] ``(L+C).apply_transpose(φ) == M(σ_C)ᵀφ`` bit-identical (adjoint teeth).

    The adjoint sibling of the teeth gate. The override
    ``InvertibleOperator.apply_transpose`` MUST realise ``M(σ_C)ᵀφ`` directly,
    bit-identical to a SEPARATE ``StreamingOperator(σ_r)``'s
    ``loss_action_transpose``. ``array_equal`` so the ≤2-ULP leaf-sum
    realisation fails (the strict-xfail teeth).

    2-D Cartesian EXCLUDED: ``ScanMarch.loss_action_transpose`` is a deferral
    raise (multi-D Cartesian adjoint not implemented — the override must not
    silently route around it). The adjoint teeth are slab/sphere/cyl only; the
    2-D adjoint is gated by the deferral-raise contract elsewhere.
    """
    sn = _REMOVAL_CASES[case]()
    sig_t, sig_r = _removal_sigmas(sn, seed=sum(map(ord, case)) + 1)
    L = StreamingOperator(sn)
    op = InvertibleOperator(L, CollisionOperator(sn, sig_r))
    if "apply_transpose" not in op.capabilities:
        pytest.fail(f"[{case}] apply_transpose not advertised on InvertibleOperator.")
    phi = _random_state(sn, seed=sum(map(ord, case)) + 1)

    composite_t = op.apply_transpose(phi).bulk.values
    # Post-#240-Step-B signature: loss_action_transpose takes the σ array (σ_r).
    L_ref = StreamingOperator(sn)
    M_t_sigma_r = L_ref.loss_representation.loss_action_transpose(
        sig_r, phi,
    ).bulk.values

    np.testing.assert_array_equal(
        composite_t, M_t_sigma_r,
        err_msg=(
            f"[{case}] (L+C).apply_transpose.bulk != M(σ_r)ᵀφ bit-identical — "
            "the composite adjoint is still the leaf sum, not its own action."
        ),
    )


# ═══════════════════════════════════════════════════════════════════════
# (b) THE VALUE GROUND — the removal-form matvec ≡ sweep round-trip and the
#     structurally-independent reference. These PASS under BOTH the leak and
#     the override (the leak is value-correct — affine-in-σ). They prove the
#     PRESERVED value is CORRECT; they do NOT discriminate the carve (that is
#     gate (a)). NOT xfail.
# ═══════════════════════════════════════════════════════════════════════


# Sphere EXCLUDED from the bare round-trip: the curvilinear M-M sweep reads the
# PREVIOUS ITERATE for the Carlson coupled-pole closure (ERR-058 family), so a
# single ``op.solve`` is NOT the one-shot inverse of ``op.apply`` — the bare
# operator round-trip DIVERGES on sphere (probe: 9.8e5 single, 2e99 iterated).
# Cylinder round-trips cleanly (its degenerate pure-azimuthal ordinate does not
# break the one-shot inverse — probe 3.1e-15). The sphere matvec≡sweep CLAIM is
# carried by the structural teeth gate (a) ``apply == M(σ_r)`` (which does NOT
# round-trip) + the existing fixed-point bridge ``TestInvertibleSolveBridgeRegression``
# (production σ). DO NOT add sphere here without the warm-started fixed-point loop.
_ROUNDTRIP_CASES = ["slab_2g", "cyl_2g", "cart2d_2g"]


@pytest.mark.parametrize("case", _ROUNDTRIP_CASES)
def test_removal_form_matvec_sweep_roundtrip(case):
    r"""[L0] ``op.apply(op.solve(q)) ≈ q`` AND ``op.solve(op.apply(ψ)) ≈ ψ`` for σ_r ≠ σ_t.

    The matvec-and-sweep-are-inverse-twins gate ON THE REMOVAL FORM
    (slab / cylinder / 2-D-Cartesian — sphere excluded, see ``_ROUNDTRIP_CASES``).
    Both directions to solver tolerance.

    HONEST SCOPE (the premise correction): this round-trip passes under the
    PRE-FIX leak too — the leaky ``apply`` is value-correct (affine in σ in the
    forward direction). So this gate is NOT the teeth; it is the consistency
    ground that the removal-form solve/apply ARE mutually inverse at all (a
    real property the carve must preserve), and that the override does not
    BREAK the round-trip. The teeth are gate (a).

    Vacuum (slab/cyl) ⟹ a single ``solve`` is the exact (lower-triangular)
    inverse — no reflective iteration. ``op.solve`` writes the outflow trace
    into ``.boundary``; ``op.apply`` reads it back as the inflow context, so the
    round-trip closes in one solve + one apply. (2-D reflective: the same
    one-shot inverse holds for the bare operator with its captured trace — probe
    1.8e-15 — because Cartesian has no iterate-dependent pole seed.)
    """
    sn = _REMOVAL_CASES[case]()
    sig_t, sig_r = _removal_sigmas(sn, seed=sum(map(ord, case)) + 2)
    op = InvertibleOperator(
        StreamingOperator(sn), CollisionOperator(sn, sig_r),
    )

    # Direction 1: apply ∘ solve = identity on a volumetric source.
    q = _random_state(sn, seed=sum(map(ord, case)) + 2)
    # Vacuum: no boundary inflow — zero the trace so q is purely volumetric.
    for face in q.boundary.layout.faces:
        q.boundary.face_view(face)[...] = 0.0
    psi = op.solve(q)
    q_back = op.apply(psi)
    np.testing.assert_allclose(
        q_back.bulk.values, q.bulk.values, rtol=1e-10, atol=1e-12,
        err_msg=(
            f"[{case}] op.apply(op.solve(q)) != q on the removal form — the "
            "matvec and sweep are NOT inverse twins for σ_r ≠ σ_t."
        ),
    )

    # Direction 2: solve ∘ apply = identity on an arbitrary ψ (with its trace).
    psi0 = _random_state(sn, seed=sum(map(ord, case)) + 3)
    # #257 S8a — the matvec leaf is a base arrow: ``op.apply`` returns a
    # timeless FullField source.  ``solve`` consumes the driver's timed rhs, so
    # wrap the source back into a TimedFullField (what the SI / Krylov driver
    # does — it carries the comonad; the operator does not).  Byte-identical.
    Lpsi_source = op.apply(psi0)
    Lpsi = TimedFullField(
        bulk=Lpsi_source.bulk,
        boundary=Lpsi_source.boundary,
        _history=(),
        history_depth=psi0.history_depth,
    )
    psi_back = op.solve(Lpsi)
    np.testing.assert_allclose(
        psi_back.bulk.values, psi0.bulk.values, rtol=1e-10, atol=1e-12,
        err_msg=(
            f"[{case}] op.solve(op.apply(ψ)) != ψ on the removal form."
        ),
    )


@pytest.mark.parametrize("case", list(_REMOVAL_CASES))
def test_removal_form_apply_value_equals_M_of_sigma_r(case):
    r"""[L0] The removal-form matvec VALUE equals ``M(σ_r)ψ`` (the affine-in-σ fact).

    Structurally-independent reference for the matvec VALUE: build a SEPARATE
    ``StreamingOperator(σ_r)`` (σ_r AS its σ_t) and take its ``loss_action`` —
    this is an UNAMBIGUOUS ``M(σ_r)ψ`` (no removal/leak ambiguity, since the
    leaf's own σ IS σ_r). The composite's ``apply`` must equal it to nULP
    REGARDLESS of leak-vs-override.

    Pillar: this is the affine-in-σ characterisation — the property that makes
    the leak value-correct. ``assert_allclose`` at nULP (NOT array_equal): the
    leaf sum and ``M(σ_r)`` re-associate ≤ 2 ULP; the override is the SAME
    ``loss_action`` call so it is bit-identical to this reference (0 ULP) but
    nULP is the regime-agnostic contract. 1G companion included as
    characterisation ONLY (no eigenvalue claim — vv H1).
    """
    sn = _REMOVAL_CASES[case]()
    sig_t, sig_r = _removal_sigmas(sn, seed=sum(map(ord, case)) + 4)
    op = InvertibleOperator(
        StreamingOperator(sn), CollisionOperator(sn, sig_r),
    )
    psi = _random_state(sn, seed=sum(map(ord, case)) + 4)

    leaky_or_override = op.apply(psi).bulk.values
    # Independent M(σ_r): a leaf whose OWN σ_t IS σ_r — unambiguously M(σ_r)ψ.
    # Post-#240-Step-B signature loss_action(sigma, psi) — pass σ_r directly.
    # Value-ground (not the teeth): holds under leak AND override.
    L_ref = StreamingOperator(sn)
    M_sigma_r = L_ref.loss_representation.loss_action(sig_r, psi).bulk.values

    np.testing.assert_allclose(
        leaky_or_override, M_sigma_r, rtol=0.0, atol=1e-12,
        err_msg=(
            f"[{case}] (L+C).apply != M(σ_r)ψ — the removal-form matvec does "
            "not realise the FULL loss at σ_r (a genuine σ-routing bug, "
            "distinct from the leak which is value-correct)."
        ),
    )


def test_removal_form_kinf_independent_reference_2g():
    r"""[L1] Structurally-independent eigenvalue ground for the removal-form sweep.

    The value gates above prove the removal-form matvec/sweep are CONSISTENT
    (inverse twins) and value-correct vs ``M(σ_r)``. Per vv-principles
    criterion 2, consistency + old-vs-new ULP is necessary but NOT sufficient —
    the PRESERVED value must be cross-checked against a structurally-independent
    reference. Here: a homogeneous reflective medium where the removal form
    folds the within-group self-scatter into the diagonal must recover the
    closed-form ``k_inf = νΣ_f / Σ_a`` (mesh-independent, ≥2G — NOT 1G, vv
    Cardinal Rule).

    STUB: requires a fissile 2G mixture whose within-group self-scatter is
    folded into σ_r and whose remaining (down/up)-scatter feeds S. This wires
    through ``solve_sn`` with the removal-form inner operator. Marked xfail
    until the removal-form solver entry exists (issue #200 — "fold within-group
    self-scatter into the diagonal collision term"); the OPERATOR-level
    round-trip (above) is the landable Step-B ground, the SOLVER-level k_inf
    is the #200 follow-on. Pillar: closed-form.
    """
    pytest.xfail(
        "removal-form k_inf recovery needs the #200 solver entry (fold "
        "within-group self-scatter into the diagonal). The Step-B operator "
        "round-trip + M(σ_r) value gate are the landable grounds; this "
        "eigenvalue cross-check lands with #200."
    )


# ═══════════════════════════════════════════════════════════════════════
# (c) THE PRODUCTION-σ INVARIANT — the override MUST NOT change the production
#     path (σ_C == σ_t). nULP-equal to the inherited leaf sum on every
#     geometry (FP-re-association of the round-trip the override removes; measured
#     ≤2 ULP slab/sphere, ≤4 cyl, ≤5 2-D). Pins the value-preserving claim.
#
#     RE-BASELINE (#240 Step B impl, 2026-06-15): the test-architect spec marked
#     slab/sphere ``bitexact=True`` from a probe that compared two override-style
#     ``loss_action`` paths (0/32). The ACTUAL production gate compares the
#     override (``loss_action(σ_t)`` direct) against the inherited leaf sum
#     (``(loss_action(σ_t) − σ_t·ψ) + σ_t·ψ``): the override REMOVES that
#     ``−σ_t·ψ + σ_t·ψ`` round-trip, which FP-re-associates by 1–2 ULP even on
#     the slab/sphere _OneDimScanWalk (re-probed: slab 1/48 @ 2 ULP, sphere 1/288
#     @ 2 ULP, cyl 16/288 @ 4 ULP, 2-D 63/960 @ 5 ULP). The override is the MORE
#     accurate path; the re-baseline to nULP is principled per ``vv-principles``
#     §"Bit-identity vs principled-equivalence" (named intermediate ``loss_action``
#     / verified against M(σ_r) by the teeth gate (a) / drift = reduction-depth ×
#     ULP). Bit-identity is NOT a property of the override-vs-leaf-sum pair.
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.parametrize(
    "case",
    ["slab_2g", "sphere_2g", "cyl_2g", "cart2d_2g"],
)
def test_production_sigma_apply_value_preserved(case):
    r"""[L0] Production σ_C == σ_t: override apply ≈ inherited-leaf-sum value (nULP).

    The value-preserving guard. With σ_C == σ_t the override
    (``loss_action(σ_t)`` direct) and the inherited leaf sum
    (``L.apply + C.apply``) compute the SAME ``M(σ_t)ψ`` to nULP.  #257 S8b:
    ``L.apply`` is now pure σ-free ``streaming_action(ψ) = loss_action(0, ψ)``,
    so the leaf sum is ``loss_action(0, ψ) + σ_t⊙ψ`` while the override is
    ``loss_action(σ_t, ψ)`` directly — the affine relation
    ``M(σ_t)ψ = streaming_action(ψ) + σ_t⊙ψ`` makes them value-equal, but the
    σ-free walk re-associates the FP reduction tree relative to the σ-bearing
    matvec, so the agreement is to FP-non-associativity (rel ≤ 1e-14, the
    per-element ULP metric spiking at near-zero cancellation values:
    measured ≤ 33 on cylinder, ≤ 1024 on 2-D Cartesian).  The bound is
    widened to ``nulp=2048`` to absorb the near-zero spikes; the
    ``rel < 1e-14`` guard below is the genuine value-preservation ground.

    This builds BOTH the production composite (σ_t == σ_t) AND an explicit
    ``L.apply(state) + C.apply(state)`` reference (the inherited leaf-sum value,
    computed via the free expression so it does NOT route through the override).
    POST-S8b the override (``loss_action(σ_t)``) and the pure-L leaf sum
    (``loss_action(0) + σ_t⊙ψ``) agree to FP — so this gate does NOT
    discriminate the carve (gate (a) does); it is the value-PRESERVATION pin
    that the override did not ALGORITHMICALLY change production.
    """
    sn = _REMOVAL_CASES[case]()
    # Production: σ_C == σ_t (the same array on both leaves).
    rng = np.random.default_rng([sum(map(ord, case)), 1])
    sig_t = rng.uniform(0.5, 3.0, size=(sn.ng, *sn.spatial_shape))
    L = StreamingOperator(sn)
    C = CollisionOperator(sn, sig_t)        # σ_C == σ_t (production)
    op = InvertibleOperator(L, C)
    state = _random_state(sn, seed=sum(map(ord, case)))

    composite = op.apply(state).bulk.values
    leaf_sum = (L.apply(state).bulk.values + C.apply(state).bulk.values)

    np.testing.assert_array_almost_equal_nulp(
        composite, leaf_sum, nulp=2048,
    )
    # And the drift is genuinely tiny (FP re-association, not algorithmic).
    # pytest.fail (NOT bare assert) so the bound fires under -O (vv Mode 8).
    rel = np.abs(composite - leaf_sum).max() / max(
        np.abs(leaf_sum).max(), 1e-300,
    )
    if not rel < 1e-14:
        pytest.fail(
            f"[{case}] production-σ override drift {rel:.2e} exceeds the "
            "FP-re-association bound — an algorithmic change masquerading "
                "as ULP noise (vv criterion 3)."
            )


# ═══════════════════════════════════════════════════════════════════════
# (d) NEGATIVE CONTROL — the σ_r > 0 construction guard. A removal form with
#     σ_r dipping non-positive MUST be rejected at construction (operator.py:
#     1244), so the gate's fixtures cannot silently pass on a malformed σ_r.
# ═══════════════════════════════════════════════════════════════════════


def test_removal_form_nonpositive_sigma_r_rejected():
    r"""[L0 negative] ``InvertibleOperator(L, C(σ_r ≤ 0))`` raises at construction.

    The positive control of the gate's own fixture: if a future cross-section
    set drives σ_r = σ_t − Σ_s0 non-positive, the WDD sweep emits NaN at those
    cells; surfacing it at construction (vs a silent NaN matvec) is the
    illegal-states-unrepresentable contract (``coding-elegance`` Pattern 4).
    Pairs with the positive removal cases above (vv L11: a contract-validation
    method needs BOTH a positive (must-not-raise — the removal cases) AND a
    negative (must-raise — here) test).
    """
    sn = _slab(ng=2)
    sig_t = np.full((sn.ng, *sn.spatial_shape), 1.0)
    sig_r = sig_t.copy()
    sig_r.flat[0] = -0.1   # one cell over-scatters → σ_r < 0
    L = StreamingOperator(sn)
    C_bad = CollisionOperator(sn, sig_r)
    with pytest.raises(ValueError, match="strictly positive"):
        InvertibleOperator(L, C_bad)
