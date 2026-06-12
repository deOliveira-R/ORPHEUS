---
name: s6-relayering-verification
description: S6 operator/representation re-layering V&V plan (capstone of the sweep-strategy carve, #222). SweepStrategy→SpatialRepresentation; residual→loss_action returns (L+C)ψ (operator does −C glue). Bit-identity anchors per stage; A2D-1 regenerate-in-commit-pin-output protocol; the "one representation instance" IS-discriminating test (concrete spec); curvilinear angular-transpose hazard (g_adjoint_reciprocity); as_scipy_linop Krylov round-trip; Mode-8/Mode-9 discipline.
metadata:
  type: project
---

# S6 — Operator / representation re-layering verification plan

S6 is the capstone of the SweepStrategy carve (issue #222). Extends the S0 parent
[[sweep-strategy-carve-verification]] and the S5 [[scan-march-verification]] — it
REUSES their L2/L3/L5/L6/L7/L8/L9 anchor vocabulary and adds the structural
"one-instance" payoff gate. Worktree `worktree-sn-nd-layout`, HEAD `716a0c1`
(verified this session). Authoritative design = `.claude/plans/sn_sweep_strategy.md`
§S6 (lines ~475–679, LOCKED, cross-domain-attacker-validated across 4 expert frames).
This is the proactive `test-architect` gate (lesson L17) before the main agent writes
S6 turn-by-turn — the carve crosses the operator↔representation subsystem boundary.

## EXECUTIVE SUMMARY (5 lines)
1. **Anchor set (all VERIFIED EXISTS @ 716a0c1):** A2D-1 source-hash, `window≡full`
   oracle (`test_2d_full_field_oracle.py` + `test_sweep_graph_window_equivalence.py`),
   `test_affine_carve_bit_identity` (3 sha256-golden cases, `-O`-safe `raise`),
   `SI≡Krylov≡k_inf` (`test_keff_2d`), ScanMarch G2.c (`test_scan_march_equivalence.py`).
2. **The ONE discriminating test (S6.5 payoff, concrete spec §4):** `assert`
   `L.sweep_strategy IS the SpatialRepresentation that InvertibleOperator.solve uses` —
   FAILS today (two doors → two instances), PASSES after S6.5. `-O`-safe via `pytest.fail`.
3. **A2D-1 protocol (§2):** S6.3 WILL regenerate the source hash (the `_apply_2d_cartesian`
   BODY moves into `loss_action`). Refresh `EXPECTED_SHA256` IN THE SAME COMMIT; assert
   OUTPUT byte-identity via `window≡full` oracle — NEVER pin a stale source hash.
4. **loss_action convention pin (§3):** `representation.loss_action(op, ψ)` returns
   `(L+C)ψ` (NOT `L·ψ`); `op.apply(ψ) == loss_action(...) − σ_t·ψ.bulk` (the `−C` glue,
   verified ONCE). Pure-z octant + zero-source probe are the cross-checks.
5. **Per-stage gate table = §8.** S6.2 rename (bit-id) · S6.3 walk-move (A2D-1 regen +
   curvilinear-transpose hazard §5) · S6.4 `_OctantWalk2D` (bit-id) · S6.5 unify doors
   (one-instance test turns green) · S6.6 deferred. S6.7 `as_scipy_linop` Krylov gate §6.

## Claim-layer / pillar gate (vv-principles §1.5 — MUST pass before writing)
Every S6 gate is a **behavior-PRESERVING carve gate**, not a new-physics claim. The
claim layers and pillars:
- **Bit-identity anchors** (A2D-1, affine-carve golden, `window≡full`): NOT a pillar
  claim — they are *free verification by inheritance* from already-verified references.
  Convention-order claim is irrelevant; the inheritance is exact (`np.array_equal`/sha256).
- **`SI≡Krylov≡k_inf` (G6):** EIGENVALUE claim → **closed-form pillar** (k_inf =
  λ_max(A⁻¹F) transfer-matrix, structurally independent). Correct: MMS is NOT used for
  the eigenvalue leg (MMS proves zero eigenvalue information). ≥2G (cardinal rule — no
  1G eigenvalue test). PASS.
- **`loss_action ≡ FullFieldWavefront.loss_action` (G2.c, the renamed ScanMarch gate):**
  flux-shape claim → **`FullFieldWavefront` oracle**, transitively pinned to k_inf (L8)
  + φ=Q/Σ_t (L7). Structurally independent (two topological linearizations of the SAME
  lower-triangular solve). PASS.
- **Curvilinear angular-transpose (§5):** the reciprocity identity `⟨Aψ,φ⟩_G=⟨ψ,A.Hφ⟩_G`
  → **closed-form pillar** (the defining adjoint property, with an INDEPENDENT g_inner
  metric — anti-R1). PASS.
Structural-independence terminates in: closed-form k_inf transfer matrix, the
`FullFieldWavefront` d-generic spine (different schedule, same operator), and the
independent `_g_inner` metric — NONE is another procedurally-different derivation of the
same identity. Gate PASSES. Matrix ready to write.

---

## 0. What S6 IS (the behaviour the gates target) — and what stays invariant

S6 re-layers WITHOUT changing the converged answer. The native frame (S6.0):
`(L+C)` is a lower-triangular `LinearOperator`; the four walks (`CumprodScan`,
`MovingFrontierWindow`, `ScanMarch`, `FullFieldWavefront`) are four *representations*
of ONE operator. The two structural defects S6 closes (both `coding-elegance` Smell #16):

1. **Half-done dispatch** — `operator.apply(ψ) → strategy.residual(op,ψ) → op._apply_<x>(ψ)`;
   the matvec impl lives back ON the operator (`_apply_1d`, `_apply_2d_cartesian`,
   `_apply_2d_cartesian_scanmarch`, `_apply_full_field`). Dispatch was relocated, not eliminated.
2. **TWO DOORS** — VERIFIED THIS SESSION at `operator.py`:
   - `apply` → `self.sweep_strategy` (a `@cached_property`, `operator.py:1566-1584` →
     `default_for(self.sn_mesh)`).
   - `solve` → `_solve_timed_full_field` → `transport_sweep` (`operator.py:2530`) →
     `default_for(sn_mesh)` (`sweep.py:239-241`) — a SECOND, independent `default_for`
     call. So L21 ("matvec ≡ sweep") is a TEST-MAINTAINED coincidence, not a construction
     theorem. The two `default_for(...)` calls return TWO instances today (§4 is the test).

**The convention that S6.3 fixes in place (confirmed THIS SESSION — pin it, §3):**
`loss_action` returns `(L+C)ψ`, the operator subtracts `C`. Three code facts re-verified:
- (a) `residual_kernel_batch` natively yields `(L+C)ψ` at zero source (`diamond.py:~363`).
- (b) the `−C` glue is LIVE at `operator.py:1483-1487` (`_apply_1d`):
  `cell_values = LpC_result.bulk.values − self.sigma_t[None] * psi.bulk.values`. The
  transpose twin mirrors it at `operator.py:1552-1556` (`_apply_1d_transpose`). This is
  the `L = (L+C) − C` Resolution-A subtraction that S6.3 keeps as the ONLY operator glue.
- (c) the pure-z octant sets `LpC = σ_t·ψ` so `L·ψ=0` after subtraction (`operator.py:~1731`).

---

## 1. Bit-identity inheritance anchors (pins S6 as behavior-preserving)

All paths VERIFIED to exist at HEAD `716a0c1`. Classify each: **bit-identity** (free
inheritance), **value-ground** (structurally-independent reference), **principled-equiv**
(nulp). For each: which S6 stage it gates + byte-identical-vs-principled call.

| Anchor | File::test (worktree-relative) | Class | Gates stage | Stays |
| ------ | ------------------------------ | ----- | ----------- | ----- |
| **A2D-1 source-hash** | `tests/sn/operators/test_streaming_operator.py::TestT4dApply2DCartesianSourceHashPin::test_apply_2d_cartesian_source_hash_unchanged` (`EXPECTED_SHA256` of `inspect.getsource(StreamingOperator._apply_2d_cartesian)`, `operator.py:1306-1310`) | bit-identity (source text) | **S6.3** (body moves) | **REGENERATES at S6.3** — see §2. Stays byte-free through S6.2/S6.4. |
| **affine-carve golden** | `tests/sn/solve/test_affine_carve_bit_identity.py::test_converged_flux_bit_identical_after_affine_carve` (3 params: `si_2d_p1_aniso_het` / `krylov_2d_p1_aniso_het` / `si_slab_2g_het`; sha256 of converged `bulk.values` + `phi`; `-O`-safe `raise AssertionError`) | bit-identity (output bytes, end-to-end) | **ALL S6 stages** | **STAYS byte-identical** — S6 does NOT change FP-association of the default path (Fork B1: window stays the d=2 default, scan stays d=1). If S6 ever flips a default, REGENERATE per the scan-march Fork-B2 discipline ([[scan-march-verification]] §G5). For S6.2–S6.5 as designed: byte-identical. |
| **window≡full SWEEP** | `tests/sn/sweep/cartesian_2d/test_2d_full_field_oracle.py::test_sweep_window_equals_full_field_end_to_end` (4 cases incl. NON-SQUARE 12×7, 5×9; `np.testing.assert_array_equal`) | bit-identity | **S6.3/S6.4** (the OUTPUT-identity oracle the A2D-1 regen leans on) | **STAYS** — pins `MovingFrontierWindow.loss_action ≡ FullFieldWavefront.loss_action` after the bodies move into representations. |
| **window≡full MATVEC** | `tests/sn/sweep/cartesian_2d/test_2d_full_field_oracle.py::test_matvec_window_equals_full_field_end_to_end` (same 4 cases; bulk residual + boundary-block residual, `assert_array_equal`) | bit-identity | **S6.3** (the matvec OUTPUT oracle) | **STAYS** — this is the test that proves the relocated `_apply_2d_cartesian` body produces byte-identical output. |
| **window≡full solve/residual (d=1/2/3)** | `tests/sn/sweep/core/test_sweep_graph_window_equivalence.py::test_solve_window_equals_full_field` + `::test_residual_window_equals_full_field` (synthetic shapes incl. non-square (12,7),(5,9) + d=3 (3,2,3),(4,3,2)) | bit-identity | **S6.3/S6.4** | **STAYS** — d-generic window≡full at the graph layer. |
| **ScanMarch G2.c (residual≡oracle)** | `tests/sn/sweep/cartesian_2d/test_scan_march_equivalence.py::test_scanmarch_residual_equals_oracle` + `::test_scanmarch_sweep_equals_oracle` (`_NULP_BOUND`, NON-SQUARE 12×7/5×9, LS-4, 2G het aniso) | principled-equiv (nulp, cross-schedule) | **S6.2/S6.3** (the renamed `loss_action`) | **STAYS as principled-equiv** — `ScanMarch.loss_action ≡ FullFieldWavefront.loss_action` at nulp. NOT bit-identity (row-march vs anti-diagonal differ at FP-association BY CONSTRUCTION). The renamed gate of S6.8. |
| **SI≡Krylov≡k_inf (value-ground)** | `tests/sn/eigenvalue/test_keff_2d.py::test_si_krylov_heterogeneous_2g_nonflat_flux` + `::test_default_entry_hits_kinf` + `::test_2g_eigenvector` + `::test_homogeneous_exact` | value-ground (closed-form k_inf transfer-matrix; ≥2G) | **ALL S6 stages** (esp. S6.5, both doors must still converge to k_inf) | **STAYS** — the eigenvalue structural ground; ≥2G non-flat (cardinal rule). |
| **L7 φ=Q/Σ_t d=2** | `tests/sn/sweep/cartesian_2d/test_2d_octant_sweep_equivalence.py::test_2d_octant_sweep_closed_form_anchor` (`np.linalg.solve` 2×2) | value-ground | underpins G2.c transitively | **STAYS** (vv §1.5: ULP-distance necessary-never-sufficient). |
| **L8 k_inf=1.875 d=1** | `tests/sn/sweep/core/test_wavefront_cumprod_equivalence.py::test_cumprod_path_hits_analytical_kinf` | value-ground | underpins G2.c transitively | **STAYS**. |

**Collection-verified this session** (`-O`, `PYTHONPATH=<worktree>`,
`pytest --collect-only`) — the A2D-1 pin, the 3-case affine golden, the 4-case window≡full
sweep+matvec oracle, the 4-case ScanMarch sweep/residual/moment, and the 5-case
g_adjoint_reciprocity all collect:

```
35 tests collected in 0.31s
```

```
PytestConfigWarning: assertions not in test modules or plugins will be ignored
because assert statements are not executed by the underlying Python interpreter
(are you using python -O?)
```

(The advisory is Mode-8-relevant — `test_g_adjoint_reciprocity.py` uses bare `assert`,
inert under `-O`. See §5 and §7.)

---

## 2. The A2D-1 regeneration protocol (S6.3 — DO NOT pin a stale source hash)

S6.3 moves the `_apply_2d_cartesian` BODY off the operator INTO
`MovingFrontierWindow.loss_action` (the walk leaves the operator; only the `−C`
subtraction remains as glue). The A2D-1 test pins
`sha256(inspect.getsource(StreamingOperator._apply_2d_cartesian))`. After S6.3 the
method either (a) is DELETED from `StreamingOperator` (the body lives in the
representation), or (b) shrinks to a thin forwarder. EITHER WAY the source text changes →
the hash MUST change. The protocol — three rules, NON-NEGOTIABLE:

1. **Refresh `EXPECTED_SHA256` IN THE SAME COMMIT as the body move.** The A2D-1 docstring
   (`test_streaming_operator.py:1155-1158`) already prescribes the update-and-bump
   ritual; follow it. Add a history-line comment recording the regen reason ("S6.3:
   `_apply_2d_cartesian` body relocated into `MovingFrontierWindow.loss_action`").
2. **The regenerated hash is NOT the correctness gate.** A source-hash change proves
   nothing about output. The correctness proof that the relocated body is behavior-
   preserving is the **`window≡full` MATVEC oracle**
   (`test_2d_full_field_oracle.py::test_matvec_window_equals_full_field_end_to_end`,
   `assert_array_equal`) PLUS the **affine-carve golden** end-to-end. Run BOTH green
   in the same commit; cite them in the commit body as the output-identity evidence.
3. **If S6.3 DELETES `_apply_2d_cartesian` entirely** (the body now lives only in the
   representation), the A2D-1 test must RE-TARGET: change `inspect.getsource(...)` to point
   at the NEW home (`inspect.getsource(MovingFrontierWindow.loss_action)` — or its
   extracted `_interior_walk` after S6.4), regenerate, and bump the docstring. NEVER leave
   the test pointing at a deleted method (it would `AttributeError` at collection — a hard
   red, not a stale-pin, but still a wrong gate). The mechanism that makes a stale pin
   IMPOSSIBLE to ship silently: the test fails loudly (hash mismatch OR AttributeError)
   unless the implementer regenerates against the NEW source — there is no "passes with
   the old hash" failure mode here, by construction. The ANTI-pattern to guard is the
   reverse — pinning the OLD `_apply_2d_cartesian` text after the body moved, which would
   pin a now-DEAD method. The retirement audit (§8 S6.3 row) catches that: grep
   `_apply_2d_cartesian` post-S6.3; if the method is gone, the A2D-1 target MUST have moved.

**Why output-not-source:** the brief's anti-recommendation. The source hash is a tripwire
for "did someone touch this body without thinking"; once we DELIBERATELY move the body,
the tripwire's job transfers to the OUTPUT oracle. Pinning source identity across a
deliberate relocation would be pinning the implementation to its pre-carve shape — the
exact thing S6.3 is undoing.

---

## 3. The `loss_action` convention pin (S6.3 — make `(L+C)ψ` a tested invariant)

NEW foundation test. Today the convention is enforced only by the `−C` subtraction
living inside `_apply_1d`/`_apply_2d_cartesian` and asserted nowhere at the representation
boundary. S6.3 exposes `loss_action` as the representation's public method returning
`(L+C)ψ`; pin BOTH halves of `op.apply(ψ) = loss_action(op, ψ) − σ_t·ψ.bulk`.

**File:** NEW `tests/sn/operators/test_loss_action_convention.py` (foundation; `-O`-safe).

**Spec (drop-in-ready shape, ≥2G het — vv §1-group-degeneracy + §H2):**

```python
@pytest.mark.foundation
def test_loss_action_returns_LpC_not_L():
    """[L0] loss_action returns (L+C)ψ; op.apply = loss_action − σ_t·ψ.bulk (the −C glue).

    Pins the S6.3 convention (sn_sweep_strategy.md §S6.1 locked-decision 2): the
    representation returns the FULL within-group loss (L+C)ψ; the operator subtracts
    the collision diagonal C = σ_t⊙. Returning L·ψ would re-couple C into every walk.
    """
    # ≥2G heterogeneous slab (H2: non-flat σ_t so C is not a uniform scalar;
    # the redistribution term is out of cancellation).
    L, sn = _build_loss_streaming_operator(ng=2, het=True)   # mixture A 2g, random σ_t(ng,nx)
    rep = L.sweep_strategy                                    # the SpatialRepresentation
    rng = np.random.default_rng(2026)
    psi = _random_composite(sn, rng)                          # non-flat ψ on bulk + boundary

    LpC_psi = rep.loss_action(L, psi)                         # (L+C)ψ  — the representation half
    L_psi   = L.apply(psi)                                    # Lψ      — the operator (with −C glue)

    # (i) the −C glue is EXACTLY σ_t·ψ.bulk (verified ONCE, here):
    glue = L_psi.bulk.values - (LpC_psi.bulk.values - sn_sigma_t(L)[None] * psi.bulk.values)
    np.testing.assert_allclose(
        glue, 0.0, atol=0.0, rtol=0.0,
        err_msg="op.apply must equal loss_action − σ_t·ψ.bulk EXACTLY (the −C glue).",
    )
    # (ii) loss_action is (L+C)ψ, NOT L·ψ — they differ by σ_t·ψ (non-trivially, het σ_t):
    diff = LpC_psi.bulk.values - L_psi.bulk.values
    np.testing.assert_allclose(diff, sn_sigma_t(L)[None] * psi.bulk.values, rtol=1e-13)
    # guard against a vacuously-zero σ_t·ψ (would make (i) and (ii) trivially pass):
    if float(np.max(np.abs(diff))) < 1e-6:
        pytest.fail("σ_t·ψ ≈ 0 — the het/non-flat config degenerated; the pin is vacuous.")
```

**Cross-check via the pure-z octant (locked-decision 2c):** add a second test that drives
a config where one octant is pure-z (`mu_x=mu_y=0` for some ordinate at d=2 LS-4). For
that ordinate `loss_action` returns `σ_t·ψ` (so `L·ψ = 0` after the `−C`). Assert
`L.apply(psi).bulk.values[pure_z_ord] ≈ 0` to `rtol=1e-13`. This is the structural
cross-check that the `(L+C)` convention is correct, not just internally consistent — a
sign-flipped or `L·ψ`-returning representation would NOT zero the pure-z streaming term.

**Zero-source probe (the `diamond.py:363` fact):** add a third assertion — with the
volumetric source set to zero and ψ a probe, `loss_action(L, psi).bulk` at the cell layer
equals `(σ_t + Σ_s,0)ψ̄ − Σ_s,0·in` term-for-term (the `residual_kernel_batch` natively-
`(L+C)ψ` fact). Pin against a hand-evaluated 1-cell value (L0 term verification, vv §10).

**Mode-9 note:** this is NOT a splitting/acceleration FP-invariance gate, so the diagonal-
cubature requirement does not bind here — BUT the config MUST be ≥2G + heterogeneous (a
1G uniform box makes `σ_t·ψ` flat and the `(L+C)`-vs-`L` distinction degenerate; vv §H1/§H2).

---

## 4. The "one representation instance" discriminating test (S6.5 — CONCRETE, drop-in)

The structural payoff. After S6.5 `InvertibleOperator.solve` and `StreamingOperator.apply`
reach ONE `SpatialRepresentation` instance → L21 becomes a TYPE FACT. The test asserts
identity (`is`). It MUST FAIL on current code (two `default_for` calls → two instances) and
PASS after S6.5 — a genuine tripwire, not a tautology.

**File (drop-in-ready, exact path):**
`tests/sn/operators/test_one_representation_instance.py`

```python
r"""S6.5 — the 'one representation instance' discriminating test (#222 capstone).

The structural payoff of the operator/representation re-layering: the
``SpatialRepresentation`` (renamed ``SweepStrategy``) that ``StreamingOperator.apply``
uses for the matvec ``(L+C)ψ`` MUST BE (object identity) the instance that
``InvertibleOperator.solve`` uses for the forward substitution ``(L+C)⁻¹q``. This makes
the L21 invariant ('matvec ≡ sweep — two actions of ONE operator') a TYPE FACT enforced
by construction, closing ``coding-elegance`` Smell #16 (two doors to one operator's
representation).

NEGATIVE PRE-CONDITION (why this is a tripwire, not a tautology): at the pre-S6.5 HEAD
the two paths each call ``default_for(sn_mesh)`` INDEPENDENTLY —
``StreamingOperator.apply`` via the ``sweep_strategy`` cached_property
(``operator.py:1566-1584``), ``InvertibleOperator.solve`` via ``transport_sweep`` →
``default_for`` (``operator.py:2530`` → ``sweep.py:239-241``). ``default_for`` constructs
a FRESH ``cls(mesh)`` each call, so the two instances are NOT identical today. This test
THEREFORE FAILS at pre-S6.5 HEAD and PASSES after S6.5 wires both doors to one instance.

``-O``-safe: uses ``pytest.fail`` (a function call that fires under ``python -O``), NEVER
a bare ``assert`` (vv Mode 8 — bare assert is stripped to a NO-OP under ``-O``).
"""
import numpy as np
import pytest

from orpheus.sn.geometry import SNMesh
from orpheus.sn.operator import CollisionOperator, StreamingOperator
# (import the project's standard 2-D Cartesian builders; mirror test_keff_2d / the
#  _build_2d helper in test_affine_carve_bit_identity.py — fuel|moderator, LS-4, 2G.)
from tests.sn._fixtures... import build_2d_cartesian_loss_operands  # ← reuse existing builder


def _build_invertible_2d():
    """Build a composite InvertibleOperator A = L + C on a 2-D Cartesian ≥2G mesh.

    The composite ``L + C`` (StreamingOperator.__add__ with a CollisionOperator) returns
    an InvertibleOperator automatically (operator.py:2251-2255). Its ``.streaming``
    operand IS the leaf ``StreamingOperator`` whose ``apply`` does the matvec.
    """
    L, C = build_2d_cartesian_loss_operands(ng=2, lvl=4)   # L: StreamingOperator, C: CollisionOperator
    A = L + C                                               # → InvertibleOperator (composite algebra)
    return A, L


def _representation_used_by_apply(L: "StreamingOperator"):
    """The SpatialRepresentation StreamingOperator.apply reaches (the leaf's cached one)."""
    return L.sweep_strategy            # S6.2 renames the property → spatial_representation; update here.


def _representation_used_by_solve(A: "InvertibleOperator"):
    """The SpatialRepresentation InvertibleOperator.solve reaches.

    PRE-S6.5: solve → transport_sweep → default_for(sn_mesh) builds a FRESH instance per
    call; there is no stored handle to compare, so the test reaches it via the SAME public
    surface solve uses (a fresh default_for) and asserts it is NOT the apply-instance
    (documents the two-doors defect). POST-S6.5: A exposes ONE shared instance (e.g.
    ``A.spatial_representation`` or ``A.streaming.spatial_representation`` — the door solve
    actually calls). Update this accessor to the S6.5 single-source handle.
    """
    # S6.5 target: solve calls representation.sweep on the SAME object apply uses.
    return A.spatial_representation     # ← the S6.5 single handle (does not exist pre-S6.5).


@pytest.mark.foundation
def test_apply_and_solve_share_one_representation_instance():
    """[L0 structural] apply's representation IS solve's representation (object identity)."""
    A, L = _build_invertible_2d()
    rep_apply = _representation_used_by_apply(L)
    rep_solve = _representation_used_by_solve(A)
    if rep_apply is not rep_solve:
        pytest.fail(
            "S6.5 NOT satisfied: StreamingOperator.apply and InvertibleOperator.solve reach "
            f"DIFFERENT SpatialRepresentation instances (apply={rep_apply!r} id={id(rep_apply)}; "
            f"solve={rep_solve!r} id={id(rep_solve)}). L21 (matvec ≡ sweep) is a coincidence, "
            "not a type fact — the two-doors smell (#222 S6.0) is still open."
        )
```

**Construction notes for the implementer (verified this session):**
- The composite `A = L + C` returns `InvertibleOperator` automatically
  (`operator.py:2251-2255`, "Operator algebra dispatch"). `A.streaming` is the leaf `L`
  (`operator.py:2360-2361`, alias for `self.a`). So `A.streaming.sweep_strategy` is exactly
  the representation `apply` uses TODAY.
- Pre-S6.5 there is NO shared handle, by design. To make the test a CLEAN
  fail-now/pass-after tripwire WITHOUT requiring the S6.5 API to exist before S6.5,
  stage it in two steps: (i) at S6.0-prime (now) land it `xfail(strict=True,
  reason="S6.5 unifies the two doors; A.spatial_representation single-handle not yet
  minted")` keyed on the `AttributeError`/`is-not` so the standing red is RECORDED; (ii)
  at S6.5, the `xfail` FLIPS to xpass — remove the marker in the S6.5 commit. This is the
  recorded-gap discipline ([[scan-march-verification]] G3.b pattern). The `xfail` reason
  MUST name the unlocking API (`A.spatial_representation`) so it is not forgotten.
- **2-D config is mandatory** (the two-doors defect is most visible where `apply` and
  `solve` route through genuinely different methods — the d=2 window matvec
  `_apply_2d_cartesian` vs the d=2 sweep). A 1-D mesh works too (both route CumprodScan),
  but d=2 is the discriminating config. Use the SAME `_build_2d` builder as the affine
  golden so the operands are a known het-aniso ≥2G case.
- `-O`-safe: `pytest.fail` fires under `-O`. The `xfail` marker also evaluates under `-O`.

**Why this is not a tautology:** it compares two INDEPENDENT selection paths that today
each `default_for(...)` a fresh object. A tautology would compare `L.sweep_strategy` to
itself. Here the LHS comes through `StreamingOperator.apply`'s door and the RHS through
`InvertibleOperator.solve`'s door — the test is exactly the "are the two doors one door"
question. Pre-S6.5 the answer is structurally NO (verified: two `default_for` call sites).

---

## 5. The curvilinear angular-transpose hazard (S6.3)

S6.3 HAZARD (sn_sweep_strategy.md §S6.3): when the `_apply_1d_transpose` body moves into
`CumprodScan.loss_action_transpose`, the curvilinear **second triangular factor** — the
angular adjoint `closure.angular_adjoint` (`operator.py:727`, reached via
`_compute_LpC_transpose` `operator.py:579`) — MUST ride with the spatial transpose. A
spatial-only reverse sweep silently DROPS the angular transpose → a wrong adjoint that
still looks plausible (the spatial part is right; only the curvilinear angular
redistribution adjoint is missing).

**EXISTING coverage (NAMED — use it, do NOT write a new numerical reciprocity test):**
`tests/sn/operators/test_g_adjoint_reciprocity.py::test_g_adjoint_reciprocity_full_block`
— parametrized over `["slab", "sphere", "cyl", "slab_2g", "sphere_2g"]` (collection-
verified this session). The **sphere** and **cyl** cases drive `A.H.apply(phi)` through
`apply_transpose → sweep_strategy.residual_transpose → _apply_1d_transpose →
_compute_LpC_transpose → closure.angular_adjoint`. The test asserts the DEFINING adjoint
property `⟨Aψ,φ⟩_G = ⟨ψ, A.Hφ⟩_G` with an INDEPENDENT `_g_inner` metric (anti-R1). If the
S6.3 move drops the angular adjoint on a curvilinear mesh, the sphere/cyl reciprocity
residual BLOWS UP (the spatial-only transpose is not the G-adjoint) — the test catches it.

**TWO actions required at S6.3 (do not skip):**
1. **Mode-8 fix (CRITICAL).** `test_g_adjoint_reciprocity.py` uses BARE `assert rel < 1e-12`
   (`test_g_adjoint_reciprocity.py:212`) — **inert under `python -O`** (the collection
   advisory above is exactly this file). As-is, the curvilinear-transpose tripwire CANNOT
   FIRE under the canonical `-O` invocation. S6.3 MUST migrate these asserts to
   `np.testing.assert_array_less` / `pytest.fail` (function calls that fire under `-O`)
   BEFORE relying on this test as the S6.3 hazard gate. Until migrated, the gate is a
   FALSE green (vv Mode 8). This is a prerequisite, not optional.
2. **Add the sphere_2g companion to the S6.3 gate set** (already parametrized) and assert
   it runs `-O`-green post-migration. The 2g curvilinear case is the strongest stressor
   (group coupling × angular redistribution × curvilinear closure).

**If the implementer prefers an explicit S6.3-local guard** (in addition to reciprocity):
a `loss_action_transpose`-vs-`_apply_1d_transpose` bit-identity pin on a sphere mesh —
`np.testing.assert_array_equal(CumprodScan(sphere).loss_action_transpose(L, phi).bulk.values,
L._apply_1d_transpose(phi).bulk.values)` immediately after the body move (the move is a
pure relocation → byte-identical). This is the cheapest "the body did not change when it
moved" pin; the reciprocity test is the "the body is mathematically the adjoint" pin. Ship
BOTH: relocation-identity (cheap, S6.3) + reciprocity (semantic, already exists, fix Mode-8).

---

## 6. The `as_scipy_linop` Krylov round-trip (S6.7 gate)

S6.7 requires `as_scipy_linop(representation/operator)` still drives the SN Krylov inner
loop after the rename. VERIFIED THIS SESSION:
- `as_scipy_linop` (`numerics/operator.py:1661-1711`) wraps `op.apply` as the scipy
  `matvec` (and `op.apply_transpose` as `rmatvec` when `CAP_APPLY_TRANSPOSE`). After S6.3
  `op.apply` routes through `representation.loss_action − σ_t·ψ` — so the adapter is
  agnostic to the rename PROVIDED `op.apply` keeps its signature + return type. The gate
  is: `op.apply` unchanged at its public boundary ⟹ `as_scipy_linop` unchanged.
- ⚠ NUANCE (verified): the SN Krylov eigenvalue/fixed-source inner does NOT call
  `as_scipy_linop` directly — `KrylovAcceleration.solve` (`numerics/iteration.py:755`)
  builds its OWN `spla.LinearOperator((n,n), matvec=A_matvec)` where `A_matvec` closes over
  `LC.apply − Σ gᵢ.apply` (`solver.py:230-257` `_within_group_krylov`). So the load-bearing
  contract is `LC.apply` (the composite `InvertibleOperator.apply` → its `.streaming`
  leaf's `apply` → `loss_action`). The `as_scipy_linop` adapter is exercised by the
  diffusion BiCGSTAB path + `tests/numerics/test_operator.py`.

**Gate (NAME the existing tests — the round-trip is covered transitively):**
- **SN Krylov inner:** `tests/sn/eigenvalue/test_keff_2d.py::test_si_krylov_heterogeneous_2g_nonflat_flux`
  (drives the GMRES inner via `_within_group_krylov` → `LC.apply` → `loss_action` post-S6).
  This is the load-bearing SN gate (the inner actually runs GMRES on `(L+C−S−B)`).
- **`as_scipy_linop` adapter itself:** `tests/numerics/test_operator.py` (the adapter's own
  unit coverage — `CAP_APPLY` requirement, `matvec`=`op.apply`, `rmatvec`=`op.apply_transpose`).
  This pins the adapter contract the rename must not break.
- **Diffusion BiCGSTAB co-tenant:** the diffusion solver round-trips through the SAME
  `as_scipy_linop` (`diffusion/solver.py`). The diffusion eigenvalue/fixed-source tests
  confirm the shared adapter still drives both consumers. (S6 touches only SN, but the
  shared adapter must keep its contract — assert no diffusion regression.)

**S6.7 action:** after the S6.2 rename + S6.3 move, run the Krylov eigenvalue gate
(`test_keff_2d::test_si_krylov_heterogeneous_2g_nonflat_flux`) + `tests/numerics/test_operator.py`
`-O`-green. If green, the round-trip survives the rename. No NEW test required — the
contract (`op.apply` stable) is already covered; STATE this in the S6.7 commit rather than
adding a redundant gate.

---

## 7. Mode-8 and Mode-9 discipline (explicit)

- **Mode-8 (`-O` strips bare `assert`).** EVERY S6 tripwire/sentinel MUST be `-O`-safe:
  `np.testing.assert_*`, `pytest.fail`, or `raise AssertionError` — NEVER a bare `assert`.
  The NEW tests (§3, §4) follow this. The EXISTING `test_g_adjoint_reciprocity.py` does
  NOT (bare `assert`, §5) — S6.3 MUST migrate it before leaning on it as the curvilinear-
  transpose gate. The collection advisory pasted in §1 is the live evidence that bare
  asserts exist in the S6 gate neighbourhood. Grep `^\s*assert ` in every new S6 test file;
  expect zero.
- **Mode-9 (FP-invariance only in a degenerate regime).** S6 is a RE-LAYERING, not a
  splitting/accelerator — the converged FP does NOT change at all (Fork B1: defaults
  unchanged). So the classic Mode-9 splitting trap does not directly apply. BUT the
  convention pin (§3) and any equivalence gate MUST still run on an ANISOTROPIC +
  HETEROGENEOUS + ≥2G config (NOT the isotropic-reflective box): the `(L+C)`-vs-`L`
  distinction (§3), the curvilinear angular adjoint (§5), and the two-doors d=2 matvec (§4)
  ALL degenerate to triviality on a 1G flat box. The §3/§4 specs already mandate ≥2G het;
  the §5 gate uses sphere/cyl (curvilinear, non-flat by geometry). If S6 ever flips a
  default (Fork B2), the FP-invariance gate inherits the full Mode-9 discipline from
  [[scan-march-verification]] §G4 (diagonal cubature + vacuum/streaming).

---

## 8. Per-stage gate table (S6.2–S6.5; S6.6 deferred)

| Stage | Gates that MUST stay green | NEW test(s) turning green | bit-id vs principled | `-O`-safety note |
| ----- | -------------------------- | ------------------------- | -------------------- | ---------------- |
| **S6.2** rename `SweepStrategy→SpatialRepresentation`, `residual→loss_action`, `residual_transpose→loss_action_transpose` (bodies still delegate to `_apply_*`) | A2D-1 (UNCHANGED — body not moved yet); affine-carve golden (byte-id); `window≡full` sweep+matvec; ScanMarch G2.c (renamed `loss_action`); SI≡Krylov≡k_inf; `test_g_adjoint_reciprocity` (after Mode-8 fix); `tests/numerics/test_operator.py` | (rename touches every call site; the §4 one-instance test stays `xfail`) | **bit-identity** (pure rename, zero behavior change) | Confirm no NEW bare assert; the renamed gate names update in `test_scan_march_equivalence.py` + `test_2d_full_field_oracle.py` |
| **S6.3** move walk OFF operator (`apply = loss_action − σ_t·ψ`); delete `_apply_*`; `loss_action_transpose` carries `closure.angular_adjoint` | affine-carve golden (byte-id, OUTPUT); **`window≡full` MATVEC oracle** (the output-identity proof for the relocated body); ScanMarch G2.c; SI≡Krylov≡k_inf; **`test_g_adjoint_reciprocity` sphere/cyl/sphere_2g** (curvilinear-transpose hazard, §5 — Mode-8-fixed) | **§3 loss_action convention pin** (`(L+C)ψ` + `−C` glue + pure-z + zero-source); **§5 relocation-identity pin** (`loss_action_transpose ≡ _apply_1d_transpose` on sphere) | **A2D-1 source-hash REGENERATES** (§2 — same commit, output-identity via oracle); everything else **bit-identity** (relocation is byte-preserving) | §3/§4/§5 tests `-O`-safe (`np.testing`/`pytest.fail`); MIGRATE `test_g_adjoint_reciprocity` off bare assert FIRST |
| **S6.4** `_OctantWalk2D` base; `MovingFrontierWindow`+`ScanMarch` provide only `_interior_walk`; retire the two Fork-B1 IOU twins | `window≡full` sweep+matvec (bit-id, the shared frame must be byte-preserving); ScanMarch G2.c; affine-carve golden | (none new — consolidation; the §4 one-instance test still `xfail`) | **bit-identity** (extract-shared-frame refactor); retirement audit = the two IOU notes (`operator.py:1828` + `sweep.py`) gone | grep the IOU comment strings → expect zero; no new asserts |

> **S6.4 has its OWN full gate plan** — see the sibling memo
> [[s6-4-unified-walk-verification]] (`s6_4_unified_walk_verification.md`). It
> supersedes this one-row summary: the one-walk SPY test, the A2D-1 RETIRE
> recommendation, the DAG-ownership structural tests (mesh-no-`sweep_graphs` +
> `geometry.py:343/350` None-slot gone), the §6 anti-boolean guard, and the
> (a)→(d) sub-step staging. Read it before writing S6.4 code.
| **S6.5** unify the two doors → ONE `SpatialRepresentation` instance reached by both `apply` and `solve` | affine-carve golden (byte-id — same convergence); SI≡Krylov≡k_inf (BOTH doors converge to k_inf); `window≡full`; ScanMarch G2.c | **§4 one-representation-instance test FLIPS `xfail`→xpass** (remove marker in this commit) — the structural payoff | **bit-identity** (wiring change, not a numeric change — `loss_action→(L+C)ψ`+`−C`, but `sweep→(L+C)⁻¹q` DIRECTLY; keep the asymmetry) | §4 uses `pytest.fail` (fires under `-O`); the `xfail` removal is the visible green flip |
| **S6.6** (DEFERRED) `ExplicitMatrix` representation | — | (future) `spsolve_triangular(M,q)==CumprodScan.sweep(q)` + `M@ψ==CumprodScan.loss_action(ψ)` to nulp; `scipy.sparse.tril(M)` exactly zero above visit-order diagonal | principled-equiv (nulp) when built | n/a now |

**Retirement audit (elegance-enforcer deliverable, gated at S6.5):** zero `_apply_*` methods
survive on `StreamingOperator` (grep `def _apply_1d\|def _apply_2d_cartesian\|def
_apply_full_field\|def _apply_2d_cartesian_scanmarch\|def _apply_1d_transpose` → expect the
bodies MOVED into representations, not lingering on the operator); ONE octant frame
(`_OctantWalk2D`); the two IOU twins retired. The §4 one-instance test green is the
construction-level proof Smell #16 is closed.

---

## 9. Deselect / standing reds (inherited)
- **DESELECT** `tests/sn/eigenvalue/test_keff_slab.py::test_heterogeneous_absolute_keff`
  (#212 `continuous_get` hang).
- Held reds untouched by the 2-D-gated re-layering: **#206** cyl-matvec (deferred
  curvilinear adjoint — note S6.3 routes the curvilinear transpose through
  `_apply_1d_transpose`, the SLAB/sphere/cyl 1-D path, NOT the deferred 2-D Cartesian
  adjoint; #206 is the 2-D cylinder-matvec and stays `xfail`); **#195** MMS@160.
- The 2-D Cartesian ADJOINT deferral (`_DAGWavefront.residual_transpose` raises
  `NotImplementedError`, `sweep_strategy.py:369-383`) MUST be PRESERVED through S6 →
  `loss_action_transpose` keeps the raise for the multi-D Cartesian case (NEVER a silent
  wrong answer). Inherit the S0 gate `test_2d_cartesian_adjoint_raises_incompatible` shape
  if it landed; if not, add it: a d=2 Cartesian `loss_action_transpose` MUST raise
  `NotImplementedError`, not return garbage.

## 10. Self-improvement triggers (fired BEFORE delivery)
- **NEW failure mode?** S6 introduces NO new failure mode beyond the `vv-principles` table.
  The "two doors → coincidence not theorem" defect is `coding-elegance` Smell #16 (already
  cataloged); the curvilinear-angular-transpose drop is a relocation hazard covered by the
  existing reciprocity test + the Mode-8 fix. NO skill-table append required. The Mode-8
  hazard in `test_g_adjoint_reciprocity.py` (bare assert under `-O`) is a CONCRETE instance
  of the already-documented vv Mode 8 — fix it in S6.3, no new row.
- **Plan rejection?** N/A (pre-implementation). If qa/user rejects, log the counter-example
  per the AGENT.md self-improvement directive.

## 11. Pre-reads cross-check (file:line, VERIFIED this session @ 716a0c1)
- `orpheus/sn/operator.py:1462` — `apply` → `self.sweep_strategy.residual(self, psi)` (door 1).
- `orpheus/sn/operator.py:1483-1487` — the `−C` glue in `_apply_1d` (`LpC − σ_t·ψ`).
- `orpheus/sn/operator.py:1552-1556` — the transpose `−C` glue in `_apply_1d_transpose`.
- `orpheus/sn/operator.py:1566-1584` — `sweep_strategy` cached_property → `default_for` (door 1 selection).
- `orpheus/sn/operator.py:579` / `:727` — `_compute_LpC_transpose` / `closure.angular_adjoint` (the §5 angular second-triangular-factor).
- `orpheus/sn/operator.py:2530` — `solve` → `transport_sweep` (door 2).
- `orpheus/sn/operator.py:2251-2255` / `:2360-2361` — `L+C → InvertibleOperator`; `.streaming` alias for `self.a` (§4 construction).
- `orpheus/sn/sweep.py:239-241` — `transport_sweep` → `default_for(sn_mesh).sweep(...)` (door 2 selection — the SECOND `default_for`).
- `orpheus/sn/sweep_strategy.py:169-221` — the `SweepStrategy` Protocol (`sweep`/`residual`/`residual_transpose`/`supports`) that S6.2 renames to `SpatialRepresentation`/`loss_action`/`loss_action_transpose`.
- `orpheus/sn/sweep_strategy.py:369-383` — `_DAGWavefront.residual_transpose` deferral raise (preserve through S6).
- `orpheus/numerics/operator.py:1661-1711` — `as_scipy_linop` wraps `op.apply`/`op.apply_transpose` (§6).
- `orpheus/numerics/iteration.py:755` — `KrylovAcceleration` builds its own `spla.LinearOperator(matvec=A_matvec→LC.apply)` (§6 nuance).
- `orpheus/sn/solver.py:230-257` — `_within_group_krylov` (`matvec = LC.apply − Σ gᵢ.apply`).
- `tests/sn/operators/test_streaming_operator.py:1306-1310` — `EXPECTED_SHA256` A2D-1 pin (§2 regen target).
- `tests/sn/operators/test_g_adjoint_reciprocity.py:191-215` — `test_g_adjoint_reciprocity_full_block` (sphere/cyl curvilinear-transpose coverage; bare assert at :212 — Mode-8, §5).
```
