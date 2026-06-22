---
name: issue-236-phase1b-curvilinear-selection
description: #236 Ph1b honest (scheme×geometry) selection — PASS-WITH-NITS; supports_curvilinear trait + _curvilinear_capability single-source gate; the two-parallel-predicate (selection vs pairing) keep-separate ruling
metadata:
  type: project
---

# #236 Phase 1b — honest curvilinear scheme selection (uncommitted review)

PASS-WITH-NITS (no blockers). Closes a dishonest-selection gap: LD IS
`is_affine_scannable` (1-D Schur), so curvilinear-LD PASSED `CumprodScan.supports`
and only failed via a mid-sweep `_require_slab` raise. Fix = `supports_curvilinear:
ClassVar[bool]` trait (DD `True`, LD `False`, base `False`) + `_curvilinear_capability(mesh)
-> Compatibility` single-source gate composed at 3 sites (`CumprodScan.supports`,
`ScanMarch.supports` 1-D arm, top-of-`default_for`).

## The load-bearing ruling (answers "should it cohabit with 1a's pairing.py?")
NO — keep `_curvilinear_capability` (loss_representation.py) and
`pair_diffusion_limit_consistent` (spatial/pairing.py) SEPARATE. They are
parallel in spirit (both per-axis-trait conjunctions) but live on DIFFERENT
axes and speak DIFFERENT return types to DIFFERENT consumers:
- `pairing.py` = (scheme ⊗ angular-closure) tensor-pair → bare `bool` → frontend
  PHYSICS pairing-validity query (diffusion-limit). LMM-1987 ∧ BMC-2010.
- `_curvilinear_capability` = (scheme × geometry) → `Compatibility(ok,reason)` →
  sweep-strategy DISPATCHER selection/capability query (gray-out reason).
Co-locating would conflate "is this pairing physically sound?" with "can this
strategy sweep this mesh?" — a real concept-merge, NOT a single-source win.
Two predicates is NOT the unify-after-two trigger (different return types argue
AGAINST premature unification). Also NOT a job for SNMesh-construction Pattern-4:
a curvilinear-LD mesh is a LEGAL object a frontend must be able to construct +
query `supports()→graceful False`; making it unconstructable denies the gray-out.
Selection-layer rejection is the right granularity.

## Defense-in-depth verdict (keep ALL THREE — none dead)
Three guards at three DISTINCT bypass boundaries:
1. selection guard (`supports`/`default_for`) — honest for FRONTENDS (gray-out).
2. `affine_scan_coefficients` raise (linear_discontinuous.py:741) — fires at
   CollisionCache build for a DIRECT-CONSTRUCT bypass; catches RUNTIME-data
   curvature (`dA_w/c_out != 0`) the selection layer can't see.
3. `_require_slab` (linear_discontinuous.py:194) — per-cell kernel fail-safe on
   `angular_upstream is not None`, innermost.
Selection guard does NOT subsume the raises — it makes the honest-selection path
honest; the raises defend the BYPASS paths. Retiring any re-opens a silent-wrong-
math route for its boundary.

## NITs (1 do-now)
- DO-NOW (Finding 5, anti-pattern-11 / rationale-invariant tell): the
  top-of-`default_for` check (loss_representation.py:1796) is JUSTIFIED (surfaces
  the SPECIFIC curvilinear reason; the loop fall-through :1804 names DIMENSIONALITY
  which misleads) and pinned by `test_curvilinear_ld_rejected_at_selection`
  (asserts `"curvilinear" in reason` ∧ `"no sweep strategy supports" not in`). BUT
  the comment doesn't warn a future DRY-deletion away or name the pin → a
  well-meaning "this is duplicated with the loop, delete it" edit silently
  regresses message quality, caught by NO structural test. Harden the comment:
  state the non-redundancy invariant + name the pinning test.
- RECORD (Finding 6): `ScanMarch.supports` 1-D arm is structurally UNREACHABLE via
  `default_for` (CumprodScan registered first, also gates `is_1d`) but CORRECT to
  keep — independently frontend-callable + future-live if registry reorders; the
  single-source `_curvilinear_capability` keeps it honest whether or not reachable.
- RECORD (Finding 7): `_curvilinear_capability` dereferences `mesh.scheme.*`
  unguarded — SAFE (SNMesh.__init__ geometry.py:271 always assigns a concrete
  scheme, default DD; illegal "mesh w/o scheme" unrepresentable by construction).

## What was nailed (praise)
Trait default polarity (`False`=slab-only-until-cited, mirrors is_affine_scannable
/transverse_coupling_is_facewise/diffusion_limit_consistent); LD's explicit-`False`-
with-citation (NOT silent inherit — declares the selection consequence); Protocol-
member append in declaration order; genuine-bool teeth list gained the 7th trait
(closes the -O truthy-non-bool footgun — `_curvilinear_capability` reads it in
boolean context); `supports()` single-line→layered guard reads like the
precondition lattice with per-failure reasons (more honest than the old conflated
message). Mocks minimal+correct (Identity→False/FakeCurvilinear→True mirror their
angular_upstream behavior; `_fake` `curv_capable=True` default keeps d=3 pins green).

Verify: 24 passed standard; 7 passed under `python -O` (pytest.fail/raises, no bare
assert). Uncommitted at review time; Phase 1a already at HEAD.
