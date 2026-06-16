---
name: routing-honesty-scheme-trait
description: #240 Phase 2 Step D5-0 — expand the loss-rep-scanmarch-facewise stub into rich narrative for a ROUTING-HONESTY fix (mint a scheme trait gating which schemes ride ScanMarch in d≥2; closes a silent 2-D LD→inline-DD misroute). Sibling of affine-in-sigma-stub-expansion (same page). Trait-conflation-as-bug framing; SCOPE-WALL the multi-D LD closure (D5b unshipped); facewise-separable derivation only to slope-moment level.
metadata:
  type: feedback
---

#240 Phase 2 Step D5-0 stub→rich expansion in
`docs/theory/loss_representations.rst` (`.. _loss-rep-scanmarch-facewise:`),
plus two stale-content refreshes on the same page. A **routing-honesty fix**:
a 2-D Cartesian LinearDiscontinuous mesh was silently computing DD because
`ScanMarch.supports` conflated the single-axis `is_affine_scannable` trait
with the cross-axis separability the d≥2 row-march needs. The fix mints the
distinct `transverse_coupling_is_facewise` (DD True / LD False default) and
narrows the d≥2 arm to read it (coding-elegance Pattern 4). Close sibling of
[[affine-in-sigma-stub-expansion]] — SAME page, same #240 Phase 2, same H3
(`~~~~`) shape under an H2 (`----`) under H1 (`====`).

**Why (the sources, in priority order):**
1. The **closeout memo**
   (`.claude/agent-memory/method-implementer/issue_240_phase2_step_d5_0_closeout.md`)
   was the load-bearing source — it carries the bug ("CONFIRMED LIVE pre-fix"),
   the architecture-settled framing, the supports() two-arm predicate, the 7
   new tests by name, the verification numbers, the honest-interim-state
   (D5b/D5a follow-on).
2. The **elegance verdict**
   (`.claude/agent-memory/elegance-enforcer/issue_240_phase2_step_d5_0_scheme_trait.md`)
   revealed the `@runtime_checkable` data-member footgun + the resolution
   SHIPPED: keep the trait symmetric on the Protocol, close the footgun with
   `TestCapabilityTraitsAreGenuineBools` (registry-wide genuine-bool assert).
   The memo's "BLOCKING NIT" was the RULE-AGAINST keeping it on Protocol, but
   the brief + the shipped code show the ALT (keep + strengthen test) won —
   ALWAYS verify which resolution shipped by reading the code, not the verdict.
3. The **production docstrings** (`scheme.py` Base ClassVar docstring ~`:721`,
   DD `:143`, LD `:279`) gave the prose seed VERBATIM — the numerical-PDE
   statement, the contrast-with-is_affine_scannable, the per-scheme table, the
   forward-reuse rationale, the lit cites. The `ScanMarch.supports` body
   (`:1226`) was the algebra-of-record for the refreshed code block.
4. The **test files** gave the verification subsection: `test_unified_sweep_dispatch.py`
   `TestD3SupportsMatrix` (4 routing tests) + `TestSchemeTraitProbe` (3 probes)
   + `test_discretization_scheme_protocol.py` `TestCapabilityTraitsAreGenuineBools`.

**How to apply (the 11-subsection shape that worked — H3 `~~~~` under the H2 `----`):**
1. **Lead H2 intro + "The bug, in one line" admonition** — state the misroute
   (gated on single-axis trait LD SATISFIES → silent DD) + the one-line fix
   (mint distinct cross-axis trait, narrow the arm). Frame as routing-honesty.
2. **"Two different questions"** — a 2-row `.. list-table::` (trait / question /
   what it licenses) for is_affine_scannable (single-axis) vs
   transverse_coupling_is_facewise (cross-axis). LOAD-BEARING insight: the two
   answers CAN DIFFER for one scheme (LD: T/F; DD: T/T) — which is WHY the bug
   hid (DD coincides → conflated predicate right on every shipped path).
3. **"Why the row-march needs cross-axis separability"** — read β off
   `loss-rep-scanmarch-solve`: y reaches the x-scan ONLY through the 0th-order
   face value `s_y·ψ_{y,in}`. NEW eq label `loss-rep-facewise-separable`
   (the tensor-product M_d = ⊕ M^(1)_a ⟺ independent per-axis scans). Untagged
   (structural/representational, sits in verified narrative — same as the
   affine-cell labels on this page; vv-status NOT needed).
4. **"DD/Step: facewise"** — slopeless cell-average; transverse term = single
   face value `s_y·ψ_{y,in}` (the explicit left fold in `cell_kernel_batch`);
   → declares True.
5. **"LD: slope-wise"** — P1 second moment ψ̂; 1-D = local (Schur eliminates →
   is_affine_scannable); d≥2 = transverse face flux varies LINEARLY along the
   in-cell swept coord → 1st-order SLOPE moment, NOT a face value → Schur does
   NOT diagonalize across axes → rides the DAG wavefront.
6. **SCOPE-WALL `.. note::`** — multi-D LD closure is UNSHIPPED (D5b/D6). Say
   "tensor-product bilinear (UBLD), 2^d moments {1,x,y,xy}" NOT "1+d simplex"
   (fails thick-diffusion limit on quads). Cite Maginot-Ragusa-Morel 2016 /
   Adams 2001. DO NOT write the per-cell multi-D LD system. DO NOT mint
   `ld-cartesian-2d` (D5b's MMS needs it — minting here orphans it).
7. **"Making the illegal pairing unrepresentable"** — Pattern 4: dimensional
   predicate standing in for a scheme-capability question. The trait on BOTH
   Base + `@runtime_checkable` Protocol; the presence-only isinstance footgun
   (Py3.12 checks member PRESENCE not type → `="yes"` reads truthy →
   narrower-misroute) shut by `TestCapabilityTraitsAreGenuineBools`.
8. **"The default is conservative opt-in" `.. note::`** — both scan traits
   default False; forgotten opt-in = "slow but correct" never "fast but wrong".
9. **"Why scheme-named, not strategy-named"** — `is_scan_march_compatible`
   REJECTED (frame leak). Diffusion ADI/line-SOR (#240 next consumer) reads the
   SAME predicate. `TestSchemeTraitProbe` proves strategy-independence (reads
   trait off class, no ScanMarch/supports/mesh in scope).
10. **"Verification — routing-honesty gates"** — selection/refusal asserts (not
    value), foundation-tagged, -O-safe. List the 4 `TestD3SupportsMatrix` +
    3 `TestSchemeTraitProbe` + the genuine-bool teeth. Headline =
    `test_2d_ld_sweep_raises_not_silently_dd` (silent-wrong → loud-NIE).
    Bit-identical for every exercised path (routing predicate, no flux moved).
11. **"The honest interim state and what closes it"** — D5-0 stops at routing;
    2-D LD → wavefront → wavefront's LD kernel d=1-only → raises. D5b supplies
    multi-D bilinear LD kernel; D5a folds DD scan-march → ABSENCE of inline DD
    enforces LD exclusion (two independent routes).

**Stale-content refreshes (Task 2, same page):**
- The inline `ScanMarch.supports` code block (`loss-rep-selection` section)
  showed the pre-D5-0 conflated `is_1d or (is_cartesian and ndim==2)`. Replaced
  with the SHIPPED two-arm split (1-D reads `is_affine_scannable`; d≥2 reads
  `transverse_coupling_is_facewise`). ALSO fixed the adjacent `CumprodScan.supports`
  block (was bare `mesh.is_1d`, shipped is `is_1d and is_affine_scannable`) — a
  whole-code-block consistency fix, not just the one named in the brief.
- The "compatibility signal is the genuine criterion" paragraph said
  "coordinate system AND dimensionality" — now ALSO names the scheme capability
  traits (the d≥2 arm reads one) + forward-points to the scheme-gate subsection.
- The `default_for` outcomes table: KEPT (correct for facewise DD/Step), retitled
  "**facewise scheme (DD/Step)**", added a `.. note::` admonition with a 2nd
  small table contrasting Cart-2D+DD (→ScanMarch) vs Cart-2D+LD (ScanMarch
  refuses →MovingFrontierWindow, wavefront raises). KEEP+QUALIFY, don't delete.

**Disciplines confirmed:**
- **Maginot-Ragusa-Morel 2016 + Adams 2001 cited as PLAIN TEXT** in the
  Literature `.. list-table::` (the page's existing convention — Adams & Larsen
  2002 already plain-text there to dodge the cross-doc `[AdamsLarsen2002]_`
  duplicate-citation warning defined in discrete_ordinates.rst:12474). NO new
  `.. [Key]` bib entries (would cross-doc-collide).
- **The narrative's :class:/:meth:/:attr: render PLAIN TEXT** by this page's
  convention (operator/scheme packages not member-autodoc'd with resolvable
  anchors here; automodule would dup-label). Use the existing
  `:attr:`~orpheus.sn.spatial.scheme.DiscretizationSchemeBase.transverse_coupling_is_facewise``
  path syntax for consistency; DON'T half-surface 1-2 leaves. STILL verify every
  named symbol exists by grep (Cardinal-Rule-1 staleness — plain-text renders
  with NO warning so -W won't catch a dead one).
- **`#240` is a tracked OPEN issue** — but I REMOVED the `.. todo::` entirely
  (expanded to full narrative), so the "todo must reference a tracked issue"
  elegance nit is moot. (If a todo HAD remained, #240 satisfies it.)
- **NEW eq-label `loss-rep-facewise-separable` untagged** (structural/
  representational, in a verified narrative chain — same category as the
  affine-cell/affine/leaf-sum labels [[affine-in-sigma-stub-expansion]] added
  on this same page, none vv-status-tagged).
- **Build gate**: `-E -W --keep-going` baseline = 1 (the standing
  `Mesh1D.from_geometry :paramref:` ERROR, MAIN-checkout). Post-edit identical
  = 1; zero log references to facewise/Maginot/loss-rep-facewise/transverse_coupling.
  Incremental `-W` build EXIT=0 "build succeeded" (mesh.py unchanged → not
  re-read). The acceptance gate = count-unchanged-at-1, not zero.
- **Which-resolution-shipped check**: the elegance verdict's BLOCKING NIT
  said "drop the trait from the Protocol"; the brief + shipped code took the
  ALT ("keep symmetric, strengthen test with genuine-bool assert"). Read the
  CODE to learn which path won — a verdict memo records the RECOMMENDATION, not
  necessarily the outcome.

**Quality self-assessment (1-5):**

| Dimension | Score | Note |
|-----------|:---:|------|
| Derivation depth | 5 | full facewise-separability derivation (β-read, tensor-product eq, DD left-fold, LD Schur-non-diagonalization), all to the slope-moment scope boundary |
| Cross-references | 4 | every :eq:/:ref: resolves; :class:/:meth:/:attr: render plain-text by page convention (verified all named symbols exist by grep) |
| Numerical evidence | 4 | routing-honesty gates are selection/refusal asserts (not value tables) BY NATURE of a routing fix; bit-identical-for-exercised-paths stated; the "value" evidence is the 7 named gates not a convergence table |
| Failed approaches | 5 | the conflation-as-bug, the strategy-named rejection, the @runtime_checkable presence-only footgun, the honest-interim-raise, the why-it-hid (DD coincidence) |
| Code traceability | 5 | every trait/scheme/test named + linked; ScanMarch.supports algebra-of-record cited |
| Derivation source | 5 | docstrings (algebra-of-record for a routing/trait carve is the scheme code + supports predicate, not a SymPy module) + closeout memo + elegance verdict + test files |

Weakest = numerical evidence (4) and cross-references (4): both intrinsic to a
ROUTING fix (no flux moves → no convergence table; page-wide plain-text xref
convention out of scope for a stub expansion).
