# P2-D — `orpheus/numerics/operator.py` code-prose rebalance

Phase 2 of #231, batch **P2-D** (the LinearOperator ABC + composite/
adjoint/inverse machinery). Pre-edit: **3119 lines, 1779 string-token
(docstring) + 325 comment lines.** Branch `docs/sn-doc-architecture`.
Pilot calibration: `.claude/plans/phase2_code_prose/P2-A_scattering.md`.

## Headline finding (calibrates the remaining operator batches)

**This ABC is CONTRACT-DOMINATED (~90%). ZERO MOVED. The rebalance is
NECESSARILY MODEST — and that is the correct outcome, not a shortfall.**
Unlike the pilot's `scattering.py` (huge MATH-derivation TWIN → 36 %
cut), `operator.py` is the algebra ITSELF: its docstrings STATE the
binding laws (closure/composition/adjoint-swap/homomorphism), the
raise-conditions, and the typing-rationale that guards plausible wrong
"simplifications". Per the brief's binding ABC rule — *"the capability-set
laws, composition/closure laws (`replace()`-routing, `.H` propagation,
sum/product algebra), and every raise-condition are CONTRACT"* — almost
none of it is cuttable. The Haiku classifier's **28 MOVED essentially ALL
collapsed** to CONTRACT (the laws) or a couple of theory-pointer swaps.
The realized cut is **−51 prose lines (−9 comment, −42 docstring),
67.5 % → 66.9 % prose share** — HISTORY-narration blocks + inline
campaign-step provenance + two use-case trims + one theory-pointer swap.

### Why the classifier's MOVED verdicts were wrong here (measured)

The classifier proposed relocating "closure laws", "composition laws",
"three-layer surface", "biproduct roles", "codomain semantics", etc. to
new theory sections. But for an ABC file the LAW **is** the contract a
modifier needs in-file: `OperatorSum.is_invertible`'s "ONLY the LEADING
term need be invertible (preconditioned-splitting)" is CONTRACT — cut it
and the next modifier "fixes" it to require both. Every one of the 28
proposed MOVED is either (a) a law that must stay at point-of-use, or
(b) categorical exposition the operator-algebra book **already** carries
(TWIN, not MOVED — the book is complete, per the pilot headline). I wrote
NO theory-page content and recorded the (near-empty) MOVED candidate set
below.

## Structural safety (identical to the pilot — cutting is safe)

- `operator.py` is **NOT `automodule`'d** anywhere (grep clean) →
  docstrings invisible to Sphinx; no `:ref:`/`:eq:`/`:class:` link
  breaks on a cut; the `-E -W` baseline (1 warning) is definitionally
  unchanged. **Sphinx build not required** (pilot precedent).
- **No `.. math:: :label:` and no `vv-status`** in any docstring →
  cutting math orphans no `:eq:` target.
- **No `verifies(...)`/`catches(...)`** test marker points at an
  `operator.py` docstring (they point at RST labels).

## BATCH SPECIAL — dual-A bridge verification (crosswalk row 8)

**Verdict: the module head CONFORMS. No fix needed** (classifier report
confirmed). `notation.rst §notation-crosswalk` row 8 promises the bridge
"is stated at both ends (the module heads of `orpheus.numerics.iteration`
and `orpheus.numerics.operator`, and the solver page)". `operator.py`
lines 9–25 state this end in full, verbatim-consistent with row 8:

- poses `(A − Σᵢ gᵢ)ψ = q` (fixed source) AND `(A − Σᵢ gᵢ)ψ = (1/k)Fψ`
  (eigenvalue) — the variadic `(A − Σᵢ gᵢ)` record ✓
- "`A` is the INVERTIBLE resolvent operand" ✓ (= row 8 "invertible
  resolvent operand")
- "for SN the binding is `A = L + C` … with gains `(S, B)`, the honest
  within-group operator `L+C−S−B`" ✓ (= row 8 "the SN binding hands
  `A = L + C` with gains `(S, B)`, composing to the same honest
  `L+C−S−B`")
- "`F` … never a gain in the eigenvalue posing: the outer loop scales it
  by `1/k`" ✓ (= row 8 "Fission is never a gain … stays on the RHS under
  `1/k`")

The bridge equations (lines 9–25) were **left untouched**; the only edit
in the module head was trimming the campaign tag "carve P4" from a
downstream `#226` reference (line 34, outside the bridge).

## Pointer convention

Phase-ratified literal greppable form `docs/theory/<part>/<file>.rst
§<label>` for NEW theory pointers (2 added — both anchors verified to
resolve). Pre-existing `:ref:`/`:class:`/`:meth:`/`:mod:` roles that
survive a trim are kept as roles (pilot rule). One pre-existing **stale**
`:ref:` was repointed (see discrepancy below).

---

## Adjudication (by region; classifier line-map re-adjudicated)

| region (symbol) | classifier verdict | MY verdict | action |
|---|---|---|---|
| module head 1–65 (bridge + three-layer) | CONTRACT | **CONTRACT** | KEEP; trimmed only "carve P4" tag @ ln34 |
| typevar block 95–152 (invariance + covariant legs) | TWIN | **CONTRACT** (guards covariant "fix") | KEEP; trimmed "W-A collapse" → "#65" |
| BlockRole comment 188–206 (biproduct 2×2 + dispatch) | MOVED | **CONTRACT** (defines the enum + dispatch fact) | KEEP; trimmed "Wave O step O.2", "/ Wave O" |
| BlockRole enum 208–237 | CONTRACT | **CONTRACT** | KEEP; **HISTORY-cut** the O.1/O.4a.1-β/-γ ship/retire narration |
| SystemRole enum 240–285 | CONTRACT | **CONTRACT** | KEEP; campaign-name → `coupled_block_operator.rst §coupled-block-operator` pointer |
| `_join_block_roles` 288–319 | TWIN | **CONTRACT** (the join LAW + collapse-trigger) | KEEP; trimmed "Wave O", cut "retiring the hand-stamped … tag" |
| `_join_system_roles` 322–340 | TWIN | **CONTRACT** | KEEP verbatim |
| `_BlockRoleMeta` + markers 343–394 | CONTRACT | **CONTRACT** | KEEP; trimmed "#290 P3" → "#290" |
| exceptions 397–476 (NotInvertible/MissingAdjoint/MissingAssembly/Incompatible/MatrixTooLarge) | CONTRACT | **CONTRACT** (raise-conditions) | KEEP; trimmed "taxonomy §12 step 6" ×2, "taxonomy §17 A2" |
| `_resolve_basis_shape` 479–505 | CONTRACT | **CONTRACT** | KEEP verbatim |
| `LinearOperator` base 508–540 | SPLIT (CONTRACT+TWIN) | **CONTRACT** | KEEP verbatim (the Notes on shape/dtype is a real constraint) |
| `block_role`/`system_role` attrs 542–576 | CONTRACT | **CONTRACT** (ClassVar/dataclass trap) | KEEP; trimmed "/ Wave O", "at O.2" |
| domain/codomain + predicates 578–659 | CONTRACT | **CONTRACT** | KEEP; **cut retired `CAP_SOLVE`/`CAP_APPLY_TRANSPOSE` lineage** |
| `apply` 661–678 | CONTRACT | **CONTRACT** | KEEP verbatim |
| algebra dunders 680–821 | SPLIT | **CONTRACT** | KEEP laws; **trimmed `__neg__`/`__truediv__` use-case lists**; `#65` tag; `__call__` retired-2-arg-BC history cut; `__and__` Grand-Report line-cite → `operator_tensor_network.rst §tensor-network-decomposition` |
| `adjoint`/`H` 827–876 | CONTRACT | **CONTRACT** | KEEP the metric law + eager raise; trimmed "(carve P4, spec §38)", Grand-Report cite |
| `as_matrix` 882–965 | SPLIT (TWIN+MOVED+CONTRACT) | **CONTRACT** | KEEP column formula/C-order/params/resource-gate/assembly-delegation; trimmed "taxonomy §2", **cut "né _as_dense" + "taxonomy §12 step 5 decided"**; "(stencil-assembly 2b, ruling R2)" tag |
| `_as_matrix_by_probing` 986–1007 | SPLIT (CONTRACT+MOVED) | **CONTRACT** (the oracle-gate rationale IS why it's named) | KEEP verbatim |
| `__repr__` 1012–1027 | COMMENT-cut | **CONTRACT** (comment explains the token format) | KEEP; trimmed "(carve P4 §44.F)" |
| narrowing-targets section 1030–1055 | MOVED | **CONTRACT** (why-NOT-runtime_checkable) | KEEP; trimmed "carve P4", "taxonomy §12 step 6" |
| SupportsInverse/Adjoint/Assembly 1058–1102 | CONTRACT | **CONTRACT** | KEEP verbatim |
| invertible/adjointable/assemblable 1105–1158 | CONTRACT | **CONTRACT** (TypeGuard-not-TypeIs) | KEEP; trimmed "spec §39.1", "spec §44.E" |
| `_AdjointOperator` 1171–1302 | SPLIT (CONTRACT+MOVED) | **CONTRACT** (swap law + eager gate + PEP-696 pin) | KEEP; trimmed "(carve P4, spec §38)", "(Wave O / O.2b)", "(Phase 2.5c, ruling R11)", "2.5c"; cut "bit-identical to the former in-line…" refactor-history |
| `OperatorSum` 1305–1524 | SPLIT (TWIN+MOVED+CONTRACT) | **CONTRACT** (all closure laws) | KEEP the laws + `#300`(SMW)/`#261`; trimmed "C4 F1", "taxonomy §12 step 4"×2, "taxonomy §11.1"; cut "transitional gated spelling retired at carve P4" + "Replaces the former hand-stamped … tag" |
| `OperatorProduct` 1527–1725 | SPLIT (TWIN+MOVED+CONTRACT) | **CONTRACT** | KEEP the laws + `#285`; trimmed "C4 F1", "carve P4"×2; **cut "before taxonomy §12 step 5 … reversed-product spelling" narration** |
| `ScaledOperator` 1728–1859 | SPLIT (TWIN+MOVED+CONTRACT) | **CONTRACT** | KEEP the laws; trimmed "C4 F1", "carve P4"×2, "carve P5", "taxonomy §13 I2", "Wave O" |
| `IdentityOperator` 1862–1887 | CONTRACT | **CONTRACT** | KEEP; trimmed "carve P4" |
| `ZeroOperator` 1890–1962 | SPLIT (MOVED+CONTRACT) | **CONTRACT** (codomain_zero contract) | KEEP the two-regime contract + `#208`; **trimmed the SN-fission-slot consumer-walkthrough + "B.5.2 fix"** |
| `_WrappedForward`/`_InvertibleForward` 1965–2018 | CONTRACT | **CONTRACT** | KEEP minimal-contract + per-sibling narrowing; **cut the "incidentally tighter … exposed the true minimum" + "carve P4 … caps-GATED retired" discovery-history** |
| `InverseWrapMixin` 2024–2119 | SPLIT (MOVED+TWIN) | **CONTRACT** (back-half laws + seeded-apply #285 + #280) | KEEP; **cut the extraction-trigger narration** (kept the defer-until-≥2 principle); trimmed "taxonomy §12/§13" ×4 |
| `InverseOperator` 2122–2190 | SPLIT (TWIN+MOVED) | **CONTRACT** (name-earning + one-realization-not-reciprocal-twin) | KEEP both (critical); trimmed "taxonomy §13/§12 step 5"; **cut the extraction narration** |
| `PermutationOperator` 2193–2303 | CONTRACT | **CONTRACT** | KEEP; trimmed "carve P4" |
| `IncomingOrdinateMaskTensor` 2306–2397 | SPLIT (MOVED+CONTRACT) | **CONTRACT** (mask vs ZeroOperator + projection + raise-surface) | KEEP verbatim (Grand-Report §16A.10 cite left — see note) |
| `PeriodicWrapOperator` 2400–2449 | SPLIT (MOVED+COMMENT-cut) | **CONTRACT** | KEEP the aliasing contract; **HISTORY-cut the "Wave 7 update" block + apply comment** |
| `TensorProductOperator` 2452–2590 | SPLIT (TWIN+MOVED) | **CONTRACT** (def + algebraic laws + numpy-relation guard) | KEEP verbatim; trimmed "carve P4"×2, "taxonomy §13 I2" |
| `SumOfTensorProductsOperator` 2593–2684 | SPLIT (MOVED+CONTRACT) | **CONTRACT** | KEEP verbatim (Grand-Report §15.2 cite left — see note) |
| `DiagonalOperator` 2687–2928 | SPLIT (TWIN+CONTRACT+MOVED) | **CONTRACT** (two regimes + params + back-compat alias + raise) | KEEP verbatim; trimmed "taxonomy §13" |
| `RankOneOperator`/`outer` 2930–3119 | SPLIT (MOVED+CONTRACT) | **CONTRACT** (dyad + transpose law + adjointable-iff-IPF) | KEEP verbatim; trimmed "campaign #276" → "#276"; **fixed stale `:ref:`**; **cut the "W2 as-built fix" aside** |

---

## Summary

### Line counts (`git show HEAD` vs edited; tokenize)

| metric | HEAD | edited | Δ |
|---|---|---|---|
| total lines | 3119 | 3068 | −51 (−1.6 %) |
| comment lines | 325 | 316 | −9 (−2.8 %) |
| string(docstring) lines | 1779 | 1737 | −42 (−2.4 %) |
| prose (cmt+str) | 2104 | 2053 | −51 |
| prose share | 67.5 % | 66.9 % | −0.6 pt |
| **code lines** | — | — | **0 (token-invariant)** |

`git diff --stat`: **155 insertions, 206 deletions**, 1 file.

### Verdict counts

- **POSING-HARMONIZE**: 0 (module head already conforms — bridge verified, see BATCH SPECIAL).
- **CONTRACT (kept)**: the whole file — every closure/composition/adjoint/
  homomorphism law, every raise-condition, every typing-rationale. ~40
  CONTRACT blocks were TRIMMED of embedded HISTORY (the pilot's "~28
  CONTRACT trimmed" pattern, at larger scale here).
- **HISTORY (cut narration)**: ~11 multi-clause narration cuts — BlockRole
  ship/retire, ZeroOperator consumer-walkthrough, PeriodicWrap Wave-7
  block, RankOne W2-aside, `_WrappedForward` discovery-story,
  `InverseWrapMixin` + `InverseOperator` extraction-story,
  `OperatorProduct.inverse` "before step 5", `_InvertibleForward`
  caps-retired, `_join_block_roles` retired-tag, `__call__` 2-arg-BC.
- **Campaign-step provenance trims** (inline, uniform rule): ~30 tags —
  ALL of `Wave O`/`Wave 7`/`carve P4`/`carve P5`/`spec §NN`/`taxonomy
  §NN step N`/`Phase 2.5c ruling RNN`/`W-A collapse`/`né _as_dense`/
  `O.2b`/`B.5.2`/`W2 as-built` removed; the LAW/constraint they annotated
  KEPT. **KEPT (rubric: live-issue one-liners + named patterns)**: bare
  `#NNN` issue anchors (#65/#208/#226/#261/#276/#280/#285/#300), named
  patterns with theory anchors (`Design C`, `coding-elegance Pattern 2`).
- **TWIN (cut → theory pointer)**: 2 — `__and__` Grand-Report line-cite →
  `operator_tensor_network.rst §tensor-network-decomposition`; SystemRole
  campaign-name → `coupled_block_operator.rst §coupled-block-operator`.
- **Use-case-list trims**: 2 (`__neg__`, `__truediv__`).
- **Stale-ref fix**: 1 (`:ref:`operator-algebra-adjoint`` →
  `:ref:`operator-adjoint``, Cardinal Rule 1 — see discrepancy).
- **MOVED**: **0**.

### Gates

- **GATE 1** `import orpheus.numerics.operator` → **OK**.
- **GATE 2** `-O -m pytest --collect-only -q` → **6652 tests collected**
  (matches expected).
- **GATE 3 token invariance**: code-token sequence **5960 == 5960,
  IDENTICAL** (sha256 pre/post; strings + comments excluded) — proves
  zero behavioral/code change. Diff audit confirms every changed line is
  a comment or docstring-prose line (no raise-message string edited; the
  lone ``TypeError`` grep hit is RST inline-literal markup in a docstring).
- **GATE 4** every cited §label resolves: both new literal pointers
  (`§coupled-block-operator`, `§tensor-network-decomposition`) resolve;
  both surviving `:ref:` targets (`operator-adjoint`, `operator-algebra`)
  resolve; stale `operator-algebra-adjoint` count now 0.
- **GATE 5 no theory-page edits**: confirmed (0 files under `docs/`
  changed).

### Hardest judgment calls (calibrate the remaining operator batches)

1. **The classifier's MOVED verdicts invert for an ABC file.** On a
   file that IS the operator algebra, "closure law" / "composition law" /
   "role classification" / "codomain semantics" are the CONTRACT a
   modifier reads in-file, NOT relocatable teaching. RULE for the
   remaining batches: on an ABC/base-class file treat the classifier's
   "MOVED closure/composition/role law" as **CONTRACT by default**; the
   theory-page TWIN is the *concept's derivation*, which the docstring
   does not re-teach (it states the law the class obeys).

2. **Discriminate CAMPAIGN-STEP codes from ISSUE anchors and NAMED
   PATTERNS.** I applied ONE uniform rule (internal consistency): trim
   `Wave O`/`carve PN`/`spec §NN`/`taxonomy §NN step N`/`Phase 2.5x`/
   `W-A collapse`/`né`/`O.2b` (landed archaeology git+theory own); KEEP
   bare `#NNN` (rubric: live-issue one-liners) and named patterns with
   theory anchors (`Design C`, `Pattern 2`). This is the "constraint
   stays / narration cuts" line drawn at *citation-vs-narration*: a bare
   `#280` is a citation (kept); a "carried as documented twins until the
   3rd sibling fired the extraction trigger" is narration (cut).

3. **A retired-system LINEAGE note is HISTORY even when terse.** "The
   RUNTIME successor to the `CAP_SOLVE` capability tag" reads like a
   one-clause aid, but `CAP_*` is fully retired — no live code references
   it, so the lineage is pure archaeology (aggressive-retirement). Cut it;
   the predicate law stands alone.

4. **A one-realization-not-reciprocal-twin / typing-choice rationale is
   CONTRACT even at length** (pilot precedent #4). `InverseOperator`'s
   "(a) `(1/c)·b ≠ b/c` in FP, (b) a units-dishonest reciprocal-XS field"
   guards a Cardinal-Rule-1 wrong change — KEEP verbatim. Same for the
   ClassVar/dataclass trap on `block_role`, the `Generic[Domain,Codomain]`
   PEP-696 pin, TypeGuard-not-TypeIs, why-NOT-runtime_checkable.

5. **Modest is correct — do NOT manufacture cuts to chase the pilot's
   36 %.** The −51-line result is the honest floor for a CONTRACT-dense
   ABC. Cardinal Rule 1 (doc = the LLM's brain) forbids stripping a law
   to hit a number. The remaining operator batches (fission/streaming/
   boundary/multiplication) will cut MORE than this ABC (they carry MATH
   TWIN the book owns) but LESS than the pilot's scattering.py.

### Reported discrepancy (found + fixed; not introduced by this pass)

**Stale cross-reference: `:ref:`operator-algebra-adjoint`` in the
`RankOneOperator` class docstring pointed at an anchor that does NOT
exist anywhere under `docs/`** (grep-verified). Since `operator.py` is
not `automodule`'d it renders plain-text with no `-W` warning (L-002), so
it was silently dangling. Live anchor is `operator-adjoint`
(`operator_adjoint.rst:1`). Repointed → `:ref:`operator-adjoint``
(Cardinal Rule 1 — a stale pointer in a docstring I was already editing).
Flagged per the pilot's discrepancy-reporting convention.

### Note — internal-report `Grand Report v3` cites left in place

`IncomingOrdinateMaskTensor` (§16A.10 line 3165), `SumOfTensorProducts`
(§15.2/§15A.2), `DiagonalOperator` (§9 `W`), `RankOneOperator` (§15) still
carry `Grand Report v3 §NN` design-doc cites. These are section
references to an internal design record, NOT campaign-STEP codes, so my
uniform trim rule left them (out of scope — a separate "internal-report
cite → theory-anchor" cleanup). The ONE I swapped (`__and__`'s
"§6.3 line 721 and §15.1 line 2044") was traded for a durable theory
pointer only because I was already editing that block and a fragile
double-line-number cite had a clean anchor replacement. **Flagged** as a
candidate follow-up if the main agent wants uniform internal-cite
retirement.
