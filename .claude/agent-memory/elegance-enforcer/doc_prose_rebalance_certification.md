---
name: doc-prose-rebalance-certification
description: Reusable method for certifying "docstring/comment-only" doc-prose-rebalance batches (#231 Phase 2) — the invariance harness, the dropped-contract net, pointer honesty
metadata:
  type: project
---

Method for certifying a **docstring/comment-only** doc-prose-rebalance batch (the
#231 Phase 2 P2-* batches: trim HISTORY/campaign-tags/TWIN-teaching from docstrings,
delegate to `docs/theory/<part>/<file>.rst §<label>` pointers, keep CONTRACT inline).
Certified P2-D/E/F/G 2026-07-22 (9 files) — ALL CERTIFIED, zero MUST-FIX. Precedent:
P2-A pilot @ `97fb2a38`.

**Why:** the ONLY two failure modes of a prose-trim are (1) a behavior change smuggled
into a "doc-only" edit, and (2) a CONTRACT (shape/unit/raise/mutation-semantic/sign-index)
trimmed away without surviving elsewhere. Everything else (which tag to cut) is the batch
agent's judgment, already recorded in the per-batch map (`.claude/plans/phase2_code_prose/P2-*.md`).

**How to apply — the four load-bearing checks (in order):**

1. **Dual invariance harness proves ZERO behavior change.** Token-stream check (tokenize,
   drop COMMENT+STRING+layout, compare count AND sha of exact-text sequence) proves code
   structure identical. BUT it drops ALL strings, so it is blind to a changed raise-message /
   dispatch-string. Add an AST check: **strip EVERY bare string-expression statement**
   (`Expr(Constant str)` anywhere, not just leading docstrings — attribute docstrings after a
   `ClassVar` annotation are bare Expr statements too), then `ast.dump` and compare. This
   KEEPS raise-messages/dict-keys/`__all__`/dispatch-strings (none are bare Expr) so a diff
   means a real code-string changed. GOTCHA: a "strip leading docstring only" stripper
   FALSE-POSITIVES on the P2-E spatial-scheme files because they edit ClassVar attribute
   docstrings — strip ALL bare-str exprs and they pass. Also diff the code-string multiset
   explicitly to prove no raise/dispatch string moved. Harness: `/tmp/tokcheck*.py` pattern.

2. **The dropped-contract net (the highest-value MUST-FIX catch).** `git diff -- <file> |
   grep '^-'` then filter for contract-signal tokens (`Shape (ng,`/`(N,`, `raise`/`:raises`,
   `Returns`/`:returns`, `in-place`/`mutat`/`overwrit`, `REQUIRED`/`must be`, units) and
   EXCLUDE campaign/wave/carve/taxonomy/retir/former narration. For each surviving deleted
   contract line, grep the CURRENT file to confirm the contract survives ELSEWHERE. In P2-*
   every hit was a dedup/relocation (a module-docstring restatement of a class/method SSOT) —
   the trim IMPROVED single-sourcing (Cardinal Rule 2). None were losses. But VERIFY each;
   do not assume.

3. **In-pass fixes verified against LIVE code, never the map's word** (Cardinal Rule 1).
   P2-* carried 4 fix classes worth re-checking every time: (a) a stale trait/Protocol
   description corrected to match a live `ClassVar` value — grep the live value, read the
   corrected prose; (b) a mislabeled adjoint — "Hilbert"→"Euclidean" for `apply_transpose`
   (Euclidean = plain `ᵀ`; Hilbert = metric `.H` applied AROUND by `_AdjointOperator`) —
   confirm summary==body==sibling-file convention; (c) a repointed dangling `:ref:` — grep
   the OLD anchor is dead tree-wide AND the NEW anchor is defined; (d) a `.. math:: :label:`
   cut from a NON-automodule'd module docstring — grep the label across docs/+tests/+orpheus/
   for ZERO `:eq:`/verifies/catches consumers before crediting orphan-safety.

4. **Pointer honesty = resolve ALL + content-spot-check ≥3/batch.** Cheap comprehensive
   check: collect every `.. _anchor:` / `:label:` in `docs/theory/**` (≈1260 defined), then
   confirm every named `§word-word` / `:ref:` in the 9 files resolves (skip `§1`/`§13` —
   those are literature/issue citations, not theory pointers). Then for ≥3/batch READ the
   landing section head to confirm it carries the DELEGATED concept, not just that the anchor
   exists. NOTE seen 2026-07-22 (honest, not a defect): a pointer can land on a section whose
   body is nested under `.. todo:: Archivist expansion needed` — the concept was still fully
   present there; and a HISTORY-narration cut correctly points at a "Superseded by #NNN"
   retained-historical section.

**Batch texture (calibration):** an ABC/operator-algebra file (operator.py) is
~90% CONTRACT — the docstrings STATE the binding laws, so the cut is NECESSARILY MODEST
(−51 ln) and that is CORRECT, not a shortfall. The leading-term preconditioned-splitting law
(`OperatorSum.is_invertible` "ONLY the LEADING term need be invertible") is the canonical
must-keep: cut it and the next modifier "fixes" it to require both. Curvilinear files
(P2-F) KEEP every sign/index convention at point of use (`c_out=α/τ`, `c_in=(1−τ)/τ·α+α`,
`faces[g,m,i]` off-by-one, τ_raw, seed `t`) — the ERR-026 family is where these bit hardest.
A "spot-check target" the task names may have NO referent in the file (the "replace()-routing
law" was never in operator.py at HEAD) — report that plainly, do not manufacture a finding.
