---
name: apply-solve-asymmetry — code-verify a draft already in the working tree
description: A docs task whose deliverable (an `apply-solve-asymmetry` subsection on operator_algebra.rst + expanded `streaming-action-pure-l` S8b stub) was ALREADY drafted as uncommitted working-tree edits at pickup; the real job was the code-verification pass the brief demanded ("verify against the code, do NOT contradict the implementation"), not fresh authoring. The thesis: forward apply distributes over (L+C); inversion does not — solve belongs to the bundled InvertibleOperator, apply to the leaves. SIBLING of the #238 decomposition-leaf note (same page, same StreamingOperator/loss-representation/fused-matvec context, same `=`/`-`/`~`/`^` marker ladder).
type: feedback
---

A second-generation sibling of [[feedback_decomposition_leaf_retirement_rationale]]
(#238) and [[feedback_affine_in_sigma_stub_expansion]] (#240 Step B) — same
`operator_algebra.rst` page, same `StreamingOperator` / `loss_representation`
/ fused-`loss_action` subject matter. The brief asked me to ADD an
`apply-solve-asymmetry` subsection near the `InvertibleOperator` /
capability-set discussion AND expand the `streaming-action-pure-l` `.. todo::
Archivist (#257 S8b)` stub.

**The decisive finding (Cardinal Rule 1 / [[feedback_drift_verification]]):
the entire deliverable already existed as UNCOMMITTED working-tree edits.**
At committed HEAD `3e4fccb` the file had (a) the `streaming-action-pure-l`
math block + a `.. todo:: Archivist expansion needed (#257 S8b)` stub, and
(b) NO `apply-solve-asymmetry` section. The working tree (`git status` =
`M docs/theory/operator_algebra.rst`) had the stub EXPANDED into full
narrative ("Why the matvec is affine in σ" + "Probe evidence and the retired
fold") AND a brand-new ~285-line `apply-solve-asymmetry` H1 section inserted
between the stub and "Capability set semantics". So the work was drafted in
a prior turn of the same session. **Diagnose this FIRST** — `git show
HEAD:<file> | grep <new-anchor>` returning 0 + `git status` = `M` proves the
content is working-tree-only, NOT a re-authoring need. The real task was the
code-verification pass the brief explicitly demanded, and confirming the
`-W` build stays clean. Do NOT rewrite a complete, correct draft.

**The code-verification recipe (every brief fact checked against
`orpheus/sn/operator.py` + `orpheus/sn/spatial/diamond.py` at HEAD):**
- `StreamingOperator.capabilities = frozenset({CAP_APPLY, CAP_APPLY_TRANSPOSE})`
  — **no `CAP_SOLVE`** (operator.py:316). The draft's "`L.solve` is not a
  live call" `.. note::` is CORRECT — frame `L⁻¹` as the MATHEMATICAL
  advection inverse ("even if L alone were inverted…"), never as an
  invocable method. ✓
- `CollisionOperator`: `{CAP_APPLY, CAP_APPLY_TRANSPOSE}` always; `CAP_SOLVE`
  iff `min|σ|>0` (multiplier spectrum law `spec(M[σ])=ess range(σ)`);
  `solve(q)=q/σ` element-wise (operator.py:531-539). ✓
- `InvertibleOperator`: `{CAP_APPLY, CAP_APPLY_TRANSPOSE, CAP_SOLVE}`
  (operator.py:682, 766); docstring itself says "no generic (A+B)⁻¹ formula
  — OperatorSum cannot solve; the WDD sweep IS the inverse algorithm (Lewis
  & Miller §3.2; Adams & Larsen §III)". ✓
- Cell-denom proof: `affine_scan_coefficients` → `geometric_streaming_term +
  collision_volume_term = denom`, `inverse_denom = 1.0/denom`
  (diamond.py:617-623); `_cartesian_streaming_diagonal` → `denom = Σ_t +
  Σ_a 2 g_a` with `g_a = |μ_a|/Δ_a` (diamond.py:329-362, so `2 g_axis =
  2|μ_axis|/Δ_axis`). The draft's "divide by the SUM, you must invert the
  sum" cell-shadow proof is exact. ✓
- ⚠ **`StreamingOperator` is NOT pure `Ω·∇` — it carries the curvilinear
  angular redistribution `(1-μ²)/r ∂ψ/∂μ` TOO** (operator.py:224 docstring
  "Pure streaming + angular-redistribution operator L"). The asymmetry
  section's physical-meaning table frames `L = Ω·∇` (the advection picture)
  as a CONCEPTUAL anchor — acceptable, NOT a contradiction, because the
  expanded `streaming-action-pure-l` block (same page) already names the
  redistribution term. Confirm the simplified-table framing is supported by
  a fuller statement elsewhere on the page rather than standing alone.

**Citation discipline ([[feedback_drift_verification]] L3 + AGENT.md
cross-ref reality):** `[LewisMiller1984]_` and `[AdamsLarsen2002]_` are
defined ONLY in `discrete_ordinates.rst` (cross-doc resolution; this
project is NOT `-n` nitpicky but undefined CITATIONS DO warn). `operator_
algebra.rst` already cites `[AdamsLarsen2002]_` (lines 3937/3979) the same
cross-doc way — so the new cites are NOT the first and create NO duplicate-
`.. [Key]` definition. The `[AdamsLarsen2002]_ §II for ρ=c` attribution
matches the bib-entry text verbatim ("§II gives the source-iteration
spectral radius ρ=c"); `[LewisMiller1984]_ §3.2 (sweep) / §4 (source
iteration)` are the standard sections AND match the InvertibleOperator
docstring's own §3.2 attribution. ALWAYS read the bib-entry text and any
existing in-code section attribution to confirm a §-number before trusting
a draft's citation.

**vv-status tagging on the new labels (matches the affine-in-σ / scanmarch
sibling pattern):** 2 CLAIM-bearing algebraic-law labels (`apply-distributes`,
`solve-does-not-distribute`) carry `.. vv-status: <label> documented` — they
are distributivity/non-distributivity IDENTITIES, not solver claims, so
`documented` is correct (lands in the matrix Documented-only bucket). The
other 7 (`apply-solve-within-group-balance`, `-cell-resolvent`,
`-denominator-inequality`, `-neumann-series`, `-neumann-expansion`,
`-parallel-identity`, `-source-iteration-series`) are DERIVATION-decomposition
labels → left UNTAGGED (orphan-list, beside the affine siblings). Grep
`verifies("<label>")` across `tests/` for ALL new labels (0 hits here) to
confirm none is a verifies-target you'd be mis-tagging.

**Build gate:** forced `-E -W` build EXIT=1 with exactly **1** warning — the
pre-existing `orpheus/geometry/mesh.py:Mesh1D.from_geometry` `:paramref:`
ERROR (the standing MAIN-checkout baseline; needs `sphinx-paramlinks`, out
of scope). Count UNCHANGED from baseline; my edit added ZERO new
WARNING/ERROR/CRITICAL. The Python `SyntaxWarning: invalid escape sequence`
lines from test files are import-time warnings (autodoc/nexus collection),
NOT docutils severities and NOT mine — pre-existing, ignore. `sphinx.ext.todo`
+ `todo_include_todos=True` are enabled, so the 3 out-of-scope sibling stubs
(S3b/S5/S6) render as todos without erroring; I left them untouched.

**File marker ladder (recurs, same as #238):** `operator_algebra.rst` uses
`=`/`-`/`~`/`^`. The new `apply-solve-asymmetry` is an H1 `===` section with
`-----` H2 subsections ("The asymmetry", "What each separate inverse would
mean physically", "The crispest proof — the WDD cell denominator", "The
Neumann series…", "Why this is the right architecture, not a limitation") —
correct level. Em-dash in the WDD-cell-denominator H2 title → underline
sized by `len(title)` code-points NOT bytes.
