---
name: protocol-retired-abc-survives docs
description: #248 sequel to type-retirement-but-concept-survives — a @runtime_checkable PoleAngularClosure Protocol DELETED but its nominal twin PoleAngularClosureBase ABC SURVIVES carrying every accessor/method; repoint is a CLASS-NAME-IN-DOTTED-PATH swap (accessors live on the ABC unchanged), collapse the dual Protocol-AND-ABC duality prose to single-ABC, retitle heading off the type KEEP anchor (:refs auto-pick new title incl cross-doc), the deleted __call__ R-output narrative repoints to cell_contribution+precompute_psi_state. KEY: module NOT automodule'd → brief's "build breaks (dangling xref)" premise is WRONG (renders plain-text, no -W warning) → repoint is a CORRECTNESS gate not a build gate.
type: feedback
---

The sibling of [[feedback_type_retirement_concept_survives]] (S6.4(f)
WavefrontFlux) but the survivor here is a NOMINAL TWIN, not a realization
in two places: #248 deleted the orphaned `@runtime_checkable
PoleAngularClosure` **Protocol** (+ its dead `__call__` bundle / `tau_mm`
arg / orphaned recurrence helpers) made redundant when #236 Phase 2 B2
retyped every consumer onto the `PoleAngularClosureBase` **ABC** and made
the three strategy methods abstract on it. The concept is unchanged; the
contract just lost its structural-typing half. So the repoint is the
SIMPLEST flavour of this family — every accessor (`c_in_per_ordinate`/
`c_out_per_ordinate`/`tau_per_ordinate`) and method (`cell_contribution`/
`precompute_psi_state`) lives on the ABC with the SAME name → the fix is a
**class-name-in-dotted-path swap** (`...PoleAngularClosure.tau_per_ordinate`
→ `...PoleAngularClosureBase.tau_per_ordinate`), bulk via `replace_all` on
each distinct dotted-path fragment. GREP-confirm the accessors are ON the
ABC first (`grep -nE "def (c_in_per_ordinate|...)" the_module.py`) so you
KNOW the swap resolves.

**The brief's "these will break the Sphinx build (dangling xref)" premise
was WRONG and you must verify it, not trust it.** The module
`orpheus.sn.spatial.pole_angular_closure` is NOT `automodule`'d anywhere
(`grep -rn "automodule:: orpheus.sn.spatial" docs/` empty) → EVERY
`:class:`/`:meth:`/`:attr:`/`:func:` to its symbols ALREADY renders
PLAIN-TEXT, and an unresolvable one emits NO `-W` warning in this
non-nitpicky project (AGENT.md "Cross-ref reality": unresolvable code-xref
= plain text, silent). So the build does NOT break on the deleted symbol;
the `-E -W` count stays at the standing baseline (1 = mesh.py paramref).
The repoint is mandatory on **Cardinal-Rule-1 correctness** grounds (a
reader following the dotted path lands on a deleted symbol), NOT because the
build fails. Report this premise-correction explicitly — it reframes the
acceptance gate from "count drops to 0" (impossible, baseline=1) to
"count UNCHANGED + zero NEW lines reference the retired symbol".

**Three narrative sites beyond the mechanical swap (the load-bearing prose):**
1. **The dual-declaration "both-sites mint" passage collapses.** A
   `replace_all` on the class name leaves an INCOHERENT sentence: "declared
   on the `@runtime_checkable` PoleAngularClosureBase **Protocol** **and**
   on the PoleAngularClosureBase **ABC**" — same class named twice, once
   wrongly as "Protocol". Rewrite to single-ABC declaration + a
   parenthetical preserving the WHY (Phase B declared on BOTH to serve
   structural-typing AND nominal-inheritance consumers; #236 B2 retyped,
   #248 deleted the Protocol). Grep the post-swap file for "Protocol"
   adjacent to the survivor class to catch these.
2. **The deleted `__call__` as an equation's R-output.** The Phase C
   `R_{n,i,g}` term was "the `:meth:`PoleAngularClosure.__call__` output".
   `__call__` is GONE → repoint to the live production verb
   `cell_contribution` (consuming the per-level state `precompute_psi_state`
   stamps once per sweep) + a parenthetical tombstone naming the retired
   `__call__` bundle. Same move for a Pattern-2 "single source of truth"
   todo that claimed the kernel routes through `__call__`: repoint the
   invariant to bind `compute_psi_half_per_level` + `precompute_psi_state`
   directly.
3. **Heading retitle off the type.** `PoleAngularClosure (Issue #168 Phase
   B)` → `The pole angular closure (Issue #168 Phase B)`. **KEEP the anchor**
   `sn-pole-angular-closure-protocol` (cross-doc — boundary_conditions.rst
   `:ref:`s it; 4 live source `:ref:` sites + an in-file `:doc:` parenthetical
   point at it). Retitling + keeping the anchor means ALL `:ref:`s (incl. the
   cross-doc one) auto-render the NEW title — PROVE it in the built HTML
   (`grep 'href="...#anchor"...><span...>NEW TITLE</span>'`). Renaming the
   anchor buys nothing and risks 5 dangling refs. Add a succession `.. note::`
   under the heading: Protocol→ABC contract-evolution, what #236 B2 did, what
   #248 deleted, "read 'strategy ABC' where the original said 'strategy
   Protocol'", anchor-retained rationale.

**History literals are intentional + must survive the swap.** Past-tense
mentions naming the retired `@runtime_checkable PoleAngularClosure`
**Protocol** stay as double-backtick literals ``PoleAngularClosure`` (my
new succession note + the collapsed-duality note both do this) — they are
NOT cross-refs, render correctly, and preserve the WHY. The FINAL gate is
two greps: (a) `grep -oE ":(class|meth|attr|func):\`~?...PoleAngularClosure[A-Za-z_.]*\`"
docs/theory/ | grep -v Base | grep -vE "MorelMontry|LegacyTau|BaileyFlat|IdentityAngular"`
EMPTY (all live cross-refs repointed; exclude the sibling strategy classes
which are NOT the deleted Protocol), and (b) `grep "\`\`PoleAngularClosure\`\`"`
returns ONLY the intentional history literals.

**SCOPE WALL — adjacent pre-existing staleness you FLAG not fix.** The same
file carried `:func:`_mm_weighted_angular_recurrence_single_level`` (×3, in
Phase D/F historical sections) — a symbol renamed away at commit 9eb5d4e
(PR-TYPED-6.5), stale LONG before #248, unrelated to the Protocol
retirement, and rendering plain-text (no build break). `git log -S` to
date the staleness; if it predates your issue, FLAG it as out-of-scope
adjacent staleness (its own follow-up), don't fold it into the #248 edit.
The IN-scope `__call__`/`_mm_psi_half_grid_single_level` refs were the ones
INSIDE the heading section the brief targeted — fix those.

**Build env**: this ran on `main` in the MAIN checkout (NOT a worktree) →
venv at the main repo root, no PYTHONPATH gymnastics, baseline=1 (the
mesh.py paramref). `rm -rf docs/_build/<throwaway>` is permission-BLOCKED —
leave the throwaway build dir (gitignored) for the orchestrator. `matrix.rst`
auto-regenerated on build (expected artifact; main agent commits it, do NOT
hand-edit).
