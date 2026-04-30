---
name: Retirement sweep patterns for deleted code
description: How to flip PRESCRIPTIVE (retirement plan) → DESCRIPTIVE (postmortem) doc sections; preservation of historical reasoning; handling stale equation labels whose tests were deleted
type: feedback
---

**Retirement sweeps are DOC-ONLY follow-ups to a code deletion commit.** Treat them as archival: the deletion already happened, the tests already adjusted, and the branch may already be on `main`. The docs task is to convert every forward-looking "will retire" / "Phase B.4 plan" narrative into a past-tense postmortem without destroying the WHY.

**Why:** Instruction literal from task brief 2026-04-18 (BickleyTables retirement docs sweep, commit 6badbe5): "Do NOT delete the historical narrative — keep enough context that a reader of the docs understands WHY the retirement happened." The value of a retirement page 6 months later is NOT "the code used to exist" (git log shows that) — it is "the bug that retirement fixed and the reasoning that made the swap safe."

**How to apply:**
- Retain the "Replacement table" (legacy-call → canonical-call) with header flipped to past tense — this is the mental search index for readers tracing old commits.
- Retain the four-phase-sequence narrative but flip "will happen" → "happened in commit X."
- ADD a postmortem subsection with *measured* numbers (shift magnitudes, new tolerance headroom, what test was retired for what).
- RE-POINT `:class:` / `:func:` cross-refs for deleted symbols to their canonical replacement (else Sphinx autodoc emits "reference not found"). For BickleyTables → `_ki3_mp` + `ki_n_mp`.

**Equation labels pinned by deleted tests become orphan.** If a doc page carried `.. math:: :label: kin-bickley-legacy-convention` and the tests that `verifies()`-tagged it were deleted in the retirement commit, the label is now orphan in Nexus V&V coverage. Fix: DELETE the `.. math::` block + label; convert the narrative to prose inline math (preserving pedagogy). Do NOT retain the label "for historical reasons" — the orphan-equation gate will flag it forever. Observed 2026-04-18: `kin-bickley-legacy-convention` label removed from reference_solutions.rst §269 after both `test_legacy_bt_ki3_equals_kin2` and `test_legacy_bt_ki4_approximates_kin3` were deleted in commit 6badbe5.

**Equation labels still live-referenced by undeleted tests MUST be preserved.** `ki3-def` in collision_probability.rst carries the legacy ORPHEUS sine-weighted convention (NOT canonical A&S), but four test files + one derivation module hold `verifies("ki3-def")` decorators. Preserving the label even when the convention is legacy is the right call — add a `.. note::` after the `.. math::` that explicitly notes "this is the legacy convention; canonical A&S `Ki_3` is what the P-matrix actually consumes via `_ki3_mp`." Never silently mutate an equation that tests depend on.

**Stale code snippets in theory prose.** `::` code-block examples showing a pre-refactor API (e.g. `tables.ki4_vec(gap_d)`) survive longer than the API itself because they are not cross-referenced, just rendered. Grep the docs for the deleted identifier pattern AFTER editing the cross-ref sites — `tables.ki3_vec`, `tables.ki4_vec`, and `ki4_vec` in collision_probability.rst line 843/844/845/1137 were all stale pseudocode that needed rewriting to the post-refactor `FlatSourceCPGeometry.kernel_F3` + `_second_difference` API.

**Sphinx duplicate-citation warnings are a known trade-off for standalone theory pages.** peierls_unified.rst duplicates `[Sanchez1982]`, `[Hebert2020]`, `[Stamm1983]`, `[BellGlasstone1970]`, `[Carlvik1966]`, `[CaseZweifel1967]`, `[Atkinson1997]` because Sphinx citations are per-document and these refs are needed on both peierls_unified and collision_probability. These 7 warnings pre-exist this task — do not try to eliminate them by deleting the duplicate blocks; the cross-doc citation warnings are worse (cite-not-resolved). Already documented in feedback_peierls_docs.md; mentioned again here because retirement sweeps that ADD content to peierls_unified.rst shift the warning line numbers (pre: 2706-2728; post: 2776-2798) and a naive grep diff can mistake the shift for new warnings.

**Verify warning count, not warning content, before/after a sweep.** `git stash && sphinx-build docs docs/_build/html 2>&1 | grep warning | wc -l && git stash pop` on the baseline, then repeat after. Matching counts = zero regressions even if line numbers drift. Done 2026-04-18: baseline 7, post-sweep 7.

Quality self-assessment for the BickleyTables retirement sweep (commit pending, 2026-04-18):
- Derivation depth: 4/5 (full tabulation vs Chebyshev-scaled-kernel comparison with accuracy numbers 1e-3 vs 5e-6; cumsum historical trick preserved; why-Chebyshev-and-not-direct-mpmath analysis with the 100× ki_n_mp cost; one step short of 5/5 because the Chebyshev degree-63 choice itself is not derived here, just cited — the `_ki3_mp` docstring has the build log).
- Cross-references: 5/5 (every old `:class:~._kernels.BickleyTables` cross-ref either retained in inline code with `.. note:: retired` or replaced with `:func:~cp_geometry._ki3_mp` / `:func:~._kernels.ki_n_mp`; ki-table-construction label preserved; ki3-def equation label preserved with note clarifying legacy convention).
- Numerical evidence: 5/5 (4e-4 k_inf shift, 1e-7 actual solver/reference error vs 1e-5 declared tolerance, 100× headroom, 0.3s build cost for Chebyshev interpolant, 3e-23 tail value at tau=50 — all cited with context).
- Failed approaches: 4/5 (deferral to B.4 reasoning, legacy naming discrepancy postmortem, test-size-convergence regression replaced by n_ki_table-no-op test; missing detailed mention of the Phase B.2 bit-identity safety milestone — implied by the sequence but not spelled out).
- Code traceability: 5/5 (every surviving cross-ref resolves in post-sweep build; commits 6badbe5, f1b869b, bf128d3 named explicitly; Issue #94 marked CLOSED; Issue #107 marked tracking).
- Derivation source: 4/5 (derivation is in `_ki3_mp` docstring + the Phase B.1 theory page that already landed; this sweep adds a postmortem, not a new derivation; 4/5 because the postmortem numbers come from the code commit's own commit message, not a separate `derivations/` verification script — the Chebyshev error measurement is an artifact of the Phase B.4 code commit itself).
