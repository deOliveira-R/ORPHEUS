# Open Fronts Audit — generated 2026-05-03

Compiled from a full read of the prior session transcript
(`/home/vscode/.claude/projects/-workspaces-ORPHEUS/72485894-14b2-41da-8f34-2902915f9aaa.jsonl`,
1266 records, ~4.3 MB), plus current `gh issue list`, `git log
main..HEAD`, branch survey, theory-page stub scan, and agent-memory
review.

Symbol legend per item:
- **CR** = correctness concern (Cardinal Rule 1).
- **CL** = cleanup / hygiene.
- Effort: S (≤1 day) / M (1–3 days) / L (1–2 weeks) / XL (multi-week).

---

## A. Active V&V gaps (concrete numerical disagreement)

- [ ] **A1. 34% Sood Table 10 disagreement** — `PUa-H2O(N)-1-0-SL` registry value vs NM 1980 / Burkart 1976. Three suspect causes: Sood typo (analogous to Eq 28 already caught), case-builder XS/thickness convention drift, or NM 1980 / Burkart 1976 limitation. Surfaced by Wave 2-A; deferred for focused investigation. Per `feedback_fix_bugs_immediately.md` this should be next-in-line. **CR**, M, depends on E. *(transcript line 1451)*
- [ ] **A2. Atalay 1.5 % systematic z_0 gap** — Wave 2-B Atalay 1997 case_method slab/sphere disagree with Atalay Tables 2–3 and Sood `Ua-1-0-SP` by 0.48–7 %. Closeout memo blames "1.5 % systematic z_0 gap." Wave 2-C agent's convention-bridge finding (Sood Σ_s1 vs DS μ̄ semantics) is a probable shared root cause with A1 — investigate jointly. **CR**, M, depends on E.
- [ ] **A3. NM 1980 reflected-slab F_N at 5e-5 (5 × soft of 1e-5 spec)** — Wave 2-A NM Table 2 cases hit 5e-5 abs at F_7 vs Burkart truth where the spec asked 1e-5. Engineering-grade, not research-grade. Investigation: F_N order convergence + integrand quadrature audit. **CR**, M.
- [ ] **A4. Path A.i ↔ Path B flux cross-check at 5e-2 (7 orders soft of 1e-9 spec)** — both KLL Fredholm and F_N projection paths share a log-singular kernel; neither has singularity-aware quadrature. Hardening pathway named (Atkinson product-Nyström or Gauss-Jacobi log-weight). Multi-day work; agents already dispatched on this front per the "let's pursue steps a/b/c" message. **CR**, L, depends on I3.
- [ ] **A5. WM-72 cylinder Sood `Ua-1-0-CY` post-hardening at 3.3e-7** — meets gate but the Phase B1 prototype was at 1 %; 2 SymPy structural errors caught by Wronskian consistency only during hardening. The pattern (prototype 1 % → hardened 3e-7 after typo catches) is now known-replicable but suggests every Phase B prototype deserves an early Wronskian-style structural cross-check. **CL**, S (process change).
- [ ] **A6. Phase 3A slab vacuum convergence** — earlier reported as "5e-4 self-consistency, structurally slower than cylinder." This was partially a misdiagnosis: ERR-034 (missing μ factor in Δx = μ·s) was hiding inside the "slow convergence." Now fixed via Phase-3A-delegates-to-Phase-3B pattern. Verify the convergence floor is truly research-grade post-delegation — currently quoted at 6.5e-9 method-of-images but the original 5e-4 bound was never re-audited. **CR**, S.
- [ ] **A7. CP production solver returns k_eff = k_inf on bare-critical sphere** — Phase A explorer found `solve_cp` returns 2.25 (= k_inf) at all refinements on bare-critical vacuum BC sphere because its monodirectional kernel `F(τ) = exp(−τ)` was designed for PWR pin cells with white BC, not bare vacuum. L1 smoke gate downgraded to "is finite & positive." This is the **production CP solver's fundamental applicability domain** — flagged as informational but it is the original-goal blocker (see G1). **CR (informational)**, depends on G1.
- [ ] **A8. Cylinder rank-1 abstraction has no MR test** — the rank-1 `compute_resolvent_T` works for sphere + cylinder + slab-symmetric, but cylinder has no multi-region cross-check. Hindsight-review qa flagged this as research-grade-blocker for heterogeneous cylinder claims.
- [ ] **A9. Sphere rank-1 V_α2 — sphere V_α2 still shares the integrand origin for one path** despite the hardening pass. QA noted "structurally independent matrix-element vs polar-escape integration paths" but the cross-check remains in the same SymPy module — the structural-independence claim deserves a paranoid re-audit. **CR**, S.

---

## B. Bugs / known correctness issues

- [ ] **B1. ERR-034: Phase 3A slab `_apply_operator_slab` first-leg parametrization** — fixed in commit (Phase 3B work) but `catches("ERR-034")` binding was missing for ~5 commits; fixed inline in `d63f654`. **DONE**.
- [ ] **B2. ERR-035: Phase 3A symmetric rank-1 closure wrong at intermediate α** — fixed by delegating Phase 3A to Phase 3B rank-2 path; -595 LoC. **DONE**.
- [ ] **B3. Sood Eq 28 typo** — caught by SymPy det(M)=0 derivation. Documented in `k_inf_derivations.py:280-330`, `multi_group/k_inf.py:155-190`, `sood_registry/la13511.py:228-235`. **DONE**.
- [ ] **B4. Westfall-Metcalf 1972 Eq 17 typo (`ν_0` → `ν_0²` in numerator)** — caught during cylinder W-M hardening. Documented in commit `6e545e4`. **DONE**.
- [ ] **B5. WM-72 cylinder V_se-cyl.4 omitted middle Fredholm-coupling term** — caught by hardening agent's Wronskian consistency check (our own SymPy module bug, not literature typo). Fixed in commit `6e545e4`. **DONE**, but confirms B5/B6 pattern: every Phase-B prototype's SymPy needs Wronskian-style cross-check.
- [ ] **B6. WM-72 cylinder V_se-cyl.5 had `R·K_1·I_0` instead of `(R/μ)·K_1·I_0`** — same hardening pass as B5. **DONE**.
- [ ] **B7. Suspected Sood Table 10 typo (Pua-H2O reflected-slab)** — see A1 above; investigation pending. **CR**, M.
- [ ] **B8. fn_method `cylinder/__init__.py:14` says "TODO: populate after literature-researcher delivers Westfall-Metcalf"** — but cylinder W-M was implemented in `singular_eigenfunction/`, NOT in `fn_method/cylinder/`. The fn_method/cylinder folder is now an empty stub that should either be deleted or re-purposed (e.g., for a future 2G cylinder F_N à la Siewert-Thomas extension). **CL**, S.
- [ ] **B9. cases/diffusion.py lines 117 + 196 say "TODO: no labels yet — docs/theory/diffusion.rst does not exist"** — orphan TODOs since diffusion theory page is yet to be written (see issue #35). **CL**, blocked by D-block.
- [ ] **B10. Existing pytest skips not on documented stub list** — `tests/derivations/test_peierls_rank_n_class_b_mr_mg.py:600` and `test_sood_registry_wide_bare_critical.py:179, 193` carry `pytest.skip` but the audit hasn't validated they trace to a documented stub registry case. Review for orphan-skip drift. **CL**, S.

---

## C. Architectural improvements

- [ ] **C1. Hindsight architectural review (delivered but not actioned)** — cross-domain-attacker + qa + explorer trio surfaced ChordOracle + BaseAtlas extraction, unified `power_iterate` driver, single `GreensResult` dataclass, Branch-1 SymPy `pass_*` schema standardisation. User said "let's hold back on architectural improvement" while V&V hardens. After A1-A4 close, this work block becomes the next major architectural arc. **CL**, L.
  - C1a. ChordOracle + BaseAtlas extraction (rule-of-3 PASSED, ~1 day; pollinates `peierls_nystrom`).
  - C1b. Unified `power_iterate` driver + single `GreensResult` dataclass (rule-of-6, OVERDUE).
  - C1c. Standardise Branch-1 SymPy `pass_*` schema rank-1 vs rank-2.
- [ ] **C2. ADR-001 — Six divergent `*Result` dataclasses** — `/workspaces/ORPHEUS/.claude/ADR-001-architectural-evaluation.md` lists 5 architectural debts not opened as issues. **CL**, L.
  - C2a. Single `TransportResult` ABC (Finding 1).
  - C2b. `_e3` exact duplication in cp/solver.py (Finding 2).
  - C2c. `solve_*` entry-point signature unification (Finding 3).
  - C2d. Per-solver Mesh wrapper proliferation (Finding 4).
  - C2e. `derivations/` subtree mixes four different things (Finding 5).
- [ ] **C3. Variant α geometry unification — already done**, BUT the family is now in two epistemic states (sphere/slab-asym research-grade; cylinder/slab-sym/hollow/annulus internally-consistent-but-externally-unanchored). Closing the unanchored set is sequenced under A8 + B-block.
- [ ] **C4. Cross-domain frame: fiber bundle abstraction** — full bundle/G-structure (Casimir invariants, isotypic decomposition) for sphere/hollow_sphere + cylinder/annulus pairings is not yet justified (only 2 rank-1 instances). Wait for instance N+1.
- [ ] **C5. Production XS protocol (`Mixture`/`Mesh1D`/`BC`) bridge to derivation registries** — Phase A delivered `mixture_to_arrays` extractor for fn_method's raw-array internals. Other registries (case_method, carlvik_galerkin, singular_eigenfunction) need analogous bridges to make CP-vs-reference cross-checks plug-and-play. **CL**, M.
- [ ] **C6. Per-paper sub-registries** — `kll_registry`, `dahl_sjostrand_registry`, `atalay_registry`, `bis_1976_registry`, `westfall_metcalf_registry` planned (Wave 3); umbrella aggregator deferred to Wave 4. **CL**, M (mostly data-entry).
- [ ] **C7. `git worktree` isolation for parallel sub-agent dispatches** — Wave 2 saw 3 parallel agents independently each create their own feature branch with cross-contaminated commits (Wave-2-A duplicating NM, Wave-2-C drifting onto case_method, etc). Lesson: pin each parallel dispatch to a worktree path or to a specific branch in the brief. Agent already requested `git worktree` adoption explicitly. **CL** (process), S.
- [ ] **C8. Shared `_e3` helper duplication in cp/solver.py** — lifted from C2b for visibility. **CL**, S.
- [ ] **C9. Variant α MR + group iteration replicated 5×** — source-build `scatter_source` einsum + fission contribution; 6 result dataclasses; 4 distinct `_scalar_flux_from_psi` impls. C1b is the explicit consolidation.
- [ ] **C10. Sanchez/Chandrasekhar Green's-function vs MoC distinction** — User asked for clarification on what we built (Variant α = bouncing-trajectory MoC with closed-form BC closure) vs canonical Sanchez/Chandrasekhar Green's function (boundary-to-boundary integral kernel). Awaiting deeper write-up. **CL**, S (write-up only).

---

## D. Documentation debts

(Theory-page TODO counts measured 2026-05-03 via `grep -c "TODO — Archivist" docs/theory/*.rst`.)

- [ ] **D1. fn_method.rst — 32 archivist TODOs**. Stub structure is in place; rich-narrative expansion is the largest single doc debt. **CL**, M.
- [ ] **D2. case_method.rst — 8 archivist TODOs.** **CL**, S.
- [ ] **D3. sood_registry.rst — 7 archivist TODOs**, including the 2G bare-critical stub block awaiting Siewert-Thomas 1986. **CL**, S.
- [ ] **D4. singular_eigenfunction.rst — 5 archivist TODOs** (cylinder W-M derivation, Mitsis-WM Fredholm narrative). **CL**, S.
- [ ] **D5. carlvik_galerkin.rst — 0 marker TODOs but stub-quality**. The Wave 2-C closeout shipped a Sphinx stub but rich narrative not written. **CL**, S.
- [ ] **D6. peierls_nystrom.rst** — multiple deferred-phase notes referencing Phase 5+ (multi-region sphere annuli, 3-surface topologies). **CL** unless fronts re-prioritise. M.
- [ ] **D7. peierls_greens.rst** — 4225 lines after archivist sweep, multi-region deferred notes still flag Phase 1b / Plan-2 / annulus-MR continuations. Audit for stale "future-tense will-deliver" wording per `feedback_phase2c_staleness_sweeps.md`. **CL**, S.
- [ ] **D8. Diffusion theory page (`docs/theory/diffusion.rst`) does not exist** — issue #35; blocks B9 cleanup. **CL**, M.
- [ ] **D9. Fuel theory page (`docs/theory/fuel_behaviour.rst`) does not exist** — issue #40. **CL**, M.
- [ ] **D10. TH module not in Sphinx autodoc** — issue #45. **CL**, S.
- [ ] **D11. Numerics theory chapter (eigenvalue methods)** — issue #76. **CL**, M.
- [ ] **D12. SN documentation gaps** — quadrature-dimension mismatch behaviour (issue #10). **CL**, S.
- [ ] **D13. MC remaining doc gaps** — derivation scripts, autodoc, convergence table (issue #19). **CL**, M.
- [ ] **D14. Sanchez/Chandrasekhar narrative** — explicit user request for write-up of how Variant α relates conceptually to the named Green's-function methods (overlaps C10). **CL**, S.
- [ ] **D15. Three new feedback memos suggested but never institutionalised** — (a) "hybrid State-1A+1B fallback when SymPy chokes," (b) "every Phase-B prototype gets a Wronskian-style structural cross-check," (c) "pin parallel dispatches to a branch / worktree." Save under `feedback_*` if pattern recurs.

---

## E. Registry / case coverage

- [ ] **E1. Sood registry coverage 42 of 75 cases** as of session end. Target was ~75-80; achieved a partial expansion this session. **CL**, M.
- [ ] **E2. 12 stubs in WIDE_SLICE_STUBS** — `PUB_1_0_CY_STUB`, `UD2O_1_0_CY_STUB`, plus the 10 2G bare-critical stubs (`PU_2_0_SL`, `PU_2_0_SP`, `U_2_0_SL`, `U_2_0_SP`, `UAL_2_0_*`, `URRA_2_0_*`, `UD2O_2_0_*`). All await Siewert-Thomas 1986 implementation. **CL**, M.
- [ ] **E3. Reflected-cylinder cases (~3 in Sood)** — none of the new papers cover; need cylinder-specific F_N reflector paper or extension. **CL**, L (literature gap).
- [ ] **E4. 2G cylinder bare-critical (~2-5 cases)** — methodologically hardest, no published 2G cylinder F_N in our local folder. Needs literature pull or original derivation. **CL**, XL.
- [ ] **E5. ISLC (Infinite Slab Lattice Cell) — 4 Sood cases** — formalism in Hébert §3.8.5 (Hebert(2009)Chapter3.pdf, in folder). Needs implementation. **CL**, M.
- [ ] **E6. Anisotropic-scattering cases (~5-8)** — Burkart-Ishiguro-Siewert 1976 covers but its multi-week implementation cost was deemed not justified; Atalay 1997 covers some but at the 1.5 % z_0 gap. **CL**, L.
- [ ] **E7. Per-paper registries (Wave 3 plan)** — kll_registry, dahl_sjostrand_registry (DS-79 11-eigenvalue spectrum tables — unique), atalay_registry, bis_1976_registry (80 reference values data-only), westfall_metcalf_registry. **CL**, M (data-entry-heavy).
- [ ] **E8. KLL Tables III + VII flux-ratio super-grid** — KLL has more (r/r_c) sample points than Sood's 4-point sampling; lift adds research-grade flux V&V anchors for slab + sphere. Sub-task of E7. **CL**, S.
- [ ] **E9. Burkart 1971 cross-check truth table for NM** — Wave 2-A used Burkart 1971 implicitly as truth source; lift the 8 cross-check cases into bis_1976_registry (or a new burkart_1971_registry) for explicit V&V provenance.
- [ ] **E10. Sood-Forster-Parsons 2003 vs 1999 alignment** — 2003 journal version confirmed has the SAME Eq 28 typo and otherwise matches. Cross-reference in the registry notes; meta-extraction memo at `.claude/agent-memory/literature-researcher/sood_2003_vs_1999_extraction.md`. **CL**, S (cross-reference only).
- [ ] **E11. Cache infrastructure** — `sood_registry/cache.py` exists, but per Phase B4 closeout the largest cache wins are sphere F_N high-N convergence sweeps and Variant α fine-quadrature cross-checks; verify these are wired and used by the relevant test suites.

---

## F. Branch / Git hygiene

(Local branches as of 2026-05-03; remote state mirrored where shown.)

- [ ] **F1. Dangling `feature/wave-2a-fn-extensions`** — fully merged (cherry-picked) onto `feature/peierls-greens-cylinder`; HEAD `74c921c`. Safe to delete. **CL**, S.
- [ ] **F2. Dangling `feature/wave-2b-case-method`** — fully merged; HEAD `144044b`. Safe to delete. **CL**, S.
- [ ] **F3. Dangling `feature/wave-2c-carlvik-galerkin`** — superseded by `feature/wave-2c-carlvik-galerkin-clean`; HEAD `b1ebc15` (a Branch-1 SymPy commit that drifted on case_method, hence the "clean" re-do). Safe to delete after audit. **CL**, S.
- [ ] **F4. Dangling `feature/wave-2c-carlvik-galerkin-clean`** — fully merged; HEAD `21200f9`. Safe to delete. **CL**, S.
- [ ] **F5. `feature/peierls-greens-cylinder` 98 commits ahead of main** — merge target for the entire session arc. Per `docs/development.rst`, `git merge --ff-only feature/peierls-greens-cylinder` once V&V cleanup (A-block) closes. **CL**, S.
- [ ] **F6. Old branches still local** — `feature/agentic-behaviour`, `feature/method-implementer-agent`, `feature/quadrature-architecture`, `feature/rank-n-cin-aware-basis`, `feature/rank-n-class-b-mr-mg`, `feature/peierls-greens-function`, `feature/peierls-specular-bc`, `refactor/peierls-facade-narrow`, `investigate/peierls-solver-bugs`, `develop`, `master`. Each needs an audit (merged? superseded? abandoned?). **CL**, S each.
- [ ] **F7. `feature/peierls-specular-bc`** — `4dc03cf [origin/feature/peierls-specular-bc: ahead 3] docs(peierls): Phase 5 RETREAT — continuous-µ K_bc is hypersingular, abandon production wiring`. Closed-out experiment; still local. Decide: delete vs keep as historical reference. **CL**, S.
- [ ] **F8. `develop` and `master` (legacy)** — both are pre-rename branches on early commits. Safe to delete. **CL**, S.
- [ ] **F9. Push-state hygiene** — `restructure/package-layout` and `feature/geometry-module` track origin; check whether origin still needs them after their work merged into main.

---

## G. Production solver work (the original goal)

- [ ] **G1. Production CP solver applicability domain** — Phase A surfaced that production `solve_cp` returns k_eff = k_inf on bare-critical vacuum-BC sphere because its monodirectional `F(τ) = exp(−τ)` kernel was designed for PWR pin cells (white BC). The user's stated goal — *"verify production CP solver against analytical references"* — is currently blocked at the architectural level: the production CP cannot solve the regime that all of our continuous references cover. Decision point: extend production CP to bare-critical vacuum (research add-on) OR build a new production CP variant for vacuum-BC OR re-anchor the V&V on PWR-pin-cell white-BC truth values. **CR**, L (architectural). Likely the highest-leverage front of all.
- [ ] **G2. CP heterogeneous independent reference** — Issue #21 (high-stats MC instead of CP proxy for cross-check). Tied to G1.
- [ ] **G3. CP interface current method for multi-cell lattice** — Issue #56. Direct production extension; ties to G1.
- [ ] **G4. CP ray-tracing for arbitrary 2D geometry** — Issue #55. Production roadmap item.
- [ ] **G5. CP vectorise group loop** — Issue #57. Performance, not correctness.
- [ ] **G6. Sanchez 1982 cylinder reference scatter/fission split** — Issue #144. Tied to disambiguating production-CP-vs-reference comparisons.
- [ ] **G7. Programmatic ingestion of Case-Zweifel 1967 bare-sphere R_c table for Peierls sphere tie-point** — Issue #143. The Peierls sphere reference is currently anchored on PS-1982; tie-point lift would add a second anchor and unblock cross-method audits.
- [ ] **G8. V&V harness retrofit Peierls cylinder/sphere tests at geometry-specific labels** — Issue #142. Wires the new continuous reference family into the harness in the form the production CP cross-checks need.
- [ ] **G9. MG benchmark Peierls multi-group residual against discrete CP** — Issue #140. Direct cross-check from continuous → production; depends on G1.
- [ ] **G10. Ki_1 retirement (issue #116)** and slab_polar retirement (#115) — both are production-side cleanup items waiting on the new continuous references reaching production-grade precision.
- [ ] **G11. Lift SymPy second-difference into `cp_geometry.derive_second_difference()`** — Issue #141. Architectural integration of derivations into production codepath.
- [ ] **G12. Unified flat-source CP in single cp_geometry.py module** — Issue #107. Production-side unification analogous to what we did for Variant α on the reference side.
- [ ] **G13. Production CP solver verified against analytical references** — the original user-stated objective. Currently blocked by G1; once unblocked, all of A-D-E feed into this. **CR**, XL.

---

## H. Tests promised but not written

- [ ] **H1. F_N projection (Path A.i) for sphere** — Wave 2-A only implemented Path A.i for slab; sphere Path A.i deferred. **CL**, M.
- [ ] **H2. Slab-sym vs Nyström at intermediate α** — Nyström supports `{white, vacuum}` BCs only; α-aware Nyström extension would close this. Documented as deferred. **CL**, M.
- [ ] **H3. Hollow / annulus published-reference benchmark** — no external benchmark found this session. Currently relies on R_in→0 reduction + V_α1 closed-shell + intermediate-α self-consistency. **CL**, L.
- [ ] **H4. Cylinder MR test** — qa hindsight flag (A8); cylinder/slab/hollow/annulus all lack MR tests. **CR**, L.
- [ ] **H5. Off-diagonal asymmetric intermediate α suite** — slab-asym 1e-5, hollow 3.5e-4, annulus 1.9e-3 (geometric-class noise). Documented as "geometric-class noise" but worth a tighter audit. **CL**, S.
- [ ] **H6. Slab-sym method-of-images analog** — qa hindsight flag (gap 5). Free given existing slab-asym infrastructure. **CL**, S.
- [ ] **H7. Cylinder external reference cross-check at MG** — single-group anchor exists (Sood `Ua-1-0-CY` 8.5e-6, then 3.3e-7 post-hardening); MG external anchor TBD.
- [ ] **H8. F_N method `verifies()` label coverage on rich-machinery flux tests** — Wave 2 added flux reconstruction but not all flux test labels are bound to `:label:` blocks in `fn_method.rst`. Audit needed. **CL**, S.
- [ ] **H9. ERR-NNN catalog audit** — `error_catalog.md` has ERR-034 + ERR-035 from this session; verify both have `@pytest.mark.catches("ERR-NNN")` bindings in `tests/_harness`. **CL**, S.

---

## I. Literature with unexploited content

- [ ] **I1. Burkart-Ishiguro-Siewert 1976** — 80 reference values across 40 Milne extrapolated endpoints + 40 reflected-slab critical half-thicknesses. Ingest as data-only registry (E7); method implementation deferred (multi-week cost). **CL**, S (data-only).
- [ ] **I2. Dahl-Sjöstrand 1979** — 11 eigenvalues per case (full spectrum, unique among the four method papers); P_1 anisotropic. Carlvik 1968 errata flag noted (Eq 4b sign). Spectrum data is genuinely unique reference material. **CL**, S (data-only) → already partially done in Wave 2-C.
- [ ] **I3. Atkinson 1972 + Atkinson book "The numerical solution of integral equations of the second kind"** — user explicitly asked for literature-researcher → numerics-investigator dispatch on Atkinson's product Nyström for the log-singular kernel hardening (front A4). User requested per session-end message. **CR**, L.
- [ ] **I4. Garcia 2017/2019/2021 spherical harmonics series** — multiple papers in folder (P_N expansion in sphere/spherical-shell, multi-region spherical, neutron-transport stable-solution). Not yet exploited; would unlock ψ(r,μ) angular-flux verification anchor for SN. **CL**, L.
- [ ] **I5. Sanchez 1977 / 2002 / Mazumdar 2015** — CP and MoC method papers in folder, surfaced during the Variant α MoC analog discussion. Untapped for MoC verification. **CL**, M.
- [ ] **I6. Pomraning 1989 / Pomraning-Siewert 1982** — sphere integral-form transport. Pomraning-Siewert is already used by `peierls_nystrom/ps1982_reference.py`; Pomraning 1989 (general geometry) untapped.
- [ ] **I7. Mitsis 1963** — used as substrate for cylinder W-M hardening; supplementary slab/sphere monoenergetic critical solutions also in scope. Partially used.
- [ ] **I8. Halsall 1980 (CACTUS) + Askew 1972 (characteristics formulation)** — historical MoC references; not exploited.
- [ ] **I9. Roy 1989 (3D heterogeneous cells) + Garcia 2006 (CP in r-θ-z)** — multi-D CP method papers; untapped.
- [ ] **I10. Sanchez 1986 (spherical integral form linearly anisotropic)** — already used as Variant α substrate but the linearly-anisotropic extension is unexploited for our isotropic-only Variant α.
- [ ] **I11. Hébert 2009 Ch. 3 + Stammler 1983 Ch. 4 + Ch. 6 + Ligou 1982 Ch. 8 + Stacey 2007 Ch. 9** — textbook chapters in folder; chapter-by-chapter extraction would surface canonical formulations not in the journal-paper corpus. ISLC (E5) sits in Hébert §3.8.5.
- [ ] **I12. Grechanuk 2016 (Semi-Analytical Benchmarks for MCNP6)** — semi-analytical benchmark catalog; unexplored as MC verification anchor.
- [ ] **I13. Makine 2018 (exact transport representations of classical and nonclassical SP_N)** — likely tangential; flag for future SN/diffusion verification.
- [ ] **I14. Sanchez 1986 + Sanchez 2002 (boundary conditions in trajectory-based deterministic transport methods)** — directly relevant to Variant α boundary-closure operator framework; untapped.
- [ ] **I15. Valougeorgis 1985 PhD thesis (F_N method in kinetic theory)** — broader F_N treatment; unexplored.
- [ ] **I16. Morel 1989 (Hybrid Collocation-Galerkin-SN)** — possible cross-pollination with case_method or SN production.
- [ ] **I17. Brockmann 1981 (anisotropic scattering in numerical neutron transport)** — supporting paper for anisotropic-scattering deferred fronts.
- [ ] **I18. Calvik (Carlvik) 1967 finite cylinders/cuboids CP** — substrate for production CP geometry extension (G4).
- [ ] **I19. Benoist 1981 (integral transport theory for diffusion coefficient calculations in Wigner-Seitz cells)** — untapped diffusion-coupling work.
- [ ] **I20. Siewert-Thomas 1986** — 2G slab + sphere F_N. The 10 2G stubs (E2) are blocked on this paper's implementation.

---

## J. Cross-domain insights to revisit

- [ ] **J1. Fiber-bundle abstraction (BaseAtlas, AngularFiber, ChordOracle)** — cross-domain-attacker memo `variant_alpha_2surface_bie_frame.md`. Top match for the 6-geometry family. Sequenced under C1a.
- [ ] **J2. Phase-space pair (sphere, hollow_sphere) and (cylinder, annulus)** sharing same angular fiber over different BaseAtlas SUPPORT — the b-partition is the pullback of the inner-boundary characteristic function through the chord map. Captured but not formalised in code.
- [ ] **J3. Operator-theoretic resolvent `T = (I − S)^{-1}`** — already actioned (Phase 2 unification + Phase 3B rank-2 introduction). The framework correctly predicted rank-2 architecture before any 2-surface code was written.
- [ ] **J4. F_N method elegance — angular pre-integration via Wiener-Hopf X-function** — captured in late-session conceptual narrative. Pollination opportunity to SN (which currently does numerical angular integration; pre-integration could close grazing-ray quadrature issues).
- [ ] **J5. Variant α as conceptually closed-trajectory MoC + analytical BC closure** — user's question; partial answer in session, deserves a Sphinx narrative section (D14). The geometric-series resolvent is the analytical closure for what discrete MoC handles by iterating multiple bounces numerically.
- [ ] **J6. "Sphere = odd-mode slab + antisymmetric BC; one sign flip + one sin↔cos"** — Atalay 1997 confirms this (already in case_method); also confirmed in Siewert-Thomas 1986 for 2G F_N. Generalisable abstraction worth lifting once a third instance appears.
- [ ] **J7. Three-pillar verification structure achieved on Sood `Ua-1-0-SP` and on bare-critical slab** — F_N (existing) vs Carlvik-Galerkin (Wave 2-C) vs Variant α (existing). When this triangulation passes at 1e-5 and disagrees at 0.5 %, the disagreement IS V&V information about the weakest pillar.
- [ ] **J8. "Build new instance standalone first; only unify after ≥2 working instances"** — vindicated by ERR-035; already saved as `feedback_unify_after_two_instances.md`.
- [ ] **J9. cross-domain-attacker `elegance_smell_rank_non_monotone.md`** — dangling memo; subject not actioned this session.
- [ ] **J10. `phase5_continuous_mu_frames.md`** — cross-domain-attacker memo from earlier work; relates to abandoned `feature/peierls-specular-bc` (F7).

---

## K. Quick-wins (1-day or smaller)

- [ ] **K1. Branch deletions** — F1, F2, F3, F4, F8 (S each, batchable).
- [ ] **K2. fn_method/cylinder/__init__.py stub deletion or repurpose** — B8.
- [ ] **K3. Slab-sym method-of-images analog test** — H6 (free given infra).
- [ ] **K4. Sphere V_α2 paranoid re-audit** — A9.
- [ ] **K5. Phase 3A slab vacuum convergence floor re-audit post-delegation** — A6.
- [ ] **K6. ERR-NNN catalog binding audit** — H9.
- [ ] **K7. F_N flux test `verifies()` label binding audit** — H8.
- [ ] **K8. peierls_greens.rst staleness sweep** — D7.
- [ ] **K9. peierls_nystrom.rst Phase 5+ deferred-phase audit** — D6.
- [ ] **K10. Sanchez/Chandrasekhar narrative section** — D14 / C10.
- [ ] **K11. KLL Tables III + VII flux super-grid lift** — E8.
- [ ] **K12. Sood 2003 vs 1999 cross-reference into registry notes** — E10.
- [ ] **K13. Burkart 1971 cross-check truth lift** — E9.

---

## L. Long fronts (multi-week)

- [ ] **L1. G1 — production CP applicability resolution** (the original goal). Architectural decision + extension or rework. Highest leverage. XL.
- [ ] **L2. A4 + I3 — Atkinson product-Nyström hardening of log-singular kernel** for Path A.i + KLL flux paths. Deep-numerics work the user explicitly green-lit. L.
- [ ] **L3. C1 — full architectural refactor pass** (ChordOracle, unified power_iterate, GreensResult dataclass). Sequenced after V&V hardening. L.
- [ ] **L4. C2 — ADR-001 implementation** (six dataclasses, three solver-entry conventions, mesh wrappers, derivations subtree split). XL.
- [ ] **L5. E2 — 2G bare-critical Siewert-Thomas implementation** (10 stubs unlock). M-L.
- [ ] **L6. E5 — ISLC implementation against Hébert §3.8.5** (4 stubs unlock). M.
- [ ] **L7. Hollow/annulus published-reference acquisition** — H3. Likely needs a fresh literature-researcher pass on radiative-transfer / heat-radiation literature where annular geometries are studied. L.
- [ ] **L8. MR for cylinder/slab/hollow/annulus** — H4. Substantial; touches every geometry. L.
- [ ] **L9. Continuous-flux hardening for KLL Path B near boundary singular kernel** — closure check `∫ψdμ = φ` currently 1e-5 (spec was 1e-10). Closing this requires ~1000-point quadrature + singularity-aware scheme. L.
- [ ] **L10. F_N flux Path A.i for sphere** + cross-method agreement triangulation. M-L.
- [ ] **L11. Fully-fledged sub-registries (E7) wired into harness** with per-paper case enumeration + cross-check tests. M.
- [ ] **L12. Anisotropic-scattering case enumeration unlock (E6)** — Burkart-Ishiguro-Siewert 1976 is heavy; Atalay covers some. L.
- [ ] **L13. Garcia 2017/19/21 spherical harmonics integration (I4)** as SN angular-flux verification anchor. L.
- [ ] **L14. F_N sphere "from-first-principles" full SymPy** — original Path-2 commitment was to derive sphere F_N from Wiener-Hopf in SymPy, not from PS-1982 wrapper. Slab achieved this; sphere via the parity-flip works but a pure first-principles SymPy derivation for sphere is still open (despite the agent reporting "real" Path 2). M-L.
- [ ] **L15. Reflected-cylinder cases (E3)** — lit gap. XL (literature acquisition + implementation).
- [ ] **L16. 2G cylinder bare-critical (E4)** — lit gap, methodologically hardest. XL.

---

## Appendix — not-yet-categorised but flagged in transcript

- **AP1.** Three feedback memos requested institutionalisation but not saved: hybrid State-1A+1B SymPy fallback pattern (line 193), Wronskian-style cross-check on every Phase-B prototype (B5–B6 lesson), and parallel-dispatch worktree pinning (Wave 2 lesson). See D15.
- **AP2.** Phase 2 unification noted "method-implementer didn't commit (misread of git protocol — conventional commits on feature branches are routine)." Workflow nudge to feed back to method-implementer's AGENT.md. Issue #148 covers method-implementer skill follow-ups; possibly already captured.
- **AP3.** Phase B3 agent "explicitly didn't commit" (closeout: "Awaiting explicit user instruction to commit") — same misread. Same workflow follow-up as AP2.
- **AP4.** Wave 2 cross-contamination caused docs/index.rst and MEMORY.md merge conflicts; resolved by union but the resolution is worth a 1-paragraph note to method-implementer / archivist about parallel writes to those files.
- **AP5.** "Sphinx -W issue, agent reports pre-existing Sphinx 9 incompatibility in `docs/conf.py`" — flagged in Wave 2-B closeout. Verified independently? Worth a clean-room re-test.
- **AP6.** B1 closeout's hardening-pathways estimate ("2-3 day work" for full Mitsis-WM Fredholm iteration) was actually achieved in a single hardening dispatch; calibration data point for future estimates.
- **AP7.** Production CP cross-check disagreement on bare-critical sphere returned **k_inf** rather than **k_eff** at all refinements (G1 / A7) — deeper symptom analysis would map which production fronts of CP, MoC, SN are similarly bare-critical-blind.
- **AP8.** The 4 dangling Wave-2 branches (F1-F4) suggest a project-level branch-cleanup pass after every parallel-dispatch wave. Could become an automated post-merge hook.

## Counts

- A: 9 V&V gaps (3 CR-tagged blockers; A4 user explicitly green-lit)
- B: 10 bug entries (6 fixed, 4 outstanding)
- C: 10 architectural items (most under hindsight-review block)
- D: 15 doc debts (52 total archivist TODOs across 4 main pages + 4 missing pages)
- E: 11 registry/coverage items
- F: 9 branch/git hygiene items (5 dangling branches confirmed for deletion)
- G: 13 production solver items (G1/G13 are the original-goal blockers)
- H: 9 missing-test items
- I: 20 literature items
- J: 10 cross-domain insights
- K: 13 quick-wins
- L: 16 long fronts
- AP: 8 appendix items

**Total: ~143 distinct items.**

The original user objective ("verify production CP solver against analytical references") sits under G13; it is currently architecturally blocked by G1 (production CP returns k_inf on bare-critical vacuum-BC sphere). Resolving G1 unblocks the entire G-block and turns the rich continuous reference family this session built into the V&V ground the project was designed for.
