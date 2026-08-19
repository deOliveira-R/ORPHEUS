# Stencil-assembly + DSA roadmap — k-estimator → spatial reification → DSA (#2)

**⏹ CAMPAIGN COMPLETE — TERMINAL (2026-07-27).** All four phases merged to
main: P1 @ `a4952c3` · P2 = the `f6079be…` chain (2026-07-04) · P2.5 = #280
@ `3f0b8c74` (2026-07-13) · **P3 (consistent DSA, #2) ff-merged
`e4c1a81c..37ffc310` (2026-07-27), branch deleted, #2 + #215 CLOSED.** This
file is now ARCHAEOLOGY — the full phase record, the rulings R1–R14 + R5/R6,
and the R6 execution log. The living record of the DSA build is the theory
page `docs/theory/methods/sn/acceleration.rst` (Key Facts + the three
consistency discoveries + Development history) + `error_catalog.md`
ERR-070/ERR-071; the evidence packs are the committed
`.claude/plans/dsa_*.md`; the open follow-ups are GitHub issues
**#312 (LD arm) · #313 (reflective lag mode) · #314 (2-D DSA) · #315
(scheduled full-space inverse) · #316 (frame-backed angular_moment(ℓ))**,
with dispositions on #200 (DSA = the first re-enabled Krylov M; issue stays
open for the block-inverse candidate) and #227 (windowed-sweep ↔ DSA
cross-link). Final gates at merge: full serial tree 6627/0 · Sphinx `-E -W`
0 · audit 0 orphans of 325, ERR 71/71 · pyright ratchet = the accepted #288
floor (`transport: 1`).

**Status (2026-07-26 — superseded by the terminal block above): Phases 1, 2,
AND 2.5 COMPLETE, merged to main (P1 @ `a4952c3`; P2 = the `f6079be…` chain
ff-merged 2026-07-04; Phase 2.5 = #280 CLOSED 2026-07-13 @ `3f0b8c74` — R9's
"keep `sn/spatial` the name" was superseded by execution: the residual package
is now `orpheus/sn/sweep/`, task #54 @ `588f2429`). Phase 3 (DSA #2) ran on
branch `feature/sn-dsa`: the 3-P0 dispatches DONE (four memos, see the
Phase-3 plan-of-record subsection below), the phase plan USER-APPROVED
2026-07-26 after the diffusion-readiness checkpoint.**

Campaign chain: **Phase 1** k-estimator unification (#259 + #291) → **Phase 2** spatial
substrate promotion + the ASSEMBLY third mode (#272 + #158 + #253, user ruling below) →
**Phase 2.5 (per R1)** the #280 orientation×kernel walk unification, its own
branch/campaign → **Phase 3** DSA for SN (#2, folding #200's Krylov posture; wires into
the UNIFIED walk). Each phase is a separate branch + ff-merge, in order; later phases
consume earlier ones.

## Provenance (how this plan was decided)

Follows directly from campaign #290 (merged @ `3a19133`; see
`.claude/plans/archive/diffusion_integration_290.md`). The user chose the sequence
warm-up → spatial extraction → DSA, then supplied the architectural key (below). The
recon facts in this file were verified against the tree at `3a19133`.

## THE USER RULING (2026-07-03) — assembly is the missing third consumption mode

> DD and LD carry only the closure for forward substitution and Krylov application
> because we don't yet have a matrix representation for them; a matrix representation
> demands the DD/LD stencil implementation, just as diffusion's does. Extracting the
> spatial discretization therefore ALSO builds the assembly machinery.

Verified against the tree — the claim is true and *conservative*:

- `orpheus/numerics/operator.py:751` — `as_matrix` is **basis-vector probing** (column j
  = apply to e_j; the "functor out of the operator category"; `max_dimension=4096`).
  NO operator in the codebase has stencil assembly. SN's only matrix route is O(N)
  sweeps — the historical 640 s matvec-twin gate (#236) is this cost made visible.
  Diffusion's `MatrixInverseOperator(FlattenedOperator(A, template))` also
  materializes by probing (tolerable only because 1-D N is small).
- The per-cell stencil coefficients ALREADY EXIST as data: `orpheus/sn/spatial/
  cell_balance.py` (`CellBalanceTerms`, `cell_balance_terms`) for DD;
  `linear_discontinuous.py` (`_LDCellTerms`) + `_ubld.py` (`assemble_ubld`,
  `assemble_inflow_axis`, `per_cell_solve` — LD already assembles per-cell block
  systems) for LD. The closure IS the stencil, consumed in-flight.

**The reframe that makes it safe (load-bearing design rule):** assembly is NOT a new
discretization — it is the third consumption mode of the ONE closure algebra:

| Mode | Consumer | Walk |
|---|---|---|
| `solve` | sweep = forward/triangular substitution | cells in walk order |
| `apply` | Krylov matvec = action | cells in any order |
| `assemble` | **NEW** — emit the SAME coefficients as (row, col, value) | cells, scatter into global sparse |

**Twin-path guardrail (Pattern 2; the Phase-F lesson):** `assemble` MUST consume the
same coefficient source (`CellBalanceTerms` / `_LDCellTerms` / closure config) that
`solve`/`apply` consume. NEVER a parallel stencil implementation. Pinned by permanent
equivalence gates: `assembled @ x ≡ apply(x)` (ULP-tight) and per-ordinate
`triangular_solve(assembled, q) ≡ sweep solve` on Cartesian.

**What one assembly layer buys (why this ordering is right):**
1. Sparse matrix representation for SN L(+C) — replaces O(N)-probe twins in tests
   (perf), enables sparse-direct/preconditioner consumers (#200), spectral analysis.
2. Diffusion assembly upgrades from probed-dense to direct-sparse (same machinery;
   optional sparse-LU resolvent).
3. **#284 becomes provable at the object level**: sweep solve ≡ triangular solve of
   the assembled matrix in walk order (Cartesian: expected EXACT). The sphere's #282
   lagged pole seed becomes VISIBLE as non-triangularity of the assembled operator in
   walk order — a characterization gate, not a fix, in this campaign.
4. **Consistent DSA by derivation**: the low-order operator can be DERIVED as the
   discrete moment reduction of the ASSEMBLED SN operator (Alcouffe's consistency,
   computed rather than hand-imposed) — Phase 3's foundation.

## Execution mode + compaction protocol

**Surgical mode throughout** (main agent writes; NO `method-implementer`; other agents
dispatched freely — the #290 pattern). **R6 RULED**: Phase 2a's mechanical relocation
is ALSO full-surgical — the delegation exception is struck.

Compaction points marked ⏸ Cn. Before each: (1) commit everything; (2) append a
`## CHECKPOINT` / per-phase `**STATUS:**` block here (done @ hash / deviations /
rulings); (3) tell the user it is safe to `/compact`. Re-anchor after compaction from
THIS file + `git log` (trust git over any summary — process-discipline rule).

Canonical gates every phase (the #290 set): serial walls
`.venv/bin/python -O -m pytest <paths> -m "not slow" -p no:xdist -p no:cacheprovider
--timeout=300`; full `-m "not slow"` tree before each merge (0 reds); CLI
`npx pyright orpheus/` = 1 (the #288 residual — diffusion + the six folders at 0);
sphinx `-W` exit 0 (generate_rst only if derivations change); `python -m
tests._harness.audit` exit 0; demo k = 1.022173 where diffusion is touched;
mutation-verify every new gate's teeth (in-process monkeypatch, never git-checkout on
uncommitted files; pytest-plugin sentinels for -O bare-assert hazards — Mode 8/11).

---

## Phase 1 — k-estimator unification (#259 root + #291 symptom) [~1 session]

The diffusion P5 discipline, transplanted while hot: derive the k-update denominator
THROUGH the loss operator (`⟨1,(Aψ)_bulk⟩_V` = absorption + leakage by the column-sum
theorem + telescoping of the conservative divergence) so no term CAN be forgotten —
then unify the fragmented copies.

**Inventory (verified @ `3a19133`) — six `compute_keff`/`compute_production_rate`
sites:**
- `orpheus/diffusion/solver.py:296,306` — the CLEAN pattern (P/⟨1,(Aψ)_bulk⟩;
  IntegratedReactionRate #270 arm). The template.
- `orpheus/sn/solver.py:941,966` — **#291: omits leakage** (biased k-update for
  vacuum-bounded eigenproblems).
- `orpheus/cp/solver.py:689`, `orpheus/moc/core.py:211` — hand-rolled copies.
- `orpheus/numerics/eigenvalue.py:123,163` — the protocol declarations.
- `orpheus/numerics/iteration.py:1217,1233` — per #259, the KEigenvalue
  estimator-injection seam is DEAD in production — retirement candidate (3-search
  audit; rewire-not-delete for behavioral tests).

**P0 findings (explorer, 2026-07-03 — durable distillation; full report ephemeral):**
- Single production consumer of ALL six sites = `power_iteration`
  (`numerics/eigenvalue.py:225,228`). SN's flux trajectory + fixed point are
  estimator-INDEPENDENT (unit-production renormalization cancels the k scaling) ⟹
  #291 shifts the REPORTED k only; `converged()` consumes dk — unaffected.
- NO tight vacuum absolute-k anchor exists in tests/sn; every SN k assertion consumes
  the `compute_keff` functional (no `direct_eigenvalue` use outside
  tests/homogeneous + tests/diffusion — the P1.4 cross-engine gate fills a real gap).
  Movers under the fix: `tests/sn/eigenvalue/test_keff_curvilinear.py:757-792`
  (P1-effect band — re-measure), `:794-827` (ordering, wide margin — re-measure);
  `tests/sn/operators/test_boundary_conditions.py:145-200` STRENGTHENS. All
  reflective anchors unmoved (leakage ≡ 0 ⟹ the fixed functional degenerates to P/A).
- n2n fork → R7: SN numerator +2Σ₂ (`sn/solver.py:959-964`) over Σa incl.
  rowsum(Σ₂) (`data/macro_xs/mixture.py:111`); MoC same + the L0 pin to re-derive
  (`tests/moc/test_verification.py:379-408`); CP/diffusion already operator-
  consistent (CP L1-pinned `tests/cp/test_verification.py:558-580`).
- Seam (site 6): DEAD in production, alive-via-tests (4 constructors; the single
  injection @ `tests/numerics/test_iteration.py:718`; default-formula pins in
  `tests/numerics/test_estimators_as_functionals.py`). Retirement blast radius:
  `iteration.py:192-206,423-472,1063-1069,1097-1098,1129-1136,1217-1237`; docs
  `discrete_ordinates.rst:13385-13394,13762`, `operator_algebra.rst:2042-2043`
  (+~30 render-silent refs — grep pass at close-out).
- ⚠ Load-bearing: SN's outer `power_iteration` iterate MUST stay a bare ndarray —
  `tests/sn/operators/test_fission_kernel_crosscheck.py:369-375` sentinel reds BY
  DESIGN on a typed-carrier retype. Leakage data comes from solver-held state
  (trace of the last inner solve), scale-consistent with the renormalized φ (check
  the renormalize-vs-report order at `eigenvalue.py:224-228`).
- CP: numerator fold mechanical-with-a-transpose ((N,ng) group-LAST arrays); NO
  loss operator exists — denominator stays the explicit additive spelling (#270
  pattern). MoC: substrate work first (no measure object; MOCMesh is not a Mesh).
- Doc drift for close-out: `docs/theory/methods/monte_carlo.rst:281` (phantom MC
  `compute_keff`), `tests/cp/test_verification.py:523-527` (stale docstring).

**Steps:**
1. **P0 dispatch — explorer**: per-site map (who calls each estimator in production;
   which tests pin each; whether SN's L1 vacuum anchors are sensitive to the #291
   bias or reflective-dominated). Output: rewire table + expected-shift table.
   **DONE — findings above.**
2. **Characterize #291 FIRST** (before fixing): run a vacuum-bounded SN anchor, log
   reported-k vs `direct_eigenvalue`-style ground truth (or balance-identity defect).
   Decide bit-identical vs principled-re-baseline per vv-principles BEFORE the carve.

   **STATUS — DONE (2026-07-03, pre-fix @ `d1daaac`).** Harness: converged-fixed-
   point map ratio k* = P(Mφ*)/P(φ*) with M = solve_fixed_source∘F — structurally
   independent of the estimator under test; noise ≤ 2e-11 (second-application
   drift); reflective control bias 1.2e-10 confirms the estimator-independent
   fixed point. Numbers (GL8, tight tols):

   | case | k_reported | k* (true) | bias |
   |---|---|---|---|
   | homog. 2G VACUUM slab w=8, 40c | 1.83767525 | 0.98163269 | **+87.2%** (L/A=0.872) |
   | het vacuum sphere P0 (`test_p1_lowers` cfg) | 0.86484694 | 0.70601977 | +22.5% |
   | het vacuum sphere P1 | 0.85080423 | 0.67876772 | +25.3% |
   | reflective control | 1.87500000 | 1.87500000 | 1.2e-10 ✓ |
   | reflective Σ₂≠0 (R7 exhibit) | 1.92857143 | 2.61278195 | **−26.2%** (zero leakage) |

   Sphere Δ(P0−P1): reported 1.404e-2 → true 2.725e-2 — STAYS inside the
   (1e-3, 5e-2) band ⟹ `test_p1_lowers` survives (docstring measured-values
   re-measured post-fix); the monotone test re-measured post-fix. The R7 exhibit
   confirms the n2n defect numerically at zero leakage (k* = P_fis/(A_xs−E_n2n)
   checks exactly: 0.78/(0.5185−0.2200) = 2.61278).

   **vv decision: PRINCIPLED RE-BASELINE.** The old reported k is a DIFFERENT
   functional from the posed problem's eigenvalue (coincides only when L=0 ∧
   Σ₂=0); the corrected value is pinned by the map-ratio harness now and the
   P1.4 cross-engine gate permanently. **Bit-identity design constraint for the
   carve:** reflective Σ₂=0 anchors (incl. the regression-snapshot suite) stay
   BITWISE unmoved — spell L as a STRUCTURAL zero over reflective faces (only
   vacuum faces contribute outflow) and keep the IRR P/A spellings (the n2n
   term is exact 0.0 on Σ₂=0, so the float arithmetic is unchanged:
   `P/(A + 0.0 − 0.0)` ≡ `P/A`). The reported-k contract change is confined to
   vacuum-bounded and Σ₂≠0 problems. Harness script ephemeral (job tmp); the
   map-ratio spelling graduates into the P1.4 gate.
3. **The carve**: one estimator discipline — production via the typed
   `IntegratedReactionRate`; denominator derived through each method's loss operator
   (per-method wiring, ONE shared spelling where the algebra permits — extract the
   helper only if ≥2 sites are literally isomorphic; no premature abstraction).
   SN gains the leakage term structurally (#291). Per R7 the n2n placement flips
   operator-consistent in SN + MoC (νΣf-only numerator; 2Σ₂ on the loss side; the
   MoC L0 expected value re-derived). CP/MoC rewired to the same discipline (mark
   OPTIONAL — drop to a follow-up if friction; they are outside the six-folder
   focus; per P0, MoC's missing measure substrate makes its drop-branch likely).
4. **Retire** the dead iteration.py estimator-injection seam (#259) — 3-search audit
   done at P0; shape per R8 (hardwire the Rayleigh defaults as plain methods; drop
   kwargs/aliases/default-factories; rewire the `:718` injection test to the
   fixed-point-agreement assertion).
5. **Gates**: per-method balance identity `P/k = absorption + leakage`; SN
   reported-k ≡ eigen-solve k on a vacuum anchor (tolerance from step 2); teeth: a
   leakage-drop mutation must go RED on a vacuum case and stay GREEN on reflective
   (that asymmetry IS the #291 catcher).
6. Close #291 + #259 via trailers; comment #270 (CP/MoC arms status).

**PHASE 1 STATUS — COMPLETE (2026-07-03, branch `refactor/k-estimator-unification`;
merged per the checkpoint below):**
- Commit chain: `f9a7171` (roadmap + R1–R6) → `d1daaac` (P0 findings + R7/R8) →
  `308f92c` (P1.1 characterization) → `1247677` (the carve; **Closes #291**) →
  `40e2528` (R8 seam retirement) → `da9c942` (gate + teeth) → `a3ada82`
  (`partial_current_metric` accessor) → the close-out commit (docs + catalog +
  markers; **Closes #259**).
- Delivered: the ONE estimator law `k = R_νΣf/(R_Σa + L − E_2n)` across
  SN/CP/MoC/diffusion; `AngularBoundaryFlux.net_current(face)` minted (the angular
  sibling of the scalar `J⁺−J⁻` — the SAME family DSA's restriction consumes later);
  SN leakage from the solver-held trace, reflective = STRUCTURAL zero (bitwise
  preservation, pinned by the gate's reflective bit-pin); R7 n2n flip (MoC L0
  re-derived 1.125→1.25); R8 hardwired KEigenvalue estimators (injection kwargs
  gone; loud-TypeError tooth); CP mesh-required (dead fallback retired);
  ERR-064 + ERR-065 catalogued with `catches()` wiring; the permanent cross-engine
  gate (`tests/sn/eigenvalue/test_keff_estimator_gate.py`, reported ≡ map-ratio k*,
  mutation teeth in-file, `verifies("sn-keff-update")`).
- Docs (archivist, reviewed): `sn-keff-update` equation + the `sn-keff-estimator`
  derivation section (divergence-telescoping balance, the R7 fork algebra, the
  characterization table); KEigenvalue injection sections rewritten to R8;
  `moc-keff-update` re-derived; `cp-keff-update` verified + linkage note;
  `monte_carlo.rst` phantom-member drift fixed; sphinx `-W` exit 0; dead-xref
  grep 0.
- Deviations from the plan, both deliberate: (1) CP's IRR-numerator rewire NOT
  done — it would mint an INTERNAL twin of the ⟨Σx,φ⟩ contraction against CP's
  measure-based per-group diagnostics; the R3 friction branch fired and #270
  (commented) remains the substrate follow-up for CP+MoC. (2) MoC got the R7
  formula flip only; IRR substrate blocked on the fiber-bundle architecture ruling.
- Walls: sn eigenvalue+operators + transport + iteration 1006; cp+moc 242;
  numerics 825; the new gate 7/7 incl. teeth; pyright CLI = 1 (the accepted #288
  residual — the `space: AngularTraceSpace` narrowing and the G_s accessor kept
  the ratchet taut); `tests._harness.audit` exit 0 (ERR-064/065 covered); FULL
  serial not-slow tree **5941 passed / 0 failed** (18 skipped, 54 xfailed —
  baseline; 50:33); sphinx `-W` exit 0.

⏸ **C1** (small phase — compaction optional; checkpoint block regardless).

---

## Phase 2 — spatial substrate promotion + the assembly mode (#272 + #158 + #253)

### 2-P0 — dispatches (before any move)
- **explorer**: the promotion split map for `orpheus/sn/spatial/` — for each module,
  method-generic vs SN-sweep-specific (see draft split below), with the full 3-search
  blast radius (graph callers + text grep code/tests/docs + direct constructors).
- **test-architect** (MUST trigger — cross-subsystem carve): the bit-identity gate
  plan for the relocation + the equivalence-gate spec for the assembly mode
  (assembled≡apply, triangular≡sweep, probed≡assembled).
- **cross-domain-attacker**: the assembly framing — global DOF numbering as the flat
  functor; block structure (per-ordinate triangular streaming blocks, moment-coupled
  scattering); whether the (row,col,value) emission wants its own algebraic type.

**2-P0 FINDINGS (2026-07-04, dispatches complete — durable distillation; full reports
ephemeral):**
- **Split VERIFIED as drafted; the carve is import-clean**: all five PROMOTE modules have
  ZERO `orpheus.sn.*` imports (not even TYPE_CHECKING). `pairing.py` = **STAY** (its
  second argument type `PoleAngularClosureBase` is SN angular machinery; zero production
  callers — its docstring's "production calls pass…" is a latent-consumer overclaim).
- `CellVisit` promotes WHOLESALE with scheme.py (its c_in/c_out/tau are plain floats with
  Cartesian-neutral defaults per #236; graph degree 132). Its "SN-specific by design"
  docstring paragraph is wording drift — archivist fixes at 2c, NOT in the bit-identical mv.
- Production imports are direct-module-path ONLY; package-path imports are tests-only
  (13 files); `sn/__init__` re-exports nothing from spatial. Production rewire list:
  `solver.py:50-51,1991,2022`; `augmented_mesh.py:57-64`;
  `loss_representation/__init__.py:138-146,198,2478,2813`; `sweep_graph.py:85-86`;
  stay-side cross-package imports `scan.py:79`, `sweep_cache.py:109`, `pairing.py:49`.
- The `_ubld` moment-helper hop (#245 re-export at `_ubld.py:84-87`): cross-package
  consumers rewire DIRECTLY to `numerics.moment_layout` for
  `AVERAGE_MOMENT`/`face_moment_tail`; `octant_moment_frame_signs` is genuinely `_ubld`'s.
- **#253 census: ~20 open-coded sites** (the issue's 5-site list is stale; its
  `loss_representation.py:802` address is dead). Canonical home
  `numerics/moment_layout.py` has `face_moment_count`; the CELL-count helper
  `moment_count(per_axis, ndim)` does NOT exist — mint at promotion. 19-site cell-count
  migration list + 2 LD face-count bypasses (`linear_discontinuous.py:557,598` hardcode
  `2**(d-1)`); Kronecker STRIDES (`_ubld.py:154`) and derivations twins exempt.
- #158 disposition: DD+LD registry-keyed = satisfied for the landed arms. REMAINING:
  ExponentialCharacteristic; Step (the first real `(False,·)` pairing occupant, needs the
  LMM 1987 Eq. 5.20 citation); the curvilinear-LD arm (`_require_slab` cites #158).
  Comment at 2c with a state refresh (its comment-2 file addresses are dead), keep open.
- Docs blast radius: **169 promote-side xrefs across 7 rst files**. PRE-EXISTING dangles
  found (fold into the 2c sweep): 23 lines → the retired `sn.spatial.boundary_face_flux`
  (6 live-role), ~15 stale `tests/sn/spatial/` path refs (files long since in `sweep/`),
  `cell_balance.py:64-65` stale test pointer (fix at mv).
- **Test-migration ruling (dependency direction decides):** MOVE to
  `tests/transport/spatial/` = the 5 solver-independent unit files (`test_affine_closure`,
  `test_ld_ubld_primitive`, `test_ld_ubld_symbolic`, `test_linear_discontinuous`,
  `test_scheme_reaction_rate_contract` — zero `orpheus.sn` imports post-rewire). STAY =
  the SN-integration surface (`test_ld_slope_frame`, `test_moment_axis_predicates`,
  `test_spatial_moment_field_space` drive `solve_sn`/`SNMesh`; `test_ordinate_scan_reset`,
  `test_pairing_diffusion_limit`). The `tests/sn/sweep/core/` unit tests STAY
  (directory-stamped `cap("sweep_core")` + `_c_surrogate` shared with genuine
  stay-tests) — import rewires only.
- **Test-architect crux (binds 2b):** NO assembly gate is 0-ULP (CSR summation order ≠
  einsum order): G1 `assembled@x ≡ apply(x)` at rtol≈1e-11, het + asymmetric-SigS +
  non-uniform-h configs, non-flat seeded x, never a scalar functional; G2
  `solve_triangular(assembled L+C) ≡ sweep` at rtol — LAPACK's structural independence
  EARNS the #284 discharge its L2 status; the triangularity leg
  `triu(PᵀMP, k=1) == 0` is the one exact structural zero; G3 anti-tautology — once
  `as_matrix` delegates, the oracle pins the RETAINED probing loop vs assembled (exactly
  ONE probed≡assembled pin per family); #282 = a POSITIVE back-edge assertion with an
  actionable message, never xfail; one-source teeth: a monkeypatched sign-flip in the
  SHARED coefficient source must red BOTH the new gates AND the existing sweep suites
  (DD = the higher twin risk — fused scalar kernel, no dense block today; LD's apply
  already consumes `assemble_ubld` blocks). Full spec: test-architect memory L16.
- **Cross-domain-attacker crux (binds 2b):** `assemble` = the SAME additive-monoidal
  functor `Op → Mat` (already named in the `as_matrix` docstring) into a sparse carrier —
  leaves override `assemble()`, composites recurse via the homomorphism laws (Sum→`+`,
  Product→`@`, Scaled→`*`, TensorProduct→`kron`); `as_matrix` delegates to
  `assemble().toarray()` when available. The C-ravel order IS the Kronecker factor order
  — express as the typed product order, don't inherit it silently. Per-ordinate block
  assembly + lift is math-demanded (triangularity is per-octant; DSA's `R·A·P` becomes a
  clean triple product on the angle block). NO new emission type (a COO-builder with laws
  would twin the operator algebra one layer down): scipy.sparse carrier + a thin
  `SparseAssembledOperator(LinearOperator)` wrapper (`FlattenedOperator`'s parallel);
  consume the carrier's ravel / `FlattenedOperator` template — NEVER re-derive the
  local-to-global map (a reified `LocalToGlobalMap` type is DEFERRED until an
  unstructured consumer — MoC rays / DG connectivity). Negative test: with every leaf
  `apply` monkeypatched to raise, `assemble()` must still succeed. Diffusion's
  `LeakageOperator` = the FIRST emitter (small, elliptic, bit-gateable vs the probed
  matrix; makes the abstraction two-consumer-justified from birth).

### 2a — Relocation: `sn/spatial/` → `transport/spatial/` (discharges #272)
Draft split (P0 verifies; verified inventory @ `3a19133`):
- **PROMOTE** (method-generic trial-space/closure layer): `scheme.py`
  (`DiscretizationSchemeBase` + registry), `diamond.py` (`DiamondDifference`),
  `linear_discontinuous.py` (`LinearDiscontinuous`, `_LDCellTerms`), `_ubld.py`
  (per-cell mass/assembly kernels), `cell_balance.py` (`CellBalanceTerms`),
  the moment-count policy (#253 — single-source cell `per_axis^d` / face
  `per_axis^{d-1}` in ONE place at promotion time).
- **STAY in sn/** (sweep-walk + angular machinery, not spatial trial space):
  `pole_angular_closure.py`, `psi_half_angle_seed.py`, `sweep_cache.py`, `scan.py`;
  `pairing.py` per P0's call.
- Bit-identical relocation (git mv + import rewires; NO behavior). Full 3-search
  retirement audit on the old paths; Sphinx xref sweep (unresolved :mod:/:class:
  refs render silently — grep the built HTML).
- #158's ask (cell updates as first-class concrete types) is largely ALREADY TRUE
  (registry-keyed classes exist) — sharpen contracts/docstrings at the new home,
  close or comment #158 honestly per what's left.

**2a STATUS — COMPLETE (2026-07-04):**
- `5b6598f` — the relocation: 5 modules `git mv`'d (rename similarity 92–99% — the
  bit-identity witness), production rewires per the P0 list, the #245 `_ubld`
  moment-helper re-export RETIRED (consumers → `numerics.moment_layout` directly),
  5 solver-independent unit-test files → `tests/transport/spatial/`, 169 docs
  role-refs swapped (history literals preserved), `matrix.rst` regenerated (the
  layer-import tripwire auto-adopted the new package: 310→311 cases, all green —
  the structural proof `transport/spatial` imports cleanly).
- `11f14f0` — #253: `cell_moment_count` minted as `face_moment_count`'s
  codimension-0 sibling (deviation from the issue's `moment_count` strawman:
  family-pattern naming); all ~20 executable spellings routed through it (2 LD
  face-count bypasses included); the `octant_moment_frame_signs` Kronecker STRIDE
  deliberately exempt (layout indexing, not a count); solver.py's function-local
  moment_layout imports hoisted (the late-import rationale died with #245).
- Deviations, both deliberate: (1) test migration ruled by DEPENDENCY DIRECTION,
  not module ownership — only the 5 solver-independent unit files moved; the
  SN-integration surface (`test_ld_slope_frame`, `test_moment_axis_predicates`,
  `test_spatial_moment_field_space`) and the `cap("sweep_core")`-stamped
  `tests/sn/sweep/core/` unit tests stayed with import rewires. (2) The `_ubld`
  re-export retirement went one hop further than the plan's "keep or point
  directly" — every consumer (incl. LD intra-package + 2 tests) now imports from
  the canonical home; the re-export block + its two `__all__` entries deleted.
- Pre-existing findings parked for 2c: the `ld-slab` verifies-marker has no
  matching equation label (build-time info, predates the carve); 23 dangling
  `boundary_face_flux` doc lines + stale `_modules/` HTML for retired modules
  (build-dir archaeology — consider a clean docs/_build rebuild at close-out).
- Walls: ring-1 549; ring-2 2012; layer tripwire 311; pyright CLI = 1 (the #288
  residual) after BOTH commits; sphinx `-W` exit 0; harness audit exit 0; import
  smoke OK. Full `tests/sn` + `tests/diffusion` serial wall at the 2a→2b seam:
  **1996 passed / 0 failed** (5 skipped, 36 xfailed — baseline shapes; 9:00).

### 2b — The assembly mode (the reification)
1. **Numerics home**: an assembled-sparse operator (scipy CSR/COO carrier) in
   `orpheus/numerics/` conforming to `LinearOperator` (apply = sparse matvec;
   `as_matrix` = densify). **R2 RULED**: structural assembly is a separate
   `assemble()` → sparse assembled operator; `as_matrix` keeps its dense
   Mat-functor semantics and DELEGATES to densified assembly when available
   (probing remains the fallback for assembly-less operators; #213 not blocked on).
2. **Global DOF numbering** = the existing flat layouts (`FullField` ravel order,
   `FlattenedOperator`'s template contract, face layouts) — do NOT mint a second
   numbering; the flat functor already exists.
3. **Emitters**: DD via `CellBalanceTerms` → per-cell (row,col,value); LD via
   `_LDCellTerms`/`assemble_ubld` blocks → block scatter (moment-blocked DOFs — the
   #253 single-sourcing is a genuine prerequisite here, not a rider). Scope:
   **Cartesian (slab + 2-D) only**; curvilinear assembly is OUT (blocked on the
   #282 pole-seed structure — see the characterization gate below).
4. **Consumers wired**:
   - SN: per-ordinate `L(+C)` blocks assemble triangular-in-walk-order; scattering
     stays the K_iso/kernel form (its dense kernel exists via `full_scatter_kernel`).
   - Diffusion: `LeakageOperator`/family assembly replaces probing for the matrix
     path — gate bit-identical vs the probed matrix, THEN optionally switch the
     resolvent to sparse LU (perf; principled-equivalent gate).
   - Tests: replace O(N)-probe twin gates with assembled gates where they exist
     (keep ONE probed≡assembled pin per family as the permanent oracle — the
     fuller-view-oracle exception, decided EXPLICITLY per the retirement rule).
5. **Gates (the teeth of the whole campaign)**:
   - `assembled @ x ≡ apply(x)` per scheme × geometry × ng (ULP-tight).
   - Per-ordinate `scipy triangular solve(assembled L+C) ≡ sweep solve` (Cartesian
     DD + LD; bit-or-ULP per test-architect's spec) — **this discharges #284's
     contract question with an object-level proof.**
   - **#282 characterization**: assembling the spherical walk order and showing the
     pole-seed row breaks triangularity = the lag made visible. Record on #282 as
     evidence; the FIX stays in #282/#280's scope.
   - Mutation teeth: a coefficient sign-flip in the emitter must red BOTH the
     equivalence gate and (via the shared source) the sweep gates — proving there is
     ONE source. If only the equivalence gate reds, a twin path exists → stop, fix.

**2b STATUS — COMPLETE (2026-07-04):**
- Commit chain: `83a0db7` (numerics home) → `3238651` (diffusion emitters) →
  `ba26644` (SN closure-generic walk + gates). Design contract:
  `.claude/plans/archive/assembly_mode_crosswalk.md` (the L16-cited convention file).
- `83a0db7`: `SparseAssembledOperator` (CSR carrier, flat-1-D apply, exact
  transpose, densifying `as_matrix`, idempotent `assemble`; COO duplicate-sum =
  the FEM scatter); the assembly axis joins inverse/adjoint as the third
  three-layer surface (`is_assemblable` / `SupportsAssembly` / `assemblable()`;
  `MissingAssembly` refusal); composer laws Sum→`+` Product→`@` Scaled→`*`;
  **R2 realized**: the probing loop extracted byte-identical to
  `_as_matrix_by_probing` (the retained fuller-view pathway), `as_matrix`
  delegates to densified assembly when assemblable; `FlattenedOperator`
  passthrough. Cast ledger +1 (the scipy `.shape → None` inline-typing gap,
  one documented boundary freeze).
- `3238651`: every diffusion loss leaf emits (L direct-structural from the
  SAME precomputed conductance/closure attributes; B via law-probing; C the
  coefficient array bit-exact; S/N2N by group-impulse probing through the
  production einsum kernels — `dense_per_material` stays the untouched L11
  oracle); the production resolvent materializes THROUGH assembly (Mode-11
  dimension-discriminating sentinel: law-sized probes expected,
  composite-sized forbidden); trace DOF numbering single-spelled
  (`_trace_dof_columns`). Probed≡assembled measured BIT-IDENTICAL on the het
  fixture (max |Δ| = 0.0, nnz 50 = 50); gates at the ruled nulp/rtol anyway.
- `ba26644`: the SN per-ordinate walk assembler — ONE generic body for DD and
  LD (**coefficient extraction by unit probes of the production
  `residual_kernel_batch`**, valid by `is_linear`; zero stencil spelling to
  drift ⟹ the one-source teeth hold by construction); the #253 moment counts
  as DOF strides; `octant_moment_frame_signs` conjugates sweep→global frame
  (the LD negative-octant lesson — raw sweep-frame blocks are wrong by slope
  signs); `ordinate_walk_order` = the gates' P.
- Gate results: **G2 = the #284 discharge** — LAPACK forward substitution ≡
  the production sweep at ~6e-16 (DD slab + 2-D; LD via LU), triangularity
  EXACT; **#282 characterized** — the spherical back edge positively asserted
  (probe-matrix `triu ≠ 0` in sweep order, loud-flip message; cylindrical
  control exactly triangular); one-source teeth green on the 2-D fixture
  (all three modes move together, equivalences persist).
- Deviations, all deliberate: (1) DD+LD landed as ONE closure-generic walk
  (plan had them as separate emitters) — the batched-kernel contract made the
  block walk subsume DD as `cm = 1`; (2) the sparse-LU resolvent switch NOT
  taken (the roadmap's \"optionally\"; probing→assembled-LU already landed via
  the delegation); (3) no O(N)-probe twin test gates existed to replace (the
  640s gate was already retired) — nothing to swap; (4) SN gates live at
  `tests/sn/sweep/test_assembly_mode.py` (dependency direction — they drive
  SNMesh/solve; supersedes L16's `tests/transport/spatial/` placement note),
  diffusion gates beside the stencil gate they extend.
- Parked for later phases: TensorProduct→kron law (no 2b consumer);
  per-octant batched emission (a vectorization for the first large-scale
  consumer); `Identity`/`Zero`/`Diagonal` leaf assembly (defer-until-consumer);
  the #242 dual-form diagonal seam noted in the teeth docstring.

### 2c — Docs + close-out
- Theory: `loss_representations.rst` / `discrete_ordinates.rst` — the three-modes
  section (solve/apply/assemble as consumption modes of one closure algebra);
  `operator_algebra.rst` — `as_matrix` probing vs structural assembly; dev-history
  rows. Close/comment #272, #253, #158, #284; comment #282 with the
  characterization; comment #200 (blocks now assemblable).

**2c STATUS — COMPLETE (2026-07-04):**
- `e066297` — the archivist-authored (main-reviewed) theory record: the
  three-modes section (`loss_representations.rst`, incl. the reconciliation of
  its own 'deferred extension point' prose — it had anticipated exactly this
  mode), the assembly axis (`operator_algebra.rst`), dev-history rows; the
  parked 2a dangles repaired (6 dead `boundary_face_flux` ROLES → literals,
  5 stale test paths); the CellVisit "SN-specific" drift + the pairing.py
  latent-consumer overclaim fixed; the `ld-slab` phantom verifies-marker
  resolved by per-test repoint to the real `ld-ubld-*` equations (ALSO
  covered two former orphan equations); `matrix.rst` regenerated.
- Issue dispositions: **#284 CLOSED** (the object-level discharge — comment
  records the G2 evidence + the honest non-source-rhs scope note); **#282
  commented** (the back-edge characterization + the loud-flip gate contract);
  **#200 commented** (blocks assemblable — the preconditioner's structural
  blocker removed); **#158 commented** (state refresh: new paths, the
  assembly-for-free property of the scheme contract, the three remaining
  arms). **#272 + #253 close via this commit's trailers** (fire at push —
  pushes held).
- Parked (pre-existing hygiene, out of campaign scope): the
  `ld-cartesian-1d` sibling phantom verifies-marker (6 tests, no equation
  home — repoint per-test like `ld-slab`, or mint the 1-D LD umbrella
  equation); a clean `docs/_build` wipe-rebuild for the stale `_modules/`
  pages of long-retired modules (build-dir archaeology only).
- Walls at close-out: sphinx `-W` exit 0 (×2 — archivist + main); harness
  audit exit 0; broad seam wall 3245/0; FULL not-slow tree + ff-merge gate
  recorded below at the merge.

⏸ **C2**.

---

## Phase 2.5 — the #280 orientation×kernel walk unification (own branch per R1)

**Branch `refactor/sn-walk-unification` off main @ `28dbaee` (opened 2026-07-04).
Scope = GitHub #280; the issue COMMENT ("Redesign onto the landed #226 taxonomy")
is the AUTHORITATIVE spelling — the a3 plan files' operator vocabulary is DEAD
(reconciliation below). Gate chain = `a3_solve_transpose_verification.md` §§9–15
(the 2026-07-04 extension, landed with this section). Full surgical mode.**

**The shape (#280 body, unchanged by the redesign):** the coherence axis is
ORIENTATION (fwd↔adj), NOT kernel; execution {scan(solve) / cell-loop(apply) /
wavefront(multi-D)} is a non-free third axis determined by (kernel,
dimensionality) — scan/loop/wavefront stay DISTINCT. Two frames, each shared
across orientation: (1) the **apply-loop frame** `{_apply_walk (fwd),
loss_action_transpose (adj)}` → ONE orientation-parametrized per-cell loop over
the one DAG; (2) the **solve-scan frame** `{_run (fwd), sweep_transpose (adj)}` —
the REVERSE-SCAN coherent with `_run`'s Blelloch scan, NOT a reverse-loop bolted
onto `loss_action_transpose`.

**The post-#226 surface (redesign comment, P0-verified vs the tree):** no CAP
tags exist to flip; the live spelling of the deliverable is —
- `sweep_transpose` lands as **`SweepOperator.apply_transpose`** — the EUCLIDEAN
  reverse-scan `(L+C)ᵀx = b` on the inverse-family sibling
  (`orpheus/sn/operators/sweep_operator.py`; `SweepOperator = (L+C).inverse()`).
- `SweepOperator.is_adjointable` flips True; `_AdjointOperator.inverse() =
  inner.inverse().H` + `is_invertible` make the swap law `A.H.inverse() ≡
  A.inverse().H` an identity of the algebra. Predicate honesty: the comment
  recommends **(b) with (a)'s spelling** — land predicate + method together with
  the sweep arm; predicate = `inner.is_invertible and
  inner.inverse().is_adjointable` (→ ruling R11).
- The METRIC adjoint-solve `A.inverse().H.apply(b) = G⁺·apply_transpose(G·b)`
  falls out of the EXISTING `_AdjointOperator.apply` for free — the a3
  "Deliverable 3" (`_AdjointOperator.solve` + CAP_SOLVE) DISSOLVES. Metric code
  never enters the sweep.
- Sibling scope: ONLY the Sweep arm this issue; GreenOperator (transposed
  splitting) / MatrixInverseOperator (`lu_solve(trans=1)`) / generic
  InverseOperator transposes defer until consumers.

### 2.5-P0 FINDINGS (2026-07-04, explorer + test-architect — durable
distillation; full memos in agent memory)

- **Inference (a) CONFIRMED + SHARPENED:** the 1-D adjoint
  `loss_action_transpose` is DD/scalar-only (buffers `(ng,N,nx)` with no moment
  tail; no `_reframe`/frame-signs; a hand-transposed DD diamond march) — AND the
  gap is **UNGUARDED**: the Protocol promises a loud NotImplementedError
  (LR:344-346) but `CumprodScan.supports` admits LD slab,
  `StreamingOperator.is_adjointable` is unconditionally True, and an LD-slab
  `.H.apply` broadcast-crashes or silently mis-computes on shape coincidence.
  S0 lands the loud guard (fix-now); the unified frame then replaces it
  STRUCTURALLY (kernel-PAIR registration: DD registers fwd+transpose, LD
  forward-only → typed deferral). CLOSING is rejected — the UBLD Schur VJP is
  new kernel math with zero consumers (A4 needs slab/sphere DD).
- **Inference (b) CONFIRMED:** `geom.chain_idx` ≡ `dag_walk_cell_indices` BY
  CONSTRUCTION — both materialize `dag_walk`'s order (sweep_cache.py:260-270
  iterates dag_walk; both iterators resolve `_representative_ordinate`;
  representative-invariance test-pinned). Residual gap: no direct pin on the
  `dag_walk_cell_indices` twin (zero test refs) — S0 adds the one-liner gate.
- **Cell-ordering materializations now SIX** (the map counted four): + the
  assembly `ordinate_walk_order` + the test-local `_curvilinear_sweep_order` —
  strengthens the map's Pattern-2 consolidation argument for the one-DAG walk.
- **The existing matvec `array_equal` canaries are SELF-REFERENTIAL**
  (removal-form gates compare against the SAME `loss_action[_transpose]` body on
  a fresh instance — override-not-leak discriminators, blind to a relocation
  that moves both paths together) ⟹ 2.5a's bit-identity claim rides FROZEN
  pre-carve `assert_regression --capture-baseline` snapshots (fwd + adj,
  slab/sphere/cyl, 2G, `_random_composite`). And the 1-D walk has NO one-walk
  spy today — 2.5a mints `tests/sn/sweep/core/test_one_dim_loop_walk.py`
  (wrap-spy proving BOTH orientations execute the ONE frame + an AST tripwire
  banning `is_adjoint`/`is_forward`/`is_transpose`/`is_reverse` identifiers —
  orientation is an OBJECT, the `_SweepEmit` discipline's sibling).
- **Assembly gives a NEW transpose oracle (Cartesian only):**
  `solve_triangular(assembled.T, b[order], lower=False)` — LAPACK back-sub,
  independent of the reverse-walk AND of the ORPHEUS scan (catches a wrong
  transposed scan coefficient a'/b'). Instantiates on the DD SLAB this phase
  (2.5's transpose scope is 1-D DD); the SPHERE keystone stays G2's dense
  forward-apply oracle (assembly refuses curvilinear until #282's structure
  fix).
- **#282 entanglement split (→ ruling R10):** the in-pass DIRECT
  starting-direction solve — route (a)'s dynamics; the direct solver
  `carlson_inward_sweep_from_source` EXISTS (psi_half_angle_seed.py:433);
  plumbing gap = `CarlsonSweepContext` carries no source field — touches exactly
  the bodies 2.5 rebuilds (`_run` seed block, `_apply_walk` seed consumption,
  the transpose seed adjoint, the two moment_frame guard rationales that cite
  the lag) = rides the carve. The literal CARRIER augmentation (per-level
  ψ(·,μ=−1) DOFs on FullField/to_flat/metrics + every flat consumer) = a genuine
  scope extension OUTSIDE the walk modules. Sequencing fact: `sweep_transpose`
  on sphere must reverse WHATEVER seed treatment the forward has — landing the
  directness first/together avoids building the reverse-scan against the lagged
  formulation twice.
- **R9 correction (the estate):** R9 keeps the `sn/spatial/` NAME; 2.5 owns the
  LAYOUT DECISION. Priced: 17 test files import `orpheus.sn.spatial.*`; ~150
  docs refs (grep-driven migration — the build does not warn on dangled roles);
  `pairing.py` has ZERO production call sites (test/V&V-facing predicate);
  hidden coupling: `_run` stashes `mesh._geom_cache`/`_coll_cache`.
- Doc-fix riders for 2.5e: `pole_angular_closure.py:963-971` stale default-seed
  docstring (says CarlsonInwardSweep; C5 default is AngularEdgeExtrapolation);
  issue-body line drifts (LR:3197→:3372, seed solver :428→:433).

### Steps

- **S0 — pre-carve scaffold:** G3 full-loss `(L+C−S−B)` G-reciprocity (S† live
  @ `15185e5`; asymmetric SigS ≥2G real mixture, per-group one-hot φ per vv L27,
  slab + sphere — extends `test_g_adjoint_reciprocity.py`, the composite
  adjoint-matvec canary hardening the surface 2.5a rebuilds); FROZEN fwd+adj
  matvec baselines (`--capture-baseline`, slab/sphere/cyl 2G); the LD-slab
  transpose loud guard (the FLAG-2 fix-now); the `dag_walk_cell_indices ≡
  dag_walk` pin.

  **S0 STATUS — COMPLETE (2026-07-04): `368cbbe` → `d9a6881` → `f869806`.**
  - `368cbbe` (S0.1, the FLAG-2 fix-now): scheme trait `has_transpose_kernel`
    (base False / DD True / LD False-with-citation);
    `StreamingOperator.is_adjointable` scheme-honest (a∧b propagates; eager
    `.H` on LD → MissingAdjoint at construction); the reverse-walk entry
    guard backstops direct Euclidean calls; the stale "S/F advertise no
    apply_transpose" domain paragraph retired; gates =
    `test_ld_adjoint_deferral.py` (trait pins / predicate flips / eager
    refusal / guard raises / DD positive controls).
  - `d9a6881` (S0.2): `ScatteringOperator.apply_transpose` gained the
    composite FullField arm (mirrors the forward lift — bulk cotangent only,
    implicit-zero trace; +1 documented cast, the #257 S8c runtime-truth
    precedent; overloads per the MultiplicationOperator pattern) — WITHOUT it
    the full-loss `.H` was predicate-reachable but crashed at the S leaf. G3
    rows: full-loss G-reciprocity rel<1e-12 on het 2-material slab+sphere ×
    {2G asym P0/P1/Sig2, 4G asym P0}; per-group one-hot rows (vv L27); the
    S-transpose-drop tooth (O(1) red). Mode-12 honesty in-file: a
    posing-sign (+S/−S) mutation is INVISIBLE to reciprocity by construction
    — never credited to this gate.
  - `f869806` (S0.3): the frozen pre-carve baselines (`walk_matvec_{slab,
    sphere,cyl}_2g.npz`, both orientations, nulp=1 + DriftWarning tripwire;
    captured + re-verified 0-ULP; curvilinear rows re-capture at 2.5d per
    R10) + the `dag_walk_cell_indices ≡ dag_walk` twin pin (inference (b)
    closed as a GATE — zero prior test refs to the twin).
  - Walls: S0 closing wall (sn operators+sweep+regression, transport/spatial,
    numerics) **2245/0** (~4 min serial); pyright CLI = 1 throughout (two
    transient overload/typing errors caught by the CLI mid-S0.2 and fixed
    pre-commit).
  - Deferred DELIBERATELY to 2.5a: the multi-D predicate-honesty question
    (`is_adjointable` stays True for DD multi-D whose transpose raises — the
    loud-raise deferral; the kernel-pair/orientation-frame registration is
    the structural place to key the predicate on dimensionality too).
- **2.5a — the apply-loop frame** (bit-identical BOTH orientations): ONE
  orientation-parametrized per-cell loop over the one DAG (orientation carries:
  cell order fwd/reversed, boundary in↔out swap, Carlson mirror routing, the
  angular_adjoint second factor). Proof = frozen-baseline array_equal both
  orientations + the new spy/AST-tripwire pins + all forward canaries green.

  **2.5a STATUS — COMPLETE (2026-07-04): `b1a4c78` → `1cf07ed` → `c439ed1`.**
  - `b1a4c78` (the fold — bit-identical BOTH orientations): the frame =
    `_WalkLeg` (one (μ-level × direction) chain: within/ordinates/abs_mu
    bundles + cells in traversal order) + `_dag_legs()` (DEPENDENCY order —
    all −1 legs then all +1; the ±eps masks + `dag_walk_cell_indices` order
    materialize ONCE for both orientations; `_MU_DIRECTION_EPS` single-
    sources the direction trichotomy) + `_reverse_traversal()` (exact
    reverse-mode order; the pole edge reverses with it) +
    `_OneDimScanWalk._loop_walk(legs, *, open_leg, visit, close_leg)` (the
    `_OctantWalk._interior_walk` sibling) + `_degenerate_positions()`.
    Kernels stayed VERBATIM closures. The adjoint's leg schedule re-nested
    per-level (+1,−1) → the exact reverse — value-identical (disjoint leg
    slots; level-local mirror handoff), held 0-ULP on the frozen baselines
    under `-W error::DriftWarning`. New pins
    `tests/sn/sweep/core/test_one_dim_loop_walk.py`: wrap-spy (both
    orientations execute the ONE `_loop_walk`, slab+sphere), AST tripwire
    (octant ∪ orientation smell sets), the reversal law, the dependency
    order, the leg/degenerate trichotomy partition. Teeth dev-verified
    in-process (monkeypatch): M1 cells-unreversed → adj moves O(10) all
    geometries, fwd untouched; M2 legs-unreversed → slab adj BIT-IDENTICAL
    while sphere/cyl move O(10) — the geometry-selective pole-handoff
    discriminator, exactly the DAG analysis's prediction.
  - `1cf07ed` (**ERR-066** fix-now, the FLAG-2 silent-garbage class, found
    while carving): `loss_action_transpose` had NO degenerate volumetric
    branch — `.H.apply` on any 4|n_phi product-quad cylinder (the DEFAULT
    `product(8,8)` included; φ = π/2, 3π/2 ⇒ |μ_x| ≈ 6e-17) silently
    DROPPED those rows; every cylinder gate rode `level_symmetric` (empty
    degenerate class — vv Mode 7 in the QUADRATURE dimension). Fix = the
    degenerate reversal block in the frame (a `_degenerate_positions`
    consumer; zero participation on empty class → all baselines bitwise).
    Gate = the `cyl_product_2g` G-reciprocity row (`_make_cyl_product`;
    RED O(1) pre-fix → GREEN rel<1e-12 post-fix, dev-verified in that
    order) + `catches("ERR-066")`; catalog entry appended.
  - `c439ed1` (the S0-deferred multi-D predicate honesty, resolved in the
    frame design): `is_adjointable = scheme.has_transpose_kernel ∧
    rep.has_transpose_walk` — the predicate FACTORIZES along the #280
    kernel×orientation axes. New protocol property `has_transpose_walk`
    (base opt-in False; CumprodScan True; ScanMarch mesh-keyed 1-D/2-D;
    the DAG-wavefront family inherits False). Eager `.H` on DD cart2d now
    raises MissingAdjoint at CONSTRUCTION; the representations' loud
    raises stay the direct-Euclidean-call backstop (the S0.1 layering).
    `TestMultiDOrientationHonesty` pins; deferral file 17/17.
  - Walls: sn+transport+numerics **3170/0** (~9 min serial); pyright CLI =
    1 throughout (the accepted #288 residual). Scope note: the angular
    SECOND factor stays method-level in both orientations (forward:
    `psi_state` precompute + in-visit `cell_contribution`; adjoint:
    in-visit `numer_bar` accumulation + ONE `angular_adjoint` close) — the
    loop frame unifies the SPATIAL march; the closure delegation was
    already single-sourced (never re-inlined, ERR-058/#195).
  - NEXT per R10: **2.5d** — the #282 direct starting-direction seed +
    carrier augmentation, landing ONCE on this unified frame; the
    curvilinear baseline rows re-capture there (Cartesian stays bitwise).
- **2.5b — `sweep_transpose` as the REVERSE-SCAN** (re-baseline, NOT bit-id):
  the transposed affine recurrence over the reversed chain, coherent with
  `_run`. Gates: G1 round-trip + G2 dense-`Mᵀ` (SPHERE keystone; iterate-thread
  the Carlson seed) + the assembled-`Mᵀ` LAPACK oracle (DD slab) + mutations
  (forward-DAG order; ±μ mirror on sphere; wrong a'/b'; the ×V/÷V two-denom
  seam). Scope: 1-D DD all geometries; LD-slab = typed deferral; multi-D +
  `ScheduledInvertibleOperator` (schedule-folded reflect-transpose) = explicit
  defers.
- **2.5c — the inverse-adjoint wiring** (the redesign comment's spec):
  `SweepOperator.apply_transpose` exposure + `is_adjointable` True;
  `_AdjointOperator.inverse()` + `is_invertible` (per R11). Gates =
  `tests/sn/operators/test_inverse_adjoint_coherence.py` (comment Gates 1–3:
  the forward-matvec G-reciprocity pin `⟨A ψ, x⟩_G = ⟨ψ, b⟩_G` for
  `x = A.H.inverse().apply(b)` — never calls the transpose path; the swap-law
  value gate rtol 1e-12; the `A.H.apply∘A.H.inverse ≈ I` round-trip) + Mode-11
  wrap sentinel (the reverse scan EXECUTES; forward solve counter 0) + M-ADJ-swap
  / M-ADJ-metric mutations (metric one reds on sphere/cyl ONLY — the
  `.H`≠Euclidean discriminator) + predicate flips (InverseOperator/GreenOperator
  STAY non-adjointable) + the `assert_type(A.H.inverse(), LinearOperator[D,C])`
  static pin.
- **2.5d — #282 per ruling R10.** If directness rides: solve/apply/transpose all
  adopt the direct starting-direction treatment (principled re-baseline where
  the fixed point moves); the characterization gate flips LOUD by design → its
  successor per the conditional spec (spherical G2: triangularity + LAPACK ≡
  sweep over the [seed-rows-first, ordinate-blocks] order; teeth = coupling-
  direction swap + Hébert 3.432–3.435 sign flip).

  **2.5d PLAN-OF-RECORD (2026-07-04, post-P0; rulings R12–R14 collected).**
  P0 corpus: explorer carrier map (`.claude/agent-memory/explorer/
  campaign_280_phase25d_carrier_map.md`, transient), cross-domain-attacker
  framing memo (`.claude/agent-memory/cross-domain-attacker/
  psi_half_seed_angular_trace_frames.md`), test-architect gate plan =
  **§16 of `a3_solve_transpose_verification.md`** (the 2.5d gate spec;
  its §16.A A1 cylinder-HAS-block leg is superseded by R12 — cylinder
  carries NO block; the gates' `# FRAMING` notes bind to Frame A).
  Diagnostics re-confirmed @ `ba16c4c`: sphere cold residual 5.18e5 /
  seedΔ 4.57e-2; slab 8.11e-16 / cyl 6.88e-16 both Δ=0.0-bit.

  The design (Frame A per R12 + the R13 corner completion):
  - **Carrier**: `FullField.starting_direction: StartingDirectionFlux | None
    = None` (Optional — all ~69 existing ctor sites compile; None ⟺ no
    seed-carrying levels). Presence per level by the STRUCTURAL predicate
    "μ_start ∉ the level's μ-nodes" (⟺ τ_raw ≠ 0): sphere-GL yes (1
    level), cylinder-product NO (its starting direction IS ψ₀ — the #229
    τ-clamp fact), Cartesian no. The leaf mirrors the
    `AngularBoundaryFlux` precedent: ONE flat backing buffer + shaped
    views — per seed level, BOTH directions (inward μ=−1 state + the
    OUTWARD starting direction, forced as state by R13: the corner
    outflow row must be linear in state), `cells(level, sign)` →
    ``(ng, nx)`` and `corner(level, sign)` → ``(ng,)`` (the (R,∓1)
    corner slots — inflow corner = given data / identity row; outflow
    corner = defect row). Role siblings: `StartingDirectionFlux` +
    `StartingDirectionSourceSink` on ONE `StartingDirectionSpace`
    (metric ZERO everywhere — the ghost treatment; the angular
    through-flux (1−μ²) vanishes at μ=±1, the SAME fact as α_{1/2}=0;
    masked-pseudo-inverse rides; the honest-scope note per §16.A A4:
    metric-invisible YET ACTIVE — constrained by B1 + C(i) + 2.5b's
    Euclidean Mᵀ, NEVER credited to G3 reciprocity).
  - **Dynamics**: the seed sub-system = per-seed-level ±1 SLAB-form legs
    on the 2.5a frame (Hébert 3.434 is the plain slab DD — the μ=±1 rays
    are straight characteristics; |μ_s|=1 on the sphere; the |μ_s|<1
    generalization documented as a seam for any future seed-carrying
    cylinder level). DAG order: seed⁻ ≺ ordinate legs (the ψ̃ chain
    consumes ψ½⁻ cells) and seed⁻ ≺ seed⁺ (pole continuation ψ½⁺(0) =
    ψ½⁻(0) — the SAME r=0 edge ordinates use); solve marches seed⁻ FIRST
    via `carlson_inward_sweep_from_source` on the TRUE q½ (the lag +
    `initial_guess`'s seed role die; kwarg survives per #285); apply
    evaluates seed rows on the GIVEN block (the extrapolation closure
    LEAVES the operator — the back edge vanishes); transpose reverses
    both legs + `angular_adjoint` STOPS at the seed cotangent (the
    strategy `seed_adjoint` delegation dies).
  - **Sources (R14)**: ONE fold helper `Q̄(μ=±1) = Σ_ℓ (2ℓ+1)/2·Q_ℓ·
    P_ℓ(±1)` (full (−1)^ℓ, B2b pin live); q_ext factories + S composite
    arm (+its S0.2 transpose arm) + F arm populate q½ on seed-carrying
    meshes; C = M[σ] gains the σ·ψ½-cells arm (the additive (L+C)
    decomposition pins force it — corners are trace-like, C skips them);
    B gains the corner arm (reflective: corner-out → corner-in via the
    seed pair; vacuum 0; albedo α·; white loud-deferred if unclear).
  - **Retirement** (aggressive-retirement + test migration per §16.D):
    the seed-strategy zoo (`PsiHalfAngleSeed` Protocol + ABC/registry +
    `ZeroSeed` + `CarlsonInwardSweep.__call__` + `AngularEdgeExtrapolation`
    + `seed_adjoint`s + `CarlsonSweepContext` + the closure's
    `psi_half_seed` field) — `carlson_inward_sweep_from_source` SURVIVES
    as the engine; 3-search audit (graph + code/tests/docs grep + direct
    ctors) before deletion.
  - **Commits**: d1 carrier + §16.A gates → d2 fold helper + source/C/B
    arms + §16.B gates → d3 the walk triple + retirement + §16.C/E/F
    gates + baselines (sphere re-captures with the 3-criteria record;
    slab AND cylinder assert-unmoved-first — a cylinder move HALTS; the
    characterization flips RED → the augmented triangularity certificate
    lands, LAPACK≡sweep leg stays 2.5b) → d4 estate rewires
    (`test_seed_threading_spy` → seed-independence;
    `test_inverse_operator_equivalence` premise rewrite; removal-form
    sphere JOINS; decomposition σ-arm update; dd_regression sphere
    re-capture; MMS/L1 absorb-gates confirmed) → d5 STATUS + memory.
    Acceptance = §16.G; keystone: sphere cold residual **5.18e5 →
    <1e-11**, seed-insensitivity **4.57e-2 → bitwise**, coarse-LS
    end-to-end **NaN → finite+positive** (both drivers).
- **2.5e — layout ruling (the R9 estate) + docs + close-out:** the sweep-layer
  layout decision from the carve's end-state; theory pages
  (`loss_representations.rst` the two-frames story + the orientation×kernel×
  execution table; `discrete_ordinates.rst` dev-history row); the doc-fix
  riders; close #280; comment #276 (A4 unblocked: the daggered posing consumes
  `A.H.inverse()`), #200 (curvilinear posture), #282 per ruling; the frozen
  baselines' disposition (keep as permanent canaries or retire post-carve —
  decide explicitly per the fuller-view rule).

### Deferral ledger (every entry must SURVIVE the unification, typed/loud)

multi-D Cartesian adjoint faces (existing raises); LD-slab adjoint faces
(kernel-pair registration — LD registers forward-only); `ScheduledInvertible
Operator` transpose; 2-D LD (pre-existing); Green/MatrixInverse/InverseOperator
`apply_transpose` (per-sibling, until consumers).

**RULINGS COLLECTED at open (2026-07-04): R10 = FULL route (a) in-phase; R11 =
(b) with (a)'s spelling — see the Rulings section. EXECUTION ORDER under R10:
S0 → 2.5a → 2.5d → 2.5b → 2.5c → 2.5e** (the #282 fix lands on the unified
apply-loop substrate — one seed-treatment change in the shared frame instead of
three bodies — and the reverse-scan is built once, against the fixed
formulation). Consequences threaded: the S0.5 frozen baselines pin the LAGGED
formulation for 2.5a's bit-identity, then the CURVILINEAR rows re-capture at
2.5d (principled re-baseline; Cartesian rows stay bitwise); 2.5b's sphere leg
gains the probed/assembled triangular-`Mᵀ` oracle the fix unblocks (the §12
successor gate built from matvec probes of the augmented system — extending the
2b assembly EMITTER to curvilinear stays an optional follow-up, not a 2.5
deliverable).

⏸ **C2.4 — COMPACTION TAKEN 2026-07-04 at the S0→2.5a seam.** Re-anchor
from: this section (shape + surface reconciliation + P0 findings + R10/R11 +
the S0 STATUS) → the gate spec `a3_solve_transpose_verification.md` §§9–15
(2.5a plan = §10, incl. the spy/AST-tripwire specs; scaffold order §13) →
the fragmentation map `a3_sweep_fragmentation_map.md` (structural claims
durable; line anchors re-derive by grep — the explorer's refreshed anchor
table lives in `.claude/agent-memory/explorer/campaign_280_phase25_p0a_map.md`,
transient until 2.5 merges). Branch `refactor/sn-walk-unification` @
`a3b8b38`, 7 ahead of main `28dbaee` (pushes held); harness tasks #22–#26
carry the sub-phase specs with R10 sequencing. **NEXT ACTION = 2.5a**: fold
`_apply_walk` (LR ~:2426) + `loss_action_transpose` (LR ~:2777, now behind
the S0.1 scheme guard) into ONE orientation-parametrized per-cell loop frame
over the one DAG — bit-identical BOTH orientations (the frozen
`walk_matvec_*` snapshots are the anchor; removal-form `array_equal` stays
the override-not-leak guard); mint
`tests/sn/sweep/core/test_one_dim_loop_walk.py` (wrap-spy proving both
orientations execute the ONE frame + the AST tripwire banning
`is_adjoint`/`is_forward`/`is_transpose`/`is_reverse` — orientation is an
OBJECT, the `_ApplyOperands`/`_SolveOperands`/`_SweepEmit` discipline's
sibling); the orientation object carries {cell order fwd/reversed, boundary
in↔out swap, Carlson mirror routing, `angular_adjoint` second factor};
resolve the S0-deferred multi-D predicate-honesty question in the frame
design. Then 2.5d → 2.5b → 2.5c → 2.5e per R10.

⏸ **C2.4b — COMPACTION TAKEN 2026-07-04 at the 2.5d plan→execution seam**
(user-called; the P0 corpus + rulings are banked, no 2.5d code yet).
Re-anchor from: the **2.5d PLAN-OF-RECORD block above** (design + commit
plan d1–d5 + acceptance numbers) + **rulings R12/R13/R14** (Rulings
section) → the gate spec **§16** of `a3_solve_transpose_verification.md`
(the full 2.5d gate plan; §16.A A1's cylinder-HAS leg superseded by R12)
→ the three P0 memos in agent memory (explorer
`campaign_280_phase25d_carrier_map.md` — the file:line consumer map,
TRANSIENT; cross-domain-attacker `psi_half_seed_angular_trace_frames.md`;
test-architect lessons). Branch `refactor/sn-walk-unification` @ the
checkpoint commit after `ba16c4c` (2.5a complete + this plan; pushes
held). Harness task #25 carries the ruled scope. Diagnostics
re-confirmed pre-fix: sphere 5.18e5 / 4.57e-2; slab+cyl machine-exact.
**d1 LANDED (the commit carrying this edit; 2026-07-04):** everything
the d1 NEXT ACTION specified, with TWO in-execution deltas ruled/forced:
(1) **R12a** — the presence predicate re-spelled to "first-ordinate raw
τ_raw ∈ (0,1) exclusive" (user-ruled; see the Rulings entry — the R12
letter mis-decides LS cylinders; `morel_montry_tau_raw_per_level` split
out of the clamped producer, bit-identical, as the single source);
(2) **`StartingDirectionDisplacement` minted alongside Flux/SourceSink**
— FORCED by the composite torsor algebra (flux−flux mints a displacement
PER LEAF through the Rep-keyed registry; without it every seeded
composite subtraction crashes). Deliverables: `StartingDirectionSpace`
(numerics/spaces/, flat layout level→sign(−1≺+1)→cells⊕corner, ghost
metric, typed views) + `StartingDirectionField` locus base (_bases.py)
+ 3 role leaves + FullField third member (REQUIRED-kwarg `_recombine`,
mixed-presence raises, mesh-keyed zeros presence, flat tail) +
TimedFullField (zeros/advance presence law/_recombine) + FullFieldSpace
third-block field-driven dispatch (illegal pairing raises) +
`SNMesh.starting_direction_levels`/`starting_direction_space` feeding
`full_field_space` + the §16.A gates (26, incl. the A2 drop-slice teeth
+ the cyl-LS discriminator row; A4 positive Mode-12 pin deferred to d3
as spec'd). Production call sites untouched — d1 is
bitwise-by-construction (fields stay 2-block until d3; the space may
carry the third block, dispatch is FIELD-driven — the one forbidden
quadrant, seeded-field-on-seedless-space, raises). Walls: pyright CLI
= 1 (#288); sn+transport+numerics serial green + DriftWarning-strict.
**d2 LANDED (the commit carrying this edit):** the R14 fold helper
`fold_moments_to_starting_direction(moments, sign)` — HOME =
`numerics/spaces/starting_direction_space.py` (NOT sn/spatial: the
transport-layer emitters consume it; sn→transport layering forbids the
reverse — the one in-execution relocation) — + the four operator seed
arms, ALL field-presence-gated (unseeded input ⇒ unseeded output —
verified dormant; d2 is bitwise): **C** = M[σ] apply σ·cells / solve
cells∕σ (post spectrum-gate) / transpose = apply (self-adjoint free),
corners zero both ways (trace-like, R13); **F** iso fold ½Q₀ both legs
through the helper; **S** forward = K_iso(φ₀) folded (ℓ=0 production
reach — honest-scope note; helper accepts ℓ≥1) + the S0.2 transpose
arm (seed-cells cotangent → ψ̄ₙ += wₙ·K_isoᵀ(½Σχ̄); output seed =
present-zero; measured Euclidean reciprocity 0.0 EXACT); **B** the R13
corner arm keyed on the law KIND tag (the realized operator is a
quadrature-row object and cannot act on the off-quadrature ray):
reflective = specular corner swap (μ=+1 ↦ μ=−1, its own mirror),
vacuum = zero, white/albedo/periodic = loud-deferred
NotImplementedError; transpose = the mirror image (reciprocity 0.0).
**Deviation (forced, banked):** the q_ext-FACTORY population moved
d2 → d3 — a 3-block q_ext against 2-block iterates raises the
mixed-presence law, so factories flip WITH the birth sites in d3's
atomic activation. §16.B gates in
`test_psi_half_angle_seed.py` (+9): B1b uniform+graded O(Δr²) to the
exponential closed form (measured ratios 3.75→3.93; smallness 2e-4 vs
the true 7.4e-5@nx64, wrong-limit plateau O(1e-1)) + the 4-mutation
matrix (index-drift: uniform-blind-by-construction GREEN, graded RED —
the Mode-5 keystone); B2a ½Q₀ exact; B2b 2-term ½Q₀−(3/2)Q₁ live +
sign-drop teeth. Walls: pyright CLI = 1; sn+transport+numerics serial
green.
**d3 IMPLEMENTED (uncommitted, 2026-07-04) — the walk triple + retirement
+ §16.C/E/F + baselines + the 3-wave test migration + the elegance fixes.**
KEYSTONE VERIFIED: **sphere cold residual 5.18e5 → 2.5e-16, seedΔ
4.57e-2 → 0.0 bitwise** (`test_282_direct_seed_fixed_point.py` 12/12 —
§16.C C(i)/C(ii)/C(iii)/C(iv) + §16.F Mode-11 wrap sentinel / Mode-12 A4
positive pin / Mode-10 activation). Delivered:
- **The walk triple.** `_run` marches seed⁻ directly per carrying level
  via `carlson_inward_sweep_from_source` on the TRUE q½ (the LR lag
  block is gone); `_apply_walk`/`loss_action` emit the seed rows
  (`_seed_rows_forward`); `loss_action_transpose` reverses them
  (`_seed_rows_transpose`) + `closure.angular_adjoint` STOPS at the
  seed cotangent (returns `(psi_bar, seed_cells_bar: dict)`); non-carrying
  cylinder levels inline `edge_extrapolated_seed` (the banked product-cyl
  ig-consumption preserved bit-exactly — cyl frozen baseline UNMOVED).
- **Strategy-zoo RETIRED.** `psi_half_angle_seed.py` 851→161 lines (ONE
  `carlson_inward_sweep_from_source` engine, now returns `(cells, face)`);
  `PsiHalfAngleSeed`/`ZeroSeed`/`CarlsonInwardSweep`/`AngularEdgeExtrapolation`/
  `CarlsonSweepContext`/`seed_adjoint` all gone; 3-search audit clean.
- **Atomic birth-site flip** (solver SI/Krylov cold-starts + q_ext
  factories + reconstruction sweep + `transport_sweep` self-sufficiency;
  `streaming.solve` allocates the ψ½ carrier). The FOURTH role member
  `StartingDirectionResidual` minted (forced by `evaluate_residual`).
- **§16.D re-baselines:** sphere `walk_matvec_sphere_2g.npz` +
  `affine_carve_baseline/*_SPH.npy` re-captured (git-verified slab+cyl
  BYTE-UNMOVED); the §16.E characterization flipped RED→the augmented
  triangularity certificate (`test_282_augmented_walk_order_is_triangular`,
  probed `triu==0` EXACT + teeth).
- **The R14 FULL (−1)^ℓ fold LANDED** (`StartingDirectionSourceSink.
  from_angular_source` now folds ALL Legendre moments — REQUIRED: an
  isotropic trial streams to a μ-linear source `q=μA'+σ_t A`, so the
  ℓ=1 term carries `−A'` at μ=−1; ℓ=0-only floored the anisotropic
  curvilinear MMS. Isotropic sources collapse to ½q₀ bit-exactly ⟹
  eigenvalue path unchanged). The anisotropic pole-cell MMS now
  converges O(h²)-to-exact — route (a) VERIFIED correct.
- **Elegance review (clean, no blockers) fixes applied:** single-source
  `_require_starting_direction` guard, `_seed_residual_march` collapses
  the forward legs ("orientation is data"), 8 `# type: ignore` dropped.
- **3-wave test migration** (test-architect ×2 + main agent): the whole
  sphere-consumer estate rewired to the 3-block composite (uniform
  `starting_direction=StartingDirectionFlux`; value tests use the
  consistent `starting_direction_edge_seed`; Euclidean reciprocity uses
  random-seed + full-space dot; G-reciprocity uses present-zero seed).
- **pyright CLI = 1** (the accepted #288 residual — the carve added zero
  new type errors).
- **The sphere EIGENVALUE re-pose — RESOLVED (numerics-investigator, ruled
  PRINCIPLED).** The decisive discriminator for a seed/angular-CLOSURE
  change is NOT h→0 (which differs — a seed IS a closure, so it changes the
  O(N) truncation) but an **angular-order N-sweep at fixed mesh**: the OLD
  edge-extrapolation and the NEW direct Carlson seed both converge to the
  SAME transport eigenvalue as N→∞ (agree ~1e-6 by GL32; the ~1.7e-3 gap at
  GL8 is pure low-N seed truncation). The MMS is BLIND to this (every
  curvilinear ansatz is ≤ linear-in-μ = the seed's EXACT regime, Mode 7),
  so MMS-O(h²) does NOT certify the seed — the N-sweep does (new gate
  `test_heterogeneous_1g_angular_order_consistency`). **HONESTY:** route (a)
  is justified STRUCTURALLY (the honest single-pass direct inverse — cold
  residual 5.18e5→2.5e-16, for #200/#280), NOT by angular accuracy — at the
  tests' GL8 the OLD seed is actually CLOSER to the N-limit; NEVER frame the
  re-baseline as "more accurate". Landed: `dd_regression[sphere_2g_3reg_dd_n40]`
  §16.D re-captured (sphere-only, git-verified); the fragile n∈[5,10,20]
  `diff_2<diff_1` ladder robustified to n∈[10,20,40] (the n=5→10
  near-coincidence Δ≈8e-7 is a REAL coarse-mesh feature; the pole cell is
  O(h^~1.4) ⟹ global rate sub-quadratic but convergent).
- **The Krylov INNER stall — FIXED (NEW route-(a) regression, ERR-053
  family, NOT #200).** Root cause: `n_dof` (→ scipy-gmres `restart`) was
  sized from the BULK (`N·ng·∏spatial`) at both solver Krylov drivers, but
  route (a) grew the composite to 3 blocks (bulk⊕trace⊕ψ½-seed), so
  `restart < to_flat` on the sphere ⟹ restarted GMRES truncated the
  augmented subspace and STALLED (info>0, 868 s, WRONG keff under an outer
  cap). Fixed: both sites size `n_dof = initial_guess.to_flat().size`.
  Regression gate promoted: `test_krylov_restart_covers_augmented_composite`
  (the trace+seed deficit pin). SI≡Krylov sphere equivalence now green.
- **pyright CLI = 1** (accepted #288; zero new). **d3 COMMITTED @ `a29ab2d`
  (2026-07-04) — full sn+transport+numerics wall GREEN (3149 passed).** The
  3-wave sphere-consumer test migration (test-architect ×2 + main agent),
  the elegance fixes, the full-fold, and the Krylov+eigenvalue resolutions
  all in that commit. d4 estate rewires were ABSORBED into the migration
  (test_invertible_operator / decomposition / dd_regression re-capture all
  landed green); d5 memory banked ([[feedback-nsweep-discriminator-closure-repose]]).
  **NEXT: theory-doc §16/#282 close-out (Cardinal Rule 3 — dispatched to
  archivist) → then 2.5b (reverse-scan sweep_transpose) per R10.** Pushes
  HELD (branch `refactor/sn-walk-unification`).

⏸ **C2.4c — COMPACTION TAKEN 2026-07-04 at the 2.5d→2.5b seam (d3 landed).**
**2.5d is COMPLETE** — d1 `9fc066d` (carrier) + d2 `5170f20` (fold helper +
arms) + **d3 `a29ab2d`** (the walk triple + retirement + §16 gates +
baselines + the 3-wave migration + the full-fold + the Krylov/eigenvalue
resolutions). Branch `refactor/sn-walk-unification` @ `a29ab2d` (through
d3), **pushes HELD**; the tree is clean EXCEPT (a) the sub-agents'
`.claude/agent-memory/*` updates (NOT main-agent's to commit) and (b)
`scratch/` (the user's). **Re-anchor after /compact from:** this C2.4b/c
block (the full d3 record above — keystones, the R12a predicate, the
full-fold gotcha, the N-sweep discriminator, the ERR-053 Krylov fix) +
`git log` (trust git, not the summary — d3 = `a29ab2d`) + the gate spec
`a3_solve_transpose_verification.md` §16 (verified) and, for 2.5b, §§11–12
(the reverse-SCAN G1/G2 + the spherical LAPACK-≡-sweep leg the augmented
triangularity certificate now unblocks). Durable memory: the N-sweep
discriminator is banked at
`~/.claude/.../memory/feedback_nsweep_discriminator_closure_repose.md`.
**d3 fully landed** — the theory-doc §282 close-out committed @ `5b03a37`
(Cardinal Rule 3: discrete_ordinates.rst +753 lines — back edge, augmented
triangular normal form + certificate, R12a trichotomy, R14 full fold,
retired-zoo rationale, the N-sweep methodology + the "structural NOT
accuracy" `.. warning::`; Sphinx `-W` clean), the investigation diagnostics
archived @ `e71c225`. Branch @ `5b03a37`, **19 ahead of main**, pushes HELD;
tree clean except agent-memory (sub-agents' own) + scratch/. **ONE follow-on
for the fresh session:** FILE the FullField construction-invariant issue the
elegance reviewer flagged (make "carrying mesh ⟹ ψ½ block present" a
`__post_init__` invariant so the `_require_starting_direction` guards can
retire) — a `module:sn` + `module:transport` `type:improvement` GitHub
issue. **NEXT ACTION = 2.5b**
(the reverse-scan `sweep_transpose`: G1 round-trip, G2 dense `(L+C)ᵀ⁻¹`
oracle, the spherical LAPACK≡sweep successor leg on the augmented system,
`SweepOperator.apply_transpose` per 2.5c/R11) → 2.5c → 2.5e (R9 layout +
docs + close #280).

### 2.5b EXECUTION STATUS (2026-07-05) — reverse-scan slab+sphere DONE; cyl reframed

**The reverse-scan `(L+C)⁻ᵀ` for SLAB + SPHERE is implemented + verified to
machine precision** (uncommitted → commit imminent). Deliverables:
- `ordinate_scan_transpose` (`sn/spatial/scan.py`) — the affine scan's adjoint,
  itself an `ordinate_scan` on reversed-shifted coefficients (inherits the
  Blelloch/pair-monoid denormal robustness; single source of truth). Unit-verified
  vs a dense transpose (2e-16, scalar + vector lanes).
- `carlson_inward_sweep_transpose` (`sn/spatial/psi_half_angle_seed.py`) — the
  Hébert §3.9.4 seed-march adjoint. Unit-verified 2e-16.
- `_OneDimScanWalk._run_transpose` + `.sweep_transpose` + `CumprodScan.sweep_transpose`
  (`sn/loss_representation/__init__.py`) — the reverse-mode adjoint of `_run`. SLAB
  arm (scan transpose + boundary) + SPHERE arm (M-M thread reverse + coupled-pole
  mirror reverse + seed-march transpose). Curvilinear-non-sphere (cyl) `raise
  NotImplementedError` — the typed deferral until the cyl-fwd fix (below).
- `InvertibleOperator.solve_transpose(b)` (`sn/operators/streaming.py`) — the
  FullField packer (mirror of `.solve`; 2.5c wires `SweepOperator.apply_transpose`
  + the swap law to it).
- Gate `tests/sn/operators/test_loss_transpose_solve.py` — G1 bulk round-trip +
  G2 augmented dense-`Mᵀ` (source-carried slots) + assembled-`Mᵀ` LAPACK (slab),
  slab + sphere. (Full §11/§14 mutation suite + cyl = the test-architect pass.)

**KEY FINDING — the #284 outflow subtlety (extends to the transpose).** `solve`
COMPUTES the outflow slots (slab boundary-outflow / sphere seed `corner_out`)
while `apply` treats them as FREE DOFs, so `solve = (L+C)⁻¹` only on the SOURCE
SUBSPACE (`M_solve·M_apply = I` on bulk rows to 5e-16; the outflow slots differ
O(1)). `sweep_transpose` is the FAITHFUL transpose of the solve — verified
`matrix(solve_transpose) = M_solveᵀ` to 2e-16 (slab AND sphere augmented) — so it
matches the structurally-independent apply-oracle `M⁻ᵀ` on EVERY source-carried
slot and deviates ONLY on the outflow slot (where `M_solve`'s outflow column is
provably 0). ⟹ every reverse-solve gate compares source-carried slots (bulk ⊕
seed cells ⊕ inflow corner), NOT the full `to_flat` (the base gate spec §3's
full-`to_flat` oracle is WRONG on the outflow — superseded by the memo's augmented
`_probe_augmented_matrix_one_group` oracle). Document in 2.5e.

**CYLINDER — user ruled "fix it FUNDAMENTALLY first" (2026-07-05).** The cyl is
NON-carrying (R12a), so route (a) never gave it a direct seed: its forward solve
seeds the M-M thread from `edge_extrapolated_seed(ITERATE)` — a lag (cold cyl
solve ≠ `(L+C)⁻¹`, 0.57 bulk err; threaded → 6.7e-16). So a cyl transpose-solve
cannot be single-pass until the FORWARD is fixed. The user chose the fundamental
fix ("route (a) for the non-carrying cylinder") over defer/thread.

**Feasibility CONFIRMED (numerics-investigator + explorer, 2026-07-05):**
- **FEASIBLE — a PURE-DIAGONAL fold.** For a PRODUCT quad, R12a `τ_raw=0` ⟹ the
  seed is `t=0` = the first-swept ordinate ψ_{m0}'s own flux (#229 "starting
  direction ≡ ψ₀"). Probed: m0's output couples to NO other ordinate, its
  self-block is triangular, and the seed's contribution is EXACTLY diagonal —
  coefficient `κ_c = dA_w[m0,c]·c_in[m0]` (the #229 τ-clamp `max(0.5,·)` DEFINES
  the live `(1−τ)/τ=1`). POC: fold κ into m0's cell diagonal → single-pass =
  `M_aug⁻¹` (5e-16). LEVEL-SYMMETRIC is the dead case (`c_in[m0]=0` — the fold is a
  no-op). Degenerate pure-azimuthal ordinates sit downstream of m0 — no interaction.
- **SOLVE-ONLY + BASELINE-NEUTRAL.** The matvec is ALREADY lag-free/triangular
  (needs NO change). Unlike the sphere route (a) (which moved keff(N)), the cyl
  converged seed IS ψ_{m0} either way, so the FIXED POINT is machine-identical:
  keff/MMS/matvec `walk_matvec_cyl_2g`/converged snapshots do NOT move. Only a cold
  single-sweep value moves (no fixture captures it) — a principled cold-correctness
  + preconditioner improvement (positive #200 cyl-half implication).
- **CORRECTS A FALSE THEORY CLAIM (Cardinal Rule 1/3).** `docs/theory/
  discrete_ordinates.rst:10809-13` states "the cylinder was already exact… α-dome
  telescopes… spherical-only" — FALSE for the product cyl (it genuinely lags). The
  claim propagated from a diagnostic (`diag_curvilinear_seed_sensitivity.py:72`)
  that measured ONLY the LS (dead-seed) rule. This mis-attribution (also in
  `test_282` docstring:629 + ERR-026 docs + the old L14 taxonomy) is WHY the carve
  was never scheduled. The archivist must correct it when the fold lands.
- Minimal change site: `_run` non-carrying `else` branch (the iterate reads at
  `sn/loss_representation/__init__.py` ~lines 4036/4039/4096-4100); fold κ into m0's
  WDD denom (like `inverse_denom` carries `dA_w·c_out`). Investigation diagnostics
  banked: `derivations/diagnostics/diag_280_cyl_{product_seed_lag,fold_feasibility}.py`
  (the feasibility proofs are permanent triangularity/routing gates → promote to
  `tests/sn/sweep/`; the strict-xfail `test_cold_single_pass_equals_matvec_inverse`
  is the acceptance certificate — flips GREEN when the fold lands). Explorer map:
  `.claude/agent-memory/explorer/campaign_280_phase25b_product_cyl_seed_map.md`.

**NEW SEQUENCE (user-ruled):** ✓ primitive+slab+sphere reverse-scan (done) →
✓ **2.5b-cyl-fwd DONE @ `ba202a1`** → ✓ **cyl reverse-scan DONE @ `f1ddeb6`** →
2.5c wiring (NEXT) → 2.5e docs + close #280.  (The #29 "full gate suite" folded
INTO `f1ddeb6`: G1–G5 + the (a)/(b)/bdry mutation matrix were built + teeth-
verified inline, so #29 is discharged — no separate test-architect pass owed.)

**2.5b-cyl-fwd LANDED @ `ba202a1` (2026-07-05):** the product-cyl direct-seed
fold. For a product quad ψ½ ≡ ψ̄_{m0} (t=0, #229 — first-swept ordinate's own
average) is a per-ordinate SELF-reference; its coupling κ=(ΔA/w)·c_in folds
into m0's cell diagonal (`c_out → c_out − c_in`, precomputed as `seed_fold` via
the single-source `affine_scan_coefficients`) with NO angular-upstream source,
and m0 reads its own average as the M-M upstream. Cold cyl solve 0.57 → 4.4e-16
(single-pass = M⁻¹). Iterate read RETIRED (both geometries direct-inverse now).
LS is a no-op (c_in=0); matvec UNCHANGED; a t≠0-and-live-c_in level raises loud.
Baseline-neutral (cold ≡ SI-converged, machine-identical). Gate
`tests/sn/sweep/test_cyl_direct_seed_fold.py` (6 gates, promoted from the retired
`diag_280_*` diagnostics; mutation-verified RED under -O on c_out+c_in). Verified:
sweep+operators -m "not slow" 1305/0; pyright 0 on the file; slow L1 cyl
cross-check `test_unified_cylinder_l1_mr_2g_trajectory_resolvent` 618s PASS (a
`@slow` test right at the 600s cap — not a regression). Theory-claim correction
LANDED (follow-on docs commit): the false "cylinder α-dome telescoping absorbs
the wrong seed / was already exact" (a level-symmetric-only mis-attribution)
corrected across 9 theory sites + 4 test docstrings + a Development-history
changelog row; the VALID scalar-blindness claim (anti-pattern #8) preserved;
Sphinx -W -E clean. The correct framing: seed-insensitivity is the LS dead
first-ordinate weight (c_in[m0]=0), NOT telescoping, and false for product.

DEFERRED to 2.5c (coherent carve, not fragmented): the now-vestigial iterate
threading — `_initial_guess_values` (0 callers) + the `_run`/`sweep`/solve
`initial_guess` param — retires WITH the SweepOperator direct-inverse predicates
(`is_invertible`, `apply_transpose`). The curvilinear solve no longer reads it.

**2.5b cyl reverse-scan LANDED @ `f1ddeb6` (2026-07-05):** `(L+C)⁻ᵀ` for the
cylinder — `_OneDimScanWalk._run_transpose` extended as the transpose of `_run`'s
unified curvilinear body (the sphere-only block generalised to sphere+cyl,
Pattern 2; sphere path bit-identical since `seed_fold` is empty ⇒ `is_seed_ord`
always False). Three cyl-specific structures: (1) the multi-level M-M thread
transpose (sphere is single-level `[None]`; cyl multi-level, each independent,
`psi_angle_bar` re-init per level); (2) the m0 direct-seed-fold transpose — the
seed ord routes `−mm_a_in` to `psi_avg` (ψ½≡ψ̄_{m0}, its OWN average), folded
coeffs `(c_out−c_in)`, no ang_contrib, `m_seed=None` (non-carrying, no carrier);
(3) the pure-azimuthal DEGENERATE ords (product quad, 8 for `product(4,8)`) as
slot-local DIAGONAL transposes — the caches are VALID at A_down=0 (`inverse_denom
= 1/(dA_w·c_out+Σ_t·V)`, probe-confirmed 0-ULP), so NO scan / NO recompute (a=−1
would wrongly couple cells — the reason the forward uses `dag_walk` not the scan).
Retired the cyl `NotImplementedError` guard.

**The G3 catch (a real bug G1/G2-bulk missed):** the degenerate ord does NOT
overwrite its `bc_outer` slot (no face march), so the forward passes
`sol.boundary[deg] = q.boundary[deg]` (identity); μ<0 degenerate ords ALSO read
that slot as the `|μ|`-weighted spatial upstream. The first cut DROPPED both —
G1 round-trip + G2 dense-Mᵀ (bulk-only for the non-carrying cyl) stayed GREEN, but
G3 full-field Euclidean solve-reciprocity (bulk⊕boundary) reddened at 2.77e-2.
Adding the passthrough + the `|μ|·A_total` spatial cotangent made it **EXACT
(0.0)** — the #284 boundary-subspace lesson, live.

Gates (`tests/sn/operators/test_loss_transpose_solve.py`, +`cyl_product`/`cyl_ls`
in `_MESHES`): G1 round-trip, G2 dense-Mᵀ oracle (structurally independent from
the forward apply), G3 full-field reciprocity (all EXACT on cyl), G4 `m_seed=None`
contract, G5 the mandatory-config activation sentinel. **Teeth mutation-verified**
(committed baseline → mutate → RED → `git checkout` revert): (a) seed misroute →
G2+G3 RED; (b) degenerate DROP → G3 RED / `cyl_ls` control GREEN; (bdry)
passthrough DROP → G3 RED / G2 bulk GREEN (blind) / `cyl_ls` GREEN — the
product-RED/ls-GREEN asymmetry IS the evidence the product config is mandatory
(the ERR-066 blindness pin). Verified: reverse-scan gate 15/15; SN
operators+sweep -m "not slow" 1315/0; pyright ratchet sn 0 (total 1 = #288).

DEFERRED to 2.5c (unchanged): the vestigial iterate threading
(`_initial_guess_values` + the `initial_guess` param) retires WITH the
SweepOperator direct-inverse predicates.

**2.5c LANDED — inverse-adjoint wiring + the `initial_guess` retirement
(2026-07-05):**

*Part A+B @ `45529c2`* — the swap law `A.H.inverse() ≡ A.inverse().H` as an
OBJECT IDENTITY of the algebra: `SweepOperator.apply_transpose(b) =
inner.solve_transpose(b)` (the 2.5b reverse-scan `(A⁻¹)ᵀ = (Aᵀ)⁻¹`) +
`is_adjointable = isinstance(inner, InvertibleOperator) and inner.is_adjointable`
(schedule-folded sibling stays deferred); `_AdjointOperator.inverse() =
inner.inverse().H` + `is_invertible = invertible(inner) and
adjointable(inner.inverse())` (R11). The metric adjoint-solve
`A.inverse().H.apply(b) = G⁺·solveᵀ(G·b)` falls out of the EXISTING
`_AdjointOperator.apply` FOR FREE (the a3 "Deliverable 3" dissolved). Gates
`test_inverse_adjoint_coherence.py` (19, test-architect): G1 forward-matvec
G-reciprocity / G2 swap-law value bit-identical / G3 round-trip / Mode-11 wrap
sentinel / M-ADJ-swap+metric mutations RED-verified / predicate flips /
assert_type. **Spec refinement:** M-ADJ-metric reds on ALL geometries (the
non-uniform slab metric `G = V·wₙ ⊕ |Ω·n|·wₙ` is non-trivial) — stronger than
§13's sphere/cyl-only prediction; the `.H`≠Euclidean claim holds universally.

*Part C @ `8cf5215`* — the `initial_guess` seed retirement, scoped by the
USER's warm-start insight (a NN-learned x₀ is the 3-vs-30-iterations lever; it
lives at the ITERATION layer, NEVER on a direct exact sweep). Post-2.5d the
curvilinear ψ½ is the DIRECT #282 seed, so the SN sweep's seed-USE is dead:
DELETE `_initial_guess_values` (0-caller orphan) + drop the dead param/threading
from `_solve_timed_full_field` (Invertible+Scheduled) + the loss_representation
`sweep`/`_run`/`transport_sweep` kernels. The inverse family's `.solve`/`.apply`
KEEP `initial_guess` as the uniform `SupportsSeededApply` kwarg (#285) but
ACCEPT-AND-DROP (joining `WindowedSweep`/`MatrixInverseOperator`/
`_SeededExactApply`). **The surviving warm-start is untouched + separately
gated** — Green's `M-GRN-SEED` (the iterative inverse's genuine warm start) +
`SourceIteration.solve(initial_guess=x₀)`. Test migration: `test_seed_threading_
spy.py` RETIRED (dead premise — "dropped seed = wrong FP" is FALSE post-2.5d);
seed-INDEPENDENCE + accept-drop path pins replace the thread pins; the extractor
+ seed-mesh-validation pins retired; `test_rhs_boundary_seeds_the_sweep_inflow`
migrated (the rhs.boundary→boundary_buf contract stays). Part-A completion
folded in (broad-run catch): `_AdjointOperator.inverse()` guard matched to
`is_invertible` (NotInvertible not MissingAdjoint when the inner's inverse is
non-adjointable); AdjointWrapper keystone row STRUCTURAL_ABSENT→VALUE_RAISE. Net
−151 LoC. Full inverse-family blast radius 2271/0; pyright transport:1.

✅ **2.5e LANDED — campaign #280 walk-unification COMPLETE (2026-07-05).** Docs +
the retired-ref cleanup + close-out committed on branch `refactor/sn-walk-
unification` (pushes HELD). Delivered:

- **Theory docs (archivist).** `loss_representations.rst` (+366): the new section
  *"The orientation axis — the adjoint completes the 2×2"* — the two frames
  (apply-loop `{loss_action, loss_action_transpose}`, bit-identical both
  orientations via 2.5a; solve-scan `{sweep, sweep_transpose}`, the coherent
  reverse-scan via 2.5b), the **orientation×kernel×execution taxonomy** table
  (coherence axis = ORIENTATION; execution the NON-free third axis pinned by
  `(kernel,dim)`), the discrete-Euclidean-transpose `.. warning::` (not μ-reversal,
  not the continuous adjoint), the deferral ledger, the **swap law**
  `A.H.inverse() ≡ A.inverse().H` as object identity (+ the free metric
  adjoint-solve, no metric code in the sweep), the `initial_guess`/warm-start
  architecture split (direct inverses accept-and-drop; the warm start lives at the
  iteration layer). `discrete_ordinates.rst` (+53): dev-history rows for 2.5b
  (reverse-scan `sweep_transpose` — the empty 2×2 cell filled) + 2.5c (swap-law
  wiring onto the #226 taxonomy). Sphinx `-W` exit 0 / 0 warnings.
- **Retired-ref cleanup (main-agent — the 3-search text-grep leg the 2.5d audit
  missed).** 3 live-framed refs to the retired seed-zoo classes repointed at the
  route-(a) treatment: `solver.py:2287` + `augmented_mesh.py:285` (→
  `carlson_inward_sweep_from_source`), `reduced_operator.py:338` (mu_start →
  `MorelMontryAngularSweep` seed treatment). Full grep now clean — only historical
  tombstones (`psi_half_angle_seed.py` dev-history, `spatial/__init__.py`,
  `pole_angular_closure.py:1361`) remain. The archivist also corrected the
  `pole_angular_closure.py` module/class docstrings (they described the retired
  proxy-source seed as live).
- **R9 estate — CONFIRMED (the layout decision).** `sn/spatial/` (`pairing`,
  `pole_angular_closure`, `psi_half_angle_seed`, `scan`, `sweep_cache`) KEEPS its
  name this campaign per R9 — renaming churns ~150 docs refs + 17 test imports, and
  Phase 2.5 already touched every module twice. The name is STALE (sweep-walk +
  angular machinery, not "spatial"); the rename is DEFERRED to a follow-up issue.
  **EXECUTED (task #54, 2026-07-13): `sn/spatial/` → `sn/sweep/`** (the honest
  name — the package is the sweep's kernel substrate; the walk executors stay in
  `loss_representation`); `sweep_cache.py` → `cache.py` rode along, the
  `tests/sn/spatial/` mirror dissolved to true homes (sweep/core ×2,
  transport/spatial ×2, numerics ×1), and the pre-existing #272 debt
  (mutation-TOML module-paths, the broken `diag_276` LD import, the dangling
  `tests.sn.spatial.*` doc roles) was settled in the same pass. This bullet is
  now archaeology.
- **Frozen baselines — KEEP (fuller-view rule, explicit).**
  `test_walk_matvec_baselines.py` + `walk_matvec_{slab,sphere,cyl}_2g.npz` stay as
  permanent regression canaries: they pin the FORWARD-walk matvec output
  byte-for-byte across all 3 geometries — distinct coverage from the 2.5b G2
  transpose oracle and the augmented-triangularity certificate; the cheapest
  forward-walk-drift tripwire; 2.5d already re-captured the sphere snapshot
  (principled). NOT procedural redundancy — a genuine structural reference.

**Issue close-out (prepared — outward `gh` writes gated on the user's go-ahead,
pushes HELD):** the 2.5e commit carries `Closes #280` + `Closes #282` (fire at
push/merge; #282's fix landed @ `a29ab2d` without a trailer, so the campaign
close-out commit carries it — go-forward, no history rewrite). Comments prepared:
**#280** (swap law + walk unification + warm-start split), **#276** (A4 unblocked —
the daggered posing consumes `A.H.inverse()`; `sweep_transpose` is A4's first
consumer), **#200** (curvilinear posture), **#282** (route (a) fix landed, keystone
5.18e5→2.5e-16 — closing as fixed). Two follow-ups to FILE: (1) the FullField
construction-invariant (make "carrying mesh ⟹ ψ½ block present" a `__post_init__`
invariant so the `_require_starting_direction` guards retire) — `module:sn` +
`module:transport`, `type:improvement`; (2) the R9 `sn/spatial/` rename —
`module:sn`, `type:improvement`.

**CAMPAIGN #280 COMPLETE.** sn+numerics wall GREEN (2809 passed / 5 skipped / 36
xfailed @ 2.5e); pyright ratchet `transport:1` (accepted #288, zero new).

**NEXT = the codim-1 `FaceField` substrate carve — BEFORE P3 DSA** (user ruling
2026-07-06). Design note: **`.claude/plans/archive/facefield_codim1_design.md`** (a full
design-only session dug the ψ½ pole-seed physics + the `BulkField`/`FaceField`
duality). It: names the one `face_streaming_normal` measure (collapsing the twin
trace-metric constructions); typed-keys `FaceLayout`; introduces the `FaceField`
codim-1 parent + re-parents the trace/starting-direction fields; makes the
composite's block-list mesh-derived (retiring the 7 `_require`/`_refuse`
starting-direction guards, subsuming the banked `FullField.__post_init__`
follow-up). Lobatto fold-in of the pole was measured **affordable but declined**
(keeps the cell-centered bulk clean). Rationale for sequencing (clean-before-
extend): DSA's R/P consume the composite + trace under the shared metric, so a
clean codim-1 substrate + unified face measure benefits the DSA boundary
restriction. Deferred-with-triggers: `SaddlePointOperator` ← diffusion mixed-form
(#294); structured *spatial* `FaceField` ← #294 / CP-2D.

**THEN P3 DSA #2** — opens with the 3-P0 dispatches (literature-researcher checks
`scratch/literature/` FIRST; test-architect FP-invariance/⚠#215 Mode-9;
cross-domain-attacker R/P Petrov–Galerkin frame). Branch
`refactor/sn-walk-unification`, pushes HELD.

---

## Phase 3 — DSA for SN (#2; folds #200's Krylov posture; runs AFTER the #280 walk
unification per R1 — the accelerator wires into the unified orientation×kernel walk)

### ⏵ Phase-3 plan-of-record (USER-APPROVED 2026-07-26; branch `feature/sn-dsa`)

**3-P0 is DONE** — the four memos are the working spec set (all in `.claude/plans/`):
`dsa_literature_memo.md` (Alcouffe 1977 / Larsen 1982 I+II / Adams-Larsen 2002 /
Morel 1982, all local + scan-verified; §6 = the build synthesis; §7–§8 Adams-Martin
1992 + WWM-McGhee-Lehoucq 2004 extraction in flight — the user added both
2026-07-26), `dsa_verification_spec.md` (13 gates D1–D13 + mutation matrix),
`dsa_rp_frame_analysis.md` (frame verdicts), `dsa_landing_zone_recon.md`
(current-tree facts @ `e4c1a81c`).

**The four structural verdicts:**
1. **(R,P) = the ℓ=0 faces of the EXISTING `Quadrature.angular_frame(0)`**
   (`GalerkinFrame`; `numerics/quadrature/directional.py:417`) — mint NOTHING;
   0-ULP pins vs `integrate_angular`.
2. **Consistency's named algebra: `A_low = Schur_{ℓ=1}(R₁·A_assembled·P₁)`** — the
   two-moment Galerkin triple product over the ASSEMBLED DD operator, then
   Schur-eliminate the current. "Consistent" = derived-by-moment-reduction (Reed's
   scheme had fixed-point compatibility and still diverges — keep lexically separate).
3. **The correction→0 partition** (spec §0): FP-invariance gates are blind to 7/8
   accelerator errors — object + rate gates carry the verification weight.
4. **Marshak = the boundary Fick** (half-range Schur under the shared `|Ω·n|w`
   metric); vacuum→Marshak(𝒜=0) is exactly the error-problem BC.

**Diffusion readiness (the approval checkpoint's finding):** the "2-group, not
operator-algebra" picture is the PRE-#290 island — today's `orpheus.diffusion` is
ng-generic operator algebra (`loss = leakage + collision − scattering − boundary`,
`diffusion/solver.py:240`; `test_three_group_is_structurally_beyond_the_island`),
TransportMethod-conformant, assembling, O(h²) MMS-gated. Open gaps #292/#293/#294
do not bite arm 1 (#292 is legacy-table data — DSA reads the SN problem's own
mixtures; #294 is off-face interfaces — DSA shares the SN mesh; diffusion is 1-D
today — matches arm-1 slab scope). Whether the landed harmonic-mean stencil IS the
DD-consistent operator is 3a's computed question (R4 rules on the D2 matrix diff);
if ≢, the derived stencil wires through the same assembly machinery (double-duty:
standalone keeps physics-chosen RT0, DSA uses consistency-derived — two defining
laws, not a twin path).

**Scope (arm 1):** slab Cartesian DD, within-group fixed-source, P0-DSA, f-form;
BOTH postures (SI+DSA = new driver construct — `SourceIteration` has no hook,
injection after the resolvent apply `iteration.py:682`; Krylov = replace the
explicit identity preconditioner at `sn/solver.py:487` with sweep + P·A_low⁻¹·R).
OUT with seams/follow-ups: curvilinear (#282-blocked; no stability theory), LD arm
(R5 decide-from-table; Larsen I §IV–V + Adams-Martin are the shelved references),
2-D Cartesian (Alcouffe corner moments + windowed-sweep interaction — file at
close), k-outer, upscatter accel, P1-DSA (gated on Σs1/Σt per A&L (7.14)).
`full_field_space` lands on `TransportMethod` (the P7b trigger fires). The
foldable/σ_r accessors get their intended consumer through the LOW-ORDER BUILD
only — never the sweep (the #215 trap stays fenced; ERR-070 files when D3 first
reds on it).

**Key gate anchors (3c):** S2 K₂=0 one-iteration exactness (the sharpest BC unit
test); ρ vs A&L (3.65) (Σ_th-independent, the primary quantitative gate);
reflective stability = D12 (the historical spike failure); the partial-consistency
negative control (WD a=0.5 + DD low-order must diverge Table-II-shaped); ρ
estimator = residual RATIO. Docs close-out = a NEW page in
`docs/theory/methods/sn/` (the old `discrete_ordinates.rst` target is dissolved);
fix the content-drifted `diffusion_1d.rst:517` xref; flip the
`field_algebra.rst:528` promise. Derivation-side watch items: Larsen is the
transcription reference (Alcouffe (17)/(23) carry PRINT SIGN ERRATA — memo §1.5);
the sources' Σw=1 vs A&L's Σw=2 vs ORPHEUS's Σw=2 slab convention maps ONCE with
numeric asserts; the ⅓ enters as the quadrature-moment property W₂/W₀, never a
transcribed constant.

⏸ COMPACTION points at 3a→3b and 3b→3c (commit + checkpoint this file first).

**3a PROGRESS (2026-07-26, branch `feature/sn-dsa`):**
- ✅ **3a.1 @ `dbdbb2b9`** — `orpheus/derivations/discrete/sn/dsa.py`: Larsen's
  four-step executed SYMBOLICALLY over a general symmetric quadrature; THE MAIN
  THEOREM proven (shared-edge f₁ continuity ≡ (27) with (23a–f), minor-based
  proportionality on the W0=1/W2=⅓ variety, weights solved linearly — no
  radicals); annihilation identities ((14b)'s "3" = 1/W₂); the TWO distinct
  ⅓-mechanisms (Legendre recursion vs W₂) kept separate; (16a–d) derived with
  the identity-verified explicit/OPAQUE-lagged split (the slot assignment IS
  the method — moment-expanding L₀[γψ] silently promotes a lagged term; caught
  by the machinery, documented); updates (28), Marshak (38a) with DERIVED
  (γ_N, W2⁺/W2), one-sided vacuum closure, reflective, DD instance. 8 pins
  green under `-O` (~89 s, shared solve cache).
- ✅ **3a.2 @ `614eee19`** — the one-sided (25)/(26) closures verified; the
  numeric builder `build_consistent_dd_system` (transcribes NOTHING; the ONE
  Σw=2→ω=w/2 convention boundary, guarded); **THE PARENT TIE**: the
  edge-eliminated DD parent ≡ `assemble_ordinate_blocks` entry-for-entry at
  1e-11 per (ordinate, group) on the het non-uniform slab — production
  convention confirmed as the DENSITY form (μ/h)Δψ + σ_Tψ̄ (diamond.py kernel).
  5 tie/structure tests green; CLI pyright 0.
- ✅ **D2 DONE 2026-07-26** — report `.claude/plans/archive/dsa_d2_characterization.md`
  (regenerate: `d2_characterization.py`, job tmp). The structural diffs
  CONFIRMED + one found in flight: (1) homes — f₀ on K+1 EDGES vs φ on K cells
  (+trace); (2) removal — consistent ¼(1,2,1)·σ̂_R h mass vs lumped diag (the
  derived off-diagonal flips sign at thick cells — the consistency
  fingerprint); (3) boundary — DISCRETE γ_S4 = 0.260634 vs continuum ¼; (4)
  **the D definition** — landed = outflow Σ_tr (TOTAL P1 out-scatter), the
  derivation (23c) = WITHIN-GROUP P1 self-scatter (coincide only iso). The
  derived operator needs NO harmonic mean anywhere (edge unknowns straddle
  material-homogeneous cells); the landed one does (cell unknowns straddle
  faces; matches `_interior_conductance` exactly). **The measured ρ scan (the
  decisive evidence)**: derived edge — ρ ≤ 0.181 over σ_t h ∈ [0.1, 30] ×
  c ∈ {0.9, 0.99} × {vac/vac, refl/vac} (inside A&L's 0.2247c); **S2 anchor
  ρ ≈ 3e-15 incl. reflective** (K₂=0 one-iteration exactness at machine zero —
  self-verifies T, G, A_edge, the discrete Marshak rows, AND the (28a) update
  in one shot); landed cell — ok ≤ 0.5, degraded at 1, **DIVERGES σ_t h ≥ 2
  (ρ up to 54.7 at c=0.99, σ_t h=30)** — the A&L (3.43)/(3.44) divergent
  class measured on OUR operator, the historical spike's regime.
- ⏸ **R4 RULED (user, 2026-07-26)**: **(a) the DSA correction operator = the
  DERIVED edge-centered consistent system** — 3b builds it in production,
  pinned entry-for-entry vs `dsa.build_consistent_dd_system`; the standalone
  diffusion solver KEEPS its cell-centered RT0/harmonic stencil (two defining
  laws — consistency theorem vs standalone accuracy — NOT a twin path). **(b)
  the production build home = SN-SIDE** (new `orpheus/sn` acceleration home
  consuming SNMesh + quadrature + the foldable σ accessors): the operator's
  coefficients are quadrature-dependent (γ_N, W₂) and scheme-dependent (ρ from
  the WD α) — a property of the SN discretization; `orpheus/diffusion` stays
  untouched. 3b wiring notes: for DD, ρ=0 ⟹ (28a) update = ½(f₀L+f₀R); d₀ =
  raw moment displacement (G carries σ̂_S·h internally); d₁ sources = 0 for
  the P0 arm; the low-order D at (23c)'s own σ_s1 = within-group self-scatter.
- **3a COMPLETE.** ⏸ COMPACTION at the 3a→3b boundary — taken.

**3b PROGRESS (2026-07-26, branch `feature/sn-dsa`):**
- ✅ **The accelerator core built + both postures verified end-to-end**:
  - `orpheus/sn/acceleration/dsa.py` — `DSALowOrderSystem` (the SN-side
    production build, R4b home: per-group (23a–f)/(27)/boundary rows
    vectorized off SNMesh + the foldable accessors' FIRST production
    consumer, admission guards for geometry/scheme/walls/Σw/D-positivity)
    + `DSACorrection` (ONE operator both postures consume; R =
    `integrate_angular`, P = normalized iso injection — NOTHING minted,
    per the anti-mint verdict).
  - `SourceIteration(corrector=…)` — the no-op extension through the
    single generic body (`iteration.py`; the stop-identity corrected-arm
    exemption documented; byte-inert when None).
  - Krylov posture: `_within_group_krylov(corrector=…)` — the #200
    identity replaced by the A&L transport-corrected M = sweep +
    correction-of-sweep (first re-enabled preconditioner).
  - Facade: `solve_sn_fixed_source(acceleration="dsa")` (additive
    opt-in; both inner solvers).
  - **Measured**: SI 204→12 (1g vac), 221→16 (1g refl), 599→33 (2g het
    aniso, both walls); Krylov 119→10 / 158→11 / 285→17 / 300→19; FP
    identity ≤ 1e-9 everywhere.
- **TWO consistency discoveries (the derivation caught real physics):**
  1. **The trace arm is load-bearing under lagged reflection** — bulk-only
     correction DIVERGES on every reflective regime (measured res ~1e+35;
     the production splitting lags B_a reading the iterate's outflow,
     while Larsen's reflecting row models a current-inflow error
     equation). Fix: the wall-EDGE solutions f0[0]/f0[K] inject
     isotropically into the boundary-trace displacement — the low-order's
     own edge unknowns ARE the trace correction. Regression-pinned
     (mutation teeth: bulk-only must diverge).
  2. **σ_s1 enters the low-order ONLY when the sweep retains ℓ≥1**
     (`scattering_order` threads into the build) — consistency is with
     the ITERATED (truncated) operator, not the mixture data.
- ✅ The typed algebra extended minimally: the correction IS a
  displacement (torsor `ψ ⊕ Δψ`); `AngularDisplacement.integrate_angular
  → ScalarDisplacement` as the tangent map over the ONE hoisted
  reduction body (`AngularField._integrate_angular_values`) — no twin.
- ✅ Gates: `tests/sn/acceleration/` 18/18 green under `-O` — the
  production tie (≡ the derivation reference builder, atol 1e-15, both
  BC variants + the solve), admission teeth, D7 (hand-posed R
  conservation), D8 (frame-row/tangent-map/round-trip pins), D3/D4
  (FP-invariance, aniso ℓ≥1 het 2G, vacuum AND reflective, both
  postures), D6 (correction→0 exact), sign-flip + trace-zero mutation
  teeth (residual-VALUE witnesses — an iteration-count bar is
  off-by-one-prone: diverging runs report max_inner−1). CLI pyright 0.
- **Plan deviations (R4-driven, for the enforcer + 3c docs):**
  `as_dsa_source` NOT minted (a third moment-0 spelling — Smell 16; the
  correction consumes displacements; `field_algebra.rst:528` rewrites to
  the landed truth in 3c). The P7b `full_field_space` trigger does NOT
  fire (no DiffusionMesh in the loop after R4a; the deferral note's
  "#2" expectation goes stale — 3c documents).
- **NEXT**: elegance-enforcer review of the 3b surface → fixes → commit
  → 3c (task #28: rate tier D11–D13 + S2/ρ/thickness/c→1/negative
  control + docs + close-out).

**3c PROGRESS (2026-07-26, branch `feature/sn-dsa`):**
- ✅ **The rate/stability tier built + measured**
  (`tests/sn/acceleration/test_dsa_rate.py`, 69 gates; the evidence
  pack `.claude/plans/archive/dsa_rate_characterization.md`, parts A–H):
  - **D11** — production S8 ρ_est 0.176–0.180 at c = 0.9 (bound
    0.2022); the spec's pre-measurement ±0.03 band replaced by the
    honest split: one-sided Fourier bound + measured floor + the
    plain-SI estimator-honesty control (ρ_est(SI) = 0.894/0.903 ≈ c).
    qa ratified the deviation (a two-sided band would be reference
    contamination — the discrete ρ sits BELOW the continuum sup).
  - **D12** — thick fully-reflective c → 1 converges FLAT (n = 21
    over σ_t·h ∈ {1,5,20} at c = 0.99); ρ band ≤ 0.35. The
    production reflective rate is LAG-limited (~0.28–0.31 — the
    Jacobi wall gain), NOT Fourier-limited; the fully-coupled
    refl/refl matrix certificates (qa IMPROVE-1) pin the operator
    healthy at the gated BC: ρ(DSA) = 0.191/0.024/0.050.
  - **D13** — CI count gates (speedup floor 225→15/249→20;
    c-independence 2110→16/2554→21; Krylov 195→11/218→12) + the
    slow grid on the SCALE-FREE 10-decade metric. Two corner
    findings: the c→1 FP FLOOR (absolute tol meets the double-
    precision floor at flux scale 1/(1−c) — callers must scale
    inner_tol) and the DOUBLE-LAG reflective mode at σ_t·h = 100
    (ρ_prod ≈ 0.745, c-independent, operator certified healthy —
    improvement issue files at close).
  - **S2 exactness** — n == 2 (vacuum, machine-zero landing 3.2e-15)
    + exactly +1 per lagged reflective wall; qa mutation-verified
    the anchor across four corruption classes (Marshak/sign/G/
    update all break it by 13+ decades) — the de-facto
    M-D12-BOUNDARY-ROW catcher.
  - **The negative control** — WD(a≠½) + the PRODUCTION low-order
    reproduces Table II (a = 0.75, c = 0.99: ρ = 1.44/1.78 at
    σ_t·h = 10/30; consistent flat ≤ 0.181); in-file independence =
    the a = ½ composite S2 machine-zero anchor.
  - **ERR-070 FILED** (the σ_r-fold, #215's class): measured 43.2%
    FP shift (11.4%/43.2% per group) on the vacuum+het anisotropic-
    flux config; identically zero on the reflective-isotropic box
    (the designed-green degeneracy, qa-verified at 1.2e-11).
    Catchers: `TestSigmaRFoldCaught` (the folded-Mixture realization)
    + the **D10 routing sentinel** (AST fence: the foldable
    accessors' production consumers ≡ {definition site, split layer,
    the DSA build}; planted-consumer tooth).
  - **D9 no-masking** landed (seeded ×1.05 scatter bug: DSA
    reproduces the same wrong FP to 1.4e-11; activation 1.05e-2).
- ⏸ **R5 RULED (user, 2026-07-26)**: **(a) LD arm DEFERS** — the
  measured evidence: the LD iterate cannot type-admit the arm-1
  restriction (moment-0 is (ng,K,2) — the pairing is UNSPELLABLE
  without the M4S reduction, which IS the arm-2 build), and the WD
  family shows what partial pairings do when spellable (Table-II
  divergence). Follow-up issue files at close. **(b) P1-DSA WIRES
  NOW** — the d₁ moment-pair arm:
- ✅ **The P1-DSA (d₁) arm landed** (`dsa.py`): `moment1_update` =
  the proven (28b) at ρ = 0 (`f₁ = −(D/h)Δf₀ + a·d₁`); the
  restriction pair (the frame's ℓ=1 row w·μ — SH slab component ≡ μ
  bit-exact); the injection = Larsen (33) `Ψ = f₀ + (μ/W₂)f₁`;
  gated on `scattering_order ≥ 1` (the SAME consistency-with-the-
  iterated-operator rule as the data row; so = 0 byte-identical).
  The trace arm stays ℓ=0 BY THEOREM ((39) wall-edge f₁ = 0).
  **Measured: the ladder FLATTENS 24/39/86 → 14/15/15 (ρ 0.175
  flat, the A&L band) and the ℓ=1 S2 system lands at machine zero
  (5.4e-15)** — S2's angular space IS span{1, μ}, so the exactness
  anchor verifies the whole d₁ convention chain in one number.
  Gates: `TestP1DSAArm` (S2-so=1 exactness, the P0-forced tooth,
  the ladder pin, the moment-pair object pins).
- ✅ **ERR-071 FOUND + FIXED** (the third consistency discovery,
  Krylov-flavored): the composite sweep inverse DROPPED the rhs's
  outflow-trace rows (buffer seeded, march clobbered) — exact on
  every physical rhs (outflow rows ≡ 0 there), SINGULAR on the full
  composite space. The P1-DSA GMRES posture excited the kernel
  (‖Mq‖/‖q‖ = 1e-15 on a pure outflow-row residual; GMRES stalled
  at O(1) true residual; **the end-of-solve certificate made the
  catch**). Fix: one post-march restore in
  `streaming.py::_solve_timed_full_field` (`−=` — the forward's row
  is `streamed − ψ_out`, sign pinned by the new round-trip gate);
  bit-inert on all physical paths. NEW gate:
  `tests/sn/operators/test_sweep_inverse_identity.py` (the
  `A∘A⁻¹ ≡ I` full-space round-trip + the pure-outflow leg +
  selector-emptying tooth; `catches("ERR-071")`). Measured after:
  2G het ℓ≥1 Krylov+DSA **287→12 / 305→16**.
- ✅ Reviews: enforcer PASS (2 IMPROVEs applied: the lru_cache/
  monkeypatch autouse guard; the lumped-mutation cross-ref) + qa
  SUPPORTED (no blockers; IMPROVE-1 = the coupled refl/refl
  certificates, applied; NITs applied). Enforcer DELTA review
  (d₁ arm + ERR-071 fix): PASS — must-fix (sign-convention doc
  harmonization to `streamed − ψ_out`, the authoritative forward
  at loss_representation/__init__.py:1164) + 2 NITs (the ℓ=1 row
  sourced from `angular_frame(1).table[:,1,1]`; the cache-guard
  docstring scoped) — ALL APPLIED 2026-07-27, re-verified green +
  pyright 0. The enforcer confirmed the restore's PLACEMENT
  (`_solve_timed_full_field` is the inverse's one assembly site)
  and the restore's sign against the forward. qa delta review of
  (d₁ + ERR-071) PENDING — runs after the root fix below.

⏸ **R6 RULED (user, 2026-07-27): the ERR-071 resolution is the ROOT
FIX — the honest full-space inverse. Campaign extension accepted.**
(The bounded preconditioner-local alternative was REJECTED.) The user
directed: update this plan, COMPACT, then execute. This block is the
post-compaction execution spec.

### R6 root-fix worklist — ✅ DONE 2026-07-27 (all 5 items + 2 discoveries)

**Executed as specified; measured outcomes:**

1. ✅ `solve_transpose` symmetric completion (`E_out` diagonal ⟹
   symmetric ⟹ `(Aᵀ)⁻¹ = S_oldᵀ − E_out`, the same one-site restore
   on the reverse-scan's output boundary) — transpose file 18/18,
   the 3 G3 reciprocity gates green.
2. ✅ The component idiom + the 2-D production path shared ONE root:
   the trace→source ROLE CONFLATION at
   `AngularBoundarySourceSink.from_mesh(<flux trace>.values.copy())`.
   The named projection ALREADY EXISTED —
   `AngularBoundarySourceSink.prescribed_inflow` (inflow slots only;
   outflow rows unrepresentable, Pattern 4) — the fix routes both
   sites through it: the `solve_sn` eigenvalue-finalize
   (solver.py:2174ff; the SOLE production conflation — audit: the
   driver folds build `zeros_on` + inflow-only B gains, measured 0.0;
   the :2747 composite arm is caller-owned source-typed data) and the
   `sweep_once` test helper. Q/Σt + 2-D reflective balance green.
3. ✅ Identity gate extended: parametrized {slab_vacuum,
   slab_reflective, cyl_product} × {random-composite, pure-outflow}
   + the sentinel tooth — 7/7.
4. ✅ **Discovery (cylinder honest scope)**: on a product quadrature
   the degenerate pure-azimuthal rows (μ_r = 0, excluded from BOTH
   selectors) are FREE DOFs — forward = structural zero row, inverse
   = identity passthrough (#284); the gate asserts both halves
   explicitly. NOT a partial-inverse regression.
5. ✅ **Discovery (the scheduled sibling)**: the full `tests/sn` net
   caught `test_w2_round_trip_machine_precision` — the W2 fixture was
   a hidden consumer of the old clobber (`LC.solve(random-full-trace)`
   only "trace-consistent" because the drop erased outflow rows).
   The SCHEDULED walk (G-S `M = (L+C)−B_lower`) is exact ONLY on the
   source subspace {y_out = 0}: the mid-march reflect consumes
   un-restored streamed values. Resolution: fixture fixed
   (inflow-only rhs), honest domain documented in
   `ScheduledInvertibleOperator._solve_timed_full_field`, and an
   off-domain characterization PIN added (measures the B(y_out) gap;
   REDS the day a per-group restore completion lands — the flip-me
   tripwire). NO value-guard (exact-zero rejects FP-dust round-trips;
   thresholds are arbitrary); the production catcher for a future
   off-domain consumer is the end-of-solve certificate itself.
   Catalog ERR-071 updated with parts (1)–(5).

**Re-run evidence (2026-07-27)**: operators+derivations+acceleration
fast 2461 passed/0 failed (66 min); tests/sn `-m "not slow"` 2318+1
(the +1 = W2, fixed after, its file now 8/8; production delta since
that run is docstring-only); slow acceleration 33/33; identity 7/7;
transpose 18/18; CLI pyright 0 errors on all touched files
(streaming.py, solver.py, scheduled_invertible.py, _test_helpers.py,
both gate files; `_test_helpers.py`'s 10 CLI errors pre-date the edit
— stash-verified identical count).

**NEW close-out addendum**: the issues step gains one more filing —
the scheduled-walk full-space completion (per-group restore ordering;
unblocks a future G-S-preconditioned Krylov posture; the W2 pin is
the tripwire) — module:sn, type:improvement.

#### The original worklist spec (as executed)

**The contract being completed**: `InvertibleOperator.solve` /
`SweepOperator.apply` treat `rhs.boundary` as the ALGEBRAIC trace rhs
of the composite system — inflow rows = given data (unchanged);
OUTFLOW rows = the defect rhs, honored as `ψ_out = streamed − rhs_out`
(the forward's row is `streamed − ψ_out`; signs harmonized everywhere
2026-07-27). The `.solve` half is LANDED (streaming.py
`_solve_timed_full_field`, the post-march `−=` restore over
`angular_trace.outflow_indices_for_face` — enforcer-verified; KEEP).
Every caller that used `rhs.boundary` as an ITERATION-STATE seed
(stale outflow rows, pre-fix harmlessly clobbered) must now pass an
inflow-only SOURCE trace — that conflation of roles is exactly what
ERR-071 documents.

**The 5 red tests ARE the worklist** (deterministic reproductions):

1. `tests/sn/operators/test_loss_transpose_solve.py::
   test_g3_full_field_solve_reciprocity[slab|cyl_product|cyl_ls]` (3)
   — **`solve_transpose` needs the SYMMETRIC completion**: writing
   A⁻¹ = S_old − E_out with E_out the diagonal partial identity
   (rhs outflow rows → outflow slots), E_out is symmetric ⟹
   (Aᵀ)⁻¹ = old_solve_transpose − E_out: the SAME one-site restore in
   `solve_transpose`'s output-boundary assembly (streaming.py:1010ff;
   the per-face membership loop handles curvilinear xmax-only faces).
   The reciprocity gates ⟨A⁻¹q,p⟩ = ⟨q,(Aᵀ)⁻¹p⟩ on random FULL
   composites are the committed catchers.
2. `tests/sn/operators/test_solver_components.py::
   TestQuadratureWeightConservation::
   test_homogeneous_scalar_flux_equals_Q_over_sigt` — the component
   iteration idiom: a persistent full-trace `boundary_flux` passed to
   `sweep_once` each iteration (stale outflow rows now subtracted ⟹
   diverges 1e+46). Fix the idiom: after
   `_reflect_outflow_into_inflow`, zero the outflow rows (or build a
   fresh inflow-only trace per iteration — what the 1-D production
   fold does; measured outflow-rhs ≡ 0.0 there). Check whether
   `sweep_once`/`_reflect_outflow_into_inflow` are test helpers or
   production surface — its comment says "mirroring production:
   _solve_fixed_source_si / solve_sn".
3. `tests/sn/operators/test_bc_extraction_2d.py::
   TestBoundaryResidual2DDrivesToZero::
   test_reflective_boundary_balance_at_convergence` — **the 2-D
   PRODUCTION path** (`solve_sn`, 2-D reflective Krylov eigenvalue):
   keff still exactly k_inf 1.875, but the converged trace balance
   defect = 8.8e-2 ⟹ somewhere the 2-D drive feeds `solve` an rhs
   whose boundary carries stale outflow rows. AUDIT METHOD (proven
   this session): monkeypatch-instrument `_solve_timed_full_field`
   recording max |outflow-row rhs| per call during THIS test → find
   the caller → make that builder emit an inflow-only source trace.
   Candidates: the 2-D SI/Krylov rhs fold, production
   `_reflect_outflow_into_inflow` usage, the windowed
   `WindowedSweep = P @ A.inverse()` path.

**Order**: (1) transpose completion → reciprocity green; (2) the
component idiom fix; (3) the 2-D audit + fix → balance green; (4)
extend `tests/sn/operators/test_sweep_inverse_identity.py` with a
reflective-wall leg and (if cheap) a curvilinear leg; (5) re-run:
the 5 tests, tests/sn/operators, tests/sn, the acceleration battery
(58 fast + 33 slow), derivations; CLI pyright on every touched file.

**Verification discipline**: keff-level greens are NOT sufficient
(the 2-D failure kept keff exact while the trace drifted — a Mode-12
lesson: the balance gate caught what the eigenvalue could not); the
trace-balance and reciprocity gates are the load-bearing catchers.

### Post-root-fix close-out sequence (unchanged, runs after R6)

1. **qa delta review** (math-bearing, NOT yet reviewed: the d₁ arm +
   the FULL ERR-071 root fix) — findings fixed in-session.
2. **The gated commits** (each green at its commit): (i) fix(sn):
   the honest full-space sweep inverse [streaming.py solve +
   solve_transpose + call-site fixes + the identity/reciprocity
   gates + ERR-070/071 catalog]; (ii) feat(sn): the P1-DSA d₁ arm
   [dsa.py]; (iii) feat(test)/docs: the 3c rate tier + the evidence
   packs + method.py note + this roadmap; (iv) chore(agents): the
   agent-memory memos.
3. **Docs close-out** (archivist dispatch): the NEW theory page
   `docs/theory/methods/sn/` (acceleration/DSA — Fourier story, the
   derived-consistency algebra, R4/R5/R6 rulings, the D2 + rate
   evidence-pack tables as teaching artifacts, the THREE consistency
   discoveries incl. ERR-071's certificate catch, honest scope);
   verification/sn.rst DSA section; refs.bib keys (Alcouffe1977,
   Larsen1982a, AdamsMartin1992 — AdamsLarsen2002 exists; standing
   Zotero back-port note); the content-drifted diffusion_1d.rst:517
   xref → field_algebra; slab_one_group.rst σ_r cross-ref; then
   verifies() markers across the DSA battery (labels minted on the
   page) + matrix regen. MIGRATE the gate-file docstring facts to
   the page (enforcer NOTE f1: migrate, don't duplicate).
4. **Issues**: `Closes #2`; #200 disposition (the DSA preconditioner
   IS the first re-enabled M); #215 disposition (fence + fenced
   consumer + ERR-070 filed); #22 comment (N-D interaction); FILE:
   the LD arm (draft: job-tmp `issue_ld_arm.md`), the reflective
   lag-mode improvement (`issue_lag_mode.md`), 2-D Cartesian DSA
   (`issue_2d_dsa.md`), and the frame-backed `angular_moment(ℓ)`
   field reduction (enforcer's architectural NOTE — module:sn,
   type:improvement; makes d₀/d₁ symmetric frame-backed reductions).
5. **Pre-merge**: full serial tree `-m "not slow"` (~52 min,
   UNCONTENDED) + Sphinx `-E -W` + `python -m tests._harness.audit`
   + pyright floor (#288 = the accepted 1). Then ff-merge to main +
   branch delete + close #2 — ON THE USER'S GO (never push without
   the explicit ask).
6. **Campaign close-out**: roadmap Phase-3 terminal block; memory
   (the stencil/DSA project memory → terminal; MEMORY.md index);
   lessons if any (the ERR-071 diagnosis chain — component-by-
   component elimination converging on the composition seam — is
   catalog-recorded).

### Tree inventory at this checkpoint (2026-07-27, all UNCOMMITTED)

- Modified: `orpheus/sn/acceleration/dsa.py` (d₁ arm; enforcer
  PASSED + NIT applied), `orpheus/sn/operators/streaming.py` (the
  solve-half restore + sign-harmonized docs; enforcer PASSED),
  `orpheus/transport/method.py` (P7b stale-note fix; PASSED),
  `.claude/skills/vv-principles/error_catalog.md` (ERR-070 +
  ERR-071), this roadmap, `.claude/agent-memory/*` (enforcer + qa
  memos).
- New: `tests/sn/acceleration/test_dsa_rate.py` (69 gates),
  `tests/sn/operators/test_sweep_inverse_identity.py` (3 gates),
  `.claude/plans/archive/dsa_rate_characterization.md` (the evidence pack).
- GREEN: acceleration 58 fast + 33 slow; identity gate 3/3;
  derivations 13; CLI pyright 0/0/0 on all touched files.
- RED (the R6 worklist): the 5 tests listed above.
- Ephemeral instruments (job tmp, `/Users/rodrigo/.claude/jobs/
  cdba0ba9/tmp/`): `rate_design_scan.py`, `rate_characterization.py`
  (the evidence-pack generator), `d2_characterization.py`,
  `issue_{ld_arm,lag_mode,2d_dsa}.md` (the filing drafts).
- Standing constraints: never push without the user's ask; scratch/
  stays untracked (never `git add -A`); commit trailers per the
  session convention; Host env `.venv/bin/python`; canonical pytest
  `python -O -m pytest`.

### 3-P0 — dispatches (MUST, before the plan-of-record for this phase) — ✅ DONE 2026-07-26
- **literature-researcher** — check `/Users/rodrigo/git/nuclear/ORPHEUS/scratch/
  literature/` FIRST (user maintains it; all NSE volumes local): Alcouffe 1977
  (consistent DSA for DD), Adams & Larsen 2002 review (fast iterative methods —
  the ρ = 0.2247c Fourier results + partial-consistency taxonomy), Adams & Martin
  1992 (M4S), Warsa–Wareing–Morel 2004 (Krylov-DSA robustness where SI+DSA
  degrades). If a paper is not local: report back, do NOT pivot to a secondary
  source without user approval.
- **test-architect** — the verification plan: FP-invariance gates, spectral-radius
  measurement protocol, mutation matrix. MUST encode the ⚠#215 / vv Mode-9 trap:
  every FP-invariance gate runs an ANISOTROPIC config (vacuum/heterogeneous) AND a
  diagonal cubature — never the isotropic reflective box alone.
- **cross-domain-attacker** — the restriction/prolongation pair as a Petrov–Galerkin
  frame on the angular axis (moment-0 analysis face, isotropic reconstruction face)
  — reuse the frame vocabulary; do not mint an ad-hoc R/P pair.

### 3a — Consistency by derivation (Phase 2's payoff)
Compute the discrete P1/diffusion moment reduction of the ASSEMBLED Cartesian DD
operator (slab first) and compare — as MATRICES — against the diffusion family's
assembly on the same mesh:
- If ≡ (or ≡ up to a characterized boundary row): the landed `LeakageOperator` IS the
  consistent partner — Alcouffe's result recovered computationally. Gate it, done.
- If ≢: **R4 RULED "decide from 3a data"** — bring the matrix diff to the user with
  its STRUCTURE characterized (boundary row only? interior stencil? a scaling?) and
  collect the derived-stencil vs partial-consistency ruling THEN. Do not pre-commit
  either way; both remain cheap once assembly exists.

### 3b — The accelerator
1. **R** (restriction): moment-0 of the typed SN residual (`evaluate_residual`,
   `sn/solver.py`) → diffusion source on the scalar composite; boundary arm = the
   ℓ=0 half-range moment under the SHARED `|Ω·n|w` metric (the #290 ruling-2 seam).
2. **Low-order solve**: the (3a) operator over `DiffusionMesh.from_material_mesh(
   sn_mesh)` — same axes/materials/BCs by construction. Correction-equation BCs:
   vacuum faces → Marshak (𝒜=0; the error problem has zero inflow — the #290
   vacuum law is EXACTLY right), reflective → reflective.
3. **P** (prolongation): isotropic correction onto the φ moments (+ the PG-frame
   shape from 3-P0).
4. **Wiring**: both postures through the existing iteration layer — SI+DSA
   (`SourceIteration` wrap) and Krylov-preconditioned (`KrylovAcceleration`; folds
   #200's intent). The Protocol trigger FIRES: declare `full_field_space` on
   `TransportMethod` (the recorded P7b deferral — first generic consumer arrives).
5. **Scope**: slab/Cartesian DD = arm 1. **LD arm 2 gated on measured need** —
   **R5 RULED "decide from measured table"**: measure LD's ρ under the arm-1
   low-order operator across the c→1 × optical-thickness sweep, present the table
   to the user, collect the build-LD-consistent ruling then (M4S / DG-diffusion
   territory if triggered) — deliberately demand-driven, per defer-until-
   consumption. Curvilinear DSA: OUT (blocked on #282; documented seam).

### 3c — Gates
- **FP-invariance (Mode 9, #215-proof)**: converged accelerated flux ≡ plain-SI
  fixed point to solver tolerance, on vacuum+heterogeneous anisotropic config AND a
  diagonal cubature; separately from any rate claim.
- **Spectral radius**: measured ρ vs the ~0.2247c Fourier bound on the classic
  homogeneous-infinite/periodic problem; degradation curve vs cell optical thickness
  (the consistency stress axis); Krylov iteration counts vs SI+DSA vs plain SI —
  the c→1 table IS the teaching artifact (demo + theory page).
- **Correction→0 gate** at convergence (the correctness-safe-by-construction pin).
- **R/P object gates** (Mode 12): conservation `⟨1, R r⟩ = ⟨1, r⟩`-family +
  adjoint-consistency of the pair; mutation-verified.
- **No-masking control**: a seeded transport-operator mutation must converge (fast)
  to the SAME wrong answer with DSA on — the accelerator must not launder bugs.
- ⚠ within-group solves must NOT route through the #215 σ_r-fold spelling.
- Full #290-style close-out: theory docs (`discrete_ordinates.rst` DSA section with
  the Fourier analysis + literature eq numbers; cross-refs to `diffusion_1d.rst`'s
  DSA-seam section), dev-history rows, `Closes #2`, dispositions for #200/#215,
  comment #22 (N-D interaction).

⏸ **C3** + campaign close-out.

---

## Rulings — COLLECTED 2026-07-03 (execution start)

- **R1 (#280 sequencing)**: **BETWEEN Phases 2 and 3** — #280's orientation×kernel
  walk unification runs as its own branch (Phase 2.5) after assembly lands and
  before DSA, so the accelerator wires into the UNIFIED walk. (Deviates from the
  recorded recommendation "after the campaign": the user chose the cleaner
  consumption order for DSA over campaign brevity.)
- **R2 (assembly spelling)**: **separate `assemble()`** → sparse assembled operator;
  `as_matrix` keeps its dense Mat-functor semantics and DELEGATES to densified
  assembly when available; probing stays the fallback for assembly-less operators.
  The parked #213 capability-morphism redesign is not blocked on.
- **R3 (Phase-1 CP/MoC arms)**: **include, drop on friction** — attempt the CP/MoC
  rewires after SN+diffusion land; on substrate friction, file a labeled follow-up
  instead and close #259 scoped to sn/numerics/diffusion.
- **R4 (3a outcome)**: **decide from 3a data** — deliberately deferred; bring the
  derived-vs-landed matrix diff (structure characterized) to the user at 3a.
- **R5 (LD arm criterion)**: **decide from measured table** — deliberately deferred;
  present the measured-ρ c→1 × optical-thickness sweep at Phase 3, rule then.
  **RULED 2026-07-26 (3c)**: (a) the LD arm DEFERS — the measurement came back
  STRUCTURAL (the LD iterate cannot type-admit the arm-1 restriction; the pairing
  is unspellable without the M4S reduction, which IS the arm-2 build; the WD family
  shows Table-II divergence when partial pairings ARE spellable) — follow-up issue
  files at close; (b) **P1-DSA wires NOW** (the d₁ moment-pair arm — landed;
  the ladder flattens 24/39/86 → 14/15/15).
- **R6 (ERR-071 resolution, ruled 2026-07-27)**: **the ROOT FIX — the honest
  full-space sweep inverse** (`rhs.boundary` IS the algebraic trace rhs; every
  iterate-state-seed caller converts to inflow-only source traces; `solve_transpose`
  gains the symmetric completion). Campaign extension accepted; the bounded
  preconditioner-local alternative REJECTED. Execution spec = the "R6 root-fix
  worklist" block in the Phase-3 3c PROGRESS section.
- **R6 (2a delegation)**: **full surgical** — the mechanical relocation is main-agent
  work too; no delegation exception.
- **R7 (n2n convention — from the P0 finding, ruled 2026-07-03)**: **operator-
  consistent everywhere** — νΣf-only numerator; the (n,2n) gain enters through the
  loss-side denominator (the CP/diffusion convention, L1-pinned vs a dense
  eigensolver). MoC's L0 expected value is RE-DERIVED to the posed problem's
  eigenvalue — principled re-baseline (the old reference encoded the estimator's
  convention, not the eigenproblem's). Rationale: every inner solve poses the
  eigenproblem with ONLY fission scaled by 1/k; the SN/MoC numerator spelling equals
  that eigenvalue only when Σ₂=0 or k=1 (same failure class as #291).
- **R9 (residual package fate, ruled 2026-07-04)**: **keep the `sn/spatial/` name for
  the residual this phase** (`pole_angular_closure`, `psi_half_angle_seed`,
  `sweep_cache`, `scan`, `pairing` — sweep-walk + angular machinery; the name goes stale
  but Phase 2.5's #280 walk unification owns the sweep-layer layout, and renaming now
  would churn the same modules twice, incl. 129 stay-side docs refs).
- **R8 (KEigenvalue seam shape, ruled 2026-07-03)**: **hardwire defaults, drop
  kwargs** — the injection kwargs/aliases/default-factories retire; KEigenvalue
  keeps `compute_keff`/`compute_production_rate` as plain methods (the leakage-
  inclusive Rayleigh spelling `Σ(Fψ)/(Σ(Aψ)−Σ(Sψ))` is its fixed estimator); the
  `test_iteration.py:718` injection test rewires to assert adapter-k ≡ solve_sn-k
  at the fixed point directly (the theorem: all consistent estimators agree at a
  converged eigenpair). The protocol is the wiring point — production self-
  implements BY DESIGN; the seam is dead by design, not by missed wiring, and
  post-unification an injection could only introduce an inconsistent functional
  (Pattern-4 retire).
- **R10 (#282 scope, ruled 2026-07-04)**: **FULL route (a) inside Phase 2.5** —
  the in-pass direct starting-direction solve AND the typed per-level ψ(·,μ=−1)
  carrier augmentation both land this phase (unblocks the spherical
  triangular-G2 gate and the curvilinear transpose oracle this campaign). The
  fix runs ON the unified substrate — execution order re-sequenced:
  S0 → 2.5a (bit-identical; the lagged treatment rides the relocation) → 2.5d
  (the #282 fix; principled re-baseline on curvilinear, Cartesian bitwise
  unmoved — the seed exists only on curvilinear) → 2.5b (the reverse-scan,
  built once against the FIXED formulation) → 2.5c → 2.5e.
- **R11 (adjoint-inverse predicate honesty, ruled 2026-07-04)**: **(b) with
  (a)'s spelling** — `_AdjointOperator.is_invertible` + `.inverse()` land
  TOGETHER with `SweepOperator.apply_transpose` (honest on arrival); the
  predicate is spelled generally: `inner.is_invertible and
  inner.inverse().is_adjointable`.
- **R12 (seed-carrier framing, ruled 2026-07-04 at 2.5d open)**: **Frame A**
  (the cross-domain-attacker recommendation) — the `starting_direction` block's
  presence is keyed PER LEVEL by the STRUCTURAL predicate "μ_start ∉ the
  level's μ-nodes" (⟺ τ_raw ≠ 0), NOT by geometry name: sphere-GL carries it
  (1 level); cylinder-product does NOT (its starting direction coincides with
  its first azimuthal ordinate — the #229 τ-clamp fact; a block would be a
  dead rank-duplicate of ψ₀); Cartesian never. Supersedes the R10 letter's
  uniform "per-level" wording and §16.A A1's cylinder-HAS-block leg. Ghost
  (all-zero) metric — a structural fact ((1−μ²)|₋₁ = 0 ≡ α_{1/2} = 0), never
  a fabricated volume weight.
- **R12a (predicate spelling refinement, ruled 2026-07-04 at d1 execution)**:
  presence per level iff the level's FIRST-ORDINATE **raw M-M weight
  τ_raw ∈ (0, 1) EXCLUSIVE** — "the recurrence consumes independent
  starting-direction state". The R12 letter's μ-node-membership spelling and
  its claimed ⟺ with τ_raw ≠ 0 are empirically FALSE on level-symmetric
  cylinder rules: their duplicate-η nodes collapse the midpoint edge onto η₀,
  giving **τ_raw,0 = 1.0 bit-exact** (μ_start ∉ nodes, gaps ~0.1) — the
  seed's only consumption path, the recurrence weight (1−τ₀), vanishes (the
  direct read is always α_{1/2} = 0), so the letter would mint DEAD blocks
  (and it is WHY the measured cyl-LS seedΔ is 0.0-bit). The bit-exact
  trichotomy: τ_raw,0 = 0 (product rules — seed ≡ ψ₀ rank-duplicate, #229) /
  = 1 (LS rules — dead) / ∈ (0,1) (sphere-GL ≈ 0.39–0.42 — genuine state).
  Single-sourced from the NEW `morel_montry_tau_raw_per_level` (the clamp's
  own input, split out in d1; clamped producer bit-identical); implemented as
  `SNMesh.starting_direction_levels` (+ `_SEED_TAU_EPS = 1e-12` FP-noise
  guard — never decides a real case); the A1 gate's cyl-LS row is the
  discriminator pin. Every named R12 instance decides identically.
- **d3 HAZARD banked at d1 (product-cylinder ig-consumption)**: TODAY the
  cylinder solve seeds its per-level thread from the extrapolation ON THE
  INITIAL GUESS (LR `_run` builds the ctx from ig for every curvilinear
  level; product rules: stencil t = 0 exactly ⟹ seed = the ig's ψ₀ row —
  formally a lag, harmless at the fixed point, and the product-cyl MATVEC is
  lag-free/triangular since apply extrapolates from the input). d3's
  strategy-zoo retirement MUST preserve the non-carrying levels' seed
  data-flow BIT-EXACTLY (solve: from ig; apply: from input — inline the t=0
  stencil read, do not "simplify" to the marched ψ₀) or the §16.D cylinder
  assert-unmoved gate HALTS the phase. The LS-cylinder seed value is dead
  (τ₀ = 1), so any spelling is bitwise there; the PRODUCT cylinder is the
  binding row.
- **R13 (seed outer BC, ruled 2026-07-04)**: **mint the corner slots NOW** —
  the (R, μ=∓1) corner pair lands this phase (inflow corner = given-data /
  identity row; outflow corner = defect row; B's reflective arm maps
  corner-out → corner-in). Forced consequence: the OUTWARD starting-direction
  leg becomes STATE (the corner outflow row must be linear in state),
  pole-continued from the inward leg — the seed leaf carries per level BOTH
  directions' cells + corners on one flat backing (the `AngularBoundaryFlux`
  precedent), all at zero metric weight.
- **R14 (q½ source fold, ruled 2026-07-04)**: **the full (−1)^ℓ fold helper
  from day one** — one helper `Q̄(μ=±1) = Σ_ℓ (2ℓ+1)/2·Q_ℓ·P_ℓ(±1)` with the
  §16 B2b 2-term Legendre sign pin LIVE (production reach may stay ℓ=0;
  the anisotropic case is manufactured before it is needed).

## Issue map

| Phase | Drives | Touches / comments | Follow-ups it may file |
|---|---|---|---|
| 1 | #259, #291 | #270 (CP/MoC arms) | CP/MoC estimator rewire (if R3 defers) |
| 2 | #272, #253, (#158) | #284 (discharged), #282 (characterized), #200, #256 | curvilinear assembly (post-#282) |
| 2.5 (R1) | #280 | #282 (structure fix candidate) | — |
| 3 | #2 | #200 (posture folded), #215 (disposition), #22, #293 | LD low-order op (per R5 table), curvilinear DSA |
