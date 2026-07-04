# Stencil-assembly + DSA roadmap — k-estimator → spatial reification → DSA (#2)

**Status: IN EXECUTION — Phase 1 COMPLETE, merged to main @ `a4952c3`. Phase 2 opened
2026-07-04 on `refactor/spatial-promotion-assembly` (off main @ `a4952c3`); the three
2-P0 dispatches are COMPLETE (distillation in the Phase-2 section); R9 ruled. Rulings
R1–R9 collected.**

Campaign chain: **Phase 1** k-estimator unification (#259 + #291) → **Phase 2** spatial
substrate promotion + the ASSEMBLY third mode (#272 + #158 + #253, user ruling below) →
**Phase 2.5 (per R1)** the #280 orientation×kernel walk unification, its own
branch/campaign → **Phase 3** DSA for SN (#2, folding #200's Krylov posture; wires into
the UNIFIED walk). Each phase is a separate branch + ff-merge, in order; later phases
consume earlier ones.

## Provenance (how this plan was decided)

Follows directly from campaign #290 (merged @ `3a19133`; see
`.claude/plans/diffusion_integration_290.md`). The user chose the sequence
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
- Doc drift for close-out: `docs/theory/monte_carlo.rst:281` (phantom MC
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

### 2c — Docs + close-out
- Theory: `loss_representations.rst` / `discrete_ordinates.rst` — the three-modes
  section (solve/apply/assemble as consumption modes of one closure algebra);
  `operator_algebra.rst` — `as_matrix` probing vs structural assembly; dev-history
  rows. Close/comment #272, #253, #158, #284; comment #282 with the
  characterization; comment #200 (blocks now assemblable).

⏸ **C2**.

---

## Phase 3 — DSA for SN (#2; folds #200's Krylov posture; runs AFTER the #280 walk
unification per R1 — the accelerator wires into the unified orientation×kernel walk)

### 3-P0 — dispatches (MUST, before the plan-of-record for this phase)
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

## Issue map

| Phase | Drives | Touches / comments | Follow-ups it may file |
|---|---|---|---|
| 1 | #259, #291 | #270 (CP/MoC arms) | CP/MoC estimator rewire (if R3 defers) |
| 2 | #272, #253, (#158) | #284 (discharged), #282 (characterized), #200, #256 | curvilinear assembly (post-#282) |
| 2.5 (R1) | #280 | #282 (structure fix candidate) | — |
| 3 | #2 | #200 (posture folded), #215 (disposition), #22, #293 | LD low-order op (per R5 table), curvilinear DSA |
