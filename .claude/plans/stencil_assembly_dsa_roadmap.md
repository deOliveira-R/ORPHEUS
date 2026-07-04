# Stencil-assembly + DSA roadmap — k-estimator → spatial reification → DSA (#2)

**Status: Phases 1 AND 2 COMPLETE, both merged to main (P1 @ `a4952c3`; P2 =
the `f6079be…` chain ff-merged 2026-07-04 — trust `git log`, this note is a
snapshot). Full-tree merge gate 5990/0 not-slow serial. Rulings R1–R9
collected. Phase 2.5 (the #280 orientation×kernel walk unification) is OPEN
on branch `refactor/sn-walk-unification` (2026-07-04; P0 findings + the
plan-of-record in its section below). Then Phase 3 (DSA #2, opening with
the 3-P0 dispatches — the literature brief checks `scratch/literature/`
FIRST).**

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

**2b STATUS — COMPLETE (2026-07-04):**
- Commit chain: `83a0db7` (numerics home) → `3238651` (diffusion emitters) →
  `ba26644` (SN closure-generic walk + gates). Design contract:
  `.claude/plans/assembly_mode_crosswalk.md` (the L16-cited convention file).
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

⏸ **C2.5** after 2.5d (post-fix, pre-reverse-scan) or at the 2.5b→2.5c seam,
whichever the session boundary hits first.

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

## Issue map

| Phase | Drives | Touches / comments | Follow-ups it may file |
|---|---|---|---|
| 1 | #259, #291 | #270 (CP/MoC arms) | CP/MoC estimator rewire (if R3 defers) |
| 2 | #272, #253, (#158) | #284 (discharged), #282 (characterized), #200, #256 | curvilinear assembly (post-#282) |
| 2.5 (R1) | #280 | #282 (structure fix candidate) | — |
| 3 | #2 | #200 (posture folded), #215 (disposition), #22, #293 | LD low-order op (per R5 table), curvilinear DSA |
