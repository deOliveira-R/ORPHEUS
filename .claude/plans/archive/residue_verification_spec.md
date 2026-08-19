# The #280 residue carve — binding verification gate spec

**Status:** PRE-CARVE verification design (mandatory before implementation).
**Carve family:** normal↔adjoint operator-algebra (surgical; main agent writes,
user steers). This spec SHAPES the implementation — it is not advisory.
**Author:** test-architect, 2026-07-23. Distills lessons L13–L19 (the
operator-taxonomy adjoint family) + the a3 reverse-scan recipe; applies the
`vv-principles` claim-taxonomy, pillar, Mode-10/11/12, and anti-pattern
instruments.

Read alongside: `a3_solve_transpose_verification.md` (§§9–16 reverse-scan),
`.claude/agent-memory/test-architect/a3_reverse_scan_transpose_verification.md`
(the reusable recipe), lessons L17/L18/L19.

---

## §0. Orientation — what is landed, what is the residue, and the ONE cardinal rule

### 0.1 What Phase 2.5 already landed (the substrate this carve completes)

On `main` today the SN adjoint surface is:

- **Two-factor scope gate** `StreamingOperator.is_adjointable =
  type(scheme).has_transpose_kernel ∧ loss_representation.has_transpose_walk`
  (`streaming.py:258-272`). `.H` eager-raises `MissingAdjoint` at construction
  when False (`operator.py:813-856`; Pattern 4 — illegal state unrepresentable).
- **Kernel factor:** DiamondDifference `has_transpose_kernel=True`,
  LinearDiscontinuous `False` (`scheme.py:691`; the base default is `False`,
  opt-in).
- **Orientation factor:** `CumprodScan.has_transpose_walk=True`,
  `ScanMarch.has_transpose_walk = mesh.is_1d`, both DAG-wavefront classes
  (`MovingFrontierWindow` / the windowed policy) `False`.
- The 1-D DD reverse walk (`loss_representation/__init__.py` `visit`/`open_leg`/
  `close_leg` closures, ~3044-3100), the 1-D reverse-SCAN solve
  (`ordinate_scan_transpose`, `scan.py:239`), the swap-law coherence
  (`A.H.inverse() ≡ A.inverse().H`), and the frozen pre-carve walk baselines
  (`walk_matvec_{slab,sphere,cyl}_2g.npz`, DD-only) all landed (2.5a–2.5d).

### 0.2 The residue (what this carve builds)

| Arm | Deliverable | Flag it flips |
|---|---|---|
| **R1a** | Relocate the hand-coded DD 1-D transpose into a scheme-registered kernel PAIR on `DiscretizationSchemeBase`. Bit-identical. | (none — DD already True; `has_transpose_kernel` becomes REGISTRATION-coupled) |
| **R1b** | LD transpose kernel from the UBLD algebra-of-record (`_ubld.py` / SymPy `ld_ubld.py`). | (paired with R1c) |
| **R1c** | LD 1-D reverse walk (moment-tailed cotangent buffers; the slope-fold transposes). | `LinearDiscontinuous.has_transpose_kernel → True` |
| **R2a** | Multi-D Cartesian **full-cochain** reverse walk (DD) — the reverse oracle. | (oracle; predicate stays False until R2b) |
| **R2b** | Multi-D Cartesian **windowed** reverse (DD) — production. | `ScanMarch(2d).has_transpose_walk → True`, `MovingFrontierWindow.has_transpose_walk → True` |
| **R2c** | LD-2D reverse (`n_face_moments=2`). | (rides the R2b flip; LD-2D kernel from R1b) |

**Explicitly OUT of scope** (stays typed + loud-deferred — DO NOT gate as
active, DO gate the refusal stays loud): the `ScheduledInvertibleOperator`
G-S schedule reverse (`sweep_operator.py:185-192`, no consumer), d≥3
interior-axis inflow interleave (`_ubld.assemble_inflow_axis` supports
`axis∈{0,d-1}` only), curvilinear LD (Issue #158/#6).

### 0.3 The claim layer and pillar (the `vv-principles` §1.5 gate) — declared ONCE for the whole carve

Every claim in this carve is an **operator-identity / flux-shape claim**
(`Lᵀ`, `L.H`, `(L+C)⁻ᵀ` reproduce the true transpose/adjoint). **NOT ONE gate
is an eigenvalue claim.** The pillars:

| Reference | Pillar | Proves | Structural independence of the reverse walk |
|---|---|---|---|
| Frozen pre-carve walk baseline (`walk_matvec_*.npz`) | closed-form (captured artefact) | bit-identity of a RELOCATION | ✓ captured BEFORE the code |
| G-adjoint reciprocity `⟨Aψ,φ⟩_G=⟨ψ,A.Hφ⟩_G` | closed-form (algebraic identity) | the DEFINING adjoint property | ✓ arithmetic = FORWARD apply + independent metric; never calls the transpose |
| Dense-Euclidean-`Mᵀ` (`M` column-probed off FORWARD apply) | closed-form | `apply_transpose ≡ Mᵀ·x` / `solve_transpose ≡ solve(Mᵀ,b)` | ✓ the reference NEVER touches the reverse code |
| Assembled-`Mᵀ` (LAPACK `solve_triangular(permutedᵀ,·)` / CSR `Mᵀ@x`) | closed-form | transposed-scan coefficient + triangularity | ✓ LAPACK/scipy shares no code with the ORPHEUS scan (Cartesian only) |
| SymPy `ld_ubld` transpose-oracle (`simplify(VJP − AᵀM⁻¹)==0`) | closed-form (Branch 1) | the LD cell VJP IS the transpose | ✓ symbolic algebra-of-record |
| G1 round-trip `solve_transpose∘apply_transpose=I` | — | necessary-NOT-sufficient (both arms share the walk) | ✗ shares the walk — corroborating smoke only |

### 0.4 THE ONE CARDINAL RULE (Mode-12) — object-level gates ONLY

**NEVER credit any transpose/adjoint gate on a spectral or norm functional.**
`eig(Kᵀ) = eig(K)`, `‖Aᵀ‖ = ‖A‖`, and `k(A) = k(Aᵀ)` by construction — every
eigenvalue/keff/spectral-norm gate is DESIGNED-GREEN on the entire adjoint
mutation class, at every tolerance, in every regime. This is memory L16's
"NEVER `keff(asm) ≡ keff(apply)`" and the #226 step-5b overclaim. **Every gate
in this spec measures the OBJECT** (full-field matvec output, per-element
matrix, bilinear reciprocity residual) — never a functional whose invariance
group contains the transpose. Full Mode-12 analysis: §7.

---

## §1. The normal↔adjoint convention crosswalk + the R1a contract-design feedback

Per the `coding-elegance` Pattern-7 crosswalk discipline, a carve crossing the
normal↔adjoint boundary writes its convention table BEFORE any code. The VJP
contract:

### 1.1 The VJP contract (verified against the landed DD hand-code)

Forward `residual_kernel_batch : (psi_bar, psi_in, s_axes, reaction_xs,
Q_cells) → (residual, psi_out)`. Its transpose kernel is the VJP:
`(res_bar, psi_out_bar) → (psi_bar_cot, psi_in_cot [, source_cot])`.

Confirmed against the landed DD `visit` closure (`__init__.py:3067-3074`) — the
VJP is exactly:

```
psi_bar_cot += 2·psi_out_bar + (S/V)·res_bar      # DD: ∂psi_out/∂psi_bar=2, ∂res/∂psi_bar=S/V
psi_in_cot  += −psi_out_bar − (|μ|A/V)·res_bar     # DD: ∂psi_out/∂psi_in=−1, ∂res/∂psi_in=−|μ|A/V
```

which is `psi_bar += 2·f_bar; f_bar = −f_bar; psi_bar += denom·ob/V; f_bar +=
−|μ|A·ob/V` in the hand-code (`ob=res_bar`, `f_bar=psi_out_bar`). The stated
signature `(residual, psi_out) → (psi_bar, psi_in)` is **correct for the
SPATIAL cell relation**.

### 1.2 CONTRACT FEEDBACK (task requirement #7) — the signature is INCOMPLETE for curvilinear DD

The landed DD `visit` closure emits a THIRD cotangent the stated signature
omits: `numer_bar[leg.mu_level_idx][...] += −ob/V` — the **angular-numerator
cotangent** (`angular_numer_upstream`, the Morel–Montry `ψ_{n-1/2}` half-flux),
threaded downstream by `angular_adjoint`. A relocation that folds the whole
`visit` into a `(psi_bar, psi_in)`-only kernel **drops `numer_bar` and
regresses the curvilinear DD adjoint** — the sphere/cyl frozen baselines
(`walk_matvec_{sphere,cyl}_2g.npz`) would red.

**Required design decision (bind it in the crosswalk):** the registered
kernel-pair VJP covers the SPATIAL relation `(psi_bar, psi_in, source)` ONLY;
the M-M **angular-thread VJP stays on the `PoleAngularClosure`** (its
`cell_contribution` transpose, i.e. `angular_adjoint`), matching the existing
spatial (`cell_balance_for_streaming`) / angular (`closure.cell_contribution`)
separation. This keeps the stated signature honest **and** keeps the sphere/cyl
DD baselines as the regression guard that a dropped angular cotangent reds.
The alternative (a third `angular_cot` return slot on the kernel) re-fuses two
closures the codebase deliberately separates — reject it unless the user
overrules. **Gate:** §3 sphere/cyl frozen baseline + the R1a Mode-11 sentinel
must fire on BOTH the spatial kernel AND the angular closure transpose.

### 1.3 The convention axes this carve crosses

| Axis | Producer (forward) | Bridge (transpose) | Consumer (reverse) |
|---|---|---|---|
| matvec normalization | DD ÷V `residual_kernel_batch`; LD `(Aψ⃗−R)/M_ii` (`linear_discontinuous.py:616-620`) | LD applies **M⁻¹ BEFORE Aᵀ** (`AᵀM⁻¹`, NOT `M⁻¹Aᵀ`) — §4 NEW-algebra (a) | reverse walk |
| moment frame | global→sweep on input (`octant_moment_frame_signs`) | **self-transpose** (diagonal ±1 involution): same vector wraps the VJP on BOTH sides — §4 NEW-algebra (c) | sweep→global on output |
| face addressing | gather at `face_in`, scatter at `face_out` | **gather at `face_out`, scatter-accumulate at `face_in`, domain-boundary in↔out swap** — §6 | reversed walk |
| level order | forward topological | **reversed** levels + reversed frontier | reverse walk |
| slope fold | `_slope_fold` → `(ψ̄, ψ̂)` (`_ubld.py:382-405`) | the slope-fold transposes (`scan_slope_face_source`ᵀ, `scan_reconstruct`ᵀ) — §5 | LD reverse |

Any bug lives at a **Bridge** cell — that is exactly where the mutation matrix
(§11) plants its teeth.

---

## §2. R1a — the DD kernel-pair RELOCATION (bit-identity refactor)

R1a moves the hand-coded DD 1-D transpose (`__init__.py:3051-3075`) into a
registered kernel PAIR on `DiscretizationSchemeBase`; `has_transpose_kernel`
becomes coupled to registration (a registered pair), not a bare declared bool.
**This is a behaviour-preserving refactor** — the `vv-principles`
§bit-identity gate applies (criterion: pure relocation, no numerical-code byte
changes → `np.array_equal` is the RIGHT contract). The task's requirement #3
(the 0-ULP pin set + the Mode-11 hazard) is the whole of this section.

### 2.1 Claim, references, tolerance

- **Claim layer:** bit-identity of a relocation (the lowest layer — no new value).
- **Reference (structurally independent):** the FROZEN pre-carve
  `walk_matvec_{slab,sphere,cyl}_2g.npz` baselines (`_generate_walk_baselines.py`
  CASES; DD-only, `(L+C)`, both orientations, blocks `{fwd,adj}×{bulk,trace}`).
  Structurally independent because captured BEFORE the relocated code (L17 —
  the `array_equal` removal-form canaries are SELF-REFERENTIAL and must NOT be
  the sole gate).
- **Tolerance:** `assert_regression(kind="direct", reduction_depth=1)` = `nulp=1`
  (a single walk pass, no iteration), LAYERED with `-W error::DriftWarning` for
  the strict byte-identity contract R1a claims. NOT `rtol`. Demonstrate the
  DriftWarning escalation with a deliberate ±1-ULP perturbation of the relocated
  coefficient (the L7 "verify WHICH invocation ran" discipline).

### 2.2 The Mode-11 hazard — prove the RELOCATED kernel is on the executed path

A green baseline does NOT prove the relocation happened: if the old hand-code
`visit` closure still runs and the registered kernel is DEAD, the baseline is
green (old code, same value). Two-part liveness proof (requirement #3):

1. **Mode-11 sentinel (NEW, load-bearing).** An in-process autouse/monkeypatch
   WRAP on the NEW registered transpose kernel (the scheme method, e.g.
   `DiamondDifference.residual_kernel_transpose` — name per the carve),
   asserting `count > 0` during the `walk_matvec_*` baseline run. This is
   strictly stronger than a green baseline (vv Mode-11 sharpening: the
   pytest-plugin wrap runs on the SAME object production constructs — a green
   baseline that routes around the new kernel leaves the counter at 0 and reds).
   Since §1.2 keeps the angular thread on the closure, ALSO wrap the closure's
   `cell_contribution` transpose and assert `count > 0` on the sphere/cyl rows
   (`count == 0` on slab — the config-asymmetry that proves the angular wrap is
   real).
2. **Retirement grep (first-class deliverable, coding-standards).** The old
   hand-code (`psi_bar += 2·f_bar`, the inline `denom·ob/V` / `−|μ|A·ob/V`
   pullback) is DELETED, replaced by a `scheme.<transpose_kernel>(...)` call.
   `grep -n "2.0 \* f_bar\|denom \* ob"` in `loss_representation/` returns ZERO
   after the carve. A lingering twin (old code alive, new kernel dead) is the
   exact Mode-11 failure the sentinel + grep jointly close.

### 2.3 Config (the frozen baselines already encode the right config)

slab_2g / sphere_2g / cyl_2g, `(L+C)`, per-group-varying σ_t (`_per_group_sigma`),
NON-FLAT seeded random bulk AND boundary (a zero boundary nulls the boundary
in↔out swap the transpose exercises), plus the seeded ψ½ legs. This config is
already correct — R1a inherits it. Do NOT add a flat-flux row (§9 config-blindness).

### 2.4 Gates (all EXTEND / stay green — no new value)

| Gate | Disposition | What it holds |
|---|---|---|
| `test_walk_matvec_baseline` (slab/sphere/cyl) | EXTEND (stays green, `-W error::DriftWarning`) | the relocation moved no byte |
| `test_removal_form_matvec_sweep::…apply_transpose…bit_identical` (slab/sphere/cyl) | stays green | override-not-leak (self-referential — L17; corroborating only) |
| `test_g_adjoint_reciprocity` (existing rows) | stays green | the composite `.H` canary the relocation rebuilds |
| **NEW** Mode-11 sentinel (§2.2) | NEW | the relocated kernel + angular closure are on the path |
| **NEW** retirement grep test | NEW | the old hand-code is gone |

### 2.5 Teeth (Mode-10 — prove the gate class can red)

- **M-R1a-DEAD:** stub the registered kernel to return zeros (simulate "not
  wired") → the Mode-11 sentinel reds (count>0 fails) EVEN IF a lingering
  hand-code keeps the baseline green — proving the sentinel catches the
  relocation-didn't-happen case the baseline can't.
- **M-R1a-VALUE:** perturb the relocated coefficient (`2·f_bar → 1.9·f_bar`)
  in-process → the frozen baseline reds (nulp) AND reciprocity reds; the
  DriftWarning escalates. (Reverting: monkeypatch, NEVER `git checkout` — the
  file is uncommitted during the surgical carve, process-discipline.md.)

---

## §3. R1b — the LD transpose kernel from the algebra-of-record

The UBLD per-cell matrix is MATERIALIZED (`assemble_ubld`, `_ubld.py:183-257`,
`A = G + F_out + Σ_t·M`) with a SymPy Branch-1 record (`ld_ubld.py`, the
`derive_*` family closing `simplify(diff)==0`). The transpose is
`einsum("...ji,...j")` / `solve(swapaxes(A,-1,-2))`. **The ONLY new algebra is
composition order** — that is where every gate plants its teeth.

### 3.1 The SymPy transpose-oracle (Branch 1 — the algebra-of-record gate)

Add to `ld_ubld.py`, in the same `derive_*` style, a symbolic proof that the
LD cell VJP IS the transpose:

- `derive_d1_transpose_equals_At_Minv()` → `simplify(VJP_matrix − (A.T @ Minv)) == 0`
  where `Minv = M⁻¹` (diagonal, `M⁻ᵀ = M⁻¹`). This pins NEW-algebra (a) (the
  mass-inverse-BEFORE-Aᵀ order) symbolically, structurally independent of the
  numpy code. **This is the R1b keystone** (algebra-of-record; `vv-principles`
  §MMS→closed-form Branch-1 discipline).
- `derive_octant_frame_sign_is_involution()` → `simplify(D_s @ D_s − I) == 0`
  AND `simplify((D_s @ A @ D_s).T − (D_s @ A.T @ D_s)) == 0` (the frame map is
  preserved under transpose — NEW-algebra (c)).
- Test gate: one `@pytest.mark.foundation` test per `derive_*` in
  `tests/derivations/test_ld_ubld_transpose_symbolic.py` (the test count = the
  `derive_*` claim count; `math-origin` promotes carry NO `verifies` label —
  lessons L9). These run under `-O` (`assert result["pass"]` inside a collected
  test file IS rewritten by pytest — Mode-8 note: bare `assert` in `tests/` is
  fine; the SymPy `simplify(...)==0` is the payload).

### 3.2 The dense numeric transpose oracle (Branch-2 cross-check)

`per_cell_solve` on `swapaxes(A,-1,-2)` (the dense d-generic reference) vs the
production LD VJP, over a batched stack: `assemble_ubld` gives `A`; the VJP's
matrix on the moment basis must equal `Aᵀ M⁻¹` (matvec direction) / `M⁻ᵀ Aᵀ`
(as required by the forward normalization). `assert_allclose(atol=1e-13)` (one
`np.linalg.solve` reduction; NOT `array_equal` — the dense solve ≠ the
production Schur reduction tree, L16 sparse-order-≠-apply-order sibling).

### 3.3 The three NEW-algebra gates (each its own mutation)

- **(a) mass-normalization order.** Forward is `res = M⁻¹(Aψ⃗−R)` → VJP is
  `AᵀM⁻¹·res_bar` (M⁻¹ applied FIRST). Gate: the SymPy oracle (§3.1) + the
  dense numeric (§3.2). **Config: NON-UNIFORM h** (so `M=diag(h,θh)` varies
  cell-to-cell and `M⁻¹Aᵀ ≠ AᵀM⁻¹` is observable). A uniform-h stack is a
  weaker (not fully blind, θ≠1 always, but less-activating) config — §9 requires
  non-uniform h. Mutation **M-R1b-MASSORDER** (`AᵀM⁻¹ → M⁻¹Aᵀ`) reds both.
- **(b) face-coupling VJP.** `assemble_inflow_axis` transposes into the upstream
  face cotangent; the outgoing-face cotangent scatters back. Gate: the dense-Mᵀ
  and (for the 1-D LD walk) the assembled-Mᵀ (§4). Mutation **M-R1b-FACECOUPLE**
  (transpose the wrong trace `B(−1)` vs `B(+1)`, or scatter to the wrong face)
  reds the dense/assembled-Mᵀ.
- **(c) `octant_moment_frame_signs` self-transpose — THE likeliest sign-error
  site (task requirement #2).** The involution `s` (`_ubld.py:96-148`) is a
  diagonal ±1 (`s[o₀,…] = ∏ octant_signs[a]^{o_a}`). Forward realizes `y = D_s·
  kernel(D_s·x)`; the VJP is `x_bar = D_s·kernelᵀ(D_s·y_bar)` (same `D_s` both
  sides). A bug applies `D_s` on ONE side / forgets it / wrong sign. **This is
  OUTSIDE the "no-backward-octant" invariance group** (see §7): on a config
  where every octant sweeps `+x` (all `octant_signs ≥ 0`), `D_s ≡ I` and the
  mutation is a no-op — DESIGNED-GREEN. So the activating config MUST have
  backward-sweeping octants (any full quadrature does) AND anisotropic
  (non-zero slope moments ψ̂) so the sign lands on a live DOF. Gate: the SymPy
  oracle (§3.1, symbolic — sign-exact) + the dense-Mᵀ (backward-octant
  anisotropic). Mutation **M-R1b-FRAMESIGN** reds them AND stays GREEN on a
  single-direction/isotropic control (the asymmetry IS the config-blindness
  proof, L11).

### 3.4 Disposition

All of §3 is NEW (`ld_ubld.py` derive-functions + `test_ld_ubld_transpose_symbolic.py`
+ the dense-Mᵀ rows that extend `test_loss_transpose_solve.py` — see §8). The
kernel-pair registration itself (`LinearDiscontinuous` overrides the transpose
kernel) lands with R1c (flip-safety §10).

---

## §4. R1c — the LD 1-D reverse walk (moment-tailed cotangents; the slope-fold transposes)

The reverse-solve slab arm (`__init__.py:3853-3912`) is already w-generic (reads
`(a, inverse_denom, w)` off the same `CollisionCache`), so LD's flat face chain
flows through unchanged. **NEW = the slope-fold transposes** (`scan_slope_face_source`ᵀ,
`scan_reconstruct`ᵀ, both riding `_slope_fold`, `_ubld.py:382-405`) + the
moment-tailed `(...,2^d)` cotangent buffers. This flips
`LinearDiscontinuous.has_transpose_kernel → True`.

### 4.1 The TWO-guard lift (flip-safety — the "predicate lie" this file kills)

Two guards spell the LD deferral and MUST lift together onto the trait:
- the trait guard (`__init__.py:2957-2973`, `if not scheme.has_transpose_kernel:
  raise "no transpose kernel"`),
- the moment-count probe (`~3836-3844`, the moment-tail width check).

**Gate M-R1c-HALFLIFT:** lift ONLY ONE guard → a raising path remains on the LD
moment-tailed field → either the LD positive control reds (raises where it
should now compute) OR a moment-tailed apply crashes on shape. The
`test_ld_adjoint_deferral.py` LD rows (`test_ld_declares_false`,
`test_ld_streaming_is_not_adjointable`, `test_ld_composites_propagate_the_refusal`,
`test_ld_H_raises_missing_adjoint_at_construction`,
`test_ld_direct_transpose_raises_typed_deferral`,
`test_ld_bare_streaming_transpose_raises_too`) **flip from deferral-negatives to
positive controls IN THE SAME COMMIT as the flag flip** (L15/L19 — a flag
flipped before its walk is the predicate lie the file was minted to kill; a
half-lift is caught because the deferral file's LD rows now assert the POSITIVE
surface and one of the two guards still raising reds them).

### 4.2 The reciprocity moment-metric — the L18 lesson applied to LD (load-bearing, task requirement on `test_g_adjoint_reciprocity`)

**This is the sharpest verification-design point in R1c.** The task states "for
LD the independent bulk metric must pick up the moment-mass factor — spec this
precisely." Here it is:

The LD bulk field carries a trailing moment axis `{1, P₁}` of size 2. The LD
mass is `M = diag(h, θh)` (`_ubld.py:151`, θ=1/3). A moment DOF's Hilbert STATE
metric is set by its OPERATOR ROLE (L18, SHARPENED): the slope moment ψ̂ enters
`φ̂ = Σ_n w_n ψ̂_n` and its self-block in `A` carries the `θh` mass, so the
correct bulk metric is

```
G_bulk_LD = (V·w_n) ⊗ diag(1, θ)          # average moment: V·w_n ; slope moment: θ·V·w_n
```

**The Mode-12 trap (must be flagged in the gate):** if the reciprocity metric
uses average-only `V·w_n` (drops the `θ` on the slope), the slope DOFs are
mis-weighted and one of two failures occurs:
1. If the production `.H` metric ALSO drops the slope mass, then `A.H = G⁻¹AᵀG`
   with a degenerate/wrong G is a WRONG adjoint (the L18 ghost-G_sd=0 family) —
   and reciprocity with the matching wrong G is internally consistent but
   BLIND to any slope-row transpose error (the error class sits in the
   stabiliser of the wrong metric). Green, and blind.
2. If the production `.H` uses the correct `θ`-metric but the test's `_g_inner`
   drops it, reciprocity FALSE-REDS even unmutated.

**Required design (both legs, L18):**
- Extend `_g_inner`'s `_bulk_measure` to carry the moment mass: `V·w_n ⊗
  diag(1,θ)` when the interior has `spatial_moments==2`. Verify it against the
  production metric via the EXISTING cross-check
  `test_full_field_space_metric_matches_independent_reference` extended to an
  LD builder — this pins that BOTH sides carry the moment mass (catches
  failure 2, and is the evidence the production metric is NOT the degenerate
  average-only one).
- **Control leg:** the unmutated LD-slab reciprocity holds `< 1e-12`.
- **Mutated-reds leg (M-R1c-SLOPEFOLD):** a sign flip in the slope-fold
  transpose moves the reciprocity residual `> 1e-6` — BUT this tooth's teeth
  DEPEND on the moment mass being in the metric. Demonstrate the dependency:
  with the average-only metric the SAME slope-flip mutation goes GREEN (blind);
  with the `θ`-metric it reds. That asymmetry IS the proof the moment mass is
  load-bearing (the L18 "encode the metric as the catcher" discipline; NOT a
  slab-vs-curvilinear split — the θ-metric is load-bearing on the LD SLOPE DOF).
- **Belt-and-braces:** the LOAD-BEARING slope-row transpose VALUE gate is the
  metric-FREE Euclidean oracle — the SymPy `ld_ubld` transpose-oracle (§3.1) +
  the dense/assembled-`Mᵀ` (§4.3). These catch M-R1c-SLOPEFOLD regardless of
  the metric, so the reciprocity moment-metric is the COMPOSITE-integration
  canary, not the sole catcher. (Lean the correctness on the Euclidean oracle;
  use reciprocity to prove the metric — and hence `.H` — is right.)

### 4.3 The assembled-Mᵀ for LD-slab (Cartesian — a second structurally-independent oracle)

Assembly (`assemble_ubld` → sparse) is Cartesian and handles LD via
`n_face_moments`. So the LD-slab reverse gets the LAPACK oracle:
`solve_triangular(permutedᵀ, b[order], lower=False)` where `permuted` = the
LD-slab assembled block in walk order (the transposed permuted block is
upper-triangular). Structurally independent (LAPACK shares no code with the LD
scan). Catches a wrong transposed-scan-coefficient. Triangularity leg:
`triu(permutedᵀ, k=1)` structure per the LD moment bandwidth. **This resolves
the L17 "LD caveat" (the a3 memo flagged the LD assembled-Mᵀ as possibly moot
if the 1-D adjoint was DD-only) — it is NOT moot; the LD moment tail lands here,
so the LD-slab assembled-Mᵀ oracle IS available.**

### 4.4 Gates for R1c (extend `test_loss_transpose_solve.py` + `test_g_adjoint_reciprocity.py`)

| Gate | Disposition | Config | Reference (independence) | Tol |
|---|---|---|---|---|
| G1 round-trip (LD-slab) | EXTEND `_MESHES` | LD-slab, het ≥2G, non-uniform h | necessary-not-sufficient (shares walk) | `rtol=1e-10` |
| G2 dense-`Mᵀ` (LD-slab) | EXTEND | LD-slab | `np.linalg.solve(Mᵀ,b)`, M off FORWARD apply | `rtol=1e-10` |
| assembled-`Mᵀ` (LD-slab) | NEW row | LD-slab | LAPACK `solve_triangular(permutedᵀ)` | `rtol` |
| SymPy transpose-oracle | NEW file | symbolic | `simplify(VJP−AᵀM⁻¹)==0` | exact |
| G-reciprocity (LD-slab) | EXTEND `_BUILDERS` + moment-metric | LD-slab het ≥2G aniso | `_g_inner` w/ moment mass = production `.H` | `rel<1e-12` |
| moment-metric cross-check | EXTEND `test_full_field_space_metric_matches_independent_reference` | LD-slab | production `space.inner_product` | `rel<1e-13` |
| `test_ld_adjoint_deferral` LD rows | REWRITE negative→positive (same commit) | LD-slab | — | raises/holds |

**Anisotropy audit (§0.6 / L18 war-story):** an all-isotropic LD suite is BLIND
to a dropped/mis-signed slope moment. The LD reciprocity + dense-Mᵀ configs MUST
carry non-zero slope moments — build the LD-slab reciprocity input with a
NON-FLAT bulk (the existing `_random_composite` gives random ψ̂) AND the σ_t /
SigS asymmetric so the slope is genuinely exercised. Manufacture the anisotropic
2-term `P_ℓ(−1)` check FIRST if any q½/source fold is on the reverse path.

---

## §5. R2a / R2b — the multi-D Cartesian DD reverse walk (matvec `apply_transpose`)

R2 delivers multi-D `StreamingOperator.apply_transpose = Lᵀ` (the matvec
transpose — the reverse of `_CellResidual`, NOT a transpose-solve; the reverse
sweep / `solve_transpose` for multi-D is the OUT-of-scope G-S deferral, §0.2).
The reverse = reversed levels + a THIRD level-op `_CellResidualTranspose`
(kernel-pair VJP) + transposed addressing (gather at `face_out`,
scatter-accumulate at `face_in`, domain-boundary in↔out swap) + reversed
frontier. **Claim layer:** operator-identity (Euclidean `Lᵀ` + metric `L.H`) —
NEVER eigenvalue (§0.4).

### 5.1 The reverse full≡window oracle discipline (task requirement #4 — mirror the forward's)

The forward pins `walk_windowed ≡ walk_full` bit-identically (the `window ≡ full`
oracles, d=1/2/3 — same cell math, same level order, different storage). The
reverse MUST mirror it EXACTLY:

- **R2a lands FIRST: the reverse `walk_full` (full-cochain) is the ORACLE.** It
  is the fuller-view reference (every face retained) — the structurally-simplest
  reverse to get right. Its correctness is pinned by the 2-D dense-`Mᵀ` +
  assembled-`Mᵀ` + reciprocity (§5.3), NOT by the windowed path.
- **R2b lands the reverse `walk_windowed` (rolling-frontier) — PRODUCTION —
  pinned bit-identical to the reverse `walk_full` oracle.** Gate
  **`test_reverse_window_equals_full`**: reverse-`walk_windowed(φ)` ≡
  reverse-`walk_full(φ)` at `np.array_equal` (same cell math, same reversed
  level order, different storage — bit-identical is the RIGHT contract, L16
  window≡full sibling), on het + asymmetric SigS + **rectangular nx≠ny** +
  (for R2c) anisotropic. Mutation **M-R2-WINDOWDRIFT** (perturb the frontier
  seed/shed order in the reverse windowed) reds it.
- The predicate `has_transpose_walk` flips True ONLY when the PRODUCTION
  (windowed) reverse works (`apply_transpose` routes through production) — so
  the flip lands in the R2b commit, NOT R2a (R2a's oracle is green-on-contact
  but the predicate stays False; the deferral file's multi-D rows stay negative
  until R2b — flip-safety §10).

### 5.2 The reverse walk spy — the multi-D sibling of `test_one_dim_loop_walk`

The 1-D shared loop frame has `test_one_dim_loop_walk.py` (runtime spy both
orientations hit `_loop_walk` + AST tripwire banning orientation booleans). The
multi-D reverse CREATES a new shared frame (the reversed level walk + the
`_CellResidualTranspose` object) → it needs the SIBLING (L17):

- **NEW `tests/sn/sweep/core/test_multi_d_reverse_walk.py`:**
  1. **Runtime spy (Mode-11 wrap):** the reverse `apply_transpose` on cart2d
     routes through the reversed `walk_full`/`walk_windowed` frame with the
     `_CellResidualTranspose` level-op; counter `> 0`. Direction is an OBJECT
     (`_CellSolve` / `_CellResidual` / `_CellResidualTranspose`), never a flag.
  2. **AST tripwire:** forbid `is_adjoint` / `is_forward` / `is_transpose` /
     `is_reverse` as real identifiers (Name/Attribute/arg) in the shared walk
     frame; demand the orientation-carrying level-op object. Verified clean
     today (L17 — none of those identifiers exist). Mutation **M-R2-SPY**
     (route the reverse via an `is_reverse` boolean) reds the tripwire.

### 5.3 The 2-D structurally-independent oracles

- **Dense-Euclidean-`Mᵀ` (matvec form — NEW 2-D artifact).** The existing
  `_probe_augmented_matrix_one_group` (imported into `test_loss_transpose_solve`
  from `test_assembly_mode`) is **1-D-only** → a **2-D forward-apply column-probe
  is a NEW oracle artifact** (task requirement #1). Build `M` by probing the
  FORWARD `apply` on the 2-D one-group unit-vector basis (bulk DOFs,
  ordinate-blocked); then `apply_transpose(x) ≡ M.T @ x` per element
  (`assert_allclose(atol=1e-11)`; matvec-transpose is exact `Mᵀ@x` up to the
  probe reduction — nulp/rtol, NOT `array_equal`, L16 sparse-order). This is the
  stabiliser-escape: it pins the OBJECT (`Mᵀ`), outside every spectral
  invariance group.
- **Assembled-`Mᵀ` (matvec, Cartesian, 2-D).** `SparseAssembledOperator.apply_transpose`
  = exact CSR `M.T @ x`; `M` = the 2-D assembled `(L+C)` block. Structurally
  independent (assembly from forward-kernel unit-probes; the assembly mode is
  Cartesian and handles 2-D). Catches the transposed-addressing bug the dense
  probe also catches, from a second angle.
- **G-adjoint reciprocity (cart2d-DD — the DEFINING property, the composite
  canary).** EXTEND `test_g_adjoint_reciprocity._BUILDERS` / `_FULL_LOSS_BUILDERS`
  with cart2d-DD rows. `⟨Lψ,φ⟩_G = ⟨ψ,L.Hφ⟩_G` on the FULL composite (bulk⊕trace
  random, `_random_composite`), `rel<1e-12`. **Full composite is mandatory, NOT
  bulk-only:** the domain-boundary in↔out swap lives in the trace metric coupling
  — a bulk-only reciprocity is blind to it (the a4/cyl G3 lesson: bulk-only
  fails to catch a boundary-cotangent bug that full-field catches at O(1)). The
  metric `G = V·w_n (bulk) ⊕ |Ω·n|·w_n (trace)` is non-trivial on cart2d (the
  L19 correction: `.H` ≠ Euclidean on EVERY geometry, incl. Cartesian — the
  trace weight `|Ω·n| ≠ V`; do NOT design a slab/Cartesian-green metric split).

### 5.4 The frozen baseline for the multi-D reverse (NEW value → re-baseline, NOT pre-carve bit-id)

The multi-D adjoint DOES NOT EXIST pre-carve (it raises) → there is NO pre-carve
value to freeze (contrast R1a). Per L17/2.5b (NEW values → re-baseline, NOT
bit-id): capture a cart2d-DD reverse baseline (fwd already exists + the NEW adj)
POST-carve, VERIFIED against the §5.3 oracles at capture time, then frozen as a
permanent regression canary. Add a `cart2d_2g` (nonsquare) CASE to
`_generate_walk_baselines.py`. The forward cart2d matvec (`residual_kernel_batch`
via `walk_full`/`walk_windowed`) MUST stay bit-identical (existing 2-D octant
snapshots + `test_cell_kernel_batch`) — R2 ADDS the transpose kernel, it does
NOT touch the forward.

### 5.5 Config (task requirement #5 — the Mode-2 axis-swap catcher)

**Rectangular nx≠ny is MANDATORY** (`cart2d_2g_nonsquare()` already exists in
`tests.sn._test_helpers`). A square nx=ny + uniform σ + symmetric SigS config is
transpose-BLIND to an x↔y axis swap AND to a DOF-transposition (L16 Mode-12:
row/col swap on `A=Aᵀ` is invisible). Require: rectangular nx≠ny + asymmetric
SigS + non-uniform h so `A≠Aᵀ` and the axis-swap/addressing bug is observable.
≥2G always (anti-#3). Mutation **M-R2-AXISSWAP** (swap x↔y in the transposed
addressing) reds the dense-`Mᵀ` on rectangular-asymmetric AND stays GREEN on a
square-uniform control (the asymmetry IS the config-blindness proof, L16).

### 5.6 Existing-gate un-exclusions (same commit as the R2b flip)

- `test_removal_form_matvec_sweep::test_invertible_apply_transpose…bit_identical`
  — UN-EXCLUDE cart2d (currently excluded "deferral raise"). Self-referential
  (override-not-leak — necessary-not-sufficient); the genuine proof is §5.3.
- `test_phase_c_gates::test_apply_apply_transpose_reciprocity_under_sweep_frame`
  (:407) — ADD the cart2d parametrization row (Euclidean composite reciprocity).
- `test_ld_adjoint_deferral::TestMultiDOrientationHonesty` — the cart2d-DD rows
  (`test_cart2d_dd_streaming_is_not_adjointable`,
  `test_cart2d_H_raises_missing_adjoint_at_construction`,
  `test_cart2d_direct_transpose_still_raises_typed_deferral`) FLIP from negative
  to positive controls in the R2b commit; `test_rep_trait_declarations`
  `ScanMarch(sn2).has_transpose_walk` / `MovingFrontierWindow(sn2)` flip
  `False→True` (flip-safety §10).

---

## §6. R2c — the LD-2D reverse (`n_face_moments=2`)

R2c rides the R2b multi-D reverse frame with the LD kernel (R1b) at
`n_face_moments = (per_axis)^{d-1} = 2` (the frontier slabs carry the trailing
transverse moment axis). The NEW risk over DD-2D is the transverse moment
coupling + the d=2 cross-moment `x̂y` frame sign.

### 6.1 Gates

- **assembled-`Mᵀ` LD-2D (the keystone).** `assemble_ubld` handles LD-2D → the
  LAPACK `solve_triangular(permutedᵀ)` / CSR `Mᵀ@x` oracle IS available and is
  the strongest structurally-independent reference (LD-2D has no dense-block
  matvec of its own that the probe could tautologically reproduce; the assembled
  block is emitted from the shared UBLD source). `rtol` (sparse order).
- **dense-`Mᵀ` LD-2D** (the 2-D probe from §5.3, moment-tailed).
- **reverse full≡window (LD-2D)** — `n_face_moments=2` frontier ≡ full-cochain,
  `array_equal`, on het + asymmetric + rectangular + **anisotropic**.
- **reciprocity (LD-2D)** with the moment-mass metric (§4.2 `G_bulk_LD` extended
  to 2-D: `V·w_n ⊗ (diag(1,θ) ⊗ diag(1,θ))` for the `2^d` moment Kronecker
  layout).

### 6.2 The d=2 cross-moment frame sign (Mode-12, the §3.3(c) generalization)

The 2-D cross moment `x̂y` flips when an ODD number of its active axes sweep
backward (`_ubld.py:120-125`). **Config: an octant sweeping backward on EXACTLY
ONE axis** (so the cross-moment sign is `−1`) AND anisotropic. Mutation
**M-R2c-FRAMESIGN-2D** (drop the cross-moment sign / use `∏` wrong) reds the
assembled-`Mᵀ` LD-2D on the one-axis-backward anisotropic octant; GREEN on a
both-forward or both-backward control (even sign) — the asymmetry proof.

### 6.3 The moment-drop config-blindness (Mode-7 — MANDATORY anisotropic)

**M-R2c-MOMENTDROP:** `n_face_moments` reverts to 1 (drops the transverse slope)
→ reds assembled-`Mᵀ` + dense-`Mᵀ` + full≡window. **This is GREEN on an
isotropic LD-2D config** (§0.6 / L1: an all-isotropic 2-D snapshot is blind to a
dropped φ_ℓ≥1). So every R2c gate MUST run an ANISOTROPIC input (non-zero slope
moments) — reserve isotropic ONLY as the config-asymmetry control that shows the
mutation goes green there. Before carving R2c, AUDIT that the existing 2-D LD
snapshots exercise the moments being reduced; if all isotropic, manufacture the
anisotropic case FIRST (L1 the moment-reduction discipline).

---

## §7. Mode-12 invariance-group analysis (task requirement #2 — done at DESIGN time)

Before any mutation is run, enumerate each measured functional's invariance
group (stabiliser) and intersect it with the carve's mutation classes. A
mutation INSIDE the stabiliser is DESIGNED-GREEN — no tolerance/refinement/regime
change can expose it through that functional; the remedy is a gate on a
functional OUTSIDE the stabiliser (canonically the OBJECT itself).

### 7.1 The measured functionals and their stabilisers

| Functional | Invariance group (blind to) | Verdict for this carve |
|---|---|---|
| **eigenvalue / keff / spectral norm** | similarity conjugation + **transpose** (`eig(Aᵀ)=eig(A)`) + adjoint | **DESIGNED-GREEN on the ENTIRE carve.** BANNED as a gate (§0.4). |
| **G-reciprocity `⟨Aψ,φ⟩_G=⟨ψ,A†φ⟩_G`** | (1) a transpose error `E` with `G·E=0` — i.e. `E` supported on the **metric null space** (zero-G-weight DOFs); (2) an `E` where the test's G ≠ production's G masks it | Catches any `E` on metric-nonzero DOFs. **The LD slope DOF is metric-nonzero ONLY if the metric carries the `θ` moment mass** (§4.2) — else the slope-row transpose is in the stabiliser (blind). The trace in↔out swap is caught ONLY by a full-composite (bulk⊕trace) ψ,φ — bulk-only reciprocity has the trace error in its stabiliser (§5.3). |
| **`octant_moment_frame_signs` involution error** | configs with NO backward-sweeping octant (`D_s≡I`) OR isotropic (ψ̂=0) | Caught by SymPy oracle + dense-`Mᵀ` on backward-octant anisotropic (§3.3c). Blind on single-direction/isotropic — hence those are the GREEN control, not the gate. |
| **G1 round-trip `solve_transpose∘apply_transpose=I`** | a shared error `E` common to both arms (the walk/kernel both ride it) — cancels | Necessary-not-sufficient; PAIR with the Euclidean oracle (which does not share the reverse code). |
| **dense-Euclidean-`Mᵀ` (`apply_transpose≡Mᵀx` / `solve_transpose≡solve(Mᵀ,b)`)** | none relevant — pins the OBJECT (the matrix) per element | **The stabiliser-escape / keystone.** Outside every spectral group. M from FORWARD apply → independent of the reverse code. |
| **assembled-`Mᵀ` (LAPACK/CSR)** | Cartesian-only (curvilinear invisible) | Object-level; the transposed-scan-coefficient catcher. Curvilinear escapes → the dense-`Mᵀ` keystone covers sphere/cyl. |
| **SymPy `simplify(VJP−AᵀM⁻¹)==0`** | none — symbolic identity | The algebra-of-record keystone (sign-exact; catches the frame-sign + mass-order). |

### 7.2 The three carve-specific Mode-12 hazards (the ones that WILL bite)

1. **The LD moment-metric hazard (§4.2).** If the `.H` metric drops the slope
   mass `θ`, the slope-row transpose error sits in reciprocity's stabiliser
   → reciprocity is blind AND `.H` is a wrong adjoint. **Closure:** repair the
   metric (put `θ` in `G_bulk_LD`), which BOTH fixes `.H` AND makes reciprocity
   a real catcher (the ERR-067 "repair the metric, not the gate" mechanism —
   the correctness fix and the Mode-12 closure are one and the same). Gate with
   BOTH legs (control-holds + mutated-reds) AND exercise the previously-nulled
   input (nonzero slope) — a still-broken average-only baseline mimics "caught"
   otherwise (the L18 false-closure trap).
2. **The frame-sign no-backward-octant hazard (§3.3c, §6.2).** The involution
   error is identically zero when `D_s≡I`. Never let a frame-sign gate run only
   a single-direction config — it is DESIGNED-GREEN there.
3. **The bulk-only reciprocity hazard (§5.3).** The domain-boundary in↔out swap
   is in a bulk-only reciprocity's stabiliser. Full composite (bulk⊕trace)
   mandatory.

### 7.3 The escape for all three: pin the OBJECT

Every one of the three hazards is closed by the SAME move — the object-level
oracle (dense-`Mᵀ`, assembled-`Mᵀ`, SymPy) pins the matrix/VJP itself, outside
every spectral/metric stabiliser. **The reciprocity gates are the composite
INTEGRATION canaries (they prove `.H` = `G⁻¹AᵀG` with the right G composes);
the object oracles are the CORRECTNESS keystones.** This is the load-bearing
split of the whole spec: lean correctness on the Euclidean/symbolic object
oracles; use reciprocity to prove the metric wrap.

---

## §8. Existing-gate disposition (EXTEND vs NEW — task requirement #1)

Verified against the code (2026-07-23). "EXTEND" = add builders/rows to a live
file; "REWRITE" = flip negatives→positives in the same commit as a flag flip;
"NEW" = a new artifact.

| File / gate | Disposition | Action |
|---|---|---|
| `test_ld_adjoint_deferral.py` | **REWRITE (R1c + R2b commits)** | LD rows → positive controls (R1c); `TestMultiDOrientationHonesty` cart2d rows + `test_rep_trait_declarations` flags → positive (R2b). The base-default-False row (`test_base_default_is_false`) STAYS (a scheme is still opt-in). **This file is the flip-safety tripwire — it MUST be rewritten in the SAME commit as each flag flip** (a flag flipped before its walk = the predicate lie the file kills). |
| `test_g_adjoint_reciprocity.py` | **EXTEND** | `_BUILDERS` += LD-slab, cart2d-DD; `_FULL_LOSS_BUILDERS` += cart2d-DD 2g. `_bulk_measure` → moment-mass `θ` for LD (§4.2). `test_full_field_space_metric_matches_independent_reference` += LD builder (pins the moment metric). New tooth `test_tooth_frame_sign_flip_reds` (LD, backward-octant aniso). |
| `test_loss_transpose_solve.py` | **EXTEND + NEW oracle** | `_MESHES` += LD-slab (G1/G2 reduce to the slab bulk-only branch). Assembled-`Mᵀ` LD-slab row. **NEW 2-D forward-apply column-probe** (`_probe_augmented_matrix_2d_one_group`) — the 1-D `_probe_augmented_matrix_one_group` is 1-D-only; cart2d + LD-2D dense-`Mᵀ` need it. |
| `test_inverse_adjoint_coherence.py` | **EXTEND (only if multi-D solve_transpose lands — it does NOT here)** | R2 is matvec-`apply_transpose`, not `solve_transpose`; multi-D `solve_transpose` is the OUT-of-scope G-S deferral (§0.2). Leave this file's meshes {slab, cyl_product} unchanged. Add a NOTE row asserting cart2d `(L+C).inverse().is_adjointable` stays consistent with the multi-D solve deferral (the swap law does not extend to cart2d in this carve). |
| `test_removal_form_matvec_sweep::…apply_transpose…bit_identical` | **EXTEND (R2b commit)** | UN-EXCLUDE cart2d. Self-referential (necessary-not-sufficient). |
| `test_phase_c_gates::…reciprocity_under_sweep_frame` (:407) | **EXTEND (R2b commit)** | add the cart2d parametrization row. |
| `test_one_dim_loop_walk.py` | **keep + NEW sibling** | unchanged; **NEW `test_multi_d_reverse_walk.py`** is its multi-D reverse sibling (§5.2 — spy + AST tripwire). |
| `walk_matvec_*.npz` + `_generate_walk_baselines.py` | **EXTEND (R2b)** | += `cart2d_2g` (nonsquare) CASE — POST-carve re-baseline (new adjoint value verified against §5.3 oracles), NOT a pre-carve bit-id freeze. slab/sphere/cyl stay as the R1a bit-id anchor. |
| `test_walk_matvec_baselines.py` | stays green (R1a) | the R1a relocation must not move the frozen values. |
| `ld_ubld.py` (SymPy) + `test_ld_ubld_transpose_symbolic.py` | **NEW** | the transpose algebra-of-record (§3.1). |
| `test_ld_ubld_primitive.py` (`_ubld` numeric) | **EXTEND** | += the dense numeric transpose oracle (§3.2). |

### 8.1 The `catches("ERR-066")` audit (do NOT inflate coverage)

`test_g_adjoint_reciprocity_full_block` carries `catches("ERR-066")` (the
degenerate-ordinate cyl-product drop). The cart2d/LD rows added to `_BUILDERS`
do NOT catch ERR-066 (different bug) — do NOT let the marker leak onto them.
Per vv-principles, a `catches` marker is a COVERAGE CLAIM, mutation-verified per
row; keep ERR-066 on the cyl_product row only. If R1b/R1c/R2 surface a NEW
caught bug (a real transpose defect found during the carve), file an ERR-NNN and
attach `catches` to the SPECIFIC catching row (mutation-verified), never the
family.

---

## §9. Heterogeneity / asymmetry config requirements (task requirement #5)

Binding config contract — a shipped gate that swaps in a convenient config is a
gate downgrade (L3), caught by qa. Per anti-#3 / anti-#4 / Mode-2:

| Requirement | Why | Applies to |
|---|---|---|
| **≥2G always; NO 1G-only** | k=νΣ_f/Σ_a is flux-shape-independent — but more to the point, a 1G group-axis nulls the group-transfer transpose; explicitly exclude 1G, ship NO 1G "structural smoke" | every value gate |
| **Asymmetric SigS** (strong downscatter, weak upscatter) | symmetric SigS ⟹ S†=S ⟹ the group-transfer transpose is invisible (ERR-002 family) | reciprocity full-loss rows, any gate touching S |
| **Non-uniform h** | uniform h nulls the per-cell mass variation → M-R1b-MASSORDER (`AᵀM⁻¹` vs `M⁻¹Aᵀ`) less-activating; and a uniform mesh masks a `dr[k]` index drift (Mode-5) | LD kernel gates, all dense-`Mᵀ` |
| **Rectangular nx≠ny** (2-D) | square nx=ny is transpose-BLIND to an x↔y axis swap + DOF-transposition on symmetric A (Mode-2/L16) — `cart2d_2g_nonsquare()` exists | all R2/R2c gates |
| **Anisotropic (non-zero ψ̂)** | isotropic nulls the LD slope moment + the frame-sign + a dropped φ_ℓ≥1 (§0.6/L1/Mode-7) | all LD gates |
| **Backward-sweeping octant present** | `D_s≡I` on single-direction → frame-sign mutation is a no-op (§7.2) | frame-sign gates |
| **Full composite (bulk⊕trace) random** | bulk-only reciprocity is blind to the domain-boundary in↔out swap (§5.3) | all reciprocity |
| **NON-FLAT random input** (fixed seed) | flat ψ nulls the DD streaming coupling AND the boundary swap (§0.6/H2) | all matvec/round-trip gates |

**The config-asymmetry CONTROL is itself a gate.** For every Mode-12 hazard, the
gate ships BOTH the activating config (mutation reds) AND the nulling control
(mutation goes GREEN) — the asymmetry positively EXHIBITS the config-blindness
(the L11/L1 discipline; e.g. M-R1b-FRAMESIGN reds backward-octant-aniso, greens
single-direction; M-R2-AXISSWAP reds rectangular-asymmetric, greens
square-uniform).

---

## §10. Same-commit sets / land-order (task requirement #6 — the flip-safety law)

**The flip-safety law (L15/L19):** a flag flip is construction-time, a walk
raise is apply-time. A flag flipped before its kernel/walk exists is the
"predicate lie" `test_ld_adjoint_deferral.py` was minted to kill. Therefore each
`has_transpose_kernel` / `has_transpose_walk` flip MUST land in the SAME commit
as (a) its kernel/walk implementation AND (b) the rewrite of the deferral file's
corresponding negatives → positive controls AND (c) the extend-gates that
`pytest.fail` when `not A.is_adjointable` (they start RUNNING their body on that
config the instant the flag flips).

| Commit | Lands | Flag flip | Same-commit gate set |
|---|---|---|---|
| **C1 — R1a** | DD transpose relocated into the registered kernel pair; angular thread stays on the closure (§1.2) | NONE (DD already True; registration-coupled) | frozen `walk_matvec_{slab,sphere,cyl}` stay green (`-W error::DriftWarning`) + NEW Mode-11 sentinel (kernel + closure) + retirement grep + M-R1a-DEAD/VALUE teeth |
| **C2 — R1b+R1c** | LD transpose kernel (`ld_ubld` VJP) + LD 1-D reverse walk + BOTH LD guards lifted onto the trait | `LinearDiscontinuous.has_transpose_kernel → True` | SymPy `test_ld_ubld_transpose_symbolic` (NEW) + dense numeric oracle + G1/G2/assembled-`Mᵀ` LD-slab (EXTEND `_MESHES`) + reciprocity LD-slab + moment-metric cross-check (EXTEND) + **REWRITE `test_ld_adjoint_deferral` LD rows → positive** + M-R1b-{MASSORDER,FACECOUPLE,FRAMESIGN} + M-R1c-{SLOPEFOLD,HALFLIFT,METRICBLIND} teeth |
| **C3 — R2a** | multi-D reverse `walk_full` (full-cochain ORACLE) + `_CellResidualTranspose` level-op + transposed addressing | NONE (oracle; predicate stays False — the reverse `walk_full` is not the production `apply_transpose` path yet) | NEW 2-D dense-`Mᵀ` probe + assembled-`Mᵀ` 2-D + reciprocity cart2d-DD (EXTEND, green-on-contact against the oracle) + M-R2-{ADDRESSING,AXISSWAP,LEVELORDER} teeth on the oracle |
| **C4 — R2b** | multi-D reverse `walk_windowed` (PRODUCTION) + `ordinate_scan_transpose` transverse chaining + wire `apply_transpose` to production | `ScanMarch(2d).has_transpose_walk → True`, `MovingFrontierWindow.has_transpose_walk → True` | `test_reverse_window_equals_full` (full≡window) + NEW `test_multi_d_reverse_walk` (spy + AST tripwire) + **REWRITE `test_ld_adjoint_deferral::TestMultiDOrientationHonesty` cart2d + flag rows → positive** + UN-EXCLUDE `test_removal_form` cart2d + `test_phase_c_gates` cart2d row + `cart2d_2g` baseline CASE + M-R2-{WINDOWDRIFT,SPY} teeth |
| **C5 — R2c** | LD-2D reverse (`n_face_moments=2`) | (rides C4 flip; LD-2D kernel from C2) | assembled-`Mᵀ` LD-2D (keystone) + dense-`Mᵀ` LD-2D + full≡window LD-2D (aniso) + reciprocity LD-2D (moment metric ⊗) + M-R2c-{MOMENTDROP,FRAMESIGN-2D} teeth |

**Ordering rationale:** C1 is a pure bit-id refactor (safest first — no
predicate change). C2 completes the 1-D adjoint capability (LD). C3 lands the
reverse ORACLE before the production it pins (the reverse full≡window discipline,
mirror of the forward). C4 flips the multi-D predicate ONLY when production
(windowed) works, pinned by C3's oracle. C5 rides C4's frame with C2's kernel.
Each Cᵢ is independently green under `python -O -m pytest`; a Cᵢ that leaves the
deferral file un-rewritten is a predicate lie (reject at review).

---

## §11. Mutation → gate teeth matrix (task requirement #1 + Mode-10)

Every NEW gate class is mutation-verified: name the mutation that reds it,
confirm it fires under `python -O`, revert by **in-process monkeypatch, NEVER
`git checkout`** (uncommitted surgical file — process-discipline.md). All value
gates use `np.testing` / `pytest.fail` (Mode-8 — bare `assert` in `orpheus/`
strips under `-O`; in collected `tests/` it is rewritten and fires).

| # | Mutation | Reds | Config-asymmetry control (goes GREEN) |
|---|---|---|---|
| **M-R1a-DEAD** | registered transpose kernel stubbed to zeros | Mode-11 sentinel (count>0 fails) | — (proves the sentinel catches what the baseline can't) |
| **M-R1a-VALUE** | `2·f_bar → 1.9·f_bar` in the relocated DD VJP | frozen baseline (nulp) + reciprocity + DriftWarning escalates | — |
| **M-R1b-MASSORDER** | `AᵀM⁻¹ → M⁻¹Aᵀ` (mass-inverse after Aᵀ) | SymPy oracle + dense-`Mᵀ` (non-uniform h) | uniform-h less-activating (require non-uniform) |
| **M-R1b-FACECOUPLE** | transpose wrong face trace (`B(−1)`↔`B(+1)`) / scatter to wrong face | dense-`Mᵀ` + assembled-`Mᵀ` (LD-slab) | — |
| **M-R1b-FRAMESIGN** | drop `octant_moment_frame_signs` on ONE VJP side | SymPy oracle + dense-`Mᵀ` (backward-octant aniso) | **single-direction / isotropic → GREEN** (the §7.2 stabiliser) |
| **M-R1c-SLOPEFOLD** | sign flip in `scan_slope_face_source`ᵀ / `scan_reconstruct`ᵀ | SymPy oracle + dense/assembled-`Mᵀ` (LD-slab) + reciprocity (θ-metric) | isotropic → GREEN |
| **M-R1c-HALFLIFT** | lift ONE of the two LD guards only | `test_ld_adjoint_deferral` LD positive control (raises where it should compute) | — |
| **M-R1c-METRICBLIND** | reciprocity uses average-only metric (drop `θ`) | control-leg FALSE-REDS unmutated OR slope-flip goes GREEN | demonstrates the θ-metric is load-bearing (L18 both-legs) |
| **M-R2-ADDRESSING** | gather at `face_in` (not `face_out`) / drop domain-boundary in↔out swap | 2-D dense-`Mᵀ` + full-composite reciprocity | — |
| **M-R2-AXISSWAP** | swap x↔y in transposed addressing | 2-D dense-`Mᵀ` (rectangular nx≠ny, asymmetric) | **square nx=ny + uniform → GREEN** (Mode-2/L16) |
| **M-R2-LEVELORDER** | don't reverse the level order | full≡window + dense-`Mᵀ` | — |
| **M-R2-WINDOWDRIFT** | perturb reverse frontier seed/shed order | `test_reverse_window_equals_full` (array_equal) | — |
| **M-R2-SPY** | route reverse via `is_reverse` boolean | AST tripwire + reverse spy | — |
| **M-R2c-MOMENTDROP** | `n_face_moments → 1` (drop transverse slope) | assembled-`Mᵀ` LD-2D + dense-`Mᵀ` + full≡window | **isotropic → GREEN** (Mode-7/§0.6) |
| **M-R2c-FRAMESIGN-2D** | drop the cross-moment `x̂y` sign (odd backward axes) | assembled-`Mᵀ` LD-2D (one-axis-backward aniso) | both-forward / both-backward (even sign) → GREEN |

**Cross-cutting teeth (the one-source proof, L16):** the LD kernel VJP is
consumed by BOTH the 1-D reverse walk (R1c) AND the 2-D reverse walk (R2c). A
sign flip in the SHARED `ld_ubld` VJP source must red BOTH the R1c AND R2c gates
(+ the SymPy oracle). If it reds only the new 2-D gate, a parallel stencil (twin
path) exists → STOP, fix, log ERR-NNN (the Cardinal-Rule-2 guardrail). DD is the
higher twin-path risk than LD (DD has no dense block; its matvec is the fused
scalar kernel — a fresh transpose stencil is tempting), so the M-R1a-VALUE /
M-R2-ADDRESSING shared-source mutations are the highest-value teeth.

---

## §12. Scaffold order + acceptance checklist

### 12.1 Per-commit scaffold (pre-carve canary vs post-carve proof)

For each Cᵢ (§10):
- **PRE:** confirm the relevant forward references are landed+green (regression
  baseline); for C1 the frozen `walk_matvec_*` exist; for C2 the LD forward
  primitive (`test_ld_ubld_primitive`) is green; for C3/C4 the forward 2-D
  octant snapshots + `test_cell_kernel_batch` are green.
- **POST (the proof):** the extend/new gates GREEN, THEN the §11 mutation matrix
  RED-verified (each under `-O`), with the config-asymmetry controls positively
  exhibiting the blindness (mutation greens on the nulling config).

### 12.2 Acceptance checklist (a gate is done only when every row is TRUE)

- [ ] **Every gate is object-level** — no keff/eigenvalue/norm gate credits any
      transpose (§0.4 / §7). Grep the new gates for `keff`/`k_inf`/`eigval` → ZERO.
- [ ] **Every NEW gate class has a RED mutation** under `python -O` (§11), reverted
      by monkeypatch (never `git checkout`).
- [ ] **Every Mode-12 hazard ships both legs** — control-holds AND mutated-reds,
      with the previously-nulled input exercised (§7.2 / §4.2).
- [ ] **The reference is structurally independent of the reverse code** — SymPy /
      dense-`Mᵀ` (M off FORWARD apply) / assembled-`Mᵀ` (LAPACK) / frozen pre-carve
      baseline. G1 round-trip is corroborating only (§0.3).
- [ ] **The regime activates the term** — non-uniform h (mass order), backward
      octant + anisotropic (frame sign), rectangular nx≠ny (axis swap), full
      composite (trace swap), non-flat input (streaming coupling) (§9).
- [ ] **Each flag flip is same-commit with its walk + the deferral rewrite + the
      is_adjointable-gated extend-gates** (§10 flip-safety).
- [ ] **The relocated kernel is proven on the path** (Mode-11 sentinel + retirement
      grep, C1) — a green baseline alone does NOT prove the relocation (§2.2).
- [ ] **The reverse full≡window pin is bit-identical** on het+asymmetric+
      rectangular+anisotropic (§5.1) — mirrors the forward.
- [ ] **The moment metric carries `θ`** and is cross-checked against the production
      `.H` metric (§4.2) — the LD reciprocity teeth depend on it.
- [ ] **V&V tagging:** the algebra-of-record SymPy tests are `@pytest.mark.foundation`
      with NO `verifies()` (math-origin, no `:label:`); the reciprocity/oracle gates
      stay `foundation` (software/algebra invariants). `catches("ERR-066")` stays on
      the cyl_product row ONLY (§8.1). New caught bug → new ERR-NNN, mutation-verified.
- [ ] **Out-of-scope stays loud** — the multi-D `solve_transpose` (G-S), d≥3
      interior-axis, curvilinear LD deferrals still RAISE typed
      `NotImplementedError` (a positive `pytest.raises` control, NOT xfail — the
      loud-flip contract; do not let R2 silently un-defer them).

---

## §13. Consolidated contract-design feedback (task requirement #7 — where verification pressure SHAPES the carve)

Three findings where the stated design must change or be pinned:

1. **The R1a VJP signature `(residual, psi_out) → (psi_bar, psi_in)` is INCOMPLETE
   for curvilinear DD (§1.2).** It omits the M-M angular-numerator cotangent
   (`numer_bar`) the landed hand-code emits. REQUIRED: keep the angular-thread
   VJP on the `PoleAngularClosure` (`cell_contribution` transpose /
   `angular_adjoint`); the registered kernel-pair covers the SPATIAL relation
   only. The sphere/cyl frozen baselines + the C1 Mode-11 sentinel (which must
   fire on BOTH the kernel AND the closure transpose) are the guard. A kernel
   that folds the whole `visit` into `(psi_bar,psi_in)` drops `numer_bar` and
   reds the sphere/cyl baselines — flag this to the implementer at C1 design.

2. **`has_transpose_kernel` should become REGISTRATION-coupled, not a bare bool
   (§2 / §0.2).** Today it is a declared `ClassVar[bool]`. The task's R1a intent
   ("coupled to registration") is correct AND verification-improving: a declared
   `True` with no registered kernel is a predicate lie (the exact class
   `test_ld_adjoint_deferral` guards). Make the trait DERIVE from "a transpose
   kernel is registered" (e.g. a class-level `_transpose_kernel is not None` or
   an override-presence check) so the flag CANNOT desync from the kernel. Pin it:
   a scheme declaring the trait True with no kernel must fail a foundation test.

3. **The LD `.H` metric MUST carry the moment mass `θ` (§4.2 / §7.2).** This is
   a CORRECTNESS requirement on the production `FullFieldSpace` metric for an
   LD field, not just a test config: `G_bulk_LD = V·w_n ⊗ diag(1,θ)`. If the
   production metric is average-only, `.H = G⁻¹AᵀG` is a WRONG adjoint on the
   slope DOF (silent — reciprocity with the matching wrong metric is green).
   The C2 `test_full_field_space_metric_matches_independent_reference` LD row is
   the gate that forces this; if it reds because the production metric is
   average-only, that is a production BUG to fix (ERR-067 "repair the metric"
   family), not a test to relax.

---

## §14. One-paragraph summary

This carve completes the SN adjoint capability along two axes (scheme kernel
+ walk orientation). Every claim is an operator-identity claim — **object-level
gates ONLY, never a spectral functional** (Mode-12; `eig(Aᵀ)=eig(A)`). The
correctness keystones are the structurally-independent OBJECT oracles: the SymPy
`ld_ubld` transpose-oracle (`simplify(VJP−AᵀM⁻¹)==0`), the dense-Euclidean-`Mᵀ`
(M column-probed off the FORWARD apply — a NEW 2-D artifact is required), the
Cartesian assembled-`Mᵀ` (LAPACK `solve_triangular` on the transposed permuted
block), and the frozen pre-carve walk baselines (for the R1a bit-id relocation).
The G-adjoint reciprocity gates are the composite-integration canaries that
prove `.H = G⁻¹AᵀG` composes with the RIGHT metric — and the LD metric MUST
carry the moment mass `θ` or the slope-row transpose sits in reciprocity's
stabiliser (blind) AND `.H` is a wrong adjoint (the load-bearing L18/ERR-067
finding). The land-order is five flip-safe commits; each flag flip is atomic
with its walk + the `test_ld_adjoint_deferral` negative→positive rewrite + the
`is_adjointable`-gated extend-gates. The three contract-design findings (§13)
are where verification pressure reshapes the carve BEFORE the ink dries.
