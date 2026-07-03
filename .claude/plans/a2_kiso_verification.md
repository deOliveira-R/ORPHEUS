# Verification plan — campaign #276 task #119: model-independent iso ENERGY operators (`IsotropicScattering` / `IsotropicN2N`), S† for free

**Branch:** `feature/sn-adjoint-transport` · **Carve:** extract the model-independent
isotropic ENERGY operators on the scalar flux into `orpheus/transport/operators/` (sn-free),
consume at three sites (SN forward / SN adjoint S† / cross-model homogeneous LHS-fold).
Surgical operator-algebra carve **crossing the SN↔model boundary** — proactive
`test-architect` dispatch per the operator-algebra-crossing-subsystems MUST trigger
(per-ordinate ↔ iso-scalar, scalar ↔ moment, normal ↔ adjoint, sn-specific ↔ model-agnostic).

This plan is governed by the standing discipline (`AGENT.md` §0.5/§0.6/§1.5): every gate is
**provably able to red** (a named mutation reddens it, fires under canonical `python -O`), its
reference is **structurally independent** of the SUT, and its regime **activates** the term the
bug lives in. Each gate carries `config / reference / tolerance+justification / claim-layer /
mutation-that-reddens / runtime-mode / vv anti-pattern guarded`.

It SUPERSEDES `.claude/plans/a2_option2_verification.md` (the reverted "bare M₀/R₀/K_iso
iso-frame" SPLIT). The design changed: the iso ℓ=0 projection is no longer decomposed into
`R₀∘K_iso∘M₀` frame faces — instead `K_iso = IsotropicScattering + IsotropicN2N` is a
standalone OperatorSum on the SCALAR flux that routes through the EXISTING `MaterialXSField`
scalar verbs (`apply_p0_in_scatter` / `apply_n2n`, material_xs_field.py:613/652) plus NEW
scalar transpose twins. The reusable parts of the prior spec (the LD-config activation
requirement, the `full_scatter_kernel` oracle, the Euclidean-not-`.H` reciprocity discipline,
the `_CONSUMPTION_TOL=1e-8` LD re-pin) carry forward; the OperatorSum-heterogeneous-space and
no-moment-tensor risks are REPLACED by the new cross-model `Mixture`-vs-`MaterialXSField`
domain crux (RISK-SUM R1).

---

## 0. The carve, restated structurally (so the gates target the right seam)

### 0.1 The SUT — verified in-tree (the symbols do NOT exist yet; this is a pre-implementation plan)

Confirmed by grep: `IsotropicScattering`, `IsotropicN2N`, `apply_p0_in_scatter_transpose`,
`apply_n2n_transpose`, `as_dense` are **NOT present** in `orpheus/transport/`. Correct
pre-implementation posture — the gates define "correct" before the code lands.

Two NEW standalone operators in `orpheus/transport/operators/` (sn-free — the `transport/`
package is already runtime-sn-free per the #261 relocation):

| Operator | `apply(φ)` | `apply_transpose(χ)` | routes through |
|---|---|---|---|
| `IsotropicScattering(mat_xs)` | `Σ_{s0}ᵀφ` (per-cell group matmul) | `Σ_{s0}χ` (group-flip) | fwd: EXISTING `apply_p0_in_scatter` (`einsum "fg,fc...->gc..."`); transpose: NEW `apply_p0_in_scatter_transpose` (`einsum "fg,gc...->fc..."`) |
| `IsotropicN2N(mat_xs)` | `2·Σ_{2n}ᵀφ` (DISTINCT multiplication channel) | `2·Σ_{2n}χ` | fwd: EXISTING `apply_n2n` (`2·einsum "fg,fc...->gc..."`); transpose: NEW `apply_n2n_transpose` |

Both: `domain == codomain ==` the scalar-flux `FunctionSpace`; capabilities
`{CAP_APPLY, CAP_APPLY_TRANSPOSE}`; a **dense/matrix view** (`as_dense` — per-material `(ng,ng)`)
for the LHS-fold consumer. The composite `K_iso = IsotropicScattering(mat_xs) +
IsotropicN2N(mat_xs)` is an `OperatorSum` (both share `mat_xs` + the same scalar space ⇒ the sum
validates natively, RISK-SUM R2).

**The NEW scalar transpose twins mirror the EXISTING moment twins** (already in-tree at
material_xs_field.py:751 `apply_legendre_scattering_moments_transpose` and :840
`apply_n2n_moments_transpose` — both flip `"...fg->...g"` to `"...gf->...f"`). The scalar
twins are the same group-axis flip on the bare `(ng, *spatial)` field:
`apply_p0_in_scatter_transpose` = `einsum("fg,gc...->fc...", sig_s0, χ)`,
`apply_n2n_transpose` = `2·einsum("fg,gc...->fc...", sig2, χ)` — LD-`...`-spectator-safe
(the trailing spatial-moment axis rides through, byte-identical at single-moment closures).

### 0.2 The three consumption sites (verified in-tree)

1. **SN forward** (P2) — `_assemble_per_ordinate_source` (scattering.py:885) currently calls
   `add_iso_source(iso,φ)` + `add_n2n_source(iso,φ)` (lines 925-926), where those delegate to
   `apply_p0_in_scatter` / `apply_n2n`. The carve replaces those two inline calls with
   `(IsotropicScattering + IsotropicN2N).apply(φ)`. **SAME `mat_xs` verbs ⇒ MUST be 0-ULP
   bit-identical** (the `OperatorSum.apply` = `a.apply(x)+b.apply(x)` is the SAME two `+=`
   accumulations into the same scalar accumulator, in the same order). Aniso
   (`build_aniso_source`, ℓ≥1) UNCHANGED. The `/W` stays OUTSIDE (the producer-side
   `(iso/sum_w)+aniso` at line 933 is untouched — Pattern 7 / L18).

2. **SN adjoint S†** (P3, closes #118) — `ScatteringOperator.apply_transpose` is wired to
   `(1/W)·full_scatter_kernel.apply_transpose` (the harmonic-frame form, ALREADY free —
   `full_scatter_kernel` = `frame.conjugate(Λ_{ℓ≥0}+N2N)`, scattering.py:851/878, whose
   transpose falls out of `OperatorProduct.apply_transpose`); `capabilities` gains
   `CAP_APPLY_TRANSPOSE`. NOTE: the brief's "ScatteringOperator.apply_transpose =
   (1/W)·full_scatter_kernel.apply_transpose" routes the S† through the OdRACLE
   (`full_scatter_kernel`), NOT through `K_iso` — `K_iso` is the FORWARD fast-path; S† rides the
   frame form. This is consistent (the perf regression that reverted A2a was on the FORWARD hot
   path; the adjoint is not the SI-sweep hot path). `K_iso`'s OWN transpose is exercised by the
   reciprocity/dense gates (#2, #7c, #9), not by the production S†.

3. **Cross-model validation** (P4) — `homogeneous/solver.py:93-95`
   `self._A = diags(mix.SigT) − SigS0_T − 2.0·Sig2_T` → `A = diags(SigT) − K_iso.as_dense()`
   (the 2nd consumer + the `as_dense` matrix view). **HEADLINE DESIGN RISK** — the homogeneous
   solver holds a bare `Mixture`, NO mesh, NO `MaterialXSField` (RISK-SUM R1).

### 0.3 The oracle (kept per the fuller-view-oracle exception, coding-standards)

`ScatteringOperator.full_scatter_kernel` (scattering.py:851) = `frame.conjugate(Λ_{ℓ≥0}+N2N)` =
the moment-tensor representation of the SAME source, forward AND transpose. It is NOT production
forward (materializes the ℓ=0 tensor — the A2a perf regression) — it is the **permanent
structural cross-check**: harmonic-moment representation vs the production scalar representation
of the ℓ=0 projection. The existing `test_scattering_adjoint.py::TestFullScatterKernel`
already pins `(1/W)·full_scatter_kernel.apply(ψ) ≡ S.apply(ψ)` (forward, rtol=1e-12) and the
full-kernel Euclidean reciprocity — this carve EXTENDS it to an LD multi-moment config and is
the cross-check for the production S† (which IS `full_scatter_kernel.apply_transpose`, so #4
becomes a SELF-equivalence — see RISK-SUM R4).

### 0.4 Structural facts that pin the architecture (verified in-tree)

1. **`OperatorSum.apply` / `.apply_transpose` are `a∘x ± b∘x`** (operator.py — the prior spec
   verified :762/:787). For `K_iso = IsotropicScattering + IsotropicN2N`: forward =
   `Σ_{s0}ᵀφ + 2Σ_{2n}ᵀφ`, transpose = `Σ_{s0}χ + 2Σ_{2n}χ`. The transpose distributes for free
   **iff BOTH summands advertise `CAP_APPLY_TRANSPOSE`** — so the "K_isoᵀ falls out free" claim
   reduces to: each leaf (`IsotropicScattering`, `IsotropicN2N`) advertises it (Gate 9, 12).

2. **Euclidean transpose, NOT the metric Hilbert adjoint `.H`** (L12). `K_isoᵀ` is the plain
   matvec group-flip. The reciprocity gate is Euclidean `⟨K_iso φ, χ⟩ = ⟨φ, K_isoᵀ χ⟩`
   (`.sum()`, **per-group, L27** — NOT weight-summed). `.H` is a DIFFERENT object (carries the
   Gram + `1/W`). The donor `test_scattering_adjoint.py` (the Λᵀ/N2N/full-kernel reciprocity
   tests :103/:222/:273) verifies exactly the Euclidean form — mirror it.

3. **The iso fast-path helpers' callers** — `add_iso_source` / `add_n2n_source` are referenced
   ONLY inside `_assemble_per_ordinate_source` + the `ScalarFlux` arm (scattering.py:1509-1510)
   + the standalone `apply(ScalarFlux)` arm + docstrings. The within-group SI/Krylov driver
   calls `S.apply(...)`, which routes through `_assemble_per_ordinate_source` ⇒ the carve's
   behavioral change flows into the solve via `S.apply` (the solve-level forward gates reach the
   re-routed code — no Mode-11 risk AT THE SOLVE LEVEL). **BUT the LD mutation gate
   monkeypatches `_assemble_per_ordinate_source` DIRECTLY** (test_mms_ld_2d.py:947-956) — the
   re-pin verification is Gate 5 (see RISK-SUM R3).

4. **The homogeneous solver has NO mesh / NO `MaterialXSField`** (verified: solver.py:87-95
   takes a bare `Mixture`, builds `_A` from `mix.SigS[0].T` / `mix.Sig2.T`). The diffusion
   solver is **matrix-FREE and 2G-hardcoded** (`_matvec` at solver.py:240-245 applies
   `(sig_a+sig_s)*fi − sig_s[::-1]*fi[::-1]` — a 2-group flip, NOT a general group-transfer
   matmul; it builds NO explicit `A = diags(Σt) − K_iso`). **So diffusion is NOT a clean
   `as_dense` consumer today** — only homogeneous is. This is the design-risk crux (RISK-SUM R1).

---

## GATE TABLE (the deliverable — paste verbatim in the reply)

Load-bearing gates are flagged **[LB]**. `donor` = the file:line the test is modeled on.
`-O` = fires under `python -O` (function-call assertions only — no bare `assert` in
production/sentinel paths; bare-asserts in collected `tests/` modules ARE rewritten and fire).

| # | Gate | What it pins | Donor (file:line) | Tol + justification | Mutation that reddens | vv anti-pattern guarded |
|---|---|---|---|---|---|---|
| **1** | **[LB] K_iso.apply ≡ current fast-path** | `(IsotropicScattering + IsotropicN2N).apply(φ)` ≡ legacy `add_iso_source(0,φ)+add_n2n_source(0,φ)` on the same φ — per-group, asymmetric SigS, Sig2≠0, P0, **+ LD multi-moment** (trailing 2^d) | `test_scattering_kernel_crosscheck.py:60` (0-ULP pattern) + `test_scattering_adjoint.py:_mix` fixture | **`assert_array_equal` (0-ULP)** — SAME `mat_xs` verbs, SAME two `+=` into the same accumulator, SAME order; an OperatorSum that reorders nothing ⇒ byte-identical, NOT principled-equiv. If the impl introduces ANY reduction reorder (it should not), fall back to `nulp(2)` with the 3 vv criteria documented | swap `Σ_{s0}`↔`Σ_{s0}ᵀ` in IsotropicScattering → RED on asymmetric SigS; drop n2n channel → RED at iso magnitude; double-count iso → RED at 2× | Mode-6 convention drift; Mode-7 (LD row breaks P0-blindness); the bit-id-inheritance anchor |
| **2** | **[LB] K_iso transpose Euclidean reciprocity** | `⟨K_iso φ, χ⟩ = ⟨φ, K_isoᵀ χ⟩` per-group (L27, NOT weight-summed), asymmetric SigS + Sig2≠0; + a dense per-material `Σ @ vec` reference | `test_scattering_adjoint.py:103/116` (Λᵀ reciprocity + dense per-material) | `assert_allclose(rtol=1e-12)` for the inner-product identity; `rtol=1e-13, atol=0` for the dense reference (hand loop ≠ einsum reduction) | swap `Σ[g',g]`→`Σ[g,g']` in `K_isoᵀ` → RED on asymmetric SigS; drop factor-2 on n2n transpose → RED at ½; omit n2n from `K_isoᵀ` → RED; + discriminator (asymmetric SigS ⇒ `K_iso ≠ K_isoᵀ`) | L27 per-group; Mode-2 group-swap; structural-independence (hand loop vs einsum); reciprocity is the broadest single transpose gate |
| **3** | **[LB] SN forward 0-ULP bit-identity** | the WHOLE `_assemble_per_ordinate_source` output (`S.apply(ψ)`) after routing through K_iso ≡ before — P0 + P1 + **LD** | `test_scattering_kernel_crosscheck.py:60` style (capture pre-carve baseline via root-conftest `--capture-baseline`) | **`assert_array_equal`** (0-ULP) — the iso piece is the same two `+=`, the aniso + `/W` are UNTOUCHED; this is a pure re-expression of two inline calls as an OperatorSum | route iso through `full_scatter_kernel` (materializes tensor, reorders) → RED (correct scope-violation signal); drop the n2n leg → RED | L12 unchanged-sibling-stays-0-ULP; the pure-refactor bit-id-inheritance gate |
| **4** | **[LB] SN adjoint S† oracle equivalence + reciprocity** | production `S.apply_transpose ≡ (1/W)·full_scatter_kernel.apply_transpose`; AND `⟨Sψ,χ⟩=⟨φ,Sᵀχ⟩` per-group (asymmetric SigS, Sig2≠0), P0+P1 | EXTEND `test_scattering_adjoint.py:273` (`test_full_kernel_euclidean_reciprocity`) + add an LD row | `assert_allclose(rtol=1e-12)` — the production S† IS `full_scatter_kernel.apply_transpose` (RISK-SUM R4: the oracle-equiv leg is a near-tautology; the LOAD-BEARING leg is the **reciprocity** of `S.apply` vs `S.apply_transpose`, where `S.apply` is the independent FORWARD fast-path) | group-flip in `apply_legendre_scattering_moments_transpose` → RED; omit-n2n from Sᵀ → RED; omit `1/W` → RED at W | L12 transpose-falls-out-free; the reciprocity (fwd fast-path vs adjoint frame-form) is the structurally-independent value pin |
| **5** | **[LB CRITICAL] LD slope-source mutation re-pin (Mode-11)** | the iso slope-source `Σ_s·φ̂` is CONSUMED + sign-correct in production, AND the gate EXECUTES the rewired line | re-pin `test_mms_ld_2d.py:921` (`test_ld_2d_scattering_slope_source_sign_mutation_reddens`) | `_CONSUMPTION_TOL=1e-8` (existing; O(1) above the 1e-12 fixed point, NOT sub-floor) | flip iso slope rows → `|Δφ|/|φ| ≈ 2.6e-3 ≫ tol` → RED; **+ in-process WRAP counter on `IsotropicScattering.apply` must be >0** in the gate's run | **Mode-11** (a re-pinned monkeypatch that no longer executes the production line is VACUOUS) |
| **6** | **[LB] as_dense ≡ apply** | `K_iso.as_dense()` (per-material `(ng,ng)` matrix) applied to φ ≡ `K_iso.apply(φ)` | NEW (the "two consumption modes, one operator" invariant; cite `test_scattering_adjoint.py:116` dense-per-material style) | `assert_allclose(rtol=1e-13, atol=0)` (the dense `M @ φ` reduction tree ≠ the per-material einsum; principled-equiv per vv crit-3) — if the impl makes `as_dense` literally the einsum coefficient, `array_equal` | **transpose the dense block** (`as_dense()` returns `Σ_{s0}` not `Σ_{s0}ᵀ`) → RED on asymmetric SigS; drop the `+2·Σ_2n` from the dense view → RED | the dual-mode invariant (a wrong dense view ships a wrong homogeneous/LHS `A`); Mode-2 |
| **7** | **[LB] homogeneous keff invariance (cross-model)** | after migrating `homogeneous/solver.py` to `K_iso.as_dense()`, `k_inf` bit-identical (or principled-equiv) on the 2eg+4eg analytical cases | `tests/homogeneous/test_homogeneous.py:54` (`test_kinf_exact`, `abs(k_inf−analytical)<1e-12`) | the EXISTING `<1e-12` vs the closed-form analytical k_inf (NOT old-vs-new — the analytical is the structurally-independent ground); the migration's value is inherited from the analytical, so `<1e-12` is BOTH the regression AND the truth gate | build the dense `A` with `Σ_{s0}` not `Σ_{s0}ᵀ` → keff O(1) move on asymmetric SigS → RED; drop n2n from `A` → RED on the n2n-bearing case | anti-#3 (≥2G — 4eg mandatory); the cross-model consumer's truth gate; reference is closed-form (NOT another solver) |
| **8** | **[LB] homogeneous VACUUM/bit-identity correctness gate** | with Sig2=0 (n2n absent), `K_iso.as_dense()` ≡ the legacy `SigS0_T` exactly (the migration's correctness floor); AND with Sig2≠0, `K_iso.as_dense() ≡ SigS0_T + 2·Sig2_T` | NEW (the snapshot-migration VACUUM-bit-id discipline, `snapshot_migration_when_production_goes_bare.md`) | **`assert_array_equal`** legacy `(SigS0_T − ... )` vs `K_iso.as_dense()` for Sig2=0; if they differ → STOP, it's a bug (the migration changed the matrix) | a `K_iso.as_dense()` that orders the group axes wrong, or omits the `.T`, diverges from the legacy `mix.SigS[0].T` → RED | the VACUUM correctness gate (bare ≡ legacy for the zero-multiplication case); Mode-2 axis order |
| **9** | **N2N transpose channel (the multiplication arm)** | `IsotropicN2N.apply_transpose` = `2·Σ_{2n}χ`; included in `K_iso.apply_transpose` | `test_scattering_adjoint.py:222` (`N2NMomentOperator` transpose) + hand loop | `assert_allclose(rtol=1e-13)` + explicit `2·sig2 @ vec` per-material loop | drop factor 2 → RED at ½; un-transpose Σ_2n (`[g,g']`→`[g',g]`) → RED; **omit n2n from `K_iso.apply_transpose` entirely** → RED | #269 non-vacuity (`Sig2≠0` MANDATORY — `make_mixture` AND every A/B/C/D mixture zero Sig2, build on the Mixture directly, lessons L1); L12 channel-in-the-sum |
| **10** | **Aniso 0-ULP canary stays 0-ULP** | the ℓ≥1 `kernel` (`build_aniso_source`) is byte-UNCHANGED — the carve touched only iso | `test_scattering_kernel_crosscheck.py:60` (existing, leave green) | `assert_array_equal` — **do NOT relax** | if the carve accidentally re-routes aniso through a new path → this REDS (correct scope-violation signal) | L12 unchanged-sibling-stays-0-ULP |
| **11** | **Capability + S†-routing sentinel** | `S` advertises `apply_transpose`; `S.apply_transpose` routes through `full_scatter_kernel` (counter>0); `K_iso` + each leaf advertise `CAP_APPLY_TRANSPOSE`; forget a leaf cap → `MissingCapability` downstream | `test_scattering_adjoint.py:96/204` (capability) + Mode-11 in-process wrap | `require(CAP_APPLY_TRANSPOSE in S.capabilities)` + `require(... in K_iso.capabilities)` + counter>0 on `full_scatter_kernel.apply_transpose` (S†) / on a leaf transpose (K_isoᵀ) | ship bespoke inline Sᵀ bypassing the kernel → counter 0 → RED; forget `CAP_APPLY_TRANSPOSE` on `IsotropicScattering` frozenset → `MissingCapability` when K_iso transpose is requested → RED; FLIP any old "S has no transpose" pin | **Mode-11** routing; L12 capability-set/impl mismatch is a distinct catcher from the value gate |
| **12** | **Forward-safety: keff + het-transport + SI≡Krylov + Q/Σ_t** | converged SN keff (2eg+4eg) UNCHANGED; het-transport unchanged; SI and Krylov fixed points still agree; flat-flux `Q/Σ_t` exact (curvilinear) | `test_kinf_homogeneous.py` + `test_heterogeneous_transport.py` (`@catches ERR-025`) + `test_sweep_vs_apply_consistency.py` + `test_streaming_equilibrium_curvilinear.py` | gates' own thresholds (they encode the claim); `Q/Σ_t` via `SAFETY×conv_tol` (L7) | drop/sign-flip iso ℓ=0 → keff O(1) move (4eg reds); an iso change that splits SI/Krylov → consistency reds | anti-#3 (≥2G); the eigenvalue layer cross-validates #1/#3's mutations; the single most powerful curvilinear iso diagnostic (`Q/Σ_t`) |

---

## §RISK-SUM — design risk in the iso-energy-operator extraction (the explicit answer to the brief's three flags)

### R1 — [HIGH, the headline] The cross-model `as_dense` consumer holds a `Mixture`, NOT a `MaterialXSField`

**The brief asks (b): "the as_dense shape contract — per-material vs per-cell vs block-diagonal —
which do homogeneous/diffusion actually need?" The in-tree answer reshapes the carve:**

1. **`homogeneous/solver.py` (verified :87-95) builds `_A` from a bare `Mixture`**
   (`mix.SigS[0].T.tocsr()`, `mix.Sig2.T.tocsr()`) — it has **NO mesh, NO `MaterialXSField`,
   NO `cells_by_material`.** But `IsotropicScattering(mat_xs)` / `IsotropicN2N(mat_xs)` are
   defined ON a `MaterialXSField` (the `mat_xs.apply_p0_in_scatter` / `mat_xs.n2n_matrix`
   surface). So `K_iso.as_dense()` at the homogeneous site requires EITHER:
   - **(a)** a `MaterialXSField` constructed from the single-material `Mixture` (a 0-D /
     single-cell phase space — `MaterialXSField.from_mesh` needs a mesh, so this means minting a
     1-cell `SNMesh` from the `Mixture`, which is a NON-TRIVIAL cross-layer dependency: the
     homogeneous solver would gain an SN-mesh import, **breaking its model-independence**); OR
   - **(b)** `IsotropicScattering` / `IsotropicN2N` gain a **`Mixture`-native constructor**
     (`from_mixture(mix)` → a single-material operator whose `as_dense()` is just
     `{0: mix.SigS[0].toarray()}` / `{0: 2·mix.Sig2.toarray()}`), so the homogeneous solver
     builds `A = diags(mix.SigT) − K_iso.as_dense()[0]` with NO mesh. **This is the
     model-independent path** — and it is the RIGHT one, because the whole POINT of the carve is
     "model-independent iso ENERGY operators" (the brief's framing). The energy operator's
     dense view is a property of the `Mixture`'s `(ng,ng)` matrices, NOT of any mesh.

   **Recommendation (state in the carve crosswalk BEFORE coding):** the iso energy operators
   MUST be constructible from the bare ENERGY data (`Mixture`'s SigS/Sig2), not only from a
   meshed `MaterialXSField`. The cleanest realization: a small `MaterialXS`-like protocol the
   operator depends on, satisfied by BOTH `MaterialXSField` (cell-resolved, SN) AND a 0-D
   `Mixture` adapter (single-material, model-agnostic). Otherwise the "sn-free, model-
   independent" claim is FALSE (the homogeneous solver would import an SN mesh). **If the
   implementer cannot make the operator `Mixture`-constructible, the P4 homogeneous migration
   should be DEFERRED, not forced through an SN-mesh import** — flag this as an
   `AskUserQuestion` checkpoint.

2. **`diffusion/solver.py` is matrix-FREE and 2G-hardcoded** (verified :234-246: `_matvec`
   applies `(sig_a+sig_s)*fi − sig_s[::-1,:]*fi[::-1,:]` — a 2-group flip, builds NO explicit
   `A = diags(Σt) − K_iso`). **The brief's "diffusion/homogeneous build `A = diags(Σt) − K_iso`"
   is FALSE for diffusion** — diffusion has no dense LHS to migrate. So the `as_dense` consumer
   in scope is **homogeneous ONLY** (one consumer). Per coding-elegance "defer abstraction until
   ≥2 instances", `as_dense` has ONE real consumer right now — that is acceptable BECAUSE its
   benefit is established (the homogeneous `A` build IS the consumer), not speculative; but do
   NOT over-build a general block-diagonal/per-cell `as_dense` for a diffusion consumer that
   does not exist. The **shape contract = per-material `dict[int, (ng,ng)]`** (mirroring
   `n2n_matrix` / `sig_s_legendre`), and the homogeneous single-material case reads `[0]`. A
   per-cell or block-diagonal `(N·ng, N·ng)` dense view is YAGNI until a real LHS-fold consumer
   (CP net-removal? — but CP is also not a dense-`A` builder today) appears.

### R2 — [LOW] Does `IsotropicScattering + IsotropicN2N` as an OperatorSum trip the domain/codomain guard? (the brief's flag (a))

**NO — they share `mat_xs` and the same scalar-flux `FunctionSpace`, so the sum validates
natively.** Unlike the reverted SPLIT (where `S_iso`'s per-ordinate outer space had to match
`S_aniso`'s `L2[S^2]/(N,)` and risked a fresh-minted name mismatch), here BOTH summands are the
SAME map on the SAME scalar space: `domain == codomain ==` the scalar-flux space, derived from
the SAME `mat_xs`. The `OperatorSum.__init__` guard (`a.domain == b.domain` AND
`a.codomain == b.codomain` by `FunctionSpace (name, shape)`, skipped when `None`) passes by
construction. **The design choice that matters:** give BOTH operators their scalar space from
ONE source (the `mat_xs`'s scalar-flux space, Pattern 2 single-source) — do NOT let each mint
its own. Add a LIGHT positive-construction gate (`require(isinstance(K_iso, OperatorSum))`
constructs without raising) — but a negative-mismatch test is LOW value here (both summands are
the same space by construction; there is no convention to drift), so it is NOT load-bearing
(unlike the SPLIT, where the heterogeneous-space mismatch was the headline risk). The risk
**moved** from the OperatorSum-space-guard (SPLIT) to the cross-model `Mixture` domain (R1).

### R3 — [MEDIUM, the Gate-5 trap] Does the LD monkeypatch target still execute after the K_iso routing? (the brief's flag (c))

**YES — the monkeypatch target SURVIVES, but the proof must be re-confirmed.** The LD gate
(test_mms_ld_2d.py:947) monkeypatches `ScatteringOperator._assemble_per_ordinate_source` and
flips `out.values[..., 1:]` (the slope rows of the COMBINED per-ordinate source). The revised
carve keeps `_assemble_per_ordinate_source` as the production combiner — it just replaces the
inline `add_iso_source`+`add_n2n_source` (lines 925-926) with `(IsotropicScattering +
IsotropicN2N).apply(φ)`. The monkeypatch wraps the OUTER method (flips its OUTPUT), so it STILL
bites regardless of the internal iso routing. **HOWEVER**: (i) the monkeypatch flips the
combined `out`, not the iso piece specifically — it remains a valid consumed-sign gate (the iso
slope IS in `out`); (ii) the Mode-11 sentinel must confirm the NEW `IsotropicScattering.apply`
line actually executes in the gate's run (an in-process WRAP counter on `IsotropicScattering.apply`
> 0) — because if the carve accidentally left the iso slope flowing through a DIFFERENT path
(e.g. the windowed `HarmonicMomentFlux` arm at scattering.py:1601 also calls
`_assemble_per_ordinate_source` but feeds `phi_moments.scalar_flux()`), the gate could be green
while `IsotropicScattering` is never entered for the slope. **Gate 5's deliverable:** keep the
monkeypatch on `_assemble_per_ordinate_source` (it still works), ADD the in-process wrap on
`IsotropicScattering.apply` (counter>0), and mutation-verify the slope flip reddens
(`|Δφ|/|φ| ≈ 2.6e-3 ≫ 1e-8`). A green re-pinned gate whose wrap counter is 0 is VACUOUS.

   **Sharper alternative (stronger, recommended):** ADD a SECOND mutation gate that monkeypatches
   `IsotropicScattering.apply` (or `apply_p0_in_scatter`) DIRECTLY to flip its slope rows — this
   pins the NEW production reader specifically (the existing `_assemble_per_ordinate_source`
   monkeypatch is now one layer ABOVE the changed code, so a stamp bug in `IsotropicScattering`
   that the combiner faithfully passes through could be masked). The direct monkeypatch +
   counter is the Mode-11 gold-standard "the gate executes the rewired line".

### R4 — [MEDIUM] The S† oracle-equivalence leg (#4) is a near-tautology — the reciprocity is the real gate

The brief wires production `S.apply_transpose = (1/W)·full_scatter_kernel.apply_transpose`. So
"`S.apply_transpose ≡ (1/W)·full_scatter_kernel.apply_transpose`" is **by construction** (the
production IS that expression) — a self-equivalence, NOT a structurally-independent check. The
LOAD-BEARING value pin for S† is therefore the **Euclidean reciprocity** `⟨S.apply(ψ), χ⟩ =
⟨ψ, S.apply_transpose(χ)⟩`, where `S.apply` is the INDEPENDENT forward fast-path (the iso piece
goes through `K_iso`/`add_iso_source`, the aniso through `build_aniso_source`) and
`S.apply_transpose` goes through the frame `full_scatter_kernel` — two structurally-different
representations of the same operator, so reciprocity genuinely cross-checks them (the forward
scalar matmul vs the adjoint moment-tensor transpose). Gate 4's reciprocity leg is the teeth;
the oracle-equiv leg is a cheap sanity (catches a typo in the wiring, e.g. a missing `1/W` or a
`full_scatter_kernel.apply` instead of `.apply_transpose`). State this in the gate so the
implementer does not mistake the tautology for the verification.

### R5 — [LOW] The production S† rides the ORACLE, but K_iso's OWN transpose still needs a gate

Because production S† = `full_scatter_kernel.apply_transpose` (the frame form), `K_iso`'s
transpose (`IsotropicScattering.apply_transpose + IsotropicN2N.apply_transpose`) is NOT on the
production S† path — it is the FORWARD operator's transpose, exercised by Gates 2/6/9 (dense
reference, reciprocity, n2n channel) but NOT by the solve. This is fine (the carve builds K_iso
as a composable operator with a transpose for COMPLETENESS / future scalar-carrier consumers,
the same way the orphan `ScalarFlux` arm is kept per #205), but it means the K_iso transpose
gates (#2/#6/#9) are UNIT-level (against a hand loop / dense reference), NOT end-to-end — declare
them so. If a future consumer wires K_isoᵀ into a scalar adjoint, an end-to-end gate is added
then; today the structurally-independent dense per-material reference is the correct verification.

---

## §ACTIVATION — the term-activation audit (Mode-7/Mode-10, AGENT.md §0.6)

The Mode-7 central risk: a P0-only fixture NEVER builds the ℓ≥1 blocks and is BLIND to the
iso/aniso split AND to the S† aniso transpose. **EVERY equivalence/correctness/transpose gate
that exercises the full operator (#1, #2, #3, #4, #9, #11) runs on the ANISOTROPIC P1 +
heterogeneous + `Sig2≠0` `solver_p1_het` fixture** (already in the donor — asymmetric `_P0_A/B`,
`_P1_A/B`, `Sig2=[[0,.03],[.01,0]]`, 2 materials, 2G). The **LD multi-moment** config
(`build_2d_cartesian_ld_stress_mms_case`, trailing 2^d) is mandatory for #1/#3/#4/#5 — the
reverted A2a's regression was LD-only, and the iso slope-source `Σ_s·φ̂` threads through
`IsotropicScattering.apply` → `apply_p0_in_scatter` (the `...` spatial-moment spectator). The
P0-only config is used ONLY where the iso-scalar arm is the point (the homogeneous cross-model
gates #7/#8 are 2eg/4eg, which exercise multi-group group-transfer even at P0 — the eigenvalue
is multi-group, so the iso ℓ=0 group-coupling is genuinely active). **`Sig2≠0` is built DIRECTLY
on the Mixture** (`_mix` in the donor sets `m.Sig2 = csr_matrix(...)`), NOT via `make_mixture`
(which zeros Sig2 AND SigL — lessons L1). For the homogeneous keff gates, build a 2eg/4eg case
with non-zero `Sig2` on the `Mixture` constructor (reuse the `_balanced_fissile_4g` hand-built
pattern or `_make_mixture_with_n2n` from `tests/cp/test_verification.py`) so the n2n channel of
the migrated `A` is non-vacuous.

**Per-group, NOT weight-summed (L27):** every reciprocity / dense-reference gate contracts the
full `(ng, *spatial)` tensor (`.sum()` over everything, or a per-group residual), NEVER a
weight-summed scalar that telescopes (a wrong per-group transpose can satisfy a global balance
by construction — vv §H3 / anti-#8).

---

## §REBASELINE — disposition of EXISTING gates (explicit per §1.5)

| Existing gate | Carve disposition | Tolerance after carve |
|---|---|---|
| `test_scattering_adjoint.py::TestLambda/Kernel/N2N/FullScatterKernel` (foundation) | **STAYS GREEN unchanged** (the moment-space leaves + the `full_scatter_kernel` oracle the S† reuses) | `rtol=1e-12..1e-13` (unchanged) |
| `test_scattering_adjoint.py::TestFullScatterKernel::test_reproduces_forward_scattering_source` + `::test_full_kernel_euclidean_reciprocity` | **EXTEND** to an LD multi-moment row (Gate #4); the production S† wiring makes the reciprocity the live S† value gate | `rtol=1e-12` (unchanged; +LD row) |
| `test_scattering_kernel_crosscheck.py` (ℓ≥1 `kernel`) | **STAYS 0-ULP** (aniso UNCHANGED) — Gate #10 | `assert_array_equal` (unchanged) |
| `test_mms_ld_2d.py::test_ld_2d_scattering_slope_source_sign_mutation_reddens` | **RE-PIN + Mode-11 wrap** on `IsotropicScattering.apply` (Gate #5); ADD a direct-on-`IsotropicScattering` mutation (R3 sharper alternative) | `_CONSUMPTION_TOL=1e-8` (unchanged; +wrap counter) |
| `tests/homogeneous/test_homogeneous.py::test_kinf_exact` (2eg+4eg) | **STAYS GREEN unchanged** after the `K_iso.as_dense()` migration (Gate #7); the analytical k_inf is the truth gate the migrated `A` must still hit | `<1e-12` vs analytical (unchanged) |
| `test_kinf_homogeneous.py` / `test_heterogeneous_transport.py` (SN forward) | **STAY GREEN unchanged** (forward-preservation, Gate #12) | gates' own thresholds |
| `test_scattering_operator.py::TestBitIdenticalExtractionP0` (iso fast-path, if present) | **KEEP** — `add_iso_source`/`add_n2n_source` survive as the verbs `IsotropicScattering`/`IsotropicN2N` call (NOT retired); the public wrappers may become thin shims over `apply_p0_in_scatter`/`apply_n2n` | `rtol=1e-13`, teeth re-confirmed |
| **NEW** Gate #1 (K_iso≡fast-path), #2 (K_iso transpose recip + dense), #6 (as_dense≡apply), #8 (homogeneous VACUUM bit-id), #9 (n2n transpose channel), #11 (cap+routing), R1 `Mixture`-constructibility | **NET-NEW** | `array_equal` / `rtol=1e-12..1e-13` / counter>0 |

**Net-new test surface:** add scalar transpose verbs `apply_p0_in_scatter_transpose` /
`apply_n2n_transpose` to `material_xs_field.py` (mirror the moment twins at :751/:840); add
`IsotropicScattering` / `IsotropicN2N` operators + their `as_dense` + (R1) a `Mixture`-native
constructor in `orpheus/transport/operators/`; extend `test_scattering_adjoint.py` (Gate
#1/#2/#4-LD/#6/#9/#11); the LD re-pin in `test_mms_ld_2d.py`; the homogeneous migration gates in
`tests/homogeneous/`. The `S` dataclass gains `CAP_APPLY_TRANSPOSE` (Gate #11 FLIPS any "S has
no transpose" pin).

**Catalog:** if the carve catches a real bug (the n2n-omission-from-K_isoᵀ, the `as_dense`
transposed-block, the homogeneous `Σ_{s0}` vs `Σ_{s0}ᵀ` axis order, the LD spatial-moment slip,
the cross-model `Mixture`-mesh dependency), log an ERR-NNN per the `vv-principles` "Log every
caught bug" directive — the cross-model dense-view axis-order and the `as_dense ≠ apply`
mismatch are new failure-mode candidates not yet catalogued.
