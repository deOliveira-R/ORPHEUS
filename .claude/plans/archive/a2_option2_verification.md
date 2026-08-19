> **SUPERSEDED by a2_kiso_verification.md — archived 2026-08-19 (plans triage).**

# Verification plan — campaign #276 A2-Option2 (#119): composable + fast + LD-correct iso scattering, S† for free

**Branch:** `feature/sn-adjoint-transport` · **Carve:** A2-Option2 (forward iso/aniso
SPLIT as two composed `LinearOperator`s; A2b transpose distributes over the sum for free).
Surgical operator-algebra carve, **main-agent-direct**. Proactive `test-architect` dispatch
per the operator-algebra-crossing-subsystems trigger (per-ordinate ↔ iso-scalar, scalar ↔
moment, normal ↔ adjoint, fast-path ↔ composed-operator).

This plan is governed by the standing discipline (`AGENT.md` §0.5/§0.6/§1.5): every gate is
**provably able to red** (a named mutation reddens it, fires under the canonical `python -O`),
its reference is **structurally independent** of the SUT, and its regime **activates** the term
the bug lives in. Each gate row carries `config / reference / tolerance / claim-layer /
mutation-that-reddens / runtime-mode / vv-anti-pattern guarded`.

It SUPERSEDES `.claude/plans/archive/a2_iso_modernization_verification.md` (the reverted A2a, which
routed the WHOLE forward through `frame.conjugate(Λ_{skip_l0=False})` and "regressed LD/P0
badly" — scattering.py:1536-1539). The A2a foundation (`Λᵀ`, `N2NMomentOperator`,
`full_scatter_kernel`, `tests/sn/operators/test_scattering_adjoint.py`) is **already committed
and green** (`e8d9f6b`/`f3f3dd2`/`0b3275d`; baseline 23 passed under `-O`). Option-2 reuses
that foundation and changes the FORWARD `apply` arms.

---

## 0. The carve, restated structurally (so the gates target the right seam)

### 0.1 Today (the committed branch state — verified in-tree, NOT from memory)

The forward `S.apply(AngularFlux)` arm (scattering.py:1540-1542) is the **fast path**:

```
S.apply(ψ) = _assemble_per_ordinate_source(ψ.integrate_angular(), build_aniso_source(ψ), ψ.mesh)
           = (iso / W) + aniso
  where iso   = add_iso_source(0,φ) + add_n2n_source(0,φ)   # φ=∫ψdΩ, scalar reaction-rate matmul
        aniso = (1/W)·kernel.apply(ψ)                        # kernel = frame.conjugate(Λ_{skip_l0=True})
```

`_assemble_per_ordinate_source` (scattering.py:885) is the SHARED `/W` combine for BOTH the
`AngularFlux` arm (1540) AND the `HarmonicMomentFlux` windowed arm (1601). `S` (the
`@dataclass`) advertises `{CAP_APPLY}` ONLY (verified: `S default caps: frozenset({'apply'})`).
`full_scatter_kernel` = `frame.conjugate(Λ_{skip_l0=False} + N2N)` is the committed FULL frame-form
ORACLE with `apply_transpose` already free (851-883), pinned by
`test_scattering_adjoint.py::TestFullScatterKernel`.

### 0.2 Option-2 (the system under test)

The forward `apply` arm is re-expressed as TWO composed `LinearOperator`s, summed, so the
transpose distributes for free:

```
S.apply(ψ)            = (1/W)·[ S_iso.apply(ψ)             + S_aniso.apply(ψ) ]
S.apply_transpose(ψ*) = (1/W)·[ S_iso.apply_transpose(ψ*) + S_aniso.apply_transpose(ψ*) ]
```

| Channel | Operator | apply | apply_transpose |
|---|---|---|---|
| `M₀` iso analysis | ℓ=0 angular integral | `φ = Σ_n w_n ψ_n` `(N,ng,*sp)→(ng,*sp)` | `(M₀ᵀφ)_n = w_n·φ` |
| `K_iso` iso in-scatter + n2n | scalar group-transfer | `Σ_{s0}ᵀφ + 2Σ_{2n}ᵀφ` | `Σ_{s0}χ + 2Σ_{2n}χ` |
| `R₀` iso broadcast | ℓ=0 reconstruction | `φ→φ` on every ordinate `(ng,*sp)→(N,ng,*sp)` | `Σ_n χ_n` (un-weighted sum) |
| `S_iso` | `OperatorProduct(R₀, OperatorProduct(K_iso, M₀))` | `(R₀∘K_iso∘M₀)ψ` (NO moment tensor) | `(M₀ᵀ∘K_isoᵀ∘R₀ᵀ)ψ*` |
| `S_aniso` | `= ScatteringOperator.kernel` (UNCHANGED) | `(R∘Λ_{ℓ≥1}∘M)ψ` | free via OperatorProduct |

For P0-only (`scattering_order=0`), `S_aniso` is absent; `S = (1/W)·S_iso` (pure scalar, no
tensor). **`M₀`/`K_iso`/`R₀` route through the EXISTING `MaterialXSField` scalar verbs**
(`apply_p0_in_scatter`/`apply_n2n` at material_xs_field.py:613/652 — the `...` spatial-moment
spectator broadcast is already there; the scalar transpose twins
`apply_p0_in_scatter_transpose`/`apply_n2n_transpose` are to be ADDED).

### 0.3 The oracle (kept per the fuller-view-oracle exception, coding-standards)

`ScatteringOperator.full_scatter_kernel` (scattering.py:851) = `frame.conjugate(Λ_{ℓ≥0}+N2N)`
= the moment-tensor representation of the SAME source. It is NOT production (materializes the
tensor for ℓ=0 — the exact perf regression that reverted A2a); it is the **permanent structural
cross-check**: harmonic-moment representation vs the production scalar representation of the ℓ=0
projection. The existing `test_scattering_adjoint.py::TestFullScatterKernel` already pins
`(1/W)·full_scatter_kernel.apply(ψ) ≡ S.apply(ψ)` (forward) — Option-2 EXTENDS it to the
transpose AND to an LD multi-moment config.

### 0.4 Structural facts that pin the architecture (verified in-tree)

1. **`OperatorSum.apply_transpose` propagates iff BOTH operands have `CAP_APPLY_TRANSPOSE`**
   (operator.py:762-763, 787-788: `a.apply_transpose(x) + b.apply_transpose(x)` — `(A+B)ᵀ=Aᵀ+Bᵀ`).
   `S_aniso = kernel` ALREADY advertises it (`test_kernel_advertises_apply_transpose` green). So
   the WHOLE "S† falls out free" claim reduces to: **`S_iso` advertises `CAP_APPLY_TRANSPOSE`**,
   which (since `S_iso = OperatorProduct(R₀, OperatorProduct(K_iso, M₀))`) requires `M₀`, `K_iso`,
   `R₀` to EACH advertise it (`OperatorProduct.apply_transpose`, operator.py:841-842/862-864). The
   load-bearing NEW leaves are `M₀ᵀ`, `K_isoᵀ`, `R₀ᵀ`.

2. **Euclidean transpose, NOT the metric Hilbert adjoint `.H`** (L12). The carve wires
   `apply_transpose = (1/W)·[S_iso.apply_transpose + S_aniso.apply_transpose]` — the plain matvec
   `Sᵀ`. The reciprocity gate is Euclidean `⟨Sψ,χ⟩=⟨φ,Sᵀχ⟩` (`.sum()`, per-group, L27). `.H` is
   a DIFFERENT object (`_AdjointOperator.apply = G⁺⊙apply_transpose(G⊙y)`, operator.py:675-695,
   carries the SH `(2ℓ+1)` Gram + `1/W`). The fission template (`test_fission_adjoint.py`) and the
   existing `test_scattering_adjoint.py` both verify exactly the Euclidean form — mirror it.

3. **The iso fast-path helpers have NO external production caller**
   (`_assemble_per_ordinate_source` / `add_iso_source` / `add_n2n_source` / `build_aniso_source`
   are referenced ONLY inside scattering.py + docstrings + the LD MMS monkeypatch + tests). The
   within-group SI/Krylov driver calls **`S.apply(...)`** (the `@singledispatchmethod` arms). ⇒
   Option-2's behavioral change flows into the solve via `S.apply` — the solve-level
   forward-preservation gates genuinely reach the re-routed code (no Mode-11 risk AT THE SOLVE
   LEVEL). **BUT the LD mutation gate monkeypatches `_assemble_per_ordinate_source` DIRECTLY**
   (test_mms_ld_2d.py:947-956) — if Option-2 retires/reshapes that method, that gate goes
   vacuous (Gate 3, the CRITICAL re-pin).

4. **The angular `measure_space` is `FunctionSpace(name="L2[S^2]", shape=(N,))`** (measure.py:263)
   — per-ordinate `(N,)` values. `S_aniso = frame.conjugate(Λ)` reports `domain == codomain ==
   measure_space == L2[S^2]/(N,)` (conjugate wraps M/R, so the OUTER spaces are the per-ordinate
   nodal space, NOT the SH coefficient space). `SphericalHarmonicSpace.from_L(L)` is
   `name="spherical_harmonic_space", shape=(L+1,2L+1)`. `FunctionSpace.__eq__` is by `(name,shape)`
   (space.py:157-160). **This is the design-risk crux — see §RISK-SUM.**

---

## GATE TABLE (the deliverable — paste verbatim in the reply)

Load-bearing gates are flagged **[LB]**. `donor` = the file:line the test is modeled on.
`-O` = fires under `python -O` (function-call assertions only — no bare `assert` in
production/sentinel paths; `tests/` bare-asserts ARE rewritten and OK).

| # | Gate | What it pins | Donor (file:line) | Tol + justification | Mutation that reddens | vv anti-pattern guarded |
|---|---|---|---|---|---|---|
| **1** | **[LB] Oracle equivalence — forward** | `(1/W)·full_scatter_kernel.apply(ψ) ≡ S.apply(ψ)` per-group (L27), P0 + P1 + **LD multi-moment** | `test_scattering_adjoint.py:249` (`test_reproduces_forward_scattering_source`) | `assert_allclose(rtol=1e-12, atol=0)` — Y₀⁰=1 ⇒ same math, different reduction tree (scalar matmul vs `R₀(K_iso(M₀ψ))` + `(2ℓ+1)|ℓ=0` broadcast); drift = reduction-depth·ULP (`vv §bit-id` crit-3), NOT 0-ULP | flip ℓ=0 sign in K_iso → RED; drop iso (S_aniso only) → RED at iso magnitude; double-count iso → RED at 2× | L12 fold-is-principled-equiv; Mode-7 (LD row breaks P0-blindness) |
| **2** | **[LB] Oracle equivalence — transpose** | `(1/W)·full_scatter_kernel.apply_transpose(χ) ≡ S.apply_transpose(χ)` per-group, P0 + P1 + LD | NEW in `test_scattering_adjoint.py` (twin of `test_full_kernel_euclidean_reciprocity:273`) | `assert_allclose(rtol=1e-12)` (same as #1; the oracle's transpose is the structurally-independent moment-tensor `Sᵀ`) | un-transpose K_iso (Mode-2) → RED; drop n2n from Sᵀ → RED; omit `1/W` → RED at W | L12 transpose-falls-out-free; the oracle pins VALUE |
| **3** | **[LB CRITICAL] LD slope-source mutation re-pin** | the iso slope-source `Σ_s·φ̂` is CONSUMED + sign-correct in production, AND the gate EXECUTES the rewired line | re-pin `test_mms_ld_2d.py:921` (`test_ld_2d_scattering_slope_source_sign_mutation_reddens`) | `_CONSUMPTION_TOL=1e-8` (existing; O(1)-above-fixed-point, NOT sub-floor) | flip iso slope rows → `|Δφ|/|φ| ≈ 2.6e-3 ≫ tol` → RED; **+ in-process WRAP counter on the rewired reader must be >0** | **Mode-11** (a re-pinned monkeypatch that no longer executes the production line is VACUOUS) |
| **4** | **[LB] Iso bit-identity (S_iso ≡ fast path)** | `S_iso.apply(ψ)/W` ≡ legacy `(add_iso+add_n2n on ∫ψ)/W` (same ops re-expressed) | `test_scattering_kernel_crosscheck.py:135` (0-ULP pattern) | `assert_array_equal` (0-ULP) where the reduction tree is preserved; if `M₀`/`R₀` reorder the sum → `nulp(reduction_depth)` with the 3 `vv` criteria documented | swap `Σ_{s0}`↔`Σ_{s0}ᵀ` in K_iso → RED on asymmetric SigS; drop n2n → RED | Mode-6 convention drift; principled-equiv ladder |
| **5** | **Aniso 0-ULP canary stays 0-ULP** | the ℓ≥1 `kernel` (S_aniso) is byte-UNCHANGED | `test_scattering_kernel_crosscheck.py:135` (existing, leave green) | `assert_array_equal` — **do NOT relax** | if Option-2 accidentally re-routes aniso through a new sum → this REDS (correct scope-violation signal) | L12 unchanged-sibling-stays-0-ULP |
| **6** | **[LB] Free-transpose Euclidean reciprocity** | `⟨Sψ,χ⟩=⟨φ,Sᵀχ⟩` per-group (NOT `.H`), the adjoint-DEFINING law | `test_fission_adjoint.py` + `test_scattering_adjoint.py:175` style | `assert_allclose(rtol=1e-12)` (fission template threshold) | group-flip K_iso → RED; omit-n2n → RED; un-transpose Λ → RED; +discriminator (asymmetric SigS ⇒ S≠Sᵀ) | L27 per-group; Mode-2; reciprocity is the broadest single transpose gate |
| **7a** | **Intrinsic law — M₀ᵀ** | `⟨M₀ψ,φ⟩=⟨ψ,M₀ᵀφ⟩` (M₀ᵀ is the Euclidean transpose of the ℓ=0 angular integral) | `test_scattering_adjoint.py:103` (`test_moment_space_transpose_identity`) | `assert_allclose(rtol=1e-12)` | `M₀ᵀφ = φ` (forget the `w_n` weight) → RED | anti-#7 (test a math type's intrinsic law); Mode-3 (missing weight) |
| **7b** | **Intrinsic law — R₀ᵀ** | `⟨R₀φ,χ⟩=⟨φ,R₀ᵀχ⟩` (R₀ᵀ is the un-weighted ordinate sum) | same donor | `assert_allclose(rtol=1e-12)` | `R₀ᵀχ = w_n·Σχ` (wrongly weight the transpose) → RED | anti-#7; the M₀/R₀ asymmetry (M₀ weights, R₀ doesn't — Pattern-7 `/W` lives OUTSIDE) |
| **7c** | **Intrinsic law — K_isoᵀ** | `⟨K_isoψ,χ⟩=⟨ψ,K_isoᵀχ⟩` (scalar group-transfer transpose) | `test_scattering_adjoint.py:116` (`test_transpose_matches_dense_per_material`) | `assert_allclose(rtol=1e-13)` + dense per-material `sig @ vec` reference | use forward index `Σ[g',g]` in transpose → RED on asymmetric SigS; +discriminator | Mode-2 group-swap; structural-independence (hand loop ≠ einsum) |
| **7d** | **Iso-frame identity** | `M₀` relates to the ℓ=0 row of the harmonic analysis M (`Y₀⁰=1`): `M₀.apply(ψ) ≡ M.apply(ψ)[0,0]` | NEW (cites the `Y₀⁰=1` convention single-sourced at `HarmonicMomentFlux.scalar_flux`) | `assert_array_equal` (0-ULP — it IS the same integral) | scale M₀ by ≠1 → RED vs the M[0,0] moment | the convention pin (Mode-4 normalization); proves S_iso is the ℓ=0 projection |
| **8a** | **[LB] P0 no-moment-tensor sentinel** | for `scattering_order=0`, production `S.apply` does NOT call `apply_legendre_scattering_moments` / `apply_n2n_moments` (the moment-tensor verbs) | `vv` Mode-11 sharpening (in-process autouse wrap) | counter — assert **0 calls** to the tensor verbs at L=0 | route P0 iso through `full_scatter_kernel` (the reverted A2a) → counter >0 → RED | **Mode-11** + the load-bearing PERF invariant (the reason A2a reverted) |
| **8b** | **P0 forward equivalence** | `S.apply(ψ)` at L=0 ≡ legacy `(iso/W)` (pure scalar path, no tensor) | `test_scattering_operator.py` P0 fixtures | `assert_array_equal` (0-ULP — same scalar ops) | (covered by #1's P0 row) | Mode-7 (P0-only is the BLIND config for ℓ≥0 — used HERE deliberately for the perf claim only) |
| **9** | **N2N transpose (the channel in the iso sum)** | n2n transpose `(2Σ_{2n})ᵀ` is included in `S_iso.apply_transpose` (n2n rides K_iso) | `test_scattering_adjoint.py:222` (`N2NMomentOperator` transpose) | `assert_allclose(rtol=1e-13)` + hand loop | drop factor 2 → RED at ½; un-transpose Σ_2n → RED; **omit n2n from Sᵀ entirely** → RED | L12 channel-in-the-fold; #269 non-vacuity (`Sig2≠0` mandatory) |
| **10** | **Forward-safety: keff** | converged keff UNCHANGED to solver tol | `test_kinf_homogeneous.py` (2eg+4eg) + `test_heterogeneous_transport.py` (`@catches ERR-025`) | gates' own thresholds (they encode the claim) | drop/sign-flip iso ℓ=0 → keff O(1) move; confirm 4eg reds | anti-#3 (≥2G); the eigenvalue layer cross-validates #1's mutations |
| **11** | **Forward-safety: SI≡Krylov + Q/Σ_t** | SI and Krylov fixed points still agree; flat-flux `Q/Σ_t` exact (curvilinear) | `test_sweep_vs_apply_consistency.py` + `test_streaming_equilibrium_curvilinear.py` | gates' own thresholds; `Q/Σ_t` via `SAFETY×conv_tol` (L7) | double-`/W` → `Q/Σ_t` off by W; iso change that splits SI/Krylov → consistency reds | R-double-`/W` (Pattern-7 crosswalk); the single most powerful curvilinear iso diagnostic |
| **12** | **Capability + S† routing sentinel** | `S` advertises `apply_transpose` AND `S.apply_transpose` ROUTES through `S_iso`/`S_aniso` leaves | `test_scattering_adjoint.py:96` (capability) + Mode-11 wrap | `require(CAP_APPLY_TRANSPOSE in S.capabilities)` + counter>0 on a leaf transpose | ship bespoke inline Sᵀ bypassing the sum → counter 0 → RED; forget `CAP_APPLY_TRANSPOSE` on a leaf frozenset → `MissingCapability` downstream → RED | **Mode-11** routing; L12 capability-set-mismatch is a distinct catcher |

---

## §RISK-SUM — design risk in the Option-2 decomposition (the explicit answer)

### R1 — [HIGH, the headline] Can `OperatorSum(S_iso, S_aniso)` trip the domain/codomain equality guard?

**YES, this is the load-bearing design risk, and it must be resolved BEFORE coding.** The
`OperatorSum.__init__` guard (operator.py:744-757) raises `IncompatibleOperatorComposition` iff
`a.domain != b.domain` OR `a.codomain != b.codomain`, **skipped only when a side is `None`**
(operator.py:742-743 read via `getattr(..., None)`). The relevant `FunctionSpace`s
(`__eq__` by `(name, shape)`, space.py:157):

- `S_aniso = ScatteringOperator.kernel = frame.conjugate(Λ_{ℓ≥1})`. Because `conjugate` wraps
  `R`/`M` around Λ, its OUTER spaces are the per-ordinate angular nodal space:
  `S_aniso.domain == S_aniso.codomain == frame.measure_space == FunctionSpace("L2[S^2]", (N,))`
  (measure.py:263; verified).

- `S_iso = OperatorProduct(R₀, OperatorProduct(K_iso, M₀))`. `OperatorProduct.domain` =
  `b.domain` = `M₀.domain`; `OperatorProduct.codomain` = `a.codomain` = `R₀.codomain`
  (operator.py:846-853). M₀ maps per-ordinate→ℓ=0-scalar; R₀ maps ℓ=0-scalar→per-ordinate. So
  `S_iso.domain = M₀.domain` (per-ordinate) and `S_iso.codomain = R₀.codomain` (per-ordinate).

**The sum constructs cleanly IFF `M₀.domain` and `R₀.codomain` are EITHER:**
(a) **EQUAL to `FunctionSpace("L2[S^2]", (N,))`** by `(name, shape)` — i.e. M₀/R₀ declare the
SAME per-ordinate angular space `S_aniso` uses; OR
(b) **`None`** — the bare/legacy contract (`ScatteringOperator.domain` is already `None` unless
`full_field_space` is threaded; the LD/test path is `None`-spaced today).

**Three concrete failure modes the implementer must avoid:**

1. **M₀/R₀ mint a NEW per-ordinate space name** (e.g. `"per_ordinate_source"` vs `"L2[S^2]"`).
   Same shape `(N,)`, DIFFERENT name ⇒ `__eq__` False ⇒ the sum RAISES. **Fix:** M₀.domain /
   R₀.codomain MUST report the SAME `frame.measure_space` instance/name `S_aniso` reports — i.e.
   build M₀/R₀ FROM the frame's measure (reuse `frame.measure.space`), do not fabricate a parallel
   space. This is the coding-elegance Pattern-2 single-source move: the per-ordinate space has ONE
   definition (the angular measure), and BOTH summands read it.

2. **K_iso's intermediate (ℓ=0 scalar) space leaks to the sum's OUTER domain.** It does not —
   `OperatorProduct.domain` reads the INNER `b.domain` (= M₀.domain, per-ordinate), and `Cmid`
   (the K_iso scalar space) is captured internally (operator.py:810-813). So the K_iso scalar
   space is invisible to the OperatorSum guard. **But** the INNER `OperatorProduct(K_iso, M₀)`
   guard checks `K_iso.domain == M₀.codomain` (the ℓ=0-scalar space) — M₀.codomain and K_iso.domain
   must agree (both the scalar ℓ=0 space). And `OperatorProduct(R₀, ·)` checks `R₀.domain ==
   K_iso.codomain` (scalar). Mint ONE scalar ℓ=0 space, used by M₀.codomain, K_iso.domain,
   K_iso.codomain, R₀.domain.

3. **Leave M₀/R₀ spaces `None` (path (b)).** The sum then constructs (guard skipped), BUT the
   composability guard NO LONGER VALIDATES the sum — and, more importantly, the inner
   `OperatorProduct(K_iso, M₀)` / `OperatorProduct(R₀, ·)` guards also skip, so a real M₀.codomain
   ≠ K_iso.domain mismatch (a wrong scalar shape) ships silently. **This is a legitimate
   transitional choice** (matches the `None`-spaced `ScatteringOperator` default and the
   `full_scatter_kernel` oracle, which works on the bare-ndarray path) — but if chosen, the
   equivalence gates (#1/#4) become the SOLE guard against a shape error, so they MUST run on the
   LD multi-moment config (trailing `2^d` axis) where a shape slip is most likely.

**Recommendation (state in the carve crosswalk):** prefer path (a) — give M₀/R₀ real spaces that
reuse `frame.measure.space` (per-ordinate) and ONE minted ℓ=0-scalar space — so the OperatorSum
AND the inner OperatorProducts VALIDATE natively (illegal-states-unrepresentable, Pattern 4). Add
a **dedicated gate**: `require(isinstance(OperatorSum(S_iso, S_aniso), OperatorSum))` constructs
without raising (a positive-construction test), PLUS a negative test that a deliberately-mismatched
M₀ space (wrong name) RAISES `IncompatibleOperatorComposition` (the guard has teeth — anti-#11
positive+negative pair). Without the negative test, a `None`-skip would make the positive test
vacuous.

### R2 — [MEDIUM] M₀/R₀ vs `frame.measure_space`: the iso faces are NOT the frame's M/R

`S_aniso`'s M/R are the FULL spherical-harmonic faces (`_FrameAnalysis`/`_FrameReconstruction`,
domain `L2[S^2]/(N,)`, codomain `spherical_harmonic_space/(L+1,2L+1)`). The Option-2 `M₀`/`R₀`
are the ℓ=0 ROW only (scalar codomain `(ng,*sp)`, no `(L+1,2L+1)` tensor). They are **mathematically
the ℓ=0 restriction** of the frame faces (Gate 7d) but are SEPARATE operators (that is the whole
point — no moment tensor at ℓ=0). **Risk:** an implementer "reuses" `frame.analysis`/`frame.reconstruction`
restricted to ℓ=0 — which would re-introduce the `(L+1,2L+1)` allocation A2a reverted for. **Mitigation:**
Gate 8a (the no-moment-tensor sentinel) directly catches this — if M₀ is `frame.analysis[ℓ=0]`,
the tensor verb fires and the counter goes >0. M₀/R₀ MUST be the cheap scalar `Σ_n w_n ψ_n` /
broadcast, routing through `apply_p0_in_scatter` (the scalar verb), NOT `apply_legendre_scattering_moments`.

### R3 — [MEDIUM] R-double-`/W` (Pattern-7 crosswalk)

The `/W` lived in `_assemble_per_ordinate_source` (`(iso/W)+aniso`, with aniso already carrying
`/W` from `build_aniso_source`). Option-2 moves it OUTSIDE the sum: `S.apply = (1/W)·[S_iso+S_aniso]`.
**Two `/W` traps:** (i) does R₀ carry `/W` (double-counting the iso `/W`)? It must NOT — R₀ is the
bare broadcast, `/W` is the single producer-side factor at the `apply` boundary (Pattern-7/L18).
(ii) S_aniso = `kernel` is PRE-`/W` (the existing `build_aniso_source` applies `/W` separately) —
in Option-2 BOTH summands are pre-`/W` and the ONE `/W` is applied to the sum. **Mitigation:** Gate
#1 (oracle equiv at 1e-12) catches a double-`/W` directly (off by W); Gate #11 (`Q/Σ_t` flat-flux)
catches a `/W` convention error O(1). The implementer MUST write the per-ordinate↔iso `/W` crosswalk
into the carve plan BEFORE coding.

### R4 — [LOW] The HarmonicMomentFlux windowed arm

`_assemble_per_ordinate_source` is consumed by BOTH the `AngularFlux` arm (1540) and the
`HarmonicMomentFlux` windowed arm (1601). Option-2 must re-express BOTH (or keep the helper and
re-route it through `S_iso`/`S_aniso`). The windowed arm's iso is `phi_moments.scalar_flux()`
(the ℓ=0 moment, `Y₀⁰=1`) — it feeds the SAME K_iso. **Mitigation:** the existing windowed
regression + the LD gate (which exercises the moment-resolved path) cover it; confirm both arms
end at the same `(1/W)·[S_iso+S_aniso]` assembly (Pattern-2: one combine, two carriers).

### R5 — [LOW, but the Gate-3 trap] The LD monkeypatch target moves

The LD slope-source mutation gate (Gate 3) monkeypatches `_assemble_per_ordinate_source`. If
Option-2 RETIRES that method (folding `/W` into the `apply` boundary directly), the monkeypatch
silently no-ops → the gate passes vacuously (Mode-11). **The Gate-3 re-pin is the CRITICAL
deliverable:** determine where the iso slope-source `Σ_s·φ̂` flows in the NEW production (it flows
through `K_iso.apply` → `apply_p0_in_scatter`, which carries the `...` spatial-moment spectator),
RE-PIN the monkeypatch to that exact production line (e.g. `K_iso.apply` or
`MaterialXSField.apply_p0_in_scatter`), SENTINEL-instrument that the gate's run actually executes
it (in-process wrap counter >0), and mutation-verify it reddens (`|Δφ|/|φ| ≈ 2.6e-3 ≫ 1e-8`).
A green re-pinned gate that never executes the rewired line is vacuous.

---

## §ACTIVATION — the term-activation audit (Mode-7/Mode-10, AGENT.md §0.6)

The Mode-7 central risk: a P0-only fixture NEVER builds the ℓ≥1 blocks and is BLIND to the
iso+aniso SPLIT (it would test only S_iso). **EVERY equivalence/correctness/transpose gate
(#1,#2,#4,#6,#7c,#9,#10-het) runs on the ANISOTROPIC P1 + heterogeneous + `Sig2≠0` `solver_p1_het`
fixture** (already in both donor files — asymmetric `_P0_A/B`, `_P1_A/B`, `Sig2=[[0,.03],[.01,0]]`,
2 materials, 2G). The P0-only fixture is used ONLY for Gate 8 (the no-tensor perf claim) and the
iso-scalar arm. The **LD multi-moment** config (`build_2d_cartesian_ld_stress_mms_case`, trailing
`2^d` axis) is mandatory for #1/#2/#3 — the reverted A2a's regression was LD-only, so the oracle
gate MUST include an LD case (the spatial-moment spectator threads through `M₀`/`K_iso`/`R₀`).
`Sig2≠0` is built DIRECTLY on the Mixture (`_mix` in the donors), NOT via `make_mixture` (which
zeros `Sig2` AND `SigL` — lessons L1).

---

## §REBASELINE — disposition of EXISTING gates (explicit per §1.5)

| Existing gate | Option-2 disposition | Tolerance after carve |
|---|---|---|
| `test_scattering_adjoint.py::TestLambda/Kernel/N2N` (foundation) | **STAYS GREEN unchanged** (the leaves Option-2 reuses) | `rtol=1e-12..1e-13` (unchanged) |
| `test_scattering_adjoint.py::TestFullScatterKernel::test_reproduces_forward_scattering_source` | **EXTEND** to LD + add the transpose twin (Gate #1/#2) | `rtol=1e-12` (unchanged; +LD row) |
| `test_scattering_kernel_crosscheck.py` (ℓ≥1 `kernel`) | **STAYS 0-ULP** (S_aniso unchanged) — Gate #5 | `assert_array_equal` (unchanged) |
| `test_mms_ld_2d.py::test_ld_2d_scattering_slope_source_sign_mutation_reddens` | **RE-PIN** the monkeypatch target + Mode-11 sentinel — Gate #3 | `_CONSUMPTION_TOL=1e-8` (unchanged; new target) |
| `test_scattering_operator.py::TestBitIdenticalExtractionP0` (iso fast-path) | **KEEP or MIGRATE** — `add_iso_source`/`add_n2n_source` survive as K_iso's verbs (they are NOT retired; K_iso CALLS `apply_p0_in_scatter`/`apply_n2n`); if the public `add_iso_source` wrapper retires, migrate to `K_iso.apply` / the scalar verb | `rtol=1e-13`, teeth re-confirmed |
| `test_kinf_homogeneous.py` / `test_heterogeneous_transport.py` | **STAY GREEN unchanged** (forward-preservation, Gate #10) | gates' own thresholds |
| **NEW** Gate #2 (transpose oracle), #7a/b/d (M₀/R₀ laws), #8a (no-tensor sentinel), #12 (S† routing), R1 construction pair | **NET-NEW** | `rtol=1e-12` / counter>0 |

**Net-new test surface:** extend `test_scattering_adjoint.py` (Gate #1-LD, #2, #6, #7d, #8, #9-in-iso,
#12, R1 pair); add scalar transpose verbs `apply_p0_in_scatter_transpose` / `apply_n2n_transpose`
to `material_xs_field.py` (mirror the moment twins at :751/:840); the LD re-pin in
`test_mms_ld_2d.py`. The `S` dataclass gains `CAP_APPLY_TRANSPOSE` (today `{apply}` only) — Gate #12
FLIPS any "S has no transpose" pin.

**Catalog:** if Option-2 catches a real bug (the iso/aniso double-`/W`, the OperatorSum
space-name mismatch, the n2n-omission-from-Sᵀ, the LD spatial-moment slip), log an ERR-NNN per
the `vv-principles` "Log every caught bug" directive — the iso/aniso-split double-`/W` and the
OperatorSum-space-mismatch are new failure-mode candidates not yet catalogued.
