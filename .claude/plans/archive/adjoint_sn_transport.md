# Adjoint SN transport solve (φ\*) — the exact discrete transpose campaign

> **APPROVED 2026-06-28** (user). Its OWN campaign (user chose "adjoint solve as its own
> campaign" over folding it into P6): own GitHub issue + own theory page; land + verify φ\*
> standalone, THEN #51/P6-B (adjoint-weighted homogenization) lands as a thin consumer.
> Branch `feature/sn-adjoint-transport` (from `main` @ `68ceb9a`). Surgical, main-agent-direct,
> user-steered; **NO `method-implementer`** (operator-algebra carve). Per-phase ff-merge;
> `main` always green. Construction route **CONFIRMED = exact discrete transpose** (user).
>
> Verification spec (test-architect, the gate chain): `.claude/plans/archive/p6_adjoint_verification_spec.md`.
> Structural ruling (cross-domain-attacker): `.claude/agent-memory/cross-domain-attacker/sn_adjoint_construction_frames.md`.
> GitHub issue: **#276** (the campaign home).
>
> ## STATUS (reconciled 2026-06-29 — read this first; reconcile vs git)
> **THE 2026-06-28 PROSE BELOW IS STALE.** A1 (F†) + A2 (S†, P3 `15185e5`, **closes #118**) are DONE
> via the **#276 K_iso reframe** — durable plan `.claude/plans/archive/kiso_cross_model_extraction.md`. The
> old "Option 2" (a lazy ℓ=0 frame to make the iso forward composable-AND-fast) was **SUPERSEDED**:
> `K_iso` (`IsotropicScattering`+`IsotropicN2N`) achieves composable-and-fast directly — the forward
> keeps the scalar fast-path through K_iso (P2 `dcea43a`, 0-ULP) and the adjoint Sᵀ rides
> `full_scatter_kernel` (P3 `15185e5`). The asymmetric forward-fast/adjoint-frame split an earlier
> note said the user "rejected" is in fact what LANDED, and is principled (forward hot-path + validated
> frame transpose, reciprocity-cross-checked). Also done en route: meshless homogeneous on `(C−K_iso)`,
> 3 eigenvalue engines, EnergyGrid geometry (P4), P5 docs (`9cf70d9`). **NEXT = A3** (`solve_transpose`
> + `.H` CAP_SOLVE) — phase spec ~line 130; verification `.claude/plans/archive/a3_solve_transpose_verification.md`
> (test-architect). The A3–A6 phase specs + per-leaf sizing below remain valid.
>
> --- stale 2026-06-28 prose follows ---
> Branch `feature/sn-adjoint-transport` from `main` @ `68ceb9a`. **3 commits landed + green:**
> - `bbd4479` **A1 — F†** (fission adjoint via the rank-1 dual-dyad transpose on `RankOneOperator`).
> - `e8d9f6b` **A2 foundation — Λᵀ** (per-ℓ group-transpose; the aniso `kernel.apply_transpose` is now LIVE, reciprocity rel 4.5e-15).
> - `f3f3dd2` **N2N** (`N2NMomentOperator` — (n,2n) as a distinct in-frame ℓ=0 operator; `apply_n2n_moments`/`_transpose` verbs).
>
> **A2a forward-migration ATTEMPTED, then user PIVOTED to Option 2 (below).** I routed the
> scattering apply arms through `full_scatter_kernel = frame.conjugate(Λ_{ℓ≥0} + N2N)` (validated:
> reproduces the forward S to **rel 1.9e-16** for non-LD; full reciprocity rel 1.5e-16). The
> forward-preservation gate surfaced THREE findings: (1) a **latent LD IndexError** in
> `apply_legendre_scattering_moments` (the `(Ellipsis, *idx)` cell-indexer mis-targets the trailing
> spatial-moment φ̂ axis when LD's `2^d` per-cell moments are present — never exposed because LD only
> ever used the P0 fast-path) — **FIXED** (→ `(slice(None), slice(None), *idx)`, bit-identical
> non-LD); (2) the LD slope-source **mutation gate lost its teeth** with the frame forward
> (`test_ld_2d_scattering_slope_source_sign_mutation_reddens`) — the frame path DIVERGES from the
> fast-path for LD φ̂-slope scattering (open Q: real frame bug vs the mutation patching fast-path rows
> the frame bypasses); (3) **PERF** — the frame conjugation builds the moment tensor + M/R matmuls
> every SI iteration, whereas the iso fast-path (`add_iso` = reaction-rate group-matmul on the scalar
> flux, NO moment tensor) is a genuine hot-path optimization. **The forward migration was REVERTED**
> (forward back to fast-path, green). The `test_ld_2d_converges_second_order_smoke` "timeout" was a
> RED HERRING — it is `@pytest.mark.slow` (deselected by the canonical `-m "not slow"`), it timed out
> only because I ran the file unfiltered; NOT a regression.
>
> **UNCOMMITTED (to be committed as the Option-2 foundation):** the `full_scatter_kernel` property
> (the validated frame form — the A2b transpose form + Option-2 forward reference) + the cells LD-fix.
> Forward arms are code-identical to `f3f3dd2` (docstring diffs only).
>
> ## NEXT (Option 2, user-chosen 2026-06-28 — pursue in a CLEAN context)
> Make the iso forward **composable (frame) AND fast AND LD-correct** — do NOT accept the asymmetric
> "forward fast-path + transpose frame" compromise (user rejected it). The crux: the ℓ=0 frame
> conjugation `R₀∘Σ∘M₀` IS mathematically integrate→group→broadcast (the fast-path) since `Y₀⁰=1` —
> so a frame form that **compiles to the fast-path at ℓ=0** (no moment-tensor materialization for the
> isotropic part) would be both composable and fast. Approaches to explore:
> 1. A **perf-specialized / lazy frame** that does not materialize the ℓ=0 moment tensor — recognize
>    `M₀`=weighted-integrate, `R₀`=broadcast and fuse them (the fast-path) while keeping the composable
>    `apply_transpose` for free. (The ℓ≥1 aniso already builds the tensor; only ℓ=0 needs the fast path.)
> 2. **Resolve the LD slope divergence**: determine whether the frame path is genuinely wrong for LD
>    φ̂-slope scattering, or the mutation gate just targets the fast-path rows. Pin the frame-LD
>    equivalence with an LD-multi-moment gate (my 1.9e-16 probe was non-LD single-moment only).
> 3. Then **A2b** (`S.apply_transpose = (1/W)·full_scatter_kernel.apply_transpose`) lands SYMMETRICALLY
>    (forward + transpose both via the fast frame). Then A3 (`solve_transpose`) → A4 (daggered posing,
>    φ\*) → A5 (φ\* carrier) → A6 (docs).
>
> **The validated facts are pinned in COMMITTED tests** (the throwaway probes were ephemeral
> `$CLAUDE_JOB_DIR/tmp` — gone post-compaction): `tests/sn/operators/test_scattering_adjoint.py` —
> `TestFullScatterKernel` (the property `op.full_scatter_kernel` reproduces the forward S to
> rel<1e-12 + Euclidean reciprocity), `TestN2NMomentOperator` (N2N reciprocity + ℓ=0-only),
> `TestLambdaTranspose` (Λᵀ identity + dense reference + group-flip), `TestKernelTranspose` (aniso
> kernel reciprocity). The committed `full_scatter_kernel` property IS the frame form to perf-specialize.

## Why this campaign exists

P6 (eigenvalue-consistent, adjoint-weighted homogenization, #51) needs an **adjoint flux φ\***
that does not exist on `main`. φ\* is independently high-value (perturbation theory, sensitivity /
GPT, point-kinetics parameters, importance), so it is built as a first-class, reusable
adjoint-transport capability — not a homogenization-only shim. `eigenvalue.py:26-29` already
names the adjoint posing row `A_loss† ψ† = λ M† ψ†` as a "future seam"; this campaign activates it.

## The construction route — exact DISCRETE transpose (not μ-reversal, not Krylov-on-transpose)

**`k_adjoint = k_forward` is exact only for the transpose of the *discretized* operator** (same
characteristic polynomial). The μ→−Ω *continuous* adjoint diverges (~1e-3 on curvilinear WDD,
where the diamond closure is not μ-symmetric) — it would break the headline identity. The discrete
route also reuses the existing reverse-walk sweep (single source of truth); μ-reversal would be a
second path to `Lᵀ` (elegance smell). The adjoint is a **daggered posing row, not a layer** (L23b):
transpose the triple, feed the *unchanged* `power_iteration`.

**Literature:** Lewis & Miller (1984) §6.4 — the SN adjoint is the forward operator with `Ω→−Ω` +
transposed scattering/fission kernels; the *discrete* adjoint is what makes `k_adj = k_fwd` exact.

**Per-leaf sizing (most is already built):**

| Piece | Verdict | New code |
|---|---|---|
| `(L+C)†` matvec | already ships — `loss_action_transpose` (reverse walk), verified by `test_g_adjoint_reciprocity` | none |
| **`(L+C)†⁻¹` inner solve** | **the ONE genuine gap** — reverse-DAG forward-substitution `solve_transpose`; `operator.py:650-651` defers exactly this | **one sweep primitive** |
| `B†` | self-adjoint under `|Ω·n|·w`; `(L+C−B).H` already reachable | none (algebra) |
| **`F†`** | **free dyad-swap** — `F = outer(χ, ReactionRateFunctional(νΣf))` (fission.py:302); transpose = `outer(νΣf, InnerProductFunctional(χ, axis=0))` = `|νΣf⟩⟨χ|`. `RankOneOperator` docstring already promises this | none (constructor re-wire) + `CAP_APPLY_TRANSPOSE` |
| `S†` | `R∘Λᵀ∘M` — angular M/R are an adjoint pair (Funk–Hecke self-adjoint); only the per-ℓ group-transfer Λ transposes | transposed kernel input + `CAP_APPLY_TRANSPOSE` |
| `.H` G-conjugation | `_AdjointOperator` already wraps `A† = G⁻¹AᵀG` | none |

**Keep the μ-reversal physical adjoint ONLY as a uniform-slab verification oracle** (the
fuller-view-oracle exception) — never the keff-bearing path.

## Phases (each: green + CLI-pyright ≤ baseline + ff-merge; surgical, main-agent-direct)

- **A0 — verification spec (DONE).** test-architect gate chain at `.claude/plans/archive/p6_adjoint_verification_spec.md`.
- **A1 — `F†` (free dyad-swap).** Add `apply_transpose` to `FissionOperator` via the dual dyad
  `outer(νΣf, InnerProductFunctional(χ, axis=0))` (reuse `RankOneOperator`'s own transpose); flip
  `capabilities` to add `CAP_APPLY_TRANSPOSE`; retire the fission.py:96-100 "deferred adjoint" note.
  GATE: **P1.0a** (F† vs explicit loop, χ ≁ νΣf, 2G+4G; mutation swap χ↔νΣf → RED — the canonical
  adjoint-fission trap). FORWARD-SAFETY: 0-ULP `test_fission_kernel_crosscheck` + `test_reaction_rate_functional`.
- **A2 — `S†` via iso-path modernization (CLEAN-BEFORE-EXTEND; user-steered 2026-06-28).** The aniso
  path is already composable (`(1/W)·frame.conjugate(Λ)` = `R∘Λ_{ℓ≥1}∘M`) — that is WHY its transpose
  is free. The iso path (P0 + n2n) is the legacy fast-path (`integrate→group-matmul→/W broadcast`),
  never routed through the frame — so its transpose was NOT free. The fix is to MODERNIZE the iso path,
  not bolt a transpose onto legacy. **Verified empirically:** because ORPHEUS uses unnormalized
  harmonics (`Y₀⁰=1`), `integrate_angular` IS the ℓ=0 analysis `M₀` and the iso broadcast IS the ℓ=0
  reconstruction `R₀`, so the ℓ=0 frame conjugation reproduces the fast-path P0 to **rel 2.2e-16**
  (probe `$CLAUDE_JOB_DIR/tmp/probe_iso_frame_equiv.py`).
  - **A2a (forward cleanup):** route the full Legendre scattering (P0 + aniso) through
    `(1/W)·frame.conjugate(Λ_{skip_l0=False})`. **n2n** stays a DISTINCT reaction (multiplication, not
    scattering — keep it visible) but joins as a separate ℓ=0 moment-space operator **summed with Λ
    before one `R∘(·)∘M`** (one M / one R, no double moment-build, physics preserved). Retire the
    bespoke `add_iso_source`/`add_n2n_source`/`build_aniso_source` assembly. GATES (test-architect spec
    `.claude/plans/archive/a2_iso_modernization_verification.md`): the **0-ULP aniso canary STAYS 0-ULP**
    (ℓ≥1 untouched); a SEPARATE principled-equiv (~1e-14) iso-through-frame gate; keff
    (`test_kinf_homogeneous`) + het-transport + SI≡Krylov **green unchanged**; perf non-regression on
    the P0-only path (L=0 moment tensor = the scalar flux).
  - **A2b (extend — transpose free):** `LegendreMomentScattering.apply_transpose` (per-ℓ group-flip
    `mfc,fg→` ⇒ `mfc,gf→`) + `CAP_APPLY_TRANSPOSE`; then `S.apply_transpose = (1/W)·full_kernel.apply_transpose`
    (covers iso+aniso+n2n in ONE transpose). GATE **P1.0b** mirrors the A1 fission suite
    (`test_fission_adjoint.py`): Sᵀ vs dense `Sᵀ` oracle built from forward applies; **Euclidean**
    reciprocity `⟨Sφ,ψ*⟩=⟨φ,Sᵀψ*⟩` (NOT the metric `.H` — test-architect correction), per-group not
    weight-summed (L27), asymmetric SigS + non-vacuous `Sig2≠0`; group-flip + omit-n2n mutations RED;
    Mode-11 leaf-wrap of `Λ.apply_transpose`. The bare transpose feeds `.H` (the metric G-adjoint, A3/A4).
- **A3 — `solve_transpose` (the one new primitive) + wire `.H` to carry `solve`.** The reverse-DAG
  forward-substitution inner adjoint solve on `(L+C)` (the analytic inverse of `(L+C)†`, the backward
  semigroup — NOT Krylov-on-transpose). Close the `operator.py:650-651` deferral so `(L+C).H`
  advertises `CAP_SOLVE` via `solve_transpose` ⟹ the daggered operator is `.solve`-identical to the
  forward for the eigenvalue driver. GATE: **P1.1** full-loss reciprocity `⟨A_loss ψ, φ⟩_G =
  ⟨ψ, A_loss† φ⟩_G` (per-group, rel<1e-12; extends `test_g_adjoint_reciprocity` from `(L+C−B)` to the
  full `A_loss = L+C−S−B` now that S† exists).
- **A4 — daggered posing row + adjoint k-eigenvalue (φ\* is computed here).** Pose
  `A_loss† = (L+C)† − S† − B†`, `M† = F†`; feed the UNCHANGED `power_iteration` (via `KEigenvalue`
  on the daggered triple — the inner solve is A3's `solve_transpose`). GATES: **P1.2**
  `⟨ψ*,q⟩=⟨ψ,Σ_d⟩` fixed-source reciprocity (asymmetric 2G, source/detector in different groups;
  donor = `tests/sn/solve/test_fixed_source_g1.py`); **P1.3** adjoint-keff = forward-keff (exact;
  ∞-medium **+ slab + sphere** — reject ∞-medium-only, flat flux hides a μ-reversal error); **P1.4**
  analytic ∞-medium adjoint spectrum (left eigenvector of `A⁻¹F` = right eigenvector of `Mᵀ`; extend
  `derivations/common/eigenvalue.py` `kinf_and_spectrum_homogeneous` to also return `eig(M.T)`; 4G —
  the sharpest F† test); **P1.5** bi-orthogonality `⟨ψ*_i, F ψ_j⟩ = 0` (i≠j; intrinsic-law test).
- **A5 — the φ\* carrier: `AdjointFlux` / `Solution.adjoint()`.** Give φ\* a typed home so consumers
  (homogenization, perturbation theory) read it cleanly. Decide carrier shape (a parallel adjoint
  `Solution`, or an `adjoint_flux` field / `AdjointFlux` type) — mint a TYPE only per the type-vs-property
  criterion (≥2 non-isomorphic realizations + a non-identity morphism). GATE: carrier intrinsic-property
  test + the `Solution.adjoint()` API smoke.
- **A6 — docs/theory (own theory page) + V&V.** New theory page (or major section): adjoint SN
  transport — the discrete-vs-continuous adjoint distinction, the daggered posing, the F† dyad-swap,
  the Funk–Hecke S† self-adjointness, the reciprocity/keff identities. Archivist. Sphinx `-W`; V&V
  matrix rows for P1.0–P1.5. Cross-link #166 (reciprocity invariants).

## Downstream — #51 / P6-B (the thin consumer; frame-projection campaign, NOT this campaign)

Once φ\* lands: **B1** `IntegratedReactionRate.evaluate(φ*, φ)` bilinear (generalize the φ†=1
degenerate `∫⟨Σx,φ⟩dV` → `∫ φ\* Σ φ dV`); **B2** `Solution.homogenize` test weight φ→φ\* (the
adjoint-weighted PG frame — swap at solution.py:428; decide forward-vs-adjoint API: parameter vs
sibling method); **B3** gates C1 (PG-vs-Galerkin discriminator where φ\*≠φ materially, + dud-guard),
C2 (keff-preservation under coarse-resolve — **first-order-stationary, error O(δφ²), a convergence-
ORDER claim**, wires the documented-but-untested `sn-homogenization-bilinear`), C3 (Mode-11 sentinel:
weight == φ\* not φ); **B4** docs (the bilinear homogenization section already drafted at
`discrete_ordinates.rst:15880-15937` — wire the gate, mark `sn-homogenization-bilinear` tested).

## Critical files
- **F† / S†:** `orpheus/transport/operators/fission.py` (dyad-swap + capabilities), `orpheus/transport/operators/scattering.py` (`R∘Λᵀ∘M` + capabilities); `orpheus/numerics/operator.py` (`outer`/`RankOneOperator` transpose; `_AdjointOperator` `.H`, the `solve`-deferral at 650-651).
- **solve_transpose:** `orpheus/sn/loss_representation/__init__.py` (`loss_action_transpose` ≈ 2602 — the reverse-walk matvec; add its `solve` twin).
- **Posing/eigenvalue:** `orpheus/numerics/eigenvalue.py` (the documented adjoint seam, 26-29), `orpheus/numerics/iteration.py` (`KEigenvalue`).
- **φ\* carrier:** `orpheus/sn/solution.py` (`Solution` + `adjoint()`); `orpheus/transport/fields/` (a possible `AdjointFlux`).
- **References:** `orpheus/numerics/functional.py` (`InnerProductFunctional` for F†'s χ-row), `orpheus/derivations/common/eigenvalue.py` (extend for the P1.4 left-eigenvector oracle).
- **Forward-safety twins:** `tests/sn/operators/test_fission_kernel_crosscheck.py`, `tests/sn/operators/test_scattering_kernel_crosscheck.py` (both 0-ULP).

## Scope / discipline
- Surgical, main-agent-direct, user-steered; **NO `method-implementer`**. Per-phase ff-merge; `main`
  always green. Branch `feature/sn-adjoint-transport` from `main`.
- **pyright should net DOWN** (retiring the F†/S† "deferred adjoint" confessions + honest transpose
  types) — the campaign closes #226-adjacent debt, not adds it. Baseline 412.
- Commit/push/merge ONLY when the user asks; stage explicitly (NO `git add -A` — `.claude/*`
  policy-uncommitted + gitignored `.h5`); trailer `Co-Authored-By: Claude Opus 4.8 (1M context)
  <noreply@anthropic.com>`. No `# type: ignore` (`cast` OK). NEVER `git checkout/restore/stash` on
  uncommitted files (L28); `git reset` (index-only) is safe.
- **Re-baseline, NOT bit-identity** where the discrete adjoint introduces new values; the
  reciprocity + exact-keff identities are the structural completeness proofs.
- Canonical gate: `.venv/bin/python -O -m pytest <paths> -m "not slow" -q -rfE --timeout=300 -p no:xdist -p no:cacheprovider`. Sphinx `.venv/bin/python -m sphinx -b html -E -W docs docs/_build/html`.
