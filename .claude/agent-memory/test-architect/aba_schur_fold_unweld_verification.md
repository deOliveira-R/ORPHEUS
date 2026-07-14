---
name: aba-schur-fold-unweld-verification
description: The reusable recipe for gating a WELDED-fold un-weld (a per-cell/per-level closure hand-rolled at N sites → ONE single source) — worked on A_BA, the ψ½ Schur fold (Coupled Block Operator campaign Step 2). Single-source routing (Mode-11 wrap-counter), bit-identity inheritance, fold-transpose contract, iso-blindness refutation, non-carrying control.
metadata:
  type: reference
---

# Verifying a welded-fold un-weld — the A_BA ψ½ Schur fold (Step 2)

**When this recipe applies.** A per-cell/per-level closure algebra is
hand-rolled (inlined) at N production sites; a refactor extracts it to ONE
single source (a factory / a block operator) and routes the sites through it.
The un-weld is *bit-identity-preserving by intent* — the gate suite proves
(a) the sites REALLY route through the one source, (b) the output is
byte-identical to the documented old loop, (c) the transpose is single-sourced
too, and (d) the term the fold most likely gets wrong is ACTIVATED. Sibling of
`snapshot_migration_when_production_goes_bare.md` and lessons L16 (assembly
third-consumption-mode).

**The object (worked instance).** `A_BA` folds a bulk isotropic cell-emission
`q₀` onto the ray q½ source at μ=±1 per carrying radial level:
`Q̄(±1) = Σ_ℓ (2ℓ+1)/2 · Q_ℓ · P_ℓ(±1) = Σ_ℓ (2ℓ+1)/2 · Q_ℓ · (±1)^ℓ`.
Single source = `fold_moments_to_radial_characteristic(moments, sign)`
(ℓ-leading). At ℓ=0 → ½·Q₀ both signs (P₀≡1). Five hand-rolled sites: the fold
itself / `RadialCharacteristicSourceSink.from_angular_source` / S-forward /
F-forward / **S-adjoint (hard-codes `0.5`, bypasses the fold entirely)**. A_BA
lives OUTSIDE the resolvent (lagged S+F/k gain) → cannot touch the #284
forward-substitution certificate → LOW risk.

## The seven gate types (order = the refutations they encode)

1. **Fold contract on a MANUFACTURED anisotropic ≥2-moment input** (refutation
   #3, THE load-bearing one). The production S/F arms feed **ℓ=0 ONLY**
   (`emission[None]`), so an S/F-only gate is STRUCTURALLY blind to `P_ℓ(±1)`
   for ℓ≥1. Manufacture ℓ=0 AND ℓ=1 (≥2G, distinct per group); assert the
   closed form `Q̄(±1) = ½Q₀ ± (3/2)Q₁`. Tooth = a local `_fold_drop_sign`
   (drop `sign^ℓ`): anisotropic reds by `3·|Q₁|` (measured 2.70), **iso-only
   input stays 0.0** — the SAME assertion, proving the anisotropic input is
   NECESSARY (the iso-null asymmetry IS the config-blindness evidence, L1/L11).
2. **Single-source routing — Mode-11 wrap-counter (CENTERPIECE).** A counter on
   the SAME object the seed arms construct; assert BOTH S-forward AND F-forward
   enter it (a green-but-unrouted arm = a divergent inlined copy = Mode 11,
   counter 0). Count is EXACT `2·n_levels` (2 signs/level; sphere-GL S4 → 1
   level → 2). **The seed arms `local-import` the fold INSIDE `apply`**, so a
   `monkeypatch.setattr(source_module, "fold…", spy)` IS picked up (verified).
   TWO tiers: (2a LIVE) wrap the shared inner fold — proves both reach the
   shared Legendre math TODAY; (2b xfail) wrap the EXTRACTED single-source
   surface — the sharper un-weld proof, needs the landed bind.
3. **Bit-identity + closed-form.** `np.array_equal(seed, old_loop_reference)`
   (byte-identity vs the documented hand-rolled loop, survives the un-weld,
   inherits verification from gate 1) AND a **structurally-independent**
   `½·emission` closed form per (level,sign) + zero corners (independent of the
   fold fn — pins the coefficient, not just the routing). Tooth (external):
   forward fold → 0.6 reds the closed form (0.26 vs 0.22).
4. **F needs FISSILE, S needs scattering** (refutation #4). F emission
   `χ·νΣf·φ` is IDENTICALLY ZERO on a non-fissile mixture → a VACUOUS 0==0 gate.
   Split configs: S = non-fissile `Q/Σ` scatterer; F = fissile `A/2g` +
   `chi=[1,0,…]`. Mandatory **non-vacuity guard** (`max|emission| > 1e-6`)
   BEFORE the bit-identity assertion, else the fissile requirement is unenforced.
5. **The fold TRANSPOSE — TWO gates** (the S-adjoint hard-codes `0.5` =
   `coeff[0]`, must single-source through the fold's transpose). (5a
   helper-level, shape-agnostic, TODAY): the Euclidean adjoint contract
   `⟨fold(m,sign),y⟩ = ⟨m, fold_T(y,sign)⟩` with
   `fold_T(y,sign)[ℓ]=coeff[ℓ]·y`, on manufactured anisotropic `m`. Teeth: a
   `0.6` ℓ=0 coeff reds (0.02–0.04); a dropped `sign^ℓ` reds only sign=−1 (the
   P₁(−1) transpose consistency). (5b operator-level, TODAY):
   `⟨A_BA·φ, χ̄⟩ = ⟨φ, A_BAᵀ·χ̄⟩` via `S.apply`/`S.apply_transpose` (which
   exist regardless of the A_BA shape) + assert the adjoint's OUTPUT ray block
   is present-zero (`∂S/∂ψ½=0`). Tooth = monkeypatch the FORWARD fold→0.6 (the
   forward local-imports it; the adjoint hand-rolls 0.5, unaffected) → **the
   coefficients DISAGREE → reds (0.09)**. KEY: the consistency gate is blind to
   a SHARED coefficient value (both legs scale together) — it catches the
   fwd↔adj MISMATCH the hand-rolled 0.5 risks; the VALUE is pinned by gates 1/3.
   Survives the un-weld (fold & fold-transpose stay distinct entry points). If
   the chosen A_BA shape is factory-only (no transpose surface), 5a+5b carry the
   whole contract — do NOT force a `.apply_transpose` gate; flag the binding.
6. **Non-carrying CONTROL** (refutation #6, NOT "other geometries"). cylinder +
   slab have `radial_characteristic_space is None` → the seed arm emits a `None`
   ray and NO fold fires. Feed a bulk-only composite (arm gated on
   `psi.radial_characteristic is not None`); assert `None` ray + spy counter 0.
7. **Mode-8 throughout.** Every sentinel/tooth is `np.testing.*` / `pytest.fail`
   — a bare `assert` is a NO-OP under the canonical `-O`.

## The BIND isolation (shape decided in parallel)

`_apply_A_BA(emission, sn)` / `_apply_A_BA_transpose` / `_wrap_extracted_A_BA`
each carry a `# BIND:` marker line the main agent flips to the chosen shape
(factory `RadialCharacteristicSourceSink.from_moments` vs operator
`A_BA.apply`). Shape-agnostic contract gates (1,3,5a,5b,6) test the fold fn +
the real S/F operators directly and run TODAY; only the DIRECT-A_BA-surface
rows (2b + the two `_apply_A_BA` rows) are `xfail(strict=False, reason=…)` —
they flip to xpass when the bind lands (L4). Exercising A_BA *through* the real
`S.apply`/`F.apply` on hand-built seed-carrying composites (the seed arm fires
on any composite with a non-None ray, even though the production DRIVER is
"dormant until d3") is shape-stable — prefer it over the BIND for the
correctness gates.

## Result (this dispatch)

File `tests/sn/operators/test_psi_half_coupling.py::TestA_BA_SchurFold`
(`pytest.mark.foundation`); reuses `_sphere`/`_mixture`/`_random_composite`.
7 LIVE green + 3 xfail-skeleton under `-O`; every tooth mutation-verified
in-process (never `git checkout` an uncommitted file); pyright delta is the
file's pre-existing accepted `radial_characteristic_space | None` carrier-access
class (test file out of the ratchet). Did NOT commit (surgical mode — main agent
commits). Two flags to main agent: (i) F non-vacuity requires the fissile
mixture built locally; (ii) the pyright `| None` baseline is file-wide (40+
accesses) — a narrowing sweep is a separate file-level item, not this dispatch's.
