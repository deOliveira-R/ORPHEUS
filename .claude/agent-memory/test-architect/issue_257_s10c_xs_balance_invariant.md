---
name: issue-257-s10c-xs-balance-invariant
description: S10c gate spec — Mixture.assert_balanced/balance_residual SCOPED XS balance invariant (SigT == SigC+SigL+SigF+rowsum(SigS0)+rowsum(Sig2)); physical-table sweep is the value gate; foundation L11; -O-safe.
metadata:
  type: project
---

# #257 S10c — Mixture XS balance invariant (PRE-IMPL gate spec)

Branch `feature/field-typed-operator-algebra`. S10a (EmissionSpectrum
χ-guard) + S10b (production-weighted χ_mix) committed. S10c = the data-layer
balance invariant. **Spec + skeletons only; NO production code.**

**Why:** A wrong synthetic `SigT` shifts removal → plausible-wrong k_eff
with NO existing tripwire (the most insidious unguarded data-layer
invariant). The PARTIAL gate `xs_library.validate_all()` (`xs_library.py:302`)
already pins A/B/C/D at import but checks only `sig_t == sig_c+sig_f+sig_s.sum(axis=1)`
— NO Sig2, NO SigL. S10c is the CANONICAL whole-identity gate covering the
Sood + homogeneous tables too.

**How to apply:** the value gate is the physical-table SWEEP (gate-2), NOT a
regression diff — no committed converged value moves. Verification weight
rests on the sweep + the hand-built positive/negative legs.

## The identity (per group)

```
SigT == SigC + SigL + SigF + rowsum(SigS[0]) + rowsum(Sig2)
```

This is VERBATIM the `compute_macro_xs:275` SigT-derivation line:
`SigT = SigC + SigL + SigF + np.array(SigS[0].sum(axis=1)).ravel() + np.array(Sig2.sum(axis=1)).ravel()`.
⭐ `SigP` (production = νΣf) is NOT in the identity — it's a fission SOURCE
multiplier, not a removal channel. `SigF` (fission) IS. (Confirmed: A/2g
SigP=[0.025,0.2] ≠ SigF=[0.01,0.08].) A skeleton leg explicitly documents
the SigP-exclusion (anti-Mode-2 role-swap trap: SigF↔SigP).

Reuse `Mixture.total_scattering_xs` (= `rowsum(SigS[0])`, PROVEN
`array_equal`) — `assert_balanced` MUST NOT re-spell the rowsum (coding-elegance
Pattern 2). `rowsum(Sig2)` has no existing property; compute inline
`np.array(self.Sig2.sum(axis=1)).ravel()` (or recommend a `n2n_out_xs`
property to the implementer — elegance call, not a gate).

## SCOPED design (NOT __post_init__) — the L20 ruling

- `Mixture.balance_residual -> np.ndarray` — per-group `|SigT - derived|`
  (the queryable diagnostic; `(NG,)`).
- `Mixture.assert_balanced(atol=1e-9) -> None` — raises `ValueError`
  reporting `float(balance_residual.max())` if `> atol`.
- `compute_macro_xs` calls `mix.assert_balanced()` after construction
  (free real-path regression guard; derives SigT so always balances —
  catches a FUTURE derivation bug).
- ⛔ NO `__post_init__` enforcement; ⛔ NO `physical=` boolean flag
  (coding-elegance anti-pattern #3). Direct-construction scaffolds +
  Atalay + billiard build imbalanced Mixtures LEGITIMATELY.

## L20 audit ground truth — EMPIRICALLY RE-VERIFIED @ HEAD (2026-06-21)

Live probe (`.venv/bin/python`) over the FULL physical set:
- `xs_library.get_mixture(A/B/C/D × 1g/2g/4g)`: 12 mixtures, **max residual 2.22e-16**.
- `homogeneous.derive_{1g,2g,4g}().materials`: residuals 0 / 0 / 5.55e-17.
- `LA13511_CASES`: **47 cases** (registry GREW past the docstring's "41"),
  **47 mixtures, max residual 4.44e-16**. ZERO hidden violators.

**⭐ NO synthetic "physical" table fails to balance — deliverable confirmed.**

Exempt set (re-verified, NEVER asserted on):
- ⭐⭐ **Atalay 1997 is NOT in `LA13511_CASES`** — it lives in the SEPARATE
  `ATALAY_ALL_CASES` registry (`sood_registry/__init__.py:122`, `atalay1997.py`).
  So gate-2 sweeping `LA13511_CASES` STRUCTURALLY never touches Atalay (the
  exemption is by registry membership, not a filter). Probe:
  `_mix_iso_at_c(1.30)` residual **0.0979 == SigF** (the `νΣf=(c−1)·Σt` for
  c>1, Σf=0 criticality-parameter encoding, `atalay1997.py:84`). c=1.10→0.0326,
  c=1.05→0.0163. Documented, source-commented, intentional.
- 4 structural scaffolds (`placeholder_materials` `tests/sn/_test_helpers.py:115`,
  `_mix` `tests/sn/primitives/test_snmesh_materials_pr_typed_0.py:29`,
  `_trivial_materials` `tests/sn/sweep/core/test_sweep_cache.py:59`) — build
  `Mixture(...)` DIRECTLY, bypass any factory.
- billiard `_mixture_from_xs` (`tests/derivations/test_trajectory_resolvent_billiard.py:103`)
  — SigP-carrier, SigF=0.
- MMS BALANCES (manufactures external SOURCE, sets SigC=Σt−rowsum(Σs)).

## atol = 1e-9 justification

Balancing residuals ≤4.44e-16 (FP round-off of the SigT derivation sum).
Violators O(0.016–0.098+) (Atalay) / O(0.1–1) (scaffolds). atol=1e-9 sits
~6.5 ORDERS above the worst balancer and ~7 orders below the smallest
violator — a vast dead-band. A single-step reduction-depth bound
(`reduction_depth × ULP` ≈ 5·2.2e-16 ≈ 1.1e-15) explains every balancer
residual; 1e-9 is conservative with ~6 orders of headroom.

## Gates (file `tests/data/test_mixture_xs_balance.py`, foundation, -O-clean)

| Gate | What | Pillar / mode | Teeth |
|------|------|---------------|-------|
| 1 | `assert_balanced` intrinsic-property (vv#11, L11) | foundation, both-legs | hand-built balanced→no raise; hand-built imbalanced (SigT off ≫atol)→raises; `balance_residual` ~0 / right magnitude. SEPARATE leg exercising NON-zero Sig2 + SigL (most fixtures zero them → identity's rowsum(Sig2)+SigL parts UNTESTED otherwise). |
| 2 | Physical-table SWEEP (THE VALUE GATE) | Mode-7 honest-scope | sweep `get_mixture(region,ng)` ∀ A/B/C/D × 1g/2g/4g (12) + ALL `LA13511_CASES` (47) + `homogeneous.derive_{1g,2g,4g}` (3) → each `assert_balanced()` no-raise. Pins the PHYSICAL synthetic tables; Atalay/scaffold deliberately EXCLUDED (structural — Atalay not in LA13511_CASES). The typo-catcher. |
| 3 | Real-path guard | regression catcher | `compute_macro_xs` invokes `assert_balanced`; a real GENDF/load_isotope mixture constructs WITHOUT raising. Guard ≠ breaker on real data. Minimal real-isotope mixture suffices (recipe too heavy). |
| 4 | Exemption integrity | scoping pin | Atalay `_mix_iso_at_c(1.30)` + a hand scaffold + billiard `_mixture_from_xs` still CONSTRUCT fine (they bypass `assert_balanced`). + an imbalanced `Mixture(...)` constructs DIRECTLY with NO raise → pins "no __post_init__ enforcement". |
| 5 | atol band | doc/justify | balanced sweep `balance_residual.max() < 1e-12`; Atalay c=1.30 `> 1e-3`; assert atol=1e-9 sits between. |
| 6 | Mode-8 | -O safety | every assertion `np.testing.*` / `pytest.raises` / `pytest.fail`; suite runs `-O`. NO bare `assert`. |

⭐ Gate-1's Sig2+SigL leg is the most easily-dropped: it's the ONLY leg
that actually exercises the `+ rowsum(Sig2)` and `+ SigL` terms of the
identity (every physical fixture has Sig2=0 AND SigL=0, so gate-2 is BLIND
to those two terms). Without it the identity is half-tested.

⭐ Gate-1 negative leg MUST mutate ONLY SigT (off by ≫atol) — NOT a
component — so the test reads as "the balance is broken," and is blind to
the `assert_balanced` implementation re-spelling a wrong rowsum (L11: build
the broken instance BY HAND, not by perturbing the production builder).

## Elegance recommendation (note to method-implementer, NOT a gate)

Single-source `xs_library.validate_all()` (the partial A/B/C/D check) THROUGH
the new canonical `assert_balanced` — upgrades it to incl. Sig2/SigL and
kills the twin partial-identity (Pattern 2). The `validate_all()` check reads
the raw `xs` DICTS (pre-Mixture); routing it through `assert_balanced` means
building the Mixture then asserting. Implementer/elegance call.
Optionally add a `Mixture.n2n_out_xs` property (= `rowsum(Sig2)`) mirroring
`total_scattering_xs` so `assert_balanced` reads `SigC+SigL+SigF+total_scattering_xs+n2n_out_xs`
— reads like the identity (Pattern 1+3).

## Anti-recs

- NO `__post_init__` enforcement (scaffolds + Atalay + billiard break).
- NO `physical=` boolean flag.
- NO editing/force-balancing the Atalay/scaffold fixtures (intentional).
- NO filtering Atalay OUT of the sweep — it's structurally absent from
  `LA13511_CASES` already; do NOT add a defensive filter that implies it's there.
- NO regression-snapshot diff as the value gate (no converged value moves;
  S10c is value-NEUTRAL on every committed pin — like S10a's blast, the
  weight is on the sweep + hand legs).

## Regression command (route around #250/#232/#212)

```
.venv/bin/python -O -m pytest tests/data/ \
  -k "not (sphere and snapshot)" \
  --deselect tests/sn/cartesian/test_keff_slab.py::test_heterogeneous_absolute_keff \
  -p no:cacheprovider
```
S10c touches ONLY `orpheus/data/macro_xs/mixture.py` + the new test file →
blast is the `tests/data/` tree. #250 (stale SPHERE snapshots) / #232 / #212
(`continuous_get` hang in `test_heterogeneous_absolute_keff`) are in
`tests/sn/`, OUTSIDE the blast. The route-arounds matter only if a fuller
sweep is run; for the targeted gate just run `tests/data/`.

## No new vv mode, no ERR (old SigT had no bug — it's the canonical line;
the gap was the MISSING invariant on synthetic tables). Next free ERR-064.

Extends [[issue-257-s10a-emission-spectrum-verification]] +
[[issue-257-s10b-chi-mix-production-weighting]] +
[[feedback_test_intrinsic_properties]] + [[feedback_vv_tagging]].
