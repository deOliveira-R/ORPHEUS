# `ReactionRateFunctional` + the `outer(Vector, Functional)` dyad — fission as a typed rank-1 reaction operator

> **Real home:** move to `.claude/plans/reaction_rate_functional_dyad.md` on the first execution turn
> (plan-mode pinned this session to `reactive-moseying-cake.md`, which held a stale COMMITTED copy of the
> #261 plan — the live #261 record is `.claude/plans/issue_261_collision_collapse_relocation.md`). Branch
> `refactor/operator-inverse-algebra`. **Surgical, main-agent-direct, user-steered; NO `method-implementer`.**
> Host env → `.venv/bin/python` (3.14). This is P4.5 task #81's resolution (supersedes "retire / wire-S6 / keep").
> Full verification spec: `.claude/plans/reactive-moseying-cake-agent-a0c71124a7a88d541.md` (test-architect).

---

## ✦ STATUS — Phases 0–3 + Phase 4 (type + SN) COMMITTED @ `bbd9c5b` (reconcile vs git, not this line)

Phases **0–3 (core dyad)** + **Phase 4 (the `IntegratedReactionRate` type + SN keff routing)** + **Phase 5
(archivist docs)** landed as ONE commit `bbd9c5b` on `refactor/operator-inverse-algebra` (27 files, +1738/−873).
Gates: 7-and-only-7 baseline reds (disjoint); pyright 412 (neutral); 0-ULP fission canary + the Mode-11
routing sentinel + the dyad / k∞ / IRR intrinsic gates; layer-imports; Sphinx `-W`. `ProductionRateFunctional`
RETIRED; the procedural twin REPLACED by the closed-form `k∞ = λ_max(A⁻¹F)` per-term oracle.

**DEFERRED (tracked GitHub issues):**
- **#270** — route CP/MoC/diffusion keff through `IntegratedReactionRate` (bespoke denominators: CP net-removal,
  diffusion leakage, MoC per-region; homogeneous 0-D = do NOT fold). SN (`sn/solver.py:1141`) is the template.
- **#269** — (n,2n) tests vacuously green (every `xs_library` mixture has `Sig2=0`); use `_make_mixture_with_n2n`.
- **#271** — `operator_algebra.rst` "C/S/F left at domain=None" reversed by W-D (+ `_integral_kernel_helpers.py:180`
  test-docstring `ProductionRateFunctional` ref).

**The PG-campaign hook:** `IntegratedReactionRate` is the **φ†=1 degenerate** of `⟨φ†, M[Σx]φ⟩` — the
homogenization PG campaign (P6, adjoint-weighted) generalises it by replacing the implicit `φ†=1` with a real
adjoint flux.

---

## Context — why

Resolving the `ProductionRateFunctional` fate (#261 P4) opened a deeper architecture, surfaced by the user and
confirmed by the cross-domain-attacker:

- **Fission is a rank-1 dyad** `F = |χ⟩⟨νΣf| = outer(χ, ⟨νΣf|)`. The row-factor `⟨νΣf|` (contract the flux
  against the production cross section over groups) **is a `Functional`**; the column `|χ⟩` is the emission
  reconstruction. `RankOneOperator.apply`'s `inner = (right*x).sum(axis=0, keepdims=True)` line
  (`operator.py:1811`) is **character-identical** to `ProductionRateFunctional.evaluate`
  (`production_rate_functional.py:151`) — i.e. the functional already lives *inside* fission's kernel, as a
  raw array. That duplication is the procedural twin (vv-principles L11) the campaign kept flagging.
- **`ProductionRateFunctional` was the right concept, wrongly seeded.** It is production-dead (sole consumer =
  its own host property `FissionOperator.production_rate`, read only by tests) and is a per-cell-DENSITY flavor
  with no live consumer — built while the *genuinely consumed* reaction-rate functionals (the seven
  `compute_keff` / `compute_group_production_rate` / `compute_group_absorption_rate` sites across SN/CP/MoC/
  diffusion/homogeneous) went **untyped** as bare methods.
- **The generic machinery the user named:** `outer(reconstruction: Vector, functional: Functional) → rank-1
  Operator`, with the matvec **routing through `functional.evaluate`** (the functional IS the contraction, not
  a parallel description). This *dissolves* the twin (one contraction, owned by the functional, used by the
  operator) and unifies fission's row-factor with the cross-method k-estimators into ONE type.
- **Polyadic framing (validation, not new machinery):** the angular scattering kernel is an *orthogonal* CP /
  spectral decomposition (`frame.conjugate(Λ) = Σ_ℓ |Y_ℓ⟩σ_ℓ⟨Y_ℓ|`, Funk–Hecke); the full scattering operator
  is a **block-term decomposition** (rank-1 angle ⊗ FULL energy transfer ⊗ diagonal space); **fission is its
  CP-rank-1 degenerate** (single term, ℓ=0, rank-1 energy). So `outer(Vector, Functional)` is the rank-1 atom
  and the `Frame` is the multi-mode manager — the same idea at rank-1 and rank-L. Guardrails: it is a **lens,
  not a type** (do NOT mint a `CPOperator`); keep the energy block dense (it is genuinely full-rank); do not
  import general-CP fitting (physics hands us the factors).

**Outcome:** fission reads `F = outer(χ, ReactionRateFunctional(νΣf))`; the twin is gone; the `Functional`
category is seated on real consumers at both layers; the cross-method `compute_keff` family folds onto one
typed source of truth; and the weak procedural-twin oracle is **replaced** by a structurally-independent
per-term `k∞` oracle — turning the carve from tidiness into a correctness upgrade.

## User decisions (this session)

- **Generalize + promote, do NOT retire the concept.** `ProductionRateFunctional` (per-cell density) retires
  *as a symbol*; its math is absorbed into the generalized `ReactionRateFunctional(Σx)`.
- **Mint the L1 generic functional (the "full" option):** numerics `InnerProductFunctional(weight, axis)` =
  `⟨weight,·⟩`; transport `ReactionRateFunctional(cross_section)` **specializes** it (mirrors
  `HarmonicFrame(GalerkinFrame)`). The category gets an L1 *and* an L2 instance.
- **The dyad routes its matvec through `functional.evaluate`** (settled — it is the whole point). The raw
  `RankOneOperator(left, right, axis)` ndarray ctor **retires**.
- Volume-measure seam: **deferred** to the Phase 4 `compute_keff` fold (lean: functional stays density-pure,
  a separate `∫·dV` spatial measure does the integral).

## Audit findings that shape the plan (explorer + test-architect, code-verified)

- **`RankOneOperator` is FISSION-ONLY in production** — one constructor (`fission.py:306`). Scattering builds
  `frame.conjugate(Λ)` (`OperatorProduct`), **never** a `RankOneOperator`. ⇒ the "scattering 0-ULP canary" is
  a **phantom**; do NOT claim a scattering bit-identity gate (Mode-11 vacuous). The real regression walls are
  **(α)** the fission crosscheck `test_fission_kernel_crosscheck.py` B.1/B.2 and **(β)** the L0 primitive suite
  `tests/numerics/test_tensor_product_operator.py::TestRankOneOperator`.
- **Layer-CLEAN by construction.** `RankOneOperator` + `Functional` + `InnerProductFunctional` are L1
  (`numerics`); `ReactionRateFunctional` + `CrossSectionField` are L2 (`transport`). Typing the dyad's
  row-factor as the **L1 `Functional` Protocol** is L1→L1; the L2 concrete satisfies it structurally. **Keep
  `outer`/`RankOneOperator` parametrized on the L1 Protocols — never reach for `CrossSectionField`/
  `ReactionRateFunctional` concretely in `numerics`** (that is the one way to introduce an L1→L2 violation).
- **The load-bearing edit:** `RankOneOperator.apply` re-homes its `inner` contraction onto
  `functional.evaluate(x)`; the `left.shape == right.shape` guard (`__init__:1782`) goes (the functional owns
  its own contraction; the reconstruction is a separate vector). Bit-identical by construction (same numpy
  primitive, now behind the functional).
- **Two estimator surfaces:** Surface A (`KEigenvalue` operator-triple `Callable` estimators) is **unused in
  production**; Surface B (per-solver bound `compute_*` via `power_iteration` → `solver.compute_keff`) is the
  real path. Phase 4 folds **Surface B**, NOT A.
- **Two verification traps:** (1) Mixture-A **2g** has a coincidentally FLAT eigenspectrum → flat-flux ratio
  equals `k∞` there → flux-shape-blind; the per-term pin **MUST use the converged φ\*** and **4g is
  MANDATORY**. (2) The functional is density-only; `compute_keff` folds `∫dV` — on a unit-volume mesh they
  coincide and a volume bug hides (Mode-10) → a **heterogeneous, non-unit-volume row** must isolate it.

---

## Target architecture (reads-like-the-math)

```
# L1 numerics — the generic co-vector and the dyad
InnerProductFunctional(weight, axis=0)            # ⟨weight, ·⟩ ;  evaluate(x) = (weight*x).sum(axis, keepdims=True)
outer(reconstruction: Vector, functional: Functional) -> RankOneOperator
RankOneOperator.apply(x) == reconstruction * functional.evaluate(x)      # routes through evaluate

# L2 transport — the domain-typed reaction-rate co-vector
class ReactionRateFunctional(InnerProductFunctional):   # carries a CrossSectionField (.values, .mesh)
    # production = ReactionRateFunctional(νΣf) ;  absorption = ReactionRateFunctional(Σa)

# fission, as a typed rank-1 reaction operator
F.kernel == outer(χ, ReactionRateFunctional(νΣf)) & IdentityOperator()
```

`RankOneOperator` keeps its name (named for its defining property, rank = 1); `outer(v, w)` is the readable
builder. Subclass-vs-compose for `ReactionRateFunctional(InnerProductFunctional)` finalizes at implementation
(lean subclass, per `HarmonicFrame(GalerkinFrame)`); the base stores `weight = cross_section.values, axis=0`
and the derived keeps `cross_section` for the domain handle (`.mesh`, units). Built fresh per access (the
`mat_xs` read-through, as today).

---

## Phases (each: serial `-O` green + CLI-pyright ≤ baseline + Sphinx `-W` + layer-imports gate)

### Phase 0 — Verification scaffold + safety net (test-architect spec)
- Lock the canary baselines: fission B.1 (`rtol=1e-13` hand-loop) + B.2 (`χ·evaluate ≡ apply`, 0-ULP) and the
  L0 `TestRankOneOperator` suite. Capture the closed-form `k∞` reference (computed, not remembered):

  | mixture | ⟨νΣf,φ\*⟩ | ⟨Σa,φ\*⟩ | k∞ |
  |---|---|---|---|
  | A 2g | 0.15909903 | 0.08485281 | 1.8750000000 |
  | A 4g (**mandatory**) | 0.19041882 | 0.12799012 | 1.4877619048 |

- Reference dot products written as **explicit Python loops, NOT `νΣf @ φ`** (a `@` re-shares numpy's reduction
  with the functional's `.sum`, re-introducing the procedural twin in the oracle).
- Mode-8 audit (`grep '^\s*assert '` on the new gates — the suite runs `-O`).

### Phase 1 — Mint the Functional types (additive, no behavior change)
- **L1** `orpheus/numerics/functional_impls.py` (or alongside `functional.py`): `InnerProductFunctional(weight,
  axis=0)`. Intrinsic-property gate (the defining law `evaluate(x) == (weight*x).sum(axis)`; linearity).
- **L2** `orpheus/transport/reaction_rate_functional.py`: `ReactionRateFunctional(cross_section: CrossSectionField)`
  specializing `InnerProductFunctional` (generalizes `ProductionRateFunctional`; νΣf → arbitrary Σx).
- **Headline correctness gate** — `tests/transport/test_reaction_rate_functional.py`: pin
  `ReactionRateFunctional(νΣf)` AND `ReactionRateFunctional(Σa)` **each independently** against the closed-form
  `k∞` at the **converged φ\*** (A 2g + A 4g + a heterogeneous non-unit-volume row). Pinning each (not just the
  ratio) catches a shared-factor error the ratio masks. `ProductionRateFunctional` still exists (retired P3).

### Phase 2 — The `outer(Vector, Functional)` dyad + fission rewire (canary-pinned surgery)
- Refactor `RankOneOperator` (`operator.py:1694`) → `(reconstruction: Vector, functional: Functional)`; `apply`
  routes through `functional.evaluate`; **retire the raw `(left, right, axis)` ctor** + the `left.shape ==
  right.shape` guard. Add `outer(v, w)` builder. Migrate the L0 `TestRankOneOperator` (9 sites) onto
  `outer(v, InnerProductFunctional(w, axis))` — absorbed into the intrinsic-property gate.
- Rewire `FissionOperator.kernel` (`fission.py:305`) → `outer(χ, ReactionRateFunctional(νΣf)) & IdentityOperator()`.
- **Gates:** fission B.1/B.2 **0-ULP** preserved; **Mode-11 sentinel** = in-process `monkeypatch` call-counter
  on `ReactionRateFunctional.evaluate`, assert `count > 0` after `F.apply` (proves the matvec ROUTES through the
  functional — a green twin routing around it leaves the counter at 0). pyright ≤ baseline; layer-imports.

### Phase 3 — Retire `ProductionRateFunctional` (the twin dissolves)
- Delete `production_rate_functional.py` + `FissionOperator.production_rate` (the per-cell-density concept now
  lives ON as the fission row-factor `ReactionRateFunctional(νΣf)`).
- **Replace the procedural-twin oracle**: the B.2 `χ·evaluate ≡ apply` row is downgraded to a demarcated
  equivalence/Mode-11 row; correctness moves to B.1 + the Phase-1 closed-form per-term oracle.
- Rewire the single-point-of-change probes `tests/transport/_functional_helpers.py`,
  `_integral_kernel_helpers.py` → `ReactionRateFunctional`; repoint the category-partition tests
  (`test_functional_category.py` POSITIVE arm, `test_integral_kernel_category.py` foil). `hand_derived_fission_emission`
  kept verbatim (structurally independent).

### Phase 4 — Cross-method `compute_keff` fold (Surface B; LATER / optional-depth)
- Fold the per-solver `compute_group_production_rate` / `compute_group_absorption_rate` / `compute_keff` (SN
  `solver.py:1070/1120/1158`, CP `661/689`, MoC `211`, diffusion `297`, homogeneous `111`) onto typed
  `ReactionRateFunctional(νΣf)` / `ReactionRateFunctional(Σa)` — one source of truth for `⟨Σx, φ⟩`.
- **Resolve the volume-measure seam:** lean density-pure functional + a separate `∫·dV` spatial measure (the SN
  fission term also carries a `(n,2n)` `2·Σ2` contribution — preserve it). Each solver volume-weights
  differently ⇒ **its own sub-phases**, gated hard (touches live `keff` math). This is the "category earns its
  keep on live consumers" payoff; depth decided with the user when Phases 0–3 land.

### Phase 5 — Docs / theory (archivist)
- Re-narrate `docs/theory/operator_algebra.rst` + `docs/api/numerics.rst` (~9 `RankOneOperator` anchors): the
  dyad `|v⟩⟨w| = outer(Vector, Functional)`; the **polyadic/BTD framing** (fission = rank-1 term; scattering =
  Frame-managed orthogonal CP sum; collision = diagonal; energy block stays dense); `ReactionRateFunctional`
  as the energy-axis co-vector; the `Functional` category seated at L1 (`InnerProductFunctional`) + L2.

**Core = Phases 0–3** (the dyad + `ReactionRateFunctional` + fission, twin dissolved). **Phase 4** is the
unification payoff (own depth decision). **Phase 5** docs.

---

## Mutation probes (Phase 0/2/3 — all in-process `monkeypatch`; NEVER `git checkout` on uncommitted files)
M1 scale νΣf → production pin RED · M2 permute χ → reconstruction RED · M3 swap production↔absorption (O(1)
move 0.159↔0.085) → ratio RED · M4 mutate contraction axis → dyad-action RED · M5 bypass `evaluate` with an
inline reduction → Mode-11 sentinel counter-0 RED.

## Critical files
- `orpheus/numerics/operator.py` (`RankOneOperator` → `(reconstruction, functional)`; `outer` builder; retire raw ctor)
- `orpheus/numerics/functional_impls.py` *(new)* — `InnerProductFunctional`
- `orpheus/transport/reaction_rate_functional.py` *(new)* — `ReactionRateFunctional(InnerProductFunctional)`
- `orpheus/transport/operators/fission.py` (kernel → `outer(χ, ReactionRateFunctional(νΣf))`; drop `production_rate`)
- `orpheus/transport/production_rate_functional.py` *(retire)*
- Phase 4: `orpheus/{sn,cp,moc,diffusion,homogeneous}/solver.py` (+ `moc/core.py`) `compute_*` methods
- tests: `tests/transport/test_reaction_rate_functional.py` *(new headline oracle)*;
  `tests/numerics/test_tensor_product_operator.py::TestRankOneOperator` (rewire → `outer`/`InnerProductFunctional`);
  `tests/sn/operators/test_fission_kernel_crosscheck.py` (B.2 → demarcated row + sentinel);
  `tests/transport/{_functional_helpers,_integral_kernel_helpers,test_functional_category,test_integral_kernel_category}.py`
  (rewire → `ReactionRateFunctional`); `test_production_rate_functional.py` (rewire/retire)
- `docs/theory/operator_algebra.rst`, `docs/api/numerics.rst` (archivist, Phase 5)

## Verification (per phase + end-to-end)
- `.venv/bin/python -O -m pytest tests/numerics tests/sn tests/transport -p no:xdist --timeout=300 -p no:cacheprovider -q -rfE`
  — SERIAL (xdist unstable); baseline = the **7-and-only-7 pre-existing reds**. Canaries: **fission B.1/B.2
  0-ULP + L0 `TestRankOneOperator`** (NO scattering canary — phantom).
- `npx --no-install pyright --outputjson orpheus/` — the ORACLE (not the `<new-diagnostics>` LSP, #226);
  baseline **412**, target ≤ (the typed dyad should not regress).
- `tests/test_layer_imports.py` — numerics must not import transport; the L1 `Functional` Protocol is the
  decoupling device.
- Sphinx `-W` clean throughout.
- **Proactive `test-architect` (done — spec banked) + `archivist` (Phase 5).**

## Scope / discipline
- Surgical, main-agent-direct; **NO `method-implementer`**; per-phase ff-merge; `main` always green. Commit only
  when the user asks; stage files explicitly (no `git add -A`); trailer
  `Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>`. NEVER `git checkout`/`restore`/`stash`
  on uncommitted files (monkeypatch to revert mutation probes). User rejects `# type: ignore` (`cast` OK).
- **Banked architecture:** a rank-1 reaction operator = `outer(reconstruction: Vector, functional: Functional)`,
  matvec through `evaluate`; `Functional` is the operator's row-factor (co-vector), `Frame` its multi-mode
  generalization; fission = CP-rank-1 atom, scattering = Frame-managed orthogonal sum (BTD), collision =
  diagonal. `ReactionRateFunctional(Σx)` is the one type for fission's row-factor AND the cross-method
  k-estimators. **No `CPOperator`/`ReactionOperator` umbrella type** (lens, not a type; L-004).
