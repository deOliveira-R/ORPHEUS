# The streaming path's objects say what they are, live where they belong, and carry only what they use

> ## ▶ RESUME STATE — written at the 2026-08-26 compaction point
>
> **This file is the resume surface. Trust it over any summary.**
>
> ### Where the tree is
>
> Branch **`refactor/unweld-phase-b`**, **8 commits, NOT MERGED**; `main` is
> at `8f05ca14`. Reconcile with git before believing any of this
> (`process-discipline`: trust git, never a frozen claim).
>
> | commit | what |
> |---|---|
> | `27576937` | `AngularRedistribution` minted; six `Optional` fields retired off `ReducedStreamingOperator`; the fused `redist_dAw` cache retired |
> | `6adf6680` | the doc blast radius of that carve (16 sites) |
> | `250fcd16` | two theory chapters — the factorization, the τ arity theorem, the seed cone risk |
> | `226cc6ca` | `vv-principles` **#32** (a positivity review cannot see a wrong limit) |
> | `6859ca05` | **Phase B.2** — the closure takes `cls(angular, gram)`; the gram carries `(n_mom, n_thread)`; the `n_mom > 1` refusal |
> | `82853743` | this plan |
> | `1f1d2cce` | this plan's **§5b** (the three factors) |
> | `bd1404ea` | `vv-principles` **#17** sharpening + agent memories |
>
> ### ⚠ THE GATE WAS IN FLIGHT AT THE COMPACTION POINT — read it before merging
>
> The canonical fast gate (`python -O -m pytest -p no:randomly -m "not slow"`,
> serial, full tree) was launched at ~09:40 against `bd1404ea` and logs to
> **`scratch/_phaseb_merge_gate.log`**. At the compaction point it was at 6 %
> with zero failures. **Read that log's summary line before merging.**
>
> ⛔ **Do NOT quote the `4661 passed / 1 failed` figure in `6859ca05`'s commit
> message as this branch's gate.** `[M]` it finished at 03:52, BEFORE the last
> three fixes in that same commit (the write-only `_angular`, the dead local,
> two doc bugs). What covers those is the narrower `1058 passed` on
> `tests/sn/sweep/curvilinear` + `tests/geometry`. The in-flight gate is the
> first number that covers `bd1404ea`.
>
> ⚠ **The one known red is NOT ours and is NOT in this gate.**
> `cyl_2g_3reg_folded_4x8_dd_n40` fails at `0.1268` vs a `0.12` tolerance —
> `[M]` **bit-identical to 16 digits** at the pre-branch commit `8f05ca14`
> (baselined in an isolated worktree), which is the third independent carve to
> confirm it (#404, itself a **duplicate of #397** — flagged, not closed). It
> is `@pytest.mark.slow` and `[M]` DESELECTED by the canonical gate.
>
> ### What is DONE, and what this plan is for
>
> ✅ **Phase B of the un-weld arc is complete** (charter §5b O-2): the angular
> closure family no longer takes a mesh — it takes its two tensor factors. What
> this plan covers is the SUCCESSOR work the four audits surfaced: the names
> that are false, the homes the filtration contradicts, the welds still held,
> and (§5b) whether `L`'s three factors become first-class.
>
> ### ✅ ALL FOUR FORKS RULED 2026-08-26 (user) — P1 is unblocked
>
> Full text + the rejected alternatives' reasons at §9. The digest:
> **`ChartConnection`** · **`sn/angular/`** · the chart sweep SPLITS
> (`curvature` dies in P1/P3, the `geometry_kind` synonym is **P6**) ·
> P4 mints the pairing, the moment measure is **P7**.
>
> ⭐ **The governing ruling, which applies past these four:** *"scheduled to
> the end of this plan, **not defer and forget**."* Two phases exist because
> of it — **P6** and **P7** — and §8's "Choosing `G`" row was AMENDED rather
> than left standing. Any future "we'll do that later" in this campaign owes
> a phase number or a named external blocker (P5 has one: O-3).
>
> ⭐⭐ **P7 grew past this plan while being scheduled.** `[M]` 2026-08-26 the
> `FunctionSpace` metric is Hadamard **by construction** at every surface
> (`space.py:538/542/570/592`), so a non-diagonal `G` is unspellable — and
> **two** campaigns need one (this plan's §6; CS4b frame-square F-0's slab
> arm), with the debt recorded in three places and owned by nobody. P7 is
> that primitive, not "choose `G`" (which really is blocked on #158).
>
> ### The four memos this rests on are in `scratch/` and are NOT tracked
>
> `ld_curvilinear_shape_derivation.md`, `tau_under_ld_dip_analysis.md`,
> `adjoint_gram_ownership_audit.md`, `specialization_audit.md`, plus their
> probes (`probe_ld_*`, `probe_adjoint_gram_0{1..6}_*`,
> `_probe_beta_spatial_independence.py`, `_probe_cone_transmission.py`,
> `_probe_o2_operand_derivability.py`). ⚠ Untracked ⟹ **a `git clean` destroys
> them.** The theory chapters at `250fcd16` carry their conclusions; the memos
> carry the derivations.
>
> ### Issues opened or corrected this session
>
> **#407** DD's blend inverse in scheme-neutral `cell_balance.py` (+ a second
> instance in the seed march) · **#408** `is_positivity_preserving` conflates
> three properties — ⚠ its body's `(0, ½)` recommendation is **REFUTED in a
> comment** (that member is inconsistent: right rate, wrong limit) · **#404**
> the pre-existing cylinder red, now triply confirmed, duplicate of **#397** ·
> **#158** carries the curvilinear-LD literature + shape record.
>
> ### Durable rules landed this session
>
> `vv-principles` **#32** (never rank scheme candidates by positivity alone —
> add a consistency leg) · `vv-principles` **#17** sharpening (read a mutation's
> verdict by the red set's IDENTITY, not its size) · `lessons.md` **L60** (a
> type SIGNATURE is evidence about its author's assumptions, never about a
> theorem).


**Status.** Proposed 2026-08-26, un-ruled. Successor work to the un-weld arc's
Phase B (`27576937`, `6859ca05`); governed by
`.claude/plans/posing_filtration_charter.md` (the posing filtration, R1–R23)
and `.claude/rules/plan-authoring.md`.

**What this rests on** — three measured artefacts, all self-contained. Read
them before designing to any phase here; this plan CITES their numbers rather
than copying them (§9 of plan-authoring):

| memo | what it settles |
|---|---|
| `scratch/ld_curvilinear_shape_derivation.md` | the 1-D spherical LD moment balances; the three matrices `M` / `G` / `R` |
| `scratch/tau_under_ld_dip_analysis.md` | the tensor factorization; τ's arity theorem; the seed cone risk |
| `scratch/adjoint_gram_ownership_audit.md` | what shapes the adjoint; what `R` IS; the ownership table |
| `scratch/specialization_audit.md` | coordinate/scheme specialization; the weld list; both smuggling directions |

---

## 0. The three findings that reorganize everything

Everything below is a consequence of these. If one is refuted, re-derive the
phase that rests on it rather than patching the step.

**F1 — the redistribution operator FACTORS, and the factorization is the
ownership map.** `[M]` `𝓡 = R_spatial ⊗ A_angular(τ, α, w)` (tau memo §3).
`R`'s **measure comes from the chart**; its **basis comes from the scheme**;
its **angular partner is the measure's weights**. No single existing object
owns it, which is why it has been homeless and kept getting fused onto
whichever neighbour was convenient.

**F2 — `R` is the Galerkin matrix of `∇·ê_r`, and `R₀₀ = ΔA` is the
divergence theorem.** `[M]` adjoint audit A3:
`R_kj = ∫ b_k b_j (∇·ê_r) dV`, `∇·ê_r = (d−1)/r`, `∫∇·ê_r dV = ∮ê_r·n̂ dA`.
⭐ The consequence that matters for naming and homes: the *"geometry-dependent
factor-of-2 in the α normalization"* this tree has carried as a per-chart quirk
is an **artefact of spelling the measure `r^(d−2)`**. Written as `∇·ê_r` the
constant is `(d−1)` and **nothing per-chart remains** — so a great deal of
apparent sphere-vs-cylinder specialization is a spelling, not a fact.

**F3 — the closure STRADDLES TWO CACHE STRATA, and that is the weld nobody
named.** `[R]` from the audits, newly assembled here: the closure's `τ`, `α`,
`c_in`, `c_out` and the angular partition are **mesh-free** (quadrature ×
chart only — `[M]` bit-identical across 5 meshes incl. a graded one, Phase B's
probe), while its `_dAw_per_level` is **mesh-dependent** (`ΔA ⊗ 1/w`). Today
one object holds both, so a re-mesh rebuilds algebra that cannot have changed.
⟹ **§5c's ruled function/evaluated-table split applies to the ANGULAR closure
verbatim.** The charter already ruled it for the SPATIAL closure; the angular
one has the identical structure and the mesh-binding hid it.

---

## 1. Renames — every one of these names is currently false

⚠ Renames are `git mv` + import sweep + the THREE-search retirement audit
(code, tests, **docs** — a Python-domain xref renders as plain text with no
warning at any severity). Each row's evidence is the reason, not decoration.

| current | proposed | why — with evidence |
|---|---|---|
| `ReducedStreamingOperator` | ✅ **`ChartConnection`** *(ruled §9.1; `StreamingCoefficients` rejected — it lends `L`'s own word to a non-operator, reproducing the very defect being retired)* | `[M]` **0 of 13** operator-surface members (`apply`/`domain`/`codomain`/`H`/`inverse`/`solve`/`apply_transpose`/`__matmul__`/`__add__`/`is_adjointable`/`block_role`/`system_role`); `SweepOperator` has 12 of 13. It maps nothing to anything. ⛔ And the name is **TAKEN** — the real reduced streaming operator is `L`, in the SN algebra, with genuine spaces. "Operator" is this codebase's most load-bearing word (the S4 amendment: *an operator is not an operator without its two spaces*); a struct wearing it teaches every reader the wrong thing. |
| `redistribution_gram` | `redistribution_pairing` | `[M]` its own `(n_mom, n_thread)` axes admit the **rectangular** ONETRAN case (Hill 1975 Eq. 32, the angular index closed on the cell average only). A rectangular object is a **pairing, not a Gram** — the word over-claims off the diagonal of its own design space. ⚠ **This name is MINE, from `6859ca05`.** |
| `GeometryCoefficients` | `ChainScanCoefficients` *(or `SweepCoefficientCache`)* | `[M]` **0 of 15 fields** are un-permuted chart data: 4 are Morel–Montry constants, 3 angular, 2 pure traversal, and the closest 3 are **chain-ordered** (i.e. already traversal artefacts). The name promises geometry and delivers a permuted sweep cache. |
| module `sn/sweep/pole_angular_closure.py` | `…/angular/closure.py` | Two lies in one path. **"pole"** names the special case (the pole cell) for a family that closes the *whole* angular axis — `IdentityAngularClosure` never sees a pole. **"sweep"** is traversal (see §2). R15 already ruled the family `AngularClosure`. |
| `SNMesh.curvature: str \| None` | **retire** — use `coord` | `[M]` a **stringly-typed duplicate** of the `CoordSystem` enum (`"cylindrical"` / `"spherical"` / `None`), with `is_cartesian` defined as `curvature is None`. Read at two matvec entries through a **defaulted `getattr`** (`vv-principles` #28's temporal twin). ⚠ And `augmented_mesh.py:124` still documents it as storing `α_n/r_i` — a numeric quantity it has never held. |
| the chart's four spellings | one | `[M]` `coord` 228 hits / `curvature` 68 / `geometry_kind` 105 across `orpheus/`, with 8 / 3 / several case-spellings respectively. One concept, four vocabularies — a grep for any one returns a confident partial answer. |

⚠ **Do NOT rename `ReducedStreamingOperator` in isolation** — §3 shows it is a
bundle of three things at three stages, and a rename would name a bundle that
is about to stop existing. Sequenced accordingly in §7.

---

## 2. Re-homes — what the filtration says, versus where things sit

**The governing test**, from the charter's own ontology: an object's home is the
stage at which **all** its inputs exist. A module that must invent a Protocol to
avoid importing one of its own inputs is in the wrong stage.

| object | today | belongs | evidence |
|---|---|---|---|
| `AngularRedistribution`, `angular_redistribution()`, `alpha_dome`, `_assert_alpha_dome_closes` | `orpheus/geometry/reduced_operator.py` | ✅ **`sn/angular/`** (ruled §9.2; `transport/angular/` rejected — 0 candidate consumers) | `[M]` `alpha_dome(cosines, weights)` takes **no geometry argument at all**; the `1/r` lives in `ΔA` (the spatial factor), not in α. Chart-dependence is a *selection*, not spatial data. `[M]` the whole object is buildable from `(quad, coord)` with **no mesh** — and is already called that way on the multi-D path (`augmented_mesh.py:417`). |
| the `AngularMeasure` Protocol | same file | **dissolves** | Its own docstring gives the reason it exists: *"the geometry layer needs no import from the quadrature package at all."* `[M]` trace which members survive moving α out — **4 of 6 go with it**, and the other 2 exist only because `StreamingTerms` bundles angular data into a geometry packet (one of those two is dead, §3). ⟹ **a boundary Protocol that would shrink 6 → 0–2 is the shadow of a misplacement, not a boundary.** |
| the `AngularClosure` family | `sn/sweep/` | `sn/angular/` | `sn/sweep/` is the **traversal** package (`scan.py`, `cache.py`, `pairing.py`). The closure is a *discretization* object the sweep consumes. §5c's hard guard says traversal is contraband in the scheme; this is the exact mirror. And it breaks R15's symmetry — `SpatialClosure` lives at `transport/spatial/`, its sibling at `sn/sweep/`. |
| `R` / the pairing | a property on the geometry object | minted by the **scheme** (§5) | `[M]` today it lives on the geometry object and **never reads the scheme** — correct only because DD's basis is the constant. F1 says the basis is the scheme's. |
| `StreamingTerms` | `geometry/reduced_operator.py` | stays in geometry — but **shed its angular fields** | `[M]` it is **not pure geometry**: `mu`, `abs_mu`, `mu_start` and the `/w` inside `delta_A_over_w` are all quadrature reads. Those four are exactly what forces the 6-member Protocol. |
| `GeometryCoefficients` / `CollisionCache` | `sn/sweep/` | **stay** (home right, name wrong — §1) | `[M]` chain-ordered ⟹ genuinely traversal-time. This is the audit's D2 verdict: mostly NO smuggling, but the names lie. |

---

## 3. Un-welds — what is held versus what is used

| weld | what is held | what is used | verdict |
|---|---|---|---|
| `requires_upstream_angular_state`, `angular_marching_axis` | two fields on `ChartConnection` | `[M]` **0 production readers**; 12 test assertions. Made RED rather than grepped: flipping both on **997** constructed operators over 4 test trees gives **6 failures — exactly the 6 assertions that name the fields**. | ⛔ **RETIRE, not wire.** The concept is respelled twice already (`upstream_state.angular_upstream is None`, `SNMesh.is_cartesian`). |
| the three-link dead chain | `AngularRedistribution.mu_start_per_level` → `StreamingTerms.mu_start` → `GeometryCoefficients.mu_start` → **∅** | `[M]` the terminal has **zero readers** of any kind (dynamic access checked). So `StreamingTerms.mu_start`'s only production consumer is *the write into it* — while its docstring claims `MorelMontryAngularSweep` consumes it, which reads **the owner** instead. | Retire the two downstream links; the owner stays. |
| `ChartConnection._quadrature` | a second reference to the measure | `[M]` `_quadrature IS angular.quadrature` → **True**. A twin reference to one object. | Retire; read through the angular factor. |
| `μ_start`'s edge-extrapolation branch | a code path + the #361 `argsort` hazard | `[M]` **structurally unreachable** on the probed fixtures (perturbation `0.0`, control `4.04e-01`): its only consumer runs on NON-carrying levels, and every level carries. ⚠ `[M]` on 2 fixtures — **a sample.** | ⛔ **BOTH THE VERDICT AND THIS ROW'S DECISION RULE ARE REFUTED — measured 2026-08-26**, see below. **KEEP** the branch. |

> ⛔⛔ **REFUTED 2026-08-26 — and the refutation is of the METHOD, not the
> number.** The row instructed: *"evaluate `consumes_independent_seed` over
> every rule `assert_carrying_quadrature` admits. All-carrying ⟹ #361 is a
> **retirement**, not a repair."* Run exactly as written it returns
> **0 of 88** — true, and it licenses the wrong conclusion.
>
> `[M]` **`assert_carrying_quadrature` has exactly ONE call site**
> (`augmented_mesh.py:347`), inside `case CoordSystem.CYLINDRICAL`. The
> `SPHERICAL` arm (`:352-356`) calls **no admission gate at all**, and
> admission is by *structure, not provenance* (the gate's own docstring). So
> the prescribed denominator was **cylinder-only** while the branch's live
> witness is a **sphere** rule: a hand-built **Gauss–Lobatto** polar rule
> builds a production `SNMesh(SPHERICAL)` and **fires the branch** at
> `n = 6, 8, 9, 10, 11, 17` — `[M]` **6 of 11 orders**.
>
> ⭐ **And the naive check is TAUTOLOGICAL**, which the census caught and
> reported rather than banking: the gate and `_carrying_levels` read the
> *same producer* on the same `(quad, coord)`, so "an admitted rule with a
> non-carrying level" is a contradiction in terms. Its `0` carries **zero**
> information. What decided the question was execution: positive control
> (forcing `_carrying_levels := frozenset()` moves state `rel 4.246e-01` /
> `3.555e-01`) plus a gate-bypass separation — **carrying 105 built / 0
> fired; non-carrying 13 built / 13 fired.**
>
> ⭐⭐ **The #361 hazard is ANNIHILATED where it is reachable, by theorem.**
> `_edge_seed_stencil:1696` computes
> `t = (μ_start − μ[m0]) / (μ[m1] − μ[m0])`, and `on_edge_node` *means*
> `μ_start == μ[m0]` — so the numerator is exactly `0` and the ambiguous
> `m1` is multiplied by `0.0`. `[M]` over 75 reachable non-carrying levels:
> tie-break live in **0 of 75**, `t == 0` exactly in **74 of 75**;
> perturbing slot `m1` by `1e3` on the production-admitted Lobatto witness
> moves the seed by **`0.000000e+00`** (control on `m0`: `1.0e+03`).
>
> ⟹ **the disposition is a THIRD state the row did not offer** — neither a
> live repair nor a clean retirement. **Keep** the branch (retiring it
> removes the only seed path a `μ=−1`-noded sphere rule has, and that rule
> class is live project work — `rules_sphere.py:876` names "a future
> Gauss–Lobatto sphere" in `MarchStart.on_edge_node`'s own docstring); write
> the measured fraction into its docstring; and land any `argsort` change
> **with a Lobatto witness row in the same step** (§6c — otherwise the change
> ships with no gate that can redden).
>
> ⚠ Two riders. The 5 Lobatto *refusals* are **round-off-gated**
> (`tau_raw[N-1] = 1 + O(1e-15)` against an exact `[0,1]` guard), not
> structural — so the 6-of-11 is a floor, not a property. And `M` is closed
> over **shipped factories** only: Gauss–Radau (same `on_edge_node` witness,
> unprobed), Gauss–Chebyshev, `quotient(g≠Mirror('y'))` and arbitrary
> hand-built `Quadrature(measure=…)` are named holes, the last genuinely
> unclosable.
| `ChartConnection` itself | three stages in one object | pure chart (`face_areas`, `delta_A` — need **no** quadrature); the head-stage angular pairing (extracted at Phase B); and `streaming_terms`' evaluated per-`(cell, ordinate)` view (**O-3's layer**) | Split three ways. The Phase B carve already thinned it toward this; the two dead flags are the residue. |
| `redistribution_pairing` allocates per access | — | `[M]` `g1 is not g2` → **True**: `delta_A[:, None, None]` re-allocates on every call. Harmless today (one call per mesh) but it is a *cache-shaped* object with no cache. ⚠ **MINE, from `6859ca05`.** | Compute once at the mint (§5). |
| `#407` — DD's blend inverse in scheme-neutral code | `2.0 = 1/w_DD` at `cell_balance.py:248/:343` | `[M]` `cell_balance_terms` has ONE production consumer (`diamond.py:194`); `linear_discontinuous.py` never imports either helper; and the matvec **declares** the specificity the producer does not (`loss_representation:3091`: *"Curvilinear (DD-only) … with DD's diamond march inlined"*). | The whole module is a DD body in a scheme-neutral name and home. Carve with O-3. |
| `#407`'s second instance | three more `2.0 = 1/w_DD` at `psi_half_angle_seed.py:180-185` | Latent, guarded only by `_require_slab` **in a different package, on the other tensor factor**. | ⭐ Sharper than the first: that march is the **angular** factor of a product `pairing.py` declares orthogonal. |

---

## 4. Strata and laziness — F3 made operational

The tree already has the right pattern and it is the model to extend: the
two-strata cache (`sn/sweep/cache.py`) separates *geometry × quadrature*
(survives Σ_t rebinds) from *× Σ_t* (rebuilt once per epoch). The question this
plan answers is **which stratum each object in §1–§3 belongs to.**

| stratum | invalidated by | members |
|---|---|---|
| **S0 — mesh-free angular algebra** | a new quadrature or chart ONLY | `α` dome, the angular cell partition, `τ`, `c_in`, `c_out`, `μ_start`, the carrying-level set. `[M]` bit-identical across 5 meshes incl. graded. |
| **S1 — chart × basis** | a re-mesh, or a scheme change | `face_areas`, `delta_A`, the **pairing** `R`, the closure's `dAw = R ⊗ 1/w`, the moment measure `M/V` |
| **S2 — × Σ_t** | a cross-section rebind | `CollisionCache` (already correct) |
| **S3 — traversal** | a sweep-order change | chain ordering, leg decomposition, the pole-mirror permutation, inflow/outflow partition |

⟹ **The laziness rule this yields:** an object is built at the stage of its
*latest* input and cached against its *earliest* invalidator. Today the angular
closure is built once per **mesh** while `[M]` its S0 half cannot change under
a re-mesh — that is the F3 weld, and it is why the split is worth doing rather
than merely tidy.

⚠ **What NOT to make lazy.** `CollisionCache`'s eagerness is load-bearing —
`_build_count` instruments a cardinal invariance gate (*exactly once per Σ_t
epoch*). Do not convert a gate-instrumented eager build into a lazy one; the
gate measures the build.

---

## 5. Where `L` and the ray-characteristic operators MEET — the scheme's mint

**The question the user posed:** if two objects share something, and the shared
thing is discretization-related, that is where the shared concept should be
constructed.

`[M]` the sharers and the shared: `StreamingOperator` (`L`) and the
radial-characteristic family (`A_BB`, `Emission`, `Reconstruction`, `Seeding`)
**both consume the angular closure** — `L` through `cell_contribution` /
`angular_adjoint` in the matvec, the RC family at
`radial_characteristic.py:886/:894/:961`. They also both need the pairing,
transitively, because the closure carries it.

⟹ **The meeting point is the scheme's minted package** (R14/R19: the scheme is
a stage-2 generator on the Frame pattern), and the shape follows from F1 —
the scheme owns the basis, the chart owns the measure, the quadrature owns the
angular side:

```
package = scheme.mint(chart)              # S1: needs the mesh's chart
    .spatial_closure                       #   the cell relation (DD / LD)
    .redistribution_pairing                #   R  — (nx, n_mom, n_thread)
    .moment_measure                        #   M/V for the space's moment axis  ⚠ §6

angular_algebra = angular_redistribution(quad, coord)      # S0: mesh-free
angular_closure = device(angular_algebra,                  # device ∈ {WDD/M-M,
                         package.redistribution_pairing)   #   plain diamond, angular-LD}
```

and **both** `L` and the RC operators receive the minted `angular_closure` —
neither builds it, and neither reaches for a mesh to get it.

⭐ **Why the device is a separate axis from the scheme** (and why the scheme
must not simply *be* the angular closure's factory): `[M]` the rank
contradiction. Adams–Martin close the angular index **per spatial moment**
(square `R`); ONETRAN closes it on the cell **average only** (rank-1 column).
Those are two different *angular devices* over the same spatial scheme. Folding
the device into the scheme would make that choice unspellable — the very error
the `(n_mom, n_thread)` axes exist to prevent.

**Where the closure's own split lands (F3):** `angular_algebra` is S0 and
shareable across meshes; the `device(…)` binding is S1 and rebuilds with the
chart. That is the function/evaluated-table split, on the angular axis.

---

## 5b. ⭐ The three factors — make the algebra the implementation

**Added 2026-08-26 after a user question this plan had NOT anticipated.** The
adjoint audit's closing observation was *"signature of a product of one diagonal
factor and two triangular ones"*, and the question is whether those three become
first-class. This section says yes, states what must be verified first, and is
deliberately scoped as a SUCCESSOR campaign, not a phase of §7.

### The measured state — they are not explicit

`[M]` 2026-08-26:

* the **composition machinery already exists**: `OperatorProduct`
  (`numerics/operator.py:1650`) composes and is a `LinearOperator`, so `.H`
  propagates through it (`:983`); `PermutationOperator` (`:2458`),
  `ScaledOperator`, `PointwiseOperator`, `TensorProductOperator` and
  `_AdjointOperator` all ship;
* there is **no triangular-factor concept anywhere** — `grep` for
  `Triangular` / `forward_substitution` over `orpheus/` returns **one** hit,
  and it is `OperatorProduct`'s unrelated line;
* and the SN loss matvec is a **fused body** in `sn/loss_representation/`, not
  a composition — so today `apply` and `apply_transpose` are **two hand-written
  implementations of one relation**.

### Why this is the highest-value item the audits surfaced

`[M]` `apply_transpose ≡ apply.T` to `3.05e-16` / `1.99e-16` / `0.0`
(sphere / cylinder / slab). That gate exists **because they are twins** — it is
Pattern 2's signature, and the project rule is to make the twin unspellable
rather than to test that the two copies agree.

If the factorization is real and declared:

```
L  =  D  @  T_ang  @  T_spatial
L.H  =  T_spatial.H  @  T_ang.H  @  D          # D diagonal ⇒ self-adjoint
```

then **the adjoint stops being written and starts being derived.** The audit's
own inventory says exactly what each factor's `.H` must do, and it is small:
`[M]` every coefficient enters the adjoint **at the same value**, and the
reversal needs only **two order reversals and one block swap**. Those are
`T_ang.H` (the reversed march), `T_spatial.H` (the reversed scan), and the
trace-block swap. ⟹ `angular_adjoint` and the hand-written reverse scan RETIRE,
and the seam where a future angular device (ONETRAN's rank-1 column,
angular-LD) plugs in becomes ONE factor instead of two hand-mirrored bodies.

This is the *"build the machinery; realize the operator algebra ALWAYS"* ruling
applied to the last welded operation in the streaming path: a welded, un-named
composition is a failure to realize the algebra.

### ⛔ P0 — VERIFY THE FACTORIZATION BEFORE DESIGNING TO IT

⚠ **The three-factor claim is `[R]`, not `[M]`.** It was *inferred* from the
symmetry characters of `∂L/∂k` (`0` skew / `2` symmetric / `1.414`
triangular-like), which is evidence of triangular STRUCTURE and **not** a proof
that `L` equals a product of exactly three named factors. Designing to an
unverified factorization is how a campaign inherits a false premise.

**Done when** — all three answered, on the sphere (⛔ never the cylinder: `[M]`
`R/ΔA` is bit-exactly `diag(1, 1/3)` there and cannot discriminate):

1. `‖D·T_ang·T_spatial − L‖ / ‖L‖` at machine precision on ≥2 charts, with the
   factors constructed independently of `L`'s assembly;
2. **the factor COUNT is confirmed, not assumed** — the audit says "two order
   reversals AND one block swap", so the honest arity may be four
   (`D · P_mirror · T_ang · T_spatial`) once the pole-mirror permutation and
   the inflow/outflow trace partition are placed. Report the count you MEASURE;
3. each factor's triangularity is checked **in its own index order** (`T_ang`
   in the ordinate index at fixed cell; `T_spatial` in chain order), because a
   matrix is triangular only with respect to an ordering, and the ordering is
   itself a traversal artefact (§4, S3).

A refutation here is the most valuable outcome, not the least — it would mean
the adjoint is genuinely irreducible to a composition and the twin must stay,
which is a thing worth knowing before anyone tries.

### The constraint that shapes the design — algebra ≠ realization

⚠ The sweep is the **hot loop**, and `CumprodScan` (Blelloch) exists to make it
fast. Applying three factors sequentially through typed carriers would allocate
three times and destroy that fusion. **Do not read this section as "apply three
operators."**

The resolution is already this tree's own vocabulary, twice over: `WindowedSweep`
is *"the FUSED typed product `P @ A.inverse()`"* (a declared composition with a
fused evaluation), and §5c's ruled **function / evaluated-table** split says the
same thing. So:

* **declare** the factors — they are the algebra, and they carry `.H`;
* **realize** the composite's `apply` through the existing fused scan;
* and **gate** that the two agree.

⭐ That gate is the `3e-16` measurement, **repurposed**: today it compares two
hand-written twins (and so can only ever say "the copies match"); afterwards it
compares the DERIVED adjoint against the FUSED realization — an
algebra-vs-implementation equivalence, which is a strictly stronger claim from
the same number.

### Scope

This is **not** a phase of §7 — §7 is names, homes and welds, all bit-identical
or nearly so, and this is an architecture change to the operator algebra. It
belongs with the operator-strategy campaign / O-3's `StreamingOperator` binding
(`posing_filtration_charter.md` §5c), and it should follow §7 so the factors are
declared in objects that already have honest names and homes.

⚠ **Dependency, stated because it is easy to miss:** `T_ang`'s content is the
angular closure's S0 algebra (§4), so F3's function/evaluated-table split is a
PREREQUISITE, not a parallel. Splitting the closure after declaring `T_ang`
would cut the same seam twice.

---

## 6. The one CORRECTNESS item — not cosmetics

⛔ **A curvilinear cell is being given the slab's moment metric.** `[M]` an LD
spherical `SNMesh` **constructs today** and installs `spatial_moment weights =
[1, 0.3333]` (the slab mass) while the true `M/V` at the pole cell is
`[[1, 0.5], [0.5, 0.4]]`. It does not bite only because `loss_representation`
refuses one layer out (`IncompatibleRepresentation`, #158).

At `n_mom > 1` there are **three different matrices on one index** — `R/ΔA`,
the true `M/V`, and the shipped `diag(1, θ)` — and **only the third is a
metric.**

⚠ The hard part is not the value, it is that the true `M/V` is **cell-dependent
and non-diagonal**, which a per-axis `Axis(weights=…)` cannot express. `[M]` a
Gram–Schmidt basis diagonalises `M` and does **not** diagonalise `R`, so a
basis change cannot make both diagonal.

⟹ **This plan does not fix it** (it needs #158's cell solve to have a
consumer). It must **refuse honestly** instead: the metric installer should
reject a curvilinear multi-moment space rather than install the slab mass. §7
P1 carries that guard.

#### ⭐ Where the slab mass ENTERS — `[M]` located 2026-08-26

The installer is `SpatialScheme.moment_axis(ndim)` →
`Axis(…, weights=self.moment_mass_diagonal(ndim))`
(`transport/spatial/scheme.py:1375-1416`), and LD's override
(`linear_discontinuous.py:623-636`) builds it as

```python
base = np.diagonal(mass_1d(np.ones(()), self.theta))   # [1, θ]
```

⭐⭐ **`mass_1d` evaluated at UNIT WIDTH is the Cartesian-only shortcut, and
it is exactly where the slab assumption enters.** For a slab, `M/V` really is
width-independent — which is *why* the shortcut is valid and why nothing has
ever caught it. On a curvilinear chart `M/V` depends on `r_i`, so unit width
is not a normalisation, it is a different matrix. The docstring's claim
*"This IS `diag(M)/V`"* is true on exactly one chart family.

⛔ **And the guard §7 P1 promises cannot be written where the plan assumed.**
`[M]` **neither** `moment_axis(self, ndim: int)` **nor**
`moment_mass_diagonal(self, ndim: int)` receives a chart — so the installer
cannot refuse a curvilinear space *even in principle*: it is never told which
one it is on. ⟹ P1's item is **not** a guard insertion; it is a signature
change (the installer must be told the chart) plus the refusal, i.e. a §6b
call-site change. Sized accordingly in P1.
`[M]` the call-site set is small — **8 references total**, 5 in `orpheus/`
(2 defs, `scheme.py:1414`, and 2 docstring mentions) and 5 in `tests/`
(`test_spatial_moment_field_space.py:211`, `test_angular_bulk_space.py:217,262`,
plus 2 prose).

⭐ **Two expressiveness failures, not one** — worth stating separately because
they need different machinery, and #409 must cover both:
1. the true metric is **non-diagonal within the moment axis** ⟹ a matrix
   metric (#409);
2. it is **cell-dependent**, i.e. it varies along a *different* axis ⟹ it is
   not a per-axis weight at all, but a coupling between the cell axis and the
   moment axis (block-diagonal in cell, dense in moment).
   `Axis(weights=…)` fails **twice over**. #409's operator-valued framing does
   cover both — the metric is an operator on the product space — but a
   "matrix metric on the moment axis" framing would cover only (1).

⛔⛔ **And the instrument that cannot adjudicate it:** `[M]` reciprocity
(`A† ≡ G⁻¹AᵀG`) is an **identity for every invertible `G`** — `1.4e-16` under
Euclidean, random, and `(V·w)³`; mismatch control `8.22`. A reciprocity gate
proves *loadedness and consistency*, **never choice of metric**. What pins `G`
is the physical functionals (scalar flux, reaction rate). Do not cite a
reciprocity gate as evidence that a metric is right — here or anywhere.

---

## 7. Phasing — bit-identical first, behaviour last

Ordering principle: every phase that can be bit-identical **is**, and says so in
its done-when, so a value move is never hidden inside a move-and-rename.

### P1 — the dead stops being carried *(bit-identical; no renames, no moves)*
**Goal.** Nothing in the streaming path is held that nothing reads.
**Means.** Retire the two flags (with their 12 test assertions **migrated**, not
deleted — the successor spellings exist); retire the three-link dead chain's two
downstream links; retire `_quadrature`; cache the pairing instead of
re-allocating. ✅ **`μ_start`'s branch is MEASURED and the answer is KEEP**
(§3's row, refuted in place 2026-08-26) — so P1's item is no longer
"retire or repair" but: write the measured fraction into
`_edge_seed_stencil`'s docstring, and **retire the stale decision rule in
the `#361` comment at `:1683-1697`**, which still says *"this path may be
dead on the cylinder arm, in which case the answer is retirement"* — the
same too-narrow denominator this plan inherited from it.
⭐ Also fix the falsified docstring the census found: `precompute_psi_state`
(`pole_angular_closure.py:1707`) advertises `psi_view` as
`(N, ng, nx, 1) canonical`; `[M]` the production caller
(`loss_representation/__init__.py:3589`) passes **3-D** `(N, ng, nx)` and a
4-D array **raises** `too many values to unpack` at `:1752`. A falsified doc
is a bug (deferred only to the end of the in-flight merge gate, so the gate
covers one tree).
Add §6's honest refusal for a curvilinear
multi-moment metric — ⚠ **re-sized 2026-08-26**: `[M]` neither
`moment_axis(ndim)` nor `moment_mass_diagonal(ndim)` receives a chart, so this
is a **signature change + refusal** (a §6b call-site change over 8 references),
NOT the guard insertion this line originally implied. See §6's
"Where the slab mass ENTERS".
**Done when.** `grep` finds no reader for each retired name; the µ_start
measurement is recorded with its denominator; the metric guard REFUSES a
constructible input (⚠ §6c — name that input in the gate, or the guard is not
gated); full fast set bit-identical.

#### P1's §6b call-site set — `[M]` audited 2026-08-26, complete BY SYMBOL

⚠ "Complete **by symbol grep**" is the method, stated per §2's FILTER clause.
The 2026-08-24 surprise (a duck-typed kwarg surrogate + an attribute read that
no symbol grep returns) says the consumers' test run is part of this
enumeration, not a check performed after it.

| target | `orpheus/` | `tests/` | `docs/` sources |
|---|---|---|---|
| `requires_upstream_angular_state` | 7 | 6 | 4 — `methods/sn/index.rst:876,884`; `foundations/structured_geometry.rst:733,741` |
| `angular_marching_axis` | 6 | 6 | 2 — `foundations/structured_geometry.rst:734,742` |
| `_quadrature` (the FIELD) | **8** | **2** | 2 (both past-tense, about a *different* object) |

⚠ `index.rst:884` is a Python-domain `:attr:` xref — it renders as plain text
with **no warning at any severity**. Grep is the only gate on it.

⭐ **The audit changed this step: `_quadrature` is NOT dead.** The plan's §3
row called it "a twin reference to one object", which is true, and reads as
*delete it*. `[M]` it has **five live production reads** —
`reduced_operator.py:667`, `:686` (`mu_x[direction_idx]`), `:717`
(`level_indices`), `:719` (`eta[global_n]`), plus the narrowing assert at
`:650` — three writes (`:1065/:1137/:1243`), one declaration (`:565`), and
**two test reads** (`sn/sweep/core/_c_surrogate.py:163,167`, one of which is
an error-message string naming the attribute — grep its SHORTEST distinctive
fragment, not your own rewording).
⟹ the retirement is a **re-point of 7 reads** onto `self.angular.quadrature`,
not a deletion.

⭐⭐ **And that re-point is worth more than tidiness — it dissolves an
`Optional`.** `_quadrature: AngularMeasure | None` is optional; `AngularRedistribution.quadrature`
is **not**. So routing through the angular factor makes the `None` state
unrepresentable and **retires the narrowing `assert` at `:650`** with it —
Pattern 4, and one fewer bare `assert` in `orpheus/` (which `python -O`
strips anyway, so it was never a guard). This is the un-weld paying out: the
twin reference existed *because* the angular factor had no owner.

⚠ **Three false-positive families in the `_quadrature` grep — triage by
MEANING before treating any hit as a site** (raw substring: 16 code / 10 test,
i.e. **2× and 5× over-counted**):
1. every `*_quadrature` symbol (`chord_quadrature`, `select_quadrature`,
   `angular_quadrature`, …) — use `(^|[^a-zA-Z0-9])_quadrature\b`;
2. a **retired BC-shim** attribute of the same name (#176 / C176.1) —
   `test_snmesh_realizer_wiring.py:372,373,433`, `test_bound_compat.py:258,266`,
   and both `docs/` hits. Past-tense history about a *different object*: it
   **STAYS**;
3. a test-local helper *function* `_quadrature(coord)` —
   `test_kinf_homogeneous_tolerance.py:105,159,208`. Unrelated.

### P2 — the angular factor comes home *(pure `git mv` + imports; bit-identical)*
**Goal.** The angular half of the factorization lives in an angular package, and
the geometry layer stops declaring a Protocol to avoid its own input.
**Means.** `sn/angular/redistribution.py` ← `AngularRedistribution`,
`angular_redistribution`, `alpha_dome`, `_assert_alpha_dome_closes`.
`sn/angular/closure.py` ← the `AngularClosure` family (out of `sn/sweep/`).
Shed `StreamingTerms`' angular fields; the `AngularMeasure` Protocol shrinks or
dies.
**Done when.** `alpha_dome` has no geometry import; the Protocol's member count
is ≤2 or the file is gone; `sn/sweep/` contains only traversal; `dead_references`
0; bit-identical.

### P3 — the names stop lying *(mechanical; bit-identical)*
**Goal.** Every §1 row reads true.
**Means.** The §1 table, one rename per commit, each with the three-search audit.
`ChartConnection` LAST — after P4 has taken `streaming_terms` out of it, so the
name describes what remains.
**Done when.** Each old name has zero live references (past-tense history stays);
`sphinx -W` clean; `dead_references` 0.

#### `SNMesh.curvature`'s audit — `[M]` 2026-08-26, and it is NOT a rename

⭐ **ORDER IT FIRST, and the reason is a dependency the plan did not state.**
`is_cartesian` is the successor spelling P1 migrates the two dead flags ONTO
— and `[M]` `augmented_mesh.py:521` implements it as `return self.curvature
is None`. So the migration target is built on the field P3 retires. Retiring
`curvature` before re-basing `is_cartesian` breaks **12 production sites**
(`solver.py:894,1105`; `dsa.py:203`; `loss_representation` ×6;
`windowing.py`). §6b's shape, with a *definition* in place of a call site.

✅ **But the fix is small, because `is_cartesian` is the CONTRACT and
`curvature` is its implementation detail.** Re-base the one-line body onto
`coord` and all 12 consumers — plus both duck-typed `SimpleNamespace`
surrogates (`test_si_gate_dispatch.py:68`,
`test_unified_sweep_dispatch.py:354`) — never notice, because the surrogates
stub the **property**, not the field. ⚠ That is luck worth naming: those two
files are the *same family* as the 2026-08-24 kwarg-surrogate surprise; had
they stubbed `curvature=` they would have been invisible to every grep here.

**The real consumer set** (prose hits excluded — `[M]` "curvature" appears in
~50 more lines of docstring where it is the *physics word*, not the field):

| site | what | becomes |
|---|---|---|
| `augmented_mesh.py:348,355,1839` | the 3 writes | deleted with the field |
| `augmented_mesh.py:521` | `is_cartesian` body | `coord` — **do this first** |
| `augmented_mesh.py:879` | `if self.curvature is None` | `coord` |
| `loss_representation:712, 2727` | 2 repr / error strings | `coord` |
| `loss_representation:3129-3134, 3169, 3228, 3324` | matvec entry **1** | see below |
| `loss_representation:3479-3481, 3610, 3684, 3764` | matvec entry **2** | see below |
| 5 test reads | `test_radial_characteristic_slot_coordination:94`, `test_native_matvec:392`, `test_axis_native_construction:214`, `test_snmesh_consumes_reduced:192`, `test_sweep_regression:225` | migrate |

⭐⭐ **The two matvec entries are where the value is, and they are TWINS.**
Each opens with the identical three lines —

```python
curvature_raw = getattr(sn_mesh, "curvature", None)
curvature = curvature_raw if curvature_raw is not None else "cartesian"
if curvature not in ("spherical", "cylindrical", "cartesian"):
    raise ValueError(f"Unknown curvature: {curvature!r}")
```

— then dispatch on `curvature == "cartesian"` / `!= "cartesian"`. That is
**stringly-typed dispatch re-derived from a value that is already an enum**,
guarded by a **runtime re-validation of that enum's own domain**, written
**twice** (Pattern 2). Routing to `coord: CoordSystem` makes the `ValueError`
branch unrepresentable (Pattern 4) and collapses the twin. ⟹ this row is a
correctness-adjacent carve wearing a rename's clothes, and it should be
costed as such rather than batched with the cosmetic rows.

### P4 — the pairing is minted where its two halves meet *(behavioural at the seam)*
**Goal.** `R` is produced by the object that owns its basis, from the object that
owns its measure — and `L` and the RC family receive one minted closure.
**Means.** §5's `scheme.mint(chart)`; the angular device becomes the separate
axis it is; the closure splits S0 algebra from S1 binding (F3).
**Done when.** The scheme is the only producer of `R`; two meshes over one
`(quad, coord)` share one `angular_algebra` object by identity; `L` and the RC
operators receive the closure and neither constructs one; ⚠ a gate pins that the
rank-1 (ONETRAN) pairing is *expressible* — §6c: it must be constructible, not
merely refused.

### P5 — `ChartConnection`'s three stages separate *(rides O-3)*
**Deferred on a DEPENDENCY, not a punt.** Its third piece is `streaming_terms`'
evaluated per-`(cell, ordinate)` view, which is exactly the layer O-3 redesigns.
Splitting it here would churn that seam twice.
**Unblocked when:** O-3 lands its per-`(cell, ordinate)` layer. That is the
whole gate — checkable, and named here so this row cannot become a
defer-and-forget (the ruling at §9.3 applies to every "later" item in this
plan, and P5 is the one whose blocker is genuinely external).

### P6 — the chart speaks ONE vocabulary *(mechanical; bit-identical)*
**Goal.** One concept, one spelling. A grep for the chart returns the whole set
rather than a confident partial answer.
**Why it is a phase and not an issue.** §9.3: scheduled, not deferred. It is
last because it is a ~400-hit sweep across files P1–P4 never open — landing it
inside a phase whose done-when is bit-identity would dominate that phase's diff
and dilute its evidence. Its own phase, its own diff.
**Means.** Collapse `geometry_kind` (105) → `coord` (228). `curvature` (68) is
already gone at P1/P3. ⚠ Three-search audit per §1's warning, and grep the
**case-spellings** (`[M]` 8 / 3 / several) — a Python-domain xref renders as
plain text with no warning at any severity.
**Done when.** `grep -rn "geometry_kind" orpheus/` is empty; `coord` is the only
live spelling; `dead_references` 0; `sphinx -W` clean; full fast set
bit-identical.

### P7 — the space metric admits a NON-DIAGONAL `G` *(behavioural; the terminal phase)*
**Goal.** A function space can carry a metric that is not a Hadamard weight, so
a cell-dependent non-diagonal `M/V` is *expressible* — and §6's honest refusal
can become an honest value.

⭐ **This is NOT "choose `G` for curvilinear LD", and the difference is what
makes it schedulable.** Choosing `G` needs #158's cell solve to have a
consumer (§6) and is genuinely blocked. What is **not** blocked is the
machinery, and `[M]` 2026-08-26 the machinery is absent by construction:

| surface | `orpheus/numerics/space.py` | shape |
|---|---|---|
| `inner_product` (dense arm) | `float(np.sum(w * x * y))` `:542` | Hadamard |
| `inner_product` (axes arm) | `np.sum(_apply_axes_weights(x) * y)` `:538` | per-axis diagonal |
| `apply_metric` | docstring spells it `G⊙x` `:570` | Hadamard **by contract** |
| `apply_inverse_metric` | `(1/G)⊙x` `:592` | Hadamard |
| the one extension point | composite spaces override `apply_metric` per BLOCK `:577` | each block still diagonal |

⟹ there is **no** path for a non-diagonal `G`, at any level. `[M]` the
`__post_init__` guard at `:214` calling `inner_product_weights` "dense" means
*a dense array of diagonal weights*, not a matrix — which is exactly the sort
of near-miss vocabulary that makes an absent capability read as present.

⭐⭐ **Cardinal Rule 2 — TWO campaigns demand this one primitive**, which is
what promotes it from this plan's loose end to a shared foundation:

| campaign | the non-diagonal `G` it needs | what it does today |
|---|---|---|
| CS4b frame square, F-0 | the **slab** frame's live discrete Gram, off-diag ~1.15 | records `parseval unavailable` as a limitation; the Parseval gate SKIPS slab |
| this plan, §6 | curvilinear multi-moment `M/V` = `[[1, 0.5], [0.5, 0.4]]`, cell-dependent | installs the **slab mass** `[1, 0.3333]`; refused one layer out (#158) |

`[M]` recorded independently in three places as a "CS4c matrix-metric debt"
(`frame_square_recarve.md:122,325`; `collapse_pair_frame_induction.md:103`;
`orpheus-operator-machinery-report-v2.md:634`) — three records, one missing
capability, no owner. This phase gives it one.

**Means** (proposed, NOT verified — design at the phase, per §1's discipline):
the metric becomes a thing that is *applied* rather than a thing that is
*multiplied* — i.e. the Hadamard weight is the diagonal SPECIAL CASE of an
operator-valued metric, and `apply_metric`/`apply_inverse_metric` are its two
faces. That is the same shape CS4c already named ("the leg-as-operator
machinery"), so the two campaigns should converge here rather than each
minting one. ⚠ Check CS4c's state before designing — if it has landed legs,
this phase CONSUMES them and shrinks to the curvilinear installer.

**Done when.** A space carrying a non-diagonal `G` satisfies
`⟨x, y⟩ = yᵀGx` and `apply_inverse_metric ∘ apply_metric == id` on `G`'s
range; the CS4b slab arm's Parseval gate STOPS skipping (⚠ §6c — that
un-skip IS this phase's witness, and it must be named in the gate); §6's
curvilinear refusal is re-derived — either it becomes a value, or it states
which of the two blockers still stands. Euclidean and diagonal paths
bit-identical.

---

## 8. Explicitly NOT in scope

- **#158's cell solve.** `n_mom > 1` stays refused (`6859ca05`'s guard).
- **#407's carve.** Both instances are recorded; the fix belongs with O-3, where
  `cell_balance.py` reorganises under the scheme family anyway.
- **Re-baselining the `cyl_2g_3reg_folded_4x8_dd_n40` red.** `[M]` triply
  confirmed pre-existing and bit-identical across three independent carves
  (#404); it is `@pytest.mark.slow` and deselected by the canonical gate.
- **Choosing `G`** — the VALUE. §6 says what CANNOT adjudicate it; choosing it
  is its own campaign, pinned by physical functionals, and it needs #158's cell
  solve to have a consumer.
  ✅ **AMENDED 2026-08-26 (user ruling, §9.4) — the row was too wide.** It read
  as excluding the *machinery* along with the *value*, and those have different
  blockers: `[M]` the value is blocked on #158, the machinery is blocked on
  nothing and has **two** waiting consumers. The machinery is now **P7**, a
  scheduled terminal phase. Only the VALUE stays out of scope.
  > ⚠ This is plan-authoring §3's REMEDIED-fact case: the row was true when
  > written and is now half-false, and nothing would have prompted its edit
  > because being true is exactly why nobody re-reads it.

---

## 9. Forks — ✅ ALL FOUR RULED 2026-08-26 (user)

The fork text stays (plan-authoring §3: a decided premise is edited in place,
never dropped — the alternative's REASON is what stops it being re-litigated).

1. ✅ **`ChartConnection`** — the first names what it is (connection
   coefficients on a chart, the differential-geometry term the α-dome already
   cites); the second is more familiar and less right.
   ⛔ **`StreamingCoefficients` REJECTED**, and the reason is structural, not
   taste: *"streaming" is what `L` DOES*. Lending `L`'s vocabulary to a
   non-operator struct is the exact defect `ReducedStreamingOperator` is being
   retired for, so the alternative reproduces the bug it is meant to fix.
   ⭐ And the positive evidence is **F2**: after P1's retirements and P4's
   extraction the residue is `face_areas` + `delta_A`, and
   `ΔA = ∮ê_r·n̂ dA = ∫∇·ê_r dV` with `∇·ê_r = (d−1)/r` — the *contracted
   connection coefficient*. The residue is connection data, not metric data,
   so the name is literally right rather than merely evocative.
2. ✅ **`sn/angular/`.** `[M]` today the closure family has
   **zero** non-`sn` importers, so `sn/` is honest. But R19 requires the *scheme*
   family to serve diffusion; if the angular family is ever expected to serve a
   second method, `transport/` is the forward-compatible home. My read: the
   angular-cell closure is genuinely SN's (diffusion has no angular axis to
   close; MoC marches space; CP is integral) ⟹ `sn/angular/`.
   ⛔ **`transport/angular/` REJECTED** — it asserts cross-method generality
   with **0** candidate consumers. A later move is a pure `git mv`.
3. ✅ **SPLIT, and the remainder is SCHEDULED — not deferred.** Retire
   `SNMesh.curvature` inside P1/P3 (it is a stringly-typed *duplicate* of the
   `CoordSystem` enum read through a defaulted `getattr` — a live defect, and
   already load-bearing for P1's flag retirements). The `geometry_kind` ↔
   `coord` *synonym* unification is **P6**, a terminal phase of THIS plan.
   `[M]` 228 + 68 + 105 hits is a large mechanical sweep that touches
   files this plan otherwise never opens — which is why it gets its own phase
   and its own diff, **not** why it gets postponed out of the plan.
   > ⭐ The user's ruling, and it governs #4 too: *"we'll do option 1 to not
   > dilute the diff, but the rest of the work should be scheduled to the end
   > of this plan, **not defer and forget**."*
4. ✅ **Mint the pairing (+ spatial closure) at P4; the moment measure is
   P7, a terminal phase of THIS plan.** Not `.moment_measure` at P4 — that
   member requires choosing `G`, and §6 measured that no gate available here
   can adjudicate the choice. But per the ruling above it is **scheduled, not
   deferred**: see **P7**, whose scope is materially different from "choose
   `G`" and is unblocked (§7).
