# The streaming path's objects say what they are, live where they belong, and carry only what they use

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
| `ReducedStreamingOperator` | `ChartConnection` *(or `StreamingCoefficients`)* | `[M]` **0 of 13** operator-surface members (`apply`/`domain`/`codomain`/`H`/`inverse`/`solve`/`apply_transpose`/`__matmul__`/`__add__`/`is_adjointable`/`block_role`/`system_role`); `SweepOperator` has 12 of 13. It maps nothing to anything. ⛔ And the name is **TAKEN** — the real reduced streaming operator is `L`, in the SN algebra, with genuine spaces. "Operator" is this codebase's most load-bearing word (the S4 amendment: *an operator is not an operator without its two spaces*); a struct wearing it teaches every reader the wrong thing. |
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
| `AngularRedistribution`, `angular_redistribution()`, `alpha_dome`, `_assert_alpha_dome_closes` | `orpheus/geometry/reduced_operator.py` | an **angular-discretization** package (`sn/angular/`) | `[M]` `alpha_dome(cosines, weights)` takes **no geometry argument at all**; the `1/r` lives in `ΔA` (the spatial factor), not in α. Chart-dependence is a *selection*, not spatial data. `[M]` the whole object is buildable from `(quad, coord)` with **no mesh** — and is already called that way on the multi-D path (`augmented_mesh.py:417`). |
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
| `μ_start`'s edge-extrapolation branch | a code path + the #361 `argsort` hazard | `[M]` **structurally unreachable** on the probed fixtures (perturbation `0.0`, control `4.04e-01`): its only consumer runs on NON-carrying levels, and every level carries. ⚠ `[M]` on 2 fixtures — **a sample.** | **Measure first** (§7 P1): evaluate `consumes_independent_seed` over every rule `assert_carrying_quadrature` admits. All-carrying ⟹ #361 is a **retirement**, not a repair. |
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
re-allocating. Measure `μ_start`'s branch over **every admitted rule** and then
retire or repair per the result. Add §6's honest refusal for a curvilinear
multi-moment metric.
**Done when.** `grep` finds no reader for each retired name; the µ_start
measurement is recorded with its denominator; the metric guard REFUSES a
constructible input (⚠ §6c — name that input in the gate, or the guard is not
gated); full fast set bit-identical.

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
**Deferred, deliberately.** Its third piece is `streaming_terms`' evaluated
per-`(cell, ordinate)` view, which is exactly the layer O-3 redesigns. Splitting
it here would churn that seam twice.

---

## 8. Explicitly NOT in scope

- **#158's cell solve.** `n_mom > 1` stays refused (`6859ca05`'s guard).
- **#407's carve.** Both instances are recorded; the fix belongs with O-3, where
  `cell_balance.py` reorganises under the scheme family anyway.
- **Re-baselining the `cyl_2g_3reg_folded_4x8_dd_n40` red.** `[M]` triply
  confirmed pre-existing and bit-identical across three independent carves
  (#404); it is `@pytest.mark.slow` and deselected by the canonical gate.
- **Choosing `G`.** §6 says what CANNOT adjudicate it; choosing it is its own
  campaign, pinned by physical functionals.

---

## 9. Forks the user must rule

1. **`ChartConnection` vs `StreamingCoefficients`** — the first names what it
   is (connection coefficients on a chart, the differential-geometry term the
   α-dome already cites); the second is more familiar and less right.
2. **`sn/angular/` vs `transport/angular/`.** `[M]` today the closure family has
   **zero** non-`sn` importers, so `sn/` is honest. But R19 requires the *scheme*
   family to serve diffusion; if the angular family is ever expected to serve a
   second method, `transport/` is the forward-compatible home. My read: the
   angular-cell closure is genuinely SN's (diffusion has no angular axis to
   close; MoC marches space; CP is integral) ⟹ `sn/angular/`.
3. **Does the chart-spelling unification (§1's last row) belong here or its own
   pass?** `[M]` 228 + 68 + 105 hits is a large mechanical sweep that touches
   files this plan otherwise never opens.
4. **P4's scope:** mint the package *whole* (spatial closure + pairing + moment
   measure), or only the pairing now and the rest at O-3?
