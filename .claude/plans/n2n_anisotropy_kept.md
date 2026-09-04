# The (n,2n) channel keeps its angular distribution — one scattering-type family, lossless in ℓ

Plan of record for GitHub **#426**'s remedy. Opened 2026-09-03 on branch `fix/n2n-anisotropy`
(from `main` `1e02f6b1`). Governed by `.claude/rules/plan-authoring.md`. Campaign-2 context:
`.claude/plans/cs4c_binding_design.md` §18.6 step 1 — the user ruled #426 and #428 are fixed
BEFORE CS4c step 5 resumes; this plan is that step 1. Ground memos (untracked, `scratch/`):
`_426_probe_ground.md` (explorer, the parser→solve path), `_426_remedy_blast_radius.md`
(explorer, the §6b set + design ground), `_426_be_reflected_results.md` (the measurement),
`_426_repro.md` (qa, independent reproduction), `_426_be_reflected_probe.py` (the instrument).

---

## 0. The measurement that chose the scope — `[M]` 2026-09-03, `main` `1e02f6b1`

**Instrument.** `scratch/_426_be_reflected_probe.py`. The ℓ = 1…6 MT=16 moments are read off the
GENDF tape with the production parser, yield-stripped with the SAME per-row diagonal production
applies at ℓ = 0, reversed to canonical order, scaled by number density, and injected as scattering
moments `SigS[ℓ] += 2·Σ₂ₙ,ℓ` (the 2 is `N2NKernel.multiplicity`; the `(2ℓ+1)` factor lives on
`LegendreBasis.reconstruct`, so a stored moment enters `SigS[ℓ]` exactly as the shipped elastic
P1/P2 do — same MF=6 record type, same parser loop). `SigS[0]`, `Sig2`, `SigT` untouched: the ℓ = 0
channel stays on `N2NOperator`, nothing is double counted. Every material's `SigS` is zero-padded
to `L_solve + 1` so the solver's silent min-over-materials clamp cannot lower the order. 1-D slab,
vacuum both sides, GL S8, **421 groups**, `keff_tol=1e-9, flux_tol=1e-8, inner_tol=1e-10`, every
arm `fully_converged`.

**Controls (all read as designed).** C0: the pipeline's ℓ = 0 product equals the shipped
`Mixture.Sig2` **bit-exactly (0.000e+00)** on every material of every fixture. C_pad: the shipped
mixtures padded to L = 6 give the baseline k to **0.000 pcm**. C_sign: the ℓ ≤ 2 moments with their
sign flipped move k the OTHER way by a comparable magnitude (linear). C_refl: the reflector alone
carries the effect (the fuel's own (n,2n) ℓ ≥ 1 is < 2 pcm — U-235's MT=16 is 13× weaker).

⚠ **Convention.** `pcm` = 10⁻⁵ and says nothing about what was divided by what. The ladder columns are **Δk·10⁵ (absolute)**; the three-convention column gives ℓ = 1's reading as Δk·10⁵ / Δk/k₀·10⁵ / Δρ·10⁵ = (1/k₀ − 1/k)·10⁵. They differ by k₀, and on the k₀ = 1.53 fixture the choice inverts the thin-vs-thick comparison (qa finding D1; every row of every convention is in `scratch/_426_be_reflected_results.md`).

| fixture (Be-9 N = 0.1236 \| core \| Be-9) | k shipped (L=2) | ℓ=1 added: Δk·10⁵ / Δk/k₀·10⁵ / Δρ·10⁵ | ℓ=2 | ℓ=4 | ℓ=6 | C_pad | C_sign (ℓ≤2) | C_refl-only |
|---|---|---|---|---|---|---|---|---|
| fast: Be 3 cm \| U-235 metal 4 cm (N 0.04894) \| Be 3 cm; 12/16/12 cells | 1.095322188 | **−413.6 / −377.6 / −346.0** | −412.3 | −413.8 | −413.8 | +0.000 | +409.8 | −412.0 |
| fast: Be 10 cm \| same core \| Be 10 cm; 40/16/40 | 1.526231521 | **−529.3 / −346.8 / −228.0** | −529.0 | −529.8 | −529.9 | +0.000 | +520.8 | −528.9 |
| thermal: Be 10 cm \| U-235 5e-4 + H 0.0669 + O 0.0334, 30 cm \| Be 10 cm; 40/60/40 | 1.745071904 | **-155.6 / -89.2 / -51.1** | (remaining arms running at plan time — `scratch/_426_be_reflected_results.md`) | | | | | |

Flux: relative L2 change of the normalised 421-group flux 1.8e-3 / 3.4e-3 / 2.2e-3; group-integrated
cell flux max relative change 5.9e-3 / 7.9e-3 / 9.3e-3 — concentrated in the reflector.

**Reading.** The ℓ = 1 moment carries the whole effect (ℓ = 2…6 add < 2 pcm); the sign is the physics
(forward-peaked (n,2n) emission in the reflector sends the emitted pair outward, less returns); the
magnitude is 2–4 hundred pcm of reactivity (3–5 hundred in Δk·10⁵) on a Be-reflected fast system and
50–150 even behind a water-moderated core. ⟹ **#426 is a DEFECT, not a documentable approximation.** The remedy is the carve.
*Independent reproduction:* `scratch/_426_repro.md` (qa agent, own instrument, own pipeline) — **REPRODUCED**
2026-09-03: all three k values to every published digit, C0 exactly 0.0, six controls as designed, the two
shared conventions (index orientation: 8195 of 8195 entries upper-triangular; Legendre normalisation:
|Σ_ℓ|/Σ_0 ≤ 1 on every entry, a stray (2ℓ+1) would read 2.9) closed against PHYSICS, not against this probe.
One shared premise stays assumed: the per-row yield scale applies to every ℓ (a yield y(E) multiplies the
whole emission distribution). **CLOSED by the format spec**: ENDF-102 Eq. (6.3) `f = Σ_ℓ [(2ℓ+1)/2] f_ℓ P_ℓ`
with Eq. (6.1) `σ_i = σ·y·f/2π`, and NJOY Eq. (242), give `σ_Xℓ = σ(E)·y(E)·f_ℓ` for EVERY ℓ — one yield per
incident energy multiplies the whole Legendre stack, so the ℓ = 0 per-row scale is the ℓ ≥ 1 scale
(`.claude/agent-memory/literature-researcher/endf6_gendf_njoy_n2n_formats.md`, "The Legendre convention").

⚠ Configuration of every number above: pure-isotope mixtures (σ₀ = 0 → the interpolation clips
to σ₀ = 1 b, identically in both arms), 294 K, the shipped elastic/inelastic/thermal stacks at
their ingest-truncated **P2** in every arm. The ladder is the (n,2n) anisotropy alone.

## 0b. Ground facts about the tree (all `[M]`, from the two explorer memos, re-derivable)

- **The drop is one subscript**: `orpheus/data/micro_xs/gendf.py:378` keeps `sig2_data[(0, 0)]` of a
  dict `_extract_mf6` already fills for ℓ = 0…6. Every MT=16 section has ONE σ₀ column (11 of 11
  tapes); NL = 7 on 10 of 11, NL = 1 on NA023 (the parser pads to 3 zero keys); absent on B_010/H_001.
- **The same site truncates every scattering channel at P2**: `_init_scattering` / `_accumulate_scattering`
  hard-code `range(3)` (`gendf.py:518`, `:545`) while the tape carries 7 elastic orders (× 6 or 10
  σ₀ columns). `compute_macro_xs(n_legendre=3)` and `interp_sig_s`'s "0, 1, or 2" are downstream copies
  of that 3. **The elastic P3…P6 effect is UNMEASURED.**
- **The moment machinery above the material field is channel-agnostic**: `HarmonicFrame.conjugate` /
  `reconstruct_after` / faces / the `(iso/W) + aniso` combine / `OperatorProduct.apply_transpose`.
  `ScatteringMaterialField._moment_blocks` is duck-typed over `kernel.moments[l]` / `.order` but has
  **no scale hook** (only the P0 verbs take `scale=`, where `N2NKernel.multiplicity` enters);
  `LegendreMomentScattering` is annotated to the scattering-named field; `ScatteringKernel` is a
  Legendre stack with implicit multiplicity 1, `N2NKernel` a rank-2 matrix with `multiplicity = 2`
  and no stack. `N2NOperator` mints its frame at `for_space(interior, 0)` (`n2n.py:200`) and is
  `ScatteringOperator`'s arms with `aniso=None`; `IsotropicN2N` is member-for-member `IsotropicScattering`.
- **Exactly ONE production read of `Mixture.Sig2` is the moment seam** (`kernels.py:239`); the other 17
  are ℓ = 0 by physics (row sums for removal / `SigT` / `absorption_xs` / the k balance
  `production/(absorption + leakage − emission_n2n)`; CP/MoC/MC P0 sources; the dense cache; the K_iso
  leaf `S.isotropic_energy + N2N.energy` for the ψ½ seed).
- **The HDF5 cache** stores one `sig2` triplet and `sigS/L{j}_S{k}` (loader derives the ℓ count from the
  keys). The `*.h5` files are UNTRACKED (`.gitignore`), regenerated by
  `orpheus/data/micro_xs/convert_gxs_to_hdf5.py`; there is **no round-trip test** (0 test sites of
  `hdf5_io`); recipe tests SKIP when the h5 is absent — a silent-green hazard for a schema change.
- **§6b surface of the `Sig2` retype**: production 18 attribute reads + 17 constructor kwargs; tests 42
  constructor sites (6 of them FACTORIES — `material_xs_from_raw`, `_config._mix`, `xs_library.make_mixture`,
  `homogeneous._make_mixture`, two MMS makers) + 38 attribute reads, every one of which breaks on a list.
  The gates that ENCODE the ℓ = 0 truncation and must flip: `test_n2n_moment_emission_touches_only_l0`,
  `TestLiftIsTheConjugation` at the L=0 frame, the tier-2 equivalence rows at `for_space(interior, 0)`.
  Every (n,2n) fixture in the tree is P0-only (#269's residue) — blind to ℓ ≥ 1 by construction.
- **Layer contract**: `data` may not import `transport` (why the yield strip is an MF=3 renormalisation
  in `data/` and the multiplicity's one home is `transport/kernels.py`); `transport` may not import `sn`.
- **Dead code found on the way**: `MaterialXSField.with_overridden_sig_s_and_n2n` (0 callers; its
  docstring names consumers that construct the field directly) and `derivations/diagnostics/
  diag_276_scattering_p0_fastpath_perf.py` (calls signatures that no longer exist).

---

## 1. Goal (the domain's terms)

**A collision channel's secondaries carry the angular distribution the evaluation gives them.** The
multigroup collision source of a channel `c` with yield `y_c` and transfer moments `Σ_{c,ℓ}(g'→g)` is
`q_c = y_c Σ_ℓ (2ℓ+1)/(4π) Σ_{g'} Σ_{c,ℓ}(g'→g) φ_ℓ(g') P_ℓ(Ω)` — scattering is the case `y = 1`,
(n,2n) the case `y = 2`; nothing else distinguishes them in the transport equation. The data layer is
therefore **lossless in ℓ** (it keeps what the tape carries), and truncation is the SOLVE's decision,
spelled once (`scattering_order`). The algebra keeps its two TERMS `(L+C) − S − N2N − B` (CS4c §14.1:
the channel's bundling is context-dependent), realised as two INSTANCES of one family.

**Done when** (checkable):
1. `Isotope.sig2` / `Mixture.Sig2` carry the tape's Legendre stack (a tape-vs-Isotope pin on Be-9's
   ℓ = 1 is non-zero; an HDF5 round-trip test exists and passes).
2. `solve_sn(..., scattering_order=2)` on the fast Be-3-cm fixture of §0 reads `k = 1.0912 ± 0.0001`
   (today `1.095322188`) with the SHIPPED library and NO probe — a tracked gate reading tracked tapes,
   pinning the sign, the ℓ-ladder's convergence at ℓ = 1, and bit-identity at `scattering_order = 0`.
3. `grep -rn "N2NKernel\|N2NMaterialField\|N2NMomentOperator\|ScatteringKernel\|ScatteringMaterialField\|LegendreMomentScattering" orpheus/`
   returns only history prose: the (n,2n) term is the `TransferOperator` core with `multiplicity = 2`,
   wearing its role name `N2NOperator`; the role subclasses carry no arithmetic (an AST gate says so).
4. `tests/sn/architecture/test_monomorphic_leaves.py`'s ledger and every gate in §0b's "must flip" list
   are re-keyed with an ℓ ≥ 1 fixture that has teeth (a mutation zeroing the ℓ = 1 (n,2n) block reddens).
5. `sphinx -E -W` clean from the repo root; `dead_references` 0; a new error-catalogue entry records the
   truncation as the defect it was; the 13-tree gate 13 of 13.

---

## 2. Forks for the ruling — proposed 2026-09-03; ALL THREE RULED the same day (§2b carries F2/F3)

**F1 — how lossless is the ingest?** ✅ **RULED (user, 2026-09-03): (b) — every channel keeps the tape.**
- (a) (n,2n) keeps its 7 orders; elastic/inelastic/thermal stay at the hard-coded P2.
- (b) **RECOMMENDED** — EVERY channel keeps what the tape carries (7 orders); the `3` disappears from
  `gendf.py`, `compute_macro_xs`, `interp_sig_s`; a mixture's order is a fact about its data, not a
  parser constant; `scattering_order` is the only truncation. Cost `[R, unmeasured]`: h5 size ×~2.3 for
  the elastic stack, `compute_macro_xs` σ₀-interpolation ×7/3; both measured at step 1 before landing.
  By-product: the elastic P3…P6 effect becomes measurable with the §0 instrument for free.
  Why (b): §0b's second bullet is the same defect one channel over, and a data layer that silently
  truncates is the lossy-return-type root cause ([[feedback-lossy-return-type-is-the-root-cause]]);
  #60's "order consistency validation" then has one honest home — the solve's clamp/refusal.

**F2 — the operator family.** ✅ RULED — (a), after the §2b analysis the user asked for; §2b carries the ruling and its reasons.
- (a) **RECOMMENDED** — ONE scattering-type family, instantiated twice: the kernel carries the Legendre
  stack AND its multiplicity as a datum (`ScatteringKernel` absorbs `N2NKernel`); the material field's
  moment verbs take the scale the P0 verbs already take; `IsotropicScattering` absorbs `IsotropicN2N`
  (CS4c step 5's parked redirect (a), corroborated by the {S, N2N} | {F} data split);
  `ScatteringOperator.from_solver_data` builds the N2N term at the SAME `scattering_order` and the same
  interned frame (`N2NOperator`, `N2NMomentOperator` dissolve; `full_n2n_kernel` becomes the family's
  `full_scatter_kernel`). The σ_r fold family stays scattering-only (it is about within-group
  self-scattering; whether (n,2n)'s within-group block folds is a later question — record, don't build).
- (b) N2N mirrors S: give `N2NKernel` a stack, `N2NOperator` an aniso arm. ⛔ A twin path of every
  member the explorer measured at 0.85–1.0 similarity; Cardinal Rule 2 rejects it.

**F3 — the family's name.** ✅ RULED — neither option as first posed: the KERNEL tier is `Transfer*`, the OPERATOR tier keeps its role names; §2b.
- (i) **RECOMMENDED** — keep `Scattering*`. (n,2n) IS a scattering-type collision (the standing ruling:
  "N2N belongs with SCATTERING; `Isotropic…` records a truncation, not a property"); the yield is a
  datum of the kernel; the algebra's terms keep their names as roles (`S`, `N2N` in
  `explicit_gains=(S, N2N, B_a)`). `[M]` a rename to `Transfer*` would touch `ScatteringOperator`
  alone in 73 files / 282 occurrences (orpheus 18, tests 33, docs 22 files) plus its four siblings —
  churn without a physics gain.
- (ii) rename the family to `Transfer*` (the NJOY/GROUPR term "transfer matrix"; scattering becomes
  the y = 1 instance). Precise, consistent, expensive.

**F4 — the ℓ = 0 consumers.** No fork: removal / `SigT` / `absorption_xs` / the k balance / the K_iso
leaf / CP / MoC / MC read the P0 moment's row sum, which IS the reaction rate; they read `Sig2[0]` (or
the named verb `n2n_out_xs`). CP's and MoC's hand-rolled `multiplicity · Sig2.T @ φ` sites (4) stay as
#279's residue, re-spelled minimally. MC keeps sampling the exit direction isotropically — an honest
per-solver statement the #428 doc split records (the ℓ ≥ 1 PDF is a separate MC feature).

**F5 — HDF5.** No fork: schema `sig2/L{j}` (and `sigS/L{j}_S{k}` up to the tape's order under F1(b)),
a format-version key, a loader that REFUSES the old layout loudly (never a skip), regeneration by the
existing script, and the missing round-trip test lands with the schema.

## 2b. What S and N2N share and where they differ — the analysis behind the F2/F3 rulings

**The tape's lens (the format is the physics).** GENDF MF=6 stores, for EVERY channel, the same object:
`σ_{c,ℓ}(g'→g) = σ_c(g')·y_c(g')·f_{c,ℓ}(g'→g)` — ENDF-102 Eq. (6.1)/(6.3) with NJOY Eq. (242) (`[M]` the
literature memo; `[M]` #427 measured the folded yield at 2.000000 on Be-9 MT=16; the qa reproduction at
1.9999999 / 2.0000000 on Be-9 / U-235). This is the **emission stack**: neutrons emitted into `(g, ℓ)` per
unit incident flux in `g'`, yield included. Elastic has `y = 1`; (n,2n) has `y = 2`; fission has `y = ν`
(GROUPR's production matrix `νσ_f χ` is the same convention). The user's recollection is exact: the
multiplicity is INSIDE the tape's transfer function. `[R, page-check owed]` Bell & Glasstone §1.1b define
the transfer function the same way — `σ f = Σ_x σ_x f_x` over every reaction type `x` with
`∫∫ f_x dΩ dE = c_x`, the mean number of neutrons emerging per collision of type `x` — so "a transfer
function with its multiplicity" is the textbook object and scattering is its `c = 1` member.

**What ORPHEUS did at #427** (`7433507d`): the ingest DIVIDES the yield out (`_strip_transfer_yield`,
per-row `σ_MF3/rowsum`), so `Isotope.sig2` / `Mixture.Sig2` are the **reaction** stack `σ·f`, and
`N2NKernel.multiplicity = 2` re-applies `y` at emission. Why: the removal side (`SigT`, `absorption_xs`)
needs `σ_c` (one neutron absorbed per event) and the gain side needs `y·σ_c·f` — two quantities from one
datum; the tree stores the reaction and lets the kernel own `y`. For scattering the two coincide. **Either
convention needs `y` beside the stack; the choice is which of the two the datum spells.** This plan keeps
the shipped one (reaction stack + `multiplicity`): `rowsum(Sig2[0])` stays `Σ₂ₙ`, the removal/balance
reads stay literal, #427 stands.

**Shared, `[M]` at HEAD** (the two explorer memos, the qa census):
1. The same tape record, parser loop and Legendre convention (`(2ℓ+1)` on reconstruction, none stored).
2. A full, non-separable `g'→g` matrix (Be-9 (n,2n): rank 50, best rank-1 error 58 %; fission is a
   rank-1 dyad by construction) — the {S, N2N} | {F} split.
3. The collision-gain operator at any order: `R·Λ·M / W` with the P0 fast path `(iso/W) + aniso`; the
   transposes are the same product chain; the K_iso leaf for the ψ½ seed is the P0 energy binding of each
   (`S.isotropic_energy + N2N.energy`, `coupled_system.py:592`).
4. Removal counted ONCE in `Σ_t` (`mixture.py:658`), for both.
5. Neither is scaled by `1/k`: both sit in the loss operator `A = (L+C) − S − N2N − B` (ERR-065 / R7,
   `coupled_system.py:551`); the eigenvalue scales fission only.
6. `[M]` AST similarity of the shipped classes 0.85–1.0 member for member (`IsotropicN2N` ≡
   `IsotropicScattering`; `N2NOperator` = `ScatteringOperator`'s arms with `aniso=None`).

**Different, and what each difference IS:**
- (a) **`y`**: 1 vs 2 — a scalar datum of the channel (exactly 2 for (n,2n) at every energy; the tape
  folds it; #427's integer-yield admission pins it).
- (b) **Production accounting** — the user's point, and a CONSEQUENCE of (a): a channel with `y > 1` is a
  net neutron source, `(y − 1)·Σ_{c,0}ᵀφ`. ORPHEUS spells it in exactly two verbs, both on the P0 energy
  binding (`N2NMaterialField.add_to_group_rate`): the k denominator's net removal `absorption + leakage −
  E₂ₙ` (`compute_keff`, `solver.py:1713`; CP `net_removal = total − scatter − 2·n2n`; MoC
  `(σ_a − 2σ₂ₙ,out)φ`) and the ERR-052 scale anchor `compute_production_rate = fission + E₂ₙ`
  (`:1637-1656`: "TOTAL physical production … the renormalisation scale anchor ONLY"). For scattering
  these are identically `(1 − 1)·… = 0`, which is why the code never spells them for `S`. ⟹ written ONCE
  as functions of `y` they are the SAME verbs for every channel and vanish for scattering — no branch, no
  `isinstance` (Pattern 4 by arithmetic).
- (c) **σ₀ columns**: elastic 6–10 (resonance self-shielding), (n,2n) 1 (threshold) — an ISOTOPE-level
  difference; at the `Mixture` level both are one Legendre stack.
- (d) **Values, not structure**: (n,2n) is strictly downscatter (upper-triangular, 8195 of 8195), threshold
  1.75 MeV; elastic + thermal upscatters. Same type, different entries.
- (e) **The σ_r fold**: the within-group block of `S` is folded into removal for the sweep
  (`foldable_part` / `residual_part`); (n,2n)'s `y·Σ₂ₙ,gg` is left explicit. `Σ_r = Σ_t − Σ_{s,gg} −
  y Σ_{2n,gg}` is equally legal — a choice the tree made for one channel, not a structural difference.
  Record; do not build.
- (f) **The one TYPE door**: `solver.py:1214` guards the seedless SI record's gains with
  `isinstance(S, ScatteringOperator) and isinstance(n2n, N2NOperator)` — it enforces the tuple ORDER
  `(S, N2N, B_a)` ("loud on drift"). Under F3's ruling the role classes survive, so the door keeps
  working; a named gains record would still be the better spelling (later, not this plan).

**Reading against the type-vs-property rule** (`coding-standards`): the two realizations are ISOMORPHIC
(one Legendre-stack shape) and the only morphism that distinguishes them is *scale by `y`* — the same
arithmetic on both. ⟹ `y` is a **property (datum) of the kernel**, not a type. The first-classness the
§14.1 ruling protects — N2N is its own TERM, bundled with S for anisotropy and with F for production
accounting — lives in the ALGEBRA (two objects in two role slots, `explicit_gains = (S, N2N, B_a)`,
`self.scattering_op` / `self.n2n_op`) and in the two `y`-verbs. Fission stays its own type by the same
rule: a separable dyad (non-isomorphic realization) under a different morphism (`1/k`).

### ✅ F2 RULED (user, 2026-09-03): one transfer family, `y` a kernel datum, two instances

Kernel `(moments, multiplicity)` — the shipped `ScatteringKernel` stack with `N2NKernel`'s multiplicity;
the field's moment verbs take the scale the P0 verbs already take; the two accounting verbs become
`emission = y·Σ₀ᵀφ` and `net_production = (y − 1)·Σ₀ᵀφ` on the field, spelled once; one P0 energy binding
core; one angular binding core built twice at the same `scattering_order` and interned frame.
`compute_keff` / `compute_production_rate` keep the ERR-065 / R7 convention (not re-posed here) and read
the N2N instance's verbs. Rejected: (b) shared substrate with two role classes carrying their own
arithmetic — no; (c) N2N mirrors S — the twin path.

### ✅ F3 RULED (user, 2026-09-03): "the kernel and the operator name different things"

*"The Kernel can be TransferKernel, but the Operator should be the ScatteringOperator or N2NOperator."*
The principle, articulated: **the kernel tier names the MATHEMATICAL OBJECT; the operator tier names the
TERM of the algebra it realises** — as `StreamingOperator`, `CollisionOperator`, `FissionOperator` already
do. Realisation:
- **Kernel / field / moment tier — `Transfer*`**: `TransferKernel(moments, multiplicity)` (absorbs
  `ScatteringKernel` + `N2NKernel`), `TransferMaterialField` (absorbs `ScatteringMaterialField` +
  `N2NMaterialField`; the moment verbs take `scale`; `emission` / `net_production` spelled once),
  `LegendreMomentTransfer` (absorbs `LegendreMomentScattering` + `N2NMomentOperator` — the `Λ` factor is
  kernel-tier math, not a term of the algebra).
- **Operator tier — role names on shared cores**: the angular binding core `TransferOperator`
  (`ScatteringOperator`'s body, channel-agnostic) with the two terms as **thin role subclasses**
  `ScatteringOperator(TransferOperator)` and `N2NOperator(TransferOperator)` whose ONLY content is the
  extraction classmethod (which `Mixture` channel `from_solver_data` / `from_material_xs` reads) and the
  role name; the P0 energy binding core `IsotropicTransfer` with role subclasses `IsotropicScattering` /
  `IsotropicN2N` the same way (the K_iso leaves, diffusion, CP and homogeneous consumers keep their term
  names; "Isotropic" is honest here — it IS the ℓ = 0 projection those consumers want).
- ⛔ **Gate the principle**: an AST test asserts the role subclasses define NOTHING but classmethods and
  `ClassVar`s — no `apply`, no `apply_transpose`, no field — so a twin path cannot regrow one override at
  a time. `isinstance(S, TransferOperator)` is `True` for both; the door at `solver.py:1214` keeps working.
- What this changes in §1's done-when 3: the grep is for the RETIRED spellings only — `N2NKernel`,
  `N2NMaterialField`, `N2NMomentOperator`, `ScatteringKernel`, `ScatteringMaterialField`,
  `LegendreMomentScattering` return history prose only; `N2NOperator` / `IsotropicN2N` survive as role
  classes with no arithmetic.

---

## 3. Step ladder — hypothesis, re-derived at each phase opener

- **Step 1 — the data layer is lossless in ℓ** (bit-identical exit: no operator reads ℓ ≥ 1 of (n,2n)
  yet, and `scattering_order ≤ 2` truncates the wider elastic stack to today's). `Isotope.sig2`,
  `Mixture.Sig2` → `list[csr_matrix]` over ℓ (single σ₀); the ingest keeps the tape (F1); the yield strip's
  diagonal applied per ℓ (derived at ℓ = 0 — the integer-yield admission is an ℓ = 0 statement);
  `_to_canonical_group_order`, `condense`, `from_dense_channels`, `compute_macro_xs`, the h5 schema +
  round-trip test (F5); the 17 + 18 production sites and the test surface. §6c witness at landing: the
  tape-vs-Isotope ℓ = 1 pin (Be-9), the round-trip test. Retire the two dead items of §0b.
- **Step 2 — the (n,2n) term is anisotropic** (F2/F3 as ruled in §2b): `TransferKernel(moments,
  multiplicity)`; `TransferMaterialField` with scaled moment verbs and the two `y`-verbs; `IsotropicTransfer`
  and `TransferOperator` cores with the four thin role subclasses; the Be-reflected gate
  flips to its measured value (§6c: the witness IS the shipped library); the "must flip" gates re-keyed
  with an ℓ ≥ 1 fixture; the K_iso leaf and the adjoint transposes unchanged in code, re-verified.
- **Step 3 — the corpus**: `adjoint.rst` (`sn-n2n-p0-truncation` → dated history; `sn-n2n-isotropic-lift`
  and `sn-n2n-adjoint-source` gain their per-ℓ form — a labelled equation is an API), `slab_multigroup.rst`
  `_n2n-reactions`, `cross_section_data.rst` (the P2 statement, the h5 schema, the nnz table; #428's split
  lands here too), the error catalogue entry, the in-code docstrings the explorer listed; `Closes #426`.

Proactive dispatches: **test-architect before step 1** (the carve crosses data ↔ transport ↔ sn);
elegance-enforcer after each step; archivist for step 3. Batteries in-process, crash-safe. The 13-tree
gate before each merge, pre-registered per tree against the 11007 baseline.

## 4. Standing constraints

Those of `.claude/plans/cs4c_binding_design.md` §18.7, unchanged. Baseline `[M]` `main` `1e02f6b1`:
11007 collected, 13 of 13 rc=0. Branch `fix/n2n-anisotropy` carries one prior commit (`8707c53a`, the
vv-principles #29 sharpenings the step-5 census owed).
