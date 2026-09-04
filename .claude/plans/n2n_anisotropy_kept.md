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
carries the effect (the fuel's own (n,2n) ℓ ≥ 1 is 0.2 pcm·10⁵ — U-235's (n,2n) emission is 13× LESS
forward-peaked than Be-9's, μ̄ +0.023 vs +0.30, and open in 22 incident groups against Be-9's 50; `[M]` its peak
reaction XS is not smaller — 0.813 b vs 0.559 b — the archivist's correction of an earlier gloss).

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
  parser constant; `scattering_order` is the only truncation. Cost `[M]` 2026-09-03 (`scratch/_426_f1_cost_probe.py`, tapes at 294 K): every shipped scattering
  section stores NL = 7 (elastic, inelastic, thermal; U-235 n_σ₀ = 10, Be/O 6, H 1); raw COO bytes over all
  scattering sections grow ×1.98 (U-235) to ×2.33 (Be, H, O) — the shipped 3-order h5 files are 6.9–50 MB
  each, untracked and regenerated by script; `interp_sig_s` over 3 orders costs 0.14 s on U-235
  (0.04 s Be/O), so 7 orders ≈ 0.3 s per isotope per mixture build. Negligible on both axes.
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

## 2c. The test-architect's plan (2026-09-03, `scratch/_426_verification_plan.md` + five drafts) — rulings on its five open points

`[M]` 41 draft rows / 22 RED today / 19 pass / 54 s; predicted per-tree deltas data +41, transport +23,
sn +21, diffusion +2, mc +1, root+harness +1 (the layer gate is parametrized over 346 modules: ±1 per
production module) ⟹ **11007 → 11096** if every design row lands as drafted. Its refutations that
change the design:

- **R1 ⛔⛔ the solver's silent order clamp** `L = min(scattering_order, min(len(m.SigS) − 1))`
  (`sn/solver.py:1359`) must NOT read the (n,2n) stack: `[M]` H-1 and B-10 carry no MT=16, so a
  two-list `min` would drop every mixture containing them to P0 — deleting the elastic P1/P2, worth
  +5787 pcm-relative on the §0 fixture. **O-1 RULED (main agent): (a) — the clamp reads `SigS` alone;
  a channel that stores fewer orders than requested is ZERO there, by the evaluation's own statement
  (absent section, or NL = 1 = declared isotropic), so the (n,2n) stack zero-pads to the requested
  order.** Padding is the exact datum, not an approximation. The elastic clamp itself stays as today
  through step 1 (bit-identical); #60 owns its future.
- **O-2 RULED: an isotope's (n,2n) stack has the TAPE's order** — `[Σ₀]` alone when the section is
  absent (H-1, B-10) or NL = 1 (NA023); the `_extract_mf6` pad-to-3 goes (the ingest invents nothing);
  the Mixture sums stacks of unequal length by zero-padding to the longest, and the SAME helper sums
  `SigS` (Pattern 2).
- **R2 ⚠ step 1's "bit-identical by design" has the denominator `scattering_order ≤ 2`**: `[M]` no
  shipped GENDF-library solve runs above P2 (130 × `=0`, 54 × `=1`, three synthetic `=3`), so F1's
  un-clamping lands with zero witnesses unless one is written — **G1.10 (a library solve at
  `scattering_order = 3` differs from 2) lands with step 1**, and the elastic P3…P6 EFFECT is measured
  at step 1's exit on the §0 fixture. **O-4 RULED: the elastic higher orders LAND with F1 (they are the
  data) and are MEASURED at step 1's exit; no separate remedy — the solve already honours
  `scattering_order`.** `[M]` BE009 elastic `max|·|` per ℓ: 48.4, 3.61, 1.70, 0.664, 0.454, 0.370,
  0.286 — ℓ = 3 is 18 % of ℓ = 1, not noise. ✅ **MEASURED at step 1's exit** (`scratch/_426_elastic_ladder.log`,
  the §0 fast/thin fixture, (n,2n) still P0): k(P2) = 1.0953221881419453 (= the pre-carve record, bit-identical),
  P3 1.0930284426239356, P4 1.0936923910519523, P5 1.093572139010408, P6 1.093592425144843 — relative to P2,
  Δk·10⁵ / Δk/k₀·10⁵ / Δρ·10⁵: P3 **−229.4 / −209.4 / −191.6**, P4 −163.0 / −148.8 / −136.0, P5 −175.0 /
  −159.8 / −146.1, P6 −173.0 / −157.9 / −144.4. The ingest's P2 cut on elastic was a truncation error of the
  SAME order as the (n,2n) anisotropy on this fixture; until step 1 a `scattering_order = 3` request was
  silently served P2. Gate: `tests/sn/solve/test_scattering_order_is_the_only_truncation.py`.
- **R3/R4**: a monotone-decay leg is a FALSE red (Be-9 MT=221 has ℓ = 6 > ℓ = 5); the `|Σ_ℓ| ≤ Σ_0`
  bound is one-sided (catches a stray (2ℓ+1) or ×2, blind to deflation; Be-9 MT=2 reaches 0.9997) —
  the two-sided catcher is **G1.2: the yield strip is ONE diagonal** (ratio invariance across ℓ).
- **R9 ⚠ the rank-2 harmonic head has no witness in a slab-only plan**: `_block_contraction` dispatches
  on the head's rank; every incumbent (n,2n) fixture is a 1-D rule. Step 2's moment gates parametrize
  over BOTH head ranks (a `product`/`lebedev` row).
- **O-3 RULED: the flagship's two RECORD rows pin k at `rtol = 10 × keff_tol = 1e-8`**, not
  `array_equal` — a record of a converged eigenvalue is a claim at the solver's tolerance; step 1's
  ingest ledger (G1.9) IS `array_equal` (stored values must not move). Step 2 states per path whether
  it is bit-identical or principled-equivalent.
- **O-5 RULED: frozen artefacts live at `tests/<tree>/data/`** (#444; the numerics tree's convention) —
  `tests/data/data/pre_426_ingest_ledger.npz` for the step-1 ledger; the flagship's two k records are
  literals in the test with their provenance.
- **H-c**: the h5 store is untracked and regeneration costs `[M]` 5.4 min today (~2.3× after F1) — the
  new loader REFUSES a pre-F1 file loudly (never a skip) so an un-regenerated checkout fails at the
  first recipe with the regeneration command in the message.
- **H-l ⛔**: the `+1e-30` epsilon the ingest adds is what makes `interp_sig_s`'s shared-sparsity
  assumption true — `[M]` BE009 elastic has genuine exact zeros at ℓ = 3/5/6 — so the widened loop keeps
  the epsilon on every order, or the σ₀ interpolation reads wrong positions silently.
- **H-i**: 18 attribute STORES `x.Sig2 = …` in tests are a §6b class the "38 reads" census did not
  contain; the migration rewriter handles Store targets too.
- **H-m**: `transport_xs` (→ D → DSA) reads `SigS[1]` only and stays (n,2n)-blind; the carve keeps it
  so and says so (a `.. warning::` at step 3).
- **G2.1 is NOT `slow`** (26.6 s of a ≥ 90-min gate; #36's lesson); **G2.9** gives MC a 0.9 s fast
  companion witness (17 σ under the mutation); **G2.7** adds `homo_2eg_n2n` to the SN k_inf ladder.
- **An ERR entry couples on two harness arms** (a `catches` marker + the regenerated
  `error_index.md`): entry + catcher + Sphinx build are ONE commit.

---

## 3. Step ladder — hypothesis, re-derived at each phase opener

- **Step 1 — the data layer is lossless in ℓ** (bit-identical exit: no operator reads ℓ ≥ 1 of (n,2n)
  yet, and `scattering_order ≤ 2` truncates the wider elastic stack to today's). ✅ LANDED `f96de34c`
  (2026-09-03; §5 carries the gate). `[M]` exit: all 13 pure-isotope mixture digests (SigT, SigP, χ,
  SigS[0..2] row sums, Sig2[0] row sums) and the fast/thin k at L = 0/1/2 reproduce to the BIT — ⚠ with ONE
  exception the "bit-identical" claim did not foresee: NA023's ℓ = 1, 2 scattering rows lost the `1e-30`
  phantoms the retired `_extract_mf6` pad had salted in at its NL = 1 MT=91 section's 5250 positions (nnz
  16809 → 13196; row sums moved ≤ 5250e-30; SigT unchanged; the honest datum is zero) — re-baselined in the
  ledger with the reason. Data tree 313 passed (237 + 76). ⚠ One §6b spelling the rewriter mis-fired on: the test
  helper `material_xs_from_raw(sig2=...)` takes a per-material DICT under the datum's lowercase name, and a
  name-keyed kwarg rule wrapped six of its call sites in a list (9 reds, one loop) — a homonym parameter is a
  member of the rename's NEGATIVE set, and only the consumers' red loop finds it. Battery (`scratch/_426_step1_battery.py`, in-process monkeypatch, one process per arm, over
  `tests/data`; baseline 0 reds): A cut every channel at three orders → 11 reds (the tape-threading rows;
  the store-reading ledger/round-trip rows are BLIND to an ingest mutation by design — they pin the store);
  B the yield scale powered per order → 5 reds (the one-diagonal rows for ℓ = 2…6; ℓ = 1 is a null arm since
  `scale¹ = scale`; the entrywise-bound rows stayed GREEN — the measured proof of R4's one-sidedness); C the
  higher orders left un-reversed → 12; D positive control (no yield strip) → 9; E the macroscopic sum drops
  short stacks → 1 (the new mixing witness). Elastic-order exit measurement: §2c O-4. Elegance review (enforcer,
  2026-09-03): 3 violations + 5 should-fix + 8 nits, ALL taken — the epsilon's rationale was FALSE on NA023
  (orders need not share P0's pattern; the invariant `interp_sig_s` needs is across σ₀ columns of ONE order)
  and is now the named `_SPARSITY_EPSILON` with the measured statement; the retype had minted a byte-identical
  gather pair in `MaterialXSField` (now ONE `_gather_stack(channel, order)`); `interp_sig_s`'s pad guard was
  an unreachable second spelling of the pad rule claiming to be the mechanism (deleted); the Legendre-stack LAW
  (non-empty, square (ng, ng), both channels) now lives in `Mixture.__post_init__` / `Isotope.__post_init__`
  as a real raise. ⏸ DEFERRED to step 2's opener (same branch): `_legendre_order` returns a COUNT the producer
  `_extract_mf6` already knew (`n_lgn`) — widen the section return to carry it (a named product, Pattern 3),
  and make the sig2/sigS paths agree on a skipped order (KeyError vs silent zero today). ⭐ Collapse trigger
  recorded: when step 2 mints `TransferKernel`, `_gather_stack`, `_macroscopic_stack`'s pad and the stack law
  all become the kernel's constructor and `moment(ℓ)`. Two residue issues filed: `SigT` derived three times;
  `transport_xs` is (n,2n)-blind now that `Sig2[1]` exists (DSA/diffusion consumer). **13-tree gate PRE-REGISTRATION** (baseline `main` `1e02f6b1`: 11007): data +78 (h5 store 12,
  ingest ledger 43, yield-convention file +21 — 1 guard row, 6 tape-physics, 13 threading, the two module
  fixtures cost none — canonical order +1, condense +1), sn +3 (the clamp/un-clamping gate), root+harness +0
  (no new production module: the layer gate's 346-module parametrization is unchanged), every other tree +0
  ⟹ **predicted 11088 collected, 13 of 13 rc=0**. `[M]` MEASURED 11087, 13 of 13 — reconciled per file: the
  yield-convention file gained +20 (I had added +21) and the ledger's 43 already held the mixing witness I
  later counted again; every row that exists is accounted for (data +77, sn +3, all other trees unchanged). `Isotope.sig2`,
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
  ✅ CARVE COMMITTED `1a3b78ec` (2026-09-04; opener `7b44ee68` = the deferred `_legendre_order` producer fix:
  `_extract_mf6` returns `_MF6Section(ifrom, ito, moments, n_legendre, n_sigma_zero)` and refuses a
  non-rectangular section — `[M]` 13 of 13 tapes rectangular, elastic NZ = header n_sig0 on 13/13; the
  sig2/sigS paths now both index the section directly inside `range(section.n_legendre)`). Gate status at
  the commit: every affected tree green (transport 738, sn/operators 1277, sn/solve 208, analytical 62+12,
  data 314, diffusion 117, homogeneous 50, numerics 3255, cp 141, moc 121, mc 40, architecture 152, mms 23);
  collected **11087 → 11129**, reconciled per file against a HEAD worktree collection (+21 transfer_kernel,
  +7 roles, +5 flagship, +3 n2n_operator, +3 material_field, +2 diffusion witness, +1 mc companion, +1 layer
  gate for `transfer.py`, −1 the retired diagnostic's registry row); the 13-tree gate, the battery, the
  elegance review and the corpus pass follow on the branch.
  `[M]` the flagship's ladder (`scratch/_426_step2_ladder_record.txt`; the §0 fast/thin fixture, elastic P2
  in every arm): the P0-only control reads **1.0953221881419453 — bit-identical** to the pre-carve record
  (measured, not argued: the y = 1 fast path is skipped and an all-zero Λ_{ℓ≥1} is skipped by `is_isotropic`);
  ℓ ≤ 1 **1.0911866898558749** (−413.55 / −377.56 / −346.01, the §0 table to every digit); ℓ ≤ 2
  **1.0911996566537725** (−412.25 / −376.38 / −344.92); ℓ ≤ 6 = ℓ ≤ 2 **to the bit** (a P2 solve never reads
  ℓ ≥ 3 — the A9 "declared null" arm confirmed by construction); k(L = 0) **1.1587120371368607** bit-identical.
  ⚠ REALIZATION REFINEMENTS that supersede §4.1's table (edited here, not there — §4.1 is the pre-carve
  design): (a) the angular core lives in a NEW module `operators/transfer.py` with `LegendreMomentTransfer`
  (the test-architect's draft placed it there; a core named `TransferOperator` in `scattering.py` misfiles
  it), the P0 core `IsotropicTransfer` stays in `isotropic_scattering.py` beside `IsotropicFission`, the roles
  keep their modules so the 107 rst xrefs to `ScatteringOperator` and the `N2NOperator` ones stay valid;
  (b) the tier-2 mints are `TransferOperator.from_field(transfer, *, scattering_order, space)` and
  `LegendreMomentTransfer.from_field(transfer, basis, *, skip_l0)` on the cores, `from_solver_data` /
  `from_material_xs` on the roles — the moment operator's `from_material_xs` (channel-ambiguous on one class)
  retired; (c) the role→P0-binding link is a `ClassVar` `isotropic_binding` on the core (default the core
  itself), overridden per role, so `S.isotropic_energy` IS an `IsotropicScattering` and `N2N.isotropic_energy`
  an `IsotropicN2N` (the terms keep their names down to the K_iso leaves); (d) `at_order` (not `truncated`)
  returns `self` at the identity, truncates below, pads exact zeros above (§4.3); (e) **`is_isotropic`** on
  the kernel/field/operator — Λ_{ℓ≥1} exactly zero — generalises the `scattering_order == 0` early return
  from SHAPE to VALUES: both bindings now run at the solve's order and the (n,2n) stack of an isotope with no
  MT=16 is exact zeros padded to L, so without it every SN solve at L ≥ 1 would run a second RΛM product on
  zeros (bit-identical either way; the skip is the "performance lazy" half of the algebra-eager lens);
  (f) ⛔ REFUTED: F2's sentence "the two accounting verbs become `emission = y·Σ₀ᵀφ` and `net_production =
  (y − 1)·Σ₀ᵀφ` on the field" — `net_production` has NO consumer (`compute_keff` keeps the ERR-065/R7
  convention `absorption + leakage − emission`, F2's own "not re-posed here"), so it is not minted
  (coding-standards defer-until-consumer); the verbs that exist are the existing set with the yield as
  `scale`, and `add_to_group_rate` on a y = 1 field is the group in-scatter rate; (g) ⛔ REFUTED: step 1's
  "collapse trigger" (`_gather_stack`, `_macroscopic_stack`'s pad and the stack law all become the kernel's
  constructor and `moment(ℓ)`) — the three pads are three TIERS' data (the isotope sum over sparse stacks, the
  σ₀-frame projection over the facade, the dense per-material kernel); one law, three data, no shared home
  possible without one tier reaching into another; the facade docstring says so now.
  ✅ BATTERY (`scratch/_426_step2_battery.py`, in-process monkeypatch, one subprocess per arm, bite check per
  arm, `^ERROR` counted apart from `^FAILED`; scope = transport, diffusion, homogeneous, sn/operators,
  sn/verification/analytical, the clamp gate, the yield-convention file; A1/A2 add cp/moc/mc; reds as
  NEW-file / PRE-EXISTING): baseline **0**; A1 ν₂ₙ 2→1 **50** (7/43); A2 2→3 **48** (6/42); A3 the (n,2n) ℓ=1
  block zeroed **10** (3/7); A4a the yield dropped in `_moment_blocks` **11** (2/9); A4b in `add_p0_source`
  **39** (4/35); A4c in `add_p0_source_transpose` **3** (0/3); A4d in `add_to_group_rate` **13** (4/9);
  A5 an `apply_transpose` override on the N2N role (via a patched `inspect.getsource`) **1** (the AST gate's
  N2N row); A6 the N2N frame minted at L=0 — the retired defect re-introduced — **9** (3/6, the flagship's
  value/Δ/ladder rows among them); A7 the MF=6 section's NL read as 1 **13** (0/13, the tape-threading rows;
  the flagship reads the STORE and is blind by design); A9 ℓ=3/ℓ=4 (n,2n) moments permuted at a P2 solve
  **0 — the declared null, confirmed**; A10 the clamp reading both stacks **2** (the B-10 witness rows);
  A11 positive control (`add_p0_source` a no-op) **166**; A12 `is_isotropic` forced True **24** (4/20);
  A13 `at_order` padding with the P0 block **7** (2/5); A14 the role's P0 binding ignored **1** (the tier-2
  row). Unspellable by admission (not run): ν₂ₙ = 0, ν₂ₙ = 2.0. ⚠ harness lesson: a bite check that reads
  `inspect.getsource` of a function exec'd from a string raises `OSError` — the four A4 arms and A10 first
  reported INSTALL FAILED for that reason and were re-run reading the recorded mutant source.
  ✅ ELEGANCE REVIEW (enforcer, `scratch/_426_step2_elegance_review.md`): 3 violations + 9 should-fix + 8
  nits + 3 reshapes; taken in the same landing: V1 the kernel gate's docstring asserting the reversed ruling;
  V2 the four extraction classmethods were a 2×2 twin of ONE recipe — the very defect this step repaired (two
  mint bodies, one minting at L=0) — now `channel: ClassVar[...]` on each role and ONE `from_solver_data` /
  `from_material_xs` on each core, a role being two class constants and NO code (the AST gate tightened to
  refuse any method, its positive control gaining the classmethod shape, its population walk made
  RECURSIVE); V3 the last production docstring naming `N2NOperator.energy`; S1 the face-binding admission
  checks BOTH faces (`[M]` an 8-ordinate reconstruction face on a 4-ordinate interior was accepted and the
  windowed arm returned an `(8, 2, 4)` source — gated now); S2 the mint's refusal names `from_field`; S3 one
  `_scalar_interior_space` (the file's last bare `assert` went with it); S4 `scattering_order` →
  **`legendre_order`** on the core (a channel name on a channel-agnostic core — the (n,2n) binding's order IS
  the elastic clamp's; the `from_field`/`from_solver_data` kwarg keeps the solve's name; 26 test reads +
  1 production read re-spelled); S5 `MaterialXSField.is_p0_diagonal_with_zero_n2n` RETIRED (zero callers; its
  ⚠ note asked THIS step to re-derive — done: `rowsum(Σ₂ₙ,0) = σ₂ₙ` and `|Σ_ℓ| ≤ Σ₀` entrywise make a zero P0
  row force every ℓ block of that row to zero, so the P0 test was sufficient and the method merely dead);
  S7 `kernel` refuses on the ONE predicate `is_isotropic` (a P2-bound (n,2n) binding over an isotope with no
  MT=16 used to build an all-zero product silently); S9 the two live message strings and the docstring
  residues; R2 the field's order is uniform BY CONSTRUCTION (`__post_init__` pads every kernel to the widest
  stored order — the pad law's fourth spelling, the lazy `order` refusal, dissolved); the `block_role` nit
  REFUTED (`[M]` pyright: the base declares it an instance attribute, so a `ClassVar` override is a typing
  violation — the unannotated class-level default stays, its comment now says why); `LegendreMomentTransfer.
  from_field` → `on_basis` (two `from_field`s of different shape in one module). Deferred with a trigger:
  S8 the module rename `isotropic_scattering.py` → `isotropic_transfer.py` (`[M]` 106 sites / 37 files, 53
  in docs — rides the corpus pass; done right after the archivist returns so the sweep runs once); the
  three architectural opportunities filed as issues — **#448** (the finalize's P0-only angular-flux
  reconstruction, measurement first), **#449** (the facade fold family y-blind), **#450** (`FissionOperator`
  adopting the derived-energy-binding shape). ✅ REVIEW ROUND COMMITTED `f52877db` (every affected tree
  re-verified green: transport 737, sn/operators 1278, sn/solve 208, analytical 74, primitives 339,
  architecture 152, mms 23, diffusion 115, homogeneous 50, root gates 405).
  ⚠ SURPRISES logged (mine, process, not plan structure): (1) an edit script that asserts occurrence counts
  SEQUENTIALLY and aborts leaves every LATER edit silently unapplied — the census and k_inf-ladder edits were
  lost twice behind a material-field count mismatch, and only the per-file collected-count reconciliation
  against a HEAD worktree caught it (the ladder showed +0 where +12 was owed). ⟹ a multi-file edit script
  writes each file as it goes and reports per file, and the reconciliation is not optional; (2) the
  name-keyed rewriter's negative set again (the 2026-09-03 row): `.n2n` field reads on `IsotropicN2N`
  instances and a second `_moment_scattering` site were outside its rules and surfaced only in the red loop.
- **Step 3 — the corpus**: `adjoint.rst` (`sn-n2n-p0-truncation` → dated history; `sn-n2n-isotropic-lift`
  and `sn-n2n-adjoint-source` gain their per-ℓ form — a labelled equation is an API), `slab_multigroup.rst`
  `_n2n-reactions`, `cross_section_data.rst` (the P2 statement, the h5 schema, the nnz table; #428's split
  lands here too), the error catalogue entry, the in-code docstrings the explorer listed; `Closes #426`.

Proactive dispatches: **test-architect before step 1** (the carve crosses data ↔ transport ↔ sn);
elegance-enforcer after each step; archivist for step 3. Batteries in-process, crash-safe. The 13-tree
gate before each merge, pre-registered per tree against the 11007 baseline.

## 4. Step 2 design — the transfer family (proposed 2026-09-03 for the checkpoint; F2/F3 ruled in §2b)

**Goal (domain).** The two collision-gain terms of `A = (L+C) − S − N2N − B` are two instances of one
object — a Legendre transfer stack with a yield, bound to an angular frame at the solve's order — and the
(n,2n) term's angular distribution reaches the solve. Done-when: §1 items 2–4.

### 4.1 Class map (kernel tier `Transfer*`, operator tier role names on shared cores)

⚠ This table is the PRE-CARVE design (2026-09-03). What LANDED differs in module placement, the mint
names, the `isotropic_binding` ClassVar, `is_isotropic`, and the un-minted `net_production` — see the
step-2 record in §3 (items a–g), which is the authority.

| tier | today | step 2 | content |
|---|---|---|---|
| kernel (`transport/kernels.py`) | `ScatteringKernel(moments)`, `N2NKernel(matrix; multiplicity ClassVar 2)` | **`TransferKernel(moments, multiplicity)`** frozen; `N2N_MULTIPLICITY: Final[int] = 2` the ONE home of the channel constant (MC's alias reads it); constructors `TransferKernel.scattering(mixture)` (`SigS`, y = 1) / `.n2n(mixture)` (`Sig2`, y = 2); `ng`, `order`, `p0`, `emission_matrix()` = `y·p0.T`; **`at_order(L)`** — truncate below the stored order, PAD with zero moments above it (see 4.3) | absorbs both; the multiplicity is a DATUM |
| field (`transport/material_field.py`) | `ScatteringMaterialField`, `N2NMaterialField` (member-identical verbs, `scale=multiplicity` only on the P0 verbs) | **`TransferMaterialField(MaterialField[TransferKernel])`**: `.scattering(mat_xs)` / `.n2n(mat_xs)`; `order`, `multiplicity` (uniform by admission), `at_order(L)`; verbs `add_p0_source(+_transpose)`, `moment_source(+_transpose)`, `add_to_group_rate` — every one with `scale = multiplicity` (the `if scale != 1.0` fast path keeps y = 1 bit-identical) | `add_emission`/`moment_emission`/`_moment_l0` retire — the ℓ = 0-only body was the truncation |
| P0 energy binding (`operators/isotropic_scattering.py`) | `IsotropicScattering`, `IsotropicN2N` (member-for-member twins) | **`IsotropicTransfer`** core with field `transfer`; role subclasses `IsotropicScattering(IsotropicTransfer)`, `IsotropicN2N(IsotropicTransfer)` carrying ONLY `from_material_xs` (which channel) and the role name; `dense_per_material` = `y·p0.T` once | the K_iso leaves keep their term names |
| moment operator (`operators/scattering.py`) | `LegendreMomentScattering`, `N2NMomentOperator` | **`LegendreMomentTransfer`** (field `transfer`; `skip_l0`, `L` derived) | kernel-tier math, not a term |
| angular binding (`operators/scattering.py`, `operators/n2n.py`) | `ScatteringOperator`, `N2NOperator` (S's arms with `aniso=None`, frame at L = 0) | **`TransferOperator`** core (`ScatteringOperator`'s body: faces, arms, `kernel`, `full_transfer_kernel`, `apply_transpose`, `isotropic_energy`, the fold verbs generic in y); role subclasses `ScatteringOperator(TransferOperator)`, `N2NOperator(TransferOperator)` carrying ONLY `from_solver_data(*, mat_xs, scattering_order, space)` — **N2N gains the order** and mints the SAME interned frame | `n2n.py` shrinks to the role subclass + its module docstring |
| gate | — | **AST role gate**: the four role subclasses define nothing but classmethods and `ClassVar`s (two filters: AST for the body, runtime `__subclasses__` for the population; a validated positive control) | F3's ⛔ |

Consumers re-spelled to ONE name: `.energy` → `.isotropic_energy` (`sn/solver.py`, `coupled_system.py:592`
`S.isotropic_energy + N2N.isotropic_energy`), `.scattering` / `.n2n` field reads → `.transfer`
(`solver.py:1601/1655/1709/2346` → `self.n2n_op.isotropic_energy.transfer.add_to_group_rate`), the
`isinstance` door at `solver.py:1214` keeps working (role classes survive). `MaterialXSField.n2n_matrix`
stays the P0 reaction matrix (removal / the fold predicate — P0 by physics). Diffusion / homogeneous keep
`IsotropicN2N.from_material_xs`. CP / MoC / MC untouched (P0 by construction).

### 4.2 Bit-identity expectations, stated per path (to be MEASURED, not argued)
- Every P0-only (n,2n) fixture (the tree's entire (n,2n) corpus, #269's residue): the padded ℓ ≥ 1 moments
  are exact zeros, so `aniso` is exactly zero and `(iso/W) + 0.0` is bit-identical; the P0 verbs take the
  same einsum with `scale` applied as today. Expect BIT-IDENTICAL; the pre-T3 snapshots and the tier-2
  equivalence rows are the pins — measure, and state per row.
- Scattering at y = 1: the `scale != 1.0` branch is skipped — bit-identical by construction of the branch.
- The flagship: `solve_sn(scattering_order=2)` on the §0 fast/thin fixture reads the A2 arm's
  `1.091199657` (`rtol` 1e-8, ruling O-3); at `scattering_order = 0` the pre-carve `1.1587120371368607`.

### 4.3 The one design question left open by F2 — `at_order` above the stored order — ✅ RULED (user, 2026-09-03): **`at_order(L)` pads**
`ScatteringKernel.truncated` REFUSES a request above the stored order ("moments beyond L do not exist and
are not invented"); ruling O-1 PADS the (n,2n) stack. Both are right about different stacks: a stack of
length 1 (NL = 1 — the evaluation declares isotropy; or an absent channel) is COMPLETE, and zeros above it
are the evaluation's statement; a stack of length 7 is GROUPR's cap (`lord = 6`), and zeros above it would be
a fabrication. The kernel cannot tell the two apart from the stack alone — but it never needs to: the SOLVER
clamps `L` to the scattering stack (`min(len(SigS) − 1)` over materials), so `at_order(L)` is only ever
asked for `L ≤` the mixture's widest scattering order, where padding a SHORTER channel is honest. ⟹
**proposed: `at_order(L)` pads; the clamp's own honesty (refuse/warn when a request exceeds the data) is
#60's, unchanged here.** Alternative: keep `truncated` refusing and give `TransferMaterialField.n2n` a
pad at construction to the mixture's scattering order — same arithmetic, the pad spelled at the field
instead of the kernel — not taken.

### 4.4 The fold verbs' home — ✅ RULED (user, 2026-09-03): **on the `TransferOperator` core, generic in y**
`foldable_part` / `residual_part` / `foldable_sigma` (the within-group self-transfer folded into removal
for the sweep, `Σ_r = Σ_t − y·Σ_{c,gg}`) are generic in y and today applied to S only (the SI splitting
reads S). Proposed: they live on the `TransferOperator` core, generic, and the splitting keeps calling them
for S — no behaviour change; folding (n,2n)'s within-group block is a later, measured decision. Alternative:
leave them on `ScatteringOperator` as role-specific methods — which breaks the "no arithmetic in a role
class" gate and keeps a twin waiting.

### 4.4b Landing — ✅ RULED (user, 2026-09-03): **step 1 merges to `main` FIRST**, then compaction, then step 2
on a fresh branch from that `main` (smaller merge units; every checkout regenerates its HDF5 store — the loader
says how). *"Also prepare for context compaction before further development."*

### 4.5 Sequencing
Step-2 opener: the deferred `_legendre_order` producer fix (§3); a fresh census of the `.energy`/`.n2n`/
`.scattering` consumer spellings (AST) and of the must-flip gates (the test-architect's list); the
test-architect's drafts `_426_draft_test_be_reflected.py`, `_426_draft_test_role_ast.py`,
`_426_draft_test_diffusion_n2n.py` land re-keyed to the final names. Then the carve as ONE landing (a
signature/type change across the five tiers cannot be split — §6b), the battery (11 arms, per arm and per
call site), the elegance review, the corpus pass (step 3 folds in: `sn-n2n-isotropic-lift` / `-adjoint-source`
gain their per-ℓ form, the ERR entry, `Closes #426`).

## 5b. ⏹ COMPACTION POINT #2 — THE CAMPAIGN IS COMPLETE (2026-09-04)

**READ THIS SECTION FIRST on pick-up. §5 below is the step-1 checkpoint (history); §3 carries every step's
landing record; §4.1 is the pre-carve design and §3's step-2 entry (items a–g) is the authority on what landed.**

### 5b.1 What landed (all on `feature/n2n-transfer-family`, ff-merged to `main` — verify with `git merge-base --is-ancestor`)

| item | commit |
|---|---|
| the MF=6 reader returns a named section carrying the header's NL / NZ (the deferred producer fix) | `7b44ee68` |
| **step 2 — the transfer family; the (n,2n) gain is anisotropic** | `1a3b78ec` |
| the elegance review round — a role is two class constants and no code; both faces admitted; `legendre_order`; `is_isotropic` the one predicate; the field at one order by construction; the dead facade predicate retired | `f52877db` |
| **step 3 — the corpus reads the family; ERR-082; `isotropic_scattering.py` → `isotropic_transfer.py`; `Closes #426`** | `9e6adf3c` |
| the fold-accessor routing sentinel follows the fold family's re-home | `ceba4f5d` |

### 5b.2 The exit, as measured
- **13-tree gate (`scratch/_426_step2_full_gate.log` + `_426_step2_sn_rerun.log` for sn on `ceba4f5d`): 13 of 13
  rc=0, 11142 collected** — numerics 3255, transport 737 (+1 skip), geometry 727, data 314, homogeneous 50,
  diffusion 115, cp 141, moc 121, mc 40, cross_method 81, sn 3415, derivations 1636 (−1 = the retired
  diagnostic's registry row), root+harness 425 (+1 = the layer gate's `transfer.py` row); the sum of
  passed + skipped + xfailed is 11142 exactly. `sphinx -E -W` from the repo root EXIT=0; `dead_references`
  0 / 65; pyright 0 on `orpheus/`; the collected count reconciled per file against a HEAD-worktree collection
  at every step.
- **The done-when (§1), all four**: (1) the ingest is lossless for every channel (step 1); (2) the Be-reflected
  slab reads **1.0911996566537725** with the shipped library and no probe (the P0 model's 1.0953221881419453 is
  the flagship's P0-only control, bit-identical); (3) the retired spellings return history prose only —
  `orpheus/` clean, `tests/` clean, `docs/` swept by the archivist; (4) ERR-082 written and caught.
- The battery table and the review round are in §3's step-2 entry; the measured ladder in
  `scratch/_426_step2_ladder_record.txt` (untracked — the numbers are ALSO in the flagship's docstring and
  ERR-082, which are tracked).

### 5b.3 What this campaign leaves behind (GitHub is the backlog; these are pointers, not a second list)
#445 (`SigT` derived three times), #446 (`transport_xs` (n,2n)-blind), #447 (the roster's σ₀ column),
**#448** (the SN finalize reconstructs the returned angular flux from a P0-only source — pre-existing,
measurement first), **#449** (the facade's σ_r fold family y-blind), **#450** (`FissionOperator` should derive
its energy binding as `TransferOperator` does), #60 (the clamp's honesty). The archivist's "weakest dimension":
no `derivations/` script for the (n,2n) Legendre stack — a SymPy `derive_*()` proving
`emission_matrix(y=2) ≡ 2·emission_matrix(y=1)` would turn `sn-n2n-transfer-binding` from a sentinel into a
foundation edge — **#451** (filed 2026-09-04; the eigenvalue-tier route is an anisotropic extension of
`derive_2g_n2n`, not an MMS).

### 5b.4 Resume pointer (per plan-authoring §1: the OUTCOME)
**This plan is DONE. The next act is Campaign 2's step 5** — `.claude/plans/cs4c_binding_design.md` §18.8
(updated 2026-09-04: step 5 UNBLOCKED, with the ⚠ that §17's N2N/S design text predates the family and must be
re-derived against `orpheus/transport/operators/transfer.py` at the opener), then §17's ▶ block. The two rulings
still owed: #429 umbrella-or-close; step 5's own.

## 5. ⏸ COMPACTION POINT #1 (2026-09-03, written pre-compaction with full context — HISTORY; §5b supersedes)

**READ THIS SECTION FIRST on pick-up; then §4 (the step-2 design, all rulings in), then §2b/§2c.**

### 5.1 What landed (all on `main` after the ff-merge; every hash verifiable with `git merge-base --is-ancestor`)

| item | commit | note |
|---|---|---|
| vv-principles #29 sharpenings the step-5 census owed | `8707c53a` | |
| the plan of record + #35/#36/#17-filter | `72aa4a59` | the measurement, F1/F2/F3 ruled |
| #428 CLOSED — the doc split, per-solver census | `f1422f24` | agent memories ride along |
| cs4c §18.8 note + F1 cost | `ef915e2d` | |
| **step 1 — the ingest is lossless in ℓ** | `f96de34c` | 13 of 13 rc=0, 11087 collected (data +77, sn +3), 61 min 33 s; sphinx −W clean; dead_references 0 / 52 |

### 5.2 Step 1's exit, as measured
- Bit-identity: 13 of 13 pure-isotope mixture digests and the fast/thin k at L = 0/1/2 reproduce to the bit;
  the one exception (NA023's ℓ = 1, 2 phantoms, 1e-30) is re-baselined in
  `tests/data/data/pre_426_ingest_ledger.json` with its reason. Elastic P3…P6 measured (§2c O-4): −229 / −209
  / −192 at P3 on the fast/thin fixture, converging to −173 / −158 / −144 by P6.
- Gates landed: `tests/data/test_n2n_yield_convention.py` (+21: tape physics, threading, one-diagonal),
  `tests/data/test_hdf5_store.py` (12), `tests/data/test_ingest_ledger.py` (44), canonical-order +1, condense +1,
  `tests/sn/solve/test_scattering_order_is_the_only_truncation.py` (3). Battery 5 arms + control, all bite.
- The HDF5 store is format 2; **every checkout must regenerate it**
  (`.venv/bin/python orpheus/data/micro_xs/convert_gxs_to_hdf5.py`, ≈ 7–8 min; the loader refuses a stale
  file with that command). Test skips that used to hide a missing store are now loud refusals.

### 5.3 What step 2 does (the transfer family) — design RULED, nothing written
§4 is the design: `TransferKernel(moments, multiplicity)` with `N2N_MULTIPLICITY = 2` as the channel constant's
one home and `at_order(L)` that pads (§4.3); `TransferMaterialField` with `scale = multiplicity` on every verb;
cores `IsotropicTransfer` / `LegendreMomentTransfer` / `TransferOperator` (the fold verbs on the core, §4.4);
thin role subclasses `ScatteringOperator`, `N2NOperator`, `IsotropicScattering`, `IsotropicN2N` with only their
extraction classmethods — `N2NOperator.from_solver_data` GAINS `scattering_order` and mints the same interned
frame; consumers re-spelled to `.isotropic_energy` / `.transfer`. Verification plan: the test-architect's
`scratch/_426_verification_plan.md` §"Step 2" (G2.1–G2.10, the 11-arm battery, the must-flip list with each
gate's ℓ ≥ 1 fixture) and its drafts `scratch/_426_draft_test_{be_reflected,role_ast,diffusion_n2n}.py`.
Step-2 OPENER owes: the deferred `_legendre_order` producer fix (§3 step 1 record); an AST census of the
`.energy` / `.n2n` / `.scattering` consumer spellings; re-derivation of the blast-radius memo's §B against the
post-step-1 tree (`scratch/_426_remedy_blast_radius.md` was written BEFORE step 1 — its data-layer half is
history now, its operator half stands). Expect the P0-only (n,2n) corpus bit-identical (§4.2) — MEASURE it.

### 5.4 Standing constraints (unchanged)
§4 of this plan's former numbering — now §6 — plus: main agent writes production, user steers at
AskUserQuestion checkpoints; test-architect BEFORE the carve (its step-2 plan exists — extend, don't
re-dispatch); elegance-enforcer after; archivist for the corpus; in-process monkeypatch batteries, one process
per arm; the 13-tree gate before the merge, pre-registered per tree; `sphinx -E -W` from the REPO ROOT +
`dead_references` at exit; never `git add -A`; `git commit -F`; L37 no source edits under a running gate.

### 5.5 Resume pointer (per plan-authoring §1: the OUTCOME, with existence-checks) — ⛔ SUPERSEDED by §5b.4 (step 2 landed 2026-09-04)
**Next act: step 2 — the (n,2n) term is anisotropic, realised as two instances of one transfer family.**
Existence-checks at resume (one grep each): `TransferKernel` (must be ABSENT — 0 hits — or step 2 has begun),
`N2NKernel.multiplicity` (present, the ClassVar to retire into `N2N_MULTIPLICITY`), `for_space(interior, 0)` in
`orpheus/transport/operators/n2n.py` (present — the truncation site step 2 removes), `Isotope.sig2`
annotated `list[csr_matrix]` (present — step 1 landed). Open a fresh branch from `main` (`feature/n2n-transfer-family`).
Campaign 2's own pointer: `.claude/plans/cs4c_binding_design.md` §18.8 — step 5 waits behind this carve.

## 6. Standing constraints

Those of `.claude/plans/cs4c_binding_design.md` §18.7, unchanged. Baseline `[M]` `main` `1e02f6b1`:
11007 collected, 13 of 13 rc=0. Branch `fix/n2n-anisotropy` carries one prior commit (`8707c53a`, the
vv-principles #29 sharpenings the step-5 census owed).
