# P2-A — `orpheus/transport/operators/scattering.py` code-prose rebalance (PILOT)

Phase 2 of #231, PILOT batch. Pre-edit: **1807 lines, 1127 docstring +
196 comment = 73% prose.** Branch `docs/sn-doc-architecture`.

## Headline pilot finding (calibrates the six follow-on batches)

**The book is exhaustively complete for this file. ZERO MOVED.** Every
teaching / derivation / design-rationale concept in `scattering.py`
already has a TWIN home in the theory corpus — nothing needed relocating.
Verified by grepping the landing chapters for every concept: the
operator-algebra content (§5.6 kernel, capability surface, Funk–Hecke,
`S=R∘Λ∘M`, `Λ=ΣP_ℓ⊗Σ_{s,ℓ}`, the foldable/σ_r-fold + #215 trap, K_iso
model-shared, the Euclidean adjoint `S^T`, the 1/W chain) is all present
across `slab_multigroup.rst`, `spherical_harmonics.rst`,
`operator_algebra.rst`, `adjoint.rst`, `slab_one_group.rst`,
`loss_representation.rst`. **Implication for the operator-file batches
(fission, streaming, boundary, multiplication, …): expect TWIN-cutting,
not MOVED-writing** — the operator-algebra book already carries the
theory. The rebalance is a docstring→pointer contraction, not a
relocation.

Two structural facts that make cutting safe here:
- `scattering.py` is **NOT automodule'd** anywhere (grep clean), so its
  `:ref:`/`:eq:`/`:class:` cross-refs already render plain-text; no
  Sphinx link breaks on a cut.
- **No `.. math:: :label:` and no `vv-status` in any docstring** → cutting
  math blocks orphans no `:eq:` target; no `verifies(...)`/`catches(...)`
  test marker points at a `scattering.py` docstring (they point at RST
  labels `pn-scatter`/`flux-moments`/`matrix-eigenvalue`).

## Pointer convention used

Literal greppable form `docs/theory/<part>/<file>.rst §<label>` (ratified
for this phase). Every cited label verified to resolve (grep gate).
Existing contract-level `:eq:`/`:class:`/`:meth:` code-refs that survive
a cut are kept as roles (unchanged).

## Posing harmonization (in-scope, text-only)

Module head opened `(L - S - F)ψ = q` with `L` = "streaming + collision"
— a dated fused-L posing. Harmonized to the ratified record
`A = L + C − S − B` (L streaming-only, C collision, S scattering gain,
B boundary; `K = A⁻¹F`), matching `notation.rst §notation-symbol-table`
row 6 verbatim. Three defects fixed: fused-L split into L+C, B added,
F moved to the RHS eigenvalue posing (not a LHS subtractive gain).

## Live-issue pointers preserved (one-line form, per brief)

`#215` (σ_r-fold DSA trap), `#2` (DSA consumer), `#205` (cross-method
field arch — the ScalarFlux keep-decision), `#118` (closed by the
adjoint). NO `#305`/`#306` exist in the current file (the brief's caution
named a `sig_s`→#306 note that is not present in this checkout; the live
`sig_s` docstring references only the retired P1.7 helpers — reported).

---

## Adjudication table (one row per block; pre-edit line ranges vs 1807-line file)

| pre-edit lines | verdict | destination / record | content id |
|---|---|---|---|
| 1–13 (mod head posing) | POSING-HARMONIZE | notation.rst §notation-symbol-table (oracle) | `(L−S−F)ψ=q` fused-L → honest `A=L+C−S−B`, `K=A⁻¹F` |
| 14–20 (S aggregates: bullets) | CONTRACT (kept, trimmed) | — | what S aggregates: P0 / Pℓ / (n,2n) — one-line each |
| 22–29 (apply-only rationale) | TWIN | operator_algebra.rst §capability-set-semantics | "apply-only on inverse axis; no tractable solve; capability-set design" |
| 31–38 (Cardinal-Rule-2 lift / Wave D bit-identical / SNSolver delegators) | HISTORY | slab_multigroup.rst §sn-scattering-fission-operators (439–444, 578–582) — verified delegators still live at solver.py:1884/1892/1921 | Wave D Issue 13 extraction narration |
| 40–60 (Math: P0 in-scatter) | TWIN | slab_multigroup.rst §mg-inscatter-source | `Q^0_g = Σ_s0^T φ` derivation |
| 61–93 (Math: Pℓ Galerkin on Y_lm) | TWIN | spherical_harmonics.rst §spherical-harmonics-eigenbasis + slab_multigroup.rst §pn-scatter/§addition-theorem | moment projection + addition theorem + (2ℓ+1) |
| 95–102 (Math: (n,2n) doubling fold) | TWIN | slab_multigroup.rst §n2n-source (incl. the "why fold into scattering" 3-reason list, 364–375) | 2·Σ_2n φ secondary emission |
| 104–117 (Capability advertisement) | CONTRACT (kept, trimmed) + TWIN pointer | operator_algebra.rst §capability-set-semantics; adjoint.rst §sn-scattering-adjoint-source | apply + apply_transpose, no solve; S^T via full_scatter_kernel |
| 119–152 (§5.6 Kernel reading) | TWIN | operator_algebra.rst §integral-kernel-category | scattering IS a nonlocal integral kernel; kernel = R∘Λ∘M |
| 204–214 (`_spatial_moments_of`) | CONTRACT (kept) | — | φ̂ width read off ScalarFlux space; ndarray→1 (#240 D5b-S3 invariant) |
| 218–229 (LegendreMomentScattering: Λ=ΣP_ℓ⊗Σ_{s,ℓ}) | TWIN | operator_algebra.rst §scattering-as-tensor-product-sum | the sum-of-tensor-products form |
| 231–252 (action math + skip_l0) | CONTRACT (kept, trimmed) | — | input shape `(L+1,2L+1,ng,*sp)`; per-material fold; skip_l0 semantics |
| 254–276 (Implementation + capability set) | TWIN (impl-narration) + CONTRACT (caps) | operator_algebra.rst §capability-set-semantics; spherical_harmonics §…eigenbasis | einsum/loop narration cut; `{apply,apply_transpose}` kept; Λ^T kept |
| 278–291 (Parameters) | CONTRACT (kept, trimmed) | — | mat_xs / L / skip_l0 |
| 313–357 (apply: role-changing edge) | CONTRACT (kept) + TWIN teaching | operator_algebra.rst §integral-kernel-category (role-edge note ~3007) | typed-in→typed-source-out; ndarray→ndarray; m-axis shift; contract kept |
| 375–400 (apply_transpose Λ^T) | CONTRACT (kept, trimmed) + TWIN | adjoint.rst §sn-scattering-adjoint-source | Euclidean transpose, role-reversing; `(RΛM)^T` fall-out is TWIN |
| 415–434 (domain/codomain) | CONTRACT (kept, trimmed) | — | endomorphic on SH space; the "cast that papered over None" teaching cut |
| 438–458 (N2NMomentOperator class) | CONTRACT (kept, trimmed) + TWIN | adjoint.rst §sn-scattering-adjoint-source (92–93 names N2NMomentOperator) | distinct (n,2n) multiplication channel; endomorphic ℓ=0; pointer |
| 463–493 (N2N methods + domain/codomain) | CONTRACT (kept) | — | short — apply/apply_transpose/domain |
| 498–531 (ScatteringOperator class) | CONTRACT (kept, trimmed) | slab_multigroup §sn-scattering-fission-operators (Wave-D history cut) | attributes: mat_xs/quadrature/scattering_order; caps line |
| 537–541 (block_role comment) | COMMENT-keep (constraint) | — | BULK operator, no boundary action (#208) |
| 550–559 (full_field_space field comment) | CONTRACT (kept, trimmed) | — | composite space threading; None default; #261 forward-note trimmed |
| 566–594 (domain/codomain/is_adjointable) | CONTRACT (kept, trimmed) | operator_algebra §carrier-grid (BULK 2-cell teaching cut) | advertises full_field_space; endomorphic; None-default contract |
| 604–638 (read-through props: n_ordinates/spatial_shape/ng/weights/Y) | CONTRACT (kept) | — | short accessors; shapes |
| 640–677 (`frame` property) | CONTRACT (kept, trimmed) + TWIN | operator_algebra.rst §… (harmonic-frame-is-galerkin ~2983) | frame = HarmonicFrame(order L), M/R faces; from_galerkin bit-id narration cut |
| 680–716 (sig_s/sig2/sig_s0/cells_by_mat TRANSIENT) | CONTRACT (kept, trimmed) | — | TRANSIENT read-throughs; retired-consumer history trimmed (NB: no #306 here) |
| 720–774 (`_aniso_source_from_moment_values`) | CONTRACT (kept, trimmed) + TWIN | operator_algebra.rst §integral-kernel-category (explicit-vs-fused ~3020) | RΛ map, pre-1/W, shape contract; the "supersedes _PerLegendre" teaching cut |
| 781–875 (`kernel` property) | TWIN | operator_algebra.rst §integral-kernel-category | §5.6 semantic reading — cut to contract (returns R∘Λ∘M, raises on order 0) |
| 877–910 (`full_scatter_kernel`) | CONTRACT (kept, trimmed) + TWIN | adjoint.rst §sn-scattering-adjoint-source (80–103) | returns R∘(Λ+N2n)∘M; used by apply_transpose; teaching cut |
| 912–952 (`isotropic_kernel` K_iso) | CONTRACT (kept, trimmed) + TWIN | adjoint.rst §…-adjoint-source (154–164) + infinite_medium.rst | Σ_s0+2Σ_2n on scalar flux; cached; space=None; model-shared teaching cut |
| 954–1004 (`_assemble_per_ordinate_source`) | CONTRACT (kept, trimmed) + TWIN | slab_multigroup.rst §…normalization-chain (377–410) | the 1/W combine SSOT; params kept; Pattern-7 teaching trimmed |
| 1006–1037 (`from_solver_data`) | CONTRACT (kept, trimmed) | — | constructor; full_field_space contract; PR-TYPED-1 history trimmed |
| 1041–1107 (`add_iso_source`) | CONTRACT (kept, trimmed) | — | raw-in-mutates vs typed-in-returns-new — critical caller semantics |
| 1109–1162 (`add_n2n_source`) | CONTRACT (kept, trimmed) | — | same typed-action overload |
| 1164–1275 (`build_aniso_source`) | CONTRACT (kept, trimmed) + TWIN | slab_multigroup.rst §pn-scatter-rlm + spherical_harmonics | (1/W)·RΛM; params/returns/None-cases kept; the M/R/Λ derivation cut |
| 1277–1328 (foldable/residual comment block + #215 trap) | COMMENT-cut (constraint kept) | slab_one_group.rst §si-sigma-r-fold-mismatch + loss_representation.rst §loss-rep-removal-sigma | split narration cut; #215/#2 one-line constraint kept |
| 1330–1443 (foldable_part/residual_part/is_foldable/foldable_sigma) | CONTRACT (kept, trimmed) | slab_one_group.rst §si-sigma-r-fold-mismatch | what each returns + the additive contract; teaching trimmed |
| 1447–1510 (`_apply_impl` dispatch table) | CONTRACT (kept, trimmed) | — | the 4-carrier dispatch surface a modifier needs; per-arm math cut |
| 1512–1578 (FullField register arm) | CONTRACT (kept, trimmed) + TWIN | operator_algebra §…; coupled_block_operator.rst (route-a) | bulk-only, implicit-zero boundary; math teaching cut; #282/#208 one-line |
| 1580–1611 (ScalarFlux arm) | CONTRACT (kept) | — | the W-F keep-decision (named-future #205 consumer) — KEEP verbatim intent |
| 1613–1642 (AngularFlux arm) | CONTRACT (kept, trimmed) | adjoint.rst §…-adjoint-source (fast-path) | the iso fast-path PERF constraint kept; math teaching trimmed |
| 1644–1703 (HarmonicMomentFlux windowed arm) | CONTRACT (kept, trimmed) + TWIN | operator_algebra.rst §… (explicit-vs-fused windowed ~3033) | moments-are-Mψ contract kept; bit-identity teaching cut |
| 1705–1725 (TYPE_CHECKING overload comment) | COMMENT-keep (constraint) | — | why apply=_apply_impl + overloads (#257 S8c) — kept, terse |
| 1731–1785 (`apply_transpose`) | CONTRACT (kept, trimmed) + TWIN | adjoint.rst §sn-scattering-adjoint-source | S^T=(1/W)fsk^T; Euclidean≠.H (L12); FullField arm; derivation cut |

---

## Summary (finalized post-edit)

### Line counts (measured, `ast` + `tokenize`)

| metric | before (HEAD) | after | Δ |
|---|---|---|---|
| total lines | 1807 | 1326 | −481 (−27%) |
| **docstring lines** | **1127** | **721** | **−406 (−36%)** |
| **comment lines** | **196** | **121** | **−75 (−38%)** |
| prose share | 73% | 63% | −10 pts |
| **code lines** | **484** | **484** | **0** |

`git diff --stat`: 418 insertions, 899 deletions, 1 file.

### Verdict counts (primary verdict per adjudicated block; 44 blocks)

- **POSING-HARMONIZE**: 1 (module head `(L−S−F)ψ=q` → `A=L+C−S−B`, `K=A⁻¹F`).
- **TWIN (block was pure teaching → cut to pointer)**: 8 (apply-only rationale;
  P0/Pℓ/(n,2n) math; §5.6 Kernel reading; Λ=ΣP⊗Σ tensor-sum; the
  implementation-narration; the `kernel` property).
- **HISTORY (cut; record owns it)**: 1 primary block (module-head Cardinal-Rule-2 /
  Wave-D-extraction narration) + HISTORY-cuts folded into 3 CONTRACT rows
  (ScatteringOperator class, `sig_s` TRANSIENT, `from_solver_data`). All verified
  already in the record (slab_multigroup.rst 439–444 / 578–582; delegators live at
  solver.py:1884/1892/1921).
- **CONTRACT (kept; ~28 trimmed of embedded TWIN/HISTORY)**: every public/consumer
  surface — the typed-action overloads, dispatch table, `apply_transpose`,
  domain/codomain, the foldable methods, the read-throughs, `full_scatter_kernel`,
  `isotropic_kernel`, `build_aniso_source`, etc.
- **COMMENT-cut**: 1 (the 51-line foldable/residual block → constraint + pointers).
- **COMMENT-keep (constraint)**: 2 (`block_role = BULK`; the TYPE_CHECKING overload note).
- **MOVED**: **0** (headline — see below).

### Gates

- `import orpheus.transport.operators.scattering` → **OK**.
- `pytest --collect-only -q` → **6652 tests collected** (unchanged).
- **Code-token invariance vs HEAD**: 2397 == 2397 non-string/non-comment tokens,
  IDENTICAL — proves zero behavioral/code change.
- All 9 pointer FILES resolve; all 14 cited § labels resolve (as `.. _` anchor OR
  `:label:` eq-label); the only `:ref:`/`:eq:` roles in added lines
  (`:eq:`flux-moments``, `:eq:`pn-scatter``) are retained-from-original and resolve.
- **Sphinx build NOT run** — zero theory-page edits (zero MOVED), and
  `scattering.py` is not `automodule`'d anywhere, so its docstrings are invisible
  to the build; the `-E -W` baseline (1 warning) is definitionally unchanged.

### The 5 hardest judgment calls (calibrate the six follow-on batches)

1. **ZERO MOVED — resist the "this rationale must be book-worthy-but-unwritten"
   instinct.** On a 73%-prose operator file the reflex is that some design rationale
   is novel and needs relocating. Three concepts FELT unique — the
   forward-fast-path-vs-adjoint-frame asymmetry, N2N-as-a-distinct-moment-operator,
   the foldable/σ_r split — and each was **fully TWIN** after grepping
   `adjoint.rst` / `operator_algebra.rst` / `slab_one_group.rst` + `loss_representation.rst`.
   **Discriminator: grep the landing chapter for the concept BEFORE assuming it is
   novel.** For the operator-file batches (fission, streaming, boundary,
   multiplication, N2N, isotropic_scattering), budget for TWIN-cutting, not
   MOVED-writing — the operator-algebra book is the SSOT and it is complete.

2. **A keep-vs-retire decision on a currently-orphan symbol is CONTRACT, not
   HISTORY** — even when phrased historically. The ScalarFlux `apply` arm reads like
   campaign narration ("Deliberately retained W-F 2026-06-26, user steered keep") but
   is a live constraint: it is an *intentional* orphan kept as the #205 future-consumer
   entry-point, and a naive retirement audit would delete it as dead. The keep-decision
   MUST stay at point of use. Kept the #205 rationale + the "no current caller / no
   sibling delegates" fact; cut only the date/steering provenance.

3. **Teaching → pointer, but a LATENT-TRAP imperative stays inline.** The 51-line
   foldable/residual comment is TWIN teaching EXCEPT the imperative "do NOT wire the
   σ_r-sweep as `A_wg.inverse()` (46–56 % silent error on anisotropic flux)". That is a
   constraint-stating anti-pattern guard a future within-group-solver editor reads *in
   this file*, not in the theory page. Kept the ⚠ imperative + the number + the
   #2/#215 pointers; cut the derivation to `§si-sigma-r-fold-mismatch`. Rule for the
   batches: a `⚠ LATENT … TRAP` / "do NOT" comment is COMMENT-keep even when its
   *explanation* is TWIN.

4. **A type-annotation-choice rationale that guards a plausible wrong "simplification"
   is CONTRACT.** `isotropic_kernel` returns `OperatorSum` (not the bare
   `LinearOperator` erasure) *so that* `apply_transpose` stays visible to the checker
   for the seed-arm consumer — a modifier who "tidies" the return type to
   `LinearOperator` silently breaks that consumer. Kept the typing rationale, cut the
   model-shared/bit-identical teaching. Same pattern on the TYPE_CHECKING `@overload`
   block (kept: it explains why `ScatteringOperator` is not the mixin's nominal
   endomorphism).

5. **Pointer form: cite the nearest LABELLED target — an eq-label is fine when no
   co-named section anchor exists; never invent.** Several landing concepts are anchored
   only by a `.. math:: :label:` (`mg-inscatter-source`, `sn-scattering-adjoint-source`,
   `loss-rep-removal-sigma`), not a `.. _` section anchor. I cited these as
   `§<eq-label>` (greppable, resolves) rather than inventing a section anchor or pointing
   at a coarser enclosing title. Because `scattering.py` is not `automodule`'d, ALL these
   pointers are plain-text greppable literals (the `§` is a human marker, not a rendered
   role) — which is exactly why the literal-path convention is the right form for
   operator-file docstrings.

### Reported discrepancy (not a defect)

The brief's caution named a `sig_s`-shim consumer note "pointing at issue #306". **No
`#305`/`#306` reference exists in `scattering.py` in this checkout** — the live `sig_s`
docstring referenced only the retired P1.7 `_build_rhs_*` helpers (cut to a terse
"pending retirement" marker). The live-issue one-line keeps that DO exist and were
preserved: `#215`, `#2`, `#205`, `#118`.
