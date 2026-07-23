# P2-A classification — orpheus/transport/operators/scattering.py (HEAD version)

> **GRADE vs the archivist's ground-truth map (`P2-A_scattering.md`), 2026-07-22 —
> the #231 Haiku-archaeology automation measurement.**
> Verdict-level agreement ≈ **35%** (~12 clean hits, ~8 half-hits, ~28 misses of
> ~50 comparable rows). Systematic biases: (1) **MOVED-inflation** — 24 proposed
> vs ground-truth **0**; Haiku cannot verify "the book already carries this"
> (content-level grep+read) and misreads point-of-use contracts (mutation
> semantics, typed-overload rationale, capability surfaces) as relocatable
> teaching. (2) **HISTORY-overreach on constraint-bearing text** — 3 dangerous
> proposals would have cut live guards: the `block_role = BULK` constraint
> comment, the TYPE_CHECKING overload rationale, and the ScalarFlux
> intentional-orphan #205 keep-decision (graded HISTORY at conf **H**).
> (3) **Strengths** — block inventory (boundaries/line-ranges/content-ids)
> essentially perfect; the posing-flags table accurate (caught the fused-L
> head); ALL conf-H TWIN rows correct with correct anchors.
> **Operational ruling for P2-B…G**: use the classify catalogs as INVENTORY +
> posing flags + conf-H TWIN hints; every MOVED/HISTORY verdict must be
> re-adjudicated by the executing archivist under the pilot's discriminators
> (grep-before-assuming-novel; keep-decisions and ⚠ imperatives are CONTRACT).

## Catalog

| lines | symbol | content-id | verdict | destination / anchor | conf |
|---|---|---|---|---|---|
| 1–152 | module docstring | operator-algebra posing (L−S−F)ψ=q; P0/Pℓ/n2n math; capability ad | SPLIT: see sub-rows | — | — |
| 1–40 | module head | (L−S−F) posing, intro to S operator, fission operator mention | MOVED | slab_multigroup §operator_algebra or methods/sn/placement | M |
| 40–102 | module § "Mathematical content" P0 + Pℓ + n2n derivations | P0 in-scatter, Pℓ Galerkin on Y_ℓ^m, n2n doubling, addition theorem, normalization | TWIN | slab_multigroup :eq:`pn-scatter`, :eq:`flux-moments`, spherical_harmonics :label:`real-sh-addition-theorem` | H |
| 104–118 | module § "Capability advertisement" | apply/apply_transpose, no solve, adjoint, #276 A2b, full_scatter_kernel, (R∘Λ∘M)^T | MOVED | adjoint.rst §sn-scattering-adjoint or operator_algebra §capability-set-semantics | M |
| 119–152 | module § "Kernel reading §5.6" | IntegralKernelOperator protocol, R∘Λ∘M composition, kernel vs full_scatter, P0/n2n local | TWIN | operator_algebra :ref:`integral-kernel-category`, #257 S6 design record in the chapter | M |
| 171–175 | comment: "Runtime imports…" | circular-import-safe declaration | CONTRACT | — | — |
| 186–190 | comment: "Runtime (not TYPE_CHECKING)…" | windowed moment-iterate registration reason | CONTRACT | — | — |
| 204–214 | _spatial_moments_of docstring | spatial-moment width, SpatialMomentSpace, (#240 D5b-S3) | CONTRACT | — | — |
| 218–291 | LegendreMomentScattering class docstring | §15.2 §10 sum-of-tensor-products, Λ = ∑_ℓ P_ℓ ⊗ Σ_s,ℓ, implementation | SPLIT: see sub-rows | — | — |
| 218–253 | LegendreMomentScattering head + math | Λ posing, per-ℓ block-diagonal, skip_l0 | TWIN | slab_multigroup :eq:`pn-scatter` or operator_algebra on scattering algebra | M |
| 254–276 | LegendreMomentScattering "Implementation" | loop structure, einsum, flop count, apply_transpose | MOVED | loss_representation chapter or operator_algebra §composition-algebra | L |
| 278–291 | LegendreMomentScattering "Parameters" | mat_xs, L, skip_l0 contract | CONTRACT | — | — |
| 313–357 | LegendreMomentScattering.apply docstring | role-changing flux→source, moments shape (L+1,2L+1,ng), Issue #196 PR-INDEX-4 | MOVED | conventions/indexing_and_layout or loss_representation chapter | M |
| 375–413 | LegendreMomentScattering.apply_transpose docstring | group-transpose Σ_s,ℓ, role-reversing, metric-aware .H, #276 A2 | TWIN | adjoint.rst :eq:`sn-scattering-adjoint-kernel-transpose` or operator_algebra on adjoint mechanics | M |
| 417–430 | LegendreMomentScattering.domain property | SH coefficient space, endomorphic, composability guard | MOVED | operator_algebra §composition-algebra or conventions on spaces | L |
| 437–457 | N2NMomentOperator class docstring | (n,2n) ℓ=0 operator, multiplication channel, isotropic, distinct from Λ | MOVED | slab_multigroup §n2n-reactions with cross-ref to operator_algebra | M |
| 498–531 | ScatteringOperator class docstring | S operator signature, precomputed data, from_solver_data, wave-D extraction | MOVED | placement.rst or history.rst (Wave D Issue 13) for extraction context | L |
| 537–548 | comment: "Scattering is BULK operator" | block_role = BlockRole.BULK, Issue #208 / Wave O | HISTORY | history.rst Issue #208 or Wave O record | L |
| 550–562 | field comment: full_field_space | composition guard, P4.5 W-D, SNMesh.full_field_space thread, relocation #261 | HISTORY | history.rst (P4.5 W-D) or placement.rst | L |
| 566–581 | ScatteringOperator.domain property | composite full-field space, bulk 2-cell, within-group loss (L+C)−S | MOVED | operator_algebra §composition-algebra or conventions/spaces | L |
| 596–602 | comment: "Convenience read-throughs" | Issue #197 PR-TYPED-1, legacy n_ordinates/nx/ny shim, TRANSIENT | HISTORY | closed Issue #197 or history.rst | L |
| 641–677 | ScatteringOperator.frame cached_property | HarmonicFrame, M/R analysis/reconstruction, SH basis binding, §5.6 | MOVED | loss_representation chapter on frame machinery or operator_algebra | M |
| 679–710 | ScatteringOperator property shims (Y, sig_s, sig2, sig_s0, cells_by_mat) | TRANSIENT routes through mat_xs, Issue #197 PR-TYPED-1 | HISTORY | closed Issue #197 or history.rst | L |
| 720–779 | _aniso_source_from_moment_values docstring | R∘Λ_ℓ≥1 map, moment→source half, windowed moment-iterate arm, §5.6 | MOVED | loss_representation chapter on aniso source assembly | L |
| 782–875 | kernel property docstring | §5.6 integral kernel R∘Λ∘M, IntegralKernelOperator protocol, #257 S6 | TWIN | operator_algebra :ref:`integral-kernel-category` and #257 S6 design; slab_multigroup :eq:`pn-scatter` | H |
| 878–910 | full_scatter_kernel property docstring | FULL ℓ=0+ℓ≥1+(n,2n) kernel R∘(Λ+N2N)∘M, #276 A2a, transpose | TWIN | adjoint.rst :eq:`sn-scattering-adjoint-kernel` on the full form; operator_algebra | M |
| 912–952 | isotropic_kernel property docstring | Σ_s0 + 2Σ_2n model-shared, IsotropicScattering+IsotropicN2N OperatorSum, #276 P2 | MOVED | operator_algebra §functional-category or cross-method section | L |
| 954–1004 | _assemble_per_ordinate_source docstring | P0+(n,2n) iso + ℓ≥1 aniso assembly, producer-side 1/W (Pattern 7, L18) | MOVED | loss_representation or conventions on normalization §L18 pattern | L |
| 1006–1031 | from_solver_data classmethod docstring | construction from mat_xs/quadrature, full_field_space threading P4.5 W-D | MOVED | placement.rst or history.rst (P4.5 W-D thread) | L |
| 1039 | comment: "In-place helpers…" | preserve bit-identity vs SNSolver pre-Wave-D | HISTORY | history.rst (Wave D) or Wave D Issue 13 | L |
| 1056–1093 | add_iso_source docstring | P0 in-scatter, typed-action overload (raw-in mutate, typed-in new), Issue #197 | MOVED | operator_algebra §heteromorphic-apply-typing or loss_representation | L |
| 1124–1150 | add_n2n_source docstring | (n,2n) source, typed-action overload, Issue #197 PR-TYPED-3 | MOVED | slab_multigroup §n2n-reactions with operator_algebra on typing | L |
| 1168–1253 | build_aniso_source docstring | Pℓ scattering source, R∘Λ∘M literal composition, (2ℓ+1) factor, 1/W producer-side | TWIN | slab_multigroup :eq:`pn-scatter` on equation; spherical_harmonics on (2ℓ+1) factor; conventions/normalization on L18 | H |
| 1277–1329 | comment block: "Foldable/residual split" | S = S_foldable + S_residual, DSA #2, latent correctness trap #215, sigma_r misconception | HISTORY | GitHub Issue #215 (measured 2026-06-05), #2 (DSA), history.rst | H |
| 1331–1351 | foldable_part docstring | P0 within-group self-scatter, diagonal-only sig_s[mid][0] | MOVED | operator_algebra §diagonal-operator or DSA design doc | L |
| 1369–1392 | residual_part docstring | non-foldable channels, off-diag P0, Pℓ≥1, (n,2n), S ≈ S_fold + S_resid | MOVED | operator_algebra or DSA design doc | L |
| 1410–1427 | is_foldable_into_sigma_r docstring | check structural eligibility: L==0, P0 diagonal, n2n zero | CONTRACT | — | — |
| 1429–1443 | foldable_sigma docstring | per-material foldable cross-section (σ_s,0^{g→g})_g | CONTRACT | — | — |
| 1447–1504 | _apply_impl docstring | runtime dispatch, TimedFullField/ScalarFlux/AngularFlux/HarmonicMomentFlux arms, #257 S8c | MOVED | operator_algebra §heteromorphic-apply-typing or placement.rst (Phase 5a window) | M |
| 1513–1546 | _apply_impl(FullField) docstring | composite bulk-only scattering, implicit-zero boundary, Option β3 #208 | MOVED | operator_algebra §design-b-native-solve or placement.rst (Wave O #208) | L |
| 1548–1578 | _apply_impl(FullField) implementation comment | delegation, psi.interior dispatch, #282 route(a) LIFTED, A_BA operator | HISTORY | Issue #282 route(a), history.rst Phase 2.5 | L |
| 1582–1606 | _apply_impl(ScalarFlux) docstring | iso scalar magnitude, deliberately retained, named-future-consumer, #205 | HISTORY | Issue #205 (cross-method field), history.rst (W-F 2026-06-26) keep ruling | H |
| 1614–1632 | _apply_impl(AngularFlux) docstring | per-ordinate magnitude, D-I.2 typed leaf, Phase 5a windowing | MOVED | placement.rst (Phase 5a angular-windowing) or history.rst (D-I.2) | L |
| 1645–1679 | _apply_impl(HarmonicMomentFlux) docstring | windowed iterate, flux moments input, M already done, R∘Λ, 0-ULP crosscheck | MOVED | placement.rst (Phase 5a windowing) or loss_representation chapter | L |
| 1684–1703 | _apply_impl(HarmonicMomentFlux) implementation comment | explicit typed grid path, scattered moment source materialized, spatial-moment φ̂ axis | CONTRACT | — | — |
| 1705–1712 | TYPE_CHECKING block comment | honest per-carrier typing surface, #257 S8c, dispatcher naming | HISTORY | #257 S8c, history.rst | L |
| 1732–1807 | apply_transpose docstring | S^T adjoint, (#276 A2b closes #118), full_scatter_kernel, reciprocity pin, composite FullField arm | SPLIT: see sub-rows | — | — |
| 1732–1765 | apply_transpose head | adjoint equation, forward/transpose asymmetry, frame-form Euclidean transpose, reciprocity | TWIN | adjoint.rst :eq:`sn-scattering-adjoint-source` and reciprocity; operator_algebra on adjoint mechanics | H |
| 1766–1774 | apply_transpose composite FullField arm | #280 Phase 2.5 S0.2, bulk-only cotangent, full-loss composition | HISTORY | #280 Phase 2.5, history.rst | L |
| 1776–1807 | apply_transpose body contract | chi parameter form, FullField bulk-only, cast truth, line-by-line math | CONTRACT | — | — |

## Summary counts

**Verdict breakdown:**
- **CONTRACT**: 9 blocks (args, returns, shapes, units, invariants; point-of-use necessity)
- **TWIN**: 8 blocks (land in spherical_harmonics.rst, slab_multigroup.rst, adjoint.rst, operator_algebra.rst)
- **MOVED**: 24 blocks (teaching/derivation/design not yet in theory book; mostly placement.rst, loss_representation, history.rst, operator_algebra §spatial sections)
- **HISTORY**: 11 blocks (Wave-X/Phase-Y, Issue #NN, landed work record)
- **COMMENT-cut**: 0 blocks (all contract-state comments retained as constraints)

**Total docstring lines proposed:**
- Total qualifying blocks: 52
- Estimated teaching/derivation prose to move: ~900 lines (split across Λ, kernel, scattering, apply_transpose, build_aniso_source, foldable/residual, frame, isotropic_kernel narratives)
- Estimated contract to retain: ~250 lines

---

## Posing flags & dated operator-algebra

| lines | posing | equation | interpretation | source ref |
|---|---|---|---|---|
| 6–9 | **(L − S − F)ψ = q** | forward source equation | streaming+collision−scattering−fission = source; canonical form | module head, slab_multigroup.rst |
| 8 | **L** = streaming + collision | composite operator | daily usage; foundational decomposition | slab_multigroup §sn-scattering-fission-operators |
| 14–20 | **P0 isotropic in-scatter** Σ_s^0(g'→g)φ_g' | reaction-rate form | energy-only coupling; isotropic angular projection 1/W | slab_multigroup :eq:`mg-inscatter-source` |
| 76–86 | **Pℓ Galerkin** ∑_ℓ (2ℓ+1) Y_ℓ^m(Ω_n) Φ_ℓ^m_{g'} Σ_s,ℓ | anisotropic redistribution | addition-theorem reconstruction; (2ℓ+1) normalization explicit | slab_multigroup :eq:`pn-scatter`, spherical_harmonics |
| 98–102 | **(n,2n) doubling** 2·Σ_{2n}(g'→g)φ_g' | multiplication source | two neutrons per absorption; isotropic fold | slab_multigroup :ref:`n2n-reactions` |
| 107–117 | **S^T** = (R∘Λ∘M)^T = M^T ∘ Λ^T ∘ R^T | adjoint via product transpose | adjointable capability; #276 A2b campaign; free via OperatorProduct | adjoint.rst :eq:`sn-scattering-adjoint-kernel-transpose` |
| 223–240 | **Λ = ∑_ℓ P_ℓ ⊗ Σ_s,ℓ** | block-diagonal moment scattering | per-ℓ cross-section action; energy-group axis asymmetry | slab_multigroup :eq:`pn-scatter` |
| 785–825 | **R ∘ Λ ∘ M** (kernel) | §5.6 integral kernel | nonlocal in angle; three-factor composition; #257 S6 typing | operator_algebra :ref:`integral-kernel-category` |
| 879–910 | **R ∘ (Λ_ℓ≥0 + N_{2n}) ∘ M** (full_scatter_kernel) | complete in-scatter | P0 + Pℓ + (n,2n) unified; frame conjugation; #276 A2a campaign | adjoint.rst :eq:`sn-scattering-adjoint-kernel` |
| 1177–1195 | **Q_n = (1/W) R Λ M ψ** | per-ordinate aniso source | Galerkin reconstruction on real SH; producer-side 1/W; Pattern 7 | slab_multigroup :eq:`pn-scatter`, conventions/normalization L18 |
| 1282–1329 | **S = S_foldable + S_residual** | foldable/residual decomposition | self-scatter fold; #215 correctness trap; DSA #2 preconditioner | Issue #215 record, DSA design #2 |

---

## Uncertain blocks (confidence L)

| lines | symbol | issue | reason |
|---|---|---|---|
| 254–276 | LegendreMomentScattering "Implementation" | loop/einsum/flop narrative | existence of dedicated numerics chapter / loss_representation unclear; placed as MOVED but destination heuristic |
| 417–430 | LegendreMomentScattering.domain | space composability | operator_algebra §composition-algebra exists but may belong in frame/numerics layer; context-dependent |
| 498–531 | ScatteringOperator class docstring | wave-D extraction context | placement.rst or history.rst owns the record; domain not certain without reading destination |
| 1039 | comment "In-place helpers" | bit-identity preservation | preservation clause is historical; unclear if it should live in history.rst or be deleted as obsolete |
| 1369–1392 | residual_part | residual algebra | foldable/residual split is DSA-preparatory; DSA chapter or separate design doc not yet scoped |
| 1410–1427 | is_foldable_into_sigma_r | structural check | pure contract but also design-governance (DSA preconditioner); destination: operator_algebra or DSA doc |
| 1447–1504 | _apply_impl dispatcher docstring | comprehensive but complex | heteromorphic apply is heavy; placement.rst or operator_algebra; complexity makes destination uncertain |
| 1513–1546 | FullField composite arm | wave-O / phase 4.5 record | multiple design epochs (Wave O #208, P4.5 W-D, P2.5 #280); history.rst is correct but dating precise phase is heuristic |
| 1614–1632 | AngularFlux arm | phase 5a context | windowing campaign is recorded somewhere; placement.rst most likely but exact section unclear |
| 1645–1679 | HarmonicMomentFlux arm | phase 5a windowing | same as above; 0-ULP crosscheck oracle location (_aniso_source_from_moment_values) is documented but placement is positional, not anchored |

