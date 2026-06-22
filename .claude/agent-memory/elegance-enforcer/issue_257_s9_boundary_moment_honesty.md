---
name: issue-257-s9-boundary-moment-honesty
description: "#257 S9 — production MMS boundary-source honesty + LD coherent-promise gate. PASS-WITH-NITS verdict on the SN2DCartesianLDStressMMSCase.prescribed_inflow moment-emission carve + the promoted test_ld_2d_boundary_promise.py gate. Verification-centric stage; design settled (no field type #263, no value gate). Branch feature/field-typed-operator-algebra, UNCOMMITTED."
metadata:
  type: project
---

# #257 S9 — LD boundary moment honesty + the coherent-promise gate (PASS-WITH-NITS)

Verdict on the final BEHAVIORAL stage of #257 (branch `feature/field-typed-operator-algebra`,
UNCOMMITTED). Verification-centric: makes the MMS case `prescribed_inflow` EMIT the moment-resolved
face slot (closing a #251 producer-blindness — the consumer already threaded slot-1, the producer
still emitted a scalar). Settled scope: NO new field type (→#263 collocation trigger), NO value gate
(boundary transverse-slope verified SUB-FLOOR), Deliverable-4 named accessor SKIPPED. The 3rd Mode-10
sub-floor-companion-unavailable instance (#240 D5b-S4 → #247 Leg A → #251 Leg B → S9).

## ⭐ LOAD-BEARING RULINGS (durable)

**1. `face_moment_count` gating is dispatch-on-structural-fact, NOT a boolean-flag anti-pattern.**
`prescribed_inflow` branches `if n_face_moments == 1 (DD/Step scalar) else (LD moment slot)` where
`n_face_moments = face_moment_count(sn_mesh.scheme.spatial_basis_per_axis, sn_mesh.ndim)` — the SAME
single-source primitive (`moment_layout.py:64`) the trace producer `SNMesh.boundary_face_layout` and
the interior cochain `_LossRepresentation._n_face_moments` key on. NOT anti-#3: discriminator derived
from the mesh's own discretization (not a caller-threaded flag), arms produce structurally distinct
SHAPES `(N,ng,n_t)` vs `(N,ng,n_t,n_face_moments)`, producer's downstream branch keys on the same rank.
STANDING TELL: a moment-vs-scalar branch keyed on `face_moment_count`/`face_moment_tail`/the rank
predicates is single-sourced dispatch — verify it reads the layout primitive, then PASS it; the smell
would be a SECOND spelling of "is this mesh moment-resolved?" introduced on one side.

**2. DD/Step arm byte-identical BY CONSTRUCTION — verify the `_drivers` call args, not just the closeout.**
Old `face_coords={"xmin":(np.array([0.0]),cy),...}`→`_drivers(xf,yf,g)`. New `face_specs` carries
`(const_axis,const_val,t_edges,t_centres)`; DD arm calls `_drivers(np.array([const_val]),t_centres,g)`
for x-faces. For xmin both reduce to `_drivers(np.array([0.0]),cy,g)` — IDENTICAL. Closeout pins
`np.array_equal(case.prescribed_inflow(DD).values, <pre-S9 build>)==True`. The refactor changed only
HOW the const/transverse axes are passed (to share the LD projector's signature), not the arithmetic.

**3. The L11 split (production projector ⊥ test-side reference) is CORRECT, NOT a duplication VIOLATION
— and collapsing it would DEFEAT the gate.** `_project_inflow_to_face_moments` (case-owned, `sn.py:1646`)
descends ONLY from `self._drivers` + `numpy leggauss` + `AVERAGE_MOMENT` (grep-VERIFIED: `_inflow_to_moments`/
`_ubld`/LD-ops/test-projectors appear ONLY in the docstring as negations). Test-side `_face_transverse_legendre`
(`test_mms_ld_2d.py:1050`) is the structurally-distinct leggauss hand-projection. GATE C
(`test_case_projector_agrees_with_test_face_projector:1152`) pins their agreement at machine precision.
⭐ If production imported the test projector (or vice versa) GATE C passes TAUTOLOGICALLY and a
double-applied transverse mass `diag(h_t,θ·h_t)` normalization drift slips through silently — that IS
the TRUE bug the bare-coeff convention (NO θ/h_t at the projector; the cochain mass adds them downstream)
defends. The split keeps the teeth. STANDING TELL: when a verification gate pins "production projector ==
hand reference", the two MUST be independent implementations; an import that single-sources them is the
laundering habitat, not the elegant move. The bulk `_project_scalar_to_tensor_legendre` is TEST-ONLY
(Leg A's reference) — there is NO production bulk projector to extract toward (Leg A consumes the
moment-resolved external source projected at solve time), so "parallel-but-independent" is forced here.

**4. The `_solve_with_boundary_slope` re-baseline is a principled test-mechanics fix.** Before S9
`slope_sign=None` meant average-only by routing through `build_nonvacuum_fixed_source`→
`case.prescribed_inflow`; after S9 production emits the real slope, so `None` would SILENTLY become
"average+real slope", collapsing the flat/mom/flip toggle + breaking the `phi_zero==phi_today`
byte-identity no-op control. FIX: build ALL branches TEST-SIDE from `_face_transverse_buffers`
(`None`→scalar `(N,ng,n_t)`; `±1/0`→moment slot). STANDING TELL: a verification HELPER that silently
inherits the very production behaviour it is the controlled baseline for = measurement-laundering
habitat (the null result goes vacuous unnoticed). The controlled toggle is the test's instrument;
production's honesty is the SEPARATE thing under test (pinned by GATE B `:1390`). Mode-11 sentinel
guards it: `zero==flat` byte-id AND `mom≠flat` consumed → toggle proven non-vacuous.

**5. Deliverable-4 named accessor SKIP is correct.** Only slope index in S9 production is `slot[...,1]`
at the ONE construction site (`sn.py:1681`), immediately below its `AVERAGE_MOMENT`-named slot-0 sibling
+ a descriptive comment. NO scattered `slot[...,1]` read across consumers (`_inflow_to_moments` threads
the WHOLE slot). A `FACE_SLOPE_MOMENT` constant for one local write = gold-plating (anti-#10). Test-side
projectors already use literal `[...,1]` without a constant → adding one only in production would
INTRODUCE asymmetry. Broad typed predicate = #246 (deferred, right).

## NITS (do-now, non-blocking — docstring/comment honesty only)

- **NIT-1 (the only CONCERN, anti-#11):** three structurally-parallel transverse-Legendre loops now
  coexist (production `_project_inflow_to_face_moments`, test `_face_transverse_legendre`, test
  `_face_transverse_buffers`→delegates to the 2nd). Two MUST stay independent (L11). Bound the dup:
  add one docstring sentence to `_project_inflow_to_face_moments` naming the collapse trigger — a 3-D
  face (`n_face_moments = per_axis^{d-1} > 2`) needs a tensor lift on BOTH sides, the rule-of-two→three
  trigger, folds into #263's collocation seam. Converts a latent twin into a tracked one.
- **NIT-2 (optional):** `test_ld_2d_boundary_promise.py:110 _amplify_boundary_slope` hand-lists every
  frozen-dataclass field to rebuild the subclass; a future field is silently dropped. One-line comment
  "field list must track SN2DCartesianLDStressMMSCase".

## PASSES (specifically called out — reinforce)

- Math/domain alignment: `slot[...,AVERAGE_MOMENT]=(wq*psi).sum/W2` reads as `⟨ψ,P₀⟩/⟨P₀,P₀⟩`,
  `slot[...,1]=(wq*xi*psi).sum/(W2*mean_p1sq)` as `⟨ψ,P₁⟩/⟨P₁,P₁⟩`; `mean_p1sq` named with `=1/3`
  comment; `W=quadrature.weights.sum()` computed not hardcoded (anti-#14).
- Promoted gate module reads clean: docstring records VERDICT+MECHANISM+promise+revisit-trigger
  (Sphinx-is-the-brain); Mode-8 conversion clean (`pytest.fail`/`np.testing`, no bare assert under -O,
  verified); import shim dropped (`from .test_mms_ld_2d import` package-relative sibling); each verdict
  pin carries its verdict-flip trigger in the fail message (anti-#17 — loose `_SUBFLOOR_FRACTION=0.30`
  band defended with semantics).
- Docstring honesty updates (HONEST SCOPE/class/external_source + symbolic test) are in-module honesty
  (code's local prose stops lying), NOT Sphinx-narrative scope creep (rich narrative owed to archivist,
  DISPATCH_REQUEST emitted — confirmed honored).

Extends [[issue-257-s9-boundary-moment-honesty-closeout]] (method-implementer). Same posture as the
#247/#251 Leg-A/B reviews.
