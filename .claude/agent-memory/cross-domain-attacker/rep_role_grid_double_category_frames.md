---
name: rep-role-grid-double-category-frames
description: Frame attack on the (Representation × Role) carrier grid type-machinery question. VERDICT the grid is a DOUBLE CATEGORY (objects=cells, horizontal 1-morphisms=Rep-morphisms M/R Role-fixing, vertical=Role-morphisms C/Λ/F Rep-fixing); scattering S=(1/W)·R∘Λ∘M IS the 2-cell frame.conjugate(Λ), naturality=interchange law, the scattering.py:587 bit-identical windowed-vs-full crosscheck IS its coherence witness. The brief's 4 carrier options resolve HARD: (d) the current MI leaves Leaf(RoleMixin,RepBase) is the UNIQUE permitted form (NOT a baseline to beat) — Role changes the __add__ SIGNATURE (torsor vs vector) ⇒ must be a CLASS (refutes phantom-param a,c); Rep changes the SHAPE ⇒ must be a CLASS (refutes b); MI is the only both-classes form. A Field[Rep,Role] phantom carrier is IMPOSSIBLE (params erased, can't specialize dunders). The (Rep,Role) parametrization the brief wants BELONGS on the OPERATOR contract LinearOperator[Din,Cout] (Din/Cout ARE the carrier's Rep,Role) where it is already HALF-BUILT — finish it so M is Role-generic [V:FluxRep,W:FluxRep], Λ is Rep-generic, retiring the 7 @overload stubs. Branch: read on main, ground truth.
metadata:
  type: project
---

# (Rep × Role) carrier grid — double-category / fibration / torsor frame attack

DESIGN VERDICT (detection feeding a "what needs to change" delta). Read on
`main`, clean tree — file:line is GROUND TRUTH (L-005 N/A). This is the
CONSOLIDATION of five prior memos into ONE structural picture: the grid the
brief crystallized is the union of the #208 affine torsor (Role axis), the #226
container endomorphism algebra, the #257 carrier layering, the
coefficient-field promotion (CoefficientRole), and the projection/reconstruction
Frame (M/R = horizontal axis). The brief asked which of 4 type-machinery options
(a phantom `Field[Rep,Role]` / b Role-outer `Flux[Rep]` / c Rep-outer
`Angular[Role]` / d current flat MI leaves) realizes it. The answer is decided
by TWO obstacles, not taste.

## Ground-truth structure (the grid IS already in normal form)
- Leaves are MI intersections `Leaf(RoleMixin, RepBase)`:
  `AngularFlux(FluxRole, AngularField)` (`angular_flux.py:63`),
  `CrossSectionField(CoefficientRole, ScalarField)`,
  `AngularDisplacement(Displacement, AngularField)` (`angular_displacement.py:32`,
  a 43-line THIN leaf — only `_SPACE_NAME`/`UNITS`),
  `AngularSourceSink(AngularField)` (NO role mixin — Source IS plain Field algebra).
- **Role arithmetic is ALREADY in ONE place per role**: `FluxRole` (`_flux_role.py`
  — torsor: `flux+flux`→TypeError, `flux−flux`→Displacement, `affine_combination`),
  `CoefficientRole` (`_coefficient_role.py` — marker, keeps plain Field algebra),
  `Displacement` (`_displacement.py` — marker + ρ/Aitken diagnostics).
- **Rep shape/methods ALREADY in ONE place per rep**: `AngularField`/`ScalarField`/
  `MomentField`/`BoundaryField` (`_bases.py` — `_phase_space_shape`,
  `integrate_angular`, `L`-match, `TraceSpace`/`face_view`).
- Grid ≈ 4 reps × 4 roles ≈ 16 leaves, each a 2-line MI binding. **NO cell-by-cell
  duplication** — the smell the brief worried about does not exist on the carrier.
- The OPERATOR contract is ALREADY two-param: `LinearOperator(Protocol[Din,Cout])`
  (`operator.py:275`, Din contra / Cout co / PEP-696 `Cout=Din`). **`Din`/`Cout`
  ARE the carrier's `(Rep, Role)`.** `M:[AngularFlux,MomentFlux]`,
  `Λ:[MomentFlux,MomentSourceSink]`.
- Scattering `S=(1/W)·R∘Λ∘M` is LITERAL: `frame.conjugate(Λ) =
  OperatorProduct(R, OperatorProduct(Λ, M))` (`scattering.py:696-700`); `1/W` lives
  OUTSIDE (`:667`); windowed `R∘Λ` arm == full `R∘Λ∘M` arm bit-identical (`:587`).

## The THREE frames (each FAIL-ABLE)

### Frame 1 — DOUBLE CATEGORY (the precise name)
Objects = cells `(Rep,Role)`. **Horizontal 1-morphisms** = Rep-morphisms (M, R,
integrate_angular) — Role-fixing, built from the `FrameBase` analysis/synthesis
pair. **Vertical 1-morphisms** = Role-morphisms (C, Λ, F) — Rep-fixing, the cross
sections. "Natural in the untouched axis" IS the interchange law. Scattering is
the **2-cell**: `frame.conjugate(A)=R∘A∘M` = canonical conjugation of a vertical
1-morphism by the horizontal M/R adjoint pair. The `:587` bit-identical
windowed-vs-full crosscheck IS the interchange COHERENCE WITNESS (a theorem, not
a numerical coincidence to re-pin per path). `1/W` = horizontal frame-operator
normalization (`S=4π·I`), structurally outside the square.
FIRST TEST (DISCRIMINATES): Path-A `R(Λ(M(ψ)))` vs Path-B direct per-ordinate
`S·ψ`, assert `array_equal` 0 ULP (NOT allclose). An impl applying `1/W` INSIDE
the square (between Λ and R) violates the 2-cell and REDs. Second: `(R∘Λ)∘M` vs
`R∘(Λ∘M)` association identical; a Λ that secretly reads N (Rep-dependence) REDs.

### Frame 2 — AFFINE TORSOR refutes the phantom-parameter options (THE decider)
Role changes the arithmetic INTERFACE: Flux = affine space A over difference
vector space V; Source/Displacement/Residual ARE V. A torsor `__add__` is `A×V→A`
(forbids `A×A`); a vector `__add__` is `V×V→V`. **A Python method signature is a
property of the CLASS, not of a `Generic[Role]` parameter** — phantom params are
erased and CANNOT specialize dunders. So `Field[Rep,Flux]` and `Field[Rep,Source]`
under a phantom `Role` share ONE `__add__` body ⇒ cannot make `flux+flux` raise
while `source+source` succeeds. **HARD REFUTATION of options (a) and (c).** A
runtime `role`-field branch is the stringly-typed anti-pattern
(illegal-state-representable: `replace(f, role=Source)` bypasses the gate).
⇒ Role MUST be a real class. It already is (`FluxRole` vs `CoefficientRole`).
FIRST TEST (DISCRIMINATES, NEGATIVE): the phantom impl must satisfy BOTH
`Field[Ang,Flux]()+Field[Ang,Flux]()` RAISES and `Field[Ang,Source]()+...`
succeeds — one shared body cannot; and `replace(f, role=Source)` type-checking +
bypassing the torsor gate REDs the runtime-branch escape.

### Frame 3 — FIBERED CATEGORY: parametrization belongs on the OPERATOR, not the carrier
Carrier fibered `p:E→B`, base B = Representation (shape/space layer), fiber =
Role. A Role-morphism (Λ,C,F) = cartesian morphism WITHIN a fiber (fixed Rep); a
Rep-morphism (M,R) = base change lifting to the total space. M is "Role-generic"
BECAUSE a base change fixes the fiber coordinate and lifts uniformly; Λ is
"Rep-generic" BECAUSE the fiber morphism is the same over every base point. This
IS `LinearOperator[Din,Cout]`: `M:[V:FluxRep,W:FluxRep]` (Role preserved via the
`Cout=Din`-role default), `Λ:[V, Source-sibling-of-V]` (Rep preserved). **The
carrier does NOT need `Field[Rep,Role]` — the operators express it through their
`[Din,Cout]` leaf types, and role-preservation is a THEOREM (base change fixes
the fiber), not a carrier parametrization.** Re-uses
[[operator_protocol_mixin_collapse_frames]] Frame-3 (ruled from the operator
side), now applied to the brief's carrier question. The 4 scattering + 3 fission
`@overload` stubs (`scattering.py:1388`, `fission.py:466`) ARE the un-derived
cartesian-morphism signatures.
FIRST TEST (DISCRIMINATES): pyright infers `M.apply(angular_flux)` →
`HarmonicMomentFlux` (Flux→Flux) and `Λ.apply(moment_flux)` →
`HarmonicMomentSourceSink` (Moment→Moment) with NO `@overload`. A composite-primary
typing `M:[FullField,FullField]` loses the Flux/Source distinction →
`Λ(M(flux))+flux` type-checks (WRONG cross-role sum); fibration-primary REDs it.

## VERDICT (ranked, with the elegance call)
1. **(d) current MI leaves WINS and is the UNIQUE permitted form** (NOT "the
   baseline to beat"). Role-changes-arithmetic ⇒ Role is a class (kills a,c);
   Rep-changes-shape ⇒ Rep is a class (kills b); MI is the only both-classes form
   with role-arithmetic once + rep-data once. The ~16 leaves are already the
   elegant normal form — no duplication to remove.
2. **No novel (e) on the carrier.** `Field[Rep,Role]` phantom is structurally
   IMPOSSIBLE (Frame 2).
3. **The (Rep,Role) parametrization belongs on `LinearOperator[Din,Cout]`** (Frame
   3) — already half-built. THE CARVE: finish operator-side `[Din,Cout]` so M is
   Role-generic, Λ Rep-generic; RETIRE the 7 `@overload` stubs (the payoff metric).
   Gated on the pyright variance spike already flagged in
   [[operator_protocol_mixin_collapse_frames]].

ELEGANCE CALL: building a (Rep×Role)-typed CARRIER is NOT worth it and mostly
impossible — all 4 criteria come out AGAINST (less structure-exposing: hides the
torsor/shape class-distinctions behind erased params; less expressive: forces
runtime branches; MORE complex: re-introduces dispatch the class identity already
does; zero algorithmic gain). WORTH building: (i) finish operator `[Din,Cout]`
(retires 7 @overloads), (ii) name `frame.conjugate` as the 2-cell + reframe the
`:587` crosscheck as interchange coherence (0-ULP `array_equal`), (iii) minor:
derive `_DISPLACEMENT_CLS` grid-edge instead of hand-setting per flux leaf.

## Cross-method (durable)
- Role axis is method-UNIVERSAL (transport-resolvent backbone L-007: SN/MoC/CP
  share C/Λ/F); only the Rep axis is method-specific. Grid's Role COLUMN shared
  (mixins in `orpheus/transport/`), Rep ROW per-method. `ScalarDisplacement`
  carries `contraction_ratio` identically on a CP fixed-point (inherits
  `Displacement`).
- The horizontal-morphism axis IS the existing `FrameBase` (analysis=M,
  reconstruction=R, conjugate=2-cell) — Frame 1 names existing code, builds nothing.

## UNEXPLORED (refuted with structural reasons — durable)
- `Field[Rep,Role]` phantom (option a) — REFUTED Frame 2 (params erased, can't
  specialize the torsor `__add__`). NOT merely unexplored.
- Rep-outer `Angular[Role]` (option c) — same Frame-2 refutation (inner Role is a
  param); also inverts the fibration then collapses back to (d).
- Role-outer `Flux[Rep]` (option b) — the ONLY option the torsor obstacle spares
  (Role is the outer class), but FAILS the SHAPE obstacle (Rep param is
  shape-bearing, can't be phantom). Trades torsor-obstacle for shape-obstacle ⇒
  confirms (d). No reformulation ⇒ listed here.
- Category-theory-as-abstract-nonsense — the CONTENT (double cat, fibration,
  torsor, frame-conjugation-2-cell) IS the concrete frame + produces the
  interchange tests; no operad/PROP/2-functor lever adds a discriminator (L-001).
- Tensor networks/MPO — fixed 4×4 product, bond-dim-1 degenerate (L-001).
- Group/rep theory — SO(3) on M/R frame + Λ kernel, not the Role axis or typing.
- Homology — `C¹=C¹_int⊕C¹_∂` present but `∂²≠0` (B∘B≠0), grid has no differential.
- Diff-geom/connection — no curvature in the flat algebraic grid (L-001/L-008).

Cross-refs: [[operator_protocol_mixin_collapse_frames]] (the OPERATOR-side dual,
Frame-3 already ruled leaf-primary/composite-derived — THE carve target),
[[issue_208_delta_psi_affine_frames]] (the torsor that refutes phantom Role),
[[issue_226_container_algebra_design]] (the endomorphism algebra = the diagonal
of the grid), [[projection_reconstruction_frame_pair]] (M/R = the horizontal axis
= FrameBase), [[coefficient_field_promotion_frames]] (CoefficientRole = the
complement role mixin), [[spatial_order_type_vs_property_criterion]] (the
class-vs-param decision criterion this attack instantiates twice).
