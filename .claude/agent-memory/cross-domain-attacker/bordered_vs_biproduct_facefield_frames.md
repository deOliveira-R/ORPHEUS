---
name: bordered-vs-biproduct-facefield-frames
description: VERDICT on the BorderedOperator (Hyp A) + FaceField-codim1-dual (Hyp B) two-primitive proposal for V=V_bulk⊕V_trace⊕V_sd — "append rows" splits into 3 distinct structures; ψ½ is an angular trace, a bad first consumer for a spatial-staggered FaceField.
metadata:
  type: project
---

# Bordered-operator vs biproduct vs FaceField — the "append rows" trichotomy

Branch `refactor/sn-walk-unification`; branch-verified reads (L-005), NOT Nexus.
Two design hypotheses pressure-tested for the augmented composite
`V = V_bulk ⊕ V_trace ⊕ V_sd`.

## The load-bearing code facts (branch-verified)
- RQI bordered system is built INLINE as a dense `(n+1)×(n+1)` with the corner
  **exactly 0** (`eigenvalue.py:545-556`: `border[:n,:n]=F−kA`, `border[:n,n]=−Av`,
  `border[n,:n]=v`, `border[n,n]` never set). Zero corner = KKT/saddle hallmark.
- The composite is ALREADY a biproduct block-operator algebra: `BlockRole
  {BULK,FULL,BOUNDARY}` on the direct sum, the 2×2 block matrix documented
  `operator.py:189-201`; `_join_block_roles` derives the sum's role;
  `_AdjointOperator` realizes `A.H = G⁻¹AᵀG` per-block. The boundary-boundary
  block `A_ss = −B` is a **nonzero realized boundary law** (vacuum/reflective).
- Composite metric `G` is block-diagonal PSD: `G_bulk=V·w>0`, `G_trace=|Ω·n|·w≥0`
  (zero on tangential ordinates, MP-pseudo-inverse), `G_sd≡0` (ghost metric —
  `(1−µ²)|_{µ=±1}=0`). `full_field_space.py` + `starting_direction_space.py`.
- Two FLAT-buffer codim-1 layouts, no structured addressing: `FaceLayout`
  (string face-name keys, `face_layout.py`) and `StartingDirectionSpace`
  (`(level,sign)` keys, `starting_direction_space.py`) — the latter's docstring
  says it copies "the exact FaceLayout discipline, key typed instead of stringly."

## VERDICT — native structure is DIFFERENT and SMALLER than the two proposed primitives

**Hyp A (single BorderedOperator) is an OVER-ABSTRACTION.** "Append/remove rows &
cols" fuses THREE structures with DIFFERENT invariants; a lone type enforces the
wrong one on ≥2 of 3. Discriminator = the CORNER block + metric definiteness:
1. **Saddle-point / KKT / bordered** (Benzi–Golub–Liesen; Peters–Wilkinson RQI;
   RT0 mixed-form div-constraint). Border row = a FUNCTIONAL (constraint vᵀ / ∫=1),
   corner **D=0**, metric INDEFINITE, well-posedness = inf-sup/LBB. The GENUINE
   border. Current instances: RQI only (mixed-form pending). Rank-1 scalar.
2. **Biproduct / dagger-biproduct trace** (the bulk⊕trace⊕sd composite). "Blocks"
   are TRACES ι*/ι_* (restriction/extension), corner `A_ss=−B` ≠0 (a boundary
   law), metric block-diag PSD, well-posed by the triangular sweep. NOT a border.
   The "border ⊣ un-border" adjoint the proposal intuits IS the biproduct
   inject/project dagger pair `π_iι_j=δ_ij, ι_i=π_i†` — ALREADY named in the
   algebra; do not re-mint it.
3. **Galerkin coarsening** (condensation/homogenization = "un-border"). This is
   the multigrid R⊣P prolongation/restriction pair, LOSSY PG `project=G⁻¹M`
   (memory: homogenization is Petrov-Galerkin), NOT the exact Schur complement
   `D−CᵀA⁻¹B`. Conflating condensation with Schur "un-bordering" is imprecise:
   Schur is exact elimination; homogenization is a lossy solution-weighted proj.

**Hyp B (cell/face staggered duality) is the RIGHT frame — but mimetic, not
cohomological, and ψ½ strains it.**
- The primal/dual staggered pairing (cell scalar ↔ face current, div⊣−grad
  adjoint) = mimetic-FD / staggered-grid / Hodge. FIRES for the SPATIAL faces:
  diffusion J±, RT0 mixed-form, CP interface currents, MoC track-faces.
- de Rham COHOMOLOGY (d²=0) does NOT fire (L-001 holds again): the transport
  trace pair is `ι_*∘ι*=id` (dagger adjoint), not a differential. Adopt
  mimetic/Hodge vocabulary; REFUSE chain-complex/homology.
- **The one genuine, correct unification: the face measure = |normal component of
  the streaming flux through the face|, vanishing at directions tangent to the
  face.** `|Ω·n|` (spatial face, vanishes at grazing Ω) and `(1−µ²)` (angular
  µ=±1 face, vanishes at the pole — the µ=±1 rays are characteristics) are TWO
  INSTANCES of this one functional. Worth naming; today it's built twice.
- **ψ½ is an ANGULAR trace, a BAD first consumer for a spatial-staggered
  FaceField.** ψ½ is full-in-space, keyed by µ-level — a restriction ι*_angle on
  the ANGULAR factor of the phase-space product `V_space⊗V_angle` (already
  `SumOfTensorProductsOperator`), NOT a spatial-incidence face. It has no
  `cell(i,j)→faces` map. Forcing it into a spatial FaceField is Smell #16 shape-2
  in REVERSE — one representation onto two genuinely different objects.

## Buildable NOW (defer nothing) vs DEFER (with trigger)
- NOW: (i) name the streaming-flux-normal face measure (unify the two trace metric
  constructions onto one function); (ii) generalize `FaceLayout` to a typed key so
  `StartingDirectionSpace` stops being a second flat-layout impl (Smell #16 shape-2).
- DEFER-with-trigger: (iii) the structured spatial `FaceSpace` (cell→face
  incidence) — first consumer = RT0 mixed-form diffusion currents / CP interface
  currents, NOT ψ½; (iv) a `SaddlePointOperator`/`BorderedSystem` (inf-sup
  invariant, zero corner) — trigger = the mixed-form (2nd genuine saddle instance);
  THEN retire RQI's inline `(n+1)×(n+1)` onto it as first consumer.

## Discriminating first tests (each REDs a naive/wrong impl)
- A-refute: assert RQI border corner `==0` (indefinite/KKT) while composite
  `A_ss` is a nonzero `BoundaryOperator` and `G` is PSD. A "same bordered
  structure" claim passes only "both are block matrices"; fails corner+inertia.
- B one-measure: ONE `face_streaming_normal(face)` must reproduce BOTH
  `G_trace==|Ω·n|·w` AND `G_sd==0` (array_equal, 0 ULP). The current two-source
  construction (AngularTraceSpace builds |Ω·n|; StartingDirectionSpace hard-codes
  zeros) cannot — twin-source reds.
- ψ½-not-a-spatial-face: a FaceField `cell_to_faces(i,j)` is well-defined for the
  spatial trace and ABSENT/raises for ψ½ (negative test).
- condensation≠Schur: exact Schur `D−CᵀA⁻¹B` of the fine DOFs ≠ production PG
  homogenization `G⁻¹M` on a 2-region problem (lossy) — reds a "condense=Schur"
  claim.
- biproduct-redundancy: a new `BorderedOperator.border/unborder` round-trip must
  equal the existing `π/ι` block dispatch bit-for-bit ⟹ redundant (retire) or it
  drifts (twin). array_equal.
