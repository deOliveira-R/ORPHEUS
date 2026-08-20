# CS1 design draft — the Energy axis, axis-composed spaces, and the space-family doctrine

**STATUS: DESIGN DISCUSSION CONCLUDED (2026-08-20, round 6) — every fork ruled or
resolved; execution protocol chartered by the user (see THE PROTOCOL FROM HERE,
below). Protocol step 2 (clear-context re-evaluation + §P grounding) EXECUTED
2026-08-20 at `273d431a`: the conceptual layers (§A/§B/§B2/§C, Appendix A) survived
the adversarial re-read unchanged; all corrections are §F execution-layer, marked
⭐ GROUNDED / ⭐ CORRECTED in place; §P results block below the checklist.** This file is the design record and the post-compaction resume surface. Governed by `plan-authoring.md`. The campaign plan of record is
`.claude/plans/space_and_kernel_binding_campaign.md` (§2 points here); when the design
is RULED, §2 absorbs the outcome and this file becomes the design record.

Provenance: proposal v1 + round 2 (quotient-family recap, trivial-axis questions,
composition story, F1–F5 critique, 6-phase skeleton) + round 3 (the MaterialMesh
challenge, the non-compact-invariant question, energy compactness/Bateman, the
partition/order conjecture; rulings on Q2/Q3/Q5) + round 4 (energy is NOT
finite-measure — the integrability correction; the MaterialMesh genuine-reason
test → the Medium finding; leaf order as a boundary convention) + round 5 (the
leak principle articulated; the geometry → {medium, mesh} → pullback hierarchy —
which DISSOLVES the rename; the order-was-conventional verdict; Q8 ruled;
Appendix A written as the archivist's source for the collapse doctrine) + round 6
(symmetry monotonicity down the lattice; #394 filed for the layout/JAX finding;
Q11 ruled by user delegation → CS1.5; the execution protocol chartered).

Markers: **RULED** = user ruling of record · **PROPOSED** = design stance awaiting
verdict · **DERIVED** = derivation produced in this discussion — check the math ·
**OPEN** = genuine fork · `[M]` = measured against the tree (`8d6d5ba3` unless dated).

---

## ▶ THE PROTOCOL FROM HERE (user-chartered, 2026-08-20 — the resume surface)

The user's verbatim charter: *"we will compact context first, and I will give you
the chance to re-evaluate the plan with a clear context and do any further
investigation to ground the plan so that execution is smoother. After that, we can
dispatch the test architect and compact again before actual execution."*

1. **COMPACT** (imminent — this round closes the design discussion).
2. **CLEAR-CONTEXT RE-EVALUATION + GROUNDING PASS** (next session, before any
   code): re-read THIS FILE top to bottom with fresh eyes — the adversarial-first
   method applies to the plan itself (Phase 1 harsh: how would I break this design;
   Phase 2 verdicts). Then run the grounding checklist (§P below) so execution is
   enumeration-complete. Fold the hardened outcome into the campaign plan §2 (it
   currently POINTS here; absorption was deliberately deferred to this pass).
3. **DISPATCH test-architect** (proactive trigger — the carve crosses
   numerics/transport/homogeneous; brief sketch in §T below). Review its plan with
   full context.
4. **COMPACT again**, then **EXECUTE** CS1 steps 1–5 (§F), surgical posture (main
   agent writes, user steers).

**Ruled campaign order (Q11 ruling folded in): CS3 ✓ → CS1 → CS1.5 (the Medium
un-weld) → CS2 → CS4.**

### §P — the grounding checklist for step 2 (each item cheap, none optional)

1. Re-verify the 4 strict-xfail ids + location at HEAD
   (`tests/sn/architecture/test_monomorphic_leaves.py:670` at `00085baf` — line
   drifts).
2. The §8 enumeration: ALL construction sites of the four operator classes
   (`MultiplicationOperator`, `IsotropicScattering`, `IsotropicN2N`,
   `FissionOperator`) — nexus impact + grep; classify each as
   threads-space-already / stays-None / CS1-edits.
3. `.H` and `as_matrix` consumers of homogeneous-path operators in `tests/`
   (the metric-flip audience).
4. Re-grep the F rename batch's call sites (`full_field_space=`) — the §6b set
   must be complete on the day of the edit, not the day of the plan.
5. The degenerate-path readers of `coord`/`axis_widths`/`areas` (CS1.5 design
   input; also confirms CS1 can leave `from_materials` internals untouched).
6. Read `_resolve_basis_shape` (`operator.py:558`) + `from_mesh`'s default chain
   (`multiplication_operator.py:299-338`) to finalize the exact lookup-extension
   edit (the `bulk_space` fallback).
7. `EnergyGrid`'s surface for `from_grid` (edges access; `eq=False`
   identity-equality — the axis's content-eq reads `edges` bytes, not grid
   object identity).
8. Confirm the diffusion threading sites won't collide with the F rename batch
   (`diffusion/solver.py:242-248`).
9. Reproduce the P1 ledger baseline (`[M]` 105 passed / 21 xfailed) + time the
   homogeneous suite (the byte-stability gate's budget).
10. Adversarial re-read of §F's step decomposition against §6b/§6c (call-site
    completeness per step; every gate lands with its witness).

#### §P results (run 2026-08-20 at `273d431a` — all 10 items closed)

1. ✅ The 4 ids confirmed at `:670` (no drift): `groups ∈ {2g,4g} × leaf ∈ {C,F}`
   under `_R1_XFAIL` (strict). The R1/R2/R6 xfail families are distinct; only R1's
   4 ids are CS1's.
2. ✅ Production construction sites of the four classes, classified —
   **CS1-edits**: `homogeneous/solver.py:143` (C via `from_mesh`, gains space
   through the extended chain — call-site text unchanged), `:146` (IsoS/N2N →
   `space=V`), `:193` (F → `space=V`). **threads-already**:
   `diffusion/solver.py:236/:242-243` (space=), `:247-248` (F, keyword renames);
   `sn/solver.py:1333-1341`, `:2804-2805`; `sn/coupled_system.py:419`,
   `:499-503`. **stays-None (deliberate, documented)**: `scattering.py:709`
   `isotropic_kernel` — "Space-anonymous (`space=None`: no composition-guard
   space at this producer-side use)"; legal until CS4. Test-side bare
   constructions (diffusion `test_operators.py:625/650/659/669/853-855`;
   `test_matrix_inverse_operator.py:198-208/:265-267`) stay on the legal
   None path — ⛔ CORRECTED at battery review (F3): the two
   `test_matrix_inverse_operator` blocks build their C via `from_mesh` on the
   DEGENERATE carrier, so after 3b they GAIN `bulk_space` (green either way —
   explicit `basis_shape` wins, and the sum's agreed space resolves through the
   None-silent iso members — but they are NO LONGER None-path witnesses; only
   the positionally-built `IsotropicScattering(mat_xs)`/`IsotropicN2N(mat_xs)`
   there stay bare). ⚠ Nexus `callers` returns 0 for all four classmethods
   (classmethod-attribute dispatch mints no edge — the documented phantom
   limitation); the grep enumeration is authoritative here (uniform spellings).
3. ✅ `.H`/`as_matrix` audience: NO existing test asserts homogeneous-op `.H`
   values — the vv#19 pair is NEW coverage, nothing demoted. The one None-path
   `.H` pin (`test_multiplication_operator.py:326-338`) builds its operator BARE
   and stays on the None path — unaffected. `as_matrix` consumers with explicit
   `(ng,1)`: ⛔ CORRECTED at battery review (F11) — `test_homogeneous.py:140`
   can drop `basis_shape` (its loss rides the chain), but `:284` and `:359`
   build F BARE, so under the ruled no-default-derivation a bare
   `as_matrix()` RAISES; the honest migration threads `space=V` into F at both
   sites FIRST (matching the new production spelling — `:359` is the
   line-for-line mirror of `solver.py:193-194` and would otherwise pin a
   RETIRED idiom, the same warrant that deletes the F xfail rows), then drops
   `basis_shape` at all three (gate D12, mutation M22).
   `test_matrix_inverse_operator.py:199/:211` (bare-op
   subject — keep explicit). The metric flip is ×1.0/÷1.0 per element
   (IEEE-exact) or a counting skip ⟹ bit-identity by construction, gate still
   owed per vv#19.
4. ✅ The complete `full_field_space` site census is in §F (the grounded rename
   batch); the round-2 list conflated F's callers with ScatteringOperator's own
   field sites, and missed `sn/solver.py:2805` + `sn/coupled_system.py:503`; the
   `test_si_gate_dispatch.py:63` hit is a mesh-attribute fake, excluded.
5. ✅ Degenerate-path geometry readers (CS1.5 input): `areas` → `MaterialMesh`
   itself RAISES on 2-D (`material_mesh.py:455` — a live leak-principle datum),
   consumed by `diffusion/augmented_mesh.py:231/:327`; `coord` set at
   `material_mesh.py:215`, read via `mesh.reduced.coord`
   (`radial_characteristic_field.py:308`); `axis_widths` internal to volumes.
   CS1 leaves `from_materials` untouched — confirmed viable (the fake axis keeps
   serving these readers until CS1.5).
6. ✅ `from_mesh` chain read (`multiplication_operator.py:336`):
   `if space is None: space = getattr(mesh, "full_field_space", None)`. The CS1
   edit appends the `bulk_space` fallback AFTER it — SN/diffusion callers
   short-circuit on `full_field_space` and are untouched; the degenerate carrier
   (no `full_field_space` attr) resolves `bulk_space`. `_resolve_basis_shape`
   (`operator.py:559-586`): explicit wins, else `domain.shape`, else the typed
   refusal naming both remedies — after 3b, `MatrixInverseOperator(loss)` and
   `K.as_matrix()` derive `(ng, 1)` from the threaded domain.
7. ✅ `EnergyGrid` (`data/energy_grid.py:106`): frozen, `eq=False`
   (identity-equality), sole field `edges` — strictly DESCENDING (fast-first),
   coerced contiguous float at construction ⟹ the axis's content-eq reads
   `edges` bytes (stable), never grid identity. Carrier→grid path exists:
   `Mixture.eg: ndarray | None` (`mixture.py:65`) + `Mixture.energy_grid`
   (raises when `eg is None`); `MaterialMesh.materials` is the stored dict
   (`material_mesh.py:167`), ng-validated. The axis docstring should carry the
   descending fast-first convention.
8. ✅ Diffusion threading sites (`diffusion/solver.py:236-248`) do not collide
   with the batch: IsoS/N2N already pass `space=` (untouched); the F site is a
   keyword-only rename member. No overlap with 3b's wiring.
9. ✅ Ledger baseline REPRODUCED at `273d431a`: `tests/sn/architecture` →
   **105 passed / 21 xfailed** in 1.9 s. Homogeneous suite: **23 passed** in
   2.42 s (wall 3.27 s) — the byte-stability gate's budget is seconds, not
   minutes.
10. ✅ Findings folded in place: the §F/Q11 micro-cleanup contradiction (⛔
    marked in §F), the done-when F-arm correction (C rows XPASS, F rows cannot),
    the step-3 split into 3a/3b, the uniform `bulk_space` formula, the rename
    batch scope fork (+S recommended), the §T name-injectivity addition.
    Witness check (§6c): both refusal witnesses and the modal-cone witness are
    constructible from shipped inputs at their landing commit (two `get_mixture`
    grids; a test-local generic MODAL `Axis` — within fences).

### §T — test-architect brief sketch (dispatch at step 3)

Scope: the CS1 gate battery — (a) axis intrinsic laws (structural eq/hash;
EnergyAxis grid-content eq incl. the synthetic ng-only arm; frozen-ness);
(b) `of_axes` laws (shape concat; deterministic names — ⭐ added 2026-08-20:
name-derivation INJECTIVE on structural content, because `(name, shape)` IS the
identity until S3, so a name collision between different axis tuples would defeat
Q2 pre-S3; the sharp witness pair is `synthetic(ng=2)`-built vs
`from_grid(2-group edges)`-built, which MUST compare unequal; axes threading
through `__mul__`; per-axis metric ≡ legacy broadcast on toy weighted axes; the
no-densification memory assertion; the uniform `bulk_space` formula on a MESHED
carrier — weights = cell volumes, distinguishing quotient point weight 1.0 from a
genuine one-cell `V ≠ 1`); (c) the refusal witnesses (2g-vs-4g sum → the
pinned "equal domains" provenance fragment; `M⁻¹(2g) @ F(4g)` → product guard) —
§6c: gates land WITH their witnesses; (d) the `.H` loaded/blind PAIR (vv#19: the
neutrality gate cites the counting theorem AND a weighted-toy control shows `.H`
MOVE); (e) byte-stability protocol for the homogeneous suite; (f) the modal-axis
cone-refusal witness; (g) the xfail deletion + positive-floor migration (marker
migration per coding-standards). Brief carries: this file §§A/B/F, the plan §2,
the vv#19/§6c obligations, and — per plan-authoring §4 — every relayed number WITH
its configuration. State what to do if Nexus is absent (sub-agents cannot
ToolSearch).

### §T-R — battery delivered + reviewed; the seven open rulings (main agent, 2026-08-20)

**The battery of record is `scratch/cs1_verification_plan.md`** (test-architect,
2026-08-20; 40 gates A1–A12/B1–B11/C1–C2/D1–D12/E1–E3 mapped to steps 1/2/3a/3b/4;
23-mutation battery M1–M23 with 3 positive controls and a **MUST-STAY-GREEN
anti-claim column** — M17/M19/M23 encode measured non-catchers, so the close-out is
written as "N gates, of which D5/D4a are provable non-catchers for the measure",
never a red count; 15 rejected alternatives; marker migration `[M]` **none owed** —
the deleted rows carry `foundation` only, no `verifies`/`catches`, no `_harness`
registry entry). Findings F1–F12 reviewed with full session context: F2/F3/F11
corrected this file in place (blast radius, §P items 2/3, done-when); F1/F5/F6/F10
folded above; F12 confirms R2 stays untouched by mechanism. Baselines re-measured
by the agent at `4e11731b` match §P item 9 exactly.

The battery's seven open questions, ruled:

1. **Q-T1 (name digest vs float edge cases)** — digest over
   `(weights + 0.0).tobytes()` (kills `-0.0`); non-finite weights REFUSED at
   construction (kills `nan`; a measure's weights must be finite — but NO
   non-negativity guard: CS2's quadrature axes legally carry negative weights).
2. **Q-T2 (unused materials)** — the energy arm reads only `mat_map`-reachable
   materials (leak principle: the mint consults exactly its defining data; the
   broader ng-validation is a GUARD, not an identity input). Ruled into §F above.
3. **Q-T3 (None vs ones)** — CANONICALIZE at construction: all-ones collapses to
   `None`, one spelling per measure, one identity (the None-vs-ones twin would be
   F10's defect shape — same measure, unequal identity). Gate A2's assertion
   sharpens to pin the canonicalization itself
   (`Axis(weights=ones).weights is None`). Ruled into §F above.
4. **Q-T4 (return type of `of_axes`)** — already RULED by §C ("an axis product is
   not a different kind of space"): plain `FunctionSpace` + `axes`;
   `find_factor` → query-by-axis is the recorded CS2 seam. Not reopened.
5. **Q-T5 (byte-gate population)** — producing-only as proposed (the 4
   `derivations` cases + the `eg`-bearing variant + `get_mixture` A/C); the solve
   entry is the k∞ eigenpair, which non-producing mixtures never reach.
6. **Q-T6 (pin the ruled absence)** — YES, gate **D11** ships:
   `from_solver_data(mat_xs=…).domain is None`, docstring naming the wrong-family
   hazard AND the deletion trigger ("WHEN CS4's Optional→mandatory flip lands:
   delete this gate") — the file's strict-xfail deletion-trigger idiom. Without
   it, a silently-added default derivation reds NOTHING (the rows that would have
   caught it are the ones CS1 deletes).
7. **Q-T7 (D5 shelf life)** — mark-for-retirement after the merge cycle, in D5's
   own docstring (a permanent migration snapshot is the artifact class
   `coding-standards` retires); the permanent claims live in `test_kinf_exact` +
   `test_as_matrix_equals_retired_as_dense_loop`.

Also adopted: **D12** (the three `test_homogeneous.py` mirror sites re-point to
the NEW production spelling — F11; edited inside that file, whose `l1 + verifies`
marks are correct for them); **D8** (the two "space-anonymous" refusal messages
re-worded in 3b — F5, `[M]` no test pins either string); the `_R1_XFAIL` reason
re-scope in 3b (F6 — the constant survives on the `:728` annotation row, whose
own non-XPASS after CS1 the agent verified); F4's twin resolution (legacy
`__mul__` keeps densifying + threads `axes`; the non-densifying metric is
selected by `axes is not None`; `of_axes`' docstring names the twin and its CS2
retirement); the CS1.5 re-point comments on B9/B10/B11/D6/D7 at write time.

## A. The formal skeleton (PROPOSED, v1 — stable since round 1)

**Axis** = (index shape, factor measure, basis kind, generator identity).

- `shape: tuple[int, ...]`, rank ≥ 1 (the harmonic axis is rank-2: `(L+1, 2L+1)`;
  a rank-d spatial axis is admissible — see the §F geometric-axis note).
- `weights: NDArray | None` over exactly `shape` — the factor measure. **For axes,
  `None` IS the counting measure, deliberately and always** — axes have no "unbound"
  state, so the legacy two-state `None` ambiguity cannot arise on the new type.
- `kind: NODAL | MODAL` — nodal ⟹ coordinate cone; modal ⟹ no coordinate cone.
- Identity: **structural per subclass**, from day one.

**Space** = ordered product of axes: `FunctionSpace` gains
`axes: tuple[Axis, ...] | None` (compare=False) and `of_axes(*axes)` — shape by
concatenation, NAME derived deterministically from the axes' structural content.
Pairing evaluated per-axis (F4 unrepresentable). `has_coordinate_cone` = all axes
nodal (None when legacy/unknown).

**The counting-measure theorem (Energy).** Multigroup flux = group INTEGRALS
(covariant, extensive); cross sections = flux-weighted group AVERAGES (contravariant,
intensive); the convention makes ∫σφ dE equal Σ_g σ_g φ_g EXACTLY, so the Energy
measure is counting **as a theorem**: metric = I, V ≅ V*, adjoint = counting
transpose. The V/V*-collapse hook (tightness (iv)): condensation = SUM on V /
weighted AVERAGE on V*; collapse adjoint-consistency = that pair being mutually
adjoint — declared on the axis docstring now, built at S7/Campaign 2.

Consequence `[M]`-backed: binding the homogeneous operators is **bit-identical by
construction** (`_AdjointOperator.apply` = G⁺∘Aᵀ∘G with identity metrics,
`operator.py:1307`; guards compare spaces, never values; `as_matrix`/
`MatrixInverseOperator` already derive `basis_shape` from `domain.shape`).

## B. The collapse doctrine — what survives into a degenerate axis, and why (DERIVED round 3, CORRECTING round 2)

> Correction trace. Round 2 proposed COMPACTNESS as the discriminator; round 3's
> energy question refuted its persist-vs-drop half and produced the consultation
> principle. Round 4 (user) refuted the remaining "energy is finite-measure" claim:
> the energy domain (0, ∞) has INFINITE Lebesgue measure — what converges is the
> INTEGRAL OF THE ADMISSIBLE FIELDS (the tail of the spectrum). The user's
> formulation — "does the integration of that infinite tail converge?" — IS the
> formally correct question (L¹ membership of the admissible class), and the
> discriminator is a property of the FIELD CLASS the physics admits, never of the
> bare domain.

**An axis survives its own collapse iff the collapse leaves data the surviving family
still CONSULTS — in the pairing, or in identity/guards. Whether the collapse
integrates or normalizes is decided by the INTEGRABILITY OF THE ADMISSIBLE FIELD
CLASS over the collapsed domain.** Three clauses:

1. **Collapse along NON-COMPACT GROUP ORBITS ⟹ NORMALIZE (quotient).** The
   symmetry that defines the problem class forces the fields CONSTANT along orbits
   of infinite Haar measure — a constant over infinite measure is never L¹, so
   integration is impossible *structurally*, whatever the physics. The only
   surviving functional is the normalized average ("per unit orbit measure"), a
   genuine convention consulted by the member's OWN pairing ⟹ the axis persists as
   the quotient point with unit weight. `[M]` live: `volume_measure` = 1.0 consumed
   by the reaction-rate functionals inside the homogeneous solve. — Homogeneous
   spatial slot (translations of ℝ^d, non-compact). In the modal family the
   surviving datum is the mode parameter (buckling's B).
2. **Collapse by PARTITION-INTEGRATION of an L¹ field class (no symmetry acting)
   ⟹ per-cell integrals (nodal mesh-axis).** Energy: no invariance group forces
   constancy; the physics makes spectra integrable (source-bounded above — χ decays
   super-exponentially; Maxwellian → 0 as E → 0), so per-interval integration
   converges INCLUDING the top interval's infinite tail. The practical grid top
   E₀ is the LIBRARY's assertion that the neglected tail is below tolerance — a
   data-layer truncation statement, not a space-layer fact. The PARTITION
   (boundaries, weighting spectra) is problem data consulted by identity/guards and
   by the V/V* pairing ⟹ the axis persists down to its terminal ONE-CELL member.
   — **EnergyGrid is a 1-D mesh in energy — the group boundaries are its FACES
   (user, RULED into Q2)**; the ONE-GROUP member is the one-cell energy mesh whose
   edges + weighting spectrum define σ̄ — exactly what the Bateman/depletion pairing
   ⟨σ̄, φ⟩ consumes. Likewise a one-CELL slab keeps its axis with weight V ≠ 1.
3. **Whole-domain integration over a COMPACT canonical orbit ⟹ the axis DROPS.**
   S² is a compact-group orbit: finite Haar measure, integration always available,
   and the total (4π) is UNIVERSAL — nothing problem-specific survives on the
   axis; the rebroadcast convention lives on the EMBEDDING ι (§I.9), consulted
   only when LEAVING the family ⟹ scalar spaces carry NO angular slot and are a
   different family.

The one-line discriminators: *can the admissible fields be integrated over the
collapsed domain?* (symmetry-forced constancy on infinite orbits ⟹ no ⟹ normalize)
— and *is the surviving convention consulted inside the family (⟹ axis) or only at
re-embedding (⟹ arrow)?* The tree's `(ng, 1)`-with-no-angular-slot, the one-group
`(1, *spatial)` convention, and the per-face partial-current slots (hemisphere
collapse parameterized by n̂ — the parameter survives on the FACE structure) are all
retrodicted.

> ⭐ **Preservation note (user, round 5): "make sure we don't lose it so the
> archivist can write it later — it was hard to steer until we got it out."**
> Appendix A of this file is the archivist's source: the full dialectic record
> (both refuted versions WITH their refuting questions, the mechanism, the
> retrodictions). The `spaces.rst` seed (§F done-when) writes FROM Appendix A.

**Answers to the round-3 questions this section owes:**

- *Is Energy non-compact?* The DOMAIN is infinite-measure; the ADMISSIBLE CLASS
  is L¹ (convergent tails — the user's formulation, confirmed) and no symmetry
  forces constancy ⟹ clause 2, integrate-side: nodal mesh-axis persisting to one
  cell. The grid top edge is a library-layer truncation assertion; the edges are
  axis identity (Q2).
- *What property does a degenerate NON-COMPACT axis store, and why useful?* The
  **density convention** (the quotient's unit normalization weight) as a
  first-class, composable object. Useful three ways: (a) the pairing consumes it —
  rates/normalizations follow automatically if the convention ever changes (per
  lattice cell of volume V ⟹ weight V, every functional follows); (b) family
  coherence — the B=0 fiber and B≠0 members share the slot, so fiberwise machinery
  (Fourier ρ(B)) reads a uniform signature; (c) boundary-of-family maps
  (homogenized constants into meshed contexts) get a home and a Jacobian (the
  measure ratio) instead of an invisible unit-convention shift — the classic
  missing-volume-factor bug class in homogenization.

**Standing consequences (from round 2, surviving the correction):** explicit point
axis, shape `(ng, 1)` (Cardinal-2 single-generic-body; B=0 fiber; the pairing
argument now sharpened by clause 1); point minted as a QUOTIENT object (not forged
from `AxisMesh(edges=[0,1])`); kind NODAL (positive 1-dim basis ⟹ nodal/modal
coincide at the B=0 fiber); buckling = size-1 MODAL spatial axis, live angle,
𝔽 = ℂ ⟹ S5/Campaign 2, CS1 must only not foreclose; the P_N hierarchy is the
angular buckling (irreps of rotations ↔ irreps of translations) and the trivial
angular slot stays ABSENT on scalar-family spaces (clause 3).

## B2. The carrier question — the genuine-reason test fails for "mesh" and passes for "medium" (round 4, superseding round 3's keep-with-micro-cleanup)

Round-3 challenge (user): if MaterialMesh survives on the homogeneous path, name the
functionality even the degenerate case genuinely needs; MaterialMesh was created so a
case could instantiate as SN or Diffusion (strong spatial concepts) and to serve
"run detailed case → condense/homogenize → instantiate diffusion". Find the genuine
concept or decompose/rename. Round 4 runs the test against the tree:

**`[M]` The inventory verdict.** MaterialMesh's public surface is `ng`, `volumes`,
`integrate_per_group`, `volume_measure`, `spatial_shape`, `ndim`,
`material_xs_field()`, `from_materials` — ALL region-carrier members — plus exactly
ONE squarely geometric property (`areas`, radial face areas, consumed by method
layers). The homogeneous path consumes `{material_xs_field, volume_measure, ng,
spatial_shape}` — none mesh-intrinsic. And `[M]` `SNMesh(MaterialMesh)` /
DiffusionMesh inherit it explicitly as the "method-agnostic DATA block"
(`augmented_mesh.py:96/:242-258`; "an SN phase space *is a* MaterialMesh").

⟹ **The genuine-reason test FAILS for the mesh concept and PASSES for the concept
the class already realizes: method-agnostic material content over a measured
spatial structure — a MEDIUM.** The mesh-ness lives in the NAME and in one
construction path's inputs (the AxisMesh tuple), not in the realized surface. The
user's workflow confirms it: the common thing solvers instantiate from is the
medium; homogenize/condense are medium → medium morphisms (meshed → coarser →
quotient); method instantiation (SN/diffusion layers) is what ADDS spatial-method
structure. The homogeneous solver is the method-ization of the TERMINAL quotient
medium (exact dense solve on Energy ⊗ point).

**Proposed shape (PROPOSED — the Medium finding):**

1. **Rename the base class to the medium concept** (name candidates, user rules:
   `Medium` — domain-native, "infinite homogeneous medium" is literally the
   degenerate member, `Medium.infinite_homogeneous(mix)` reads as the physics;
   alternatives `MaterialRegions`, `MaterialDistribution`). This also realizes the
   user's original intent at the API: the homogeneous constructor receives ONE
   Mixture — while the region/view machinery stays the single generic body.
2. **Type the spatial structure inside it as a small union** —
   `spatial_structure: MeshBackedRegions | QuotientPoint` — instead of optional
   geometry attributes: mesh-backed extras (`axis_widths`, `coord`, `areas`) live
   ONLY on the mesh-backed member (illegal states unrepresentable — reading
   geometry off the quotient member becomes a type error, not a runtime surprise);
   `volumes`/`spatial_shape` forward from either. The fake
   `AxisMesh(edges=[0,1])` mint dies BY CONSTRUCTION (the quotient constructor
   builds `QuotientPoint`), subsuming round 3's micro-cleanup. The SpatialAxis
   generator story aligns exactly: `SpatialAxis` forgets from `spatial_structure`
   (mesh-backed ⟹ `from_mesh`; quotient ⟹ `quotient_point`).
3. **Keep the inheritance as-is** (`SNMesh(Medium)`): the is-a was chosen
   deliberately for bit-identical data-block reuse; re-litigating is-a → has-a is
   real churn with no charter. The method layers read geometry through the typed
   `spatial_structure` member.

**Scope honesty:** this is its OWN carve (base class of the whole mesh hierarchy +
tree-wide rename of a ubiquitous name + method-layer attribute reads) — NOT CS1.
Charter it as a pre-CS2 phase ("CS1.5") or CS2's opener; the 3-search rename audit
applies. **CS1 is de-risked either way:** `bulk_space` hangs on the base class and
RIDES the rename unchanged (the property belongs to the concept, not the name).
Round 3's micro-cleanup is DROPPED if the Medium carve is chartered (avoid touching
`_init_data` twice) and retained as the fallback if the carve is deferred. Q10's
Protocol variant (structural typing instead of the union) folds into the same
charter decision.

**⭐ ROUND-5 CORRECTION — the user's "why does Medium feel PARALLEL to
geometry/mesh instead of hierarchical?" exposed that round 4's shape was still a
weld, and the true hierarchy DISSOLVES the rename:**

The round-4 proposal renamed MaterialMesh → Medium, i.e. still one class carrying
materials AND computational cells. But the physical medium — *what fills space* —
exists independent of any mesh: material interfaces are PHYSICS (zone boundaries);
cells are NUMERICS (a discretization choice). The honest lattice:

```
Geometry (domain + symmetry [+ coords])
   ├── Medium  = materials assigned over the geometry's PHYSICAL regions
   │            (mesh-free; the infinite homogeneous medium = full translation
   │             symmetry + one mixture; 1-D realization = interface positions +
   │             material per segment)
   └── Mesh    = a discretization of the geometry into computational cells
                 (conforming to the medium's interfaces — an ASSERTABLE guard,
                  today a silent user obligation: Mesh1D takes edges + mat_ids
                  coupled by hand)
Medium × Mesh ──named constructor──▶ MaterialMesh  (the PULLBACK: material-per-
                 CELL views, mat_map, per-cell XS — a DERIVED view, not a
                 primitive; its NAME is honest for exactly this)
MaterialMesh × method data ──▶ SNMesh / DiffusionMesh (method layers, as today)
```

Consequences: (1) **No tree-wide rename** — `MaterialMesh` is the *right* name for
the pullback ("materials on a mesh"); the round-4 rename is WITHDRAWN. The
dishonesty was never the name — it was making the pullback the HOME of the
homogeneous case, which has no mesh to pull back through. (2) **`Medium` is minted
ABOVE**, new; the homogeneous solver re-homes to it
(`Medium.infinite_homogeneous(mix)`); the degenerate carrier is revealed as the
pullback of the infinite medium through the QUOTIENT spatial structure — today's
`from_materials` wearing the wrong constructor story. (3) The round-4 typed union
(`spatial_structure: MeshBackedRegions | QuotientPoint`) survives INSIDE the
pullback class, so face/area surface is type-absent on the quotient member — the
user's leak test ("face area must be *not there*, not an error") passes at the
type level. (4) The mesh-conformity guard (every medium interface ∈ mesh.edges)
becomes constructible at the named binding — a latent misalignment bug class
closed by construction. (5) The axis generators now sit at distinct pipeline
stages: EnergyAxis ← medium (ng/grid); SpatialAxis ← mesh OR the geometry's
quotient; QuadratureAxis ← method — which is exactly WHY welding carriers welded
spaces. **Scope**: new class + one named constructor + homogeneous re-home +
conformity guard; substantially SMALLER than round 4's rename. Staging: CS1
untouched (`bulk_space` interim on the carrier; migrates to Medium when Q11
lands).

**Round 6 (user) — the SYMMETRY MONOTONICITY LAW down the lattice.** Symmetry
lives at Geometry, and every arrow DOWN the lattice at best preserves and usually
REDUCES it — never increases: G_medium = G_geom ∩ Stab(assignment) (a material
assignment breaks symmetry unless invariant); G_mesh = G_geom ∩ Stab(cells) (a
uniform mesh keeps discrete translations, an unstructured one nothing);
G_pullback ⊆ G_medium ∩ G_mesh. Consequence for the collapse doctrine: **clause
1's quotient constructor consumes the MEDIUM's surviving symmetry, not raw
geometry's** — the infinite homogeneous medium is exactly the member whose
assignment stabilizer is EVERYTHING (uniform), so the full translation group
survives to be quotiented; buckling restricts to an irrep OF the surviving group.
This is also why the medium (not the mesh) is the homogeneous solver's honest
home: the quotient is licensed one level ABOVE the mesh, by symmetry the mesh
would only ever shrink.

## C. Composability and the composition end-state

The four family members, against the live tree (`[M]` all four layouts):

| solver | axes | shape | live layout |
|---|---|---|---|
| homogeneous | (Energy, point) | `(ng, 1)` | the `(ng, *spatial)` contract, `spatial_shape == (1,)` |
| diffusion | (Energy, Spatial) | `(ng, *spatial)` | `scalar_bulk` mint, `diffusion/augmented_mesh.py:361` |
| SN | (Quadrature, Energy, Spatial) | `(N, ng, *spatial)` | `sn_bulk` mint (memo A F2b) |
| PN / harmonic | (Harmonic(L), Energy, Spatial) | `(L+1, 2L+1, ng, *spatial)` | SH production layout, `space.py:234` |

- Same ROLE, two REALIZATIONS on angle: `QuadratureAxis` (nodal, w_n, coordinate
  cone) vs `HarmonicAxis` (modal, per-ℓ mass metric, none) — connected by frames
  (`Moments`/`Ordinates`, B/A tightness), never by equality.
- **One composition mechanism (end-state):** a space IS its axis tuple; `of_axes`
  on the `FunctionSpace` BASE (type-vs-property: an axis product is not a different
  kind of space). CS2 collapses the three mechanisms at the live mints
  (`augmented_mesh.py:1099`, `_bases.py:381/:503`; the space.py factories are DEAD
  `[M]`); `TensorProductSpace` retires onto axis concatenation (`find_factor` →
  query-by-axis); `CoupledSpace` stays the UNIQUE ⊕ that
  `partition_along(axis, blocks)` targets; blocks by retract conjugation
  `A_ij = π_i ∘ A ∘ ι_j` (S6/2.10). `SNMethodSpace` = generator record,
  rename-not-redesign (#19/0.9).

**Axis order (Q1) and the partition/order question (round 3).** The live mints share
ONE canonical order — **(angular-method, energy, spatial, [trailing
spatial-moment])** — so "storage order = axis order" and "one canonical order" agree
on every real space (Q1 RESOLVED-proposed, round 2). Round 3's conjecture (the
pre-partition order becomes essentially meaningless; partition rearranges indexes;
the composed space carries the actual order) decomposes into three claims:

1. ✅ **Partition changes the composite's effective storage order** — TRUE. The
   partitioned composite is block-major: the partitioned axis's block index becomes
   the OUTERMOST level of the flat composite layout, wherever that axis sat in the
   dense member layout. (This is what a block-tridiagonal-in-energy materialization
   needs, and the ⊕-composite provides it without touching member arrays.)
2. ✅ **The composed space carries the actual-order information** — TRUE, and
   ALREADY SHIPPED: `CoupledSpace` offsets-as-structure / `system_slices` `[M]` ARE
   that record. Repeated partitions build a ⊕-tree; the tree is the order record.
3. ⛔ **"The standard order becomes essentially meaningless" — REFUTED.** Partition
   is not permutation: ⊕-lifting axis k restricts k per member and touches no other
   axis — `V ≅ ⊕_b (A_1 ⊗ … ⊗ A_k|_b ⊗ … ⊗ A_n)` — so every MEMBER inherits the
   member-internal order unchanged, whatever k is. The leaf order is what
   kernels/einsums/sweeps consume (the user's "advantage in multiplying / setting
   up cross sections" — exactly; that is its job, and it keeps it at every leaf of
   the partition tree). ⟹ **Two ordering records, two owners, no conflict:** leaf
   dense order = the axes tuple (canonical family convention); composite order =
   the ⊕-tree with offsets (CoupledSpace). A leaf-layout CHANGE (e.g. JAX
   vectorization regrouping) is Realization freedom exercised only when a second
   layout exists, recorded as an explicit layout map — never a silent reorder.

**Round 4 sharpening — WHEN leaf order is consulted (user: "the leaf dense order is
the last thing to be realized… only needed after pencil → resolvent → partition,
finally lowering").** Right for the algebra, with one caveat on the data side. The
algebraic middle — pencil, resolvent, partition, schedule — manipulates operators
and spaces symbolically and consults INDEX SETS and ⊕-structure, never dense
layout. Dense layout is consulted exactly where ARRAYS touch the algebra, and there
are two such boundaries: (a) the DATA-IN boundary — XS libraries and fields
materialize as arrays from step one (`[M]` the mat_xs per-cell views are
`(ng, *spatial)`; every field mint follows the family convention), so the
convention binds early on the data side whether we like it or not; (b) the
LOWERED-SOLVE boundary — materialization (`as_matrix` C-order), kernels' einsums at
apply time. ⟹ **leaf order is a BOUNDARY convention, not a pervasive one**: the
axes tuple fixes it at the two array boundaries; everything between is layout-free,
so a future lowering that wants a different layout (JAX batching) renegotiates only
the boundaries — which is the precise content of "layout freed for Realization",
and confirms the user's instinct on the solve side while the data side explains why
the convention must exist from day one anyway.

**Round 5 — the order itself is CONVENTIONAL, and the tree's data layer already
agrees (user's early-binding review confirmed).** The (method, energy, spatial)
order was chosen for the angular solve's advantage — a Realization motivation that
leaked into the universal convention; any fixed order would have served, and the
math is order-free (⊗ symmetric up to natural iso). Two sharpenings keep this from
licensing drift: (a) the requirement is UNIFORMITY, not any particular order — F2b's
implicit-transpose hazard is the cost of per-site freedom, so the convention stays
one and explicit (the axes tuple), which is what turns "arbitrary early choice"
into "recorded, renegotiable convention"; (b) `[M]` the DATA layer is already
layout-free where it matters — MaterialXSField's per-material (ng, ng) matrices
"carry the source of truth"; the (ng, *spatial) per-cell views are CACHED DERIVED
views — so late/lazy layout binding is genuinely available at the view/field layer
when a consumer (JAX lowering) arrives, without touching primitive data. The
early-bound order is a VIEW convention, revisable at bounded cost.

⟹ **Filed as #394** (round 6, user ruling: no advantage in acting now; the JAX
lowering point is strong and a port is likely after the entire plan lands) —
cross-cited to #295 (layout Protocol sibling) and #197 (the Pint-at-JAX
milestone marker). The issue carries the two-boundary analysis and the
suggested shape for the eventual layout map.

## D. The structural-flaw register (round 2, adopted; `[M]` corrections)

F1 live + wider (2 unread-file weak-equality sites; #369 measure twin; CS1's
deterministic names = the hand-mint discipline done by machine, bridging to S3).
F2 headline live; the two-methods divergence is HISTORY (fixed 2026-08-04) — do not
hunt a live divergence. F3 confirmed and wider (4–5 doctrines / 6+ sites; direction
RATIFIED bulk→trace-doctrine, volumes onto the space, JAX argument; owner CS2/S2
with per-move ERR-067 checks; CS1 neither worsens nor fixes it). F4 measured at
state size; unrepresentable for axis-built spaces; densifier retires in CS2. F5 5/5
absent; CS1 delivers cone/basis-kind + EnergyAxis + declared V/V* hook; CS2/S4 owns
WithTrace + parity (stress case: the `_RadialCharacteristicSubSpace` split pair).
Severity: F1/F2 correctness; F3 correctness-adjacent; F4 structure; F5 capability.

## E. Phase mapping and judgment calls (round 2, RULED-agreed round 3)

Phase 0 → CS1 as item 0.8 in MIGRATION FORM (chartered doctrine + still-nominal
behavior, both stated). Phase 1 → SPLIT: primitive/of_axes/per-axis-metric/instances
= CS1; SN-mint rewiring + `_broadcast_metric` deletion = CS2; **⊕ NOT in CS1 —
CS2's opening step** (user agreed round 3). Phases 2–4 → CS2 (S2/S3/S4;
logged-then-flip already plan-of-record). Phase 5 → S7/Campaign 2. Judgment calls:
interning + content-eq LAYERED (construction discipline vs semantics — user agreed);
F3 → volumes onto the space; Phase-3 breakage → logged-then-flip.

## F. CS1 concrete scope

**Ships** — `orpheus/numerics/axis.py`: `BasisKind`, `Axis` (frozen; structural
eq/hash; ⚠ naming coordinated with the transport `Axis1D`/`AxisMesh` family — see
the rename issue), `EnergyAxis` (`from_grid(eg)` / `synthetic(ng)`; **identity =
ng + edges content — RULED Q2**, synthetic ng-only; docstring: the counting
theorem, the V/V* hook, the faces reading (EnergyGrid = 1-D energy mesh; groups =
cells; boundaries = faces; condensation = mesh-overlap map), the declared morphism
family). ⭐ Construction rulings from the battery review (2026-08-20, §T-R below):
weights are stored CANONICALLY — an all-ones array collapses to `None` at
construction, so the counting measure has ONE spelling and one identity (Q-T3;
kills the None-vs-ones twin AND keeps the no-densification memory claim clean);
non-finite weights REFUSED at construction (no non-negativity guard — CS2's
quadrature axes legally carry negative weights); the derived-name digest runs over
`(weights + 0.0).tobytes()` so `-0.0`/`+0.0` cannot mint two names for one measure
(Q-T1). ⚠ F10 trap, `[M]`: `mix.energy_grid` mints a FRESH `eq=False` object per
access (`is` and `==` both False on one mixture) — every grid comparison anywhere
in CS1 is on `edges` content, never on `EnergyGrid` identity, or two mints from
one mixture disagree inside a legitimate solve. `FunctionSpace`: `axes` (compare=False), `of_axes` (deterministic names;
weights never densified), per-axis metric path (Q4: build now, toy-tested),
`has_coordinate_cone`, `__mul__` axis-threading. `MaterialMesh.bulk_space` — the
cached mint, ⭐ GROUNDED 2026-08-20 as the **UNIFORM formula** (the property is
inherited by `SNMesh`/`DiffusionMesh`, so it must be honest on EVERY member — and
it can be, with a single generic body, Cardinal 2):
`of_axes(energy_axis, Axis("spatial", spatial_shape, weights=volumes, NODAL))`.
The degenerate carrier falls out as the quotient point (`[M]` its volumes are
`[1.0]` ⟹ weight 1.0, the normalized density convention); a genuine one-cell mesh
keeps `V ≠ 1` BY THE DATA (retrodiction A.5 row 4 realized mechanically); a meshed
carrier gets the honest scalar bulk `(ng, *spatial)` with cell-volume weights —
the SEED of CS2's single scalar-bulk mint (diffusion's `scalar_bulk` unification
stays CS2; in CS1 the only consumer is homogeneous). Energy arm (`[M]`
`Mixture.eg: ndarray | None`, `mixture.py:65`; the carrier stores
`materials: dict[int, Mixture]`, ng-validated at construction): all materials'
`eg` present AND edges content-equal ⟹ `EnergyAxis.from_grid`; else
`EnergyAxis.synthetic(ng)` — deterministic per carrier, and NO new
construction-time refusal in CS1 (grid coherence across materials is a
MEDIUM-level invariant; noted as CS1.5 design input). ⭐ RULED at battery review
(Q-T2/F8): the arm reads only materials REACHABLE from `mat_map` — `[M]`
`from_materials`' own docstring retains spectator entries "unused by the single
cell", and a spectator with `eg=None` must not flip the axis identity of a problem
it does not touch (the leak principle: the mint consults exactly its defining
data; the ng-validation guard may stay broader — a guard is not an identity
input). ⚠ F1, `[M]`: ALL 12 shipped `get_mixture` pairs have `eg is None`, so the
`from_grid` arm's witness is MANUFACTURED (`dataclasses.replace(mix, eg=…)`, the
`test_homogeneous.py:415-417` idiom) — every energy-arm gate parametrizes over
BOTH arms and the battery mutates each arm separately (M15/M16, vv#17). ~~+ the §B2 micro-cleanup
(from_materials stops minting `AxisMesh(edges=[0,1])`)~~ ⛔ SUPERSEDED — the Q11
ruling (round 6) dropped it and this paragraph was not reconciled until the
2026-08-20 grounding pass: CS1 does NOT touch `from_materials` internals; the
fake mint dies in CS1.5, nowhere else. Item 0.8 migration-form docstring.
`Field.cone_violations` consults `has_coordinate_cone` (False ⟹ refuse; None ⟹
legacy behavior, documented).

**Homogeneous rewiring** (unchanged from round 2): thread `V = mat_xs.mesh.bulk_space`
through C (via the extended `from_mesh` mesh-default chain), IsoS/IsoN2N
(`space=V`), F (`from_solver_data(…, space=V)` — widen FORCED, `FullFieldSpace | None` →
`FunctionSpace | None` on field + param + the docstrings that say "composite
full-field space"; **rename `full_field_space` → `space` RULED-provisional Q5**:
proceed, judge the end result at review; one §6b batch). No default-derivation on
`from_solver_data` (wrong-family hazard). Both production `basis_shape=(ng, 1)`
spellings deleted (`solver.py:194/:202`; test-side spellings are judged per-test —
mirror tests migrate to domain-derivation, bare-op tests KEEP the explicit shape
because the bare path is their subject until CS4).

⭐ **Rename batch GROUNDED 2026-08-20 — supersedes the round-2 site list above,
whose `scattering.py` rows are `ScatteringOperator`'s OWN field (a second class
carrying the same field name), not F callers.** The §6b batch must pick its scope
explicitly:

- **F-scope (the Q5 ruling's letter)** — `fission.py` (field :182, readers
  :233/:238, classmethod :243/:255, docstrings :174/:250) + keyword callers `[M]`
  `sn/solver.py:1341`, `sn/solver.py:2805`, `diffusion/solver.py:248`,
  `tests/sn/architecture/test_monomorphic_leaves.py:519`.
- **+S-scope (⭐ RULED IN — user, 2026-08-20: "space and F+S in scope" —
  naming-consistency: fix the off-pattern family member in the same change; S's
  TYPE stays the narrow `FullFieldSpace | None`, only the NAME moves)** —
  `scattering.py` (field :410, readers :429/:434, classmethod :766/:785,
  self-forwards :1013/:1054, docstrings :403/:770-771) + keyword callers `[M]`
  `sn/solver.py:1337`, `sn/coupled_system.py:503`,
  `tests/sn/operators/test_psi_half_coupling.py:2873`,
  `tests/sn/architecture/test_monomorphic_leaves.py:516`.
- ⚠ **NOT in the batch**: `tests/sn/solve/test_si_gate_dispatch.py:63` — its
  `full_field_space=None` fakes the MESH attribute (a `SimpleNamespace` mesh
  mimic), not the operator field; and every `*.full_field_space` MESH-property
  read (`SNMesh`/`DiffusionMesh`/consumers) keeps its name. The batch grep must
  discriminate operator-field sites from mesh-attribute sites.
- Bare callers (no keyword) are untouched by the rename:
  `tests/homogeneous/test_homogeneous.py:285/:359`,
  `tests/diffusion/test_operators.py:650`.

**Geometric-axis note (round 3):** the transport `Axis1D` family already factors
GEOMETRY per axis (`spatial_shape`/`face_labels` as pure functions on axis tuples
`[M]` `transport/mesh/axis.py`). CS2's `SpatialAxis.from_mesh` may therefore choose
rank-d single axis vs per-dimension rank-1 axes (grand report §15.1's
`L = Σ_axis D_axis ⊗ Ω_axis ⊗ I_g` leans per-dimension) — recorded as a CS2 design
point, not decided here.

**Blast radius (§8):** values NONE by the theorem (homogeneous suite byte-stable);
guards newly active on `C − K_iso` + inner iso-sum + `M⁻¹ @ F` (one shared instance
⟹ trivially agreeing; new capability = refusal); `.H` flip bit-identical, gated per
vv#19 with the loaded/blind PAIR — ⛔ SHARPENED at battery review (F2, `[M]` table
in `scratch/cs1_verification_plan.md` §1): a SCALAR metric commutes with every
operator (vv#12's commutator exactly — `G = cI ⟹ G⁻¹AᵀG = Aᵀ`), and EVERY
production-reachable homogeneous space has a scalar metric (rank-1 spatial axis ⟹
scalar weight; energy counting BY THE THEOREM), so the quotient-point-vs-one-cell
`V ≠ 1` difference can NEVER move `.H` — `[M]` `w=[2.,2.]` is bit-identical to
the None path; only per-GROUP `w=[2.,5.]` moves it (`[-0.08,0.2] → [-0.38,0.2]`,
and only component 0 moves, so the control asserts the WHOLE vector). The control
(gate D4b) is therefore a deliberately NON-PHYSICAL per-group-weighted toy — its
docstring must say the counting theorem forbids this on a real EnergyAxis — and
the quotient-vs-one-cell distinction is guarded by SPACE IDENTITY alone (gate
B9, plus the M17 MUST-STAY-GREEN proof that no value gate can see the measure);
carve-time nexus enumeration of the four classes' construction sites + the
degenerate-path geometry readers (§B2).

**Done-when:** (1) the 4 strict-xfails
`test_model_generic_leaf_declares_a_space[C-2g,C-4g,F-2g,F-4g]` DELETED, replaced
by the positive homogeneous floor (all five operators + K report the SAME space;
IsoS/N2N included) — ⭐ CORRECTED 2026-08-20 (grounding pass; the §1
done-when-predicate class): only the **C** rows XPASS mechanically (the test's own
`from_mesh(field, mat_xs.mesh)` call rides the extended chain; strict ⟹ forced
deletion), while the **F** rows CANNOT XPASS — `[M]` the test's F arm calls
`from_solver_data(mat_xs=mat_xs)` BARE, and the ruled no-default-derivation keeps
that `None`. The F rows are deleted in the SAME §6b commit on a different warrant:
their mirrored production line (`solver.py:193`) now passes `space=V`, so the bare
mirror would pin a RETIRED idiom, and the floor is the successor gate for all
four; (2) both production `(ng, 1)` spellings gone — ⚠ battery finding M23: this
item has NO possible runtime witness (a leftover explicit shape is
value-identical), so it is a GREP obligation on the 3b commit, never implied
coverage; (3) refusal witnesses
RED (2g-vs-4g sum → the pinned "equal domains" fragment; `M⁻¹(2g) @ F(4g)`) + the
modal-axis cone-refusal witness; (4) homogeneous suite byte-stable; pyright
`orpheus` = 1; (5) `spaces.rst` SEEDED — Key Facts; axis taxonomy; the counting
theorem (covariant/contravariant derivation); **the collapse doctrine §B (three
clauses, the consultation principle, the Bateman one-cell reading)**; the quotient
family incl. B=0 fiber + P_N-as-angular-buckling; dev-history + the owed
`field_algebra.rst` micro-edit; Sphinx `-W` clean.

**Fences:** no Spatial/Quadrature/Harmonic axis classes (CS2; quotient point =
generic instance, re-homed as `SpatialAxis.quotient_point()`); no ⊕ (CS2 opening);
no identity flip (S3); no metric relocation / mint rewiring / `_broadcast_metric`
deletion (CS2); no Optional→mandatory flip or dispatch collapse (CS4); no 𝔽/ℂ (S5);
condensation/partition declared only; the RegionCarrier Protocol (Q10) recorded,
not built.

**Execution shape (post-ruling; surgical posture; ⭐ step 3 SPLIT 2026-08-20 by
the grounding pass — the rename/widen batch is behavior-neutral and separable,
while the `from_mesh` chain extension ALONE XPASSes the C xfail rows, §6b-fusing
the chain with the xfail handling):** 0. test-architect dispatch + enumeration
memo → 1. `axis.py` + intrinsic tests → 2. `of_axes` + per-axis metric + cone
metadata + `__mul__` threading + 0.8 docstring → 3a. F widen +
`full_field_space` → `space` rename batch (scope per the grounded list above;
suite green, zero behavior change) → 3b. `bulk_space` + `from_mesh` chain
extension + solver rewiring + xfail deletion + floor + witnesses + byte-stability
run (~~micro-cleanup~~ ⛔ CS1.5 per Q11) → 4. `cone_violations` wiring + modal
witness → 5. `spaces.rst` seed + dev-history + campaign-plan §2 ledger row.

## H. The design lens this campaign runs under (articulated round 5, user + main agent)

The recurring defect across kernel/operator (CS4), the space layer (CS1/CS2), and
the carrier (Q11) is ONE defect, and the user named its test: **a concept has a
LEAK when it can be asked a question that is meaningless for it** — the answer must
be "the method is not there" (absence by construction), never "Error: not
applicable". The principle:

> **Every concept carries exactly its defining data — everything it needs, nothing
> else. Every dependency between concepts is a NAMED CONSTRUCTOR: the moment of
> binding is a spellable object in the code. Polymorphism and derivation happen at
> that construction, never at application.**

A weld = two stages' data fused in one class ⟹ the binding moment between them
becomes unspellable ⟹ downstream, the missing name resurfaces as apply-time
dispatch, fake degenerate data, or guards that cannot be stated. Instances:
`ScatteringOperator` welding kernel data to apply-time carrier dispatch (the
binding kernel × frame × space → operator was unspellable — CS4 un-welds);
monolithic shape-tuple spaces welding axes (partitions/strategies unspellable —
CS1/CS2 un-weld); `MaterialMesh` welding medium to mesh (the pullback unspellable;
the homogeneous member forced to fake geometry — Q11 un-welds). The detector is
the report Part IV stage-inversion table; this lens is its constructive dual.
(Promotion into `coding-elegance` ACCEPTED by the user, round 6 — with the caveat
"it might even be there already, spelled in a slightly different way": at
distillation time, READ the skill for an existing spelling first and sharpen it
rather than minting a duplicate clause. Owed at the campaign's distillation pass,
not before.)

## G. Fork register

| # | fork | state |
|---|---|---|
| Q1 | axis order | **RESOLVED-proposed** (round 2): storage order = axis order = the uniform (method, energy, spatial, [moment]) convention; §C round-3 partition analysis confirms leaf-order and ⊕-order are separate records with separate owners |
| Q2 | EnergyAxis identity | **RULED** (round 3): ng + edges content. Sharpened by the user: the boundaries are the energy mesh's FACES — EnergyGrid is a 1-D mesh in energy; `from_grid` is axis-from-mesh, symmetric with `SpatialAxis.from_mesh` |
| Q3 | names | **RULED** (round 3): `EnergyAxis` (not `GroupAxis`). Geometry `AxisMesh` ruled a MISNOMER (declares geometry, creates no mesh) — rename issue filed (86 construction sites `[M]`; coordinate with numerics `Axis` naming so space-factor vs geometric-axis vocabularies stay distinct) |
| Q4 | per-axis metric now | **RESOLVED-proposed** (round 2): build in CS1, toy-tested |
| Q5 | F field rename → `space` | **RULED-provisional** (round 3): proceed; judge the end result at review. **⭐ RULED FINAL (user, 2026-08-20, post-grounding): name = `space`, scope = F+S** — ~~S's type stays the narrow `FullFieldSpace \| None` (only the NAME moves)~~ ⛔ SUPERSEDED at 3a review (user, 2026-08-21): **S's type ALSO widens to `FunctionSpace \| None`** — `[M]` the family census (run answering the user's challenge) showed C/IsoS/IsoN2N all already carry the WIDE type, so the narrow-S proposal left S the family's only off-pattern member on the type axis (the exact shape the naming-consistency rule prevents), and the narrow annotation was load-bearing for nothing (no test pins it; the real guard is instance AGREEMENT in OperatorSum, not the annotation's family). Landed with the amend commit after `e8769897`. Mesh properties (`SNMesh`/`DiffusionMesh.full_field_space`) and the `FullFieldSpace` class keep their names (they name objects, not roles) |
| Q6 | `(ng, 1)` explicit point | **RESOLVED-proposed**, principle upgraded round 3 (§B clause 1: pairing-consulted normalization) |
| Q7 | trivial angular slot | **RESOLVED-proposed**: ABSENT on scalar-family spaces (§B clause 3 — the convention lives on ι) |
| Q8 | 0.8 doctrine rewrite in CS1 | **RULED (round 5): CS1 executes it, migration form.** (Round-4 explanation retained: `space.py:145-150` still TEACHES the overturned doctrine ("two copies of ℝⁿ are 'the same' space regardless of which inner product is installed" — report correction #16). Rewriting it to state ONLY the target doctrine would lie forward (behavior stays nominal until S3); leaving it lies backward. MIGRATION FORM = state both: the chartered doctrine (identity = axes' structural content + tags; metric differences imply space differences) AND the current nominal realization AND the S3 flip plan AND the axis-built bridge (derived names). Proposal: CS1 executes it (CS1 edits space.py anyway) |
| Q9 | meshless homogeneous | **SUPERSEDED round 4 → Q11**: the axis is NOT welded to Mesh (multi-generator constructors — stands); the carrier resolution upgraded from keep-with-micro-cleanup to the Medium finding (§B2) |
| Q10 | `RegionCarrier` Protocol | **FOLDED into Q11** (round 4): the Protocol vs union-in-base choice is internal to the Medium carve |
| Q11 | **the Medium carve** | **RULED round 6 (delegated to main agent by the user; decision + reasoning):** chartered as **CS1.5 — its own phase between CS1 and CS2**. Why there and not elsewhere: (a) CS2's axis constructors must bind to the RIGHT generators from day one — `SpatialAxis.from_mesh` forgets from the MESH, the quotient axis from the MEDIUM's symmetry — so Q11 lands before CS2 or CS2 binds to the pullback and re-points later (double-touch); (b) not BEFORE CS1: CS1 is small, xfail-pinned, and the ruled order (CS3→CS1→…) stands — and Medium's quotient member wants CS1's axis machinery already present to mint its space directly; (c) not a CS2 opener: CS2 already opens with ⊕ (ruled round 3), and stacking two openers merges blast radii against §6b. Content as §B2-round-5: mint `Medium` ABOVE (minimal start: infinite member = full-symmetry marker + one mixture; 1-D member = interface positions + per-segment materials; nothing more until a consumer); `MaterialMesh` keeps its name as the PULLBACK; named binding + mesh-conformity guard; homogeneous re-homes to `Medium.infinite_homogeneous(mix)`; `bulk_space` MIGRATES carrier→Medium in CS1.5 (one planned move). Round-3's micro-cleanup is DROPPED as subsumed — ⚠ explicitly: CS1 does NOT touch `from_materials`' internals; the fake `AxisMesh(edges=[0,1])` mint dies in CS1.5, nowhere else (no double-touch, nobody "cleans" it twice) |

---

## Appendix A — the collapse doctrine: full derivation record (the archivist's source)

**Purpose.** The user's preservation instruction (round 5): *"Make sure we don't
lose it so the archivist can write it later. It was hard to steer until we got it
out."* This appendix is the self-contained source for the `spaces.rst` section —
per Sphinx-as-brain it records the refuted versions WITH their refuting questions
(what was tried and failed is part of the knowledge), the final doctrine, and the
retrodictions. Derived dialectically 2026-08-19/20 (user ⇄ main agent, rounds 2–4
of the CS1 design discussion).

### A.1 The question

The homogeneous infinite-medium solver, diffusion, SN, and PN differ in which
tensor factors their spaces carry. Two prior doctrines in the operator-machinery
report were in TENSION: §I.11 says the homogeneous quotient KEEPS a one-point
spatial axis with unit measure ("per unit volume" = the normalized quotient
measure); §I.9 says a retract's codomain DROPS the axis (scalar space is NOT
angular space with a trivial axis — conflating them is the 4π-error source). Both
collapses are rank-one. Which collapses leave an axis, and why?

### A.2 Version 1 — compactness (REFUTED)

*Proposed (round 2):* a collapse over a COMPACT domain integrates (measure
consumed, axis drops — angle, S²); over a NON-COMPACT domain it must normalize
(axis persists as the quotient point — space, ℝ^d).

*Refuted (round 3, by the energy/Bateman question):* energy would sit on the
"integrate ⟹ drop" side, yet the one-group flux manifestly KEEPS its axis
(`(1, *spatial)`), and must — σ̄'s defining interval is problem data. Compactness
decides integrate-vs-normalize at best; it cannot decide persist-vs-drop.

### A.3 Version 2 — "effectively finite measure" for energy (REFUTED)

*Proposed (round 3):* patch energy in as "effectively finite-measure (grid-topped,
flux integrable)".

*Refuted (round 4, by the user):* the energy domain (0, ∞) has INFINITE Lebesgue
measure — there is no upper limit; the grid top is practicality. The correct
question is the user's own formulation: *does the integration of the infinite tail
converge?* — i.e. L¹ membership of the ADMISSIBLE FIELD CLASS. The discriminator
is a property of the fields the physics admits, never of the bare domain.

### A.4 The doctrine (v3, standing)

**Fork 1 — integrate or normalize? Decided by integrability of the admissible
field class over the collapsed domain.** The structural mechanism:

- A quotient BY A SYMMETRY GROUP forces the fields CONSTANT along group orbits.
  If the group is non-compact (translations of ℝ^d), orbits have infinite Haar
  measure and a constant on infinite measure is never L¹ — integration is
  impossible STRUCTURALLY, whatever the physics. Only the normalized average
  ("per unit orbit measure") survives.
- If the group is compact (rotations, S² an orbit), Haar measure is finite,
  integration is always available, and the total (4π) is canonical.
- With NO symmetry acting (energy), nothing forces constancy; the physics decides
  integrability. For neutron spectra it holds: χ decays super-exponentially above
  source energies; the thermal Maxwellian → 0 as E → 0; the tail integral
  converges. The practical grid top E₀ is the LIBRARY's truncation assertion
  (data-layer), not a space-layer fact.

**Fork 2 — does the axis persist? Decided by CONSULTATION: the axis survives iff
the collapse leaves data the surviving family still consults — in its own pairing,
or in identity/guards. A collapse that leaves only re-embedding conventions puts
them on the (π, ι) arrows and drops the axis (the family changes — a §I.9
retract).** Three clauses:

1. **Non-compact-orbit quotient ⟹ normalize ⟹ axis persists** as the quotient
   point with unit weight — the density convention ("per unit volume") is
   consulted by the member's OWN pairing. `[M]` live: `volume_measure` = 1.0 is
   consumed by the reaction-rate functionals inside the homogeneous solve
   (`material_mesh.py:253/:408`, `reaction_rate_functional.py:224-228`). In the
   modal family the surviving datum is the mode parameter (buckling's B).
2. **Partition-integration of an L¹ class (no symmetry) ⟹ per-cell integrals ⟹
   nodal mesh-axis persists to ONE cell** — the partition (boundaries, weighting
   spectra) is problem data consulted by identity/guards and the V/V* pairing.
   EnergyGrid is a 1-D mesh in energy: boundaries are FACES, groups are cells
   (user's insight, ruled into EnergyAxis identity), condensation is the
   mesh-overlap map. The one-group member is the one-CELL energy mesh — edges
   retained — and its ⟨σ̄, φ⟩ pairing is exactly what Bateman/depletion consumes.
   Likewise a genuine one-cell slab keeps its axis with weight V ≠ 1
   (distinguished from the quotient point by measure).
3. **Whole-domain integration over a compact canonical orbit ⟹ axis DROPS** —
   S² is universal; nothing problem-specific survives on the axis; the 4π
   rebroadcast convention ("uniform in flux" vs "uniform in current") lives on
   the EMBEDDING ι, consulted only when LEAVING the family. Scalar spaces carry
   no angular slot and are a different family.

The one-line tests: *can the admissible fields be integrated over the collapsed
domain?* — and — *is the surviving convention consulted inside the family (axis)
or only at re-embedding (arrow)?*

### A.5 Retrodictions (the doctrine against the tree, all `[M]`)

| tree fact | clause |
|---|---|
| homogeneous shape `(ng, 1)` — spatial point present, unit volume consumed by functionals | 1 |
| scalar family `(ng, *spatial)` — NO angular slot; 4π conventions live at embeddings | 3 |
| one-group `(1, *spatial)` — energy axis persists with its edges | 2 |
| one-cell meshes keep weight V ≠ 1 (≠ quotient point) | 2 |
| partial currents: hemisphere collapse parameterized by n̂ — the parameter survives on the FACE structure, the angular content on the arrow | 3 (+ face axis) |
| buckling member: size-1 MODAL spatial axis carrying B; live angle; 𝔽 = ℂ | 1 (modal branch) |

### A.6 Consequences carried elsewhere in this draft

The explicit `(ng, 1)` point axis (§B standing consequences); the B=0-fiber
continuity of the buckling family; the P_N hierarchy as the ANGULAR buckling
(irreps of rotations ↔ irreps of translations) with the trivial angular slot
absent on scalar members; the "what a degenerate non-compact axis stores" answer —
the density convention as a first-class composable object (pairing-automatic
rates; family coherence; homogenization-boundary Jacobians); the Medium hierarchy
(§B2) — the quotient spatial structure is GEOMETRY's symmetry data, which is why
the medium (not the mesh) is the homogeneous solver's honest home. Round 6
sharpening (user): symmetry is MONOTONE down the lattice (arrows preserve or
reduce, never increase) — clause 1's quotient consumes the MEDIUM's surviving
symmetry, G_medium = G_geom ∩ Stab(assignment); the infinite homogeneous medium
is the member whose assignment stabilizer is everything.

### A.7 For the archivist

Target: `docs/theory/foundations/spaces.rst` (the CS1 seed, §F done-when). Carry:
the tension (A.1), BOTH refuted versions with refuting questions (A.2/A.3 — the
dialectic is the pedagogy), the two forks + mechanism (A.4), the retrodiction
table (A.5), and cross-references to §I.9/§I.11 of the operator-machinery report
(this doctrine RESOLVES their tension — say so explicitly). Key Facts entries:
the two one-line tests; "EnergyGrid is a 1-D mesh in energy"; "the quotient point
records the density convention".
