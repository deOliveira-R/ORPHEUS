---
name: campaign-290-diffusion-integration-map
description: Ground-truth architecture map for the #290 diffusion-integration campaign — claim audit, typing shape, SN mirror pattern, D data seam, retirement blast radius. TRANSIENT.
metadata:
  type: project
---

# #290 diffusion integration — ground-truth map

**Line numbers current at HEAD `d2a2a0c` (2026-07-03, branch `refactor/pyright-burndown` ≡ main).
Re-derive via Nexus if drifted. DELETE this file when #290 merges** (per L-003: transient
carve-maps are deletable archaeology; durable rulings graduate to the theory page / AGENT.md).

**Why:** #290 is THE active campaign; this memo saves the next dispatch the ~40-call re-derivation.
**How to apply:** trust the SHAPE claims; re-verify any `file:line` before acting.

## A. Claim audit — every #290 claim CONFIRMED at d2a2a0c

No drifted claims. Highlights: island = 3 files/443 lines exactly; scipy A_op `solver.py:192-196`;
flip trick `:244`; `print` in converged `:323`; ONE real pyright error `solver.py:194` ("No
parameter named matvec", pyright 1.1.410 — 2 additional cross-root import artifacts to discount);
`KEigenvalue` names 5 families `iteration.py:997/:1013`; stub realizer registered "diffusion"
(`diffusion/boundary_realizer.py:31`, #182 OPEN); scalar role triple complete (`ScalarFlux`
fields/scalar_flux.py:108, `ScalarSourceSink` source_sinks/scalar_source_sink.py:89,
`ScalarResidual.from_balance` residuals/scalar_residual.py:69/:114 — note: sibling subpackages,
not one package); `FullField.bulk` admits ScalarFlux (full_field.py:134-136, isinstance-BulkField
only :211); #289 role-erasure guard live in `evaluate_residual` (sn/solver.py:277-281);
BoundaryField SN-bound (`_bases.py:685` mesh:SNMesh, `:694` TraceSpace+layout required, factories
read `mesh.trace` :798/:818/:843); TraceSpace quadrature-coupled (omega_dot_n, trace_space.py:340);
MaterialMesh has NO `.trace`; Mixture has NO transport channel; #215 trap verbatim
scattering.py:1305-1320; DSA architecture operator_algebra.rst:4210-4311; baseline json = {transport:1},
diffusion absent; `"orpheus/diffusion"` in pyproject ignore :119.

## B. Typing shape (post-C5) for the scalar composite + A_diff

- `LinearOperator(Protocol[Domain, Codomain])` operator.py:403. Mandatory `apply` :526. NO
  capability tags — recursive predicates `is_invertible`/`is_adjointable` (:495/:513, default
  False) + `SupportsInverse`/`SupportsAdjoint` (:889/:904) + `NotInvertible`/`MissingAdjoint`.
  Dunders have real bodies on the Protocol (sum/product/scaled/pow/H/as_matrix).
- KEigenvalue contract: A.is_invertible must be True (iteration.py:1113 raises); the built
  `A.inverse().apply` MUST accept `initial_guess=` kwarg (SourceIteration calls it :677;
  InverseOperator.apply operator.py:1822). Dense route: MatrixInverseOperator (homogeneous
  precedent, eigenvalue.py:270-276).
- `FunctionSpace(Generic[Carrier])` space.py:104 (PEP-696 default Any); identity=(name,shape);
  carrier-generic metric surfaces apply_metric/apply_inverse_metric/inner_product with
  structured override (FullFieldSpace :267-305 per-block dispatch).
- `FullFieldSpace.from_blocks(bulk_space, trace_space: TraceSpace)` full_field_space.py:196 —
  the trace-block ANNOTATION is TraceSpace; TraceSpace requires layout+omega_dot_n kw_only
  (:332-333). Scalar-trace decision cascades here + `FullField.mesh` (reads boundary.mesh,
  typed SNMesh, full_field.py:238-254) + `to_flat` flat-backing note (:356-378).
- Field leaves: Field ABC (values,space,UNITS) field.py:143; ScalarField gives rank-adaptive
  `(ng, *spatial)` + `_SPACE_NAME` factories `_bases.py:392-477`; BulkField mesh-identity gate
  :151-166; role transitions ONLY via from_balance (named composition).
- Two driver anchors: `EigenvalueSolver` 5-method protocol (eigenvalue.py:64-143; NOTE
  DiffusionSolver is NOT ProductionRateSolver — legacy un-normalised window :221-227, adopt
  `compute_production_rate` in integration) OR the KEigenvalue (A,S,F) triple route.

## C. SN mirror pattern (mesh → spaces → operators → realizer → solver)

1. Declaration: `Mesh1D(edges, mat_ids, coord, bc_left/right)` or `AxisMesh(edges, bc_low/high)`
   (transport/mesh/axis.py:222 — NOT geometry/); BCs are `BoundaryTraceLaw` descriptors
   (geometry/boundary/).
2. `MaterialMesh(mesh, materials)` material_mesh.py:122 — ng/spatial_shape/volumes/
   volume_measure/`material_xs_field()`.
3. `SNMesh.from_material_mesh(mm, quadrature, scheme=…)` augmented_mesh.py:727 (data/behavior join).
4. Trace built ONCE in mesh init: `TraceSpace.from_quadrature_and_layout` augmented_mesh.py:448,
   property `.trace` :779.
5. BC realization AT MESH CONSTRUCTION: per face `SNMethodSpace.for_face(mesh, quad, face, trace)`
   :502-506 → `SNBoundaryRealizer().realize(law, method_space)` :508 (realizer
   sn/boundary/realizer.py:126; `_as_boundary` stamps BlockRole.BOUNDARY :107-123).
6. `SNMesh.full_field_space` :800 = from_blocks(G_bulk=V·w_n, G_trace=|Ω·n|·w_n).
7. Operators off the mesh: SNSolver.__init__ solver.py:590-763 — `mat_xs = mesh.material_xs_field()`
   :658; S/F `.from_solver_data(mat_xs, …, full_field_space)` :708-717; L+C = StreamingOperator +
   MultiplicationOperator(total_cross_section_field) :725-731; `_within_group_triple` :160-220
   returns variadic (L+C, S, B).
8. Drivers: SourceIteration/Krylov, rhs = q + Σ gains.apply (iteration.py:666-668), resolvent
   apply :677 (THE DSA plug point, before convergence test :706).
9. `solve_sn(materials, mesh, quadrature)` :1466 → `_as_sn_mesh` → SNSolver → power_iteration :1573.

## D. Data seam — D = 1/(3Σ_tr)

- Mixture channels mixture.py:57-65: SigC/SigL/SigF/SigP/SigT (NG,), **SigS = list per Legendre
  order** (P1 = SigS[1]), Sig2, chi, eg. `from_dense_channels` :228. NO transport channel anywhere
  in data/ (grep-verified).
- ENDF path: compute_macro_xs n_legendre=3 → SigS[0..2] (:466-495). Synthetic: make_mixture has
  `sig_s1` param; xs_library regions A–D all set it (μ̄ 0.05/0.60/0.10/0.30).
- Seam = derived property: `Σ_tr = SigT − rowsum(SigS[1])` (outflow transport approx), with
  len(SigS)==1 ⇒ Σ_s1≡0 ⇒ Σ_tr=Σ_t EXACTLY (correct isotropic limit, not a fallback). No new
  stored channel needed.
- Island today: hand-entered MATLAB constants `_default_xs` solver.py:88-106; face-arithmetic-mean
  of Σ_tr then D=1/(3σ_face) :180-186 (vs harmonic-mean-of-D — a discretization decision to make
  consciously); `e_per_fission` hardcoded :368 (power normalization post-processing).

## E. Retirement blast radius (graph + grep + constructors + docs)

- **4 consumer packages**: orpheus/diffusion itself; tests/diffusion (4 files, all via
  CoreGeometry/TwoGroupXS/solve_diffusion_1d); examples/diffusion/demo_diffusion_1d.py:14,24
  (zero-arg call + MATLAB keff 1.022173 pin); **orpheus/derivations/continuous/cases/diffusion.py:152-165**
  (derive_2rg Richardson RUNS the island solver on `_richardson_cache` miss — #290 does NOT list
  this consumer; `solver_cases()` :208-220 already ".. deprecated:: Phase 1.2 … Phase 2 deletes it";
  rewire test_spatial_convergence_reflected → derive_2rg_continuous, then delete).
- Docs nodes: theory/diffusion_1d.rst:62-93 (BC_REGISTRY/CoreGeometry — the :63 cross-ref target
  `~diffusion.solver.CoreGeometry` is ALREADY malformed, missing `orpheus.`), :604/:622/:632;
  api/diffusion_1d.rst (automodule page); api/geometry.rst:76-103 (BC_REGISTRY family doc names
  DiffusionSolver); api/numerics.rst:114; theory/homogeneous.rst:1635; theory/verification.rst:87
  (wrong name `orpheus.diffusion_1d`); docs/conf.py:56 (nitpick comment).
- BC_REGISTRY is a 4-member legacy family (cp/solver.py:133,420; mc/solver.py:162,209;
  moc/geometry.py:298,453) — diffusion's exit leaves 3.
- scipy A_op/bicgstab: confined to solver.py:16,192-196,285-291 — no external consumer.
- Stale TODOs to sweep: tests/diffusion/test_diffusion.py:9-13 + derivations diffusion.py:117-119
  & :196-198 (same "no theory page (#35)" claim; page EXISTS, 737 lines).

## F. Surprises / decision points for the user

1. **Derivations coupling** (E above) — the un-listed 4th consumer.
2. **Marshak-vs-zero-flux**: #182 stub specs Marshak/Robin; the island + ENTIRE 581-line analytic
   suite intentionally pin ZERO-FLUX Dirichlet (diffusion_1d.rst:74-87 "pure sinusoids without an
   extrapolation-length fudge"). Realizing Marshak changes every reference. User decision.
3. Island's reflective-BC arm has ZERO test coverage (all tests vacuum; "reflected" = reflector
   REGION).
4. scalar_flux.py staleness beyond #290's list: self-contradictory module docstring (:32-38 says
   MaterialMesh since #267 vs :51-62 SN-bound), prose-only TYPE_CHECKING SNMesh import :100-101,
   false `_check_partner`-override claim :40-45 (lives on BulkField now).
5. #242 (diamond.py:645 deferral target for the generic advection–reaction diagonal) CLOSED as
   documented-split — no open tracker; folds into #290's discretization choice. Live hook =
   scheme.py:706-710 (schemes take arbitrary reaction_xs; names diffusion-removal + DSA).
6. DiffusionSolver not ProductionRateSolver; conditions by `fi /= max|fi|` inside
   solve_fixed_source (solver.py:293-294) — adopt production-rate contract (ties to #270).
7. **AGENT.md durable-shape correction needed** (I cannot edit it): "`_reflect_outflow_into_inflow`
   RETIRED" is imprecise — the S+B FOLD is retired and the DRIVER route no longer uses it, but the
   helper SURVIVES (sn/solver.py:1408-1459) for the final eigenvalue reconstruction sweep + the
   octant-group G-S resolvent, delegating to `SNBoundaryOperator.reflect_inflow_inplace` (same B).
