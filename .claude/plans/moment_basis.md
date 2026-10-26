# The spatial-moment basis: from a layout policy to a numerics basis (#423, layer DAG ruling 19)

Living plan, opened 2026-10-25. The ontology is being searched (`plan-authoring` §0), so nothing here is implemented until the user rules the plan polished.

Parent: `layer_dag.md`, ruling 19. The context memo is `scratch/layer_dag/moment_basis/context.md` (the explorer, 2026-10-25), with its AST census in `census.out` beside it.

## The finding

`orpheus/numerics/moment_layout.py` encodes a polynomial cell basis that numerics does not hold.
- `[M]` 14 production import statements in 11 files: transport 6, sn 7, derivations 1. None is inside numerics.
- The spatial schemes spell the basis as two scheme constants: `spatial_basis_per_axis` (1 for DD and Step, 2 for LD) and `theta` (LD's 1/3).
- The scheme mints the moment axis as `Declared(Modal(e_0))` (`transport/spatial/scheme.py:1726-1740`), so the axis points at no basis.
- Four more sites spell the same shape without importing the module:
  - `scheme.py:518` writes `p**(ndim-1)`;
  - `fields/_bases.py:533` recovers `per_axis` with a float root;
  - `_ubld.py:144` decodes the Kronecker index;
  - the LD static methods hard-code `2`.

## The object, mathematically `[R]`

The LD cell expansion uses the Legendre polynomials on the reference interval ξ ∈ [−1, 1]:
- they are normalised by `P_k(1) = 1`, with degree k < p per axis;
- the per-volume Gram is `diag(1/(2k+1))`, which is `[1, 1/3]` for p = 2;
- `[M]` θ = 1/3 is that Gram: a 4-point Gauss–Legendre sum gives `[1, 0.3333]`, equal to `mass_1d`.

The d-dimensional cell basis is the tensor product of d copies, Kronecker ordered with axis 0 outermost. What `moment_layout` and the schemes carry today is each a property of that object:

| today | the basis fact |
|---|---|
| `AVERAGE_MOMENT = 0` | the constant mode is the unit vector e₀ of the tensor basis |
| `cell_moment_count(p, d) = p**d` | the dimension of the tensor basis |
| `face_moment_count(p, d) = p**(d−1)` | the dimension of the face trace: the tensor product with the face-normal factor removed |
| the face trace `B(+1) = [1, 1]`, the inflow weight `B(−1) = [1, −1]` | the 1-D basis evaluated at the interval's endpoints |
| the octant sign map `(−1)^k` | the action of the mirror ξ ↦ −ξ on the basis (its parity) |
| `theta = 1/3` | the 1-D Gram |
| `SPATIAL_MOMENT_AXIS_LABEL` | the label of the axis the basis heads |

The candidate primitives are a 1-D interval Legendre basis and a tensor product of bases, with the cell basis as their composition (`coding-elegance` Pattern 5).
- The existing `LegendreBasis` is the angular `P_ℓ(Ω·e_a)` on `S²/O(2)`. Reusing it is refuted: at d = 1 its shape `(2,)` equals the spatial one, so an angular head would accept the spatial axis.
- No tensor product of bases exists today. `TensorProductSpace`, `TensorProductOperator`, `Product` and `ProductMeasure` do.

## Not basis facts

- **The tail rule** (`face_moment_tail`, `spatial_moment_tail`: "append a trailing axis iff its width > 1") is a space-composition rule. It is applied to cell and face counts alike, and it keeps a degree-0 space byte-identical to the scalar one.
- **The rank predicates** (`is_moment_valued_by_rank`, `is_moment_valued_by_flat_rank`; 5 production uses) test a buffer's shape against the space. They are not a property of the basis.

## Hazards (`[M]` the explorer)

1. Moving the axis from `Declared` to a basis- or frame-built origin adds `support` and `layout` to its identity. Every LD-widened space then gets a new content digest, and the traced memos regenerate once.
2. `FromBasis.kind` is a closed `match` over four layouts, so a new layout joins it and the `CoefficientLayout` union. `axis_marginal.py:455` handles only `CellLayout`.
3. No gate pins the LD axis as `Declared`, although `spaces.rst:1249` and `axis.py:47` and `:268` say it is. The new origin needs its own first red.
4. No `automodule` covers `moment_layout`, `basis` or `coefficient_layout`, so their doc roles render as plain text, and grep is the only retirement gate.
5. Pins on the name or the path: `test_concept_order.py:43`, the test `..._policy_lives_in_the_layout_module`, `.claude/hooks/session-health.sh:36`, and `cartesian_multid.rst:2099` ("it did not move").
6. The bit-identity of the weights when the Gram is derived rather than typed is not measured. 1/(2k+1) computed as an exact rational and rounded gives the same float as the literal 1/3 `[R]`; the measurement is owed.

## Rulings (the user, 2026-10-25)

1. **The primitives are an interval Legendre basis and a tensor product of bases.** The cell basis is their composition. The face trace (one factor dropped, evaluated at an endpoint), the per-axis parity and the counts `p**d` and `p**(d−1)` follow from the factorisation.
2. **The moment axis is the cell basis's continuum head (`FromBasis`).** Ruled first as `FromFrame` (frame-dressed), then reversed the same day by the user once the frame machinery was read:
   - The frame face `frame_faces.basis_space` produces COVARIANT moments: the analysis outputs φ = Gc, with metric G⁻¹ (Parseval) and unit Φᵀw.
   - The LD spatial moments are CONTRAVARIANT expansion coefficients (ψ ≈ ψ̄ + ξψ̂), with metric G = `[1, 1/3]` and unit e₀.
   - Through `basis_space` the slope-row metric would invert, to `[1, 3]`. The test-architect found the same, independently.
   - So the axis is `FunctionSpace.of_basis(cell_basis)`, origin `FromBasis`: weights = the basis's continuum Gram, unit e₀.
   - The curvilinear case (#612) keeps its refusal guard, tagged with the machinery that retires it: a frame face that does not exist yet, the trial coefficient space dressed by the cell Gram (the ♯-image of `basis_space`).
3. **θ is read from the Gram.** The scheme constant retires, so there is one definition.
4. **The tail rule goes to numerics space composition, beside `FunctionSpace.of_axes`. The rank predicates go to the transport field layer.** `moment_layout.py` then retires.

## Steps (draft, for the user to rule polished)

Mode: W3 surgical carve. The carve crosses numerics, transport and sn, so the main agent writes the code. test-architect specifies the gates first; qa and elegance-enforcer review in parallel; the archivist writes the docs.

- **S0. The verification spec** (test-architect, `scratch/layer_dag/moment_basis/gate_spec.md`).
  - Interval basis:
    - the Gram `1/(2k+1)` against an independent rule (SymPy or mpmath);
    - `P_k(±1) = (±1)^k`;
    - the parity under ξ ↦ −ξ.
  - Tensor product:
    - its dimension;
    - its Kronecker order against `spatial-moment-kronecker-order`;
    - the face trace as one factor dropped.
  - The Cartesian frame's Gram equals the constant Legendre Gram.
  - The axis origin's first red, since no gate pins it today.
  - A bit-identity measurement of θ and of the LD solutions: which pins may move, and what the traced-memo regeneration costs.
  - For every re-spelled site, a mutation that reddens its gate.
- **S1. `numerics/basis/`.** Add the interval Legendre basis and the tensor product of bases. A coefficient layout joins the `CoefficientLayout` union and `FromBasis.kind`'s match, and `axis_marginal.py:455` is checked.
- **S2. Withdrawn by ruling 2's reversal.** No cell frame on this carve. The curvilinear refusal guard is re-tagged with its machinery: the frame's contravariant face (#612).
- **S3. The scheme reads the cell basis** (the continuum head, ruling 2 as reversed; no frame).
  - It mints its moment axis as the cell basis's continuum head (`Axis.of_basis`).
  - θ is read from the Gram.
  - `spatial_basis_per_axis` becomes the basis's degree.
  - The four hand spellings re-point: `scheme.py:518`, the float root at `fields/_bases.py:533`, the Kronecker decode at `_ubld.py:144`, and the hard-coded `2` in the LD statics.
- **S4. Retire `moment_layout.py`.**
  - The tail rule moves to numerics space composition, and the rank predicates to the transport field layer.
  - The `retirement-audit` A–G pass, with the dependency table at `.claude/plans/moment_basis_dependency_audit.md`, leaves first.
  - The path pins: `test_concept_order.py:43`, the test `..._policy_lives_in_the_layout_module`, `.claude/hooks/session-health.sh:36`.
- **S5. Docs** (archivist):
  - `cartesian_multid.rst`: §`spatial-moment-space`, its five labels, the "it did not move" line, and the stale module count near line 2306;
  - `spaces.rst:1249`, and `axis.py:47` and `:268`;
  - `verification/sn.rst` and `layering.rst`;
  - the History sections.
  #423 closes with the carve.
- **Bar before the merge:**
  - the full suite `-m "not slow"`, run per directory, detached;
  - the pyright ratchet;
  - Sphinx `-W`;
  - `dead_references` with a positive control;
  - CI watched after the push.

## ⏸ Progress (2026-10-25)

- **S0** done: `scratch/layer_dag/moment_basis/gate_spec.md` (the test-architect). Its findings that changed the plan:
  - `[M]` The continuum head is bit-identical to today's axis on weights and unit for d = 1, 2, 3 (8 Gram spellings, both Kronecker orders); a Gauss–Legendre-built Gram is 5–22 ulp off, a second reason for no frame.
  - `[M]` No existing pin moves under the identity change or under an S3 prototype (0 of 1262 rows in the 43-file LD scope).
  - Hazard 1 was the wrong mechanism: the traced memos go stale through module skeletons (`basis/__init__.py`, `space.py`), not LD identity, and `[M]` 0 of 8064 entries validate on the tree already, so the regeneration is owed independently of this carve.
  - Three rows demote at S3 (each compares the axis weights with `moment_mass_diagonal`, now one object); the value pin is the new `_HAND_MASS` row.
- **S1** done (uncommitted at this point): `IntervalLegendreBasis`, `TensorProductBasis`, `spatial_moment_basis`, `PolynomialLayout`, `TensorLayout`; `TensorProductBasis.face_dimension(axis)` added at S3 (total: 1 on a one-factor product, where `face` refuses).
- **S3** done (uncommitted at this point):
  - `DiscretizationSchemeBase.cell_basis(ndim)` (classmethod, minted per call so a Gram read is a read at use); `moment_axis` = `Axis.of_basis(cell_basis)` after the chart guard; the base `moment_mass_diagonal` reads the basis Gram and LD's override is deleted.
  - `LinearDiscontinuous.theta` is a property reading `IntervalLegendreBasis(2).coefficient_weights[1]` (`[M]` 2.9 µs per read).
  - Re-spelled: `scheme.py:518` (face count), `fields/_bases.py` (the float root becomes the layout's `factor_shape`, refusing a factor count other than the mesh dimension), `_ubld.octant_moment_frame_signs` (`reflection_signs`), `_FOUT_1D`/`_FIN_TRACE` (`endpoint_values`), `assemble_ubld`'s size (the mass's dimension), the LD statics' literal `2` (the basis's `factor_shape`, `face_dimension`, `endpoint_values(+1)`).
  - Gates landed: `tests/gates/numerics/test_interval_legendre_basis.py` (64), `tests/gates/transport/spatial/test_cell_basis_head.py` (16). The face-table row and the reflection row's leg 2 were re-pinned on the hand tables they replaced (kept in the test), so they did not demote. The three demoted rows carry docstrings naming the value pin; `g_bulk_measure` types θ = 1/3.
  - The curvilinear refusal names #612.
- **S4 done** (rulings of the user, 2026-10-25: the tail rule is `factor_tail`; the counts are read off the basis; `AVERAGE_MOMENT` is a constant in the basis module). The dependency audit is `.claude/plans/moment_basis_dependency_audit.md`.
  - `moment_layout.py` is deleted.
  - `factor_tail` is in `orpheus/numerics/space.py`, and `TensorProductBasis.dimension` is added.
  - The rank predicates are in `transport/fields/_bases.py`.
  - Two more hand spellings of the tail rule were re-pointed (`assembly.py`, `scheme.py`).
  - The retirement gates landed (`test_moment_layout_retired.py`, the stale-rank row in `test_concept_order.py`).
  - The pyright canary of `session-health.sh` is now `coefficient_layout.py`.
- **Review** (qa: every re-spelled site bit-identical, `[M]` 1200 of 1200 face-trace comparisons; elegance: two blockers). Fixed in the carve:
  - B1, the cell mass spelled twice. The UBLD kernels (`assemble_ubld`, `assemble_inflow_axis`, its transpose, `mass_1d`) take the cell basis instead of θ: the mass is `h·diag(Gram)` and the face traces are `endpoint_values`. `_FOUT_1D` and `_FIN_TRACE` retire. `_GRAD_1D` stays hand-typed, pinned by integration by parts.
    - Refuted FOR the gate design: the reviewer's decoy `[0.5, 0.25]`. It changes the degree-0 Gram, which is 1 for every interval Legendre factor, so it is outside the algebraic class.
    - The FACT the review established: the two spellings agreed only through that theorem. They are now one spelling, and the route row requires the kernel mass to move with the decoy.
  - B2: `moment_mass_diagonal` retired (0 production callers). The axis weights are the one spelling of the mass.
  - Scheme reads: `moment_tail(ndim)` and `face_moment_count(ndim)`, the one site stating that every face of the cell basis has one dimension.
  - `face`/`face_dimension` refuse an out-of-range axis, and the factors must be interval bases.
  - Leg 1 of the reflection row runs over n = 2, 3, 4: the parity rule and the old "flip every nonzero degree" rule differ from n = 3.
  - The width read matches `FromBasis(TensorProductBasis)` and carries `ELEGANCE-DEBT[guard]` #662 (filed: thread the cell basis, not the integer, through `spatial_moments=`).
  - Left as is: the `DiscretizationScheme` Protocol does not declare `cell_basis`. It is `runtime_checkable`, so a new member would change what `isinstance` admits for the duck-typed stubs. Consumers are typed against `DiscretizationSchemeBase`, and the pyright ratchet passes.
- **S5 done** (archivist): a new section `spaces-cell-basis` in `docs/theory/foundations/spaces.rst`, and every `moment_layout` and `moment_mass_diagonal` reference re-pointed across 7 pages. `[M]` `-E -W` clean; `dead_references` 12 → 3 (the 3 left predate the carve). The archivist's NEEDS: a `derive_interval_legendre_gram()` in `derivations/`, and the merge hash for the two history rows.
- `[M]` The scope after the review fixes: 1972 passed, 0 failed. The pyright ratchet passes.
