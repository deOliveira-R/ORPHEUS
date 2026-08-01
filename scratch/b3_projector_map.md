# B3 — the inflow/outflow projector surface: a map

**Scope**: campaign phase B3 of `.claude/plans/boundary_machinery_review.md` (§4.2, §5, B3 block).
**Tree**: `/Users/rodrigo/git/nuclear/ORPHEUS`, branch `refactor/operator-strategy-layers`,
HEAD `acb46245` (`refactor(sn,diffusion): boundary B2.2 — delete the stringly-typed dispatch`).
**Evidence tags**: `[R]` read-in-code · `[M]` measured-by-running · `[G]` grep/graph absence.
Scripts under `/Users/rodrigo/.claude/jobs/c30e4f25/tmp/`; all runs `PYTHONPATH=<repo> .venv/bin/python`.

> Line numbers current at HEAD `acb46245`. Tree is under concurrent edit
> (`docs/theory/verification/matrix.rst`, `.claude/skills/vv-principles/*` modified);
> nothing under `orpheus/` or `tests/` was dirty during this audit.

---

## HEADLINE (the durable claims)

1. **The plan's §5 measurement REPRODUCES — with one material sharpening.**
   `_reflect_trace`'s slice-write and `IncomingSourceOperator._inflow_mask` ARE the same
   map `P_in` (bit-identical). `IncomingOrdinateMaskTensor` IS `I − P_in`.
   **But `I − P_in ≠ P_out` in general** — it equals `P_out` only when the quadrature
   has NO tangential ordinates. On the **production cylinder mesh** (`CYL` +
   `Quadrature.product(n_mu=2, n_phi=4)`) 4 of 8 ordinates are tangential, and
   `IncomingOrdinateMaskTensor` preserves 12 of 24 dofs that `P_out` kills.
   So the honest statement is **two spellings of `P_in`, one spelling of `I − P_in`,
   and `P_out` is a FOURTH map that no operator type spells at all.**

2. **`[G]` No restriction/injection/selection operator exists anywhere in
   `orpheus/numerics/`.** There is no `RestrictionOperator`, no `SubspaceInjection`,
   no `SelectionOperator`, no index-set-parameterised projector. `IncomingOrdinateMaskTensor`
   is the closest thing and it is the wrong sign, wrong name, and (with tangential rows)
   the wrong set. Minting the projector is genuinely new construction.

3. **B3's bit-identity crux HOLDS for reflective and vacuum — the two SN-BC-reachable
   laws — and FAILS for albedo(α≠0) and periodic.** `[M]` For reflective the inflow-row
   output reads only outflow columns exactly (`perm[inflow] ⊆ outflow`, measured on slab /
   sphere / cylinder); for vacuum the inflow-row output is identically zero. For
   `AlbedoBoundary(α)` → `α·I` and `PeriodicBoundary` → identity-copy the inflow-row output
   reads the **inflow** rows, so narrowing the domain zeroes them (`|D|∞ = α`). That is not a
   regression B3 introduces — it EXPOSES that those two realizations were never
   `Γ₊ → Γ₋` maps.

4. **The trace partition is DISJOINT but NOT EXHAUSTIVE**, and inflow/outflow are
   **NOT contiguous blocks** except in 1-D Gauss-Legendre. A narrowed domain must be a
   **gather**, not a view/slice.

---

## Q1 — Every spelling of the inflow/outflow projector

### Q1.0 The §5 claim, re-measured independently `[M]`

Script `q1_projectors.py`, slab S4 + `ng=6` (so the face slot is `(4,6)` → the plan's 24×24):

```
face xmin slot shape (4, 6) N 4 ng 6
mu_x           = [-0.861136 -0.339981  0.339981  0.861136]
omega_dot_n[f] = [ 0.861136  0.339981 -0.339981 -0.861136]
inflow  indices = [2, 3]
outflow indices = [0, 1]

rank(P_in from slice-write)          = 12
rank(P_out from slice-write)         = 12
rank(IncomingOrdinateMaskTensor)     = 12
rank(IncomingSourceOperator mask)    = 12

np.array_equal(A, C)          (slice-write P_in  ==  IncomingSourceOperator mask) : True
np.array_equal(B, I - A)      (IncomingOrdinateMaskTensor  ==  I - P_in)          : True
np.array_equal(B, AO)         (IncomingOrdinateMaskTensor  ==  P_out)             : True
np.array_equal(A @ B, 0)      (P_in . (I-P_in) == 0)                              : True
np.array_equal(A + AO, I)     (P_in + P_out == I   <=> no tangential ordinate)    : True
all diagonal?  A: True  B: True  C: True  AO: True
idempotent?    A: True  B: True  AO: True
symmetric?     A: True  B: True  AO: True
```

**Verdict: the plan's §5 is CORRECT as stated for a slab.** All three are diagonal,
idempotent, symmetric, derived from the one index set.

### Q1.1 ⭐ The sharpening — `I − P_in == P_out` is NOT an identity `[M]`

Script `q1_cyl.py`, real `CYL` SNMesh, `Quadrature.product(n_mu=2, n_phi=4)`, `ng=3`:

```
CYL face xmax slot (8, 3) N 8 in [2, 6] out [0, 4]
rank P_in = 6  rank P_out = 6  rank IncomingOrdinateMaskTensor = 18
array_equal(A, C)        (slice-write P_in == IncomingSourceOperator mask) : True
array_equal(B, I - A)    (mask tensor == I - P_in)                          : True
array_equal(B, AO)       (mask tensor == P_out)                             : False
array_equal(A + AO, I)   (P_in + P_out == I)                                : False
rank(I - A - AO) [= #tangential dofs] : 12  of n_flat 24
```

`IncomingOrdinateMaskTensor` has **rank 18**, `P_out` has **rank 6**. The docstring's
claim (`orpheus/numerics/operator.py:2319-2321`) that it is "projection onto the outflow
subspace" is **FALSE on any mesh with tangential ordinates** — it is projection onto
`outflow ⊕ tangential`. `[R]` The implementation is honest (`out[self.inflow_indices] = 0`,
`operator.py:2374`); only the prose is wrong.

### Q1.2 Full census — production

Legend for **applied to**: **IN** = masks the operator's input · **OUT** = masks/injects the
operator's output · **SELF** = the mask IS the operator's own body · **FIELD** = writes a
restricted region of a stored field · **SET** = pure index-set arithmetic, no array touched.

| # | site (`file:line`) | map | applied to | note |
|---|---|---|---|---|
| P1 | `orpheus/sn/operators/boundary.py:359` — `out_boundary.face_view(face)[sel] = full[sel]` | `P_in` (or `P_rows ⊆ P_in` when `rows` given) | **OUT** | ⭐ the B3 site. Codomain projection `B_face = P_in ∘ law` |
| P2 | `orpheus/sn/operators/boundary.py:374-375` — `masked = np.zeros_like(face_in); masked[sel] = face_in[sel]` | `P_in` | **IN** | the transpose leg `(P_in ∘ law)ᵀ = lawᵀ ∘ P_in`. **NOTE: `P_in`, not `P_out`** |
| P3 | `orpheus/sn/operators/boundary.py:353-356` — `sel = rows[face] if rows is not None else trace.inflow_indices_for_face(face)` | index-set for `P_in` | **SET** | the single `sel` both P1 and P2 read |
| P4 | `orpheus/sn/operators/boundary.py:497-500` — `boundary_flux.face_view(face)[inflow] = reflected.face_view(face)[inflow]` | `P_in` | **FIELD** | `reflect_inflow_inplace`; ASSIGNMENT semantics |
| P5 | `orpheus/sn/operators/boundary.py:524-532` — `np.setdiff1d(trace.inflow_indices_for_face(face), lower_rows...)` | partition WITHIN `P_in` | **SET** | `split()` → `B_lower`/`B_upper` |
| P6 | `orpheus/sn/operators/boundary.py:906-908` — `boundary_flux.face_view(face)[rows] += reflected.face_view(face)[rows]` | `P_rows ⊆ P_in` | **FIELD** | `reflect_rows_inplace`; **ADDITIVE** (deliberately ≠ P4) |
| P7 | `orpheus/numerics/operator.py:2369-2379` — `out = x.copy(); out[self.inflow_indices] = 0.0` | `I − P_in` | **SELF** | ⭐ the vacuum law's entire body. NOT `P_out` (Q1.1) |
| P8 | `orpheus/sn/boundary/angular.py:331-333` (build) + `:358-360` (multiply `q * mask.reshape(...)`) | `P_in` dense diagonal | **SELF/OUT** | masks the affine source `q`, not any input |
| P9 | `orpheus/sn/boundary/angular.py:194-197` — `outgoing_mask = (outward_sign*mu_n) > 0.0; cos_w = np.where(outgoing_mask, cos_w, 0.0)` | a **weighted** `P_out'` | **SELF** | ⭐ a FOURTH predicate: strict `> 0.0`, NOT `> TANGENTIAL_EPS`. See Q5.3 |
| P10 | `orpheus/sn/boundary/realizer.py:239-240` — `outgoing = omega_dot_n[claimed] > +TANGENTIAL_EPS` | asserts `claimed ∩ outflow = ∅` | **SET** | ERR-041 orientation guard |
| P11 | `orpheus/sn/mesh/method_space.py:151` — `inflow_indices = trace.inflow_indices_for_face(face)` | index-set for `P_in` | **SET** | the realizer's plumbing source |
| P12 | `orpheus/sn/operators/streaming.py:975-978` and `:1098-1101` — `face_view(f)[out_rows] -= seed[...][out_rows]` | `P_out` | **OUT** | `L`'s outflow-definition face residual (fwd + transpose) |
| P13 | `orpheus/sn/loss_representation/__init__.py:1160-1168` — `[out_idx] = streamed−given`, `[in_idx] = given` | `P_out` **and** `P_in` | **OUT** | the 1-D `loss_action` boundary block. The only site that spells BOTH halves in one body |
| P14 | `orpheus/sn/loss_representation/__init__.py:1269-1276` | `P_out`/`P_in` | **OUT** | `loss_action_transpose`, the mirror |
| P15 | `orpheus/sn/loss_representation/__init__.py:3351-3369`, `:3520-3530` | `P_out`/`P_in` index sets | **SET/FIELD** | the 1-D scan's outer/inner face rows |
| P16 | `orpheus/sn/loss_representation/sweep_schedule.py:212-214` — `np.intersect1d(trace.inflow_indices_for_face(face), fresh_ordinates)` | sub-partition of `P_in` | **SET** | `lower_inflow_rows` |
| P17 | `orpheus/transport/source_sinks/angular_boundary_source_sink.py:272-285` — `view[inflow] = arr[inflow]` | `P_in` | **FIELD** | `prescribed_inflow` factory |
| P18 | `orpheus/sn/mesh/augmented_mesh.py:760` (docstring ref to `outflow_indices_for_face`) | — | doc | |
| **D1** | `orpheus/diffusion/operators.py:612-613` — `out_boundary.face_view(face)[INFLOW_ROW] = law.apply(trace.outflow_view(face))` | `P_out` on **IN** + `P_in` injection on **OUT** | both | ⭐ **the honest two-stage shape B3 copies** |
| **D2** | `orpheus/diffusion/operators.py:647` (`assemble`) | same pair, structural | both | the assembled twin of D1 |
| **D3** | `orpheus/transport/fields/scalar_boundary_flux.py:110-125` — `outflow_view` / `inflow_view` = `face_view(face)[ScalarTraceSpace.OUTFLOW_ROW/INFLOW_ROW]` | named `P_out` / `P_in` as **VIEWS** | **FIELD** | ⭐ the diffusion precedent for a NAMED restriction. `[R]` These are *views* because the scalar trace has exactly 2 rows — a constant index, not a computed set |

**The asymmetry that IS the defect (§4.2, re-confirmed `[R]`):** SN's P1 applies `P_in` to
the law's **OUTPUT** and never restricts the **INPUT**; diffusion's D1 applies `P_out` to
the **INPUT** and injects into the inflow row of the output. `[G]` **There is no SN
counterpart to `outflow_view`** — no `AngularBoundaryFlux.outflow_view`/`inflow_view` exists.
Grep of `orpheus/transport/fields/angular_boundary_flux.py` and `_bases.py` finds only the
generic `face_view` (`_bases.py:873`) and `face_views` (`:891`). The SN field surface stops
at the whole face slot.

### Q1.3 Full census — tests

| site | map | applied to | what it pins |
|---|---|---|---|
| `tests/sn/operators/test_sn_boundary_operator.py:150-157` | `P_in` on the law output + `P_out`-is-zero assert | OUT | ⭐ **C-1 re-pose target #1**: `expected_full = bc.apply(psi.boundary.face_view(face))` passes the **WHOLE face slot** to the law |
| `tests/sn/operators/test_sn_boundary_operator.py:208-217` | `P_in` on the law INPUT | IN | ⭐ **C-1 re-pose target #2**: pins `Bᵀ = lawᵀ ∘ P_inflow` explicitly |
| `tests/geometry/test_boundary.py:97-115` (`test_vacuum_bc_realizer_zeros_only_inflow_per_section_16A5`) | `I − P_in` | SELF | ⭐ **C-1 re-pose target #3**: pins the OUTFLOW PASS-THROUGH that B3 removes from the picture |
| `tests/geometry/test_boundary.py:462+` (`test_albedo_zero_and_vacuum_agree_on_inflow_rows`) | `P_in` | OUT | agreement on the inflow rows only |
| `tests/geometry/test_bc_equivalence_snapshot.py:141,152-160` | `I − P_in` | SELF | ⭐ **C-1 re-pose target #4**: Lebedev-17 snapshot; "outflow rows pass through unchanged" |
| `tests/sn/operators/test_snmesh_realizer_wiring.py:196,309,337` | `I − P_in` | SELF | "zeroes ONLY the inflow rows, leaving outflow rows untouched" |
| `tests/sn/operators/test_bc_extraction_2d.py:448-470` | `P_in` | OUT | ⭐ **the 2-D balance gate**: `coupled = B.apply(psi)`, compares `psi.boundary.face_view(face)[in_idx]` to `coupled...[in_idx]` |
| `tests/sn/operators/test_bc_extraction_2d.py:536,584-585` | `P_in`, `P_out` | FIELD | the inflow/outflow perturbation controls on `L` |
| `tests/numerics/test_incoming_ordinate_mask_tensor.py:29,37,53,62,73,83,95,105,128-165` | `I − P_in` | SELF | the mask's own law suite (10 tests) — incl. `test_leaves_non_inflow_untouched:37` |
| `tests/numerics/test_angular_trace_space.py:147-266` | index sets | SET | ⭐ `:190-191`/`:222-223` pin disjointness explicitly |
| `tests/sn/solve/test_gauss_seidel_reification.py:106,185,191,256,358` | `P_out`/`P_in` + `setdiff1d` | FIELD/SET | the W2 split |
| `tests/sn/operators/test_sweep_inverse_identity.py:169-170,213,253` | both | SET | `:253` names `"outflow_indices_for_face"` as a STRING |
| `tests/sn/operators/test_bc_extraction_matvec.py:649,715,720,775,784,787` | `P_in`/`P_out` | FIELD | 1-D matvec |
| `tests/sn/solve/test_b1pp_verification.py:276-277`, `tests/sn/sweep/core/test_phase_c_gates.py:793`, `tests/sn/verification/analytical/test_prescribed_inflow_consistency.py:97,202`, `tests/sn/verification/analytical/test_mms_prescribed_inflow.py:189,292`, `tests/sn/verification/mms/test_mms_ld_2d.py:1407`, `tests/transport/fields/test_angular_boundary_source_sink_residual.py:264,279` | `P_in`/`P_out` | FIELD | seeding / checking inflow slots |
| `tests/geometry/test_bc_universal_invariants.py:721-740` | `P_in` mask guards | SELF | `IncomingSourceOperator` construction guards |
| `tests/sn/primitives/test_method_space.py:114-146` | index set | SET | the delegation |
| `tests/diffusion/test_augmented_mesh.py:66-142` | — | — | `dm.bc[face].apply(_PROBE)` — probes the **already-narrowed** diffusion law |

---

## Q2 — Does a restriction/injection primitive already exist?

**`[G]` NO. Nothing in `orpheus/numerics/` restricts to or injects from an index-set
subspace.** Full inventory of `orpheus/numerics/` operator classes (`grep "class .*Operator"`,
all modules):

`SparseAssembledOperator` · `CoupledOperator` · `CoupledSubstitutionOperator` ·
`FlattenedOperator` · `_FrameAnalysis` · `_FrameReconstruction` · `GreenOperator` ·
`MatrixInverseOperator` · `DiagonalOperator` · `IdentityOperator` ·
**`IncomingOrdinateMaskTensor`** · `InverseOperator` · `PeriodicWrapOperator` ·
`PermutationOperator` · `RankOneOperator` · `ScaledOperator` · `SumOfTensorProductsOperator` ·
`TensorProductOperator` · `ZeroOperator` · `_AdjointOperator` · `OperatorSum` ·
`OperatorProduct` · ABCs `AnalysisOperator` / `ReconstructionOperator`.

Assessed candidates:

- **`IncomingOrdinateMaskTensor`** (`orpheus/numerics/operator.py:2284`). The only
  index-set-parameterised operator in the tree. API:
  `__init__(inflow_indices: np.ndarray, n_ordinates: int, axis: int = 0)`;
  `apply(x) -> x.copy() with x[idx] = 0`; `apply_transpose = apply`;
  `is_adjointable = True`; **no `inverse()`**, **no `as_matrix`**, **no `domain`/`codomain`**
  (inherits `None` from the base). It is `I − P_in`, i.e. the wrong sign for `P_in` and —
  with tangential rows — not `P_out` either. It **zeroes in place on the same space**; it
  does NOT change space. Not a restriction operator.
- **`AnalysisOperator` / `ReconstructionOperator`** (`orpheus/numerics/projection.py:88,136`).
  ABCs only — `[R]` "The ABC carries no implementation — the concrete realisation is the
  `analysis` face of a frame". Their only concretes are `_FrameAnalysis`/`_FrameReconstruction`
  (`orpheus/numerics/frame.py:412,446`), which are **basis** projections (fine→coarse via a
  `Basis`), not index selections. A `Γ → Γ₊` restriction is an *indicator* selection, and
  `orpheus/numerics/basis/indicator_basis.py` exists — but it is an energy/spatial group
  basis, not a trace-row selector. Using a frame here would be an abstraction mismatch
  (a frame carries an analysis/reconstruction PAIR + a discipline; a restriction has no
  meaningful reconstruction weight).
- **`PermutationOperator`** (`:2171`) — square, bijective, index-based. Cannot express a
  non-square restriction.
- **`orpheus/numerics/registry.py`** — despite `grep` hits for `VacuumBoundaryOperator` /
  `SpecularBoundaryOperator`, those names appear **only inside the module docstring** as a
  usage example. The module defines exactly one symbol: `RegistryMixin`.
- **`orpheus/numerics/measure.py:628`** has a `restrict` (indicator-multiplication on a
  `DiscreteMeasure`) — a *measure* restriction, not an operator, and not on a trace.

**Does `OperatorProduct` exist and would it accept a projector?** `[R]` **Yes** —
`orpheus/numerics/operator.py:1504`. Composition is `A @ B`, `apply(x) = a.apply(b.apply(x))`
(`:1533`). Its **only** construction-time guard is
`A.domain == B.codomain` (`:1560-1565`, raising `IncompatibleOperatorComposition`), and
`getattr(op, "domain", None)` is used throughout — so an operator advertising
`domain = codomain = None` (which every realized boundary law does today) composes freely.
`domain`/`codomain` of the product are read off `b.domain` / `a.codomain` (`:1573,1578`).
A new `Γ → Γ₊` restriction advertising a genuine narrowed `codomain` **would** be caught by
that guard if paired with a law whose `domain` is still the full trace — which is exactly
the type-check B3 wants, and exactly why the laws' `domain`/`codomain` being `None` today
means the guard is currently vacuous. `[G]` `OperatorProduct` is used by **no boundary
operator today** (the plan's B4 `[G]` re-confirmed: only `TensorProductOperator` `&` is used
in `orpheus/sn/boundary/realizer.py`).

**Bottom line `[G]`: the projector B3 needs has no name, no type, and no home. Minting it
is new construction, not a repoint.**

---

## Q3 — Every caller of a realized boundary law's `.apply` / `.apply_transpose`

### Q3.1 Graph leg

`mcp__nexus__callers(SNBoundaryOperator._reflect_trace)` → **exactly 2**, both internal:
`reflect_into_inflow` (`boundary.py:432`) and `_apply_faces` (`boundary.py:379`).
`mcp__nexus__context(IncomingOrdinateMaskTensor)` → the only production `calls` edge in is
`SNBoundaryRealizer.realize`; everything else is `tests/numerics/test_incoming_ordinate_mask_tensor.py`
plus 3 doc pages.

### Q3.2 Text leg — `grep -rn "law\.apply"` across `orpheus/` + `tests/` (worktree excluded)

**Production, real call sites — exactly THREE:**

```
orpheus/diffusion/operators.py:613:    law.apply(trace.outflow_view(face))
orpheus/sn/operators/boundary.py:358:  full = law.apply(face_in)
orpheus/sn/operators/boundary.py:376:  out_boundary.face_view(face)[...] = law.apply_transpose(masked)
```

(the remaining hits — `orpheus/geometry/boundary/_base.py:26,121`, `__init__.py:198` — are
docstrings stating that a *descriptor* `law.apply` does NOT exist.)

**`[R]` The plan's §4.1 claim is CONFIRMED and now precise:** the sole SN caller of a realized
law's forward `.apply` is `_reflect_trace` (`boundary.py:358`, was `:301` when the plan was
written — drift), and the sole caller of `.apply_transpose` is `boundary.py:376`.

| caller | passes IN | does WHAT with the output |
|---|---|---|
| `boundary.py:358` (fwd) | **the WHOLE face slot** `boundary.face_view(face)` | writes ONLY `[sel]` (inflow rows) into a zero sink; discards the rest |
| `boundary.py:376` (transpose) | the face slot **masked to `P_in`** (`masked[sel] = face_in[sel]`) | writes the **FULL** `lawᵀ` image (`[...] =`) |
| `diffusion/operators.py:613` | **the outflow half only** `trace.outflow_view(face)` | writes the whole result into `[INFLOW_ROW]` |

Note the SN transpose leg already masks its **INPUT** — but to `P_in`, not `P_out`. Under a
`Γ₊ → Γ₋` law the transpose becomes `Γ₋ → Γ₊`, so this leg is re-posed too, and the
`⚠` comment at `boundary.py:301-308` (the ERR "output-projecting `lawᵀ` onto the OUTFLOW
rows extracts a law's DIAGONAL block") is exactly the hazard B3 must not re-introduce.

### Q3.3 Indirect leg — who reaches `_reflect_trace`

`[R]` Complete chain (all in `orpheus/sn/operators/boundary.py`):

```
SNBoundaryOperator.apply            :430  -> _apply_faces(psi, "apply")       :425 -> _reflect_trace
SNBoundaryOperator.apply_transpose  :546  -> _apply_faces(psi, "apply_transpose")
SNBoundaryOperator.reflect_into_inflow :458 -> _reflect_trace(boundary, "apply", faces=)
SNBoundaryOperator.reflect_inflow_inplace :491 -> reflect_into_inflow
SNMaskedBoundaryOperator.apply      :866  -> inner._apply_faces(psi, "apply", rows=)
SNMaskedBoundaryOperator.reflect_rows_inplace :902 -> inner._reflect_trace(..., rows=)
```

`[R]` The realized law object reached is `sn_mesh.bc[face]`, a
`_BoundBoundaryOperator` wrapper (`orpheus/geometry/boundary/_bound_compat.py:95`) whose
`apply` is a pure pass-through (`:218-219`) and whose `apply_transpose` guards on
`adjointable(self.inner)` (`:221-229`). It has NO `domain`/`codomain` of its own —
`getattr(self.inner, "domain", None)` (`:171,175`) — which for every realized law is `None`.

### Q3.4 ⭐ Duck-typed / class-name-bypass surrogate leg (the FOURTH search)

`grep -rn "SimpleNamespace(\|class _NotALaw\|class _Fake\|class _Stub\|Mock("` across
`tests/`, filtered to the boundary layer. **Six surrogates the graph and a symbol grep
both miss:**

| surrogate | `file:line` | what it stands in for | does it call `law.apply`? | B3 exposure |
|---|---|---|---|---|
| `_NoTransposeLaw` + `_BWithStubFace(SNBoundaryOperator)` overriding `_face_laws` | `tests/sn/operators/test_sn_boundary_operator.py:237-253` | **a realized LAW** — `apply(x) -> x`, `is_adjointable=False` | **YES**, reached via the production `_reflect_trace` path | ⭐ **HIGHEST RISK.** A `_face_laws` **property override on a real subclass** — invisible to `callers()`. Its `apply` will receive whatever B3 hands the law. It currently only asserts the `MissingAdjoint` raise, so it should stay green, but any signature/shape change reaches it |
| `SimpleNamespace(law=…)` mesh (`bc` dict), `angular_trace=SimpleNamespace(layout=SimpleNamespace(faces=…))` | `tests/sn/sweep/core/test_sweep_schedule_nd.py:133-140` | a 3-axis **mesh** | no — reads `.law` factors only | ⭐ its `angular_trace` stub carries **`layout` ONLY** — no `inflow_indices_for_face`, no `omega_dot_n`. If B3 makes the schedule (or anything it calls) touch a new trace method, this dies with `AttributeError` |
| `SimpleNamespace(bc={face: SimpleNamespace(law=…)}, scheme=SimpleNamespace(key=…))` | `tests/sn/acceleration/test_dsa_low_order.py:164-171` | a **mesh** for DSA admission | no — `.law` only | low; the surrogate that bit B2 |
| `SimpleNamespace(is_1d=…, scheme=SimpleNamespace(...))` | `tests/sn/sweep/core/test_unified_sweep_dispatch.py:353-374` | a mesh for dispatch | no | none (no trace at all) |
| `SimpleNamespace(...)` ×2 | `tests/sn/solve/test_si_gate_dispatch.py:62,70` | solver gate objects | no | none |
| `class _NotALaw` | `tests/diffusion/test_boundary_realizer.py:331`, `tests/geometry/test_law_composition.py:348` | a non-law refusal probe | no | none |
| `class _StubLaw(BoundaryTraceLaw, key="_stub_for_test")` | `tests/geometry/test_boundary_trace_law.py:53` | a **descriptor** law | no (descriptors aren't callable) | low, but it registers into the PRODUCTION registry at import (documented at `:63-75`) |
| `class _FakeQuad` | `tests/geometry/test_bc_universal_invariants.py:117` | a quadrature | n/a | if B3 asks the trace/quadrature for a new datum |

**Operational consequence:** B3's gate must include a **full-suite run before the commit**,
not after — for exactly the reason B2 learned (plan lines 1011-1038). The
`_BWithStubFace._face_laws` override is a shape the graph cannot see at all.

### Q3.5 Doc leg (Operating Principle 4 — the fourth search)

`grep -rn "IncomingOrdinateMaskTensor|inflow_indices_for_face|outflow_indices_for_face|_reflect_trace|IncomingSourceOperator" docs/` — **6 pages, ~35 references.** The Python-domain roles
here are the SILENT class (no `-W` warning without `-n`):

- `docs/theory/foundations/boundary_conditions.rst` — the heaviest: `:105`, `:291`, `:296`,
  `:319`, `:353`, `:413-414`, **`:664-665`** ("`apply_transpose` writes the
  `outflow_indices_for_face` slots" — ⚠ this contradicts `boundary.py:376`, which writes the
  FULL image; a B0 survivor to check), **`:891`, `:908`, `:938-940`** (the `_reflect_trace`
  projection prose), `:1604-1606`, `:2064-2068`, `:2096-2097`, `:2260`, `:2365`, `:2371-2373`,
  `:2419-2426`, `:2541`, `:3096`, `:3202`, `:3225`
- `docs/theory/foundations/operator_algebra.rst:3882-3903`
- `docs/theory/foundations/operator_adjoint.rst:533`
- `docs/theory/foundations/operator_tensor_network.rst:177` (`IncomingOrdinateMaskTensor(...) & …`)
- `docs/theory/methods/sn/boundary_conditions.rst:305`, `:344`
- `docs/theory/methods/sn/cartesian_multid.rst:3931`
- `docs/theory/methods/sn/curvilinear_one_group.rst:2515`

---

## Q4 — The trace-space ordering contract

### Q4.1 Who decides `[R]`

`orpheus/numerics/spaces/angular_trace_space.py`:

- `inflow_indices_for_face` (`:432-439`): `np.flatnonzero(row < -TANGENTIAL_EPS)`
- `outflow_indices_for_face` (`:441-448`): `np.flatnonzero(row > +TANGENTIAL_EPS)`
- `row = self.omega_dot_n[self._face_row(face)]`, built by `build_omega_dot_n` (`:204-259`)
  as `sign_f * mu_axis(f)` from the face NAME alone.
- `TANGENTIAL_EPS = 4.0 * np.finfo(np.float64).eps ≈ 8.88e-16` (`:158`).

The module docstring (`:66-71`) states the three-way classification explicitly:
inflow `< -ε`, outflow `> +ε`, **tangential `|·| ≤ ε` — in NEITHER selector.**

### Q4.2 Is the partition exact? `[M]` — DISJOINT ALWAYS, EXHAUSTIVE ONLY SOMETIMES

Script `q4_tangential.py` / `q4_real_mesh.py`:

```
TANGENTIAL_EPS = 8.881784197001252e-16

gauss_legendre(4)       face=xmin N=  4 |in|=  2 |out|=  2 |tan|= 0 disjoint=True exhaustive=True  contig_in=True  contig_out=True
gauss_legendre(8)       face=xmin N=  8 |in|=  4 |out|=  4 |tan|= 0 disjoint=True exhaustive=True  contig_in=True  contig_out=True
gauss_legendre(5) odd   face=xmin N=  5 |in|=  2 |out|=  2 |tan|= 1 disjoint=True exhaustive=False contig_in=True  contig_out=True
      TANGENTIAL rows [2]  values=[-0.]
product(n_mu=2,n_phi=4) face=xmin N=  8 |in|=  2 |out|=  2 |tan|= 4 disjoint=True exhaustive=False contig_in=False contig_out=False
      TANGENTIAL rows [1, 3, 5, 7]  values=[-4.99959962e-17  1.49987989e-16]
level_symmetric(4)      face=xmin N= 24 |in|= 12 |out|= 12 |tan|= 0 disjoint=True exhaustive=True  contig_in=False contig_out=False
level_symmetric(6)      face=xmin N= 48 |in|= 24 |out|= 24 |tan|= 0 disjoint=True exhaustive=True  contig_in=False contig_out=False
level_symmetric(8)      face=xmin N= 80 |in|= 40 |out|= 40 |tan|= 0 disjoint=True exhaustive=True  contig_in=False contig_out=False
```

On **real production SNMeshes**:

```
--- SLB: faces=('xmin', 'xmax')  N=4  layout total=16
    face=xmin  slot shape=(4, 2)  in=[2, 3] out=[0, 1] TANGENTIAL=[]
    face=xmax  slot shape=(4, 2)  in=[0, 1] out=[2, 3] TANGENTIAL=[]
--- SPH: faces=('xmax',)  N=4  layout total=8
    face=xmax  slot shape=(4, 2)  in=[0, 1] out=[2, 3] TANGENTIAL=[]
--- CYL: faces=('xmax',)  N=8  layout total=16
    face=xmax  slot shape=(8, 2)  in=[2, 6] out=[0, 4] TANGENTIAL=[1, 3, 5, 7]
lebedev(3) N=   6 face=xmin: |in|=1 |out|=1 |tangential|=4  exhaustive=False
lebedev(5) N=  14 face=xmin: |in|=5 |out|=5 |tangential|=4  exhaustive=False
lebedev(7) N=  26 face=xmin: |in|=9 |out|=9 |tangential|=8  exhaustive=False
```

**Verdict:**

- **Disjoint: ALWAYS**, structurally (`< -ε` vs `> +ε` cannot both hold). Pinned by
  `tests/numerics/test_angular_trace_space.py:190-191` and `:222-223`.
- **Exhaustive: NO.** Tangential ordinates are **reachable on the production cylinder**
  (4/8 at `xmax`) and on every Lebedev set. `inflow ∪ outflow ⊊ all rows` in general.
- The tangential values on the cylinder are **NOT exactly zero** (`±4.9e-17`, `±1.5e-16`) —
  they are floating-point residue from the product-quadrature construction, sitting just
  under `TANGENTIAL_EPS`. The eps is load-bearing, not decorative.

### Q4.3 Contiguity `[M]` — NO, except 1-D Gauss-Legendre

- 1-D Gauss-Legendre: **contiguous** (`in=[2,3]`, `out=[0,1]`) because `mu_x` is sorted.
- 2-D/3-D level-symmetric (S4/S6/S8): **NOT contiguous** at any face, including `xmin`/`xmax`.
- Cylinder product: **NOT contiguous** (`in=[2,6]`, `out=[0,4]` — stride-4 interleave).

**⭐ Consequence for B3:** a narrowed `Γ₊` domain **cannot be a cheap view/slice**. It must be
a **gather** (`x[out_idx]`, a copy) plus a **scatter** on the way back
(`y[in_idx] = ...`). The diffusion precedent (`outflow_view` = `face_view(face)[OUTFLOW_ROW]`,
a genuine memory-sharing view) does **not** transfer: the scalar trace's row index is a
CONSTANT (`ScalarTraceSpace.INFLOW_ROW`/`OUTFLOW_ROW`), while the angular trace's is a
quadrature-dependent fancy index. Any SN `outflow_view` will return a COPY, and writes to it
will NOT propagate — a semantic difference from the diffusion sibling that must be spelled
out (or the operation must be posed as an explicit restriction OPERATOR rather than a `_view`
accessor whose diffusion namesake is a writable view).

---

## Q5 — The reflective law's read pattern (the bit-identity crux)

### Q5.1 The construction `[R]`

- `SNBoundaryRealizer.realize`, reflective arm (`orpheus/sn/boundary/realizer.py:270-293`):
  `base = PermutationOperator(quad.reflection_index(law.axis), axis=0) & IdentityOperator()`
  (a 2-factor `TensorProductOperator`), scaled by `α` when `α ≠ 1.0`.
- `PermutationOperator.apply` is `np.take(x, perm, axis=0)`, so
  `law.apply(x)[i] = x[perm[i]]` — **row `i` of the output reads row `perm[i]` of the input,
  and nothing else.**
- Vacuum arm (`:253-268`): `IncomingOrdinateMaskTensor(inflow_indices, n_ordinates, axis=0)
  & IdentityOperator()`. Output inflow rows are identically `0`.

### Q5.2 Measurement — the gate `P_in · law · (I − P_out) == 0`

That product is the exact statement "masking the input to the outflow half is a no-op for the
inflow-row output". Script `q5_read_pattern.py` (dense matrices over the whole face slot;
the end-to-end column re-applies the law to a random slot vs its outflow-masked copy and
compares the inflow rows):

```
=== SLAB S4 (no tangential)  face=xmin slot=(4, 3) in=[2, 3] out=[0, 1] tangential_dofs=0
  VacuumInflow           P_in.law.(I-P_out)==0 : True    |D|inf=0.000e+00   end-to-end inflow-rows bit-identical: True
  ReflectiveBoundary     P_in.law.(I-P_out)==0 : True    |D|inf=0.000e+00   end-to-end inflow-rows bit-identical: True
  Reflective(a=0.5)      P_in.law.(I-P_out)==0 : True    |D|inf=0.000e+00   end-to-end inflow-rows bit-identical: True
  WhiteBoundary(xmin)    P_in.law.(I-P_out)==0 : True    |D|inf=0.000e+00   end-to-end inflow-rows bit-identical: True
  AlbedoBoundary(0.0)    P_in.law.(I-P_out)==0 : True    |D|inf=0.000e+00   end-to-end inflow-rows bit-identical: True
  AlbedoBoundary(0.4)    P_in.law.(I-P_out)==0 : False   |D|inf=4.000e-01   end-to-end inflow-rows bit-identical: False
  AlbedoBoundary(1.0)    P_in.law.(I-P_out)==0 : False   |D|inf=1.000e+00   end-to-end inflow-rows bit-identical: False
  PeriodicBoundary       P_in.law.(I-P_out)==0 : False   |D|inf=1.000e+00   end-to-end inflow-rows bit-identical: False
  PrescribedInflow(0)    P_in.law.(I-P_out)==0 : True    |D|inf=0.000e+00   end-to-end inflow-rows bit-identical: True
  PrescribedInflow(q)    P_in.law.(I-P_out)==0 : False   |D|inf=2.500e+00   end-to-end inflow-rows bit-identical: True

=== CYL product(2,4) (4 tangential ordinates)  face=xmax slot=(8, 3) in=[2, 6] out=[0, 4] tangential_dofs=12
  VacuumInflow           P_in.law.(I-P_out)==0 : True    |D|inf=0.000e+00   end-to-end inflow-rows bit-identical: True
  ReflectiveBoundary     P_in.law.(I-P_out)==0 : True    |D|inf=0.000e+00   end-to-end inflow-rows bit-identical: True
  Reflective(a=0.5)      P_in.law.(I-P_out)==0 : True    |D|inf=0.000e+00   end-to-end inflow-rows bit-identical: True
  AlbedoBoundary(0.0)    P_in.law.(I-P_out)==0 : True    |D|inf=0.000e+00   end-to-end inflow-rows bit-identical: True
  AlbedoBoundary(0.4)    P_in.law.(I-P_out)==0 : False   |D|inf=4.000e-01   end-to-end inflow-rows bit-identical: False
  AlbedoBoundary(1.0)    P_in.law.(I-P_out)==0 : False   |D|inf=1.000e+00   end-to-end inflow-rows bit-identical: False
  PeriodicBoundary       P_in.law.(I-P_out)==0 : False   |D|inf=1.000e+00   end-to-end inflow-rows bit-identical: False
  PrescribedInflow(0)    P_in.law.(I-P_out)==0 : True    |D|inf=0.000e+00   end-to-end inflow-rows bit-identical: True
  PrescribedInflow(q)    P_in.law.(I-P_out)==0 : False   |D|inf=2.500e+00   end-to-end inflow-rows bit-identical: True
  WhiteBoundary(xmax)    P_in.law.(I-P_out)==0 : False   |D|inf=3.062e-17   end-to-end inflow-rows bit-identical: True
```

### Q5.3 Reading the table — five distinct verdicts

**(a) REFLECTIVE — the narrowing is an EXACT no-op `[M]`.** Confirmed both by the matrix
identity and directly on the permutation, on every geometry:

```
### (i) reflective permutation: does perm[inflow] land in outflow?
  SLB xmin: perm=[3, 2, 1, 0]  inflow=[2, 3] -> reads rows [1, 0]  outflow=[0, 1]  SUBSET_OF_OUTFLOW=True
  SLB xmax: perm=[3, 2, 1, 0]  inflow=[0, 1] -> reads rows [3, 2]  outflow=[2, 3]  SUBSET_OF_OUTFLOW=True
  SPH xmax: perm=[3, 2, 1, 0]  inflow=[0, 1] -> reads rows [3, 2]  outflow=[2, 3]  SUBSET_OF_OUTFLOW=True
  CYL xmax: perm=[2, 1, 0, 3, 6, 5, 4, 7]  inflow=[2, 6] -> reads rows [0, 4]  outflow=[0, 4]  SUBSET_OF_OUTFLOW=True
```

This is the ERR-045 "inflow→outflow table" invariant that `ReflectiveBoundary.assert_realizable`
already certifies at realization — B3 is *consuming* an invariant the law layer already
enforces, not adding a new assumption. Holds for `α = 0.5` too (scaling is diagonal).

**(b) VACUUM — no-op, trivially `[M]`.** The inflow-row output is identically zero; it reads
nothing. `|D|∞ = 0` on both meshes.

**(c) ALBEDO(α≠0) and PERIODIC — the narrowing CHANGES VALUES `[M]`.** `|D|∞ = α` / `1.0`.
`[R]` The cause: `AlbedoBoundary(α)` realizes to `α·(I & I)` (`realizer.py:312-326`) and
`PeriodicBoundary` to `PeriodicWrapOperator() & IdentityOperator()` (`:328-337`), whose body
is `np.asarray(x).copy()` (`operator.py:2421-2424`). Both are **ordinate-diagonal identity
maps**: output row `n` reads input row `n`. Since inflow rows are disjoint from outflow rows,
narrowing the domain to `Γ₊` makes their inflow-row output **identically zero**.

  This is not a B3 bug — it is B3 *revealing* that these two realizations are `Γ → Γ` maps
  where the physics is `Γ₊ → Γ₋`. A physical albedo at a slab face is
  `γ₋ψ(μ) = α·γ₊ψ(−μ)`, which needs the **permutation** `α·P`; `α·I` cannot be right on
  either domain. Same for periodic (which needs the partner-face pushforward that
  `PeriodicWrapOperator`'s docstring `operator.py:2402-2408` admits it does not implement).

  **Reachability `[R]`:** `orpheus/sn/mesh/augmented_mesh.py` admits only
  `{vacuum, reflective}` from a `BC` declaration (per the plan §6, pinned by
  `tests/sn/operators/test_snmesh_realizer_wiring.py:403`), so **neither is production-reachable
  on the SN side today**. But both are *realizable* by `SNBoundaryRealizer`, so a
  directly-constructed law reaches them, and `tests/geometry/test_boundary.py` /
  `test_bc_universal_invariants.py` exercise them at the law layer. B3 must decide
  explicitly: refuse them at realization, or fix them to carry the permutation (which is
  B4's `R ∘ G` work), or accept a documented value change on an unreachable path.

**(d) PRESCRIBED_INFLOW — the matrix column is MEANINGLESS; end-to-end it is a no-op `[M]`.**
`IncomingSourceOperator.apply` **ignores its input** (`angular.py:336-360`: `q =
self.source.evaluate(shape)`). It is therefore **not a linear operator** — `apply(0) = q ≠ 0` —
so the unit-probe densification produces a constant "matrix" and `|D|∞ = 2.5` is an artefact.
The honest end-to-end column shows `True`: masking the input cannot change an ignored input.
⚠ **Corollary worth carrying into B4:** any `as_matrix`-style probing extraction (as diffusion's
`assemble` does at `orpheus/diffusion/operators.py:626-637`) is **wrong for this operator**.

**(e) ⭐ WHITE — the narrowing is NOT strictly value-neutral, by a genuinely different
predicate `[M]`.** On the cylinder `|D|∞ = 3.06e-17` — small, but not round-off. The cause:

```
### (ii) white on CYL: AngularAverageOperator's outgoing predicate vs the trace's
  mu_x                     = [ 8.1649658092772603e-01  4.9995996217394874e-17 -8.1649658092772603e-01
 -1.4998798865218462e-16  8.1649658092772603e-01  4.9995996217394874e-17 -8.1649658092772603e-01 -1.4998798865218462e-16]
  operator cos_w           = [1.282550e+00 7.853353e-17 0.000000e+00 0.000000e+00 1.282550e+00 7.853353e-17 0.000000e+00 0.000000e+00]
  cos_w nonzero rows       = [0, 1, 4, 5]
  trace outflow rows       = [0, 4]
  trace tangential rows    = [1, 3, 5, 7]
  => operator reads rows the trace calls TANGENTIAL: [1, 5]
  probe: non-outflow rows = 1e12, outflow rows = 1.0
  full   inflow-row value = np.float64(1.0000612323399574)
  narrow inflow-row value = np.float64(1.0)
  bit-identical = False  abs diff = 6.123234e-05
```

`[R]` `AngularAverageOperator.from_quadrature` uses `outgoing_mask = (outward_sign * mu_n) > 0.0`
(`orpheus/sn/boundary/angular.py:194`) — a **strict-zero** test — while the trace uses
`> TANGENTIAL_EPS` (`angular_trace_space.py:448`). Two independent encodings of "outgoing",
and on the product quadrature they **disagree on rows 1 and 5**. The magnitude is
`~1e-16`-weighted, so it is invisible on physical data, but it is a **real twin-source-of-truth
defect** (exactly the ERR-041 shape the vacuum arm already has a cross-check for). White is
not SN-BC-reachable, so this does not bite today.

### Q5.4 Bottom line for the B3 gate

| law | narrowing is value-neutral? | SN-BC-reachable? |
|---|---|---|
| `ReflectiveBoundary(α)` | **YES, exactly** (all geometries measured) | **yes** |
| `VacuumInflow` | **YES, trivially** | **yes** |
| `AlbedoBoundary(0.0)` | YES (it is `ZeroOperator`) | no |
| `PrescribedInflow(q)` | YES (input ignored) — but the law is AFFINE, not linear | no |
| `WhiteBoundary` | **NO** — `~1e-16`-weighted tangential leak, from a second `outgoing` predicate | no |
| `AlbedoBoundary(α≠0)` | **NO** — `|D|∞ = α`; the realization is `α·I`, a `Γ→Γ` map | no |
| `PeriodicBoundary` | **NO** — `|D|∞ = 1`; the realization is identity-copy | no |

**The B3 bit-identity gate holds on every production-reachable path.** The three failures are
all on paths a `BC` declaration cannot reach, and each is a pre-existing typing defect that
the narrowing *surfaces* rather than causes.

---

## Gaps / things I did NOT find

- **`[G]` No SN `outflow_view` / `inflow_view`** on `AngularBoundaryFlux` or any
  `_bases.py` boundary family. The diffusion accessor pair has no angular sibling.
- **`[G]` No `domain` / `codomain` on any realized boundary law.** Every one inherits the
  base `None`, and `_BoundBoundaryOperator` forwards that `None`
  (`_bound_compat.py:171,175`). So `OperatorProduct`'s space guard is currently **vacuous**
  for boundary composition — B3 narrowing the domain is also the first opportunity to make
  that guard bite.
- **`[R]` A docstring/code discrepancy to check in B3**:
  `docs/theory/foundations/boundary_conditions.rst:664-665` says "`apply_transpose` writes the
  `outflow_indices_for_face` slots", but `boundary.py:376` writes the **full** `lawᵀ` image
  (`face_view(face)[...] = ...`) after masking the INPUT to inflow. The `⚠` block at
  `boundary.py:301-308` explicitly says output-projecting onto outflow is the WRONG spelling.
  The doc line looks like a pre-B.2d survivor.
- **`[R]` `IncomingOrdinateMaskTensor`'s docstring (`operator.py:2319-2321`) claims
  "projection onto the outflow subspace"** — measured false on any tangential-carrying
  quadrature (Q1.1). It is `I − P_in`.
- I did **not** run the pytest suite (read-only brief); the surrogate exposure in Q3.4 is a
  static assessment.
