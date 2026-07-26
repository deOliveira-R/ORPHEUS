# DSA (#2) landing-zone recon — Phase-3 facts re-verified against main @ `e4c1a81c`

**Verified 2026-07-26** (explorer dispatch). Baseline being re-verified: the facts
`.claude/plans/stencil_assembly_dsa_roadmap.md` recorded at `3a19133` (2026-07-03).
Landed since: Phase 2.5 / #280 walk unification (`3f0b8c74`), `sn/spatial` → `sn/sweep`
(task #54, `588f2429`), the #276/#310 adjoint campaign, the #281 P6 collapse campaign,
and the #231/#304 doc part-tree restructure. Nexus graph built at `e301b675` (code
identical to HEAD — the two commits since are plan/agent files only). Line numbers are
**at this HEAD**; re-derive via Nexus after any merge.

---

## Drift table — roadmap claim → current fact

| # | Roadmap claim (@ `3a19133`) | Current fact (@ `e4c1a81c`) |
|---|---|---|
| D1 | `numerics/operator.py:751` — `as_matrix` is probing-only; NO operator has stencil assembly | Phase 2b landed: `as_matrix` @ `orpheus/numerics/operator.py:868` **delegates to `assemble().as_matrix()`** when `assemblable(self)` (:961-967); probing retained as `_as_matrix_by_probing` (:970). Assembly axis: `is_assemblable` :637 (default False), `SupportsAssembly` :1072, `assemblable()` :1128, `MissingAssembly` :424; composer laws Sum :1478 / Product :1676 / Scaled :1812 |
| D2 | Stencil data in `orpheus/sn/spatial/{cell_balance,linear_discontinuous,_ubld}.py` | Promoted (2a): `orpheus/transport/spatial/{scheme,diamond,linear_discontinuous,_ubld,cell_balance}.py`. Residual package renamed (task #54): `orpheus/sn/sweep/` (`cache.py`, `scan.py`, `pole_angular_closure.py`, `psi_half_angle_seed.py`, `pairing.py`). Roadmap ruling R9 ("keep the name") is superseded-by-execution |
| D3 | Within-group factory = `_within_group_triple` (sn/solver.py) returning `(L+C, S, B)` | **RETIRED** → `build_within_group_system` @ `orpheus/sn/coupled_system.py:380` returning the `WithinGroupSystem` record (:283): `{loss, space, resolvent, gains}`; 1×1 seedless / 2×2 carrying grid; one LC spelling `build_streaming_collision` :356 |
| D4 | `evaluate_residual` returns the 2-block typed residual | Signature widened @ `orpheus/sn/solver.py:231`: `(loss_op, psi: TimedFullField \| CoupledField, q_ext) → FullField \| CoupledField` — a **coupled arm** (r_A ⊕ split-ψ½ r_B) on carrying meshes, with a **loud guard** (:324-329) refusing a bare System-A call on a carrying mesh (explicitly a DSA-blindness guard). Composite member spelling is now `FullField(interior=…, boundary=…)` — **not** `bulk=` |
| D5 | `ScatteringOperator` in the SN package | Moved: `orpheus/transport/operators/scattering.py:356`. `foldable_part` :978, `residual_part` :1016, `is_foldable_into_sigma_r` :1057 — **zero production callers** (Nexus-confirmed); the #215 σ_r-trap warning is inline @ :968-976 |
| D6 | `MaterialXSField.foldable_sigma()` | `orpheus/transport/mesh/material_xs_field.py:1213`; siblings `foldable_sig_s` :1148, `residual_sig_s` :1167, `is_p0_diagonal_with_zero_n2n` :1189 |
| D7 | `as_dsa_source` lands with #2 (does not exist) | Confirmed: **absent from `orpheus/` and `tests/`**; the "lands WITH DSA #2" doc note lives at `docs/theory/foundations/field_algebra.rst:528` (moved out of operator_algebra.rst) |
| D8 | KEigenvalue estimator-injection seam @ `iteration.py:1217,1233` (retirement candidate) | R8 executed: injection kwargs gone. `KEigenvalue` @ `orpheus/numerics/iteration.py:1000`, `__init__(A, S, F, *, …)` :1089 (requires `A.is_invertible`; inner = `SourceIteration(seeded_inverse(A), S)` :1143); plain methods `compute_production_rate` :1212 / `compute_keff` :1237 |
| D9 | Krylov preconditioner posture unknown / #200 | Seam **exists and is live**: `KrylovAcceleration(…, preconditioner=…)` @ `iteration.py:810-818`; default = sweep (`seeded_inverse(A).apply`, :846-849) when invertible; scipy left-M @ :912-916, :943-950. Production **overrides with explicit identity** @ `orpheus/sn/solver.py:487` (`preconditioner=lambda q: q` — "#200 tracks re-enable"); `restart=n_dof` full size :489 (ERR-053) |
| D10 | P7b deferral: `full_field_space` not on `TransportMethod`; declare with DSA | Verbatim intact @ `orpheus/transport/method.py:162-172`: "Declare it when the first generic consumer arrives (the DSA driver, #2)". Both witnesses implement it as `cached_property` (`SNMesh`, `DiffusionMesh` @ `diffusion/augmented_mesh.py:335`) |
| D11 | Diffusion resolvent materializes by **probing** (tolerable, small N) | Composition unchanged (`MatrixInverseOperator(FlattenedOperator(A, template))` @ `diffusion/solver.py:250-252`) but it now materializes **through structural assembly** — all four loss leaves emit (L :441, B :626 in `diffusion/operators.py`; C `multiplication_operator.py:258`; S `isotropic_scattering.py:305` + N2N :414), the Sum derives `is_assemblable`, and `as_matrix` delegates (D1). Probing survives only as the pinned oracle (`tests/diffusion/test_operators.py`) |
| D12 | DSA doc anchors: `operator_algebra.rst` ~:4210-4311 + `diffusion_1d.rst` DSA-seam | **operator_algebra.rst** (now `docs/theory/foundations/`, 4708 lines) has **no DSA section** — one passing mention :1594. The residual→DSA substrate discussion moved to `docs/theory/foundations/field_algebra.rst:519-529`. "The DSA seam" = `docs/theory/methods/diffusion_1d.rst:494-518` (label `diffusion-dsa-seam`). `diffusion_1d.rst:517` still xrefs operator_algebra "for the DSA-consumer discussion" — **content-drifted pointer** (fix at close-out) |
| D13 | 2.5b: LD-slab transpose = typed deferral | Extended by #310 C2: `has_transpose_kernel` is **derived True for LD** (`transport/spatial/linear_discontinuous.py:283`; registration machinery `scheme.py:694,726`). Transpose coverage today: DD + LD slab kernels; walk-side `has_transpose_walk` mesh-keyed (`loss_representation/__init__.py:410,674,1432,1489`) |
| D14 | 3c close-out targets "`discrete_ordinates.rst` DSA section" | That monolith is **dissolved** (#231): the SN theory is the part tree `docs/theory/methods/sn/` (14 pages). No acceleration/DSA page exists yet; scattered DSA mentions listed in §8 |

---

## 1. Assembly machinery (Phase 2's deliverable)

**The carrier.** `SparseAssembledOperator` @ `orpheus/numerics/assembled_operator.py:77`
— CSR carrier (COO duplicate-sum construction = the FEM scatter), flat-1-D `apply`
:140 (CSR matvec), **exact `apply_transpose`** :151, `shape` :136, densifying
`as_matrix` :182, idempotent `assemble` :176.

**The axis (three-layer surface).** `LinearOperator.is_assemblable` @
`numerics/operator.py:637` (default False — explicit override only);
`SupportsAssembly` Protocol :1072; `assemblable()` narrowing bridge :1128;
`MissingAssembly` refusal :424. Composites derive recursively: OperatorSum
`assemble = a+b` :1478, OperatorProduct `@` :1676, Scaled `*` :1812.
`FlattenedOperator.assemble` passthrough @ `numerics/flat_operator.py:145`.
`CoupledOperator.assemble` @ `numerics/coupled_system.py:1158` (block grid).
**R2 verified live**: `as_matrix` @ :868 delegates to `assemble().as_matrix()` when
assemblable (:961-967); the probing loop is `_as_matrix_by_probing` :970 — the
retained fuller-view oracle.

**Who implements `assemble()` (production):**
- Diffusion leaves: `LeakageOperator.assemble` @ `orpheus/diffusion/operators.py:441`
  (direct-structural from the conductance/closure attributes),
  `DiffusionBoundaryOperator.assemble` :626 (law-probing).
- Shared transport leaves: `MultiplicationOperator.assemble` @
  `transport/operators/multiplication_operator.py:258` (C, coefficient bit-exact);
  `IsotropicScattering.assemble` @ `transport/operators/isotropic_scattering.py:305`
  and `IsotropicN2N.assemble` :414 (group-impulse probing through the production
  einsum, `_assemble_iso_energy_operator` :161) — these are the **scalar-composite**
  (iso-energy) S; the **angular** `ScatteringOperator` does NOT assemble.
- SN: **no operator-surface assembly** (nothing in `orpheus/sn/operators/` overrides
  `assemble`; `is_assemblable` inherits False). The SN route is the standalone
  per-ordinate walk assembler:

**The 3a call.** `assemble_ordinate_blocks(sn_mesh, ordinate, *,
include_collision=True) → tuple[SparseAssembledOperator, ...]` @
`orpheus/sn/loss_representation/assembly.py:276`. One `(n_cells·cm × n_cells·cm)`
block **per energy group** for ONE ordinate; `include_collision=True` = the
sweep-invertible **L+C** (reaction_xs = Σ_t), `False` = pure L (honestly singular for
degenerate ordinates). Scope gates: **Cartesian only** (slab + 2-D; the
`SNMesh.streaming` accessor is the gate — curvilinear refuses), **linear schemes only**
(`is_linear` — DD and LD both declare it; `MissingAssembly` otherwise :317-324).
Coefficients extracted by **unit probes of the production `residual_kernel_batch`**
(`_probe_coefficient_blocks` :187) — zero stencil re-spelling; sweep-frame → global
frame conjugation via the single-sourced `octant_moment_frame_signs` (:366-369).
**The blocks are BULK-only at zero boundary inflow** — trace/boundary rows are
structurally excluded (the returns note :303-311).

So for 3a TODAY: the assembled Cartesian DD (L+C) on a slab mesh = loop ordinates ×
groups over `assemble_ordinate_blocks(sn_mesh, m)`. **−S is NOT available as an
assembled angular operator** — for the ℓ=0 moment reduction the scattering
contribution is the Σ_s0 energy matrix directly off `mat_xs` (nothing to assemble);
the dense angular kernel for cross-checks is `ScatteringOperator.full_scatter_kernel`
(`transport/operators/scattering.py:656`, an OperatorProduct R∘(Λ+N2N)∘M) via
probing `as_matrix` at small sizes.

**Equivalence gates (all pinned, `tests/sn/sweep/test_assembly_mode.py`):** G1
`assembled@x ≡ apply(x)` :153 (DD) / :401 (LD); G2 walk-order triangularity exact
:187 + LAPACK forward-substitution ≡ sweep :211 (DD) / :438 (LD, block-triangular
via LU); G3 probed-column pin :242; streaming-vs-collision diagonal :264; one-source
teeth (shared-kernel sign flip moves all three modes) :288; curvilinear refusal
:502. Transpose-side oracles: assembled-Mᵀ per-ordinate blocks @
`tests/sn/sweep/core/test_multi_d_reverse_walk.py:507` (DD 2-D) / :904 (LD 2-D).
Diffusion probed≡assembled gates: `tests/diffusion/test_operators.py`.

## 2. The unified walk (#280) — current shape

The 1-D walk executors live in `orpheus/sn/loss_representation/__init__.py`
(4983 lines): `_OneDimScanWalk` :2769 with the one orientation-parametrized frame —
`_dag_legs` :2800, `_loop_walk` :2850, `_reverse_traversal` :2753 (module level),
`sweep` :2923, `sweep_transpose` :2975, `_run` :3768, `_run_transpose` :4345. The
2-D wavefront sibling `_OctantWalk` :921. Capability predicates factorize along the
kernel×orientation axes: scheme `has_transpose_kernel`
(`transport/spatial/scheme.py:694`, derived-registration :726; DD True, **LD now
True per #310 C2** — `linear_discontinuous.py:283`) ∧ representation
`has_transpose_walk` (LR :410 protocol; mesh-keyed leaves :1432/:1489).

The inverse-adjoint wiring (2.5c): `SweepOperator` @
`orpheus/sn/operators/sweep_operator.py:83` — `apply` :104 (= inner.solve),
`is_adjointable` :130, `apply_transpose` :158 (= `inner.solve_transpose`, the
reverse-scan); `InvertibleOperator.solve_transpose` @
`orpheus/sn/operators/streaming.py:981`. The swap law `A.H.inverse() ≡
A.inverse().H` is an object identity: `_AdjointOperator` @
`numerics/operator.py:1155` with `is_invertible` :1237 + `inverse()` :1251. Gates:
`tests/sn/operators/test_inverse_adjoint_coherence.py` (swap law, G-reciprocity,
round-trip) and `tests/sn/operators/test_loss_transpose_solve.py` (reverse-scan
G1/G2/G3, all geometries). A_BB (the curvilinear ψ½ march) is
`RadialCharacteristicOperator` @ `orpheus/sn/operators/radial_characteristic.py`.
The DSA accelerator "wires into the unified walk" in the sense that both drivers
consume the ONE `WithinGroupSystem` splitting (§5) whose resolvent `.solve`/
`.apply`/`.apply_transpose` all run the same walk frames.

## 3. The diffusion family (the A_diff DSA consumes)

`orpheus/diffusion/` = `{augmented_mesh, boundary_realizer, method_space,
operators, solver}.py`.

- **Mesh**: `DiffusionMesh` @ `diffusion/augmented_mesh.py:122`;
  `from_material_mesh(material_mesh: MaterialMesh) → DiffusionMesh` :254 — NO extra
  parameters; BCs ride the axes; the docstring names the DSA path explicitly
  (:268-271: "An SNMesh promotes directly (it IS a MaterialMesh) — the DSA
  construction path (#2)"). `realize_boundary_law` :294 (the TransportMethod arm →
  `DiffusionMethodSpace.for_face` + `DiffusionBoundaryRealizer().realize`);
  `scalar_trace` :320 (per-face (J⁺,J⁻), quadrature-free); `full_field_space` :335
  (cached_property — the P7b pyright rationale).
- **Operators** (`diffusion/operators.py`): `LeakageOperator` :272 (`__init__(mesh)`
  :314, `apply` :407, `face_currents` :380, `assemble` :441);
  `DiffusionBoundaryOperator` :544 (`apply` :593, `assemble` :626). C/S/N2N come
  from the shared transport layer (§1).
- **Composite + resolvent** (`diffusion/solver.py`): `DiffusionSolver` :174;
  `self.loss = leakage + collision − scattering − boundary` :229-240 with
  `scattering = IsotropicScattering + IsotropicN2N` :235-238;
  `template = FullField.zeros(interior=ScalarFlux, boundary=ScalarBoundaryFlux,
  mesh=mesh)` :246-248; **the resolvent ruling holds verbatim**:
  `self.resolvent = MatrixInverseOperator(FlattenedOperator(self.loss,
  self.template))` :250-252 — eager LU (`matrix_inverse_operator.py:142` calls
  `inner.as_matrix`, `lu_factor` :160), NOT a structure-keyed `.inverse()`.
  Materialization now routes **through assembly** (D11).
- **Solve entry points**: `DiffusionSolver.solve_fixed_source` :287 (exact one-LU
  backsolve; EigenvalueSolver protocol), `compute_production_rate` :296,
  `compute_keff` :306; `DiffusionResult` :336; `solve_diffusion_1d` :355 →
  `power_iteration` :407.

## 4. The residual substrate

`evaluate_residual` @ `orpheus/sn/solver.py:231`:

```python
def evaluate_residual(
    loss_op: "LinearOperator",
    psi: "TimedFullField | CoupledField",
    q_ext: "FullField | CoupledField",
) -> "FullField | CoupledField":
```

Seedless (every Cartesian mesh): returns
`FullField(interior=AngularResidual.from_balance(...),
boundary=AngularBoundaryResidual.from_balance(...))` (:223-228 via the shared
`_system_a_residual` helper). Carrying mesh (R12a curvilinear): requires the coupled
pair + the 2×2 loss grid and returns `CoupledField[r_A, r_B]` with r_B =
`RadialCharacteristicInteriorResidual ⊕ RadialCharacteristicBoundaryResidual`
(:303-313); a bare System-A call on a carrying mesh **raises** (:324-329) — the
guard message names the DSA consumer. The docstring (:257-261) marks it "the
substrate the consistent-DSA low-order correction (#2) will consume … NOT in the
convergence path — additive." Sibling diagnostic: `boundary_vs_interior_split` @
`sn/solver.py:335`.

Foldable split (data API, zero production callers — Nexus-confirmed):
`ScatteringOperator.foldable_part` @ `transport/operators/scattering.py:978`,
`residual_part` :1016, `is_foldable_into_sigma_r` :1057; the ⚠#215 trap comment
:968-976 (pointers to the two theory sections, §8). `MaterialXSField` accessors @
`transport/mesh/material_xs_field.py`: `foldable_sig_s` :1148, `residual_sig_s`
:1167, `is_p0_diagonal_with_zero_n2n` :1189, `foldable_sigma` :1213 (returns
`{mid: (ng,) diag}`).

**`as_dsa_source` does NOT exist** in code (grep clean over `orpheus/` + `tests/`);
the only mentions are the docs note "lands WITH DSA #2"
(`docs/theory/foundations/field_algebra.rst:528` + stale `_build` copies).

## 5. The iteration layer (the plug point)

`orpheus/numerics/iteration.py` (1351 lines): `SupportsSeededApply` :166,
`Preconditioner = Callable[[np.ndarray], np.ndarray]` :198, `seeded_inverse` :239,
`_CarrierMatvecOperator(spla.LinearOperator)` :325 (the typed→scipy adapter),
`SourceIteration` :440, `KrylovAcceleration` :707, `KEigenvalue` :1000.

**`SourceIteration.solve`** (:592-699), the loop body per iteration:
1. `rhs = q_ext + Σ g.apply(psi)` (:657-659) — the lagged gains;
2. the ρ-honest STOP via the free identity `r_n = rhs_{n−1} − rhs_n`, checked
   BEFORE the next inverse apply (:668-673);
3. `psi = self.A_inv.apply(rhs, initial_guess=psi_prev)` (:682) — the resolvent;
4. displacement diagnostics (:689-697; not the stop).

**The DSA correction step injects between :682 and the loop end** (compute r,
low-order solve, correct ψ before the next rhs build). ⚠ `SourceIteration` has **no
hook/callback seam today** — the "SI+DSA wrap" is a new construct (wrapper or
sibling driver), not a parameter. Note the stop measures rhs-deltas; with a
corrected iterate the identity still measures the true equation residual of the
corrected ψ (it is Σg(ψ_{n−1}−ψ_n) under exact-M), but the wrap design should
re-derive this consciously.

**`KrylovAcceleration`** — see §10.

**`KEigenvalue`** (post-R8): `__init__(A, S, F, *, max_outer=500, keff_tol=1e-7,
flux_tol=1e-6, max_inner=1000, inner_tol=1e-8, eigenvalue_method="power")` :1089;
requires `A.is_invertible` :1115; builds `self._inner =
SourceIteration(seeded_inverse(A), S)` once :1143-1146; hardwired estimators
`compute_production_rate` :1212 / `compute_keff` :1237 (no injection kwargs).

**Who constructs what (the solve_sn chain):**
- `solve_sn` @ `sn/solver.py:1944` → `_as_sn_mesh` :2038 → `SNSolver(...)` :2040
  (class :837, EigenvalueSolver protocol; estimators :1173/:1204) →
  `power_iteration(solver, max_iter=max_outer)` :2049
  (`numerics/eigenvalue.py:203`; `direct_eigenvalue` :363).
- Inner: `SNSolver.solve_fixed_source` :1088 → `_solve_source_iteration` :1409 /
  `_solve_krylov` :1632, both consuming the ONE
  `build_within_group_system(sn_mesh, mat_xs, scattering_op=…)`
  (`sn/coupled_system.py:380`) → `WithinGroupSystem{loss, space, resolvent, gains}`
  :283 (seedless: resolvent = fused (L+C) `InvertibleOperator`, gains = `(S, B_a)`
  with B_a LAST; carrying: the 2×2 grids). Driver factories:
  `_within_group_si` @ `sn/solver.py:~806-830` (seedless arm: `_select_si_resolvent`
  G-S re-split + `_maybe_window` :493 → `SourceIteration(step, *gains)` :829;
  coupled arm `SourceIteration(system.resolvent.inverse(), *system.gains)` :812) and
  `_within_group_krylov` @ :458-490 → `KrylovAcceleration(LC, *gains, …)` :485.
- Fixed-source siblings: `solve_sn_fixed_source` :2857 → `_solve_fixed_source_si`
  :3001 / `_solve_fixed_source_krylov` :3210.
- Adjoint: `solve_sn_adjoint` :2348 builds `KEigenvalue(resolvent.H, gain.H,
  F_posed.H, …)` :2416.

## 6. TransportMethod protocol

`orpheus/transport/method.py`: `@runtime_checkable class
TransportMethod(Protocol[OpT_co])` :145-146. Members: `axes:
tuple[Axis1D, ...]` :178; `BOUNDARY_OPERATOR_REGISTRY: ClassVar[dict[str,
type[BoundaryTraceLaw]]]` :187; `bc` property :189-197; `realize_boundary_law(law,
face)` :199-212. The ONE generic body `resolve_boundary_conditions(method)` :215.

**The P7b deferral note is intact verbatim** @ :162-172: `full_field_space` is
deliberately NOT declared — both witnesses implement it as
`functools.cached_property` (which pyright rejects against a Protocol property
member), and it has no method-generic consumer yet; "**Declare it when the first
generic consumer arrives (the DSA driver, #2)**". The trigger fires exactly as the
roadmap planned. Witnesses: `SNMesh.full_field_space` and
`DiffusionMesh.full_field_space` (`diffusion/augmented_mesh.py:335`).

## 7. Sweep/scheme layout post-rename (supersedes R9)

| Symbol | Home @ `e4c1a81c` |
|---|---|
| `GeometryCoefficients` | `orpheus/sn/sweep/cache.py:118` |
| `CollisionCache` | `orpheus/sn/sweep/cache.py:365` |
| `DiamondDifference` | `orpheus/transport/spatial/diamond.py:81` (`key="diamond_difference"`) |
| `CellVisit` | `orpheus/transport/spatial/scheme.py:81` |
| `DiscretizationSchemeBase` | `orpheus/transport/spatial/scheme.py:527` (+ `has_transpose_kernel` :694/:726) |
| `LinearDiscontinuous` / `_LDCellTerms` | `orpheus/transport/spatial/linear_discontinuous.py` |
| `assemble_ubld` / inflow-axis kernels | `orpheus/transport/spatial/_ubld.py:183/:260/:316` |
| `CellBalanceTerms` | `orpheus/transport/spatial/cell_balance.py` |
| walk executors + assembler | `orpheus/sn/loss_representation/__init__.py` + `assembly.py` |
| scan + seed engine + closure | `orpheus/sn/sweep/{scan,psi_half_angle_seed,pole_angular_closure,pairing}.py` |

`tests/sn/spatial/` is dissolved; assembly gates live at `tests/sn/sweep/`.

## 8. Doc anchors (post-#231 part tree)

- **`docs/theory/foundations/operator_algebra.rst`** (4708 lines): the assembly
  axis (label `operator-algebra-assembly-axis`; dev-history row :4552-4580; §
  around :3843-3944). **The roadmap's ~:4210-4311 DSA block no longer exists on
  this page** — single passing DSA mention :1594 (the operator→scipy adapter
  aside).
- **`docs/theory/foundations/field_algebra.rst`**: the typed-residual → DSA
  substrate paragraph :519-529 ("the literal substrate the consistent DSA
  low-order correction (Issue #2) will consume … `as_dsa_source` lands WITH DSA
  #2"); residual diagnostics :495-518.
- **`docs/theory/methods/diffusion_1d.rst`**: `.. _diffusion-dsa-seam:` :494,
  section "**The DSA seam**" :496-518 — names `from_material_mesh`, the ℓ=0
  half-range boundary restriction under the shared |Ω·n|w metric, Marshak-for-free.
  ⚠ :517 xrefs operator_algebra "for the DSA-consumer discussion" — that
  discussion moved to field_algebra (content-drifted pointer; fix in 3's doc pass).
- **`docs/theory/methods/sn/`** (14 pages; **no acceleration/DSA page exists**).
  DSA-relevant sections: `slab_one_group.rst` — the σ_r-fold mismatch (eq label
  `si-sigma-r-fold-mismatch` :742; Key-Facts row :51; §~694-800: "DSA needs a
  consistent diffusion operator", foldable_sigma named as DSA input :795-797);
  `loss_representation.rst` — "The three consumption modes — solve, apply,
  assemble" :601, removal-σ eq (`loss-rep-removal-sigma`) :461, DSA's R·A·P moment
  reduction :777; `cartesian_multid.rst` :4010 + :4124-4128 (consistent-DSA
  forward-ref; the naive-FD DSA divergence spike); `placement.rst` :291-296;
  `slab_multigroup.rst` :453 (apply_transpose for adjoint/DSA posing);
  `solver.rst` = the driver/iteration page (KrylovAcceleration :42/:80/:85) — the
  natural home for the DSA theory section, or a new part-tree page beside it.
  **The roadmap 3c target "discrete_ordinates.rst DSA section" is stale** — close
  out into the part tree.

## 9. Anisotropic scattering + the isotropic projector inventory

`ScatteringOperator` (`transport/operators/scattering.py:356`) holds per-material
Legendre matrices Σ_{s,ℓ} up to `scattering_order` (0 = P0), (n,2n) matrices, and
precomputed real SH Y_ℓ^m at the quadrature directions; `BlockRole.BULK` :393;
`is_adjointable` True :437 (Sᵀ via the kernel). Companions in the same file:
`LegendreMomentScattering` :115, `N2NMomentOperator` :300.

**What the SI actually lags: the FULL S.** `build_within_group_system` places the
whole `ScatteringOperator` (all retained ℓ, cross-group P0, n2n) in `gains`; the SI
loop applies it to ψ every iteration (`iteration.py:657-659`). The foldable split
is wired NOWHERE in production (zero callers; "Data API only — the intended
consumer is a consistent DSA preconditioner (#2)", scattering.py:963-965). On 2-D
Cartesian the SI iterate is harmonic moments via `WindowedSweep`
(`_maybe_window` @ `sn/solver.py:493`; step = `P @ A.inverse()`, P sourced from
**S's own frame** so stored moments match S's projection term-for-term) — arm-1
slab DSA sees the full-angular iterate; the 2-D interaction is a design note.

**Moment-0-over-quadrature reductions — the complete spelling inventory** (DSA's R
must reuse, never re-spell):
1. `AngularFlux.integrate_angular()` @ `transport/fields/angular_flux.py:111` —
   THE canonical Σ wₙ ψₙ → `ScalarFlux`.
2. The scattering harmonic frame: `full_scatter_kernel` :656 (OperatorProduct
   R∘(Λ+N2N)∘M — M is the analysis face; its ℓ=0 row is pinned ==
   `integrate_angular` bit-exactly, :1222/:1245); `isotropic_kernel` :682
   (OperatorSum Σ_s0 + 2Σ_2n — the model-shared K_iso, an ENERGY kernel applied to
   φ, not itself an angular reduction).
3. `FissionOperator` = χ⊗νΣf ∘ kernel ∘ `integrate_angular`
   (`transport/operators/fission.py:401`, read :561; transpose = wₙ-weighted
   broadcast :438).
4. Boundary half-range family: `AngularBoundaryFlux.net_current(face)` @
   `transport/fields/angular_boundary_flux.py:128` and
   `ScalarBoundaryFlux.net_current(face)` @
   `transport/fields/scalar_boundary_flux.py:127` — the (J⁺,J⁻) reductions under
   the |Ω·n|w trace metric (the #290 ruling-2 seam the DSA boundary restriction
   consumes).
There is no named `P_iso` production object — "P_iso" appears only in prose
(scattering.py:970, streaming.py:233). For P (prolongation): the producer-side /W
isotropic factory is `AngularSourceSink.from_isotropic` (used at the solve_sn final
sweep, `sn/solver.py:~2058`).

## 10. Krylov preconditioner seam (#200 context)

`KrylovAcceleration` (`numerics/iteration.py:707-992`): GMRES on
`(A − Σᵢ gᵢ)ψ = q_ext` over whatever carrier the operators consume — the matvec is
the honest `A.apply(psi) − Σ g.apply(psi)` per call (:899-908), raveled at the
scipy boundary through `_as_scipy_linop`/`_CarrierMatvecOperator` (:325; typed
composites ravel via `to_flat`, so the flat vector carries bulk ⊕ trace (⊕ ψ½)).
`scipy.sparse.linalg.gmres(A_scipy, b, x0, M=M_scipy, rtol=tol, atol=0.0,
maxiter, restart=min(restart, n), callback_type='pr_norm')` :943-950; ERR-053
non-convergence warning :978-990.

**The preconditioner seam exists TODAY**: `preconditioner: Preconditioner | None`
kwarg (:814; type = `Callable` :198). Semantics (docstring :740-764): a LEFT
preconditioner M ≈ (A − Σg)⁻¹ — a typed callable (carrier → carrier; the adapter
handles ravel, :912-916). Default when `None`: `seeded_inverse(A).apply` (the
sweep — the Adams & Larsen transport-corrected M) if `A.is_invertible`, else
identity (:844-849). **The SN production path overrides with an explicit
identity**: `_within_group_krylov` @ `sn/solver.py:458-490` passes
`preconditioner=lambda q: q` (:487, "explicit identity — issue #200 tracks
re-enable") and `restart=n_dof` (:489). The DSA-preconditioned-Krylov posture (3b
wiring 4) is literally: replace that lambda with the DSA callable
(e.g. sweep + P·A_diff⁻¹·R correction) — no structural work needed at this seam.

---

## Blockers / loud flags for the 3a–3b shape

**No hard blocker.** Five scope facts to carry into the plan-of-record:

1. **3a boundary row is guaranteed-differ, not maybe-differ.** The assembled SN
   blocks are bulk-only at zero inflow (boundary/trace rows structurally excluded,
   `assembly.py:303-311`), while A_diff's assembly INCLUDES its boundary closure
   rows (B assembles). The R4 comparison is interior-stencil vs interior-stencil;
   the boundary row lands in the "characterized boundary row" branch by
   construction. Plan 3a's matrix diff accordingly.
2. **−S is not assembled in angular form** — and doesn't need to be for 3a: the
   ℓ=0 reduction of S is the Σ_s0 energy matrix straight off `mat_xs`
   (`foldable_sigma`/`sig_s_legendre`); the dense angular kernel
   (`full_scatter_kernel` + probing `as_matrix`) exists for cross-checks.
3. **SI has no correction hook.** The "SourceIteration wrap" (3b wiring) is a new
   driver construct; the injection point is after the resolvent apply
   (`iteration.py:682`) and before the next rhs build. The Krylov side needs zero
   structural work (§10).
4. **`evaluate_residual` on carrying meshes demands the coupled pair** (loud
   ValueError, `sn/solver.py:324-329`). Cartesian DSA scope = seedless = clean
   2-block; the guard only bites if curvilinear DSA is attempted (already OUT per
   3b.5).
5. **Spelling drifts for any plan text**: `FullField(interior=…, boundary=…)` (not
   `bulk=`); the within-group factory is `build_within_group_system` →
   `WithinGroupSystem` (not the `(L+C, S, B)` triple); the 3c doc target is the
   `docs/theory/methods/sn/` part tree (not `discrete_ordinates.rst`); and
   `diffusion_1d.rst:517`'s operator_algebra xref is content-drifted (the DSA
   discussion lives in `field_algebra.rst`) — fold the fix into the Phase-3 doc
   pass.
