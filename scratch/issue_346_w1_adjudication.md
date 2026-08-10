# #346 W1 — adjudication table for the 91 `.rst` dead candidates

`[M]` 2026-08-10 @ `9e550eb5`. Produced by `$CLAUDE_JOB_DIR/tmp/w1_list.py`
(MRO-aware: it resolves `Cls.attr` through the real class object, so a member
inherited from a base or set as `self.attr: T` is **not** reported dead).

**91 sites / 64 distinct**, split into five groups. Every ruling below was
measured against the tree, not inferred from the prose.

⚠ **The instrument that produced the first draft of this list was wrong three
times** (a bare suffix index over the graph). It said the `BC` enum was retired
— it is not, `orpheus.geometry.mesh.BC.vacuum` resolves — and it said five
`SNMesh` members were gone when they are inherited. Trust the MRO column, not
a tail match.

---

## Group 0 — ALIVE, merely unqualified (23 sites). Do NOT "fix"; QUALIFY.

These resolve today. They render as unlinked plain text only because the role is
relative. The exact target is known, so the edit is mechanical: replace the
relative target with `~<full.path>` (the `~` keeps the displayed text short).

| site(s) | role + target | full path |
|---|---|---|
| `docs/api/geometry.rst:163` ×3 | `:obj:`BC.vacuum`` / `.reflective` / `.white` | `orpheus.geometry.mesh.BC.<member>` |
| `docs/theory/methods/sn/boundary_conditions.rst:26,27,28` | `:attr:`geometry.mesh.BC.<m>`` | `orpheus.geometry.mesh.BC.<m>` |
| `docs/theory/verification/sn.rst:98` | `:attr:`geometry.mesh.BC.vacuum`` | same |
| `boundary_conditions.rst:6968,7197,7368`; `methods/sn/boundary_conditions.rst:617`; `methods/sn/index.rst:176` | `:attr:`SNMesh.axes`` | `orpheus.sn.mesh.augmented_mesh.SNMesh.axes` |
| `boundary_conditions.rst:7458` | `:attr:`SNMesh.coord`` | `…SNMesh.coord` |
| `boundary_conditions.rst:7494,7509,7518,7533,7544` | `:attr:`SNMesh.spatial_shape`` | `…SNMesh.spatial_shape` |
| `boundary_conditions.rst:7510,7514,7531` | `:attr:`SNMesh.axis_widths`` | `…SNMesh.axis_widths` |
| `boundary_conditions.rst:7546,7700` | `:attr:`SNMesh.volume_measure`` | `…SNMesh.volume_measure` |

## Group 1 — ASPIRATIONAL: the type does not exist (16 sites)

`[M]` all ABSENT from the graph (977 modules / 1326 classes / 5720 functions).
These pages describe a **design**, and a `:class:` role asserts existence — the
same defect class as the five reserved-placeholder `:mod:` roles fixed in
`7d7661b2`. ⚠ In `layering.rst` they sit in a table **beside a real one**
(`:class:`~orpheus.numerics.operator.LinearOperator``), so the false ones are
camouflaged by a true neighbour.

**Fix: double-backtick literal + an inline status marker** making it
unmistakable that the type is a design intent, not shipped code. Follow the
wording pattern landed for the reserved placeholders in `7d7661b2`.

| site | targets |
|---|---|
| `docs/architecture/layering.rst:122` | `Problem`, `Solver` |
| `:135` | `Eigenproblem` |
| `:136` | `PowerIteration`, `Arnoldi`, `TimeStepper` |
| `:138` | `CriticalityProblem`, `AlphaEigenproblem` |
| `:139` | `FixedSourceProblem`, `InitialValueProblem` |
| `:143` | `SweepPreconditionedSolver` |
| `:145` | `Problem` |
| `docs/theory/methods/sn/index.rst:694,697` | `Step`, `ExponentialCharacteristic` |
| `docs/theory/methods/sn/index.rst:846,847` | `Step`, `ExponentialCharacteristic` — the prose already says **"Wave C-extension *will ship*"**, so the tense is right and only the role is wrong |

⚠ `LinearDiscontinuous` on `:846` is the exception — it **EXISTS** at
`orpheus.transport.spatial.linear_discontinuous.LinearDiscontinuous`. Qualify
that one; literal-ise only its two absent siblings. A sentence listing three
names of which one is real is exactly where a blanket edit goes wrong.

## Group 2 — REPOINTABLE: the symbol exists elsewhere (14 sites)

| site | current target | ✅ repoint to |
|---|---|---|
| `references/peierls.rst:1895` | `peierls_geometry.build_volume_kernel` | `orpheus.derivations.continuous.peierls_nystrom.geometry.build_volume_kernel` |
| `references/peierls.rst:2182` | `peierls_geometry.build_white_bc_correction` | `…peierls_nystrom.geometry.build_white_bc_correction` |
| `references/fn_method.rst:1107`, `references/sood_registry.rst:453` | `...slab.solve_fn_slab_bare_critical` | `orpheus.derivations.continuous.fn_method.slab.one_group.solve_fn_slab_bare_critical` |
| `references/sood_registry.rst:442` | `...fn_method.multi_group.compute_kinf_1g` | `…fn_method.multi_group.k_inf.compute_kinf_1g` (note the extra `k_inf` segment) |
| `references/sood_registry.rst:464` | `...fn_method.sphere.solve_fn_sphere_bare_critical` | `…fn_method.sphere.one_group.solve_fn_sphere_bare_critical` |
| `references/singular_eigenfunction.rst:1370` | `...slab.solve_case_method_slab_critical` | `orpheus.derivations.continuous.singular_eigenfunction.slab.one_group.solve_case_method_slab_critical` |
| ⭐ `references/fn_method.rst:1559` | `...fn_method.core.x_function` | `orpheus.derivations.continuous.singular_eigenfunction.core.x_function` — **a DIFFERENT package** than the prose claims. Verify the sentence still reads true after repointing; if it asserts the evaluator lives in `fn_method`, that claim is false and the prose changes too. |
| `collision_probability.rst:1328` | `~geometry.factories.pwr_slab_half_cell` | `orpheus.geometry.structured_geometry.StructuredGeometry.pwr_slab_half_cell` — it is a **method now**, so the role becomes `:meth:` |
| `foundations/operator_tensor_network.rst:945` | `TensorProductOperator.assert_separable` | `orpheus.numerics.operator.SumOfTensorProductsOperator.assert_separable` — ⚠ READ the sentence first: if the claim is about the *tensor-product* operator and the method now lives on the *sum*, the claim may have shifted scope, not just moved. |
| `foundations/boundary_conditions.rst:3558`, `operator_tensor_network.rst:905` | `np.array_equal` | `numpy.array_equal` — `np` is an import **alias**; Sphinx cannot resolve an alias |
| `operator_tensor_network.rst:912` | `np.add` | `numpy.add` |
| `foundations/boundary_conditions.rst:6123` | `:class:`ReflectiveBoundary.apply`` / `:class:`WhiteBoundary.apply`` | role misuse (`:class:` on a method) **and** unqualified: `:meth:`~orpheus.geometry.boundary.reflective.ReflectiveBoundary.apply`` and the `white` sibling |

## Group 3 — ⛔ PRESENT-TENSE-FALSE: prose that MISINFORMS (5 sites) — the real payload

Not literals. The prose asserts a **current** API that does not exist, so per
the articulation-enforcement rule (a falsified doc is a BUG) the sentence must
be corrected against the tree, not merely de-linked.

1. **`docs/theory/conventions/indexing_and_layout.rst:1332–1341`** — a section
   headed *"Factory methods (`SNMesh`)"* states *"The `SNMesh` carries factory
   methods that allocate zero-initialised instances"* and lists
   `zeros_angular_flux` / `zeros_scalar_flux` / `zeros_boundary_flux`.
   `[M]` `[n for n in dir(SNMesh) if "zero" in n.lower()] == []` — **none of the
   three exists, and neither does the capability as described.** Investigate what
   replaced zero-allocation (start from `orpheus/transport/fields/`,
   `FullField`/`AngularFlux` constructors and any `zeros`/`empty_like` classmethod)
   and rewrite the section to the real mechanism. If allocation moved off the mesh
   entirely, say so — that is the architecturally interesting fact.
2. **`docs/theory/foundations/boundary_conditions.rst:7461`** and
   **`docs/theory/methods/sn/index.rst:180`** — both say, present tense, *"After
   C5.1, `SNMesh.mesh` **is** inbound provenance only … `None` when built from
   axes"*. `[M]` `hasattr(SNMesh, "mesh") is False`. Determine whether the
   provenance field was renamed or removed outright; if removed, the whole
   "inbound provenance" paragraph is stale and must go past-tense or be deleted.
   ⚠ `boundary_conditions.rst:7459` has a companion literal ``` ``mesh.coord`` ```
   in the line above — check whether that is also stale.

## Group 4 — RETIRED: past-tense history → literals (33 sites)

`[M]` all ABSENT. Per `.claude/rules/coding-standards.md`, discriminate **by
tense**: past-tense history is correct and STAYS as a literal; any
**present-tense** claim or **imperative** joins Group 3 instead. Several of
these read like history already ("Before C4 …", "Issue #188 then collapsed …",
"The former `BickleyTables` …") — but **check each one**, do not assume.

| target | sites |
|---|---|
| `SNMesh._resolve_one` | `boundary_conditions.rst:6761, 6797, 7008, 7088, 7896` |
| `_resolve_one` | `boundary_conditions.rst:7224` |
| `peierls_slab` | `peierls_nystrom.rst:138, 763, 1102, 1704, 1977, 2448` |
| `peierls_cylinder` | `peierls_nystrom.rst:2448` |
| `_quadrature_recipes` | `peierls_nystrom.rst:6851` |
| `BickleyTables` | `peierls_nystrom.rst:8514` |
| `_chord_half_lengths` | `collision_probability.rst:1183, 1194` |
| `_composite_gauss_legendre` | `collision_probability.rst:1483` |
| `_ki4_lookup` | `collision_probability.rst:1593` |
| `_build_ki_tables` | `collision_probability.rst:2007` |
| `Delta2` | `collision_probability.rst:870` |
| `_outflow_at_boundary_for_sphere` | `curvilinear_numerics.rst:988` |
| `_sweep_1d_cylindrical` | `curvilinear_numerics.rst:1756` |
| `_sweep_1d_cumprod` | `slab_one_group.rst:413` |
| `StreamingOperator._apply_2d_cartesian` | `operator_tensor_network.rst:867` |
| `Spectrum.from_problem` | `singular_eigenfunction.rst:573` |

⚠ **`peierls_slab` is 6 sites and several are IMPERATIVE**, e.g.
`peierls_nystrom.rst:1102` — *"modifications to `peierls_slab` should preserve
…"*. An imperative naming a retired module is a MUST-FIX (a contributor
following it hits `ImportError`), so those belong in Group 3, not here. The
module's successor is in `orpheus/derivations/continuous/peierls_nystrom/`
(`slab.py` / `cylinder.py`) — repoint where the instruction is still meant to
be followed.

## Group 5 — role misuse on a local VALUE, not a symbol (5 sites)

`:meth:`B.apply``, `:meth:`bc.apply``, `:meth:`bc_outer.apply``,
`:meth:`A.apply``, `:meth:`F.apply``, `:attr:`psi.boundary.face_view``. These
name **instances** in the surrounding prose (`B` the boundary operator, `A` the
loss operator, `psi` a field), not importable symbols — so no qualification can
ever make them resolve.

**Fix: double-backtick literals.** Do NOT invent a class to point them at: the
sentence is about a particular value, and repointing would change its meaning.

| site | target |
|---|---|
| `foundations/boundary_conditions.rst:3322` | `B.apply` |
| `foundations/boundary_conditions.rst:6592` | `bc.apply` |
| `methods/sn/curvilinear_one_group.rst:2018` | `bc_outer.apply` |
| `methods/sn/slab_multigroup.rst:585` | `A.apply` |
| `methods/sn/solver.rst:137` | `F.apply` |
| `foundations/operator_tensor_network.rst:872` | `psi.boundary.face_view` |

---

## Acceptance

1. `.venv/bin/python tools/check_docstring_xrefs.py` → `DEAD TARGETS: 0`.
2. `.venv/bin/python -m sphinx -E -W --keep-going docs docs/_build/html` → exit 0,
   **0 warnings** (`-E` is mandatory, see `[[lessons-L36]]`).
3. `.venv/bin/python -O -m pytest tests/test_docstring_xrefs.py -q` → all green.
4. The Group-0 count of gate-decided `docs/` roles goes **up** (it was 4690 of
   6032 at `9e550eb5`); re-run `$CLAUDE_JOB_DIR/tmp/w1_list.py` and confirm the
   91-site list shrinks to the Group-1/4/5 literals only.
