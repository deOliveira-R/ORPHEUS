# refgen_layers: the architecture constraint for Q-R1 (#405 step 2), HEAD 901f64ca, 2026-09-24

## 1. The contract as the linter enforces it [M]
- tests/gates/test_layer_imports.py:71 `"derivations": L2_PACKAGES | L3_PACKAGES` is the ONLY forbidden set. numerics (L1), geometry and data (input) are NOT forbidden: derivations may import them.
- Relative imports and the bare `from orpheus import X` form are invisible to `_check_module` (it keys on `orpheus.` prefixed absolute names). A grep for relative imports climbing out of derivations found 0 (4 hits with `....` all stay inside derivations).
- TYPE_CHECKING tolerance applies to L1/L2 sources only, not derivations (a TC import of transport in derivations is still a violation; mms/sn.py:85 is covered by the whitelist).
- WHITELIST (test lines 74-...): 4 entries; the doc (layering.rst "At the time of P3.1's landing") lists 3 (the mms/sn.py -> transport entry came later).
- Census, predicate `grep -rlE "(from|import) orpheus\.<pkg>\b"` over orpheus/derivations/**/*.py (137 files, __pycache__ excluded), files importing:
  numerics 3, geometry 9, data 10, transport 1 (mms/sn.py, whitelisted, lazy), moc 1 (mms/moc.py, whitelisted, module level), cp 1 (sood_registry/builders.py, whitelisted, lazy), sn 1 (angular_differencing.py:286, a COMMENT, not an import), diffusion 0, others 0.
  Positive control: the census finds the three whitelisted edges it should.
- STALE whitelist entry: ("derivations/continuous/cases/diffusion.py", "diffusion") -- the file has 0 `orpheus.` imports (grep `orpheus\.` hits only docstring roles at :35, :38, :727). The entry exempts nothing.

## Rationale text (quoted)
- layering.rst:103 "**L0** sits BELOW **L1** in the import hierarchy. ... They have *less* algorithmic knowledge than the production primitives in orpheus.numerics (a SymPy expression is structurally simpler than a numpy iteration). Production code that needs a structurally independent reference imports L0 — the L3-uses-L0 pattern is documented in /theory/verification/index."
- layering.rst WHITELIST para: "Each entry is a Branch-1-uses-production-as-a-black-box benchmark — the reference imports a production solver to *cross-check* a reference value, NOT to share algebra. These are categorically different from algebra-sharing imports (which would be structurally contaminating ...). The retirement trigger for each is the module's migration to a method-side test or to an external benchmark harness."
- Contradiction: layering.rst:61 "each layer may import only the layers below it" + the table places L0 lowest, so the prose implies derivations import nothing; the linter dict (reproduced verbatim in the same page) allows numerics+geometry+data. The input-layer note ("every layer (including L1) may consume") is what licenses it.

## 3. MaterialMesh
- orpheus/transport/mesh/material_mesh.py:122, layer L2 (transport). Forbidden to derivations.
- ctor `MaterialMesh(mesh: Mesh1D | Mesh2D, materials: Materials | Mapping[int, Mixture], sigma_t_cell=None)`; everything funnels into `_init_data(axes, mesh, mat_map, materials, sigma_t_cell)`. No general `from_axes` (deliberately deferred, :290 "NOTE: a GENERAL axis-native ..."). No to_dict/spec/serialisation. Content identity exists: `_identity_key` / `_contractibility_key`, `__eq__`/`__hash__` by content (usable as a cache key).
- Minimal plain data: a Mesh1D (edges via StructuredGeometry + RegionMesh refinement, mat_ids, coord, BCs) or Mesh2D (+ mat_map), and {id: Mixture}. ng derives from materials. Groups are not a separate datum.
- Uplifts: SNProblem.from_material_mesh(mm, quadrature, scheme=None) (sn/problem.py:796), DiffusionMesh.from_material_mesh (diffusion/augmented_mesh.py:258). None for CP/MoC/MC.
- MaterialMesh is DISCRETISED (cell edges): a reference can only produce one given a refinement parameter, which is the test's choice.

## 4. Precedents at L0 / input
- The method- AND mesh-independent problem is ALREADY an input-layer type: geometry/structured_geometry.py StructuredGeometry(geometry: "SLB"|"CYL"|"SPH", regions: tuple[Region(mat_id, outer_thickness_cm)], bcs: tuple[BC | BoundaryTraceLaw]) + {id: Mixture}. 1-D only. Consumed by 8 derivation files legally. Pipeline: StructuredGeometry + RegionMesh -> Mesh1D.from_geometry -> MaterialMesh(mesh, materials) -> SNProblem.from_material_mesh(mm, quad).
- sood_registry/builders.py:35 `build_mesh(case, n_cells=64) -> Mesh1D` + build_materials: this is option (b) already shipped (the test calls MaterialMesh itself).
- ProblemSpec (derivations/common/continuous_reference.py:130): materials dict[int, Any], geometry_type Literal (incl. "cartesian-2d", "homogeneous"), geometry_params: untyped dict, boundary_conditions: dict[str, Literal BC], is_eigenvalue, n_groups. A weaker parallel twin of StructuredGeometry+BC: stringly BCs duplicate geometry.BC; geometry_params keys seen in consumers: fuel_height 5, length 3, mms_case 3, refl_height 2. 19 `ProblemSpec(` sites in 9 derivation files; 55 test/orpheus sites read `.problem.<field>`. No adapter ProblemSpec -> Mesh exists.
