# P1 monomorphic-domains blast radius — read-only audit

**Audit HEAD**: `3d1a9b39` on `refactor/operator-strategy-layers`, working tree
clean under `orpheus/` and `tests/` (uncommitted files are all `.claude/` + `scratch/`).

**Nexus graph provenance**: built at `8d996f53` on branch
`refactor/operator-naming-honesty`, `git_dirty: true`. `8d996f53` IS an ancestor
of HEAD. Delta `8d996f53..HEAD` touches exactly ONE file under `orpheus/`
(`orpheus/derivations/discrete/sn/sweep_acyclicity.py`) — every other commit is
`docs(plan)` or `tests/`. **The graph is therefore trustworthy for `orpheus/`
production structure and slightly stale for `tests/`** (16 test commits landed
after the build). All test-side counts below come from grep, not the graph.

**Method**: three searches per the retirement-audit rule — (i) Nexus
`callers`/`impact`/`graph_query`, (ii) text grep over `orpheus/ tests/ scratch/
docs/ examples/`, (iii) an **AST census** (`ast.Call` walk over every `.py` in
those trees) that catches direct `T(...)` constructors AND multi-line calls AND
tells whether the space kwarg is present. The AST census is the authority for
counts; grep and the graph are the cross-checks.

---

# (1) Construction-site census

## 1.0 The three routes, and why "bare" splits in two

`MultiplicationOperator.from_mesh` is the ONLY factory with an implicit space
fallback: `multiplication_operator.py:335-337`

```python
if space is None:
    space = getattr(mesh, "full_field_space", None)
```

So a `from_mesh(σ, sn_mesh)` call with NO `space=` kwarg still ends up
space-CARRYING (an `SNMesh` always has `full_field_space`). `from_solver_data`
(Scattering, Fission) has **no** such fallback — omitting `full_field_space=`
gives a hard `None`. The plain dataclass ctor likewise gives a hard `None`.

Bucket (a) below therefore = *explicit kwarg* **or** *from_mesh on a
space-bearing mesh*; bucket (b) = the space really is `None` at runtime.

## 1.1 Production (`orpheus/`) — every site, exhaustive

| bucket | operator | site | route | space source |
|---|---|---|---|---|
| **a** | MultiplicationOperator | `orpheus/diffusion/solver.py:230` | ctor | `space=space` |
| **a** | MultiplicationOperator | `orpheus/sn/coupled_system.py:376` | ctor | `space=sn_mesh.full_field_space` |
| **b** | MultiplicationOperator | `orpheus/homogeneous/solver.py:143` | `from_mesh` | `mat_xs.mesh` = a `MaterialMesh` → `getattr` → **None** |
| **a** | ScatteringOperator | `orpheus/sn/solver.py:1035` | `from_solver_data` | `sn_mesh.full_field_space` |
| **a** | ScatteringOperator | `orpheus/sn/coupled_system.py:456` | `from_solver_data` | `sn_mesh.full_field_space` (via local at `:452`) |
| **a** | ScatteringOperator | `orpheus/transport/operators/scattering.py:1009` | ctor (`foldable_part`) | `self.full_field_space` (propagates) |
| **a** | ScatteringOperator | `orpheus/transport/operators/scattering.py:1050` | ctor (`residual_part`) | `self.full_field_space` (propagates) |
| **a** | FissionOperator | `orpheus/diffusion/solver.py:241` | `from_solver_data` | `full_field_space=space` |
| **a** | FissionOperator | `orpheus/sn/solver.py:1041` | `from_solver_data` | `sn_mesh.full_field_space` |
| **a** | FissionOperator | `orpheus/sn/solver.py:2363` | `from_solver_data` | `sn_mesh.full_field_space` |
| **b** | FissionOperator | `orpheus/homogeneous/solver.py:193` | `from_solver_data` | omitted → **None** |
| **a** | IsotropicScattering | `orpheus/diffusion/solver.py:236` | ctor | `space=space` |
| **b** | IsotropicScattering | `orpheus/homogeneous/solver.py:146` | ctor | omitted → **None** |
| **b** | IsotropicScattering | `orpheus/transport/operators/scattering.py:709` | ctor | omitted → **None** |
| **a** | IsotropicN2N | `orpheus/diffusion/solver.py:237` | ctor | `space=space` |
| **b** | IsotropicN2N | `orpheus/homogeneous/solver.py:146` | ctor | omitted → **None** |
| **b** | IsotropicN2N | `orpheus/transport/operators/scattering.py:709` | ctor | omitted → **None** |

**Production totals** — 17 sites: **11 in bucket (a), 6 in bucket (b)**, and the
6 live in exactly **TWO** files:

* `orpheus/homogeneous/solver.py` — `:143` (C), `:146` (K_iso pair), `:193` (F).
  The entire meshless path. See §2.
* `orpheus/transport/operators/scattering.py:709` —
  `ScatteringOperator.isotropic_kernel`, inside the composite that DOES carry
  `self.full_field_space`. **This one is free**: the space is in scope and simply
  not threaded (see §5, surprise S-1).

## 1.2 Tests (`tests/`)

AST totals — bucket (a) counts the `from_mesh`-on-`SNMesh` implicit thread:

| operator | (a) explicit kwarg | (a) implicit via `from_mesh(σ, SNMesh)` | (b) genuinely `None` | total |
|---|---|---|---|---|
| MultiplicationOperator | 12 | 114 | **9** | 135 |
| ScatteringOperator | 2 | — | **11** | 13 |
| FissionOperator | 1 | — | **11** | 12 |
| IsotropicScattering | 4 | — | **12** | 16 |
| IsotropicN2N | 1 | — | **10** | 11 |
| **total** | **20** | **114** | **53** | **187** |

The 114 implicit sites are `MultiplicationOperator.from_mesh(σ, <sn|sn_mesh|sn1|sn2|sph>)`
— the SN test idiom. They need **no edit** if `from_mesh` keeps deriving the space
from the mesh; they need an edit only if the fallback `getattr` is removed.

### bucket (b) test sites, per file (53 total)

MultiplicationOperator (9):
`tests/diffusion/test_operators.py:641` (mesh is a `DiffusionMesh` fixture — verify),
`tests/numerics/test_matrix_inverse_operator.py:206,265` (`mat_xs.mesh`, MaterialMesh),
`tests/sn/architecture/test_monomorphic_leaves.py:692` (`mat_xs.mesh`),
`tests/sn/sweep/cartesian_2d/test_2d_l2_matvec_correctness.py:211,312` (`mesh`),
plus 3 plain-ctor sites in `tests/transport/test_multiplication_operator.py` and
`tests/sn/_test_helpers.py` — see the JSON dump for exact lines.

ScatteringOperator (11): `tests/sn/architecture/test_monomorphic_leaves.py:1090`,
`tests/sn/operators/test_operators_apply_typed.py:187,286,321,429`,
`tests/sn/operators/test_scattering_operator.py` (4 ctor sites),
`tests/sn/operators/test_capability_survival.py`, `tests/sn/operators/test_green_operator_sn.py`.

FissionOperator (11): `tests/diffusion/test_operators.py:650`,
`tests/homogeneous/test_homogeneous.py:285,359`,
`tests/sn/architecture/test_monomorphic_leaves.py:696,1094`,
`tests/sn/operators/test_operators_apply_typed.py:140,156,291,326,428`,
`tests/sn/operators/test_capability_survival.py` (1 direct ctor).

IsotropicScattering (12) / IsotropicN2N (10):
`tests/diffusion/test_operators.py` (3/2),
`tests/numerics/test_matrix_inverse_operator.py` (3/3),
`tests/sn/operators/test_capability_survival.py` (1/1),
`tests/sn/operators/test_isotropic_scattering.py` (5/4).

## 1.3 `scratch/`, `docs/`, `examples/`

* `examples/` — **0** sites (no `.py` mentions at all).
* `scratch/` — **0** executable sites. The 21/41/31/18/19 grep mentions are all
  prose in the three uncommitted review maps
  (`scratch/review_map_reaction_operators.md`, `review_map_agnostic_diffusion.md`,
  `review_map_sn_assembly.md`) — they quote the same sites this audit measures.
* `docs/` — **ONE** executable-looking snippet:
  **`docs/theory/foundations/infinite_medium.rst:1061`** —
  `production = FissionOperator.from_solver_data(mat_xs=mat_xs)` — a bucket-(b)
  literal inside a `.. code-block::`. It is the doc twin of
  `orpheus/homogeneous/solver.py:193` and moves with it.
  (Grep also hits `docs/theory/foundations/field_algebra.rst:329` and
  `docs/theory/verification/sn.rst:3278` — both are `from_mesh` on FIELD types,
  not on these operators; not in scope.)

## 1.4 Graph vs grep — do they disagree?

**Yes, in one systematic way, and grep wins.** `mcp__nexus__callers` on
`py:class:…MultiplicationOperator` returns 16 (2 production + 14 test), missing
`orpheus/homogeneous/solver.py:143` — because that site goes through the
**classmethod** `from_mesh`, which is a SEPARATE graph node, and `callers` on the
class only sees direct-ctor edges. Same shape for `ScatteringOperator` /
`FissionOperator`, whose production sites are almost all `from_solver_data`.
A class-level `callers` query therefore **systematically under-reports exactly the
bucket-(b) sites** (they are the ones that route through a factory). The AST
census found 204 sites total vs the graph's per-class handfuls.
No case was found where the graph saw a site grep missed.


---

# (2) The meshless / homogeneous path

## 2.1 Mesh hierarchy — where `full_field_space` actually lives

`orpheus/transport/mesh/material_mesh.py:82` `class MaterialMesh` is the BASE;
the two augmented meshes subclass it:

* `orpheus/sn/mesh/augmented_mesh.py:94` `class SNMesh(MaterialMesh)` —
  `full_field_space` at `:919` (`cached_property`, always non-None)
* `orpheus/diffusion/augmented_mesh.py:122` `class DiffusionMesh(MaterialMesh)` —
  `full_field_space` at `:335` (`cached_property`, always non-None)

`grep -rn "def full_field_space" orpheus/` returns EXACTLY those two. The base
`MaterialMesh` has none — which is precisely why
`multiplication_operator.py:336` must use `getattr(mesh, "full_field_space", None)`.

## 2.2 Does `MaterialMesh` expose ANY `FunctionSpace`-valued attribute today?

**No.** Its entire public surface (`material_mesh.py`):

| member | line | kind |
|---|---|---|
| `axes`, `mesh`, `mat_map`, `materials` | `:125` ctor / `:143` `_init_data` | data |
| `nx` | `:211` (set in `_init_data`) | `int` sugar for `spatial_shape[0]` |
| `ng` | `:330` | `int` |
| `volumes` | `:365` | ndarray `(*spatial)` |
| `volume_measure` | `:370` | a `DiscreteMeasure`, **not** a `FunctionSpace` |
| `areas` | `:402` | ndarray `(nx+1,)`, **raises on 2-D** |
| `ndim` | `:423` | `int` |
| `spatial_shape` | `:428` | `tuple[int, ...]` — "the canonical dim-agnostic shape descriptor" |
| `material_xs_field()` | `:439` | `MaterialXSField` |

The only near-space object is `volume_measure` (a `DiscreteMeasure`). Grepping
`Space|space` in the whole file yields **only prose** (`:7`, `:32`, `:241`) — no
`FunctionSpace` import, no space-valued member.

**Shape-valued attributes DO exist and are sufficient**: `ng` + `spatial_shape`.
For the meshless carrier `MaterialMesh.from_materials` (`:234-284`) sets
`axes=(AxisMesh(edges=[0.0, 1.0]),)`, so `spatial_shape == (1,)` and
`(ng, *spatial_shape) == (ng, 1)` — **exactly the `basis_shape` the homogeneous
solver hand-writes.** `volumes` is `[1.0]` (unit cell), so the metric weights are
available too. **A `FunctionSpace(name=…, shape=(ng, *spatial_shape),
inner_product_weights=volumes)` is constructible on a bare `MaterialMesh` today
with zero new data** — cf. how `DiffusionMesh.full_field_space`
(`augmented_mesh.py:361-367`) builds its bulk block from exactly `self.volumes`,
`self.ng`, `self.spatial_shape`. The only thing it adds that `MaterialMesh`
cannot supply is the **trace** block (`self.scalar_trace`), and
`FullFieldSpace.interior_space` / `.trace_space` are both `Optional` fields
(`full_field_space.py:191`, `:194`) — a trace-free composite is representable.

## 2.3 `MaterialXSField`

Also **no** `FunctionSpace`. Shape metadata only, both read-through to the mesh:

* `orpheus/transport/mesh/material_xs_field.py:798-801` `ng`
* `orpheus/transport/mesh/material_xs_field.py:803-811` `spatial_shape`

It does NOT expose a single `(ng, ncells)` tuple; the caller composes
`(ng, *spatial_shape)`.

## 2.4 Other production entry points building transport operators on a `MaterialMesh`

**There are none besides `orpheus/homogeneous/solver.py`.** Evidence:

* `MaterialMesh` construction sites in production (`grep MaterialMesh(`):
  `orpheus/diffusion/solver.py:402` (immediately promoted —
  `DiffusionMesh.from_material_mesh(MaterialMesh(mesh, materials))`, so the
  operators downstream see a `DiffusionMesh`),
  `orpheus/homogeneous/solver.py:180` (**the meshless path**),
  `orpheus/sn/solution.py:743` (`Solution.homogenize` RETURNS a `MaterialMesh` —
  a data product; no operator is built on it inside `orpheus/`).
* Which packages even import the five operators:
  `orpheus/{diffusion, homogeneous, numerics, sn, sn/operators,
  transport/{fields, mesh, operators}}`.
  **`orpheus/cp/`, `orpheus/moc/`, `orpheus/mc/`, `orpheus/kinetics/`,
  `orpheus/fuel/`, `orpheus/thermal_hydraulics/` and `orpheus/derivations/`
  import NONE of them** — zero hits for all five symbols.

So "the meshless / model-generic path" is a single file:
`orpheus/homogeneous/solver.py`, three call sites (`:143`, `:146`, `:193`).

## 2.5 Production `basis_shape=` sites — the latent consumers

`orpheus/numerics/operator.py:474-500` `_resolve_basis_shape` is the single
resolution rule: **explicit `basis_shape` wins; else `op.domain.shape`; else
`ValueError`** whose message literally says *"Either construct the operator with a
space, or pass an explicit basis_shape="*. So every explicit `basis_shape=` in
production is a domain-less operator.

Production sites (exhaustive; `orpheus/numerics/*` definitions excluded):

| site | argument | operator | in P1 scope? |
|---|---|---|---|
| `orpheus/homogeneous/solver.py:194` | `basis_shape=(ng, 1)` | `MatrixInverseOperator(loss, …)` where `loss = C − K_iso` (the three bucket-(b) leaves) | **YES — deletable by P1** |
| `orpheus/homogeneous/solver.py:202` | `basis_shape=(ng, 1)` | `K.as_matrix(...)`, `K = MatrixInverse @ FissionOperator` | **YES — deletable by P1** |
| `orpheus/diffusion/operators.py:648` | `basis_shape=(ng,)` | `law.as_matrix(...)` where `law ∈ {ZeroOperator, IdentityOperator, ScaledOperator(IdentityOperator)}` from `DiffusionBoundaryRealizer.realize` (`boundary_realizer.py:142-157`, `:195`) | **NO** — space-polymorphic numerics primitives (deliverable-4 category (ii)); P1 does not touch them |

The code AGREES with the plan's claim, in its own words —
`orpheus/homogeneous/solver.py:136-141`:

> ``basis_shape=(ng, 1)`` remains the operators' group-leading ``(ng, *spatial)``
> bare-ndarray contract on the meshless single cell — **passed explicitly by
> consumers because the meshless operators carry no** ``domain`` **space to derive
> it from.**

Test-side latent consumers of the same argument (would also become deletable):
`tests/homogeneous/test_homogeneous.py:140,284,286,360,362,363,365`,
`tests/numerics/test_matrix_inverse_operator.py` (several).

---

# (3) The `None`-skip consumption sites

The house convention is stated at the base itself,
`orpheus/numerics/operator.py:573-592`:

> *"When either operand of a composition has ``None`` for ``domain`` or
> ``codomain``, the composability check is skipped — preserving the legacy
> duck-typed behaviour for code paths that do not track spaces."*

Exhaustive inventory (`grep` for `domain is (not )?None`,
`getattr(…, "domain"/"codomain", None)` over `orpheus/`). **Every single one is
in `orpheus/numerics/`** — no transport/sn/diffusion file consumes a `None`
domain. 6 sites CHANGE BEHAVIOUR; 10 merely PROPAGATE the `None`.

## 3.A Behaviour-changing (the guards P1 turns on)

| # | site | what the `None` branch does |
|---|---|---|
| A1 | `orpheus/numerics/operator.py:1348-1363` (`OperatorSum.__init__`) | **skips** the equal-domains AND equal-codomains guard when EITHER operand is `None`-spaced — a mismatched sum builds silently. *(known)* |
| A2 | `orpheus/numerics/operator.py:1556-1564` (`OperatorProduct.__init__`) | **skips** the `A.domain == B.codomain` guard for `A @ B` when either side is `None`. |
| A3 | `orpheus/numerics/operator.py:1223` (`_AdjointOperator.apply`) | `z = inner_codomain.apply_metric(y) if inner_codomain is not None else y` — with `None` the adjoint applies **no G metric**, silently degrading `A† = G⁻¹AᵀG` to the bare Euclidean `Aᵀ`. *(known)* |
| A4 | `orpheus/numerics/operator.py:1225-1226` (`_AdjointOperator.apply`) | the mirror half: `if inner_domain is not None: result = inner_domain.apply_inverse_metric(result)` — `None` skips `G⁻¹`. *(known)* |
| A5 | `orpheus/numerics/operator.py:489-500` (`_resolve_basis_shape`) | `domain is None` ⟹ **`ValueError`** ("the operator carries no domain FunctionSpace…"). The one site where `None` is LOUD rather than silent — and the reason §2.5's explicit `basis_shape=` arguments exist. |
| A6 | `orpheus/numerics/coupled_system.py:711-733` (`CoupledOperator.__init__`) | per-block placement typing: `block_domain is not None and block_domain != domain_members[j]` (and the codomain twin at `:721-724`) — a `None`-spaced block can be placed at the **wrong grid position** with no error. Comment at `:705-707` names it "the house convention for legacy space-less operators". |

## 3.B Propagation-only (`None` in ⟹ `None` out, no behaviour change)

| site | note |
|---|---|
| `orpheus/numerics/operator.py:998-1001` (`LinearOperator.__repr__`) | prints `'?'` instead of the space name — cosmetic only. |
| `orpheus/numerics/operator.py:1198`, `:1202` (`_AdjointOperator.domain/.codomain`) | swap-and-forward. |
| `orpheus/numerics/operator.py:1391-1392`, `:1396-1397` (`OperatorSum.domain/.codomain`) | **`a` else `b` fallback** — an interesting asymmetry: a sum of (`None`-spaced, spaced) reports the SPACED one, so `C − K_iso` in the homogeneous path would inherit a space from either summand if one had it. |
| `orpheus/numerics/operator.py:1581`, `:1586` (`OperatorProduct.domain/.codomain`) | `B.domain` / `A.codomain`. |
| `orpheus/numerics/operator.py:1752`, `:1756` (`ScaledOperator.domain/.codomain`) | passthrough. |
| `orpheus/numerics/flat_operator.py:112`, `:116` (`FlatOperator.domain/.codomain`) | passthrough to inner. |

## 3.C The same-`None`, different-spelling family (`space` / `full_field_space`)

These read the SAME field P1 targets, but through the operator's own attribute
rather than through `domain`. They are the assembly-mode refusals — P1 flips them
from "sometimes refuses" to "always emits":

| site | behaviour on `None` |
|---|---|
| `orpheus/transport/operators/multiplication_operator.py:253-256` (`is_assemblable`) | `False` — refuses to assemble. |
| `orpheus/transport/operators/multiplication_operator.py:278-283` (`assemble`) | raises `MissingAssembly` ("this multiplier is space-anonymous"). |
| `orpheus/transport/operators/isotropic_scattering.py:152-158` (`_iso_is_assemblable`) | `False` for BOTH iso operators. |
| `orpheus/transport/operators/isotropic_scattering.py:191-198` (`_assemble_iso_energy_operator`) | raises `MissingAssembly`. |

(`FullFieldSpace.interior_space` / `.trace_space` are independently `Optional` —
`orpheus/numerics/spaces/full_field_space.py:191`, `:194`, with the two-block
requirement enforced at `:259`. So "has a domain" and "has a *block-bearing*
domain" are two different predicates; P1 fixes only the first.)

---

# (4) The base-class population

`orpheus/numerics/operator.py:573-584` — `LinearOperator.domain` returns `None`
with the docstring naming the convention ("Operators that pre-date Issue 9.6 …
the composability check is skipped").

## 4.0 The corrected census (AST transitive closure over `orpheus/`)

| | count |
|---|---|
| `LinearOperator` descendants (transitive, `orpheus/` only) | **55** |
| … that define their own `def domain` | **27** |
| … that do NOT define `def domain` | **28** |
| of those 28: false positive (subclasses scipy's `spla.LinearOperator`, not ours) | 1 |
| of those 28: **inherit a WORKING derived `domain` from a composer/wrapper base** | **8** |
| of those 28: **actually resolve to the `None`-returning base** | **19** |

(Your "44 / ~24" figures were direct-subclass counts; the transitive closure adds
the `OperatorSum`/`OperatorProduct`/`InverseWrapMixin` descendants. Verified by a
runtime MRO probe — `next(k for k in cls.__mro__ if "domain" in k.__dict__)`.)

### The 8 that DERIVE a domain and would NOT break

These already have honest arrows and are silent beneficiaries of P1 (their
derived domain becomes non-`None` the moment their operands do):

| class | site | derives `domain` from |
|---|---|---|
| `InverseOperator` | `orpheus/numerics/operator.py:2103` | `InverseWrapMixin` (`:2056-2062`) — `inner.codomain`/`inner.domain` **swap** |
| `MatrixInverseOperator` | `orpheus/numerics/matrix_inverse_operator.py:95` | `InverseWrapMixin` swap |
| `GreenOperator` | `orpheus/numerics/green_operator.py:192` | `InverseWrapMixin` swap |
| `SweepOperator` | `orpheus/sn/operators/sweep_operator.py:83` | `InverseWrapMixin` swap |
| `CoupledSubstitutionOperator` | `orpheus/numerics/coupled_system.py:1212` | `InverseWrapMixin` swap |
| `StreamingCollisionOperator` | `orpheus/sn/operators/streaming.py:445` | `OperatorSum` a-else-b |
| `ScheduledInvertibleOperator` | `orpheus/sn/operators/scheduled_invertible.py:87` | `OperatorSum` a-else-b |
| `WindowedSweep` | `orpheus/sn/operators/windowing.py:163` | `OperatorProduct` (`B.domain`/`A.codomain`) |

False positive excluded: `_CarrierMatvecOperator`
(`orpheus/numerics/iteration.py:325`) subclasses **scipy**'s
`spla.LinearOperator`, not ORPHEUS's.

## 4.1 THE DELIVERABLE — the 19 that land on the `None` base, three-way split

### (i) Protocol / ABC marker — not really an operator (6)

Making the base `domain` abstract is HARMLESS here; these are never instantiated.

| class | site | one-line |
|---|---|---|
| `SupportsInverse` | `orpheus/numerics/operator.py:1042` | `Protocol` — the static narrowing bridge for `invertible()`; a capability marker, no instances. |
| `SupportsAdjoint` | `orpheus/numerics/operator.py:1057` | `Protocol` — same, for `adjointable()`. |
| `SupportsAssembly` | `orpheus/numerics/operator.py:1072` | `Protocol` — same, for `assemblable()`. |
| `IntegralKernelOperator` | `orpheus/transport/operators/integral_kernel_operator.py:164` | `LinearOperator[V], Protocol[V]` — the kernel-category interface; conformers supply their own arrows. |
| `AnalysisOperator` | `orpheus/numerics/projection.py:88` | `ABC` — the frame's `M : V → W` role; its docstring EXPLICITLY says "the V/W distinction lives in the runtime `domain`/`codomain`", and its concrete faces (`_FrameAnalysis:412`, `BulkAnalysisOperator:71`) DO override. |
| `ReconstructionOperator` | `orpheus/numerics/projection.py:136` | `ABC` — the `R : W → V` sibling; concrete face `_FrameReconstruction:446` overrides. |

### (ii) Space-polymorphic numerics primitive — the space IS whatever the consumer uses (9)

These are the honest hard cases: the same instance is legitimately reused across
spaces. `IdentityOperator()` / `ZeroOperator()` are constructed with **no
arguments at all** and are the structure-keyed collapse targets of the boundary
realizers (`orpheus/diffusion/boundary_realizer.py:153-157`,
`orpheus/sn/boundary/realizer.py`). Forcing a domain here means either a
mandatory-arg break at every realizer, or an opt-out.

| class | site | one-line |
|---|---|---|
| `IdentityOperator` | `orpheus/numerics/operator.py:1834` | `I` on any space; nullary ctor; realizer collapse target for albedo `𝒜=1`. |
| `ZeroOperator` | `orpheus/numerics/operator.py:1862` | `0` on any space; nullary ctor; realizer collapse target for `𝒜=0` (vacuum). |
| `DiagonalOperator` | `orpheus/numerics/operator.py:2655` | the N-D broadcast engine `MultiplicationOperator` delegates to; carries a coefficient array, not a space. |
| `PermutationOperator` | `orpheus/numerics/operator.py:2171` | index permutation on any conforming axis. |
| `TensorProductOperator` | `orpheus/numerics/operator.py:2420` | `A ⊗ B` on the factor spaces; could derive, does not. |
| `SumOfTensorProductsOperator` | `orpheus/numerics/operator.py:2561` | `Σ Aᵢ ⊗ Bᵢ` — the §15.2 white-BC form. |
| `RankOneOperator` | `orpheus/numerics/operator.py:2898` | `u ⊗ vᵀ` dyad. |
| `IncomingOrdinateMaskTensor` | `orpheus/numerics/operator.py:2284` | the sparse inflow-ordinate mask (vacuum trace); self-adjoint idempotent projector — acts on whatever ordinate axis it is handed. |
| `PeriodicWrapOperator` | `orpheus/numerics/operator.py:2378` | angular-identity periodic pushforward; docstring says the type is *reserved* for a future spatial extension. |

### (iii) Real domain-specific operator that simply never got tagged (4)

**These are the true finds.** Each has exactly one meaningful function space and
carries mesh/quadrature/face data — nothing generic about them.

| class | site | one-line |
|---|---|---|
| `AngularAverageOperator` | `orpheus/sn/boundary/angular.py:36` | `G_diff`, the cosine-weighted Lambertian hemisphere average of the white BC — stores per-quadrature cosine-weighted mask values + a normalization scalar. Its space is the face's angular trace; untagged. |
| `IncomingSourceOperator` | `orpheus/sn/boundary/angular.py:228` | the affine `q` of `γ₋ψ = RGγ₊ψ + q` (prescribed inflow) — stores an `InflowSourceSpec` **and the face's inflow-ordinate indices**. Face-bound; untagged. |
| `DSACorrection` | `orpheus/sn/acceleration/dsa.py:453` | `LinearOperator["FullField","FullField"]` — the within-group DSA correction `Δψ ↦ P A_low⁻¹ G R Δψ`. Statically typed on `FullField` in the class parameters but with **no runtime `domain`**; it is exactly the SN `full_field_space`. |
| `_BoundBoundaryOperator` | `orpheus/geometry/boundary/_bound_compat.py:77` | the 1-arg BC passthrough shim built at `orpheus/sn/mesh/augmented_mesh.py:427`. It delegates `is_invertible`, `is_adjointable`, `inverse`, `block_role`, `apply`, `apply_transpose`, `__eq__`, `__hash__` to `inner` — but **NOT `domain`/`codomain`**. A wrapper that drops the space it wraps (see §5, surprise S-2). |

## 3.D LATE ADDITION — a 7th behaviour-changing site (found on the final sweep)

| # | site | what the `None` branch does |
|---|---|---|
| A7 | `orpheus/numerics/matrix_inverse_operator.py:181-208` (`_returned_as`), reached from `apply` at `:225` (`self.inner.domain`) and `solve` at `:249` (`self.inner.codomain`) | `if space is None or not callable(zeros_factory): raise ValueError(…)` — a **TYPED (ravellable) operand cannot be re-typed** when the inner operator has no domain: the backsolve result has no zero exemplar to mint from. An ndarray operand silently takes the `reshape(self._basis_shape)` legacy path instead. **This is why the homogeneous solver must feed bare `(ng,1)` arrays, not typed carriers, through `MatrixInverseOperator`.** LOUD, like A5. |

So the corrected count is **7 behaviour-changing** + 10 propagation-only + 4
assembly-refusal (§3.C). Both LOUD sites (A5 `_resolve_basis_shape`, A7
`_returned_as`) are on the SAME consumer — `MatrixInverseOperator` — and both are
hit by the SAME production caller, `orpheus/homogeneous/solver.py:194`.

---

# (5) What surprised me — for the strategy discussion

## S-1. Production bucket (b) is SIX sites in TWO files, and one of them is free

The whole production blast radius of "C/S/F must state their arrow" is
`orpheus/homogeneous/solver.py` (3 sites) + one line in
`orpheus/transport/operators/scattering.py`. Everything else in `orpheus/`
already threads a space. That is a much smaller carve than the 220/91/69
grep-mention counts suggest.

And **`scattering.py:709` is free** — `ScatteringOperator.isotropic_kernel` builds
`IsotropicScattering(self.mat_xs) + IsotropicN2N(self.mat_xs)` while
`self.full_field_space` sits right there in scope. The docstring at `:701-702`
even declares the choice: *"Space-anonymous (`space=None`: no composition-guard
space at this producer-side use)."* But that object is NOT purely producer-side —
`orpheus/sn/coupled_system.py:498` passes it into
`RadialCharacteristicEmission(sn_mesh, S.isotropic_kernel)`, i.e. **an anonymous
sum ends up inside a `CoupledOperator` grid**, where the §3.A6 per-block
placement guard skips it. Threading `self.full_field_space` there is a
one-line change with no meshless obstacle.

## S-2. `_BoundBoundaryOperator` is a passthrough that drops the space

`orpheus/geometry/boundary/_bound_compat.py:77` forwards `is_invertible`,
`is_adjointable`, `inverse`, `block_role`, `apply`, `apply_transpose`, `__eq__`
and `__hash__` to `inner` — **but not `domain`/`codomain`.** Every other wrapper
family in the codebase does forward them (`ScaledOperator:1752`,
`_AdjointOperator:1198`, `FlatOperator:112`, `InverseWrapMixin:2056`). So every
realized SN face law reached through `sn_mesh.bc[...]`
(`orpheus/sn/mesh/augmented_mesh.py:427`) is space-blind *even when its inner
operator is not*. That looks like an oversight rather than a decision, and it is
a 4-line fix independent of P1's scope.

## S-3. L and B are ANNOTATION-only changes; the runtime is already monomorphic

`orpheus/sn/operators/streaming.py:279` and `orpheus/sn/operators/boundary.py:208`
both annotate `-> Optional["FunctionSpace"]` while unconditionally returning
`self.sn_mesh.full_field_space`. There is no runtime path on which they are
`None`. So 2 of the 5 leaves cost a type annotation and nothing else — the real
work is C/S/F, and within those, only the meshless path. `LegendreMomentScattering`
(`scattering.py:279`) and `N2NMomentOperator` (`:343`) are the worked precedent
for how a non-Optional domain is spelled here: **build the space on demand from
owned data** (`SphericalHarmonicSpace.from_L(self.L)`), not store it.

## S-4. The meshless space is already derivable — nothing is missing

`MaterialMesh` has no `FunctionSpace`, but it has everything one is made of:
`ng` (`material_mesh.py:330`), `spatial_shape` (`:428`, `(1,)` for the meshless
carrier), `volumes` (`:365`, `[1.0]`). `(ng, *spatial_shape) == (ng, 1)` is
**literally the `basis_shape` the homogeneous solver hand-writes** at
`homogeneous/solver.py:194,202`. `DiffusionMesh.full_field_space`
(`diffusion/augmented_mesh.py:361-367`) builds its bulk block from exactly those
three members. So the C/S/F domain on the meshless path is a derivation, not new
data — the only genuinely absent piece is the **trace** block, and
`FullFieldSpace.interior_space`/`.trace_space` are both `Optional`
(`full_field_space.py:191,194`). Whether an infinite medium *should* have a trace
block at all is the actual open design question (it has no boundary — that IS the
physics).

## S-5. The 28 "no `def domain`" classes are really 19 — 8 already derive one

The MRO probe (§4.0) shows `InverseWrapMixin` (5 classes),
`OperatorSum` (2) and `OperatorProduct` (1) all supply a *derived* domain. So
"make the base abstract" breaks 19, not 28 — and 6 of those 19 are Protocols/ABCs
that are never instantiated. The genuine decision surface is **9 space-polymorphic
primitives + 4 untagged domain-specific operators**.

## S-6. The space-polymorphic 9 are load-bearing in the BC realizers

`IdentityOperator()` and `ZeroOperator()` are constructed **nullary** as the
structure-keyed collapse targets of BOTH realizers —
`orpheus/diffusion/boundary_realizer.py:153-157`,
`orpheus/sn/boundary/realizer.py:262,280,301,309,311,320`. A mandatory-domain
base would break every one of those lines. This is the one place where "monomorphic
domains" collides with a deliberate design (`𝒜=0 ⟹ ZeroOperator()`,
`𝒜=1 ⟹ IdentityOperator()`), so the P1 rule almost certainly has to be scoped to
the transport/SN *leaves* rather than pushed onto `LinearOperator` itself — which
is exactly how the phase is named.

## S-7. The gate for this phase is already committed, red, on this branch

`tests/sn/architecture/test_monomorphic_leaves.py` (landed `cd7b17cd`, refined
`97b14e45`) already encodes G1.1/G1.3/G1.4/G1.5 with **13 `xfail` markers**, and
its own measured-state table names the same findings this audit reaches
independently:

* *"G1.1 value — GREEN on the SN ladder, RED off it (**R1**): the model-generic
  construction (`homogeneous/solver.py:193`) reports `domain is None`"* — matches
  §1.1 bucket (b) exactly.
* *"G1.1 ann. — RED (R1's static face): every leaf annotates
  `FunctionSpace | None` / `Optional[FunctionSpace]`"* — matches S-3.
* *"G1.5 — RED for C/S/F (**R1**/**R2**): all three construct happily with no
  space, and the `.H` they then build is **bit-identical** to the bare Euclidean
  transpose."* — that is §3.A3/A4 measured, not argued.
* Also measured there and worth carrying into the strategy: **`C` and `B` are
  metric-BLIND by structure**, not by fixture accident (`C` is diagonal so
  `G⁻¹CG = C` for any diagonal metric; `B` specular is a signed permutation that
  commutes with the trace metric). So the A3/A4 silent-degradation risk is real
  for `L`/`S`/`F` and *structurally unobservable* for `C`/`B`.

The practical consequence: the P1 change-set is a small production carve (6 sites,
2 files) with an oversized **test** migration (53 bucket-(b) test sites + the 114
`from_mesh`-implicit ones if the `getattr` fallback is retired), and the gate that
proves it is already in the tree waiting to XPASS.

---

# CORRIGENDUM to §1.2 — the exact test-side split

The §1.2 table's implicit/bare columns were computed before the mesh-argument
resolution. **Corrected, exact numbers** (`from_mesh` counted as implicit-(a) only
when the mesh argument is an SN-named local — `sn`/`sn_mesh`/`sn1`/`sn2`/`sph`):

| tree | (a) explicit kwarg | (a) implicit `from_mesh(σ, SNMesh)` | (b) genuinely `None` | total |
|---|---|---|---|---|
| `orpheus/` | 11 | 0 | **6** | 17 |
| `tests/` | 20 | 111 | **56** | 187 |

Per operator in `tests/`: MultiplicationOperator 12 / 111 / **12**;
ScatteringOperator 2 / 0 / **11**; FissionOperator 1 / 0 / **11**;
IsotropicScattering 4 / 0 / **12**; IsotropicN2N 1 / 0 / **10**.

`MultiplicationOperator` bucket-(b) test sites, exhaustive (12):

```
tests/diffusion/test_operators.py:625                            ctor
tests/diffusion/test_operators.py:641                            from_mesh(…, mesh)
tests/diffusion/test_operators.py:853                            ctor
tests/numerics/test_matrix_inverse_operator.py:206               from_mesh(…, mat_xs.mesh)
tests/numerics/test_matrix_inverse_operator.py:265               from_mesh(…, mat_xs.mesh)
tests/sn/architecture/test_monomorphic_leaves.py:692             from_mesh(…, mat_xs.mesh)
tests/sn/architecture/test_monomorphic_leaves.py:1087            ctor
tests/sn/sweep/cartesian_2d/test_2d_l2_matvec_correctness.py:211 from_mesh(…, mesh)
tests/sn/sweep/cartesian_2d/test_2d_l2_matvec_correctness.py:312 from_mesh(…, mesh)
tests/transport/test_functional_category.py:85                   ctor
tests/transport/test_integral_kernel_category.py:166             ctor
tests/transport/test_multiplication_operator.py:156              ctor
```

(`tests/diffusion/test_operators.py:641`'s `mesh` fixture is a `DiffusionMesh`,
which DOES carry `full_field_space` — so it is implicit-(a) in practice and
listed here only because the argument name is not SN-shaped. Confirm the fixture
before editing it. The two `test_2d_l2_matvec_correctness.py` sites' `mesh` is an
`SNMesh` under a non-SN name — same caveat, opposite direction.)

Census stability (lessons L-007): re-run immediately before writing this report —
204 call sites both times, `git status` clean under `orpheus/` and `tests/`,
HEAD unchanged at `3d1a9b39`.
