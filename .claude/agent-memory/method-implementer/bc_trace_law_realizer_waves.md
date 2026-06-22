---
name: bc-trace-law-realizer-waves
description: Consolidated record of the landed Boundary Operator refactor (Waves 0/1/5/7/8/9/10/11 + C188) — BoundaryTraceLaw + SNBoundaryRealizer + Wave-0 operator primitives + MixedBoundaryOperator retirement + realize_recursively walker. All merged to origin/main. Folds 9 per-wave diary entries into one durable architecture + lessons note.
metadata:
  type: project
---

# Boundary Operator refactor — BoundaryTraceLaw + SNBoundaryRealizer (landed)

The 12-wave Boundary Operator refactor (plan was `transient-giggling-cake.md`)
moved SN boundary conditions from a stringly-keyed `psi_bc: dict` + 2-arg
`bc.apply(psi, quad)` legacy ABC to a three-layer architecture:
**trace structure (Wave 2) → physical law (`BoundaryTraceLaw`, Wave 3/7) →
method realisation (`SNBoundaryRealizer`, Wave 5)**. All work is **LANDED on
origin/main**. This note consolidates the durable architecture + lessons from
waves 0, 1, 5, 7, 8, 9, 10, 11 and the C188 curvilinear-realizer follow-up
(the per-wave diary entries are retired in favour of this). Waves 2/2a/2b/2c/3
have their own committed notes ([[wave_2_trace_spaces]], [[wave_3_bc_trace_law_abc]]).

## The landed architecture (verify against code)

- `orpheus/numerics/operator.py` — Wave-0 generic primitives:
  `PermutationOperator` (reflective ordinate swap via `quad.reflection_index`,
  `inverse_perm = argsort`, `is_involution` flag), `IncomingOrdinateMaskTensor`
  (sparse inflow-only mask — distinct from `ZeroOperator` which zeroes
  EVERYTHING; this preserves the outflow trace per §16A.10), `PeriodicWrapOperator`.
  Capability `{CAP_APPLY, CAP_APPLY_TRANSPOSE}`; plain classes (match
  `IdentityOperator`/`ZeroOperator` style, NOT `@dataclass`).
- `orpheus/sn/angular_operator.py` — `AngularAverageOperator` (cosine-weighted
  Lambertian average lifted verbatim from `WhiteBoundaryOperator.apply`;
  `{CAP_APPLY}` only — the BC is self-adjoint under the COSINE-weighted inner
  product, NOT the unweighted one) + `IncomingSourceOperator` (Wave 7, for
  `PrescribedInflow`: ignores `psi_out`, returns `source.evaluate(probe_trace)`).
- `orpheus/geometry/boundary/_realizer.py` — `BoundaryRealizer` Protocol +
  `BoundaryRealizerRegistry` (stand-alone string-keyed registry, NOT
  `RegistryMixin`; `@register("X")` raises on collision at import time).
- `orpheus/sn/boundary_realizer.py` + `orpheus/sn/method_space.py` —
  `SNBoundaryRealizer` (isinstance dispatch on `BoundaryTraceLaw` subclasses →
  Wave-0/1 primitive) + `SNMethodSpace` (mesh + face + traces; `.minimal(quad)`
  / `.for_face(...)`).
- `orpheus/sn/boundary_realize.py` — `realize_recursively(op, method_space)`
  tree-walker for the 5 Wave-0 composers + `BoundaryTraceLaw` leaf (Wave 11
  replacement for `MixedBoundaryOperator`).
- `orpheus/geometry/boundary/` — `BoundaryTraceLaw` is the single registry root;
  concretes `VacuumInflow` / `ReflectiveBoundary` / `WhiteBoundary` /
  `PeriodicBoundary` / `AlbedoBoundary` / `PrescribedInflow`. `mixed.py`
  DELETED (Wave 11). MoC/MC/CP/diffusion ship stub realizers (auto-registered
  via `from . import boundary_realizer` in each `__init__`).

The realisation map (legacy class → realized primitive at α=1):
`VacuumInflow → IncomingOrdinateMaskTensor`; `ReflectiveBoundary →
PermutationOperator(reflection_index)`; `WhiteBoundary →
AngularAverageOperator.from_quadrature`; `AlbedoBoundary(0/1/α) →
Zero/Identity/ScaledOperator(α,Identity)`; `PeriodicBoundary →
PeriodicWrapOperator`.

## Durable lesson 1 — the α=1.0 fast-path is what makes the rename bit-identity-preserving

For specular and white, `SNBoundaryRealizer.realize` returns the BARE primitive
(no `ScaledOperator` wrap) when α=1.0. This is load-bearing: without it, the
legacy "perfect reflection / perfect white" case would shift by one ULP under
the realizer, and the Wave-7 isinstance-update rename would NOT be a
bit-identity-preserving change. Non-vacuum (α≠1) realizations are pinned at
`nulp ≤ 4` vs legacy. **Generalises**: whenever you replace a code path with a
composed-operator path, the identity/scalar-one case must return the
unwrapped primitive so the common case stays bit-exact.

## Durable lesson 2 — the §16A.5 vacuum semantic correction (zeros INFLOW only, not ALL)

Legacy `VacuumBoundaryOperator.apply` zeroed ALL ordinates; the new
`IncomingOrdinateMaskTensor` zeroes ONLY inflow ordinates and preserves the
outflow trace (§16A.10 trace-correct representation). This is an INTENTIONAL
semantic change, gated by routing through `SNBoundaryRealizer.realize`. The
risk was "does any consumer read the outflow rows of a vacuum BC result?" A
Wave-9 grep audit of all 13 2-arg `bc.apply` call sites in `sweep.py` +
`operator.py` confirmed **every consumer reads INFLOW rows only** (each has an
explicit `if mu_n < 0` / `if mu_x[n] > 1e-15` inflow gate). C188 then verified
this empirically end-to-end through the production curvilinear sweep (spherical
26/26 + cylindrical 25/25 green). **Durable rule**: a semantic correction that
narrows what a BC zeroes is safe IFF a call-site audit proves no consumer reads
the now-preserved rows — and the audit must be backed by an end-to-end test, not
just the grep.

## Durable lesson 3 — compatibility-shim semantics extend beyond the apply signature

The Wave-8 `_BoundBoundaryOperator` transitional shim (2-arg → 1-arg) was
sketched as "swallow the extra `quad` arg". In practice the legacy BC object's
API ALSO included `__eq__` against strings (`sn_mesh.bc_xmin == "vacuum"`), and
4 tests relied on it. The clean sketch missed this. **Fix**: the shim carries a
`kind` tag (`kind=bc.kind`) and `__eq__` returns `self.kind == other` for string
comparands, with `__hash__` falling back to `id(self)` (Python auto-disables
`__hash__` when `__eq__` is defined). **Durable rule**: when wrapping a legacy
object whose API includes `__eq__`-against-strings (or any non-`apply` surface),
the wrapper must replicate that surface OR you migrate every assertion site —
the kind-tag replication is the cheaper path.

## Durable lesson 4 — curvilinear needs an escape hatch in any trace-space-dependent design (then C188 closes it)

Wave 2's `InflowTraceSpace.from_mesh_and_quadrature` raised `NotImplementedError`
on curvilinear meshes. Wave 8 therefore needed an explicit `if _inflow_trace is
None: return bare law instance` branch in `_resolve_one` (curvilinear kept the
legacy 2-arg path; bit-identical). Without it, every curvilinear `SNMesh`
construction would crash — exactly the broad blast radius the refactor was
avoiding. **C188 then closed the hatch**: it added curvilinear `InflowTraceSpace`
support (C188.1/.2) and DELETED the bypass (C188.3), so every supported mesh
routes uniformly through `SNBoundaryRealizer`. After C188 the bound-quadrature
dual-mode backing (Wave 9) has NO production consumer (every production-issued
`bc_*` shim has `_quadrature is None`), unblocking Issue #176 Phase 2 shim
removal. **Durable rule**: a deferral (`NotImplementedError` for the hard
geometry) MUST be paired with an explicit `is None` escape-hatch branch at every
downstream consumer; the unification commit later deletes both.

## Durable lesson 5 — sum-of-stubs over a shared ABC; collision-loud registration

Per "unify after two instances": no `MethodBoundaryRealizerBase` ABC was
factored out across MoC/MC/CP/diffusion — each ships its own ~30-line stub
class. The shared abstraction emerges only when ≥2 methods ship FUNCTIONAL
realizers. Similarly, `realize_recursively` lives in `orpheus.sn.boundary_realize`
(SN-specific); cross-method generalisation to `orpheus.geometry.boundary.realize`
(via `BoundaryRealizerRegistry.get(method_name)`) is deferred until a second
functional realizer ships. Registration is collision-loud (`@register("X")`
raises at import time, no silent override). **The 1-D y-face placeholder
subtlety (C188)**: `bc_ymin`/`bc_ymax` on a 1-D `Mesh1D` route through
`SNMethodSpace.minimal(quad)` (NOT `.for_face`, which would raise — 1-D trace
spaces have only `("left","right")` faces). The realized op is an identity
`PermutationOperator` because `GaussLegendre1D.reflection_index("y")` returns
`arange(N)` (since `mu_y == 0` on 1-D). Production consumers only invoke
`bc_ymin.apply` when `curvature is None` and the inner `if mu_y[n] > 1e-15`
filter is false for every GL1D ordinate — so the no-op is observably correct.

## Durable lesson 6 — the `PrescribedInflow` field-vs-property collision (Wave 7)

`BoundaryTraceLaw.source` is a `@property`. Declaring `source: BoundarySource`
as a dataclass field on `PrescribedInflow` collides with the inherited
descriptor (no setter ⇒ `AttributeError` from dataclass `__init__`). **Resolved**
by naming the field `_source` + a custom `__init__(source=...)` forwarding via
`object.__setattr__` + overriding the inherited `source` property to read
`self._source` (same idiom `MixedBoundaryOperator` used). Public UX unchanged
(`PrescribedInflow(source=ConstantInflowSource(2.5))`). **Durable rule**: a
dataclass field that shadows an inherited read-only `@property` fails at
construction; use a `_private` field + custom `__init__` + property override.

## Retirement accounting (Wave 11)

`MixedBoundaryOperator` (class + file + tests + re-exports) DELETED; replaced by
`realize_recursively` walking the Wave-0 composer algebra
(`ScaledOperator`/`OperatorSum`/`OperatorProduct`/`TensorProductOperator`/
`SumOfTensorProductsOperator`). The Marshak example `0.3*spec + 0.7*white`
becomes Wave-0 algebra. **Snapshot bit-identity**: the deleted realizer's
internal mixed-BC path already composed via `OperatorSum`-of-`ScaledOperator`,
so the regenerated BC-equivalence snapshot was BIT-IDENTICAL to pre-Wave-11
(`np.array_equal` confirmed). Net +457 LoC but the implementation surface
NARROWED (one generic walker for all 5 composers vs a single-purpose class) and
the V&V surface WIDENED (structural tests per composer type vs only the
linear-combination path). No new ERR-NNN across the whole refactor (structural;
the realizer-path equivalence tests on composed `apply` output catch any numeric
drift).

## No Sphinx narrative per wave

Every wave deferred Sphinx (mechanical refactor; no new equation labels). The
unified BC theory page covering Waves 0-12 shipped once at Wave 12 (commit
`80112e6`) — incremental per-wave narrative would have fragmented the
architectural story. **Durable pattern**: for a multi-wave structural refactor,
defer the rich narrative to ONE archivist dispatch at the end, not per-commit.
