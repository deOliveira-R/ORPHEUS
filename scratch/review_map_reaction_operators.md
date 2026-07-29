# Structural map — the two reaction-operator facades (`orpheus/transport/operators/`)

**Grounding**: branch `main`, HEAD `8654d348`, Nexus graph fresh. Files read IN FULL:
`fission.py` (674 L), `scattering.py` (1326 L), `isotropic_scattering.py` (438 L),
`multiplication_operator.py` (530 L), `frames/harmonic_frame.py` (208 L),
`numerics/frame.py` (479 L). Every claim below is verified against the code BODY, not a
docstring (L33). Docstring claims that do NOT match the body are flagged `⚠ DOC-CLAIM`.

`file:line` current at this read; re-derive via Nexus if the tree moves.

---

## 0. TL;DR — the two facades are NOT symmetric

| Axis | `FissionOperator` (`fission.py:145`) | `ScatteringOperator` (`scattering.py:356`) |
|---|---|---|
| dataclass | `@dataclass` (mutable, `eq` default) | `@dataclass` (mutable) |
| stored state | **1 field** + 1 optional space | **3 fields** + 1 optional space + 1 lazy cache |
| quadrature | **NONE** — never touched at construction | **REQUIRED positional field** `quadrature` (`scattering.py:386`) |
| quadrature at apply | reached ONLY via `bulk.mesh.quad` inside the composite-transpose angular arm (`fission.py:436`) | reached in `weights`, `Y`, `frame`, `_assemble_per_ordinate_source`, `apply_transpose` |
| `apply` arms | 3 (`FullField`, `ScalarFlux`, `ndarray`) | 4 (`FullField`, `ScalarFlux`, `AngularFlux`, `HarmonicMomentFlux`) |
| distinct in→out arrows | **5** (3 forward + 2 transpose), *7* counting the composite sub-arms | **6** (4 forward + 2 transpose), *8* counting the composite sub-arm |
| sibling leaf operators | `RankOneOperator` via `outer` (numerics) | `LegendreMomentScattering` + `N2NMomentOperator`, both defined IN `scattering.py` |
| `.kernel` type | `TensorProductOperator` (concrete) | `LinearOperator` (erased `OperatorProduct`), **RAISES** at `L=0` |
| frame | none | `HarmonicFrame` built INSIDE (`cached_property`, `scattering.py:489`) |

**Verdict on the plan's claim** ("`ScatteringOperator` is method-specific in fact because it
needs a quadrature"): **CONFIRMED, and stronger than stated.** `quadrature: "Quadrature"` is a
**mandatory positional dataclass field** with no default (`scattering.py:386`) and no
`None` branch anywhere in the class — you cannot construct a `ScatteringOperator` without a
discrete angular quadrature, even for the pure-P0 scalar path. See §2b for the precise
list of quadrature-dependent members and §2c for the counter-nuance (the *cross-section*
half is genuinely method-neutral and is ALREADY factored out into
`isotropic_scattering.py`).

---

## 1. The exact `apply` dispatch surface

### 1a. Mechanism — identical in both facades

Both use **`functools.singledispatchmethod` on a PRIVATE `_apply_impl`, with the public
`apply` name runtime-ALIASED to it**, guarded by a `TYPE_CHECKING` split that supplies
`@overload` stubs:

- `fission.py:467` `@singledispatchmethod def _apply_impl(self, phi)` (fallback raises
  `TypeError`, `fission.py:510-515`)
- `fission.py:658-674` — `if TYPE_CHECKING:` → 3 `@overload def apply` stubs + a
  body-less `def apply`; `else: apply = _apply_impl` (`fission.py:674`)
- `scattering.py:1093` `@singledispatchmethod def _apply_impl(self, psi)` (fallback raises,
  `scattering.py:1127-1132`)
- `scattering.py:1250-1270` — same `TYPE_CHECKING` / `else: apply = _apply_impl` shape
  (`scattering.py:1270`)

**Consequence for exploration (L-001 sharpening, confirmed live here):** the Nexus graph
attributes calls to the ALIAS, and the three/four `@_apply_impl.register` arms are all named
`_` (Nexus shows them as one node `…ScatteringOperator._` with degree 55, `scattering.py:1210`,
and `…FissionOperator._`, `fission.py:643`). `callers()` on a per-arm node is meaningless
here — the consumer inventory in §9 is grep-primary.

**Nested `isinstance` INSIDE the `FullField` arm** — both facades' composite arm is *not* a
single arrow; it re-parses `psi.interior` by field FAMILY:

- `fission.py:560` `isinstance(bulk, AngularFlux)` / `fission.py:595` `isinstance(bulk, ScalarFlux)` / `fission.py:608` raise
- `scattering.py:1153` — no isinstance; it `cast(...)`s and RE-DISPATCHES through the public
  `self.apply` on `psi.interior` (so the family parse is delegated to the singledispatch table)

That is a real structural asymmetry: **fission hand-parses the composite bulk; scattering
re-enters its own dispatcher.**

### 1b. `FissionOperator` — every (input → output) arrow

| # | Input type | Output type | Mechanism site | Body site | Math |
|---|---|---|---|---|---|
| F1 | `FullField` (and `TimedFullField` via MRO) with `AngularFlux` bulk | `FullField(interior=AngularSourceSink, boundary=AngularBoundarySourceSink.zeros_on)` | `fission.py:517` register | `fission.py:559-594` | `φ = ∫ψdΩ` (`integrate_angular`), `Fφ`, then `AngularSourceSink.from_isotropic(…, mesh)` |
| F2 | `FullField` with `ScalarFlux` bulk | `FullField(interior=ScalarSourceSink, boundary=ScalarBoundarySourceSink.zeros_on)` | same register | `fission.py:595-607` | identity reduction; bare `Fφ` |
| F3 | `ScalarFlux` | `ScalarSourceSink` | `fission.py:615` register | `fission.py:639-640` | `self.kernel.apply(phi.values)` → `ScalarSourceSink.from_mesh` |
| F4 | `np.ndarray` | `np.ndarray` | `fission.py:642` register | `fission.py:656` | `self.kernel.apply(phi_arr)` |
| F5 | anything else | `TypeError` | — | `fission.py:510` | — |
| FT1 | `FullField` w/ `AngularFlux` bulk (transpose) | `FullField(AngularSourceSink ⊕ zero AngularBoundarySourceSink)` | plain `isinstance`, `fission.py:426`+`428` | `fission.py:429-445` | `w_n ⊗ kernelᵀ(Σ_m ψ*_m / W)` |
| FT2 | `FullField` w/ `ScalarFlux` bulk (transpose) | `FullField(ScalarSourceSink ⊕ zero ScalarBoundarySourceSink)` | `fission.py:446` | `fission.py:447-457` | bare `kernelᵀ` |
| FT3 | `np.ndarray` / any `.values`-carrier (transpose) | `np.ndarray` | `fission.py:464` `getattr(phi_star, "values", phi_star)` | `fission.py:465` | `kernel.apply_transpose(arr)` |

**Count for `FissionOperator`: 4 forward arms (F1/F2 share one register, so 3 registered
arms → 4 distinct arrows) + 3 transpose arrows = 7 distinct (in-type → out-type) arrows one
instance stands for.** Note F3/F4 differ only in wrapping; F1's AND FT1's angular arms are the
only ones that touch a quadrature (via `bulk.mesh.quad`, see §2).

`apply_transpose` is **NOT** singledispatch — it is a hand-rolled `isinstance(phi_star, FullField)`
chain + a duck-typed `getattr(..., "values", ...)` fallthrough (`fission.py:426`, `fission.py:464`).
Two different dispatch mechanisms in one class.

### 1c. `ScatteringOperator` — every (input → output) arrow

| # | Input type | Output type | Mechanism site | Body site | Math |
|---|---|---|---|---|---|
| S1 | `FullField` / `TimedFullField` | `FullField(interior=AngularSourceSink, boundary=AngularBoundarySourceSink.zeros_on(psi.mesh))` | `scattering.py:1134` register | `scattering.py:1148-1161` | re-dispatch `self.apply(psi.interior)` under `cast("AngularFlux \| HarmonicMomentFlux", …)` |
| S2 | `ScalarFlux` | `ScalarSourceSink` | `scattering.py:1163` register | `scattering.py:1183-1187` | `add_iso_source` then `add_n2n_source` on a `ScalarSourceSink.zeros_on(mesh)`. **No Pℓ, no 1/W.** |
| S3 | `AngularFlux` | `AngularSourceSink` | `scattering.py:1189` register | `scattering.py:1205-1207` | `_assemble_per_ordinate_source(psi.integrate_angular(), build_aniso_source(psi), psi.mesh)` |
| S4 | `HarmonicMomentFlux` | `AngularSourceSink` | `scattering.py:1209` register | `scattering.py:1231-1248` | `Λ(skip_l0=True).apply(moments)` → `frame.reconstruct(…)/W`, iso from `phi_moments.scalar_flux()` |
| S5 | anything else | `TypeError` | — | `scattering.py:1127` | — |
| ST1 | `FullField` (transpose) | `FullField(AngularSourceSink ⊕ zero AngularBoundarySourceSink)` | `isinstance`, `scattering.py:1311` | `scattering.py:1312-1322` | recursion on `chi.interior.values`, then re-wrap; `cast("SNMesh", chi.mesh)` at `1315` |
| ST2 | `np.ndarray` / `.values`-carrier (transpose) | `np.ndarray` | `scattering.py:1323` | `scattering.py:1324-1326` | `full_scatter_kernel.apply_transpose(chi)/W` |

**Count for `ScatteringOperator`: 4 forward arrows + 2 transpose arrows = 6 distinct
(in-type → out-type) arrows one instance stands for** (S1 is a delegating arrow, so counting
the delegated leaf it is 6; if you count S1's two possible interior families as separate
arrows it is 7).

Same asymmetry as fission: `apply` is singledispatch, `apply_transpose` is an `isinstance`
chain (`scattering.py:1311`) with a duck-typed fallthrough (`scattering.py:1323`).

**⚠ NOTE — S2 is a self-declared orphan.** `scattering.py:1172-1181` states in the body's own
docstring: *"Deliberately retained — a named-future-consumer surface, NOT dead weight. This arm
has no current production caller."* Verified by grep in §9: no `orpheus/` site feeds a bare
`ScalarFlux` to `S.apply`. Kept for #205.

---

## 2. Every attribute / property — method-NEUTRAL data vs quadrature-dependent

Legend: **N** = method-neutral (cross sections / group transfer / topology only) ·
**Q** = requires the discrete angular quadrature · **S** = requires a `FullFieldSpace` ·
**M** = requires an `SNMesh`-shaped carrier at call time.

### 2a. `FissionOperator` — the complete member list (`fission.py`)

| Member | Kind | Site | Class | What it needs |
|---|---|---|---|---|
| `mat_xs: MaterialXSField` | dataclass field (required, positional) | `168` | **N** | cross sections only |
| `full_field_space: FullFieldSpace \| None = None` | dataclass field (`repr=False, compare=False`) | `180-182` | **S**, optional | `None` by default |
| `block_role = BlockRole.BULK` | class constant (unannotated ⇒ not a field) | `188` | **N** | — |
| `is_adjointable` → `True` | property | `190-195` | **N** | hardcoded `True` |
| `chi` → `mat_xs.emission_spectrum` | property (read-through) | `199-205` | **N** | — |
| `sig_p` → `mat_xs.fission_production` | property (read-through) | `207-213` | **N** | — |
| `domain` → `self.full_field_space` | property | `216-231` | **S**, `None`-able | — |
| `codomain` → `self.full_field_space` | property | `233-236` | **S**, `None`-able | — |
| `from_solver_data(*, mat_xs, full_field_space=None)` | classmethod | `238-253` | **N** | quadrature NOT a parameter |
| `kernel` → `TensorProductOperator` | property, built FRESH each access | `255-309` | **N** | `outer(chi, production_rate) & IdentityOperator()` |
| `production_rate` → `ReactionRateFunctional` | property, built FRESH each access | `311-359` | **N** | `mat_xs.fission_production_field` |
| `apply_transpose` | method (`isinstance` chain) | `365-465` | **N** except the angular arm | angular arm reads `bulk.mesh.quad.weights` (`436-437`) — **Q at call time, via the CARRIER** |
| `_apply_impl` + 3 registered arms | singledispatchmethod | `467-656` | **N** except F1 | F1 uses `bulk.integrate_angular()` + `AngularSourceSink.from_isotropic` — the quadrature is inside the CARRIER, not the operator |
| `apply` | runtime alias / `TYPE_CHECKING` overloads | `658-674` | — | — |

**`FissionOperator` holds ZERO quadrature state.** It is method-neutral by construction: the
only angular knowledge it uses is read *off the carrier* (`bulk.mesh.quad` at `fission.py:436`;
`bulk.integrate_angular()` at `fission.py:561`). This is exactly why `diffusion` and
`homogeneous` can and DO construct it (§9).

### 2b. `ScatteringOperator` — the complete member list (`scattering.py`)

| Member | Kind | Site | Class | What it needs |
|---|---|---|---|---|
| `mat_xs: MaterialXSField` | dataclass field (required, positional) | `385` | **N** | cross sections + cell→material map |
| **`quadrature: Quadrature`** | **dataclass field (REQUIRED, positional, NO default)** | **`386`** | **Q — MANDATORY** | there is NO `None` branch anywhere in the class |
| `scattering_order: int` | dataclass field (required, positional) | `387` | **N** | Legendre truncation `L` |
| `block_role = BlockRole.BULK` | class constant | `393` | **N** | — |
| `_Y_cached: np.ndarray \| None` | dataclass field, `init=False, repr=False` — **mutable lazy cache on a mutable dataclass** | `398-400` | **Q** | populated in `Y` |
| `full_field_space: FullFieldSpace \| None = None` | dataclass field (`repr=False, compare=False`) | `410-412` | **S**, optional | `None` by default |
| `domain` / `codomain` → `full_field_space` | properties | `415-434` | **S**, `None`-able | — |
| `is_adjointable` → `True` | property | `436-442` | **N** | hardcoded |
| `n_ordinates` → `quadrature.N` | property | `452-455` | **Q** | |
| `spatial_shape` → `mat_xs.spatial_shape` | property | `457-464` | **N** | replaced the retired `nx`/`ny` (C5.2 #225) |
| `ng` → `mat_xs.ng` | property | `466-469` | **N** | |
| `weights` → `quadrature.weights` | property | `471-474` | **Q** | |
| `Y` → `quadrature.spherical_harmonics(L)` or `None` at `L==0` | property, MUTATES `_Y_cached` | `476-486` | **Q** (returns `None` at `L==0`) | |
| **`frame` → `HarmonicFrame`** | **`cached_property`** | **`488-515`** | **Q** | `HarmonicFrame.from_galerkin(quadrature.angular_frame(L))` |
| `sig_s` → `{mid: mat_xs.sig_s_legendre(mid)}` | property, marked **TRANSIENT** (retirement tracked as **#306**) | `517-529` | **N** | |
| `sig2` → `{mid: mat_xs.n2n_matrix(mid)}` | property, **TRANSIENT** | `531-537` | **N** | |
| `sig_s0` → `{mid: sig_s_legendre(mid)[0]}` | property, **TRANSIENT** | `539-546` | **N** | |
| `cells_by_mat` → `mat_xs.cells_by_material` | property, **TRANSIENT** | `548-552` | **N** | |
| `_aniso_source_from_moment_values` | private method | `556-596` | **Q** | `self.frame.reconstruct_after(Λ)` |
| `kernel` → `frame.conjugate(Λ_{ℓ≥1})` | property, built fresh; **RAISES `ValueError` at `L==0`** | `598-653` | **Q** | see §6 |
| `full_scatter_kernel` → `OperatorProduct` | property, built fresh | `655-679` | **Q** | `frame.conjugate(Λ_{ℓ≥0} + N₂ₙ)` |
| **`isotropic_kernel` → `OperatorSum`** | **`cached_property`** | **`681-709`** | **N** | `IsotropicScattering(mat_xs) + IsotropicN2N(mat_xs)`, space-anonymous |
| `_assemble_per_ordinate_source` | private method | `711-757` | **Q** + **M** | `quadrature.weights.sum()` at `756`; `mesh` sizes the zero accumulators |
| `from_solver_data(*, mat_xs, quadrature, scattering_order, full_field_space=None)` | classmethod | `759-786` | **Q — quadrature is a required kwarg** | |
| `add_iso_source(Q, phi)` | method + 2 overloads | `790-854` | **N** | `mat_xs.apply_p0_in_scatter`; in-place for ndarray, return-new for typed |
| `add_n2n_source(Q, phi)` | method + 2 overloads | `856-908` | **N** | `mat_xs.apply_n2n` |
| `build_aniso_source(angular_flux)` | method | `910-955` | **Q** | `self.kernel.apply(...)/weights.sum()` |
| `foldable_part()` | method | `978-1014` | **N** (but propagates `self.quadrature`) | see §4 |
| `residual_part()` | method | `1016-1055` | **N** (propagates `self.quadrature`) | see §4 |
| `is_foldable_into_sigma_r()` | method | `1057-1074` | **N** | `mat_xs.is_p0_diagonal_with_zero_n2n()` |
| `foldable_sigma()` | method | `1076-1089` | **N** | delegates to `mat_xs.foldable_sigma()` |
| `_apply_impl` + 4 arms / `apply` | singledispatch / alias | `1093-1270` | mixed | S2 is **N**; S1/S3/S4 are **Q** |
| `apply_transpose` | method (`isinstance` chain) | `1272-1326` | **Q** | `full_scatter_kernel` + `weights.sum()` |

### 2c. Verdict on the plan's claim (question 2 — the load-bearing one)

**The plan's statement is CORRECT and the code is stronger than the plan's phrasing.**

Evidence FOR "method-specific in fact":
1. `quadrature: "Quadrature"` at `scattering.py:386` is a **required positional dataclass field
   with no default**. There is no `quadrature: Quadrature | None` and no `if self.quadrature is
   None` branch anywhere in `scattering.py`. Constructing a `ScatteringOperator` **requires** a
   discrete angular quadrature, unconditionally — even for `scattering_order=0`.
2. `from_solver_data` (`scattering.py:759-767`) makes `quadrature` a **required keyword**.
   Contrast `FissionOperator.from_solver_data` (`fission.py:238-242`) which has no such parameter.
3. Even the pure-P0 per-ordinate path needs it: `_assemble_per_ordinate_source` divides by
   `float(self.quadrature.weights.sum())` at `scattering.py:756`.
4. `foldable_part()` / `residual_part()` — the *cross-section-only* split — still thread
   `quadrature=self.quadrature` into the derived sibling (`scattering.py:1011`, `1052`) even
   though the sibling has `scattering_order=0` and therefore never builds a frame. The quadrature
   rides along as **dead weight on the σ-split path**.
5. `frame` (`scattering.py:488-515`) constructs a `HarmonicFrame` from `quadrature.angular_frame(L)`
   — i.e. an SN-specific angular discretisation object lives *inside* the facade (§8).
6. **No non-SN production consumer exists.** `orpheus/cp`, `orpheus/moc`, `orpheus/mc`,
   `orpheus/kinetics`, `orpheus/fuel`, `orpheus/thermal_hydraulics` contain **ZERO** references
   to `orpheus.transport.operators` (grep over all six packages returned nothing).
   `orpheus/diffusion` and `orpheus/homogeneous` deliberately **bypass** it (§9).

Nuance the plan should NOT lose:
- The **method-neutral half is ALREADY extracted.** `isotropic_kernel` (`scattering.py:681-709`)
  is `IsotropicScattering(mat_xs) + IsotropicN2N(mat_xs)` from `isotropic_scattering.py` — bare
  `np.ndarray` in/out, no quadrature, no mesh, `space` optional (`isotropic_scattering.py:252`,
  `365`). That module's own docstring records the cross-model motivation (the 6× duplication,
  `isotropic_scattering.py:19-21`) and it is what diffusion/homogeneous consume directly.
- The **σ-split family** (`foldable_part`, `residual_part`, `is_foldable_into_sigma_r`,
  `foldable_sigma`, `sig_s`, `sig2`, `sig_s0`, `cells_by_mat`, `ng`, `spatial_shape`) is
  method-neutral *data*, currently welded onto a quadrature-bearing class.
- `LegendreMomentScattering` (`scattering.py:114-296`) and `N2NMomentOperator`
  (`scattering.py:299-352`) are **BOTH quadrature-free** — they hold `(mat_xs, L[, skip_l0])` and
  act on moment space. Their `domain`/`codomain` are `SphericalHarmonicSpace.from_L(self.L)`
  (`scattering.py:289-291`, `345-347`) — a *coefficient* space, not a quadrature. The quadrature
  enters only where the frame conjugates them.

So the honest structural statement is: **`ScatteringOperator` is a quadrature-bearing SN facade
wrapped around three already-method-neutral leaves (`IsotropicScattering`, `IsotropicN2N`,
`LegendreMomentScattering`) plus a σ-split data API that needs no quadrature at all.**

---

## 3. `full_field_space` — optional on BOTH; consumed only by the composition guard

| | `FissionOperator` | `ScatteringOperator` |
|---|---|---|
| declaration | `fission.py:180-182` | `scattering.py:410-412` |
| default | `None` | `None` |
| dataclass flags | `repr=False, compare=False` | `repr=False, compare=False` |
| read by | `domain` (`fission.py:231`), `codomain` (`fission.py:236`) — nothing else | `domain` (`scattering.py:429`), `codomain` (`scattering.py:434`); propagated into `foldable_part`/`residual_part` siblings (`scattering.py:1013`, `1054`) |
| threaded by | `from_solver_data` (`fission.py:253`) | `from_solver_data` (`scattering.py:785`) |

**When `None`:** `domain`/`codomain` report `None`, and `OperatorSum`/`OperatorProduct`'s
space-compatibility guard **silently skips** the operand — documented at `scattering.py:404-406`
and `fission.py:177-178`. Nothing raises; nothing else in either class reads the field.

Second, less obvious consequence (**not** stated in either docstring): `_AdjointOperator.apply`
(`orpheus/numerics/operator.py:1204-1227`) applies the metric ONLY when the inner operator's
`domain`/`codomain` are non-`None`:

```
z = inner_codomain.apply_metric(y) if inner_codomain is not None else y   # operator.py:1223
...
if inner_domain is not None: result = inner_domain.apply_inverse_metric(result)  # 1225-1226
```

So **`F.H` / `S.H` on a space-anonymous instance silently degrade to the bare Euclidean
transpose** (no `G`-weighting). See §7.

Threading status in production:
- SN: `sn/solver.py:1029` (S) and `sn/solver.py:1033` (F) both pass `sn_mesh.full_field_space`.
- SN coupled: `sn/coupled_system.py:458` (S fresh-construct fallback).
- SN adjoint pencil: `sn/solver.py:2354-2356` (F).
- diffusion: `diffusion/solver.py:241-243` (F, `full_field_space=space`).
- homogeneous: `homogeneous/solver.py:193` — `FissionOperator.from_solver_data(mat_xs=mat_xs)`
  **with NO space** (meshless single cell) ⇒ space-anonymous by design.

---

## 4. `foldable_part` / `residual_part` — STILL under those names, ZERO production callers

**They exist, unrenamed, unmoved:**

| Method | Site | Returns | Math |
|---|---|---|---|
| `foldable_part()` | `scattering.py:978-1014` | a NEW `ScatteringOperator` (`scattering.py:1009`) with `scattering_order=0`, diagonal-only P0, zeroed (n,2n), and `mat_xs.with_overridden_sig_s_and_n2n(...)` | the P0 within-group self-scatter `Σ_s0^{g→g}` — the piece foldable into `σ_r = σ_t − Σ_s0^{g→g}` |
| `residual_part()` | `scattering.py:1016-1055` | a NEW `ScatteringOperator` (`scattering.py:1050`) with `scattering_order = self.scattering_order`, zero-diagonal P0, Pℓ≥1 and (n,2n) verbatim | everything `foldable_part` is not |
| `is_foldable_into_sigma_r()` | `scattering.py:1057-1074` | `bool` | `L==0 AND mat_xs.is_p0_diagonal_with_zero_n2n()` |
| `foldable_sigma()` | `scattering.py:1076-1089` | `dict[int, np.ndarray]` `{mid: (ng,)}` | delegates to `mat_xs.foldable_sigma()` |

The stated algebraic contract (`scattering.py:1033-1035`):
`S.apply(ψ) ≈ S.foldable_part().apply(ψ) + S.residual_part().apply(ψ)` at `rtol=1e-14`.

### 4a. Did the DSA campaign (`e4c1a81c..37ffc310`, #215 CLOSED) rename or re-pose these? — NO.

Verified by grep across `orpheus/`, `tests/`, `docs/`:

- **The DSA implementation went AROUND them, straight to the XS field.**
  `orpheus/sn/acceleration/dsa.py:233` reads `mat_xs.foldable_sigma()` and
  `orpheus/sn/acceleration/dsa.py:236` reads `mat_xs.residual_sig_s()` — the
  `MaterialXSField` helpers (`material_xs_field.py:1148`, `1167`, `1213`), **not** the
  `ScatteringOperator` methods. So the operator-level split is bypassed, not consumed.
- **Zero production call sites of `S.foldable_part()` / `S.residual_part()`.** The only
  `orpheus/` occurrence of the token `foldable_part` outside `scattering.py` itself is a
  **docstring** in `orpheus/sn/operators/streaming.py:180` (`"(``L + C - S.foldable_part()``); no
  ``WithinGroupOperator`` wrapper"`) — verified prose inside the `StreamingOperator` class
  docstring, NOT a call. The actual within-group build
  (`sn/coupled_system.py:464`) is `A_AA = LC - S - B_a` — the **full** `S`, no split.
- **Nexus agrees**: the `foldable_part` node (`scattering.py:978`, degree 15) has its only
  `calls` edges from `tests/sn/operators/test_scattering_operator.py`.
- **The DSA campaign instead planted a SENTINEL against wiring them.**
  `tests/sn/acceleration/test_dsa_rate.py:1063` defines
  `_FOLD_ACCESSORS = frozenset({"foldable_sigma", "residual_sig_s"})` and reds if a sweep
  builder reaches for `foldable_sigma` (`test_dsa_rate.py:1100`, mutation tooth at
  `docs/theory/methods/sn/acceleration.rst:1024`). This is the codified counterpart of the
  in-file warning at `scattering.py:968-977`:
  *"⚠ LATENT CORRECTNESS TRAP (#215): do NOT wire the σ_r-SWEEP as the within-group
  `A_wg.inverse()` … ships 46–56 % silent flux errors on anisotropic problems."*

**⚠ FLAG for the plan.** If the plan asserts (or assumes) that #215's closure renamed or
re-posed `foldable_part`/`residual_part`, that is **wrong**: #215 was closed by *building
consistent DSA elsewhere* (`orpheus/sn/acceleration/dsa.py`) and by *fencing off* the
σ_r-fold wiring with a sentinel. The two methods are untouched, unrenamed, and remain a
**test-only data API** — a genuine `coding-standards` "aggressive retirement" candidate (or a
documented keep-as-anchor; that is the user's call per lessons L-004, not mine).

### 4b. Consumer inventory for the σ-split family

| Consumer class | Count | Representative sites |
|---|---|---|
| production (`orpheus/`) — **operator methods** | **0** | — (only the docstring mention at `sn/operators/streaming.py:180`) |
| production (`orpheus/`) — the `mat_xs.*` helpers instead | 2 | `sn/acceleration/dsa.py:233`, `:236` |
| tests — behavioural | ~30 assertions | `tests/sn/operators/test_scattering_operator.py:749-1231` (`TestFoldablePart`, `TestResidualPart`, `TestFoldableSigma`, `TestIsFoldableIntoSigmaR`, additivity at `:942-953`) |
| tests — DSA fold sentinel | 1 module | `tests/sn/acceleration/test_dsa_rate.py:1063-1134` |
| tests — other fold consumers | 3 | `tests/derivations/test_dsa_production_tie.py:132`, `tests/sn/acceleration/test_dsa_low_order.py:66,68` (all via `xs.foldable_sigma()`) |
| docs | 4 | `docs/theory/conventions/indexing_and_layout.rst:1055`, `docs/theory/methods/sn/slab_one_group.rst:792-793`, `docs/theory/methods/sn/acceleration.rst:988-989,1024` |

---

## 5. `(n,2n)` — secondary-neutron production lives in THREE places, at two layers

It is **both** a separate operator **and** a channel inside `ScatteringOperator`. Precisely:

| Realisation | Type | Site | Carrier | Consumed by |
|---|---|---|---|---|
| `IsotropicN2N` | `@dataclass(frozen=True)` `LinearOperator`, `2Σ₂ₙ` on the scalar flux | `isotropic_scattering.py:344-438` | bare `ndarray` (or scalar `FullField` via `_scalar_composite_source`, `isotropic_scattering.py:382-383`) | `S.isotropic_kernel` (`scattering.py:709`); `diffusion/solver.py:237`; `homogeneous/solver.py:146` |
| `N2NMomentOperator` | `@dataclass(frozen=True)` `LinearOperator`, `2Σ₂ₙ` on the **ℓ=0 moment block** | `scattering.py:299-352` | bare `ndarray` moment tensor | ONLY `full_scatter_kernel` (`scattering.py:678`) — i.e. the adjoint path |
| `ScatteringOperator.add_n2n_source(Q, phi)` | in-place / return-new accumulator | `scattering.py:856-908` | `ndarray` (mutate) or `ScalarSourceSink` (return-new) | the `ScalarFlux` arm S2 (`scattering.py:1186`); `sn/solver.py:1891` delegator `_add_n2n_source`, itself called from `sn/solver.py:2110` |

All three route to the SAME per-material kernels on `MaterialXSField`:
`apply_n2n` / `apply_n2n_transpose` (`isotropic_scattering.py:386`, `404`;
`scattering.py:902`) and `apply_n2n_moments` / `apply_n2n_moments_transpose`
(`scattering.py:336`, `340`). Single source of XS truth, three carriers.

**Physics posture recorded in the code (consistent across all three sites):** (n,2n) is a
*multiplication* channel, NOT scattering, kept as its own named operator on purpose —
`scattering.py:301-311` and `isotropic_scattering.py:27-30`. It feeds the in-scatter source
when summed with `Σ_s0`, but `isotropic_scattering.py:29` also notes it "feeds the keff
*production* numerator, not the loss". `homogeneous/solver.py` (module docstring, cited at
`:176`) and `diffusion/solver.py:236-238` both fold it **loss-side** (`A = C − K_iso`), and
`FissionOperator` carries **no** (n,2n) at all — `F` is purely `χ ⊗ νΣ_f`.

**⚠ ASYMMETRY to note:** the *forward* production hot path uses the SCALAR
`IsotropicN2N` (through `isotropic_kernel`, `scattering.py:749`), while the *adjoint*
uses the MOMENT `N2NMomentOperator` (through `full_scatter_kernel`, `scattering.py:678`,
consumed at `scattering.py:1324`). Two realisations of `2Σ₂ₙ` on two carriers; the
crosscheck oracle relationship is asserted at `isotropic_scattering.py:40-41`
("the harmonic-moment `full_scatter_kernel` is the permanent verification oracle for this
scalar form").

---

## 6. `.kernel` on `FissionOperator` — a concrete `TensorProductOperator`, already a `LinearOperator`

`fission.py:255-309`:

```python
@property
def kernel(self) -> "TensorProductOperator":
    return outer(self.chi, self.production_rate) & IdentityOperator()
```

- **Return type**: `TensorProductOperator` (annotated concretely, `fission.py:256`) — a 2-factor
  TP built by the `&` operator. Factor 1 = `outer(chi, production_rate)` (a `RankOneOperator`,
  per `orpheus/numerics/operator.py:3059` `outer`); factor 2 = `IdentityOperator()` (the
  spatial-axis broadcast).
- **Is it already a `LinearOperator`?** **Yes.** `TensorProductOperator` is a numerics composer
  in `orpheus/numerics/operator.py`; `FissionOperator` calls `kernel.apply` (`fission.py:639`,
  `656`) and `kernel.apply_transpose` (`fission.py:440`, `451`, `465`) directly. No adapter.
- **Domain / codomain**: `TensorProductOperator` carries whatever its factors declare;
  `outer(...)` and `IdentityOperator()` here are space-anonymous — the kernel's operational
  contract is the bare `(ng, *spatial)` array (input scalar flux → output `(ng, *spatial)` iso
  fission source). Confirmed by the two arms that feed it raw `.values` /raw ndarray.
- **Built fresh on every access** (`fission.py:280-284` documents this, body confirms) to honour
  the `mat_xs` read-through for depletion / thermal feedback.
- `production_rate` (`fission.py:311-359`) is the row co-vector — a `ReactionRateFunctional` over
  `mat_xs.fission_production_field`. **`apply` genuinely ROUTES THROUGH it** (verified: `kernel`
  at `fission.py:309` embeds `self.production_rate`, and F3/F4 both call `self.kernel.apply`).
  The docstring's "the procedural twin is dissolved" claim (`fission.py:333-334`) **checks out** —
  there is no second einsum in the file.

### 6b. Contrast — `ScatteringOperator.kernel` is a different animal

`scattering.py:598-653`. Return annotation is the **erased** `LinearOperator` (not
`OperatorProduct`), and it **RAISES `ValueError` when `scattering_order == 0`**
(`scattering.py:638-644`). Body: `self.frame.conjugate(LegendreMomentScattering(..., skip_l0=True))`
⇒ `OperatorProduct(R, OperatorProduct(Λ, M))` (per `numerics/frame.py:236-238`). It is the
ℓ≥1 **anisotropic subcomponent only**; the FULL in-scatter is the separate
`full_scatter_kernel` (`scattering.py:655-679`, annotated concretely as `OperatorProduct`).

**⚠ DOC DRIFT (silent Python-domain xref).** `docs/api/numerics.rst` references
`ScatteringOperator.kernel_summands` at **three** sites (`:225`, `:252`, `:273`). **No such
member exists** anywhere in the repo (grep for `kernel_summands` returns only those 3 doc
lines). Per `.claude/rules/coding-standards.md`, a Python-domain `:attr:` role renders as plain
text with **no `-W` warning** unless the build runs `-n`, so the Sphinx gate does not catch this.

---

## 7. Transpose / adjoint — both are BARE EUCLIDEAN; `.H` is the metric wrapper, and it can silently no-op

### 7a. `FissionOperator.apply_transpose` (`fission.py:361-465`)

- Mechanism: **`isinstance` chain**, not singledispatch. `@overload`s at `fission.py:362-364`.
- Math: `F† = |νΣf⟩⟨χ|` — the χ↔νΣf dyad swap, obtained by routing through
  `self.kernel.apply_transpose` (`fission.py:440`, `451`, `465`). No new kernel code; the swap
  lives on `RankOneOperator.apply_transpose`.
- Composite angular arm (`fission.py:429-445`) is the ONLY place `FissionOperator` reads a
  quadrature: `w = np.asarray(mesh.quad.weights)` (`437`), `iso_star = ψ*.sum(axis=0)/w.sum()`
  (`437`), then `np.multiply.outer(w, kernel.apply_transpose(iso_star))` (`439-441`). This is
  the honest pullback of `from_isotropic ∘ kernel ∘ integrate_angular`.
- **Explicitly Euclidean**: `fission.py:385-387` — *"This is the Euclidean transpose `F^T`; the
  metric-correct Hilbert adjoint `F† = G⁻¹F^T G` is the `.H` wrapper's job."* Body confirms —
  no metric application anywhere in the method.

### 7b. `ScatteringOperator.apply_transpose` (`scattering.py:1272-1326`)

- Mechanism: **`isinstance` chain** + duck-typed fallthrough. `@overload`s at `1273-1275`.
- Math: `S^T χ = (1/W)·full_scatter_kernel.apply_transpose(χ)` (`scattering.py:1324-1326`) —
  ONE expression covering iso ℓ=0 + aniso ℓ≥1 + (n,2n), because all three are conjugated by the
  SAME frame in `full_scatter_kernel`.
- Composite arm (`1311-1322`) recurses on `chi.interior.values`, then re-wraps as
  `AngularSourceSink ⊕ zero AngularBoundarySourceSink`, with `cast("SNMesh", chi.mesh)` at `1315`
  — an **unchecked SN-mesh cast**, the transpose path's one SN-specific weld.
- **Explicitly Euclidean**: `scattering.py:1290-1291` — *"This is the **Euclidean** transpose
  (L12) — NOT the metric Hilbert adjoint `.H` (which would carry the angular Gram)."*
- **Forward/adjoint structural asymmetry is DELIBERATE**: forward uses the scalar fast path
  (`isotropic_kernel` + `build_aniso_source`, `scattering.py:1198-1204`), the adjoint uses the
  moment frame form. `scattering.py:1200-1204` records the reason (routing the iso forward through
  the frame "regresses LD/P0 badly"), and `1292-1295` claims this asymmetry is what makes the
  reciprocity check a genuine cross-check rather than a tautology.

### 7c. `.H` — where the metric actually lives, and the `None`-space trap

`.H` is NOT defined on either facade. It is inherited: `LinearOperator.H` →
`LinearOperator.adjoint()` → `_AdjointOperator(self)` (`orpheus/numerics/operator.py:857-862`,
`847-854`). `adjoint()` **eagerly** raises `MissingAdjoint` unless `is_adjointable` — both facades
hardcode `True` (`fission.py:191-195`, `scattering.py:436-442`).

`_AdjointOperator.apply` (`operator.py:1204-1227`) is the metric-correct
`A† = G_V⁺ ⊙ apply_transpose(G_W ⊙ y)` — **but conditional on the spaces**:

```
inner_codomain = getattr(self.inner, "codomain", None)      # operator.py:1221
z = inner_codomain.apply_metric(y) if inner_codomain is not None else y   # 1223
result = self.inner.apply_transpose(z)                       # 1224
if inner_domain is not None: result = inner_domain.apply_inverse_metric(result)  # 1225-1226
```

⇒ **`F.H` / `S.H` built from a `full_field_space=None` instance apply NO metric** and reduce
exactly to `apply_transpose`. This is not documented on either facade. It bites the
space-anonymous constructions — notably `homogeneous/solver.py:193`
(`FissionOperator.from_solver_data(mat_xs=mat_xs)`, no space) and every bare/test constructor.

Also: `_AdjointOperator.apply_transpose` raises `NotImplementedError` (`operator.py:1229-1234`) —
`(S.H).H` is not available.

Production `.H` consumers on these facades (grep): `orpheus/numerics/eigenvalue.py:29`
(docstring: `KEigenvalue((L+C).H, (S+B).H, F.H)`), `orpheus/numerics/iteration.py:1312`
(docstring: `(A.H, S.H, F.H)`), `orpheus/sn/solver.py:2429` (docstring for the daggered pencil);
the live wiring is `tests/numerics/test_iteration.py:864` — `LC.H, S_total.H, F.H`.

---

## 8. Where the harmonic frame is built — INSIDE the scattering facade

| Layer | Class | Site | Role |
|---|---|---|---|
| numerics | `FrameBase` (ABC, `@dataclass(frozen=True)`) | `numerics/frame.py:113-330` | discipline-free mechanics: `table`, `basis_space`, `test_table`, `test_space`, `measure_space`, `analysis`, `reconstruction`, `conjugate`, `reconstruct_after`, `gram`, `project` |
| numerics | `PetrovGalerkinFrame` | `numerics/frame.py:333-370` | explicit `test_basis` |
| numerics | `GalerkinFrame` (`frozen=True, init=False`) | `numerics/frame.py:373-408` | `test IS trial`; `__init__(basis, measure)` sets 3 fields via `object.__setattr__` |
| numerics | `_FrameAnalysis` / `_FrameReconstruction` | `numerics/frame.py:411-479` | the two `LinearOperator` faces (`M`, `R`), both with a working `apply_transpose` |
| transport | **`HarmonicFrame(GalerkinFrame)`** | `frames/harmonic_frame.py:83-208` | adds only the **carrier-typed** verbs `analyse` (M) / `reconstruct` (R); narrows `basis` to `SphericalHarmonicBasis` (`:102`, enforced at `:124-128`) |

**Which class builds a quadrature-bearing frame?** The **quadrature** does:
`Quadrature.angular_frame(L)` (`orpheus/numerics/quadrature/directional.py:417`) returns the
generic `GalerkinFrame`. `HarmonicFrame.from_galerkin(frame)` (`harmonic_frame.py:109-129`)
*upgrades* it, **reusing** its basis + measure — no rebuild, table bit-identical.

**Built INSIDE or injected?** **Built INSIDE the scattering facade**, as a `cached_property`:

```python
@cached_property
def frame(self) -> HarmonicFrame:                                  # scattering.py:488-489
    return HarmonicFrame.from_galerkin(
        self.quadrature.angular_frame(self.scattering_order),      # scattering.py:513-515
    )
```

There is **no** `frame=` constructor parameter and no injection seam on `ScatteringOperator`.
The frame is derived from `(self.quadrature, self.scattering_order)`.

**And the frame LEAKS OUT of the facade to a non-scattering consumer** — a load-bearing coupling
the plan should know about: `orpheus/sn/solver.py:552` does

```
BulkAnalysisOperator(scattering_op.frame, sn_mesh) @ sweep
```

i.e. the 2-D Cartesian **angular-windowing** path reaches into `S` for its frame. Also
`orpheus/sn/solver.py:573` reads `scattering_op.scattering_order`. So `ScatteringOperator` is
currently the **de facto owner and distributor of the SN angular frame**, not just a scatter
source. Design rationale for that ownership is recorded in
`docs/theory/foundations/frame.rst:3145` and `:3283-3284` ("an operator owns its frame IFF the
frame is its eigenbasis" — Funk–Hecke); `frame.rst:3328` discusses the one-`ScatteringOperator`-
per-method consequence.

Internal frame consumers within `scattering.py`: `_aniso_source_from_moment_values`
(`:596` — `frame.reconstruct_after(scatter)`), `kernel` (`:649` — `frame.conjugate`),
`full_scatter_kernel` (`:674` — `frame.conjugate`), the `HarmonicMomentFlux` arm
(`:1243` — `frame.reconstruct`).

---

## 9. Consumer inventory

Method: Nexus `context` on both classes **plus** repo-wide grep. Nexus alone is
**insufficient** here for two confirmed reasons (lessons L-001 sharpening): (a) `apply` is a
runtime alias (`apply = _apply_impl`) so per-arm arrows are invisible; (b) the four registered
arms are all named `_`, collapsing into one graph node. Grep is primary; Nexus corroborates.

### 9a. `FissionOperator`

| Bucket | Count | Sites |
|---|---|---|
| **`orpheus/sn/`** | 4 | `sn/solver.py:48` (import), `:1031-1034` (construct, `self.fission_op`), `:1042` (`self.F` alias), `:1132` (`self.fission_op.apply(flux_distribution) / keff` — the `compute_fission_source` delegator), `:2345`+`:2354` (adjoint pencil: local import + construct) |
| **`orpheus/diffusion/`** | 2 | `diffusion/solver.py:155` (import), `:241-243` (construct, `self.fission`); consumed at `:285` `self.fission.apply(psi).to_flat()/keff` |
| **`orpheus/homogeneous/`** | 2 | `homogeneous/solver.py:46` (import), `:193` (`production = FissionOperator.from_solver_data(mat_xs=mat_xs)`); consumed at `:195` `MatrixInverseOperator(loss, basis_shape=(ng,1)) @ production` |
| **`orpheus/cp/`, `moc/`, `mc/`, `kinetics/`, `fuel/`, `thermal_hydraulics/`** | **0** | grep over all six returned nothing |
| elsewhere in `orpheus/` (docstring refs only) | 5 | `sn/operators/radial_characteristic.py:1130`+`:1186` (accepts `fission_op.kernel`), `sn/operators/streaming.py:325`, `numerics/iteration.py:71`, `numerics/operator.py:3081`, `diffusion/operators.py:19`, `diffusion/__init__.py:45` |
| **tests** | 9 modules | `tests/sn/operators/test_fission_operator.py`, `test_fission_kernel_crosscheck.py`, `test_fission_adjoint.py`, `test_capability_survival.py:248` (direct ctor `FissionOperator(mat_xs=mat)`), `test_operators_apply_typed.py:140,156,291,326,428`, `tests/diffusion/test_operators.py:650`, `tests/homogeneous/test_homogeneous.py:285,359`, `tests/numerics/test_iteration.py:601`, `tests/data/test_emission_spectrum.py:102` |
| **docs** | 25 doc lines / 16 pages | `docs/api/{transport,numerics,homogeneous}.rst`, `docs/theory/foundations/{operator_algebra,operator_adjoint,operator_tensor_network,infinite_medium,path_integral,index,boundary_conditions,cross_section_data}.rst`, `docs/theory/methods/{sn/index,sn/slab_multigroup,sn/solver,sn/adjoint,diffusion_1d}.rst`, `docs/theory/conventions/indexing_and_layout.rst` |

**Structural conclusion: `FissionOperator` IS genuinely cross-method** — three production
solver families construct it (sn / diffusion / homogeneous), which is exactly what its
quadrature-free state permits.

### 9b. `ScatteringOperator`

| Bucket | Count | Sites |
|---|---|---|
| **`orpheus/sn/`** | ~15 live sites | `sn/solver.py:71` (import), `:1025-1030` (construct `self.scattering_op`), `:1041` (`self.S`), `:522` + `:552` (`scattering_op.frame` → `BulkAnalysisOperator`), `:558`+`:573` (`scattering_op.scattering_order`), `:708`+`:857-860` (`_select_si_resolvent` type-guard `isinstance(S, ScatteringOperator)`), `:1569`, `:1757`, `:2141`, `:3172`, `:3196`, `:3403` (threaded as `scattering_op=` into `build_within_group_system` / windowed cold-start), `:1854`/`:1880`/`:1891` (the three thin delegators), `:2110` (`solver._add_n2n_source`); `sn/coupled_system.py:145` (import), `:384`+`:451-460` (inject-or-construct), `:464` (`A_AA = LC - S - B_a`), `:496` (`RadialCharacteristicEmission(sn_mesh, S.isotropic_kernel)`), `:519` (`N_AA = S + B_a`) |
| **`orpheus/diffusion/`** | **0** | **bypasses it**: `diffusion/solver.py:236-237` builds `IsotropicScattering(mat_xs, space) + IsotropicN2N(mat_xs, space)` directly |
| **`orpheus/homogeneous/`** | **0** | **bypasses it**: `homogeneous/solver.py:146` `k_iso = IsotropicScattering(mat_xs) + IsotropicN2N(mat_xs)` |
| **`orpheus/cp/`, `moc/`, `mc/`, `kinetics/`, `fuel/`, `thermal_hydraulics/`** | **0** | grep over all six returned nothing |
| elsewhere in `orpheus/` (docstring refs only) | 6 | `sn/operators/radial_characteristic.py:153`+`:1127`+`:1183`, `sn/operators/streaming.py:324`, `numerics/iteration.py:70`, `numerics/operator.py:561`, `numerics/quadrature/directional.py:435`, `transport/mesh/material_xs_field.py:8-15,226,307-308,825,861` |
| **tests** | 12+ modules | `tests/sn/operators/test_scattering_operator.py` (the big one — foldable/residual/apply/is_foldable), `test_scattering_adjoint.py`, `test_scattering_kernel_crosscheck.py`, `test_legendre_moment_scattering.py`, `test_frame_conjugate_carve.py`, `test_capability_survival.py:251` (direct ctor), `test_green_operator_sn.py:124` (direct ctor), `test_operators_apply_typed.py:187,286,321,429`, `test_psi_half_coupling.py:2867`, `tests/sn/solve/test_sn_adjoint_certification.py`, `tests/sn/verification/mms/test_curvilinear_aniso_scattering_p1.py`, `tests/sn/verification/mms/test_mms_ld_2d.py`, `tests/sn/primitives/test_harmonic_moment_flux.py`, `tests/numerics/test_iteration.py` |
| **docs** | ~60 doc lines / 20 pages | heaviest: `docs/theory/foundations/{operator_algebra,boundary_conditions,frame,operator_tensor_network,operator_adjoint,spherical_harmonics}.rst`, `docs/theory/conventions/{indexing_and_layout,normalization}.rst`, `docs/theory/methods/sn/{slab_multigroup,slab_one_group,adjoint,cartesian_multid,solver,history}.rst`, `docs/api/{transport,numerics}.rst`, `docs/theory/verification/sn.rst:2296` |
| derivations | 1 | `derivations/diagnostics/diag_276_scattering_p0_fastpath_perf.py:92-93` (direct `LegendreMomentScattering` / `N2NMomentOperator` ctors) |

**Structural conclusion: `ScatteringOperator` has ZERO non-SN production consumers.** The
non-SN families reach for the method-neutral leaves instead. That is the empirical form of the
plan's "method-specific in fact".

### 9c. Direct-constructor audit (bypassing `from_solver_data`)

`ScatteringOperator(...)` direct: `scattering.py:1009` + `:1050` (its own `foldable_part` /
`residual_part`), `tests/sn/operators/test_capability_survival.py:251`,
`test_green_operator_sn.py:124`, `test_scattering_operator.py:1147,1172,1198,1226`.
`FissionOperator(...)` direct: `tests/sn/operators/test_capability_survival.py:248`.
`IsotropicScattering(...)` / `IsotropicN2N(...)` direct: `diffusion/solver.py:236-237`,
`homogeneous/solver.py:146`, `scattering.py:709`, plus 9 test sites.
`MultiplicationOperator(...)` direct: `diffusion/solver.py:230`, `sn/coupled_system.py:374`,
`homogeneous/solver.py:144` (`.from_mesh`), plus ~10 test sites.

---

## 10. Sibling facades — `isotropic_scattering.py` and `multiplication_operator.py` (for contrast)

### `isotropic_scattering.py` — the model-portable leaves

Two `@dataclass(frozen=True)` `LinearOperator`s: `IsotropicScattering` (`:229`) and
`IsotropicN2N` (`:345`). Fields on each: `mat_xs: MaterialXSField` + `space: FunctionSpace |
None = None` (`:251-252`, `:364-365`); `block_role = BlockRole.BULK` class constant
(`:255`, `:366`).

Surface per class: `is_adjointable → True`; `apply` (bare `ndarray` **or** scalar `FullField`
→ scalar source composite, via the shared `_scalar_composite_source`, `:101-136`);
`apply_transpose` (**bare ndarray ONLY** — it RAISES `TypeError` on a `FullField`, `:286-292`
and `:395-401`, "composite FullField transpose is not yet wired (lands with the adjoint
diffusion consumer, #281)"); `is_assemblable` / `assemble` (the DSA-campaign stencil-assembly
mode, shared bodies `_iso_is_assemblable` `:139-158` and `_assemble_iso_energy_operator`
`:161-225`); `dense_per_material` (the **storage-side oracle**, `:310-333`, `:419-430`, kept
deliberately realization-independent per vv L11); `domain`/`codomain` → `self.space`.

**Dispatch mechanism: plain `isinstance` at the top of `apply`** (`:272`, `:382`) — a THIRD
mechanism, different from both facades' singledispatch. Zero quadrature, zero frame, zero mesh
in state (mesh read off the parsed bulk at `:132`).

**⚠ Asymmetry worth flagging:** these two leaves have `assemble()` (a live DSA consumer), while
neither `ScatteringOperator` nor `FissionOperator` implements `is_assemblable`/`assemble` at all.
`MultiplicationOperator` does (`:244-296`).

### `multiplication_operator.py` — `M[f]`, the promotion

`@dataclass(eq=False)` `MultiplicationOperator(LinearOperator["FullField"])` (`:129-130`).
Fields: `coefficient: CrossSectionField` (`:163`), `space: FunctionSpace | None = None`
(`:175`), `engine: DiagonalOperator` (`init=False`, built once in `__post_init__`, `:183`,
`:205-207`); `block_role = BlockRole.BULK` (`:191`).

Surface: `is_invertible` **delegated to the engine's spectrum gate** (`:210-215` — `spec(M[f]) =
ess-range(f)`, invertible iff `min|f| > 0`); `inverse()` → `InverseOperator(self)`, raising
`NotInvertible` otherwise (`:217-234`); `is_adjointable` delegated (`:237-240`);
`is_assemblable`/`assemble` (`:244-296`); `from_mesh` classmethod (`:298-337`, defaults `space`
from `getattr(mesh, "full_field_space", None)` at `:336`); `domain`/`codomain` → `space`;
`_apply_impl` **singledispatch** with 2 arms (`FullField` `:384-443` with an angular/scalar
family parse; bare `ndarray` `:445-458`) + the `apply = _apply_impl` alias (`:470`);
`solve(q: FullField)` (`:472-516`, same two-family parse, returns to a FLUX role);
`apply_transpose` → **literally `return self.apply(psi)`** (`:522-530`) because a real
multiplier is self-adjoint.

This is the ONLY one of the four with a `solve` / `inverse` — `#261` retired the former
`CollisionOperator` thin subclass into it (`:50-55`, `:318-321`), so `C = M[σ_t]` is now a plain
`MultiplicationOperator`.

---

## 11. Every place the plan's description does NOT match the code

| # | Claim | Reality | Evidence |
|---|---|---|---|
| P1 | "`ScatteringOperator` is method-specific in fact because it needs a quadrature" | **CONFIRMED, stronger.** `quadrature` is a required positional dataclass field with no default and no `None` branch; it is also required by `from_solver_data`; and ZERO non-SN production consumers exist | `scattering.py:386`, `:766`; grep over cp/moc/mc/kinetics/fuel/TH = 0 hits |
| P2 | (implied) `foldable_part`/`residual_part` may have been renamed / re-posed by the DSA campaign that closed #215 | **NOT renamed, NOT re-posed, NOT consumed.** Both exist verbatim; DSA reads `mat_xs.foldable_sigma()` / `mat_xs.residual_sig_s()` instead; production callers of the operator methods = **0**; a sentinel test now FORBIDS the σ_r-fold wiring | `scattering.py:978`, `:1016`; `sn/acceleration/dsa.py:233,236`; `tests/sn/acceleration/test_dsa_rate.py:1063-1134`; the in-file trap note `scattering.py:968-977` |

### Docstring-vs-body drifts found while verifying (L33 discipline)

| # | Drift | Site |
|---|---|---|
| D1 | **`ScatteringOperator.kernel_summands` does not exist** but is `:attr:`-referenced 3× in the API docs. Silent (no `-W` warning for Python-domain roles) | `docs/api/numerics.rst:225`, `:252`, `:273` |
| D2 | `fission.py`'s module docstring writes the algebra as `(L − S − F)ψ = q` and `(L − S)ψ = (1/k)Fψ` — the OLD C-folded-into-L spelling. The ratified project-wide spelling is `A = L + C − S − B` (used correctly in `scattering.py`'s own module docstring) | `fission.py:8`, `:13`, `:16-17` vs `scattering.py:4-6` |
| D3 | `FissionOperator`'s `Attributes` block glues the capability-set prose INTO the `mat_xs` parameter description (mis-indented — reads as if it documents `mat_xs`) | `fission.py:158-166` |
| D4 | The composite `FullField` arm's docstring says the output boundary is an "implicit zero `AngularBoundaryFlux`" — the body emits `AngularBoundarySourceSink` (a source role, not a flux role) | claim `fission.py:533-535`; body `fission.py:593` |
| D5 | A runtime comment claims the "legacy `orpheus.sn.angular_flux.AngularFlux` … still rides on the operator-algebra path until D-H.1c". **`orpheus/sn/angular_flux.py` does not exist.** Six PRODUCTION `orpheus.sn.angular_flux` xrefs are dangling (silent) | comment `fission.py:125-127`; dangling xrefs `numerics/iteration.py:579`, `sn/solution.py:255`, `transport/timed_full_field.py:40`, `sn/solver.py:1469,1692,2070` |
| D6 | `MaterialXSField` cites `scattering.py:405` / `scattering.py:426` as the `add_iso_source` / `add_n2n_source` sites — the real sites are `800` / `866` | `transport/mesh/material_xs_field.py:825`, `:861` |
| D7 | `scattering.py:448-451` says the TRANSIENT read-throughs are pending "a follow-up PR"; `sig_s`'s own docstring now names **#306** and states its original `_build_rhs_*` consumers "are gone". Four read-throughs (`sig_s`, `sig2`, `sig_s0`, `cells_by_mat`) survive with **no production caller** | `scattering.py:444-552`; retirement tracked as #306 |

### Structural asymmetries the plan should carry forward (facts, not proposals)

1. **Three different dispatch mechanisms** across the four files for the same job:
   `singledispatchmethod` + alias (fission, scattering, multiplication `apply`);
   hand-rolled `isinstance` chain (all four `apply_transpose`, fission's composite bulk parse);
   plain top-of-method `isinstance` (isotropic_scattering `apply`).
2. **`apply` is singledispatch, `apply_transpose` is not** — on BOTH facades.
3. **Fission hand-parses its composite bulk; scattering re-enters its own dispatcher**
   (`fission.py:560,595` vs `scattering.py:1153`).
4. **Forward and adjoint use different realisations of the same math** in scattering (scalar
   fast path vs moment frame form) — deliberate, perf-motivated, documented.
5. **The frame leaks out of the facade**: `sn/solver.py:552` consumes `scattering_op.frame` for
   angular windowing — `ScatteringOperator` is the de facto owner/distributor of the SN angular
   frame.
6. **`assemble()` exists on the leaves but not the facades**: `IsotropicScattering`,
   `IsotropicN2N`, `MultiplicationOperator` are assemblable; `ScatteringOperator` and
   `FissionOperator` are not.
7. **`ScatteringOperator.kernel` raises at `L=0`** while `FissionOperator.kernel` always
   returns; the two `.kernel` properties are not the same contract (ℓ≥1 subcomponent vs the
   whole operator).
8. **`_Y_cached` is mutable state on a mutable dataclass** (`scattering.py:398-400`), whereas
   `LegendreMomentScattering`, `N2NMomentOperator`, `IsotropicScattering`, `IsotropicN2N` are
   all `frozen=True`. `ScatteringOperator` and `FissionOperator` are the only non-frozen ones.
