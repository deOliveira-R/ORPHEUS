---
name: issue-196-phase-g-replan-sn-strategies
description: Strategy-Protocol map for the SN spatial module — CellUpdate, PsiHalfAngleSeed, PoleAngularClosure, BoundaryRealizer. The Phase G replan extends these strategies; it does NOT wrap them as LinearOperator subclasses (the failure mode of the previous Step 1).
metadata:
  type: project
---


**Purpose.** Map every strategy Protocol the SN module already ships,
so the Phase G replan extends them rather than reinventing them as
`LinearOperator` subclasses. The previous Step 1 attempt
(uncommitted `orpheus/sn/spatial/operators.py`) promoted two
strategies to `LinearOperatorMixin` subclasses — wrong architecture
per user correction. This memo is the corrected baseline.

**Status of the working tree (2026-05-13).**
`git status` shows the SN spatial module mid-refactor:

* Modified: `orpheus/sn/spatial/__init__.py`, `diamond.py`,
  `operators.py`, `sweep.py`.
* Untracked: `cell_balance.py` (Step 2's correct twin-path fix) and
  `streaming.py`.

The Step 1 misplacement (`SNCellOperator` /
`AngularRedistribution` as `LinearOperatorMixin` subclasses) is still
on disk. The Step 2 helper (`cell_balance.py`) is the **salvageable**
piece — see §4 below.

---

## 1. `CellUpdate` Protocol

### Protocol declaration

* File: `orpheus/sn/spatial/cell_update.py:418-502`.
* Decoration: `@runtime_checkable` (line 418).
* Inherits: `Protocol` only (no operator/algebra inheritance).
* Class-level traits (mandatory attribute *contracts*, not methods):
  * `is_linear: bool` (line 454).
  * `is_positivity_preserving: bool` (line 455).
* Method signature (lines 457-502):
  ```python
  def update(
      self,
      visit: CellVisit,
      total_xs: np.ndarray,            # (ng,)
      source: np.ndarray,              # (ng,), already weight-normalised
      upstream_state: UpstreamState,
  ) -> CellResult: ...
  ```
* Docstring claims (load-bearing for downstream architecture decisions):
  * Source arrives **already weight-normalised** by the sweep
    (line 484: "the sweep applies its 1/Σ_n w_n factor at the call
    site").
  * Sweep-direction is **pre-resolved** in `visit.face_area_downstream`
    so the strategy "sees no sign-of-μ branching" (line 489).
  * Slab discrimination: `visit.streaming_terms.alpha_in is None`
    (line 474).
  * Cylindrical-degenerate flag:
    `visit.streaming_terms.abs_mu < 1e-15` (line 476).

### Concrete ABC: `CellUpdateBase`

* File: `orpheus/sn/spatial/cell_update.py:510-601`.
* Inherits: `RegistryMixin, ABC` (line 510).
* Registry shape:
  * `registry: ClassVar[dict[str, type["CellUpdateBase"]]] = {}` (line 550).
  * `_registry_base(cls) -> type` returns `CellUpdateBase` (lines 555-557).
* Registration key convention: `key="..."` class-creation kwarg.
  Example (from line 539): `class Step(CellUpdateBase, key="step"): ...`.
  Lookup: `CellUpdateBase.create("step")` (line 547).
* Abstract method: `update(...)` (lines 559-567) — same signature as
  the Protocol.
* Optional capability hook: `update_batch(slice_args: SweepCellSlice)
  -> np.ndarray` (lines 569-601). Default raises
  `NotImplementedError`. Strategies that participate in the 2-D
  Cartesian wavefront sweep override it (DD does — see §1.4 below).

### Companion dataclasses

* `CellVisit` (lines 278-347):
  * `cell_idx: int`.
  * `streaming_terms: StreamingTerms`.
  * `face_area_downstream: float | None` — sweep-direction-resolved
    outgoing face area. `None` for slab and for cylindrical-degenerate.
* `UpstreamState` (lines 354-376):
  * `spatial_upstream: np.ndarray` (ng,).
  * `angular_upstream: np.ndarray | None` (ng,) — `None` for slab.
* `CellResult` (lines 383-411):
  * `cell_average_flux: np.ndarray` (ng,).
  * `outgoing_spatial_flux: np.ndarray | None` — `None` for the
    cylindrical-degenerate case where no radial face flow exists.
  * `outgoing_angular_state: np.ndarray | None` — `None` for slab.
* `SweepCellSlice` (lines 172-271) — per-topological-level batched
  packet for the 2-D wavefront sweep. Fields documented at line 193.

### Concrete implementations

| Class | File | Line | Key | Notes |
|---|---|---|---|---|
| `DiamondDifference` | `orpheus/sn/spatial/diamond.py` | 280 | `"diamond_difference"` | Single geometry-polymorphic class covering slab + sphere + cylinder + cylindrical-degenerate; bit-identical to the inlined sweep math; also implements `update_batch` (2-D wavefront) at line 381. |

A `grep -rn 'class .*\(CellUpdate'` across `orpheus/` and `tests/`
returned only the Protocol/ABC themselves and `DiamondDifference`.
Test-only mocks: `IdentityCellUpdate` and `BadCellUpdate` in
`tests/sn/spatial/test_cell_update_protocol.py`. No other production
concretes — Step / LinearDiscontinuous / ExponentialCharacteristic are
documented in the diamond module docstring (lines 208-219) as
*planned* Wave-C-extension work but unimplemented.

### What `DiamondDifference` carries vs derives

* Carries (class-level): `is_linear=True`, `is_positivity_preserving=False`
  (lines 297, 303).
* Reads per-call from `visit.streaming_terms`: `alpha_in`, `abs_mu`,
  `chord_length`, `face_area_inner`/`_outer`, `delta_A_over_w`,
  `alpha_out`, `tau_mm`, `volume`.
* Three internal branches (each is a `@staticmethod`):
  * `_update_slab` (line 470) — chord-length DD recurrence.
  * `_update_curvilinear` (line 528) — sphere / non-degenerate
    cylinder via the M-M weighted DD balance. **Post-Step-2 uses
    `cell_balance_terms` from `cell_balance.py`** (line 567 imports
    it). Returns `outgoing_angular_state` from the M-M closure
    `(psi_avg - (1-tau)·psi_angle_in) / tau` (line 588).
  * `_update_cylindrical_degenerate` (line 598) — also routed
    through `cell_balance_terms_degenerate` post Step 2 (line 618).
* All branches mirror operation order from
  `orpheus/sn/sweep.py:_sweep_1d_{cumprod,spherical,cylindrical}`
  — see §6 for the procedural twin-path problem this creates.

---

## 2. The per-cell / per-level helpers

### `PsiHalfAngleSeed` Protocol (Phase D)

* File: `orpheus/sn/spatial/psi_half_angle_seed.py:228-271`.
* Decoration: `@runtime_checkable` (line 227).
* Class-level trait: `is_linear: bool` (line 264).
* Method signature (lines 266-271):
  ```python
  def __call__(
      self,
      psi_level: np.ndarray,             # (ng, M, nx)
      context: CarlsonSweepContext,
  ) -> np.ndarray:                       # (ng, nx)
      ...
  ```
* What it produces: the **M-M recurrence half-angle seed
  φ_{1/2,i,g}** — the input to
  `_mm_weighted_angular_recurrence_single_level`
  (`pole_angular_closure.py:340`) at `m=0`.
* Companion dataclass `CarlsonSweepContext` (lines 167-220):
  bundles `sigma_t`, `dr`, `mu_quad`, `weights`, `bc_outer_value`.

### `PsiHalfAngleSeedBase` ABC

* File: `orpheus/sn/spatial/psi_half_angle_seed.py:279-312`.
* Inherits: `RegistryMixin, ABC`.
* Same self-registration pattern via `key=...`.

### Concrete implementations

| Class | File | Line | Key | Behaviour |
|---|---|---|---|---|
| `ZeroSeed` | `psi_half_angle_seed.py` | 321 | `"zero"` | Phase B regression-safety ablation — returns `zeros((ng, nx))`. Carries the ERR-026 anti-pattern by design. |
| `CarlsonInwardSweep` | `psi_half_angle_seed.py` | 422 | `"carlson_inward_sweep"` | Phase D canonical — runs Hébert §3.9.4 Eqs. (3.432)-(3.435) inward μ=−1 sweep. Delegates the recurrence to `carlson_inward_sweep_from_source` (line 363). |

### `carlson_inward_sweep_from_source` (free function)

* File: `orpheus/sn/spatial/psi_half_angle_seed.py:363-419`.
* Signature: `(Q_bar, sigma_t, dr, bc_outer_value) -> ndarray` of
  shape `(ng, nx)`.
* Purpose: factored out so the SI/sweep path (Phase F backport) and
  the Krylov/matvec path consume the **same** Hébert (3.434)-(3.435)
  recurrence kernel. This is one of the few primitives the codebase
  has already correctly de-duplicated.

### `MorelMontryAngularSweep` (PoleAngularClosure concrete)

* File: `orpheus/sn/spatial/pole_angular_closure.py:466-661`.
* Self-registers via `key="morel_montry_angular_sweep"` (line 468).
* Field: `psi_half_seed: PsiHalfAngleSeed =
  field(default_factory=CarlsonInwardSweep)` (line 534). **This is
  the canonical composition pattern**: a concrete strategy carries a
  Protocol-typed field as an inner strategy, default-instantiated to
  the canonical Phase D choice.
* Class-level trait: `is_linear: ClassVar[bool] = True` (line 529).
* Method signature (`__call__`, lines 551-559):
  ```python
  def __call__(
      self,
      psi_cells: np.ndarray,
      alpha_half: "np.ndarray | list[np.ndarray]",
      redist_dAw: "np.ndarray | list[np.ndarray]",
      tau_mm: "np.ndarray | list[np.ndarray]",
      volume: np.ndarray,
      level_indices: "list[np.ndarray] | None" = None,
      carlson_context: "CarlsonSweepContext | list[CarlsonSweepContext] | None" = None,
  ) -> np.ndarray:
      ...
  ```
* What it computes: the per-(group, ordinate, cell) angular
  redistribution term `R_{n,i,g}` that the matvec adds to streaming +
  collision. Dispatches by `level_indices`: `None` = sphere (single
  level), populated = cylindrical (per-level loop). Delegates to
  `_mm_weighted_angular_recurrence_single_level` (line 340) which is
  also the algorithmic core consumed indirectly by the SI sweep
  (through the τ=½ branch of the cell update).
* Who calls it: the apply-matvec path
  (`transport_operator_matvec_spherical`/`_cylindrical`); the
  uncommitted `AngularRedistribution.apply` wrapper at
  `operators.py:433` (which is the architecture-incorrect Step 1
  promotion — see §6).

### `PoleAngularClosure` Protocol

* File: `orpheus/sn/spatial/pole_angular_closure.py:179-276`.
* Decoration: `@runtime_checkable`.
* Class-level trait: `is_linear: bool` (line 264).
* Method signature (lines 266-276): same as
  `MorelMontryAngularSweep.__call__` above (the strategy is the
  default registry entry for the Protocol).
* Three concrete subclasses:
  * `MorelMontryAngularSweep` (canonical; line 466) —
    `key="morel_montry_angular_sweep"`.
  * `BaileyFlatFluxRedist` (legacy ablation; line 669) —
    `key="bailey_flat_flux_redist"`. Reproduces the pre-Phase-B
    flat-flux collapse. Carries Defect 3 by design.
  * `LegacyTauSymmetricInterpolation` (line 812) — `key="legacy_tau_symmetric"`.
    Pre-Phase-B inlined τ-symmetric form. Bit-identical regression
    preservation.

### `_mm_weighted_angular_recurrence_single_level` (private helper)

* File: `orpheus/sn/spatial/pole_angular_closure.py:340-458`.
* Signature: `(psi_level, alpha_level, dAw_level, tau_level, volume,
  psi_half_seed=None) -> ndarray`.
* The algorithmic core of `MorelMontryAngularSweep` AND the cell
  update's M-M closure (via `DiamondDifference` at τ=½). When the
  apply-matvec and the SI sweep "solve the same discrete fixed
  point" (line 357), this is the load-bearing function — but it is
  invoked from **two unrelated call sites** today, which is the root
  of Manifestation #7.

---

## 3. The angular recurrence "Protocol" structure recap

There is **one Protocol** for the angular recurrence work:

* `PsiHalfAngleSeed` (the seed strategy, called at `m=0`).

The recurrence itself is **not** a Protocol — it is a free function
(`_mm_weighted_angular_recurrence_single_level`). The thing dressed
as a strategy at the level above it is `PoleAngularClosure`, which
COMPOSES the seed via a `field(default_factory=CarlsonInwardSweep)`
on `MorelMontryAngularSweep`. The composition is documented in the
`psi_half_angle_seed` module docstring at lines 103-117 as
"Option α (composition, not sibling Protocol)". This is the design
the user has already validated — Phase G should extend it, not
break it.

The Protocol shape `PsiHalfAngleSeed` is shared by `ZeroSeed`
(regression-safety; carries the ERR-026 anti-pattern) and
`CarlsonInwardSweep` (Phase D canonical). Both are linear; the
Protocol's `is_linear: bool` class-level trait is the load-bearing
property the operator algebra above (matvec linearity, Krylov
correctness) relies on.

---

## 4. `cell_balance.py` — Step 2 partial work (the salvageable piece)

* File: `orpheus/sn/spatial/cell_balance.py` (untracked / current
  working tree — 235 lines).
* Module docstring frames the intent (lines 1-40): the DD
  `_update_curvilinear` solve branch and `SNCellOperator`'s apply
  residual branch were duplicating the same algebra; Pattern 2
  (Cardinal Rule 2, "no twin paths") required factoring it to one
  place.

### Surface

| Symbol | Lines | Role |
|---|---|---|
| `CellBalanceTerms` (frozen dataclass) | 53-87 | Algebraic intermediates `denom`, `numer_upstream`, `c_in`, `c_out` |
| `cell_balance_terms(st, A_downstream, total_xs, upstream_state)` | 90-182 | Non-degenerate curvilinear branch. Mirrors `sweep.py:328-355` operation order for bit-identity (line 109). |
| `cell_balance_terms_degenerate(st, total_xs, upstream_state)` | 185-227 | Cylindrical pure-azimuthal degenerate variant. |

### Who consumes it (today, in the uncommitted working tree)

* `orpheus/sn/spatial/diamond.py:567` — `_update_curvilinear` imports
  `cell_balance_terms` and rebuilds DD's cell-average solve as
  `psi_avg = (source + terms.numer_upstream) / terms.denom`. Lines
  583-588 then run the WDD spatial closure + M-M angular closure on
  top.
* `orpheus/sn/spatial/diamond.py:618` — `_update_cylindrical_degenerate`
  imports `cell_balance_terms_degenerate` and does the same.
* `orpheus/sn/spatial/operators.py:319-341` —
  `_apply_curvilinear_residual` consumes the same helper to produce
  `terms.denom · cell_avg − (source + terms.numer_upstream)` as the
  per-cell residual.
* `orpheus/sn/spatial/operators.py:344-362` —
  `_apply_cylindrical_degenerate_residual` analogously.

**Verdict: keep `cell_balance.py` as is.** It is the **correct**
Step 2 fix — one place for the cell-balance algebra, two consumers
(solve via DD, apply via residual). The architecturally wrong piece
sits **above** it, not in `cell_balance.py` itself.

---

## 5. `BoundaryOperator` / `BoundaryRealizer` strategy Protocol

### What is `BoundaryOperator` today?

* **Deprecated alias** for `BoundaryTraceLaw`. See
  `orpheus/geometry/boundary/__init__.py:420`:
  `BoundaryOperator = BoundaryTraceLaw`.
* `BoundaryTraceLaw` (`_base.py:74-298`) is a `RegistryMixin, ABC`
  descriptor — **NOT a `LinearOperator`**. The `_base.py` module
  docstring is explicit (lines 18-22): "It is NOT a callable
  operator — there is no `apply` method and no
  `LinearOperatorMixin` inheritance. The realizer is the sole
  bridge from descriptor to operator."
* Carries three first-class properties: `geometry_map: G`,
  `response_kernel: R`, `source: q` for the affine boundary
  ψ_− = R·G·ψ_+ + q.
* Carries `creates_sweep_cycle: ClassVar[bool]` (line 143).
* Has an algebra (lines 178-214): `+`, `-`, `*`, `/`, unary `-`
  build `LawSum`/`LawScaled` *descriptor-tree* nodes. The
  *operator-tree* algebra (`OperatorSum`/`ScaledOperator`) is a
  separate type family — composing the two is a static-type error.
* Registry key convention: ``key="vacuum"``, `"reflective"`,
  `"white"`, `"periodic"`, `"albedo"`, `"prescribed_inflow"`.

### `BoundaryRealizer` Protocol (the law → operator bridge)

* File: `orpheus/geometry/boundary/_realizer.py:84-126`.
* Decoration: `@runtime_checkable`.
* Class-level: `method_name: str`.
* Method: `realize(law, method_space) -> LinearOperator`.
* Registry: `BoundaryRealizerRegistry` (lines 129-198) — stand-alone
  (NOT mounted on `RegistryMixin`) because realizers are
  per-transport-method strategies, not a subclass hierarchy.

### `SNBoundaryRealizer` (the only functional realizer)

* File: `orpheus/sn/boundary_realizer.py:106-200`.
* Registered via `@BoundaryRealizerRegistry.register("SN")`
  (line 105).
* Stateless: `method_name: str = "SN"` is the only attribute.
* `realize(law, method_space) -> LinearOperator` (line 122)
  dispatches by `isinstance(law, ...)` to map each BC subclass onto
  a Wave-0 / Wave-1 primitive:
  * `VacuumBoundaryOperator` → `IncomingOrdinateMaskTensor` (line 143).
  * `SpecularBoundaryOperator` → `PermutationOperator` (α=1) /
    `ScaledOperator(α, PermutationOperator)` (line 149).
  * `WhiteBoundaryOperator` → `AngularAverageOperator.from_quadrature`
    (line 158).
  * `AlbedoBoundaryOperator` → `ZeroOperator` (α=0) /
    `IdentityOperator` (α=1) / `ScaledOperator(...)` (line 166).
  * `PeriodicBoundaryOperator` → `PeriodicWrapOperator` (line 173).
  * `PrescribedInflow` → `IncomingSourceOperator(law.source)` (line 176).
* The realized output IS a `LinearOperator` — it satisfies the
  Protocol at `orpheus/numerics/operator.py:115-183` and inherits
  the algebra dunders via `LinearOperatorMixin`
  (`orpheus/numerics/operator.py:202`).

### How does SN consume a BC today?

The descriptor is in the law registry; the realized
`LinearOperator` is bound onto `_BoundBoundaryOperator`
(`_bound_compat.py`, internal) and called as `bc_op.apply(psi)` on
the boundary trace — **once per sweep**, on the trace vector, not
per-ordinate per-cell. The per-ordinate slicing inside the cell
loop went away in Wave 9 (single-arg migration) and Wave 11 (mixed
removal). So the BC composition into the SN sweep is **already at
the trace-edge layer**, not inside the per-cell loop — the
realizer pattern works because the operator and the trace match in
type at the call site.

---

## 6. The uncommitted Step 1 work — what to keep, what to delete

### `orpheus/sn/spatial/operators.py` (511 lines, uncommitted)

The module docstring (lines 1-56) is explicit: "**pure type-system
promotion**" — wrapping two existing strategies as
`LinearOperatorMixin` subclasses to "make the operator algebra (+,
-, @, .H, capability sets, dunder composition) available at the
per-cell and per-level layers".

#### `SNCellOperator(LinearOperatorMixin)` — lines 84-313

* Field: `cell_update: DiamondDifference = field(default_factory=DiamondDifference)`
  (line 148).
* Capabilities: `frozenset({CAP_APPLY, CAP_SOLVE})` (line 139).
* `solve(...)` (line 150) — **pure delegation** to
  `self.cell_update.update(visit, total_xs, source, upstream_state)`
  (line 187).
* `apply(cell_avg, *, visit, total_xs, upstream_state, source)`
  (line 189) — the genuinely new piece. Computes the per-cell
  residual `L_cell · ψ̄ − q` by:
  * Slab branch (lines 276-290): closed-form residual
    `2|μ|·(cell_avg − ψ_in) + chord·Σ_t·cell_avg − source`.
  * Curvilinear non-degenerate (lines 301-313): delegates to
    `_apply_curvilinear_residual` (line 319) which itself consumes
    `cell_balance_terms` from `cell_balance.py`.
  * Cylindrical-degenerate (lines 293-299): delegates to
    `_apply_cylindrical_degenerate_residual` (line 344) which
    consumes `cell_balance_terms_degenerate`.

**Verdict.** The wrapping IS architecturally wrong: a per-cell
strategy is being lifted to `LinearOperator` for the sake of
acquiring algebra dunders it does not need (no SN cell-update ever
participates in `A + B` or `A @ B` at the cell level — the algebra
the codebase actually relies on lives **above** at the
streaming/collision/scattering/fission level). But the `.apply`
**body** — computing the residual via `cell_balance_terms` — is
load-bearing and correct. **Move that body into `CellUpdate`** as
a new method, e.g. `residual(visit, total_xs, source,
upstream_state, cell_avg) -> ndarray`. The Protocol grows by one
method; the wrapper class disappears entirely. The per-branch
private helpers `_apply_curvilinear_residual` /
`_apply_cylindrical_degenerate_residual` (lines 319, 344) become
either methods on `DiamondDifference` (mirroring the existing
`_update_curvilinear` / `_update_cylindrical_degenerate`
@staticmethods) or free functions in `cell_balance.py` that the
`residual` method calls.

#### `AngularRedistribution(LinearOperatorMixin)` — lines 370-505

* Field: `angular_sweep: MorelMontryAngularSweep =
  field(default_factory=MorelMontryAngularSweep)` (line 429).
* Capabilities: `frozenset({CAP_APPLY})` (line 425).
* `apply(...)` (line 433) — **pure delegation** to
  `self.angular_sweep(...)` (line 497).

**Verdict.** This wrapper adds **literally nothing**. No new
behaviour, no new method, no new field, no new capability —
`MorelMontryAngularSweep` already exposes the same callable
surface and already advertises `is_linear=True`. Promoting it to a
`LinearOperator` adds the `@`/`+` dunders but those compose only
with other operators on the same vector space; the M-M
redistribution doesn't share a vector space with streaming or
scattering (it operates on `(ng, M, nx)` flat-angular-flux, not on
the trace or the cell-centre flux). **Delete this class entirely.**

#### Private residual helpers

`_apply_curvilinear_residual` (line 319) and
`_apply_cylindrical_degenerate_residual` (line 344) — these are
the right helpers, just in the wrong module. They consume
`cell_balance.cell_balance_terms` correctly. **Move them into
`diamond.py`** as static methods on `DiamondDifference`, or into
`cell_balance.py` next to the term-builders. Either way, the
algebra is correct and shouldn't be re-derived.

---

## 7. Anti-recommendations for the Step 1 / Step 2 agent

1. **Do NOT create `SNCellOperator(LinearOperatorMixin)`.** The
   `CellUpdate` Protocol at `orpheus/sn/spatial/cell_update.py:418`
   is the right home for the per-cell strategy. Add a `.residual`
   method to that Protocol (and to `CellUpdateBase` as an
   `abstractmethod`); have `DiamondDifference` implement it by
   consuming `cell_balance.cell_balance_terms`. Do **NOT** introduce
   a wrapper class to "gain operator algebra" — `+`, `@`, `.H` are
   not used at the cell level anywhere in the SN module.

2. **Do NOT create `AngularRedistribution(LinearOperatorMixin)`.**
   The `PoleAngularClosure` Protocol at
   `orpheus/sn/spatial/pole_angular_closure.py:179` already gives
   `MorelMontryAngularSweep` the exact callable surface the matvec
   needs. No wrapper is required. If the apply-matvec needs an
   `is_linear` guarantee, read it off the strategy's
   `ClassVar[bool]` directly (line 529).

3. **Do NOT introduce a `BoundaryOperator` realiser parallel to
   `SNBoundaryRealizer`.** `orpheus/sn/boundary_realizer.py:106`
   already maps every `BoundaryTraceLaw` subclass to a Wave-0 /
   Wave-1 `LinearOperator`. Per Cardinal Rule 2, if BC composition
   into the cell loop is needed, route through
   `SNBoundaryRealizer().realize(law, method_space)` and consume the
   returned `LinearOperator` — do NOT bypass the realizer.

4. **Do NOT touch `cell_balance.py`.** The untracked
   `orpheus/sn/spatial/cell_balance.py` (235 lines) is the correct
   Step 2 deliverable: one place for the curvilinear cell-balance
   algebra (`CellBalanceTerms`, `cell_balance_terms`,
   `cell_balance_terms_degenerate`), two consumers (DD's
   `_update_curvilinear` / `_update_cylindrical_degenerate` solve
   branches, and the apply-residual sites). Step 2's twin-path fix
   IS done; the architectural problem sits in `operators.py`'s
   wrapper layer, not in `cell_balance.py`.

5. **Do NOT replace the `PsiHalfAngleSeed` Protocol with a
   `LinearOperator` Protocol.** The Protocol at
   `orpheus/sn/spatial/psi_half_angle_seed.py:228` is intentionally
   narrow — `__call__(psi_level, context) -> ndarray` — and
   compositional (the `MorelMontryAngularSweep` carries one via
   `field(default_factory=CarlsonInwardSweep)` at
   `pole_angular_closure.py:534`). This composition pattern is the
   correct extension model; emulate it for any new inner-strategy
   slot.

6. **Do NOT inline the M-M recurrence kernel.** The function
   `_mm_weighted_angular_recurrence_single_level` at
   `orpheus/sn/spatial/pole_angular_closure.py:340` is the single
   source of truth for the per-cell M-M weighted DD recurrence.
   Both `MorelMontryAngularSweep` and (transitively, via the τ=½
   substitution and the cell-update closure) `DiamondDifference`'s
   curvilinear branch must remain bit-consistent with it. Phase G's
   Manifestation #7 (SI-vs-Krylov O(h) drift) is exactly the cost of
   having two procedural twins of these primitives — DO NOT add a
   third.

7. **Do NOT promote `BoundaryTraceLaw` to a `LinearOperator`.** The
   §16A.3 three-layer split (descriptor / realizer / operator) is
   enforced by the type system at
   `orpheus/geometry/boundary/_base.py:18-29`. The deprecated
   `BoundaryOperator = BoundaryTraceLaw` alias at
   `orpheus/geometry/boundary/__init__.py:420` exists for backward
   compatibility — it is **not** an invitation to re-merge the
   layers.

8. **Do NOT re-derive `carlson_inward_sweep_from_source`.** The free
   function at `orpheus/sn/spatial/psi_half_angle_seed.py:363` is
   the SI-and-Krylov-shared kernel for Hébert (3.434)-(3.435). It
   IS the correct architecture for de-duplicating the seed
   recurrence; any new path that needs an inward μ=−1 sweep must
   call this function, not transcribe the recurrence.

9. **Do NOT add an SI vs Krylov branching switch.** The Phase G
   goal (per the agent-memory issue-196 audit) is to compose ONE
   operator stack out of LinearOperator types so SI sweep, Krylov
   matvec, and BC face-trace algebra all flow through the same
   primitives. Adding `if mode == "si": ... else: ...` is the
   anti-pattern Manifestation #7 already exhibits.

10. **Do NOT advertise a capability the operator cannot deliver.**
    The `LinearOperator` Protocol at
    `orpheus/numerics/operator.py:115` uses
    `capabilities: frozenset[str]` as the **single source of
    truth** for `apply` / `solve` / `apply_transpose` availability.
    The Step 1 wrappers correctly omit `CAP_APPLY_TRANSPOSE` and
    `AngularRedistribution` correctly omits `CAP_SOLVE`. Phase G
    extensions MUST keep this discipline.

---

## 8. Inventory summary (cross-reference table)

| Protocol | File:line | Concrete ABC (file:line) | Concrete strategies (file:line, key) |
|---|---|---|---|
| `CellUpdate` | `orpheus/sn/spatial/cell_update.py:418` | `CellUpdateBase` (`cell_update.py:510`) | `DiamondDifference` (`diamond.py:280`, `"diamond_difference"`) |
| `PsiHalfAngleSeed` | `orpheus/sn/spatial/psi_half_angle_seed.py:228` | `PsiHalfAngleSeedBase` (`psi_half_angle_seed.py:279`) | `ZeroSeed` (`psi_half_angle_seed.py:321`, `"zero"`), `CarlsonInwardSweep` (`psi_half_angle_seed.py:422`, `"carlson_inward_sweep"`) |
| `PoleAngularClosure` | `orpheus/sn/spatial/pole_angular_closure.py:179` | `PoleAngularClosureBase` (`pole_angular_closure.py:284`) | `MorelMontryAngularSweep` (`pole_angular_closure.py:466`, `"morel_montry_angular_sweep"`), `BaileyFlatFluxRedist` (`pole_angular_closure.py:669`, `"bailey_flat_flux_redist"`), `LegacyTauSymmetricInterpolation` (`pole_angular_closure.py:812`, `"legacy_tau_symmetric"`) |
| `BoundaryRealizer` | `orpheus/geometry/boundary/_realizer.py:84` | (no ABC; stand-alone registry `BoundaryRealizerRegistry` at `_realizer.py:129`) | `SNBoundaryRealizer` (`orpheus/sn/boundary_realizer.py:106`, registered as `"SN"`); 4 stub realizers in MoC/MC/CP/Diffusion |
| `BoundaryTraceLaw` (descriptor, not strategy) | `orpheus/geometry/boundary/_base.py:74` | (IS the ABC) | `VacuumInflow`, `ReflectiveBoundary`, `WhiteBoundary`, `PeriodicBoundary`, `AlbedoBoundary`, `PrescribedInflow` (each in its own submodule of `orpheus/geometry/boundary/`) |
| `LinearOperator` | `orpheus/numerics/operator.py:115` | `LinearOperatorMixin` (`operator.py:202`) | Many (see `operator.py`); `SNCellOperator` / `AngularRedistribution` at `orpheus/sn/spatial/operators.py:84,370` are the misplacements to delete. |

Note that `BoundaryTraceLaw` is in this table for completeness, but
it is **NOT a strategy Protocol** in the same sense as the SN
strategies — it is a descriptor that the realizer turns into a
`LinearOperator`. The strategy Protocol on that side is
`BoundaryRealizer`.

---

## 9. Orphan / dead / partial abstractions to flag for replan

* **`SNCellOperator` and `AngularRedistribution`** at
  `orpheus/sn/spatial/operators.py:84,370` — uncommitted Step 1
  misplacements. `AngularRedistribution` is pure delegation and adds
  nothing. `SNCellOperator.apply` carries a load-bearing residual
  body that must be ported to `CellUpdate.residual`; `solve` is pure
  delegation and disappears. The whole file `operators.py` likely
  goes away after Step 2 cleanup; the private helpers
  `_apply_curvilinear_residual` /
  `_apply_cylindrical_degenerate_residual` (lines 319, 344) belong
  in `cell_balance.py` or as DD static methods.

* **`BoundaryOperator` deprecated alias** at
  `orpheus/geometry/boundary/__init__.py:420` — flagged by the
  comment "remove in a future cleanup wave". Phase G should not
  expand its usage; new sites must consume `BoundaryTraceLaw`
  directly.

* **`PoleAngularClosure` `_registry_base` confusion** — the Nexus
  query returned an `unresolved` node
  `py:method:PoleAngularClosureBase.create("morel_montry_angular_sweep")`
  with degree 1, suggesting a documentation reference that the
  graph parser couldn't resolve to a concrete factory call. The
  `create` classmethod itself comes from `RegistryMixin` (not
  `_realizer.py`) and is exposed through every `*Base` ABC. Not a
  bug, but worth confirming the Sphinx cross-refs build cleanly
  before Phase G.

* **`update_batch` only on `DiamondDifference`** — the `CellUpdateBase`
  ABC ships a default `update_batch` that raises
  `NotImplementedError` (line 597). Only `DiamondDifference`
  overrides it (line 381). The 2-D Cartesian wavefront sweep
  pathway will silently break if a non-DD strategy is dispatched to
  `_sweep_2d_wavefront`. Phase G shouldn't surface this — but if
  it ships a second `CellUpdate` concrete, the `update_batch`
  status must be addressed.

* **The `psi_half_angle_seed` Phase D module** is a model
  implementation of the user's "unify after two instances"
  preference: `ZeroSeed` (regression-safety) + `CarlsonInwardSweep`
  (canonical), unified through Protocol + ABC + RegistryMixin only
  AFTER both ships shipped. Phase G should preserve this pattern.

---

## 10. Bottom line for the Phase G replanner

* **Existing strategy Protocols** (extend, don't replace):
  `CellUpdate`, `PsiHalfAngleSeed`, `PoleAngularClosure`,
  `BoundaryRealizer`. Each ships a `runtime_checkable` Protocol + a
  `*Base(RegistryMixin, ABC)` concrete ABC + at least one
  registered concrete.
* **Salvageable from Step 1 / Step 2**: `cell_balance.py` (keep,
  correct); the `apply` residual body inside `SNCellOperator` (move
  to a new `CellUpdate.residual` method); the private residual
  helpers (move to `cell_balance.py` or `DiamondDifference`).
* **Delete entirely**: `AngularRedistribution` (pure delegation, no
  value); `SNCellOperator` wrapper (the `solve` is pure delegation;
  the `apply` body migrates to `CellUpdate.residual`).
* **Architectural Frame Phase G must obey**:
  1. Strategies are NOT LinearOperators. Strategies live at the
     per-cell / per-level layer where vector-space operator algebra
     is not the right abstraction.
  2. The LinearOperator algebra lives ABOVE strategies, at the
     streaming/collision/scattering/fission layer where the
     operands DO share a vector space.
  3. Strategies compose via Protocol-typed fields with default
     factories (as `MorelMontryAngularSweep.psi_half_seed`
     demonstrates), not by inheriting algebra mixins.
  4. Sweep / matvec twins (Manifestation #7) close by having BOTH
     paths consume the same underlying free function or
     strategy — as `carlson_inward_sweep_from_source` and
     `cell_balance_terms` already demonstrate at the layers below.

---

*References*

* Architecture audit:
  `.claude/agent-memory/explorer/issue_196_sn_operator_architecture_audit.md`
* Phase G plan:
  `.claude/plans/issue_196_phase_g_four_operator_unification.md`
* Issue #196 on GitHub.
* `coding-elegance` skill — `.claude/skills/coding-elegance/SKILL.md`.
* `subagent-handoff-protocol` skill — used by future replanner.
