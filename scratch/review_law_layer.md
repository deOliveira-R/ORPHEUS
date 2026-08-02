# Boundary review — quadrant 3: the LAW layer

**Scope:** `orpheus/geometry/boundary/` — descriptors, invariants, admission,
composition. Companion to `.claude/plans/boundary_machinery_review.md` §8.3.
**Branch:** `refactor/operator-strategy-layers` @ `73627b71`.
**Evidence key:** `[M]` measured · `[R]` read (quoted, `file:line`) · `[G]`
exhaustive grep · `[U]` unverified (must not enter the plan).

---

## 1. Complete law inventory

**[G] The set of seven is COMPLETE and CORRECT.** Exhaustive:
`grep -rn --include='*.py' -E "^\s*class\s+\w+\(.*BoundaryTraceLaw" .`
returns exactly 8 hits outside `.claude/worktrees/` — the 7 production laws
below plus ONE test stub (`tests/geometry/test_boundary_trace_law.py:53`
`class _StubLaw(BoundaryTraceLaw, key="_stub_for_test")`). There is no
abstract intermediate: every concrete law inherits `BoundaryTraceLaw`
directly, and every one carries a `key=`.

| # | class | `file:line` | registry `key` | fields (defaults) |
|---|---|---|---|---|
| 1 | `VacuumInflow` | `vacuum.py:20` | `"vacuum"` | `kind: str = "vacuum"` |
| 2 | `ReflectiveBoundary` | `reflective.py:41` | `"reflective"` | `axis: str = "x"`, `albedo: float = 1.0` |
| 3 | `WhiteBoundary` | `white.py:26` | `"white"` | `axis: str = "x"`, `outward_sign: int = +1`, `albedo: float = 1.0` |
| 4 | `AlbedoBoundary` | `albedo.py:26` | `"albedo"` | `albedo: float = 0.0` |
| 5 | `PeriodicBoundary` | `periodic.py:20` | `"periodic"` | **none** (stateless) |
| 6 | `PrescribedInflow` | `prescribed_inflow.py:46` | `"prescribed_inflow"` | `_source: InflowSourceSpec = NoSource()` (public kwarg `source=`) |
| 7 | `ZeroFluxBoundary` | `zero_flux.py:24` | `"zero_flux"` | **none** (stateless) |

### 1.1 What each parameter carries, physically

**`VacuumInflow`** — **[R]** `vacuum.py:21` *"Vacuum boundary: :math:`R = 0`,
:math:`q = 0`. The empty sum in the tensor decomposition."* The `kind` field is
**not physics** — **[R]** `vacuum.py:59-61`: *"String tag for legacy
string-kind comparisons. Preserved across the Issue #186 descriptor cleanup so
the SN-side `sn_mesh.bc["xmax"] == "vacuum"` test contract still holds."* It is
a dataclass FIELD here (mutable-by-construction, `kind: str = "vacuum"`),
unlike the other laws which expose `kind` as a read-only `@property`. That is
an inconsistency inside the family: `VacuumInflow(kind="banana")` constructs
without complaint. **[M]** verified — see §1.3.

**`ReflectiveBoundary.axis`** — the reflection axis; selects which column of
`Quadrature.axis_cosines` the mirror `Ω ↦ Ω − 2(Ω·n̂)n̂` acts on. **[R]**
`reflective.py:88-92`. Mapped to a column index by the module-private
`_AXIS_INDEX = {"x": 0, "y": 1, "z": 2}` (`reflective.py:30`).
**`.albedo`** — the SPECULAR albedo `α ∈ [0,1]`, the response amplitude on the
mirrored ordinate (**[R]** `reflective.py:47-49`). Note `albedo` here also
drives the `kind` tag: **[R]** `reflective.py:105-107`
`return "reflective" if self.albedo == 1.0 else "partial"` — an exact float
comparison, so `albedo=0.9999999999` tags as `"partial"`.

**`WhiteBoundary.axis`** — the boundary NORMAL axis. **`.outward_sign`** — `+1`
upper face / `−1` lower face; **[R]** `white.py:75-79`: *"Selects which
ordinates are *outgoing* at this face."* This is the one law that carries its
own face-side orientation as a field rather than deriving it from the face
label. **`.albedo`** — the DIFFUSE albedo scaling the returned Lambertian
current. **[R]** `white.py:44-48`: the `1/π` Lambertian constant is replaced in
implementation by *"normalis[ing] by the outgoing cosine-weighted weight sum so
the BC is conservative for any quadrature."*

**`AlbedoBoundary.albedo`** — a pure scalar attenuation `ψ_in(Ω) = α ψ_out(Ω)`
with NO angular redistribution: the tensor pair is `(I, α)`. **[R]**
`albedo.py:29-34`. Default `0.0` — i.e. the default `AlbedoBoundary()` IS
vacuum, which is a defaulting choice worth noting (the other albedo-bearing
laws default to `1.0`).

**`PeriodicBoundary`** — stateless. Its geometric content is the SPATIAL wrap
`ψ_in(Ω, x_this) = ψ_out(Ω, x_partner)` (**[R]** `periodic.py:29-33`), and the
partner face is **not a field on the law** — **[R]** `periodic.py:42-43`: *"The
two-face plumbing is handled by whoever instantiates `PeriodicBoundary` and
orchestrates the sweep."* This is the law whose defining parameter is missing.

**`PrescribedInflow.source`** — an `InflowSourceSpec` (a `runtime_checkable`
Protocol with one method, `evaluate(shape) -> np.ndarray`; `_source.py:45`).
Two shipped implementations: `NoSource` (`_source.py:68`, zeros) and
`ConstantInflowSource(value: float)` (`_source.py:83`). The field is named
`_source` with a `source` property override to dodge a property/field collision
— **[R]** `prescribed_inflow.py:86-94` documents exactly this, and the class
hand-writes `__init__` (`prescribed_inflow.py:99`) rather than using the
dataclass-generated one.

**`ZeroFluxBoundary`** — stateless; `φ_Γ = 0`. Its physics is the P1 closure
`φ_Γ = 2(J⁺+J⁻) ⟹ 𝒜 = −1` — **[R]** `zero_flux.py:41-47`. It is
**deliberately outside the sub-Markov family** (**[R]** `zero_flux.py:49-58`):
*"zero-flux requires a NEGATIVE incoming partial current, which no non-negative
angular flux can produce… this class does not override
`assert_response_positive_if_declared` to raise, because the negative response
IS its definition."* Diffusion-only; SN refuses it (`zero_flux.py:64-71`).

### 1.2 The `kind` inconsistency across the family

**[G]** Three different spellings of the same concept coexist:

| law | how `kind` is provided |
|---|---|
| `VacuumInflow` | dataclass **field** `kind: str = "vacuum"` (`vacuum.py:61`) |
| `ReflectiveBoundary` | `@property` derived from `albedo` (`reflective.py:105`) |
| `PrescribedInflow` | `@property`, constant (`prescribed_inflow.py:114`) |
| `WhiteBoundary`, `AlbedoBoundary`, `PeriodicBoundary`, `ZeroFluxBoundary` | **NO `kind` at all** — they fall back to `RegistryMixin.key` (`registry.py:85`, a `ClassVar`) |

So `law.kind` is an `AttributeError` on four of seven laws, while `law.key` is
universal. Verified by probe in §1.3.

Two laws also hand-write `__eq__`/`__hash__` to compare equal to a STRING
(`vacuum.py:63-68`, `reflective.py:109-116`) — `VacuumInflow() == "vacuum"` is
`True`. The other five do not. This is the string-kind coupling the review's
§3.1 flags, embedded in the law layer's own `__eq__`.

### 1.3 [M] Runtime probe of the law surface

Probe: `.venv/bin/python` importing `orpheus.geometry.boundary`, then
`orpheus.sn.mesh.augmented_mesh` + `orpheus.diffusion.augmented_mesh`.

```
kind / key / properties per law                     (measured 2026-07-30)
  VacuumInflow          kind='vacuum'             key='vacuum'            gm=None rk=None src=NoSource
  ReflectiveBoundary    kind='reflective'         key='reflective'        gm=None rk=None src=NoSource
  ReflectiveBoundary(α=0.5) kind='partial'        key='reflective'        gm=None rk=None src=NoSource
  WhiteBoundary         kind=<<AttributeError>>   key='white'             gm=None rk=None src=NoSource
  AlbedoBoundary        kind=<<AttributeError>>   key='albedo'            gm=None rk=None src=NoSource
  PeriodicBoundary      kind=<<AttributeError>>   key='periodic'          gm=None rk=None src=NoSource
  PrescribedInflow      kind='prescribed_inflow'  key='prescribed_inflow' gm=None rk=None src=NoSource
  ZeroFluxBoundary      kind=<<AttributeError>>   key='zero_flux'         gm=None rk=None src=NoSource
```

**[M]** `VacuumInflow(kind="banana")` **constructs without complaint**, and then
`v == "vacuum"` is `False` while `v == "banana"` is `True`. The `kind` field is
a settable constructor parameter on a frozen dataclass whose whole purpose is
to be the immutable tag. This is an illegal state that IS representable.

**[M]** `geometry_map` and `response_kernel` are `None` on **all seven** laws —
independently confirming the review's §3 `[G]`.

---

## 2. Does each law factor as `R ∘ G`?

**The grand report's definitions** (**[R]** `.claude/plans/neutron_transport_grand_report_v3.md:2780-2786`):

```python
class BoundaryGeometryMap(Protocol):
    def map_outgoing_to_incoming_geometry(self, outgoing_trace_point): ...
    def as_tensor(self, trace_space): ...

class BoundaryResponseKernel(Protocol):
    def apply(self, mapped_outgoing_trace): ...
    def as_tensor(self, trace_space): ...
```

So `G` is the **relabeling** (which outgoing trace point feeds which incoming
one) and `R` is the **amplitude/kernel applied after relabeling**. Both
Protocols carry `as_tensor(trace_space)` — i.e. **the design already
acknowledges that the discretized form is produced by handing the spec a trace
space**, not by the law knowing it. That is the decisive architectural fact for
this question: the law owns the *spec*, the trace space owns the *matrix*.

### 2.1 The per-law factorization table

| law | (a) geometric content `G` | (b) response content `R` | (c) derivable from the law's own fields alone? |
|---|---|---|---|
| `VacuumInflow` | identity (annihilated by `R=0`) | `0` | **BOTH YES.** Zero-parameter. Grand report spells it exactly: `geometry_map=IdentityBoundaryMap(), response=ZeroBoundaryResponse()` (`grand_report:2810-2814`) |
| `ReflectiveBoundary` | specular mirror `Ω ↦ Ω − 2(Ω·n̂)n̂` across `axis` | scalar `α = albedo` | **R: YES.** **G: as a typed spec YES** (`SpecularMirror(axis)` — `axis` is the only content); **as a realized matrix NO** — needs `quadrature.reflection_index(axis)`, an N-vector of the angular discretization |
| `WhiteBoundary` | cosine-weighted average over the `outward_sign` hemisphere of `axis`, broadcast isotropically | scalar `α = albedo` | **R: YES.** **G: as a typed spec YES** (`HemisphericalAverage(axis, outward_sign)`); **as a realized matrix NO** — needs the quadrature weights AND axis cosines, because the normalizer is `Σ_{outgoing} w_n |μ_n|` |
| `AlbedoBoundary` | **identity** | scalar `α = albedo` | **BOTH YES — even in realized form.** `IdentityOperator()` requires no discretization data at all. The unique law where the factorization is complete with zero external input |
| `PeriodicBoundary` | spatial wrap to the **partner face** | `1` | **R: YES** (constant). **G: NO — not even as a typed spec.** The partner face is *not a field on the law* (`periodic.py:51` — the class body is empty). `FaceWrap(partner)` is unconstructible from `PeriodicBoundary()` |
| `PrescribedInflow` | `0` (rank-0) | `0` | **BOTH YES** (both are zero by the law's own definition — **[R]** `prescribed_inflow.py:51`: *"The rank-0 case `R = G = 0`"*). The content is in `q` |
| `ZeroFluxBoundary` | identity | `−1` **in the P1 partial-current basis only** | **R: BASIS-DEPENDENT — see §2.3.** No angular `R` exists |

### 2.2 The precise statement of where discretization is needed

Two of seven laws (`ReflectiveBoundary`, `WhiteBoundary`) have a `G` that is
**law-intrinsic as a specification but quadrature-dependent as a matrix**. The
split is exact and it is the split the grand report's Protocol already encodes:

* `SpecularMirror(axis="x")` — constructible from `ReflectiveBoundary.axis`
  alone. Its `as_tensor(trace_space)` is `PermutationOperator(quad.reflection_index("x"))`.
* `HemisphericalAverage(axis="x", outward_sign=+1)` — constructible from
  `WhiteBoundary.axis` + `.outward_sign` alone. Its `as_tensor(trace_space)`
  is `AngularAverageOperator.from_quadrature(quad, axis, sign)`.

**No law needs mesh data for `R`.** Every `R` in the family is a SCALAR
(`0`, `α`, `1`, `−1`) — **[R]** confirmed at both realizer call sites, which
multiply by a bare float: `sn/boundary/realizer.py:288`
`return stamp_boundary_role(float(law.albedo) * base)` and
`diffusion/boundary_realizer.py:214` `return float(law.albedo)`.

> Consequence for the empty-property decision: `response_kernel` is
> **already fully populatable today** — it is `law.albedo` under a different
> name, and both production realizers already read it that way.
> `geometry_map` is populatable as a typed spec for **6 of 7** laws; the one
> blocker is `PeriodicBoundary`, and the blocker is a **missing field**, not a
> missing discretization.

### 2.3 `ZeroFluxBoundary` is not an affine trace law at all

`φ_Γ = 0` is not of the form `γ₋ψ = R G γ₊ψ + q` in the angular basis: it
constrains a functional of BOTH traces at once. In the grand report's own
taxonomy (**[R]** `grand_report:2749-2775`) that is the tier ABOVE:

> ```
> BoundaryTraceLaw          # affine law gamma_- = R G gamma_+ + q
> BoundaryRelation          # more general constraint A_- gamma_- + A_+ gamma_+ = g
> ```
> …
> ```
> BoundaryCondition
>     BoundaryRelation
>         AffineTraceLaw
>             ZeroInflowTrace / VacuumInflow / PrescribedInflow /
>             ReflectiveBoundary / PeriodicBoundary / AlbedoBoundary /
>             SymmetryBoundary
>         WeakBoundaryForm
>             MarshakBoundaryForm / MarkBoundaryForm / ...
> ```

`φ_Γ = 0` is `A₋γ₋ + A₊γ₊ = 0` with `A₋ = A₊ = 1` (via `φ = 2(J⁺+J⁻)`) — a
`BoundaryRelation`, not an `AffineTraceLaw`. It becomes affine ONLY after the P1
closure supplies the dictionary. **[G]** The code ships **no `BoundaryRelation`
tier**: `grep -rn "BoundaryRelation" orpheus/` returns nothing. Everything is
collapsed onto the single `BoundaryTraceLaw` ABC.

This is the structural reason SN must refuse `ZeroFluxBoundary` — and the
refusal is currently expressed as a hand-written `isinstance` guard in the SN
realizer (`sn/boundary/realizer.py:164-180`) rather than as a type distinction.
The tier that would make the refusal *unspellable* is the one that is missing.

### 2.4 The vacuum mismatch: `R = 0` but the realized operator is a projector

**[R]** Grand report `:2810-2814` — vacuum is `response=ZeroBoundaryResponse()`,
so `R·G = 0` and `B_vacuum` should be the **zero map on the trace block**.
**[R]** The SN realizer instead builds
`IncomingOrdinateMaskTensor(inflow_indices, …) & IdentityOperator()`
(`sn/boundary/realizer.py:256-263`) — a **projector** that zeroes the inflow
rows and *passes the outflow rows through*.

The two disagree on the outflow rows. The realizer's justification
(**[R]** `sn/boundary/realizer.py:24-28`) is *"the §16A.10 trace-correct
implementation zeroes ONLY the inflow ordinates so the outflow trace survives.
The new behaviour is the right one; downstream consumers (sensitivity adjoints)
need the outflow trace preserved."* — the promise the review's §4 already
measured as having **no such caller**. Under an honest `R∘G` reading with
`R = 0`, the projector is not vacuum's `B`; it is `B` *plus* an unrequested
identity on the unconstrained block. Cf. `assert_outgoing_leakage_unconstrained`
(§3): the BC is DEFINED not to touch the outflow rows, so a `B` that writes them
(even as identity) is over-specified.

---

## 3. The `assert_*` invariant suite

### 3.1 The five as shipped (`_base.py`)

**[M]** AST measurement of the base class bodies (docstring excluded):

| # | invariant | `_base.py` line | base-body statements | what it checks |
|---|---|---|---|---|
| 1 | `assert_inflow_outflow_classification(quadrature)` | `:216` | **0 — NO-OP** | intended: every ordinate is strictly inflow or outflow (no tangential) |
| 2 | `assert_outgoing_leakage_unconstrained(quadrature)` | `:227` | **0 — NO-OP** | intended: the BC does not constrain the outgoing trace |
| 3 | `assert_geometry_map_measure_preserving(quadrature)` | `:236` | **0 — NO-OP** | intended: `G` preserves `w(Ω)·|Ω·n|` |
| 4 | `assert_response_positive_if_declared()` | `:246` | **0 — NO-OP** | intended: a declared `R` yields non-negative output |
| 5 | `assert_source_lives_on_incoming_trace(quadrature, inflow_indices=None)` | `:253` | **3 — REAL** | the only one with a base implementation: evaluates `self.source` on shape `(N,)`; if any entry is nonzero AND `inflow_indices is None`, raises `BoundarySourceNotOnIncomingTraceError` (ERR-047) |

**[M]** The override matrix (measured by identity-comparing each law's bound
method against the base's):

| law | 1 classification | 2 leakage | 3 measure | 4 response⁺ | 5 source | law-specific extras |
|---|:---:|:---:|:---:|:---:|:---:|---|
| `VacuumInflow` | – | – | – | – | – | — |
| `ReflectiveBoundary` | – | – | **OVERRIDE** | – | – | `assert_is_involutive`, `assert_reflection_maps_inflow_to_outflow`, `assert_realizable` |
| `WhiteBoundary` | – | – | – | **OVERRIDE** | – | `assert_submarkov`, `assert_realizable` |
| `AlbedoBoundary` | – | – | – | **OVERRIDE** | – | `assert_submarkov`, `assert_realizable` |
| `PeriodicBoundary` | – | – | – | – | – | — |
| `PrescribedInflow` | – | – | – | – | – | — |
| `ZeroFluxBoundary` | – | – | – | – | – | — |

**Three findings from the matrix:**

1. **Invariants 1 and 2 are overridden by NOBODY.** `assert_inflow_outflow_classification`
   and `assert_outgoing_leakage_unconstrained` are permanent no-ops across the
   whole family — they have never once been given a body. `[G]`
   `grep -rn "def assert_inflow_outflow_classification\|def assert_outgoing_leakage_unconstrained" orpheus/`
   returns only the two `_base.py` definitions. Their typed error
   (`IncomingOutgoingTraceClassificationError`, `_errors.py:68`) is defined and
   **never raised anywhere in production**. `[G]`
2. **Four of seven laws override nothing at all.** `VacuumInflow`,
   `PeriodicBoundary`, `PrescribedInflow`, `ZeroFluxBoundary` inherit five
   no-ops and one real check. For them `assert_realizable` fires exactly one
   test (the ERR-047 source probe), and for the three with `NoSource` that
   test short-circuits on `if not np.any(probe): return`.
3. **The grand report specifies SIX universal invariants, not five.** **[R]**
   `.claude/plans/neutron_transport_grand_report_v3.md:3263-3272`:

   ```python
   boundary.assert_inflow_outflow_classification()
   boundary.assert_vacuum_sets_only_inflow_trace()      # <-- MISSING from the code
   boundary.assert_outgoing_leakage_unconstrained()
   boundary.assert_geometry_map_measure_preserving()
   boundary.assert_response_positive_if_declared()
   boundary.assert_source_lives_on_incoming_trace()
   ```

   `[G]` `assert_vacuum_sets_only_inflow_trace` — **0 hits in `orpheus/` and
   `tests/`.** The check itself is not absent in spirit: the SN realizer's
   vacuum arm hand-implements the equivalent (ERR-041 orientation guard,
   `sn/boundary/realizer.py:227-247`, raising `VacuumAppliedToOutgoingTraceError`).
   But it lives in the SN realizer, not on the law, so it is **method-specific
   code doing a law's job** — and the diffusion realizer therefore does not get
   it. The "five universal invariants" phrasing in `_base.py:11-12` is thus a
   *silent narrowing* of the source it cites.

### 3.2 The invariants are `-O` safe

**[G]** `grep -rn -E "^\s+assert\s" orpheus/geometry/boundary/` → **zero hits**.
Every `assert_*` method raises a typed `BoundaryError` subclass; the `assert_`
prefix is naming, not the Python statement. The canonical `python -O -m pytest`
invocation does **not** strip them. (The Mode-8 `-O` exposure the review's §8.4
flags is a *tests/* concern, not a law-layer one.)

### 3.3 `assert_realizable` — the seam call, and who actually calls it

**[R]** `_base.py:354-360` — the base template fires exactly the five:

```python
self.assert_inflow_outflow_classification(quadrature)
self.assert_outgoing_leakage_unconstrained(quadrature)
self.assert_geometry_map_measure_preserving(quadrature)
self.assert_response_positive_if_declared()
self.assert_source_lives_on_incoming_trace(quadrature, inflow_indices)
```

Three laws EXTEND it via `super()`: `ReflectiveBoundary` adds involution +
inflow→outflow (`reflective.py:276-278`); `WhiteBoundary` and `AlbedoBoundary`
each add `assert_submarkov` (`white.py:135-136`, `albedo.py:111-112`).

**Per-law, what `assert_realizable` ACTUALLY checks today:**

| law | effective checks fired |
|---|---|
| `VacuumInflow` | ERR-047 source probe only — and `NoSource` ⇒ it returns immediately. **Effectively nothing.** |
| `ReflectiveBoundary` | measure-preservation of the reflection table (ERR-042) + involution (ERR-044) + inflow→outflow image (ERR-045). **Three real checks — the only law that is genuinely certified.** |
| `WhiteBoundary` | `α ≥ 0` (ERR-043) + `α ≤ 1` (ERR-046). Two real scalar-range checks. `G_diff` is **not** checked for measure-preservation or for conservation. |
| `AlbedoBoundary` | `α ≥ 0` + `α ≤ 1`. Same two. |
| `PeriodicBoundary` | **nothing.** (Grand report `:3283-3289` specifies three periodic tests — `assert_pairing_bijective`, `assert_normals_opposite`, `assert_measure_preserving`; **[G]** all three are **0 hits** in `orpheus/` + `tests/`.) |
| `PrescribedInflow` | ERR-047 source probe — the **one law for which it bites**: a `ConstantInflowSource` with `inflow_indices=None` raises. |
| `ZeroFluxBoundary` | nothing (and by design it must NOT be response-positive-checked — `zero_flux.py:49-58`). |

**[G] Call sites of `assert_realizable` in production: exactly ONE.**
`orpheus/sn/boundary/realizer.py:194-197`:

```python
if isinstance(law, BoundaryTraceLaw):
    law.assert_realizable(
        quad, inflow_indices=method_space.inflow_indices
    )
```

(Line ~195, not ~292 as the brief estimated.) **[R]** `DiffusionBoundaryRealizer.realize`
(`orpheus/diffusion/boundary_realizer.py:176-195`) is two statements —
`albedo = self._partial_current_albedo(law); return stamp_boundary_role(_albedo_operator(albedo))`
— and calls **no invariant at all**.

> **This contradicts the docstring's own framing.** **[R]** `_base.py:328-333`:
> *"This template method is the production seam that honors those lessons: **a
> method realizer** (e.g. `SNBoundaryRealizer.realize`) calls it ONCE at entry,
> so **every law arrives at its primitive construction already certified**."*
> Half the production realizers do not. A `WhiteBoundary(albedo=5.0)` realized
> through diffusion silently produces `𝒜 = 5.0` — sub-Markov violated, no
> error. The invariant suite is SN-only in practice.

**[G] `assert_realizable` in tests:** one file only —
`tests/geometry/test_bc_universal_invariants.py`.

### 3.4 The grand report's reflection tests vs. what shipped

**[R]** `grand_report:3274-3281` specifies four:
`assert_is_involutive`, `assert_preserves_boundary_measure`,
`assert_maps_inflow_to_outflow`, `assert_direction_norm_preserved`.

The code ships three (two renamed): `assert_is_involutive` ✓;
`assert_geometry_map_measure_preserving` ✓ (renamed, and promoted from a
reflection-specific test to a universal-slot override);
`assert_reflection_maps_inflow_to_outflow` ✓ (renamed).
**[G]** `assert_direction_norm_preserved` — **0 hits**. `‖Ω_π(n)‖ = ‖Ω_n‖` is
never checked; a reflection table that pairs an ordinate with one of different
norm would pass all three shipped checks if the cosine measures happened to
match.

---

## 4. Registry / admission — and the reachability question

### 4.1 How `BoundaryTraceLaw.registry` self-populates

**[R]** `_base.py:142-146` — the root re-declares its own dict and points
`_registry_base` at itself:

```python
registry: ClassVar[dict[str, type["BoundaryTraceLaw"]]] = {}

@classmethod
def _registry_base(cls) -> type:
    return BoundaryTraceLaw
```

**[R]** `orpheus/numerics/registry.py:87-130` — `RegistryMixin.__init_subclass__`
consumes the `key=` class-creation kwarg, sets `cls.key = key`, and inserts into
`base.registry[key]`, raising `KeyError` on a duplicate. Registration is
therefore *at class-creation time*, i.e. **at import of the defining module**.
Two mechanical properties follow, both real: a new law without a `key=` is
syntactically obvious (the kwarg sits next to the base class), and a key
collision is an import-time crash rather than a silent shadow.

One subtlety worth recording: `__init_subclass__` carries a
`@dataclass(slots=True)` re-entry branch (`registry.py:95-119`) that re-points
the registry entry at the replacement class. **[G]** No boundary law uses
`slots=True` today (`grep -rn "slots=True" orpheus/geometry/boundary/` → 0 hits),
so that branch is dormant here — though the grand report's skeleton
(`grand_report:2791`) does specify `@dataclass(frozen=True, slots=True)`.

**[M]** All seven laws are registered by importing `orpheus.geometry.boundary`
alone; importing the SN and diffusion meshes adds nothing:

```
BoundaryTraceLaw.registry = {albedo, periodic, prescribed_inflow,
                             reflective, vacuum, white, zero_flux}   (7)
```

There are therefore **two registries in play** and they are NOT the same thing:

| registry | scope | contents |
|---|---|---|
| `BoundaryTraceLaw.registry` | the **physical vocabulary** — every law that exists | 7 keys |
| `<Mesh>.BOUNDARY_OPERATOR_REGISTRY` | the **per-method admission table** | SN: 2 · diffusion: 4 |

This is exactly the separation the grand report prescribes (**[R]**
`grand_report:3250-3257`): *"Keep method-specific realizations in a separate
registry keyed by method: `BoundaryLawRegistry -> physical boundary vocabulary`
/ `BoundaryRealizerRegistry -> method-specific enforcement`. This separation
prevents `VacuumBoundary` from becoming an overloaded method-specific object."*
The implementation replaced the second registry with a per-mesh ClassVar
(#290 P7b dissolved the realizer registry), which preserves the separation.

### 4.2 The admission tables — VERIFIED, and the prior claim is CORRECT

**[R]** `orpheus/sn/mesh/augmented_mesh.py:171-174`:

```python
BOUNDARY_OPERATOR_REGISTRY: ClassVar[dict[str, type[BoundaryTraceLaw]]] = {
    "vacuum": VacuumInflow,
    "reflective": ReflectiveBoundary,
}
```

with the comment (**[R]** `:182-187`) stating the intent explicitly: *"The 5
other kinds the realizer handles today (`white`, `periodic`, `albedo`,
`prescribed_inflow`, `mixed`) are NOT registered here — adding them requires
SN-sweep-side wiring (sweep cycles for periodic, etc.)."*

**[R]** `orpheus/diffusion/augmented_mesh.py:158-163`:

```python
BOUNDARY_OPERATOR_REGISTRY: "ClassVar[dict[str, type[BoundaryTraceLaw]]]" = {
    "vacuum": VacuumInflow,
    "reflective": ReflectiveBoundary,
    "albedo": AlbedoBoundary,
    "zero_flux": ZeroFluxBoundary,
}
```

**The prior review's §6 claim is VERIFIED** — SN admits only `vacuum` and
`reflective`. **[M]** End-to-end, constructing a real mesh per tag:

| `BC(kind)` | SNMesh | DiffusionMesh |
|---|---|---|
| `vacuum` | **ADMITTED** | **ADMITTED** |
| `reflective` | **ADMITTED** | **ADMITTED** |
| `white` | `ValueError: SNMesh does not support boundary condition 'white' … Supported: 'reflective', 'vacuum'.` | `ValueError: … Supported: 'albedo', 'reflective', 'vacuum', 'zero_flux'.` |
| `periodic` | `ValueError` (same) | `ValueError` (same) |
| `albedo` | `ValueError` | **ADMITTED** |
| `prescribed_inflow` | `ValueError` | `ValueError` |
| `zero_flux` | `ValueError` | **ADMITTED** |

> **Stronger than the prior claim.** Not just white and periodic: **three laws
> — `WhiteBoundary`, `PeriodicBoundary`, `PrescribedInflow` — are unreachable
> from a `BC` declaration under EVERY method in the codebase.** 3 of 7 laws
> (43%) have no declarative entry point at all.

### 4.3 Is white/periodic reachable by ANY other production path? — NO

Four candidate paths, each audited:

1. **`BC` declaration → `resolve_boundary_conditions` → `_law_from_tag`.**
   **[M]** Refuses, per the table above. **[M]** AST: `realize_boundary_law` has
   exactly **one** call site in `orpheus/` — `transport/method.py:259`, inside
   `resolve_boundary_conditions`, whose `law` comes from `_law_from_tag`, which
   goes through the admission table. No bypass.
2. **Direct law construction in production code.** **[G]**
   `grep -rn -E "\b(WhiteBoundary|PeriodicBoundary|PrescribedInflow|ZeroFluxBoundary)\(" orpheus/`
   — every hit is a **docstring**. The ONE real non-{vacuum,reflective}
   construction anywhere in `orpheus/` is `AlbedoBoundary(albedo=…)` at
   `transport/method.py:305`, inside `_law_from_tag` itself — i.e. gated by the
   admission table.
3. **`realize_recursively` (the rank-N walker).** **[M]** Zero production
   callers — see §5.2.
4. **`BoundaryTraceLaw.create("white")`.** **[G]** No production caller; the
   only mention is a comment in `zero_flux.py:94`.

**[M]** AST: the realizer seam `<X>.realize(...)` has exactly **three** call
sites in `orpheus/` — `sn/mesh/augmented_mesh.py:426`,
`diffusion/augmented_mesh.py:315`, and the walker's own
`_realizer.py:303`. Both mesh sites are inside `realize_boundary_law`.

**Conclusion [M/G]: `WhiteBoundary`, `PeriodicBoundary` and `PrescribedInflow`
are reachable ONLY from test code.** Their realizer arms
(`sn/boundary/realizer.py:290-305`, `:323-332`, `:334-349`) are dead in
production. `ZeroFluxBoundary` is reachable — but only under diffusion.

### 4.4 Two admission-layer defects found while measuring

**[R] `BC.white` is a public convenience instance that no method admits.**
`orpheus/geometry/mesh.py:124-126`:

```python
BC.vacuum = BC("vacuum")
BC.reflective = BC("reflective")
BC.white = BC("white")          # <-- tab-completable, refused everywhere
```

**[M]** `BC.white` → `ValueError` on both `SNMesh` and `DiffusionMesh`.
The `BC` docstring even advertises it: **[R]** `mesh.py:56-58` *"Boundary
condition identifier (e.g. `"vacuum"`, `"reflective"`, `"white"`)."*

**[G] The `"partial"` tag is unregistered.** `ReflectiveBoundary.kind` returns
`"partial"` when `albedo ≠ 1` (`reflective.py:107`), `BC.to_alpha` handles
`BC("partial", {"albedo": x})` (`mesh.py:109-116`), and a derivation module
tells users to write it (`orpheus/derivations/continuous/singular_eigenfunction/spectrum.py:192`).
But `"partial"` is **not** a `BoundaryTraceLaw.registry` key and not in either
admission table, so `BC("partial", {"albedo": 0.7})` cannot construct a mesh.
The path that WOULD work is `BC("albedo", {"albedo": 0.7})` — and that one is
diffusion-only. Note also that `realize_boundary_law` tags the shim with
`law.key` (`augmented_mesh.py:427`), **not** `law.kind`, so a partial-albedo
reflective law would produce a shim tagged `"reflective"` while the law itself
reports `kind == "partial"` — two disagreeing tag sources on one object.

---

## 5. `LawSum` / `LawScaled` — the descriptor-tree algebra

### 5.1 What is closed, and what the algebra guarantees

**[R]** `_composition.py:65-71` — the algebra is over
`LawNode = Union[BoundaryTraceLaw, LawSum, LawScaled]`. Both composers are
`@dataclass(frozen=True)`. Every one of `__add__`, `__radd__`, `__sub__`,
`__rsub__`, `__mul__`, `__rmul__`, `__truediv__`, `__neg__` is implemented on
**all three** node kinds (`_base.py:177-209`, `_composition.py:103-128`,
`:150-175`) — so the algebra is genuinely closed: no operation on `LawNode`s can
escape the `LawNode` family.

Guarantees, verified by reading:

| guarantee | evidence |
|---|---|
| **Never returns an operator.** Descriptor-space is sealed from operator-space | `_base.py:174` comment; every dunder returns `LawSum`/`LawScaled` |
| **Constant folding.** `LawScaled(α, LawScaled(β, x))` never exists at rest | `_composition.py:118` `return LawScaled(float(scalar) * self.scalar, self.inner)` |
| **Not flattened.** `(a+b)+c` is `LawSum(LawSum(a,b), c)` — deliberately | `_composition.py:31-37`, which also warns tests not to assert tree shape |
| **Subtraction is sugar.** `a − b ≡ LawSum(a, LawScaled(−1, b))` | `_composition.py:38`, `_base.py:185-187` |
| **Scalar guard.** non-numeric `*` returns `NotImplemented` | `_base.py:195-196`, `_composition.py:116-117`, `:163-164` |
| **No mixing with operators.** `LawNode + LinearOperator` is unsupported by contract | `_composition.py:41-45` |

The structure-preserving walk is `LawSum → OperatorSum`, `LawScaled →
ScaledOperator`, leaf → `realizer.realize(leaf, space)` (**[R]**
`_realizer.py:301-314`), with a `TypeError` naming the offending type for
anything else (`:316-323`).

**One gap in the closure:** the algebra has **no product**. `__matmul__` /
`__mul__(LawNode, LawNode)` do not exist, so `R ∘ G` — the very composition the
affine form is built from — **is not expressible in the law algebra**. The
algebra realizes `Σ_α c_α G_α` (§15.2) and nothing else. That is consistent with
the review's `[G]` finding that `OperatorProduct` is never used by any boundary
operator: the descriptor layer cannot even *spell* the composition its own
docstring says it implements.

### 5.2 Is rank-N composition exercised in PRODUCTION? — **NO. Tests only.**

**[M]** AST walk of every `*.py` under `orpheus/`, counting real `Call` nodes
(docstrings excluded by construction):

```
AST call sites of realize_recursively in orpheus/ :  3
    orpheus/geometry/boundary/_realizer.py:307   <- its own recursion (LawScaled arm)
    orpheus/geometry/boundary/_realizer.py:312   <- its own recursion (LawSum arm, a)
    orpheus/geometry/boundary/_realizer.py:313   <- its own recursion (LawSum arm, b)
AST call sites in tests/ : 11
    tests/diffusion/test_boundary_realizer.py:325, 334, 347
    tests/geometry/test_law_composition.py:214, 228, 248, 272, 296, 322, 352, 362
```

**Every production "call site" the text-grep reports is a docstring.** The one
that looks most like production — `orpheus/diffusion/boundary_realizer.py:93-94`
— sits inside the **module docstring** (which runs `:1-114`); the module's first
executable line is `from __future__ import annotations` at `:116`.

**[G]** `LawSum` / `LawScaled` are imported in `orpheus/` **only** by
`geometry/boundary/_realizer.py:79` (to `isinstance`-dispatch on them) and
`_base.py:65` (TYPE_CHECKING). No solver, mesh, or realizer outside the boundary
package references them. The only other production mention is a *comment*
(`sn/boundary/realizer.py:191`).

> **Plainly: the entire rank-N composition path — `LawSum`, `LawScaled`,
> `LawNode`, `realize_recursively`, and the eight algebra dunders on
> `BoundaryTraceLaw` — has ZERO production consumers.** It is exercised
> exclusively by `tests/geometry/test_law_composition.py` (+3 diffusion cases).
> `realize_recursively`'s own comment says it (`_realizer.py:200-204`:
> *"Production realizes single BCs directly … this walker is the rank-N
> composition entry point (the Marshak `0.3·Reflective + 0.7·White` partial-current
> BC), exercised by the law-composition wall"*) — and that self-description is
> accurate. What the comment does not say is that the Marshak mix it names is
> **unreachable anyway**: `WhiteBoundary` is admitted by no method (§4.2), so
> even wiring the walker into production would not make that example work.
>
> This is a **fourth instance of the dead-capability pattern** the review's §4
> identifies — and the largest one by code volume.

---

## 6. The affine form's provenance — VERIFIED against the repo

### 6.1 The Grand Report IS in the repo

**[G]** `.claude/plans/neutron_transport_grand_report_v3.md` — 6600+ lines,
tracked. Not an external document. Both of the prior review's `[U]` claims are
now raisable to `[R]`, and **both are CORRECT, verbatim**:

**Claim A — §16A.2 spells it `AffineBoundaryOperator(linear=R @ G, source=q)`.**
**VERIFIED [R]** `neutron_transport_grand_report_v3.md:2791-2801`:

```python
@dataclass(frozen=True, slots=True)
class BoundaryTraceLaw:
    geometry_map: BoundaryGeometryMap
    response: BoundaryResponseKernel
    source: BoundarySource

    def as_operator(self, trace_space):
        G = self.geometry_map.as_operator(trace_space)
        R = self.response.as_operator(trace_space)
        q = self.source.as_field(trace_space.incoming)
        return AffineBoundaryOperator(linear=R @ G, source=q)
```

**Claim B — §16A.5's realizer sketch dispatches on `isinstance(law.response, ZeroBoundaryResponse)`.**
**VERIFIED [R]** `neutron_transport_grand_report_v3.md:2946-2960`:

```python
@dataclass(frozen=True, slots=True)
class DiscreteOrdinatesBoundaryRealizer:
    def realize(self, law: BoundaryTraceLaw, space: DiscreteOrdinatesPhaseSpace):
        trace = DiscreteOrdinatesBoundaryTrace(space)

        if isinstance(law.response, ZeroBoundaryResponse):
            return ZeroInflowOrdinateConstraint(
                incoming_mask=trace.incoming_mask()
            )

        G = trace.discretize_geometry_map(law.geometry_map)
        R = trace.discretize_response(law.response)
        q = trace.discretize_source(law.source)
        return SNBoundaryOperator(linear=R @ G, source=q)
```

**[G]** `AffineBoundaryOperator`, `ZeroBoundaryResponse`,
`BoundaryResponseKernel`, `BoundaryGeometryMap`, `ZeroInflowTrace`,
`BoundaryRelation` — **all 0 hits in `orpheus/`**. The design's entire
factor-object vocabulary is unimplemented.

### 6.2 TWO factors, not three — and the design intends `R` AFTER `G`

The composition is unambiguously **two linear factors plus an affine offset**:
`R ∘ G` then `+ q`. Three independent statements agree:

1. **[R]** `grand_report:2749`: `BoundaryTraceLaw   # affine law gamma_- = R G gamma_+ + q`
2. **[R]** `grand_report:2801` / `:2960`: `linear = R @ G` — a single `@`.
3. **[R]** `docs/theory/foundations/boundary_conditions.rst:191-220` — the
   theory page's own equation `:label: affine-bc-form`, which is SHARPER than
   the grand report because it types the factors:

   > `G : Γ₊ → Γ₊` is the **geometric map** — a measure-preserving permutation,
   > pushforward, spatial wrap-around, or hemispheric cosine-weighted average.
   > It carries pure geometry (it changes nothing about the physical interaction
   > at the boundary; **it just relabels the angular fluxes that meet there**).
   >
   > `R : Γ₊ → Γ₋` is the **response kernel** — a scalar amplitude in `[0,1]`
   > for the standard sub-Markov BCs (albedo, white, partial-current) or a full
   > angular kernel in general weak-form BCs (deferred …).

   So `G` is an **endomorphism of the outflow trace** (relabeling) and `R` is
   the map that **crosses to the inflow trace**. That typing is the crispest
   statement of the factorization anywhere in the repo, and it is what makes
   §2.1's table decidable: `G` is "which outflow DOF", `R` is "how much".

**[R]** The same page also carries the master-form docstring `_base.py:6-8` and
`_source.py:5-9`, and the residual definition `boundary_conditions.rst:1280`
(`r_Γ = γ₋ψ − (R G γ₊ψ + q)`), all consistent.

### 6.3 Does any design text say the realizer should dispatch on the RESPONSE?

**Yes — one, and it is the grand report itself** (§16A.5, quoted above). Nothing
else in the repo mandates it. The theory page hedges, and is careful to describe
the class dispatch as the CURRENT state, not the target — **[R]**
`boundary_conditions.rst:250-256`:

> The split lets cross-method realizers introspect the law's geometric and
> response components separately — **the SN realizer dispatches on the law's
> class today, but a future weak-form realizer might dispatch on the geometry /
> response / source independently.**

So: response-dispatch is a **grand-report design directive** and a **theory-page
"might"**, not a landed contract anyone has broken. The load-bearing sentence
for the empty-property decision is the one immediately before it, which IS
stated as fact and is FALSE:

> **[R]** `boundary_conditions.rst:245-252` — *"**The three operators are
> first-class.** The `BoundaryTraceLaw` ABC exposes `geometry_map`,
> `response_kernel`, and `source` as Python properties on every concrete
> subclass. The properties default to `None / None / NoSource()`; **concrete
> laws override when applicable.**"*

**[M]** No concrete law overrides `geometry_map` or `response_kernel` (§1.3).
Only `PrescribedInflow` overrides `source`. "Override when applicable" is true
of one property out of three, on one law out of seven.

---

## 7. Promises not kept — the law layer's own list

Ten, all verified. (The review's §4 already holds three subsystem-wide; these
are additional and specific to `orpheus/geometry/boundary/` + its theory page.)

### P1. "The Wave-7 deprecated aliases … are re-exported below"

**[R]** `orpheus/geometry/boundary/__init__.py:354-357`:

> *"The Wave-7 deprecated aliases (``VacuumBoundaryOperator`` = `VacuumInflow`,
> etc.) are **re-exported below** for backward compatibility with consumers that
> import the old names. They will be removed in a future cleanup wave once every
> consumer migrates."*

**[M]** None is importable:
```
names NOT importable from orpheus.geometry.boundary:
  ['VacuumBoundaryOperator', 'SpecularBoundaryOperator', 'WhiteBoundaryOperator',
   'PeriodicBoundaryOperator', 'AlbedoBoundaryOperator', 'MixedBoundaryOperator']
```
They were retired in Wave O step O.4a.1 — every law module's own docstring says
so (`vacuum.py:48-51`, `reflective.py:80-84`, `white.py:64-68`,
`albedo.py:49-53`, `periodic.py:45-50`). The package `__init__` was not updated
with them, and it **still advertises five of them inline** at `:91`, `:102`,
`:110`, `:117`, `:124` as "deprecated alias `X`".

### P2. "each carrying the three first-class properties"

**[R]** `__init__.py:42-46` — *"`BoundaryTraceLaw` ABC + **six** concrete
subclasses, **each carrying** the three first-class properties (`geometry_map`,
`response_kernel`, `source`)."* **[M]** Seven subclasses, and none carries the
first two. The same docstring says "seven" 40 lines later (`:82`) — the count is
internally inconsistent too.

### P3. "concrete laws override when applicable" (theory page)

`boundary_conditions.rst:250-252` — see §6.3. Same defect, in the doc of record.

### P4. `_base.py:89-91` — "Concrete subclasses … populate these"

**[R]** *"Concrete subclasses (`VacuumInflow`, `ReflectiveBoundary`,
`WhiteBoundary`, `AlbedoBoundary`, `PeriodicBoundary`, `PrescribedInflow`)
**populate these**."* **[M]** None does. (This is the review's §4 row, confirmed
independently at its source line.) The list also omits `ZeroFluxBoundary`.

### P5. `BoundaryTraceLaw` "combines `LinearOperator` … and `RegistryMixin`"

**[R]** `boundary_conditions.rst:1661-1666`:

> *"The base class `BoundaryTraceLaw` is an `abc.ABC` that combines
> `LinearOperator` (for operator-algebra dunders like `+`, `*`, `@`) and
> `RegistryMixin`…"*

**[M]** `BoundaryTraceLaw.__mro__ == [BoundaryTraceLaw, RegistryMixin, ABC, object]`.
**The same page contradicts itself at `:3157-3158`** (*"`BoundaryTraceLaw` no
longer inherits `LinearOperator`"*) — the §1661 block is pre-#186 text that
survived the descriptor-model rewrite. It also promises an `@` dunder, which the
law algebra does not have at all (§5.1).

### P6. The universal-invariants table documents behaviour that isn't there

**[R]** `boundary_conditions.rst:2500-2557`, four defects:

| doc claim (`file:line`) | reality |
|---|---|
| `:2512` *"reflective overrides to require strict partition"* (`assert_inflow_outflow_classification`) | **[M] FALSE** — `ReflectiveBoundary` does not override it. Nobody does. |
| `:2523` *"Reflective overrides (**delegates to the involution check**)"* (`assert_geometry_map_measure_preserving`) | Override is real; the delegation is **[R] FALSE** — `reflective.py:166-172` states the opposite: *"This is checked **DIRECTLY, independent of the involution property**… exactly the hole the pre-#52 delegation ('weight equality is implied by construction') left open."* The doc still describes the pre-#52 code. |
| `:2535-2536` *"Prescribed-inflow + future user-source classes override"* (`assert_source_lives_on_incoming_trace`) | **[M] FALSE** — `PrescribedInflow` does not override it; the base body handles every law. |
| `:2555-2557` documents `ReflectiveBoundary.assert_direction_norm_preserved` | **[G] 0 hits in `orpheus/` + `tests/`.** A documented method that does not exist. `:2550` also names `assert_maps_inflow_to_outflow`, whose real name is `assert_reflection_maps_inflow_to_outflow`. Both are Python-domain `:meth:` roles — per `coding-standards.md` they render as plain text with **no `-W` warning**, so the Sphinx gate cannot catch either. |

### P7. "Each `BoundaryTraceLaw` subclass overrides the five universal `assert_*` invariants"

**[R]** `__init__.py:255-257`. **[M]** Four of seven laws (`VacuumInflow`,
`PeriodicBoundary`, `PrescribedInflow`, `ZeroFluxBoundary`) override **none**;
two of the five universal slots are overridden by **nobody**.

### P8. `assert_realizable` is "the production seam … every law arrives already certified"

**[R]** `_base.py:328-333`. **[R]** `DiffusionBoundaryRealizer.realize`
(`orpheus/diffusion/boundary_realizer.py:176-195`) never calls it — see §3.3.
The guarantee holds for SN only. A `WhiteBoundary(albedo=5.0)` or an
`AlbedoBoundary(albedo=-2.0)` realized through diffusion produces `𝒜 = 5.0` /
`𝒜 = -2.0` with no sub-Markov or positivity error.

### P9. `IncomingSourceOperator` "returns `source.evaluate(psi_out.shape)`"

**[R]** `__init__.py:132-137` describes the pre-#52 unmasked behaviour.
**[R]** The operator now masks to the inflow ordinates when the realizer supplies
them (`sn/boundary/realizer.py:343-348`; `sn/boundary/angular.py` docstring:
*"with a mask, the evaluation is zeroed on every non-inflow ordinate (outflow
AND tangential …)"*). Minor, but it is the package `__init__` — the first thing
a fresh reader loads.

### P10. Two doc listings show the PRE-P7b 2-arg `realize_recursively`

**[R]** `docs/theory/foundations/boundary_conditions.rst:2831-2854` introduces
the listing as current — *"The dispatch is exhaustive on the descriptor-tree
node types"* — and then shows:

```python
def realize_recursively(
    node: BoundaryTraceLaw | LawSum | LawScaled,
    method_space: SNMethodSpace,                 # <-- geometry annotated with the SN space
) -> LinearOperator:
    if isinstance(node, BoundaryTraceLaw):
        return SNBoundaryRealizer().realize(node, method_space)   # <-- hardcoded SN realizer
```

Three defects in one block: the `realizer` parameter is missing (it is
**REQUIRED** today — `_realizer.py:212-216`), the SN realizer is hardcoded
inside the method-BLIND walker, and `method_space` is annotated `SNMethodSpace`
— i.e. the listing depicts exactly the coupling #290 P7b removed. The SAME page
uses the correct 3-arg call 20 lines later (`:2874`).

**[R]** `docs/theory/methods/sn/boundary_conditions.rst:240` repeats it in a
user-facing snippet: `op = realize_recursively(tree, method_space)`. Copy-pasted,
that raises `TypeError`.

### Notational defect (not a broken promise, but a correctness hazard)

**[R]** `__init__.py:103`, `:111`, `:117`, `:124` write the reflective law as
`R = G_refl · α`, white as `R = G_diff · α`, albedo as `R = I · α`, periodic as
*"`R` is a spatial pushforward"*. Here `R` denotes the **composite** `R·G`.
`_base.py:8` and `boundary_conditions.rst:210-220` use `R` for the **response
factor alone**. The two spellings are in the same package and collide exactly on
the symbol whose factorization this review is deciding. Any plan that populates
`response_kernel` must first pick one meaning of `R`.

---

## 8. Two additional measurements worth carrying

### 8.1 ERR-040 is a typed error that production never raises

**[G]** `IncomingOutgoingTraceClassificationError` appears in `orpheus/` only in
`_errors.py` (its definition), `_base.py:223` (a docstring `Raises` mention on a
no-op), and `__init__.py:260/456/509` (prose + re-export). **The only
`pytest.raises` on it is a synthetic hand-construction**
(`tests/geometry/test_bc_errors.py:66,73` — the test builds the error object and
raises it itself). It is a catalog entry with no production raiser, because its
invariant (§3.1 #1) has no body on any law.

### 8.2 The test stub registers into the PRODUCTION registry

**[M]** Importing `tests/geometry/test_boundary_trace_law.py` mutates the
process-global law vocabulary:

```
registry before: [albedo, periodic, prescribed_inflow, reflective, vacuum, white, zero_flux]
registry after : [_stub_for_test, albedo, periodic, prescribed_inflow, reflective, vacuum, white, zero_flux]
```

`_StubLaw(BoundaryTraceLaw, key="_stub_for_test")` (`:53`) is registered at
class-creation time into `BoundaryTraceLaw.registry` with no teardown. Any test
that asserts on the registry's exact contents is order-dependent on whether that
module has been imported. (Not a production defect — the admission tables are
per-mesh ClassVars and never consult the law registry — but it is a real
registration-timing hazard of the same family the #290 P7b note cites as the
reason the realizer registry was dissolved.)

---

## 9. Synthesis — what this quadrant decides

### 9.1 On `geometry_map` / `response_kernel`: populate, but not symmetrically

The two properties are **not** the same case, and the review's §3.2
"same symptom, different disease" framing needs one more split:

* **`response_kernel` is already real, under a different name.** Every law's
  `R` is a scalar; both production realizers already read it as `law.albedo`
  (`sn/boundary/realizer.py:288`, `diffusion/boundary_realizer.py:214`).
  Populating the property is a **rename-and-lift**, not a new capability — and
  it makes `_law_from_tag`'s `if law_cls is AlbedoBoundary` special case
  (`transport/method.py:303-310`) and the diffusion realizer's five-arm
  `isinstance` ladder (`:209-216`) collapse to one read. **[R]** all cited.
* **`geometry_map` is populatable as a typed SPEC for 6 of 7 laws** (§2.1). The
  blocker on the 7th is `PeriodicBoundary`'s **missing partner-face field**, not
  a discretization wall. The `-> Any` return signature (`_base.py:149`) is the
  spec defect the review's §2.3 correction identified — confirmed at the line.
* **`ZeroFluxBoundary` should not be forced into the pair at all** (§2.3): it is
  a `BoundaryRelation`, not an `AffineTraceLaw`, and the grand report has a tier
  for that which the code never built.

### 9.2 On rank-N composition: this is a decision, not a detail

**[M]** `LawSum` + `LawScaled` + `LawNode` + `realize_recursively` + eight
dunders on the ABC = a complete, closed, well-tested algebra with **zero
production consumers**, whose flagship example (`0.3·Reflective + 0.7·White`) is
**unreachable even in principle** because no method admits `white` (§4.2). Under
`coding-standards.md`'s aggressive-retirement rule this is a retirement
candidate; under the fuller-view-oracle exception it is not (it is not an oracle
for anything). The honest options are (a) wire it — which requires admitting
`white` into `SNMesh.BOUNDARY_OPERATOR_REGISTRY` first, or (b) retire it. What
it must not remain is half-alive.

### 9.3 The `∘` gap is the same gap in two places

The law algebra has `+` and `*` but **no product** (§5.1); the operator side has
`OperatorProduct` but **no boundary operator uses it** (review §2.3). So `R ∘ G`
— the composition the affine form is *defined by* — is unspellable at the
descriptor layer and unused at the operator layer. Populating `geometry_map` /
`response_kernel` without also giving the descriptor layer a way to compose them
would reproduce the same half-measure one level down.

---

## 10. Gaps / explicitly NOT verified

* I did **not** verify the runtime behaviour of `AngularAverageOperator`,
  `IncomingOrdinateMaskTensor`, `PermutationOperator`, or `PeriodicWrapOperator`
  — those are the operator quadrant.
* I did **not** measure test teeth (the Mode-8 `-O` question) for
  `tests/geometry/test_bc_universal_invariants.py` — quadrant 4. I DID verify
  that the *production* invariants are `raise`-based and therefore `-O` safe
  (§3.2).
* I did **not** sweep `docs/theory/methods/sn/boundary_conditions.rst` (the
  second, SN-side page) for the full P1–P9 set; I only confirmed its stale
  `realize_recursively` call (P10). A doc sweep of that page is likely to find
  more.
* The grand report's `SymmetryBoundary`, `white_reflection`,
  `diffuse_reflection` semantic names (`grand_report:2767`, `:3246-3247`) have
  no code counterpart. **[G]** 0 hits. Not pursued — no evidence anyone intends
  them.
