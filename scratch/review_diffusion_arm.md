# The diffusion boundary arm — complete review

Companion quadrant to `.claude/plans/boundary_machinery_review.md` §8 item 1.
Branch `refactor/operator-strategy-layers` @ `73627b71`. Evidence key identical
to the parent report: **[M]** measured, **[R]** read (quoted + `file:line`),
**[G]** grep/exhaustive, **[U]** unverified, **[X]** corrected.

---

## 1. The realizer IS two named stages — CONFIRMED

**[R] The prior claim is CORRECT.** `DiffusionBoundaryRealizer.realize`
(`orpheus/diffusion/boundary_realizer.py:176-195`) is *not* an isinstance chain
producing an operator. Its entire body is two calls:

```python
    def realize(
        self,
        law: "BoundaryTraceLaw",
        method_space: DiffusionMethodSpace,
    ) -> LinearOperator:
        ...
        albedo = self._partial_current_albedo(law)
        return stamp_boundary_role(_albedo_operator(albedo))
```
— `orpheus/diffusion/boundary_realizer.py:176-195` (the body is `:194-195`).

### The two stages, named precisely

| stage | symbol | function | `file:line` | signature |
|---|---|---|---|---|
| **1. law → 𝒜** | `law → float` | `DiffusionBoundaryRealizer._partial_current_albedo` (a `@staticmethod`) | `orpheus/diffusion/boundary_realizer.py:197-249` | `(law: BoundaryTraceLaw) -> float` |
| **2. 𝒜 → operator** | `float → LinearOperator` | `_albedo_operator` (module-level, private) | `orpheus/diffusion/boundary_realizer.py:142-157` | `(albedo: float) -> LinearOperator` |
| (3. stamp) | role tag | `stamp_boundary_role` (shared, geometry layer) | `orpheus/geometry/boundary/_realizer.py:109-129` | mutates `op.block_role = BlockRole.BOUNDARY`, returns `op` |

**[R]** Stage 1 is where the isinstance chain lives — but it produces a **float**,
not an operator:

```python
        if isinstance(law, VacuumInflow):
            return 0.0  # Marshak: zero incoming current (ruling 3).
        if isinstance(law, ZeroFluxBoundary):
            return -1.0  # Dirichlet φ_Γ = 0 ⟺ J⁻ = −J⁺ under P1.
        if isinstance(law, (ReflectiveBoundary, WhiteBoundary)):
            return float(law.albedo)  # P1-coincident (see module docstring).
        if isinstance(law, AlbedoBoundary):
            return float(law.albedo)
```
— `orpheus/diffusion/boundary_realizer.py:209-216`.

**[R]** Stage 2 is a **structure-keyed collapse** — it dispatches on the *value*
of 𝒜, not on the law type, and never sees a law at all:

```python
def _albedo_operator(albedo: float) -> LinearOperator:
    ...
    if albedo == 0.0:
        return ZeroOperator()
    if albedo == 1.0:
        return IdentityOperator()
    return float(albedo) * IdentityOperator()
```
— `orpheus/diffusion/boundary_realizer.py:142-157`.

**[R]** The design is stated in the module docstring, so it is intentional, not
emergent: *":meth:`DiffusionBoundaryRealizer.realize` is the composition of two
named maps: *law → albedo scalar* 𝒜 (the physics table below), then *albedo
scalar → operator* (the structure-keyed collapse mirroring the SN realizer's
``AlbedoBoundary`` branch)."* — `orpheus/diffusion/boundary_realizer.py:26-29`.

### Why this shape matters as a template

The two stages have **different natures**, and that is the whole point:

- Stage 1 is **physics** — per-law, method-specific, an isinstance chain by
  necessity (it is the law→response table). Its codomain is a *number*, so it is
  trivially testable and trivially transportable: nothing about a `float` knows
  what a mesh, a quadrature, or a face shape is.
- Stage 2 is **structure** — zero physics, no `isinstance(law, ...)`, and
  **method-blind**: `_albedo_operator` reads only a scalar. It is the same
  `α ↦ α·I` collapse SN's `AlbedoBoundary` branch performs.

**[R]** Stage 2 reads NOTHING from `method_space` — the parameter is threaded
purely for Protocol conformance: *"The rank-1 albedo family reads nothing from
``method_space`` (mirroring SN's albedo branch — a
:meth:`DiffusionMethodSpace.minimal` suffices); the parameter is the uniform
:class:`BoundaryRealizer` Protocol surface the P4/P5 wiring threads per-face
spaces through."* — `orpheus/diffusion/boundary_realizer.py:187-192`.

**[M]** Verified by direct measurement — `DiffusionMethodSpace.minimal()` (all
three fields `None`: `{'mesh': None, 'face': None, 'trace': None}`) vs
`for_face(mesh, "xmin")` vs `for_face(mesh, "xmax")`:

| law | `minimal()` | `for_face(xmin)` | `for_face(xmax)` | same class | same value |
|---|---|---|---|---|---|
| `VacuumInflow` | `ZeroOperator` | `ZeroOperator` | `ZeroOperator` | ✓ | ✓ `[0. 0.]` |
| `ZeroFluxBoundary` | `ScaledOperator` | `ScaledOperator` | `ScaledOperator` | ✓ | ✓ `[-1.1 -1.2]` |
| `ReflectiveBoundary(1)` | `IdentityOperator` | `IdentityOperator` | `IdentityOperator` | ✓ | ✓ `[1.1 1.2]` |
| `AlbedoBoundary(0.37)` | `ScaledOperator` | `ScaledOperator` | `ScaledOperator` | ✓ | ✓ `[0.407 0.444]` |

The `method_space` parameter is **inert** for every law diffusion realizes.

> **Template reading.** In the parent report's `R ∘ G` vocabulary, stage 1
> produces **R** (the response amplitude) and stage 2 produces **the realized
> operator on the method's trace carrier**. The reason diffusion could split them
> cleanly is that its **G is degenerate** — see §6. The template is real but its
> cheapness is a consequence of the degeneracy, not independent of it.

---

## 2. Which laws diffusion admits — TWO gates, not one

There are **two independent admission surfaces**, and they disagree. This is a
finding in itself.

### Gate A — the registry (what a `BC` tag can name)

**[R]** `orpheus/diffusion/augmented_mesh.py:158-163`:

```python
    BOUNDARY_OPERATOR_REGISTRY: "ClassVar[dict[str, type[BoundaryTraceLaw]]]" = {
        "vacuum": VacuumInflow,
        "reflective": ReflectiveBoundary,
        "albedo": AlbedoBoundary,
        "zero_flux": ZeroFluxBoundary,
    }
```

**[M]** Measured end-to-end through `DiffusionMesh(Mesh1D(..., bc_left=tag), mats)`:

| `BC` tag | verdict | realized operator |
|---|---|---|
| `BC('vacuum')` | ADMITTED | `ZeroOperator` |
| `BC('reflective')` | ADMITTED | `IdentityOperator` |
| `BC('albedo', {'albedo': 0.4})` | ADMITTED | `ScaledOperator(scalar=0.4)` |
| `BC('zero_flux')` | ADMITTED | `ScaledOperator(scalar=-1.0)` |
| `BC('white')` | **REFUSED** | `ValueError: DiffusionMesh does not support boundary condition 'white' on face 'xmin'. Supported: 'albedo', 'reflective', 'vacuum', 'zero_flux'.` |
| `BC('periodic')` | **REFUSED** | same message, `'periodic'` |
| `BC('prescribed_inflow')` | **REFUSED** | same message, `'prescribed_inflow'` |
| `BC('mixed')` | **REFUSED** | same message, `'mixed'` |
| `BC('albedo')` (no param) | **REFUSED** | `ValueError: BC('albedo') on face 'xmin' requires an 'albedo' parameter; got params={}.` |

**[R]** The registry-miss message is generic, from the shared
`TransportMethod` body — `orpheus/transport/method.py:289-298`; the missing-param
message from `orpheus/transport/method.py:303-309`. Neither is diffusion-authored.

### Gate B — the realizer's own stage-1 chain (what `realize` accepts)

**[M]** Calling `DiffusionBoundaryRealizer().realize(law, DiffusionMethodSpace.minimal())`
directly, bypassing the registry:

| law | verdict |
|---|---|
| `VacuumInflow()` | realized → `ZeroOperator` |
| `ZeroFluxBoundary()` | realized → `ScaledOperator(-1.0)` |
| `ReflectiveBoundary(axis="x", albedo=α)` | realized → `IdentityOperator` (α=1) / `ScaledOperator(α)` |
| **`WhiteBoundary(axis="x", albedo=α)`** | **realized** → `IdentityOperator` (α=1) / `ScaledOperator(α)` |
| `AlbedoBoundary(albedo=α)` | realized → `Zero`/`Identity`/`Scaled` by α |
| `PeriodicBoundary()` | `BoundaryError` (quoted below) |
| `PrescribedInflow()` | `BoundaryError` (quoted below) |

**[M]** Refusal messages, verbatim from the run:

> `BoundaryError: DiffusionBoundaryRealizer cannot realize PeriodicBoundary: the
> periodic law couples a face to its OPPOSITE face (a trace-block permutation),
> not a per-face albedo J⁻ = 𝒜·J⁺. The wrap lands with the boundary-operator
> assembly when a diffusion consumer exists (#290 P4 seam).`

> `BoundaryError: DiffusionBoundaryRealizer cannot realize PrescribedInflow: J⁻ =
> q is the rank-0 AFFINE law — its realization is the boundary source q.boundary,
> not a linear boundary operator B (the same operator/source split SN keeps). The
> diffusion fixed-source arm lands with the solver wiring (#290 P5).`

**[R]** Plus a catch-all for any unrecognized law class,
`orpheus/diffusion/boundary_realizer.py:239-249`, which *enumerates the realizable
set in the message* — including `WhiteBoundary`.

### The disagreement — `white` is realizable but unreachable

**[G]** `WhiteBoundary` is handled by stage 1
(`boundary_realizer.py:213`, sharing the `ReflectiveBoundary` branch) and named in
the catch-all's "Realizable laws:" list (`:242-245`), yet `"white"` is **absent**
from the registry — so no `BC("white")` declaration can reach it. **[R]** This is
deliberate and documented at the registry:

> *"``"white"`` is deliberately absent: at the P1 level white coincides with
> reflective (the P3 realizer's coincidence note) — declare ``reflective`` or
> ``albedo``."* — `orpheus/diffusion/augmented_mesh.py:154-157`

**Assessment:** this is defensible physics (§6) but it means the diffusion arm has
**a live-but-unreachable branch**: `isinstance(law, (ReflectiveBoundary, WhiteBoundary))`
where the `WhiteBoundary` arm is production-dead unless someone calls the realizer
directly. It is *not* the parent report's §4 "dead capability" pattern — the code
is exercised (**[G]** `tests/diffusion/test_boundary_realizer.py:119, 131, 157,
176, 200, 264` construct `WhiteBoundary(...)` and realize it directly), and it
costs one tuple element —
but it is a second instance of the subsystem's habit of admitting more at the
realizer than the registry lets through. Note the **SN registry has the same
shape**: **[R]** `orpheus/sn/mesh/augmented_mesh.py:171-174` admits only
`vacuum`/`reflective`, with a comment (`:182-185`) that *"The 5 other kinds the
realizer handles today (``white``, ``periodic``, ``albedo``,
``prescribed_inflow``, ``mixed``) are NOT registered here."*
**[M]** confirms `set(SNMesh.BOUNDARY_OPERATOR_REGISTRY) == {"vacuum", "reflective"}`
is also pinned by a test (`tests/sn/operators/test_snmesh_realizer_wiring.py:403`).
→ **This raises the parent report's §6 [U] to [R]/[M]: SN admits exactly
`vacuum` + `reflective`; diffusion admits those two plus `albedo` + `zero_flux`.**

### Law universe: 7 laws

**[G]** `grep -rn "^class " orpheus/geometry/boundary/*.py` — the complete
`BoundaryTraceLaw` subclass set is 7: `VacuumInflow`, `ReflectiveBoundary`,
`WhiteBoundary`, `AlbedoBoundary`, `ZeroFluxBoundary`, `PeriodicBoundary`,
`PrescribedInflow` (plus the non-law composition carriers `LawSum` / `LawScaled`).
So: **diffusion registry 4/7, diffusion realizer 5/7, SN registry 2/7.**

---

## 3. The full per-law table — law → 𝒜 → operator → what `apply` computes

Every `apply` body below is quoted, not paraphrased. The realized operator's
carrier is a **bare `np.ndarray` of shape `(ng, *face_spatial)`** — for 1-D
diffusion, `(ng,)` (§5).

| law | 𝒜 (stage 1) | realized class | `apply(x)` in array terms |
|---|---|---|---|
| `VacuumInflow` | `0.0` | `ZeroOperator` | `0.0 * x` |
| `ZeroFluxBoundary` | `-1.0` | `ScaledOperator(-1.0, IdentityOperator())` | `-1.0 * x` |
| `ReflectiveBoundary(α=1)` | `1.0` | `IdentityOperator` | `x` |
| `ReflectiveBoundary(α≠1)` | `α` | `ScaledOperator(α, IdentityOperator())` | `α * x` |
| `WhiteBoundary(α=1)` | `1.0` | `IdentityOperator` | `x` |
| `WhiteBoundary(α≠1)` | `α` | `ScaledOperator(α, IdentityOperator())` | `α * x` |
| `AlbedoBoundary(α=0)` | `0.0` | `ZeroOperator` | `0.0 * x` |
| `AlbedoBoundary(α=1)` | `1.0` | `IdentityOperator` | `x` |
| `AlbedoBoundary(α∉{0,1})` | `α` | `ScaledOperator(α, IdentityOperator())` | `α * x` |
| `PeriodicBoundary` | — | REFUSED | — |
| `PrescribedInflow` | — | REFUSED | — |

**[R]** The three `apply` bodies, verbatim from `orpheus/numerics/operator.py`:

```python
class ZeroOperator(LinearOperator[Domain, Codomain]):
    def apply(self, x: Domain, /) -> Codomain:
        if self._codomain_zero is not None:
            return self._codomain_zero(x)
        ...
        return cast("Codomain", 0.0 * x)
```
— `:1920-1928`. The realizer passes **no** `codomain_zero`
(`boundary_realizer.py:146`: bare `ZeroOperator()`), so the endomorphic
`0.0 * x` branch is the one taken. On an `np.ndarray` input this is
`np.zeros_like` bit-exact.

```python
class IdentityOperator(LinearOperator[Domain]):
    def apply(self, x: Domain, /) -> Domain:
        return x
```
— `:1843-1844`. **Note: returns the input object itself, not a copy.**

```python
class ScaledOperator(...):
    def apply(self, x: Domain, /, *extra, **kwextra) -> Codomain:
        return self.scalar * self.op.apply(x, *extra, **kwextra)
```
— `:1758-1759`. With `op = IdentityOperator()` this is `α * x` — one fresh array.

**[R]** `ScaledOperator.__init__` **refuses** `scalar == 0.0`
(`:1729-1737`, `ValueError("ScaledOperator with zero scalar is degenerate; use
ZeroOperator explicitly.")`) — which is exactly why stage 2 has to branch on
`albedo == 0.0` first. The `α == 1.0` branch is a pure optimization/structural
collapse (an `IdentityOperator` advertises `is_invertible=True` and needs no
multiply); the `α == 0.0` branch is **load-bearing** — without it, a vacuum face
would raise at construction.

### The `IdentityOperator` aliasing hazard — measured

**[M]** `IdentityOperator.apply` returns the input **by reference**. In
`DiffusionBoundaryOperator.apply` the argument is `trace.outflow_view(face)`,
which **[R]** is documented as *"a VIEW. Writes propagate to the backing buffer"*
(`orpheus/transport/fields/scalar_boundary_flux.py:110-117`). The write target is
a **different** buffer (`ScalarBoundarySourceSink.zeros_on(self.mesh)`), so the
composite's `slot[...] = law.apply(...)` copies on assignment and no aliasing bug
occurs today. **[M]** Verified: after `B.apply(psi)`, the input trace buffer is
unchanged (`[1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8]` in, same out — §8).
This is safe **by the caller's discipline**, not by the operator's contract — a
future consumer doing `y = law.apply(view)` then `y *= 2` would silently mutate
the input field on a reflective face. Worth a note, not a defect today.

**[M]** Direct probe: `law.apply(view) is view` → **True**;
`np.shares_memory(law.apply(view), view)` → **True** for the reflective
(`IdentityOperator`) face. Input trace before/after `B.apply(psi)`:
`[1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8]` → `[1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8]`,
unchanged.

---

## 4. `DiffusionBoundaryOperator` — the `:544` vs `:407` question, RESOLVED

**[G]** `orpheus/diffusion/operators.py` (668 lines, `git status` clean at
`73627b71`) contains **exactly two `LinearOperator` classes** and one private
dataclass:

| line | symbol | kind |
|---|---|---|
| `:254` | `_FaceClosure` | `@dataclass(frozen=True)` — per-face geometry + P1 closure coefficients; **not** an operator |
| `:272` | `LeakageOperator(LinearOperator["FullField", "FullField"])` | `block_role = BlockRole.FULL` |
| `:544` | `DiffusionBoundaryOperator(LinearOperator["FullField", "FullField"])` | `block_role = BlockRole.BOUNDARY` |

**[R]** `__all__ = ["DiffusionBoundaryOperator", "LeakageOperator"]`
(`orpheus/diffusion/operators.py:160`).

### Verdict on the two prior claims

- **The `:544` claim is CORRECT.** `class DiffusionBoundaryOperator` begins at
  line 544.
- **The `:407` claim is WRONG — and instructively so.** **[R]** Line 407 is
  `    def apply(self, psi: "FullField", /) -> "FullField":` — the **`apply` of
  `LeakageOperator`**, not a class definition at all. The parent report's §1 line
  `DiffusionBoundaryOperator (diffusion/operators.py:407)` should read `:544`.
- **The `:613` claim is CORRECT.** **[R]** Line 613 is verbatim
  `                law.apply(trace.outflow_view(face))`.

### There is no "second boundary-ish class" — there is a second boundary-*touching* one

The confusion is real but structural, not a duplication. **`LeakageOperator` also
writes boundary rows** — it is the only diffusion leaf that couples bulk ↔ trace:

**[R]** `orpheus/diffusion/operators.py:419-428` (inside `LeakageOperator.apply`):

```python
        # ── Trace block: outflow-definition defect + inflow identity ──
        out_boundary = ScalarBoundarySourceSink.zeros_on(self.mesh)
        for face, c in self._face_closures.items():
            j_plus = trace.outflow_view(face)             # (ng,)
            j_minus = trace.inflow_view(face)             # (ng,)
            slot = out_boundary.face_view(face)
            slot[ScalarTraceSpace.OUTFLOW_ROW] = (
                j_plus - c.c_phi * phi[:, c.edge] - c.c_inflow * j_minus
            )
            slot[ScalarTraceSpace.INFLOW_ROW] = j_minus
```

So the trace rows are **split across the two operators by design**, per the
block table **[R]** `orpheus/diffusion/operators.py:34-56`:

```
    L        →  [ A_bb  A_bs ] (FULL — the elliptic sibling of streaming)
                [ A_sb  A_ss^L]
    B        →  [ 0     0    ] (BOUNDARY — the realized albedo block)
                [ 0     A_ss^B]
```

`L` writes **both** trace rows (the P1 outflow-definition defect **and** the bare
inflow identity `J⁻`); `B` writes **only** the inflow row (`𝒜 J⁺`); the loss's
inflow row `(L − B)` reads `J⁻ − 𝒜 J⁺ = 0`. **This is the exact SN streaming
convention mirrored** — `L` supplies the inflow identity, `B` supplies the
response.

**The composite that holds per-face laws is `DiffusionBoundaryOperator` (`:544`),
unambiguously**: **[R]** `self.face_laws: Mapping[str, LinearOperator] = dict(mesh.bc)`
(`:580`). `LeakageOperator` holds `_FaceClosure` bundles, never a law.

---

## 5. The composite's action — and the "discarded preservation" question

### How it walks faces

**[R]** `DiffusionBoundaryOperator.apply`, the entire post-guard body,
`orpheus/diffusion/operators.py:610-618`:

```python
        out_boundary = ScalarBoundarySourceSink.zeros_on(self.mesh)
        for face, law in self.face_laws.items():
            out_boundary.face_view(face)[ScalarTraceSpace.INFLOW_ROW] = (
                law.apply(trace.outflow_view(face))
            )
        return FullField(
            interior=ScalarSourceSink.zeros_on(self.mesh),
            boundary=out_boundary,
        )
```

Preceded by two guards (`:596-609`): the boundary must be a `ScalarBoundaryFlux`
(`TypeError` otherwise) and `psi.interior.mesh is self.mesh` (`ValueError` — the
mesh-identity invariant).

### The structural answer: diffusion does NOT have the SN discard

This is the sharpest finding of the quadrant. **The two composites differ in the
LAW'S DOMAIN, and that difference is the whole thing.**

| | SN | diffusion |
|---|---|---|
| composite | `SNBoundaryOperator._reflect_trace` | `DiffusionBoundaryOperator.apply` |
| what is handed to the law | **[R]** `face_in = boundary.face_view(face)` — the **WHOLE face slot** (all ordinates, in AND out) — `sn/operators/boundary.py:295, 301` | **[R]** `trace.outflow_view(face)` — the **outflow half only**, `(ng,)` — `diffusion/operators.py:613` |
| what the law returns | a full face slot (`full`) | the inflow value, `(ng,)` |
| what the composite writes | **[R]** `out_boundary.face_view(face)[sel] = full[sel]` with `sel = inflow_indices_for_face` — `boundary.py:302` | the whole return, onto `INFLOW_ROW` |
| **is anything discarded?** | **YES** — `full`'s outflow rows are computed then dropped | **NO** — nothing outside the inflow row is ever produced |

**Conclusion: the parent report's §4.1 finding ("the law preserves something the
composite discards, `P_in ∘ P_out = 0` by construction") has NO diffusion
counterpart.** The diffusion law is *already* typed `J⁺ ↦ J⁻`: outflow-in,
inflow-out. There is no complement subspace to preserve and nothing for the
composite to project away.

**[M]** Confirmed by measurement (§8): the composite's OUTFLOW row is `0` on
every face, and so is the law's output shape — the law never had an outflow
component. **[G]** Nothing in `orpheus/diffusion/` reads `B(ψ)`'s outflow rows;
the only consumer of `DiffusionBoundaryOperator.apply` is the `OperatorSum` in
`A = L + C − S − B`.

> **Read this as evidence for the parent report's §4.1 diagnosis, not against
> it.** SN's discard exists *because* SN's realized law has the wrong domain — it
> is typed `FullField → FullField` when the physics is
> `outflow-trace → inflow-trace`. Diffusion picked the honest domain and the
> discard vanished. **The SN fix is a domain narrowing, and diffusion is the
> worked precedent for it.**

### Consumers of the composite

**[G]** `grep -rn "DiffusionBoundaryOperator" orpheus/` — the only production
construction site is `orpheus/diffusion/solver.py:239`:

```python
        self.boundary = DiffusionBoundaryOperator(mesh)
        self.loss = self.leakage + collision - scattering - self.boundary
```
— `orpheus/diffusion/solver.py:239-240`. Everything else is docstring
cross-references (`__init__.py`, `augmented_mesh.py:312`, `method_space.py:31`,
`boundary_realizer.py:34`, `scalar_boundary_source_sink.py:20`) and two Sphinx
references (`docs/theory/methods/diffusion_1d.rst:288`, `:1348`). The loss then
feeds `MatrixInverseOperator(FlattenedOperator(self.loss, self.template))`
(`solver.py:249-251`) — i.e. `apply` reached by basis probing, **not**
`assemble()`.

---

## 6. Does diffusion factor as `R ∘ G`? — YES, and it proves it by NEGATION

### The affine form the whole subsystem is posed in

**[R]** `orpheus/geometry/boundary/_base.py:77-87` — `BoundaryTraceLaw` is
declared as *"Method-agnostic boundary law in the affine form
:math:`\gamma_- \psi = R\,G\,\gamma_+ \psi + q`"* with the three first-class
properties `geometry_map` (`G`), `response_kernel` (`R`), `source` (`q`).
(Line-number note for the parent report: at `73627b71` the `geometry_map`
property is `_base.py:148-155` and `response_kernel` is `:157-164`, not
`:157`/`:166`.)

### Diffusion's answer, factor by factor

| factor | diffusion's realization | evidence |
|---|---|---|
| **`R`** | **the `𝒜` stage — YES, it is exactly `R`.** Stage 1 (`_partial_current_albedo`) returns the scalar amplitude and nothing else. | **[R]** `boundary_realizer.py:197-249` |
| **`G`** | **the identity — DEGENERATE, and forced.** | **[M]** every realized operator is `Zero`/`Identity`/`α·Identity`; **[R]** `_albedo_operator` never sees a law and never reads `method_space` |
| **`q`** | **structurally excluded** — `PrescribedInflow` refuses. | **[M]** refusal message §2 |

### What "degenerate" means here, precisely

Not "diffusion is simple" — a *dimension* statement:

**[M]** A diffusion boundary face carries a `(2, ng)` slot: **one** outflow DOF
and **one** inflow DOF per group. Measured on the 2-group fixture,
`trace.layout.faces == {'xmin': (offset 0, shape (2, 2)), 'xmax': (offset 4,
shape (2, 2))}`, with `face_shape(axes, label) == ()` — no spatial extent on a
1-D face.

`G` is a linear map from the outflow space to the inflow space. Per (face,
group) both are **one-dimensional**, so `GL(1) ≅ ℝ*` — *every* candidate `G` is a
scalar, and any scalar is by convention absorbed into `R`. **`G = I` is not a
modelling choice; it is the only thing a 1-D→1-D geometric map can be once the
amplitude is factored out.** Contrast SN, where `G` acts on the ordinate axis:
specular reflection is a genuine non-trivial **permutation** of ordinates, white
is a genuine **rank-one contract-and-broadcast**. Those are the same `G` slot
with ≥2 dimensions to move around in.

Two corollaries fall straight out, both already in the code:

1. **`white ≡ reflective` at P1** is exactly "the two `G`s differ, but their
   difference lives entirely in the angular axis that the half-range ℓ=0 moment
   integrates away." **[R]** `boundary_realizer.py:60-66` says this in physics
   words: *"specular and Lambertian return differ only in the ANGULAR
   redistribution of the returned particles, which the half-range ℓ=0 moments
   integrate out — both preserve the returned current."* **[M]** Confirmed
   numerically: `WhiteBoundary(α)` and `ReflectiveBoundary(α)` realize to the
   *identical* operator class and value for α ∈ {1.0, 0.3/0.6}.
2. **`G` is NOT degenerate for the spatial axis** — and that is precisely where
   diffusion refuses.

### The negation proof — the two refusals ARE the two non-`R` factors

This is the strongest structural evidence in the quadrant. The realizer refuses
exactly two laws, and they are exactly the two that step outside `R`:

| refused law | which factor it needs | the refusal's own words **[M]** |
|---|---|---|
| `PeriodicBoundary` | a **non-identity `G`** | *"couples a face to its OPPOSITE face (**a trace-block permutation**), not a per-face albedo J⁻ = 𝒜·J⁺"* |
| `PrescribedInflow` | a **nonzero `q`** | *"J⁻ = q is the **rank-0 AFFINE** law — its realization is the boundary source q.boundary, not a linear boundary operator B"* |

`DiffusionBoundaryRealizer` therefore realizes **exactly the `{G = I, q = 0}`
corner of `γ₋ψ = R G γ₊ψ + q`**, and it refuses, with named errors, at both
walls of that corner. It has independently rediscovered the `R ∘ G ∔ q`
factorization by drawing its own boundary along the factorization's seams.

> **The answer to the driving architectural question:** the `R ∘ G` factorization
> is *real* and *load-bearing*, and diffusion is its cleanest witness — precisely
> because diffusion's `G` collapses, letting the `R` factor stand alone in a
> named, method-blind, `float`-valued function. The parent report's §2.3 finding
> that SN "fuses `R` as `α * base`" is the same factorization performed
> **inline**; diffusion performs it **as two functions**. Same math, one has a
> name for each factor.
>
> **The template for the whole subsystem** is therefore not "copy diffusion's
> two-stage realizer" (it works only because `G = I`), but: **give `R` and `G`
> each a named producer whose codomain is a specification, and let a
> method-blind collapse turn the pair into an operator.** Diffusion has the
> second half already; it never needed the first because its `G` is `I`.

### What is NOT there

**[G]** Zero uses of `OperatorProduct` / `@` / `TensorProductOperator` / `&` /
`SumOfTensorProductsOperator` / `RankOneOperator` anywhere in
`orpheus/diffusion/`. The `∘` is never spelled — consistent with the parent
report's §2.3 finding for SN. **[G]** `geometry_map` and `response_kernel` are
never read in `orpheus/diffusion/` either: diffusion re-derives `R` from
`law.albedo` inside its own isinstance chain rather than reading the declared
`response_kernel` property. This is the **sixth** production site answering the
parent report's §3.1 question ("is my `R` …?") in a non-property way, and the
strongest one — because here the answer literally IS a `float` named `albedo`,
i.e. the `response_kernel` property's value, computed by a private static method
instead of read off the law.

---

## 7. Adjoint / transpose — the leaves have it, the composite throws it away

**[M]** Measured on every realizable law
(`realize(law, DiffusionMethodSpace.minimal())`, probe input `x = [1.1, 1.2]`):

| law | realized class | `is_adjointable` | `apply_transpose(x)` | `.H` | `is_invertible` | `is_assemblable` |
|---|---|---|---|---|---|---|
| `VacuumInflow` | `ZeroOperator` | **True** | `[0. 0.]` | `_AdjointOperator` | False | False |
| `ZeroFluxBoundary` | `ScaledOperator(-1)` | **True** | `[-1.1 -1.2]` | `_AdjointOperator` | True | False |
| `ReflectiveBoundary(1.0)` | `IdentityOperator` | **True** | `[1.1 1.2]` | `_AdjointOperator` | True | False |
| `ReflectiveBoundary(0.3)` | `ScaledOperator(0.3)` | **True** | `[0.33 0.36]` | `_AdjointOperator` | True | False |
| `WhiteBoundary(1.0)` | `IdentityOperator` | **True** | `[1.1 1.2]` | `_AdjointOperator` | True | False |
| `WhiteBoundary(0.6)` | `ScaledOperator(0.6)` | **True** | `[0.66 0.72]` | `_AdjointOperator` | True | False |
| `AlbedoBoundary(0.4)` | `ScaledOperator(0.4)` | **True** | `[0.44 0.48]` | `_AdjointOperator` | True | False |
| `AlbedoBoundary(0.0)` | `ZeroOperator` | **True** | `[0. 0.]` | `_AdjointOperator` | False | False |
| `AlbedoBoundary(1.0)` | `IdentityOperator` | **True** | `[1.1 1.2]` | `_AdjointOperator` | True | False |

**Every realized diffusion boundary law is adjointable, unconditionally** — it
comes free because the leaves are generic numerics primitives
(**[R]** `ZeroOperator.is_adjointable → True` `operator.py:1945-1947`;
`IdentityOperator.is_adjointable → True` `:1857-1859`;
`ScaledOperator.is_adjointable → self.op.is_adjointable` `:1801-1804`).

### Then the composite drops it

**[M]** Measured:

| operator | `is_adjointable` | has `apply_transpose`? | `.H` | `is_assemblable` | `is_invertible` |
|---|---|---|---|---|---|
| `DiffusionBoundaryOperator` | **False** | **No** | raises `MissingAdjoint` | True | False |
| `LeakageOperator` | **False** | **No** | raises `MissingAdjoint` | True | False |

**[M]** verbatim: `MissingAdjoint: DiffusionBoundaryOperator is not adjointable —
.H/.adjoint() requires is_adjointable=True (a working apply_transpose on every
constituent).`

**[G]** `grep -rn "\.H\b|apply_transpose|adjoint" orpheus/diffusion/` returns
**exactly one line**, and it is prose: `operators.py:280`
*"elliptic-self-adjoint vs characteristic-triangular"*. **There is no adjoint
implementation anywhere in the diffusion package.**

### What is gated on what

- **Leaf level:** nothing is gated — `is_adjointable` is `True` by inheritance
  from the primitive.
- **Composite level:** `DiffusionBoundaryOperator` simply **does not declare**
  `is_adjointable` or `apply_transpose`, so it inherits the base default
  **[R]** `return False` (`orpheus/numerics/operator.py:622-634`). The gate is
  *absence*, not a condition.
- **Contrast SN:** the SN composite DOES implement `apply_transpose` and gates it
  on `_RULED_CORNER_KINDS` (**[R]** `sn/operators/boundary.py:303-315`: the
  per-face `if not adjointable(law): raise MissingAdjoint(...)` bridge, described
  as *"unreachable-in-practice because :attr:`is_adjointable` gates the composite
  eagerly"*). So SN gates on a **condition**; diffusion gates on a **hole**.

### The adjoint would be free — measured

The one hard part of a boundary adjoint (SN's cosine-weighted metric vs the
Euclidean transpose — parent report §2.4) **does not arise in diffusion**:

**[M]** The trace metric is a single scalar per face, identical on both the
outflow and inflow rows and across groups:

| geometry | faces | trace metric `G` entries |
|---|---|---|
| Cartesian slab (R=5) | `xmin`, `xmax` | `1.0` × 8 |
| Cylinder (R=5) | `xmax` only (pole is not a face) | `31.41592654` × 4 = 2πR |
| Sphere (R=5) | `xmax` only | `314.15926536` × 4 = 4πR² |

Because the metric block is a scalar within a face, the metric-correct adjoint
and the Euclidean transpose **coincide exactly**. **[M]** With
`G_full = diag(V_cell ⊗ 1_ng, A_face)` built by basis-probing
`full_field_space.apply_metric`:

```
||G⁻¹ [B]ᵀ G  −  [B]ᵀ||_max  =  0.0     (Cartesian, cylinder, AND sphere)
```

So `DiffusionBoundaryOperator.apply_transpose` is a **three-line addition** —
write `law.apply_transpose(...)` from the inflow row back onto the outflow row —
with **no metric question to resolve**. It is missing, not blocked.

**[G]** No production consumer of a diffusion adjoint exists today (the solver
uses `MatrixInverseOperator` on the flattened forward loss). Per the parent
report's §4 "declared-capability-ahead-of-consumption" caution, this is the
mirror-image defect: **capability that is trivially available and simply absent**,
with a latent consumer (adjoint diffusion for perturbation theory / DSA-adjoint /
the `AdjointSolution` family) already built on the SN side.

---

## 8. The measured end-to-end — mirror of the parent report's §4.1

**Setup [M]:** homogeneous 2-group, 4-cell non-uniform slab
(`edges = [0, 0.5, 1.5, 3.0, 5.0]`, `mat_ids = [0,0,0,0]`, the standard test
mixture); boundary trace **seeded strictly positive on ALL rows**,
flat buffer `[1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8]`
(`xmin` slot = `J⁺=[1.1,1.2]`, `J⁻=[1.3,1.4]`; `xmax` slot = `J⁺=[1.5,1.6]`,
`J⁻=[1.7,1.8]`); bulk `φ ≡ 2.0`.

### Config A — `xmin = albedo(0.4)`, `xmax = vacuum`

| face | law | `law.apply(J⁺)` | composite OUT **inflow** | composite OUT **outflow** |
|---|---|---|---|---|
| `xmin` | `ScaledOperator(0.4)` | `[0.44 0.48]` | **`[0.44 0.48]`** | `[0. 0.]` |
| `xmax` | `ZeroOperator` | `[0. 0.]` | **`[0. 0.]`** | `[0. 0.]` |

Composite OUT bulk = `[0 0 0 0 0 0 0 0]` — **[M]** exactly zero, as `BlockRole.BOUNDARY` requires.

### Config B — `xmin = zero_flux`, `xmax = reflective`

| face | law | `law.apply(J⁺)` | composite OUT **inflow** | composite OUT **outflow** |
|---|---|---|---|---|
| `xmin` | `ScaledOperator(-1.0)` | `[-1.1 -1.2]` | **`[-1.1 -1.2]`** | `[0. 0.]` |
| `xmax` | `IdentityOperator` | `[1.5 1.6]` | **`[1.5 1.6]`** | `[0. 0.]` |

Composite OUT bulk = all zero.

### Reading it against the SN table

The SN §4.1 table has a column *"law alone: outflow"* showing `1.653` preserved by
the law and `0.000` after the composite — the discard. **The diffusion table
cannot have that column.** `law.apply` receives `(ng,)` — the outflow half only —
and emits `(ng,)`, which IS the inflow value. `composite OUT inflow == law.apply(J⁺)`
**bit-identically on every row of both configs**. Nothing is computed and dropped.

The `composite OUT outflow = [0. 0.]` column is **not** a discard: it is the
`A_ss^B` block being zero on that row *by design* — `L` owns the outflow-definition
row (§4). **[M]** Confirmed by applying the composed loss leaf:

```
(L − B) INFLOW  row, xmin: expected J⁻ − 𝒜J⁺ = [0.86 0.92]  got [0.86 0.92]  ✓
(L − B) INFLOW  row, xmax: expected J⁻ − 0·J⁺ = [1.7  1.8 ]  got [1.7  1.8 ]  ✓
(L − B) OUTFLOW row, xmin: [1.279034 1.190633]  ==  L alone: [1.279034 1.190633]  ✓ (B contributes 0)
(L − B) OUTFLOW row, xmax: [1.608386 0.994030]  ==  L alone: [1.608386 0.994030]  ✓
```

### The assembled form — the same statement in 2 nonzeros

**[M]** `DiffusionBoundaryOperator(mesh).assemble()` on config A, `(16, 16)`:

```
nnz = 2
  (10,  8) = +0.4000       # xmin: J⁻(g0) row  ←  J⁺(g0) col
  (11,  9) = +0.4000       # xmin: J⁻(g1) row  ←  J⁺(g1) col
```
Flat layout: bulk `0..7` (`g*nx + i`); `xmin` slot at `8` = `[J⁺g0, J⁺g1, J⁻g0, J⁻g1]`;
`xmax` slot at `12`. The **vacuum face contributes literally zero entries** (the
`ZeroOperator`'s `as_matrix` is all-zero, so `np.nonzero` yields nothing). `[B]`
is group-diagonal, face-block-diagonal, and confined to the trace block — the
`A_ss` claim of the module docstring, measured.

**[M]** `max| as_matrix(B.apply) − B.assemble() | = 0.0` on Cartesian, cylinder,
and sphere — the two consumption modes agree exactly.

**[M]** The metric adjoint check (§7): `||G⁻¹[B]ᵀG − [B]ᵀ||_max = 0.0` in all three
geometries.

### Blockers encountered — none

The measurement ran without obstruction. Two constructor details worth recording
for whoever reproduces it:

- `DiffusionBoundaryOperator.apply` returns a **`ScalarBoundarySourceSink`**, not a
  `ScalarBoundaryFlux`. **[M]** `hasattr(out.boundary, 'inflow_view')` → `False`.
  The P1-dictionary accessors (`outflow_view` / `inflow_view` / `net_current` /
  `p1_boundary_scalar_flux`) are **flux-role-only** — a probe must read the output
  through `face_view(face)[ScalarTraceSpace.INFLOW_ROW]`. This is the role-typed
  algebra working correctly (an operator's output is a source/sink defect, not a
  flux), and it is a genuinely nice property of this arm.
- A curvilinear diffusion mesh has **one** boundary face (`xmax`); the pole `r=0`
  is not a face. **[M]** `list(mesh.bc) == ['xmax']` for cylinder and sphere.

---

## 9. Summary — the eight answers

1. **Two named stages: CONFIRMED.** `_partial_current_albedo` (law → `float`,
   `boundary_realizer.py:197-249`) then `_albedo_operator` (`float` → operator,
   `:142-157`), plus the shared `stamp_boundary_role`. The isinstance chain exists
   but its codomain is a **number**, and the operator-producing half is
   **law-blind and method-blind**.
2. **Registry admits 4 of 7 laws** (`vacuum`, `reflective`, `albedo`, `zero_flux`);
   refuses `white`/`periodic`/`prescribed_inflow`/`mixed` with the shared
   `TransportMethod` message. The **realizer** admits 5 (adds `WhiteBoundary`) and
   refuses 2 with its own named `BoundaryError`s. → parent §6 [U] raised to [M]:
   SN admits `{vacuum, reflective}` only.
3. **Full table in §3**, every `apply` body quoted: `0.0*x` / `x` / `α*x`. Three
   distinct classes, chosen by the *value* of 𝒜.
4. **`:544` is right, `:407` is wrong** (`:407` = `LeakageOperator.apply`). Two
   operator classes in the file; the law-holding composite is
   `DiffusionBoundaryOperator` at `:544`. `LeakageOperator` also writes trace rows
   — the outflow-definition defect and the inflow identity — which is the source
   of the "two boundary-ish classes" impression.
5. **`:613` is right**, and diffusion has **NO** discard structure: the law's
   domain is already the outflow half only. §5 has the side-by-side.
6. **Diffusion factors, and proves it by negation.** `𝒜` IS `R`; `G` is degenerate
   (1-D→1-D per face-group ⟹ `G = I` is forced, not chosen); `q` is excluded. The
   two refusals are exactly the non-identity-`G` law and the nonzero-`q` law.
7. **Leaves are all adjointable; the composite declares nothing** and `.H` raises
   `MissingAdjoint`. Zero adjoint code in the package. The adjoint is *free* here
   — measured, the metric adjoint equals the Euclidean transpose exactly, in all
   three geometries.
8. **Measured end-to-end in §8**, both configs, plus the assembled 2-nonzero
   matrix and the `(L − B)` row algebra.

### The one-line takeaway for the plan

> Diffusion is the subsystem's **only fully-factored boundary arm**, and it is
> fully factored because its `G` collapsed to `I`. It therefore validates the
> `R ∘ G` decomposition as *the* right frame — but it does **not** demonstrate
> that the two-stage realizer shape survives a non-trivial `G`. The transferable
> lesson is narrower and sharper than "copy this": **the stage that produces the
> operator must not see the law**, and **the stage that reads the law must
> produce a specification, not an operator**. Diffusion satisfies both because its
> specification is a `float`. SN's would be a pair (`G`-spec, `R`-scalar).

---

## 10. Corrections owed to the parent report

| parent report claim | status | correction |
|---|---|---|
| §1: `DiffusionBoundaryOperator` (`diffusion/operators.py:407`) | **[X] WRONG line** | `:544`. `:407` is `LeakageOperator.apply`. |
| §3: `geometry_map`/`response_kernel` declared at `_base.py:157`, `:166` | **[X] drifted** | at `73627b71`: `geometry_map` `:148-155`, `response_kernel` `:157-164`, `source` `:166-170`. |
| §4.1: *"diffusion has its own path feeding `trace.outflow_view(face)` (`diffusion/operators.py:613`)"* | **[R] CORRECT verbatim** | line 613 is exactly `law.apply(trace.outflow_view(face))`. |
| §6 [U]: *"SN's registry admits only vacuum and reflective"* | **RAISED to [R]/[M]** | `orpheus/sn/mesh/augmented_mesh.py:171-174`; pinned by `tests/sn/operators/test_snmesh_realizer_wiring.py:403`. Diffusion additionally admits `albedo` + `zero_flux`. |
| §1: realized operators are monomorphic | **[G] holds for diffusion** | zero `singledispatch` / `@overload` in `orpheus/diffusion/`. |

### Incidental findings outside this quadrant (offered, not claimed as scope)

- **[G] A stale comment in the SN registry.** `orpheus/sn/mesh/augmented_mesh.py:182-185`
  reads *"The 5 other kinds the realizer handles today (``white``, ``periodic``,
  ``albedo``, ``prescribed_inflow``, ``mixed``)"* — but **`mixed` names a law that
  no longer exists**: **[G]** `grep -rn 'key="mixed"|MixedBoundary' orpheus/`
  returns only three *deletion notices*
  (`sn/boundary/realizer.py:69`, `geometry/boundary/__init__.py:209`, `:350`),
  all saying `MixedBoundaryOperator` was *"deleted in Wave 11"*. So the comment
  advertises 5 kinds where only 4 laws exist. **[M]** `BC("mixed")` refuses at the
  registry, so nothing is broken — but the comment is a stale capability claim of
  exactly the family the parent report's §4 is cataloguing.
- **[M] 4 conflicting-V&V-marker warnings** in the diffusion boundary suite
  (`tests/diffusion/test_boundary_realizer.py:217, 228, 242, 254` — each
  `PytestUnknownMarkWarning: ... has conflicting V&V level markers
  ['foundation', 'l0']; using 'l0'`). Relevant to the parent's §8 item 4.
- **[G] The `-O` assertion-stripping exposure is NOT present in this arm's
  production code**: `grep -rn "^\s*assert " orpheus/geometry/boundary/
  orpheus/diffusion/` returns **nothing**. (This says nothing about the *test*
  side, which is the other quadrant's question.)
- **[M] Suite status:** `python -O -m pytest tests/diffusion/ -q` →
  **112 passed, 6 warnings, 2.19 s**. `test_boundary_realizer.py +
  test_operators.py` alone → **70 passed, 0.60 s**.
