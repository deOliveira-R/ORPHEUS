# P6 — promote the MMS's hand-rolled spec into production

Phase P6 of `.claude/plans/affine_boundary_source_channel.md`, the last one.
Requirements list: `scratch/p4_mms_design.md` §10 (items 1–6; **item 7 is VOID**
— the declaration channel `985497b5` dissolved the delivery problem).

> **Status: DESIGN IN PROGRESS.** Measurements below are mine and live; the
> cross-method survey (`scratch/p6_spec_survey.md`) is what decides WHERE the
> type lives, and until it lands this file does not propose a home.

---

## 1. The contract today

`orpheus/geometry/boundary/_source.py`:

```python
@runtime_checkable
class InflowSourceSpec(Protocol):
    def evaluate(self, shape: tuple[int, ...]) -> np.ndarray: ...
```

Two shipped implementors, `NoSource` (zeros) and `ConstantInflowSource` (a
scalar broadcast). Both are shape-only by nature, which is exactly why the
signature looked sufficient for five years of these two.

The Protocol's own docstring states the design intent that P6 has to overturn:

> "The source is a property of the boundary, not of any trace space, so it takes
> the *shape* it must fill rather than a trace object."

That sentence is the thing under review. It is defensible for a constant; it is
false for anything whose value depends on **which ordinate** or **which face**.

---

## 2. `[M]` The five gaps, from the P4 stopgap's constructor

`_ManufacturedFaceInflow` (`tests/sn/verification/analytical/test_mms_declared_inflow.py`)
needs five things `evaluate(shape)` cannot give it. This is not a wish list —
each is a constructor argument that exists solely to work around the signature:

| # | what it needs | why the shape cannot give it |
|---|---|---|
| 1 | the per-row **μ in trace order** | `shape[0]` is a COUNT; row `i` is ordinate `inflow_indices_for_face(face)[i]` and nothing says so |
| 2 | **which face** it is on | the trace value depends on `x_face`; `evaluate` receives no face identity, so one instance must be built per face |
| 3 | **which index space** `shape[0]` lives in | the SAME object is called at `(N,)` and at `(\|Γ₋\|, ng)`. Today a spec disambiguates **by rank** — an accident, and it breaks outright for any face where `\|Γ₋\| == N` |
| 4 | what the **trailing axes mean** | slab `(\|Γ₋\|, ng)`, 2-D `(\|Γ₋\|, ng, ny)`. A spec that filled `(12, 2, 6)` as `(ordinate, y, group)` passes every shape check and is silently wrong |
| 5 | the **`1/W` convention** | `γ₋ψ` is a flux, so it carries `1/W`. Nothing in the Protocol says whether `evaluate` returns a flux or a `W`-normalised density, and `[M]` getting it wrong is a **×2652** error — numerically indistinguishable from a double delivery on Gauss–Legendre |

**Every one of the five is information the face's incoming trace space already
holds.** That is the design, and it is not a coincidence: the spec's job is to
produce a vector *on* Γ₋, so Γ₋ is its natural argument.

---

## 3. ⛔ `[M]` Item 6 is not a separate defect — it is a SYMPTOM of item 3

Measured live today (`$CLAUDE_JOB_DIR/tmp/p6_probe_hole.py`), het slab, GL-8, 2G:

```
declines the probe     realized OK; q_d delivered = 7.0   shapes asked for: [(8,), (4, 2)]
honest                 realized OK; q_d delivered = 7.0
```

`assert_source_lives_on_incoming_trace` (`_base.py:440`) opens

```python
probe = self.source.evaluate((int(quadrature.N),))
if not np.any(probe):
    return
```

so a source that answers the rank-1 probe with zeros and the rank-2 delivery
with `7.0` **skips the ERR-047 certification entirely** and delivers anyway. A
presence predicate the subject can decline is not a guard.

The guard's own docstring justifies the probe like this:

> "the `InflowSourceSpec` contract ('fill exactly the requested shape') makes
> the probe representative of any trailing-axis block the realizer later
> requests."

⭐ **That is a claim about a contract the Protocol does not enforce.** Nothing
stops a source varying with rank — and the P4 stopgap *has* to vary with rank,
because rank is the only signal distinguishing the probe from the delivery
(item 3). The guard and the defect share one root cause.

⟹ **Do not patch the probe.** If `evaluate` receives Γ₋:

* there is only ONE call site and ONE shape, so item 3 disappears and with it
  the rank-sniffing that made declining possible;
* `q ∈ Γ₋` holds **by typing** rather than by a presence check, so the guard has
  nothing left to certify and is *retired*, not fixed.

This is `coding-elegance` Pattern 4 (illegal states unrepresentable) and the
retirement rule in one move: the elegant fix deletes the guard instead of
hardening it.

⚠ **Before deleting, answer one question the survey is checking:** is the
guard's RAISE arm (`inflow_indices is None`) reachable at all post-B3.4a? If it
is already dead, this is a retirement of dead code plus a hole; if it is live,
the retirement owes a replacement for whatever still reaches it.

---

## 3b. `[M]` The survey — done by hand (the dispatched explorer died on a 529
## having written **zero bytes**, despite an incremental-write brief)

**Q1 — the call sites of `spec.evaluate(...)` number exactly TWO:**

| site | shape | what it is |
|---|---|---|
| `geometry/boundary/_base.py:440` | `(N,)` | the ERR-047 presence probe |
| `transport/source_sinks/angular_boundary_source_sink.py:402` | `(\|Γ₋\|,) + slot_shape[1:]` | the real delivery |

⟹ **deleting the probe leaves ONE call site and ONE shape.** Item 3 (rank as the
disambiguator) does not need fixing; it ceases to exist. Every other
`.evaluate(` hit in `orpheus/` is a different concept (basis, functional,
density).

Note what the delivery site already computes:
`wanted = (|Γ₋|,) + slot_shape[1:]`, with `inflow = trace.inflow_indices_for_face(face)`.
**The caller already holds the trace and the face** — it just doesn't pass them.

**Q2 — the concept IS cross-method, and diffusion has already written down
why it refuses it.** `orpheus/diffusion/boundary_realizer.py` refuses
`PrescribedInflow` with an unusually honest comment:

> "nothing about P1 makes a prescribed inflow unrepresentable: q is a VECTOR in
> Γ₋, and Γ₋ here is one number per face, which the trace carries perfectly
> well. This refusal disappears the day #290 P5 wires the fixed-source arm."

⟹ the future second consumer is **named and dated**, and its Γ₋ is *one number
per face* — no ordinates. So whatever `evaluate` receives must not presuppose an
angular trace. (Also there: the deliberate TYPE test rather than a source-VALUE
test, "measured while writing B2 — it is the trap in this collapse".)

**Q4 — the guard has ONE production caller**, `sn/boundary/realizer.py:757`,
inside `assert_realizable`, passing `inflow_indices=method_space.inflow_indices`.
The raise arm therefore fires only when a method space carries no inflow set —
reachable on hand-built spaces (the realizer's `VacuumInflow` arm checks for
exactly that six lines later), so the arm is **live, not dead**. A retirement
owes it a replacement.

**Q5 — the object exists and is `AngularFaceTraceSpace`** (G6.1's Γ ladder),
built by `trace.inflow_space(face)`; P3 already binds `ZeroOperator`'s codomain
to it. Diffusion's counterpart is `ScalarTraceSpace`. `[M]` its public surface:

```
face='xmin'  role='inflow'  ordinate_indices=[4 5 6 7]
inner_product_weights=[0.0665 0.1649 0.1772 0.0972]   shape  name  dual
apply_metric  apply_inverse_metric  inner_product  norm
```

⛔ **It is necessary but NOT sufficient.** It answers items 2 (face *name*), 3
and 4 outright — but it carries **no μ**, **no quadrature**, and **no face
COORDINATE**. So "just pass the space" is half a fix: a μ-dependent or
`x_face`-dependent source would still need those baked in at construction, which
is the stopgap's tax (R1) surviving under a new signature.

That gap is the design decision, and it is recorded in §4 rather than resolved
here.

---

## 3c. ⭐⭐ THE REFRAME (user, 2026-08-06) — and the collapse it produced

The three options this file was about to offer were all over-built, because the
question was posed wrong: *"which metadata bag do we hand the source?"* The user
rejected all three and asked the right question — *does it need the quadrature
at all, or can it project at instantiation time?*

**`q` is an element of Γ₋(f). There are two ways to specify an element of a
discrete space: extensionally (the numbers, in row order) or intensionally (a
function on the continuum, plus a rule for discretizing it).**
`evaluate(shape) -> ndarray` is a **half-extensional** interface: it demands the
numbers while refusing to say which row is which. That is the single root of
every gap in §2. `ConstantInflowSource` survives it only because a constant is
row-independent.

⟹ the five gaps are not five metadata items. They are **the discretization map,
itemized** — and a discretization map is an OPERATOR:

.. math::  q_{\Gamma_-} \;=\; \gamma_S \circ E[f]

with `f` the user's function, `E` evaluation at the face's angular nodes, and
`γ_S` the restriction — **`TraceRestrictionOperator`, already minted at B3.1,
whose index set IS `space.ordinate_indices`.** `[M]` Γ₋ even has a second,
independent spelling in the measure algebra that agrees exactly:
`measure.restrict(λ μ: μ > 0)` gives nodes/weights bit-identical to
`quad.nodes[S]` / `quad.weights[S]`.

### Re-scoring the five gaps under the reframe — only ONE is real

| # | verdict |
|---|---|
| **1 μ in trace order** | ⭐ **the only genuine signature gap.** Depends on the quadrature AND the trace layout; neither exists at declaration time |
| 2 face coordinate | **not a gap.** Declaration is already per face (`bcs=(law_xmin, law_xmax)`), so `x_face` is declaration-time knowledge. One instance per face is the right granularity, not a tax |
| 3 which index space | **dissolves** — with the probe deleted there is one call site and one shape |
| 4 axis semantics | `[M]` `space.shape` is already `(12, 2, 6)` on 2-D — exactly the `wanted` the delivery site recomputes by hand |
| 5 `1/W` | **not a gap.** `1/W` belongs to the *ansatz*, which the author chose. It is an undocumented CONTRACT, fixable in prose. ⚠ `[M]` full `W = 2.0` vs inflow-restricted `W = 1.0` — handing a source only the restricted measure would make it wrong by 2×, silently, in the direction that mimics a double delivery |

⟹ **the spec should carry neither the quadrature nor μ.** It carries a function;
the realizer supplies the directions.

---

## 4. The plan (steps 2–4 open; step 1 LANDED)

**Step 1 ✅** — `AngularTraceSpace.directions` (= `quadrature.nodes`) and
`AngularFaceTraceSpace.directions` (the parent's, restricted by the tier's own
`ordinate_indices` — the SAME expression the metric already uses). `[M]` zero
direct `AngularTraceSpace(...)` constructions in the tree, so a required field
was a single-site change. Dimension-generic: `(m,)` in 1-D, `(m, 3)` in 3-D.
Four gates, mutation-verified: reversing the row order reds 9/9 row-order rows
and correctly leaves the inward-pointing row green (reversal preserves the set);
dropping the restriction reds 15.

**Steps 2 + 3 + 4 ✅ — and they were NOT separable.** ⚠ A planning surprise
worth its own note: the probe is one of the two `evaluate` call sites, so
changing the signature *forces* the guard's retirement in the same change.
Landing "step 2" alone would have left the probe calling `evaluate((N,))` on a
source expecting a space. **A step boundary that cuts across a signature's call
sites is not a step boundary** — the unit of work is the call-site set, not the
conceptual tidiness of the description.

* `InflowSourceSpec.evaluate(space)`. `NoSource` → `np.zeros(space.shape)`;
  `ConstantInflowSource` → `np.full(space.shape, value)` (unchanged in
  behaviour, and therefore the migration's regression control).
* The delivery site's hand-rolled `wanted = (|Γ₋|,) + slot_shape[1:]` deleted —
  it is `space.shape`.
* The probe deleted. `assert_source_lives_on_incoming_trace` →
  **`assert_source_is_placeable`**, because the old name outlived its claim: it
  no longer certifies where the source's values live, only that Γ₋ is
  nameable. The `quadrature` parameter went with the probe (nothing read it).
* Contract written down on the Protocol: `evaluate` returns the angular flux on
  Γ₋ in ψ's units, and the `1/W` is the author's — the FULL weight sum, `[M]`
  2.0 vs the inflow-restricted 1.0 on GL-8.

⭐ **Behaviour change, deliberate and strictly safer.**
`PrescribedInflow(ConstantInflowSource(0.0))` with no inflow set now RAISES,
where the probe let it through. The check can no longer ask "is this value
zero?" — in the arm where the answer matters there is no Γ₋ to evaluate
against, which is the thing being complained about. So the discriminator is the
source's IDENTITY (`NoSource` or not). Value-at-check-time was never
value-at-delivery-time, and that gap *is* ERR-047's mechanism. The typed escape
for a genuinely sourceless prescribed law is unchanged: `PrescribedInflow()`
defaults to `NoSource`.

### `[M]` The payoff, measured against the phase's own done-when

`_ManufacturedFaceInflow`'s constructor went **5 args → 2**, and both survivors
are declaration-time:

| before | after |
|---|---|
| `mu_inflow` (per-row μ in trace order) | `space.directions` |
| `n_ordinates` (to recognise the rank-1 probe) | *gone — there is no probe* |
| `A_g`, `B_g`, `W` | `case`, `x_face` (the ansatz and where to evaluate it) |

and the driver's throwaway **probe `SNMesh` is deleted** — the two-step that
existed only so the spec could learn the row order. The done-when predicate
("no `mu_inflow=`, no `n_ordinates=`") greps clean; the only residual mentions
are past-tense history.

**G8 re-posed, and it now asserts the opposite of what it did.** It pinned that
the spec is called at TWO shapes and that its rank-1 answer is non-zero — both
properties of the defect, not of a contract. It now pins **exactly one call,
against `Γ₋(f)` by name**, plus that the delivered `q_∂` is what the spec
returned for that space (so the single call is the real one, not a discarded
rehearsal).

⚠ **Diffusion.** `evaluate(space)` typed on the shared `FunctionSpace` means
`ScalarTraceSpace` fits the slot but has no `directions`, so a directional
source there would fail with a bare `AttributeError`. Note it; let #290 P5 add
the typed refusal when diffusion's fixed-source arm actually lands, rather than
build a guard for a caller that does not exist.

---

## 4b. Superseded — the three options this file was going to offer

1. **Where does the type live?** `InflowSourceSpec` is in
   `orpheus/geometry/boundary/` — the SHARED tier, consumed by the diffusion
   realizer too. Whatever object `evaluate` receives must therefore be
   method-agnostic, and diffusion's trace has no ordinates. Either the argument
   is a common supertype of the two traces, or the Protocol splits. **This is
   the decision; it is not mine to guess from the SN side.**
   (`feedback-extract-model-independent-primitive`: grep the sibling families
   BEFORE minting.)
2. **What is the object called?** G6.1 minted "the three-tier Γ ladder as
   spaces". If a per-face inflow space already bundles (face identity + inflow
   rows + the ordinates), the spec receives THAT and nothing new is minted.
   Reuse beats minting (`coding-standards`, type-vs-property).
3. **Blast radius.** Every implementor and every call site of `.evaluate(...)`,
   in `orpheus/` and `tests/`.

---

## 5. Refuted candidates (with the structural reason each fails)

| # | candidate | why it fails |
|---|---|---|
| **R1** | **Keep `evaluate(shape)` and pass the extras at CONSTRUCTION** (what the P4 stopgap does) | It works — and it is why the stopgap exists — but it pushes trace knowledge onto every *caller* that builds a spec, so the "one spec instance per face" tax is permanent and the trace-order contract stays unwritten. The stopgap is evidence of the gap, not a design. |
| **R2** | **Widen `evaluate` with optional keyword args** (`mu=None, face=None, …`) | A Protocol whose arguments are all optional pins nothing: an implementor may ignore them and a caller may omit them, so every gap above stays *expressible*. Optional context is not a contract. |
| **R3** | **Patch the probe to evaluate at BOTH shapes and OR the results** | Hardens item 6 while leaving item 3 (rank as the disambiguator) in place — it makes the accident load-bearing by writing a second consumer of it. Treats the symptom. |
| **R4** | **Drop the guard with no replacement, keep `evaluate(shape)`** | Removes the check without removing what it checks. `q ∈ Γ₋` would then rest on nothing at all. |
