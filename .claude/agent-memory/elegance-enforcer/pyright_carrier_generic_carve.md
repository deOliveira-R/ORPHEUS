---
name: pyright-carrier-generic-carve
description: Review rulings for the #226 pyright carrier-generic carves (C2 FunctionSpace, C3 transport/sn) — generic-surface/Any-realization split, default=Any, ratchet-baseline gate, stale-ignore-after-required-field smell
metadata:
  type: project
---

# C4 (composition-leg generics + typed-field role parses + `_CellSolve` split) — reviewed 2026-07-03, PASS-WITH-NITS

Uncommitted diff on `refactor/pyright-burndown`; ratchet sn 82→60, transport 1 (the
60 residual = the two C5 pole files only). ZERO VIOLATIONS, zero code defects, 6 casts
retired, 0 casts/ignores added. Reusable rulings (the C5 pole-file review + any future
`Generic[..., LegType]` composition retype inherit these):

- **The covariant-leg keystone is PRINCIPLED, not ceremony — and the chain is FORCED.**
  `OperatorSum` → `Generic[Domain, Codomain, SummandA, SummandB]` (also OperatorProduct
  `FactorA/FactorB`, ScaledOperator `ScaledOperand`), COVARIANT legs, PEP-696
  `default="LinearOperator[Domain, Codomain]"`, legs = read-only `@property` over
  `self._a: Final`. The forced chain: retire accessor casts → pin legs as type params → a
  pinned subclass (`InvertibleOperator = OperatorSum[FF,FF,Streaming,Mult]`) must upcast to
  the defaulted base because `StreamingOperator.__add__`'s overload returns the base spelling
  → upcast requires COVARIANCE → covariance sound ONLY without a setter → read-only property
  over `Final`. VERIFY soundness, don't trust: (a) grep NO writers to `.a`/`.b`/`.op`
  tree-wide (read-only conversion runtime-safe); (b) the leg typevar escapes only via the
  property (output) + `__init__` (the standard pyright/mypy covariant-in-ctor exception) →
  genuinely sound, not silenced. **`default=covariant-leg-type` is CORRECT here — contrast
  C2's forced `default=Any`: C2's `FunctionSpace` was INVARIANT; covariance lets the narrow
  default work.** Concept-count DOWN (6 cast-assertions → structural typing). The keystone
  GENERALIZES: #289 (FullField leaf-generics) cites it as the resolution template. Payoff
  visible: `cast("X", cast(ScaledOperator, self.b).op)` → `self.b.op`.

- **Operator-leg `__init__` guard vs static leg-pinning = NOT a twin (parse-once/propagate).**
  `InvertibleOperator.__init__`'s `isinstance(streaming, StreamingOperator)` (pre-existing) is
  a boundary PARSE that runs for untyped/`Any` callers and is ERASED at runtime →
  irreplaceable by the static type; the static pinning propagates the invariant cast-free →
  irreplaceable by the guard. Complementary (Pattern 4), two distinct jobs. The discriminator
  when asked "static type + runtime guard = duplication?": the runtime guard survives
  type-erasure and serves untyped callers; the static type gives compile-time knowledge. Each
  loses something if removed.

- **The typed-field ROLE parses (`FullField.bulk`/`.boundary` role) are a DIFFERENT family at
  a type-erasure seam — honest, #289-tracked, NOT validate-twice.** `FullField.bulk: BulkField`
  is role-generic (role varies by composite → unlike C3's mesh it CANNOT get a typed accessor).
  So C4 reifies each demand with a local `isinstance` parse (loud TypeError) that pins the
  runtime ravellable-template-echo contract pyright can't see. 4× `isinstance(bulk, AngularFlux)`
  + AngularField/HarmonicMomentFlux/TimedFullField/BoundaryFlux variants across solver/solution/
  boundary. Correctly left un-extracted (a helper would be DELETED by #289's generics — the
  keystone above is #289's template) per Pattern-6/L-002 (defer, cross-ref, track). VERIFY the
  tracker exists: #289 OPEN, inventories the sites exactly. Do-now = cite `#289` by NUMBER in
  each parse comment (greppability). CORRECTNESS CHECK I ran: `AngularFlux` is imported LOCALLY
  in each fixed-source fn (2083/2270) before use (2179/2322) — no NameError; the parses fire safe.

- **`_CellSolve` Optionals+XOR-`__post_init__` → ABC + 2 kw_only subclasses = the RIGHT cut
  (Pattern 4, tagged-union would be WORSE).** 4 `Optional` buffers (16 combos, 2 legal) + runtime
  XOR guard → `_CellSolve(ABC)` shared `cell()` kernel + abstract `_emit`; `_CellSolveAngular`/
  `_CellSolveMoment` each carry 2 REQUIRED buffers. 14 illegal combos + guard DISSOLVE
  (unrepresentable by construction). "A repeated conditional is a missing type" (`if moment_buf
  is None` was in BOTH `cell()` and `_SweepEmit.pure_z`). A tagged-union `mode` slot pushes the
  branch back INSIDE `_emit` → re-introduces the conditional. Expression-problem: cases stable +
  shared kernel + forked emit → polymorphism (ABC over Protocol because `cell()` is a real ~40-line
  shared IMPLEMENTATION). `kw_only=True` REQUIRED (base has a defaulted field, subclasses add
  required → else "non-default follows default"). Verify no bare `_CellSolve(` construction
  survives (abstract→TypeError); test rewires `_CellSolve(`→`_CellSolveAngular(` = retirement-as-
  test-migration. ASYMMETRY NOTE: `_SweepEmit` kept Optionals+guard + got `angular_buffers()`/
  `moment_buffers()` accessors — the natural next split, but pre-existing + coupled through the
  ONE walker dispatch (no drift habitat) → forward-note, not required.

- **The `range(d)`-loop → `(det free + 1 determined)` restructure removes a slice/ndarray UNION
  access AND false generality.** `write` built per-axis in a loop that re-branched `box_contiguous`
  → now each branch builds `read` AND its `+1`-`shifted` counterpart whole; `write = (*(...) for a
  in range(det)), read)`. Semantically identical IFF `det == d-1` (one determined axis, LAST) —
  VERIFIED: `det = d - 1` is a single-sourced invariant (sweep_graph.py:647, documented SSOT with
  the non-last generalization path noted). The old loop's `a >= det` arm only ever fired once given
  det=d-1 → the new form is honest about the real invariant, not a narrowing. Bit-id anchored by the
  window≡full oracle (test_sweep_graph_window_equivalence).

- **THE recurring hard gate (verified C4 too): the ratchet baseline companion.** `test_pyright_
  ratchet.py:63,74` reds on IMPROVEMENT ("DECREASED ... sn: 82 -> 60"). C4's
  `tests/_harness/pyright_baseline.json` was UNMODIFIED (`"sn": 82`) in the working tree → the carve
  reds `main` until `python -m tests._harness.pyright_ratchet --update` re-tightens it (sn 82→60)
  and the baseline commits WITH the code. ALWAYS `grep git status` on the baseline; the diff's
  "82→60" prose runs ahead of the tree until it lands.

- **`full_field.py __add__/__sub__: other: T → other: "FullField"` = anti-#20 fixed HONEST
  direction.** `other: T` was the LIE (claimed same-flavor partner); `_check_partner` (read it —
  it accepts ANY FullField, load-bearing for the timed−timeless time-derivative stencil per its own
  docstring) proves "any FullField flavor" is the real contract. The type now matches runtime. When
  a param annotation is NARROWER than the runtime guard, the annotation is the bug, not the guard.

- **NITs that recur on these carves:** sibling builders typed asymmetrically (`_within_group_krylov`
  got a return type, `_within_group_si` didn't); a narrowed member accessed twice inline instead of
  bound-to-a-local-once (fragile narrowing-persistence); dead-fn code tombstone duplicating a Rule-3
  docs retirement note (trim to one line). All NIT-tier.

---

# C3 (transport slot-typing + F2 layering carve) — reviewed 2026-07-02, PASS-WITH-NITS

Uncommitted diff on `refactor/pyright-burndown`; ratchet transport 24→1, sn 86→82.
ZERO VIOLATIONS, zero code defects. All 4 user-flagged scrutiny points cleared:

- **(a) `FullField.mesh` reads the mesh off the BOUNDARY leaf — HONEST, not a trick.**
  The composite's `__post_init__` enforces `bulk.mesh IS boundary.mesh` (identity),
  so reading either leaf yields the same object; the property picks `boundary` because
  `BoundaryField.mesh: SNMesh` is a HARD declaration (line 685) while `BulkField.mesh`
  was widened to `MaterialMesh` (#267) and `AngularField` re-narrows it. So `psi.mesh`
  is the SNMesh the operators need, NO cast — retires 4 `cast("SNMesh", psi.bulk.mesh)`
  sites. The REUSABLE ruling: a composite that carries a mesh-identity invariant should
  expose it as a property reading the leaf whose STATIC type is narrowest — that is
  single-source (don't store a duplicate mesh field) AND cast-free. Forward-coupling
  (when a 2nd method widens `BoundaryField.mesh`, this widens too) is pyright-visible at
  the operator consumers → defer, not a blocker.
- **(b) Closed-grid `@overload`+isinstance BEATS `singledispatchmethod` — not lateral.**
  `singledispatchmethod`'s return type is the base method's (can't vary per registered
  type → callers got `object`); typed `@overload` faces give PRECISE per-carrier returns
  AND make an off-grid carrier a call-site type error. For an architecturally CLOSED set
  (the 2×2 carrier grid) the isinstance chain is the coding-elegance-blessed "exhaustive
  match over a closed sum type." The singledispatch form ALREADY had the same per-arm body
  duplication, so isinstance adds no duplication and gains static precision. Ruling: prefer
  overload+isinstance over singledispatch when the dispatch set is closed.
- **(c) `#288` accept-untyped is the RIGHT defer.** The #207 cross-class source-sink
  `__add__` (iso + per-ordinate → the LARGER class) is statically unspellable against
  `Field.__add__`'s `(T,T)->T`. The three candidate resolutions (named-composition rewire /
  NotImplemented protocol / accept-untyped) are DIFFERENT ARCHITECTURES with call-site
  ripple → a user decision (#288, OPEN, verified), NOT a mechanical typing fix a burn-down
  cluster should unilaterally pick. Done honestly: dunders left untyped (NOT `# type: ignore`),
  WHY documented in both docstrings, the one residual consumer error left VISIBLE in the
  baseline (not silenced). This is anti-pattern #17 done right (a pinned, tracked gap).
- **(d) Fission composite-arm narrowing to AngularFlux-only is HONEST.** Grep proved NO
  production path couples fission to a windowed `HarmonicMomentFlux` bulk (windowing is a
  within-group SI/scattering construct); `HarmonicMomentFlux` has `.scalar_flux`, not
  `.integrate_angular`, so the OLD code AttributeError'd on that (unreachable) path — the
  new guard converts it to a clear parse-boundary TypeError (Pattern 4). No production
  behavior change (AngularFlux path identical). The asymmetry with scattering's arm (which
  DOES accept the moment bulk, scattering.py:1553) is genuine and documented. A test would
  pin an unreachable defensive guard → optional NIT, not required.

## THE novel reusable smell: stale `# type: ignore` surviving a required-field carve

When a carve makes a formerly-`Optional` field REQUIRED (here `TraceSpace.omega_dot_n:
NDArray`, retiring the `_face_row`/`face_names` None-guards — good Pattern-4 work), it must
SWEEP the `# type: ignore[...]` comments that existed ONLY because the field was Optional.
C3 left two behind (`trace_space.py:406,415` — `self.omega_dot_n[...]  # type: ignore[index]`).
- **Non-destructive proof they're stale:** the SAME index op on a LOCAL `NDArray` at line 281
  carries NO ignore. Same type now → the field-access ignores are dead.
- **Why it MATTERS / bug habitat:** a dead ignore is a suppressed static check (anti-pattern
  #19) — it will swallow a REAL future index error at that line. And it is INVISIBLE to the
  ratchet: `reportUnnecessaryTypeIgnoreComment` is NOT in `[tool.pyright]`, so pyright never
  flags it and the burn-down can't self-heal it. This is precisely the elegance-review
  value-add over the mechanical gate. Recurring check for EVERY Optional→required carve in
  this campaign: grep the file for `type: ignore` and prove each survivor still bites.

## The retired-vs-kept cast split (the C3 F2 ledger) is architecturally sound

4 mesh-casts RETIRED via `psi.mesh` (composite has a narrow-typed boundary leaf to read).
1 KEPT: `ScalarSourceSink.as_per_ordinate` (scalar_source_sink.py:201) casts its own
`MaterialMesh` up to `SNMesh` — a BARE scalar source (no composite leaf), method-agnostic by
#267/#276 (narrowing its `.mesh` would break the homogeneous meshless solver). Irreducible
true-boundary cast → KEEP-and-document is correct ("casts only at true boundaries"). Ruling:
the F2 "read mesh off the boundary leaf" technique only applies to COMPOSITES; a bare
method-agnostic field projecting into SN space keeps its boundary cast.

## Doc-drift the carve introduced (NIT): retire-the-mechanism → retire-its-docstring

`harmonic_frame.py:40` module docstring still cites `:func:`~functools.singledispatchmethod``
as the dispatch mechanism the diff just RETIRED. The `:func:` xref still RESOLVES (real stdlib
symbol) so Sphinx `-W` won't warn — it just misdescribes the code (anti-pattern #20). The
method-level comment WAS updated; the module docstring was missed. Recurring: when a carve
swaps a named mechanism, grep the whole file's prose for the old mechanism name.

---

# Carrier-generic carve (pyright #226, cluster C2 — `FunctionSpace[Carrier]`)

Reviewed 2026-07-02 @ `refactor/pyright-burndown` (uncommitted). Verdict
PASS-WITH-NITS, one blocking approval condition (baseline re-tighten). Files:
`orpheus/numerics/space.py`, `orpheus/numerics/spaces/full_field_space.py`.
Zero code defects, zero VIOLATIONS. These rulings are REUSABLE for the sibling
generic carves (C3 operator-slot `FunctionSpace[Domain]`, C4 sn, and any future
`Generic[T, default=Any]` two-param-house-pattern retype).

## The pattern that recurs: generic public surface + `Any`-typed private realization

pyright (1.1.410+) runs the **override compatibility check against the BASE
signature regardless of self-specialization** — a `self: "FunctionSpace[NDArray]"`
scoped body does NOT make a composite override compatible. The clean fix:
- Make the base `Generic[Carrier]`; type the PUBLIC methods on `Carrier`
  (`inner_product(x: Carrier, y: Carrier)`). A leaf `FunctionSpace[CompositeField]`
  subclass then has overrides that are EXACTLY the base at `Carrier=CompositeField`
  → LSP-clean by construction.
- Move the bare-array numpy BODY to a PRIVATE realization method typed `(x: Any)`
  (`_diagonal_inner_product`, etc.); the public surface forwards to it.

**Why the `Any` is HONEST, not a checking-erasure trick** (this is the review
axis to re-run on every such carve):
1. The `_diagonal_*` realizations have ZERO callers outside the surface that
   forwards to them — verify with grep; they are pure private hooks.
2. The bodies are verbatim moves of already-correct, test-covered numpy.
3. The alternative — typing the realization `NDArray` — FORCES a `cast()` in the
   surface (passing `Carrier` to an `NDArray` param does not type-check), and a
   cast is exactly what the campaign burns down. So `Any` on a documented private
   seam is MORE aligned with the campaign thesis than the cast, and strictly
   better than `# type: ignore` (anti-pattern #19). Erasure confined to the
   narrowest private surface, precise checked contract on the public methods.
4. Apply the split PRECISELY: a method that only COMPOSES the surface (`norm` =
   `sqrt(inner_product(x,x))`) needs NO realization and stays genuinely
   carrier-generic (and thereby becomes valid for composite carriers — a latent
   type-bug fix). Blanket-splitting every method would be the smell.

## `default=Any` vs `default=NDArray` — Any is REQUIRED, not lazy

`FunctionSpace` is INVARIANT in `Carrier` (it appears in both param and return of
`apply_metric: (Carrier)->Carrier`). Bare `FunctionSpace` slots (operator
`domain`/`codomain`) genuinely hold EITHER realization today. Under
`default=NDArray`, `FunctionSpace[CompositeField]` (the composite space) is NOT
assignable to a bare `FunctionSpace[NDArray]` slot → 300+ consumer sites break,
pulling the C3+ slot-specialization into C2. `default=Any` keeps both realizations
assignable AND is honest for today's architecture. The slot-specialization
(`FunctionSpace[Domain]`) is the declared C3+ follow-on (TypeVar doc-comment).
Mild CONCERN: until C3+ lands, `Carrier` buys precision only at the composite
override site (`Any` at the 300+ bare slots) — but the delivered win (numerics
5→0, clean composite overrides) is real and standalone. Staged deferral, not debt.

## Protocol home + underscore-in-public-Protocol

- A structural `Protocol` for the composite-field contract MUST live in `numerics`
  (layering forbids `numerics`→`transport` import); belongs beside its sole
  carrier class. Correct home.
- `__dataclass_fields__: ClassVar[dict[str,Any]]` on the leaf Protocol makes it
  satisfy typeshed `DataclassInstance` so `dataclasses.replace()` type-checks
  without a concrete import — precise, honest. NB: this is a Protocol, so the
  institutional-smell-#5 "stringized ClassVar becomes a dataclass field" trap does
  NOT fire (that trap is frozen-dataclass-only).
- `_recombine` (underscore) as a required member of a PUBLIC Protocol is a MILD
  smell (underscore says "private" but it's a load-bearing cross-layer contract
  member). Downgraded to NIT because it FAITHFULLY MIRRORS the pre-existing
  underscore method on the concrete carriers (`full_field.py`, `timed_full_field.py`)
  — the Protocol describes reality, it does not mint a new bad name. Renaming is
  out of the numerics-only cluster scope → file/defer, don't block.

## `Generic[Carrier]` on a frozen dataclass — no hazard (verified)

`Generic` contributes no annotation → adds no dataclass field; custom
`__eq__`/`__hash__` on `(name, shape)` untouched; `Carrier` erased at runtime so
identity is unchanged; frozen→frozen inheritance chain intact (`Generic` is a
non-dataclass base); PEP-696 `default=` runs natively on the Py-3.14 host with
`pythonVersion=3.13` pinned. 171 SN-adjoint + 1133 numerics tests green.

## THE APPROVAL CONDITION that recurs on EVERY ratchet-improving carve

`tests/test_pyright_ratchet.py` FAILS on IMPROVEMENT (by design — forces lock-in):
`pyright error count DECREASED (module: numerics: 5 -> 0)`. The carve is NOT
landable until `python -m tests._harness.pyright_ratchet --update` re-tightens the
baseline (numerics 5→0, total 115→110) and that baseline change is committed WITH
the code. Check `git status tests/_harness/pyright_baseline.json` — if unmodified,
the carve is incomplete and `main` goes RED. Always verify the baseline companion
landed; the diff's "ratchet 115→110" prose runs AHEAD of the working tree until it does.
