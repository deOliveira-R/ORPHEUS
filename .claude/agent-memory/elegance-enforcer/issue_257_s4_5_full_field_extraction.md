---
name: issue-257-s4-5-full-field-extraction
description: S4.5 FullField base extraction + TransportState Protocol retirement — PASS-WITH-NITS; the _recombine polymorphic-ctor-hook ruling + from_flat-via-hook + check_partner widening soundness
metadata:
  type: project
---

# #257 S4.5 — FullField timeless base extraction (the #217 cofree split)

PASS-WITH-NITS (`feature/field-typed-operator-algebra`, HEAD `93aa016`, UNCOMMITTED, behavioral-NEUTRAL;
51 file-local + 237 transport tests ✓; pyright 0/0/0 on both files; cross-method rejection probed LIVE).
Extracts timeless concrete `FullField` (`bulk ⊕ boundary`) base out of `TimedFullField`; retires the
S4 `TransportState` Protocol. `TimedFullField(FullField)` = `Cofree(FullField, depth=d)` adds only
`_history`/`history_depth`/`advance`/`at_lag` + the `_recombine` override.

## The load-bearing ruling: `_recombine` = the polymorphic-CONSTRUCTOR hook (Pattern 2 done right)

A frozen-dataclass base whose subclass differs in a CONSTRUCTION-time field (history) should route ALL
vector-space dunders + `from_flat` + `copy` through ONE `def _recombine(self: T, *, bulk, boundary) -> T`
hook. Base spells `replace(self, ...)` (provably T, nothing to drop); subclass OVERRIDES to rebuild the
subclass with the construction-field threaded + the iteration-metadata field RESET (`_history=()`,
`history_depth` preserved). This is the STANDING REMEDY for "base + subclass share an algebra but differ
in a ctor field" — algebra defined ONCE, correct concrete return type per subclass, Liskov-correct.
Generalises the emergent-gate-via-`replace` standing pattern from [[issue-257-s1-coefficient-fields]] J3
(there: re-run post_init; here: re-dispatch the ctor polymorphically).

## 4 judgment calls ALL PASS

- **(a) `_recombine` hook**: right abstraction, well-named, base-vs-override split clean (base can't know
  history; override owns exactly history). `copy` routing through it (drops history on a timed copy)
  MATCHES documented prior behavior — old `copy` body was `TimedFullField(_history=(),...)`, bit-id.
- **(b) `from_flat(cls, flat, template: T) -> T` via `_recombine`**: REMOVING the `TimedFullField.from_flat`
  override + its `# type:ignore[override]` is a real win — the override had NARROWED `template` (genuine
  Liskov violation the ignore papered over). Generic-template + route-through-the-same-hook is the
  elegant terminal state; does NOT over-couple (`from_flat` legitimately needs template's concrete type).
  Duck-typed Krylov adapter (`iteration.py:206` `type(template).from_flat(flat, template)`) consumes
  transparently. **STANDING TELL**: a subclass `.from_flat`/factory override that NARROWS the template
  param is a Liskov violation hiding behind `# type:ignore[override]` — collapse it to a base-generic
  routed through the polymorphic ctor hook.
- **(c) `zeros` delegation**: `base = FullField.zeros(...)` then re-wrap with `history_depth` keeps the
  duck-typed `zeros_on` `# type:ignore[attr-defined]` on the BASE ONLY (count stayed 2 not 4, verified
  both at `full_field.py:200-201`). Cannot route through `_recombine` (allocates from leaf TYPES not an
  instance) → explicit delegate-then-rewrap is the honest shape.
- **(d) `_check_partner` widening** `type(other) is not TimedFullField` → `isinstance(other, FullField)`:
  SOUND + necessary. Load-bearing for BDF stencil `at_lag(0) - at_lag(1)` (timed − timeless snapshot;
  `test_timed_full_field.py:385` ✓). Result type governed by `_recombine` (timed−timeless→timed).
  Opens NO cross-method hole — member-type rejection is SSOT on the LEAVES (probed: SN FullField −
  foreign ScalarFlux-bulk FullField rejects at `AngularFlux − ScalarFlux → TypeError`). Container
  deliberately does NOT pre-check member types → preserves the legit composite torsor flux+displacement.

## NOT-a-twin verified (the role-grid discriminator again)

`multiplication_operator.py:208,233` builds `TimedFullField(...)` INLINE not via `_recombine` — CORRECT,
not a DRY violation: those are OPERATOR OUTPUTS with a ROLE CHANGE (apply: flux→`AngularSourceSink`;
solve: source→`AngularFlux`). `_recombine` PRESERVES self's leaf types → would mint the WRONG codomain.
Operator `.apply`/`.solve` mint a new leaf type and CANNOT route through the type-preserving hook — same
"operator-output vs solve-trace" discriminator as [[issue-257-s3b-multiplication-operator]] / the role-grid
in AGENT.md §4. Out of scope for any container-algebra carve.

## Retirement clean + the nominal-check upgrade

`TransportState` Protocol fully gone (file deleted, zero live refs; survivors = 1 docstring history note
`timed_full_field.py:124` + memory slug + unrelated "transport state space" prose in operator/boundary).
Test migration structural→nominal (`isinstance` against the concrete `FullField` base) is a STRENGTHENING
— base is now a concrete class so the hierarchy check is genuine, and the suite ADDS a timeless/timed
distinction the old Protocol couldn't express. `-O`-firing `_require` + `_is_a` object-boundary helper
rigorous (vv Mode-8 aware).

## NITS (all do-now / follow-up, documentation+convention only, NONE blocking)

1. `timed_full_field.py:124` stale-by-omission: the only surviving `TransportState` token in prod; note
   the S4-era Protocol is now ALSO retired so a future grepper doesn't read it as live.
2. `full_field.py:341` function-local `import numpy as np` inside `to_flat` — carried forward verbatim
   but sibling leaf `angular_flux.py:48` imports numpy module-level; prefer module-level on a fresh write.
3. `full_field.py:319-327` `to_flat` 1-D-boundary assumption is SN-specific (`BoundaryFlux` flat storage);
   one-clause caveat future-proofs the flat protocol against a non-flat CP/MoC boundary leaf (Pattern-6
   defer applies — no 2nd method yet).
