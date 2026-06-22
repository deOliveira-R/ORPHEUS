---
name: issue-257-s8c-overload-fibration
description: S8c — PASS-WITH-NITS verdict ratifying "Pattern M" (private _apply_impl dispatcher + TYPE_CHECKING @overload alias) as the FIRST @overload precedent in orpheus/ for heteromorphic multi-carrier apply
metadata:
  type: project
---

# #257 S8c — honest @overload typing of the heteromorphic `apply` fibration

**Verdict: PASS-WITH-NITS** (`feature/field-typed-operator-algebra`, UNCOMMITTED working tree,
behavioral-NEUTRAL). Branch `refactor/pyright-signal` is where I reviewed from; the S8c edits sit
on the working tree of `feature/field-typed-operator-algebra`. Files: `orpheus/sn/fission.py`
(328-482), `orpheus/sn/scattering.py` (1128-1358), `tests/sn/operators/test_operators_apply_typed.py`
(C6 332-407).

## What S8c does
`FissionOperator.apply` / `ScatteringOperator.apply` are HETEROMORPHIC multi-carrier operators:
a `@singledispatchmethod` dispatches on input *carrier* → DISTINCT output carrier
(Scalar→ScalarSourceSink, Angular→AngularSourceSink, TimedFullField→FullField, ndarray→ndarray,
HarmonicMomentField→AngularSourceSink). They inherit `LinearOperatorMixin`'s nominal endomorphism
`apply(x:V)->V`, and the singledispatch base returned `NoReturn` → callers statically saw poison +
every register arm pyright-errored. The "fibration" (apply = a family of arrows indexed by carrier)
was a TYPING LIE. S8c makes it honest.

## ⭐ "Pattern M" — RATIFIED as the project's first @overload precedent
For each operator: rename dispatcher `apply`→`_apply_impl` (`-> "Any"` base so `.register` accepts
real-typed arms), keep all register arms at NATURAL class-body indentation, append:
```python
    if TYPE_CHECKING:
        @overload
        def apply(self, ...: <Carrier>, /) -> "<Output>": ...   # one per arm
        def apply(self, x: Any, /) -> Any: ...                   # overload-impl marker
    else:
        apply = _apply_impl   # public apply IS the same singledispatchmethod descriptor
```
**Pattern M chosen over Pattern C** (`if TYPE_CHECKING:<overloads> else:<ENTIRE dispatcher+~150-215
lines of source-assembly bodies>`). The DISPOSITIVE argument = the master standard's Beck ordering:
"reveals intention" (rule 2) OUTRANKS "fewest elements" (rule 4). Pattern C buries the real R·Λ·M
Galerkin / (n,2n) / 1/W math in a DOUBLE-INDENTED `else:` whose purpose is to hide typing scaffold —
inverts the hierarchy. Pattern M keeps the math where the eye expects it; typing lives in a labeled
epilogue. ⭐ KEY: Pattern M is NOT a foreign idiom — `LinearOperatorMixin` ALREADY uses
`if TYPE_CHECKING: def apply(self,x:V,/)->V: ...` (`operator.py:400-408`) as a type-only declaration;
S8c extends the SAME idiom (type-surface-in-TYPE_CHECKING, runtime-reality-elsewhere) with overloads.
The 3 costs are benign: `_apply_impl` private name (it IS the public apply), `else: apply=_apply_impl`
(descriptor alias = single-source, NOT a 2nd path), `-> Any` base (shadowed by overloads for all
callers; only the `object` fallback reaches it → raises TypeError fail-loud).

## STANDING TELLS for reviewing @overload-on-singledispatch carves
1. ⭐ VERIFICATION METHOD — stash-to-baseline pyright delta: `git stash push -- <the N files>`,
   `pyright --outputjson | count errors by (file,line)`, `git stash pop`. S8c: **33→17 errors, net
   −16, ZERO new**. The 16 eliminated were exactly the `@apply.register`-on-`NoReturn`-base errors.
   The 17 survivors are ALL pre-existing root gaps S8c unmasks (NOT in the new overload blocks).
   ALWAYS do the stash-baseline when a brief claims "errors all pre-existing" — line numbers SHIFT
   under a diff, so categorize SEMANTICALLY (read the surviving-error lines), not by line number.
2. ⭐ RUNTIME ALIAS PROOF — `Cls.__dict__['apply'] is Cls.__dict__['_apply_impl']` is the right check
   (returns True = same singledispatchmethod descriptor). Do NOT use `Cls.apply is Cls._apply_impl`
   (returns False — `descriptor.__get__` makes fresh bound wrappers per access; a test artifact, not
   a code defect). Confirm registry parity via `desc.dispatcher.registry` (the `object` key = the
   base/fallback arm).
3. ⭐ STATIC GATE PROOF — the pyright-only uncalled `_c6_static_typing_pins(...)` of `assert_type`
   pins is CLEAN (no `test_` prefix → pytest never collects; only the type checker reads it). The
   gate passes iff ZERO pyright errors land INSIDE that function's line band. This is the right idiom
   for a static-typing gate; pair it with a runtime `pytest.fail`-not-`assert` parity test (Mode-8).
4. CORRECTED DISPATCH-TABLE DOCSTRINGS must be verified against the LIVE registry, not trusted. S8c's
   docstrings now match exactly (fission: TimedFullField/ScalarFlux/ndarray — NO AngularFlux arm;
   scattering: TimedFullField/ScalarFlux/AngularFlux/HarmonicMomentField — NO ndarray arm). A stale
   dispatch table in a docstring IS the "comment asserts a false invariant" bug-habitat tell.

## cast() bridges — HONEST not laundering
`scattering.py:1237` `self.apply(cast("AngularFlux | HarmonicMomentField", psi.bulk))` + test double-
cast `cast(ScalarFlux, cast(AngularFlux, state.bulk).integrate_angular())`. Honest because: cast
target is the EXACT runtime invariant AND both carriers → `AngularSourceSink` (output correct either
way); a 3rd carrier hits the `object` fallback → TypeError fail-loud. They bridge a PRE-EXISTING root
gap S8c merely unmasks — DEFER-not-fix was the right scope (don't balloon the carve into
transport/fields/).

## ⭐ ARCHITECTURAL OPPORTUNITY (issue recommended under #256)
Root under-typing: `FullField.bulk: BulkField` (`full_field.py:164`, over-broad base) +
`AngularFlux.integrate_angular(self) -> "object"` (`angular_flux.py:116` — Pattern-3/13 violation: the
docstring KNOWS it returns ScalarFlux, a named 1/(cm²·s) quantity, yet annotates `object`). These force
the casts. Fix at root retires them. File under #256 ("close the ndarray architecture gaps").

## NITS (non-blocking)
- N1 (do-now, precedent): document "Pattern M" as a named idiom (one para in operator-algebra doc /
  coding-elegance skill) BEFORE a 3rd multi-carrier operator copies the block by eye — #261 will
  relocate these cores → a 3rd copier is imminent. Pattern-6 "document-the-idiom-after-two".
- N2 (optional): C6 static/runtime asymmetry — `_c6_static_typing_pins` pins `S.apply(moments)` but the
  runtime `cases` list omits the HarmonicMomentField arm (its runtime dispatch IS covered by Phase-5a
  windowing + test_scattering_operator.py, so NOT a coverage hole — just a self-containedness gap).

Gates: 17✓ typed-apply (incl C6) + 85✓ fission/scattering value-level.
