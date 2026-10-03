# #405 step 2, phase P2: the solutions and the certificates (verification specification)

Written by the test-architect on 2026-10-03 for workflow W3 (surgical carve), phase P2, against `main` `317418da`. The design it specifies is "P2, the design as ruled (2026-10-02)" at the end of `.claude/plans/reference_cache.md`, with the census, the rulings G1 to G6 and the error-ontology draft above it. The form follows `.claude/plans/reference_p1_spec.md` §1.7 and §1.8: per step, a table of gates with their kind, their first red on the tree as it stands, the mutation that must redden each, and the rungs each rests on.

Gate ids here are `R<step>.<n>` (R for the reference phase); an `S<step>.<n>` id is P1's, in `reference_p1_spec.md`, and the test functions are named `test_r<step>_<n>_…`. Step 1's gates are written, run and mutated (§1.1). The other steps are specified to the level the main agent needs to write the code; each is re-specified against its prototype when it opens, as P1's steps were.

## 0. Rulings this specification rests on (2026-10-03, beyond the design of 2026-10-02)

The orchestrator answered four questions the design left open; the user ruled the second.

1. **`Printed(text, citation)`** (Q1, ruled as proposed). The printed decimal string is the one source: it is canonicalised by `str(Decimal(text))`, so `"1.0E0"`, `"+1.0"` and `"10E-1"` are one claim while `"1.0"` and `"1.00"` are two; the value (the correctly rounded double) and the half unit in the last printed digit are derived from it; the citation must carry a locator, because a value is printed at a place. The design's spelling `Printed(value, digits, citation)` is superseded: "digits" was ambiguous between significant figures and decimal places, a value and a digit count could disagree, and a float drops the trailing zeros that ARE the precision claim.
2. **The observables** (Q2, ruled by the user). The primitive is `FluxIntegral(weight)`, a linear functional of the scalar flux with the weight a mesh-free function over position and group. `Rate(cells, weight)` is not a member of the closed sum: it is a spelling the SPECIFICATION resolves into a `FluxIntegral` (the weight times the cells' coefficient field from the materials), as question keys resolve to explicit cells, so one functional has one canonical spelling after resolution. The closed sum is `Observable = FluxIntegral | Ratio | Eigenvalue | PointValue`.
3. **The algebraic-error floor** (Q3, ruled with a distinction). The `VerificationCertificate` takes production's algebraic-error `Evidence` as an explicit argument with no default, and records it. For an AGREEMENT verdict (`|v_prod − v_ref| + b_ref ≤ tol`) the comparison bounds the total error whatever its split, so `NotYet` or `NotApplicable` is recorded, not refused. For an ORDER verdict (a ladder attributing error to the discretisation) the algebraic error must be `Measured` or `Certified` at or below a tenth of the tolerance, or the verdict refuses as unestablished. No member of today's `ExitCertificate` is a stage-A error estimate of a functional, so none is mapped; the missing estimator is #564, and every fixture that has no estimate passes `NotYet(564, …)`.
4. **The homes** (Q4, ruled as proposed): `Enclosure` in `orpheus/numerics/enclosure.py`; everything that names a specification or a citation in a new input-tier package `orpheus/reference/` (it may import `specification`, `numerics`, `data`, `geometry`; it may not import `mesh`, `transport`, an L3 package or `derivations`; `derivations`, which writes reference solutions, may import it), registered in `tests/gates/test_layer_imports.py`.

## 1. The step order, the sizes, and the gates

Sizes are sizing only, never a reason (`plan-authoring` EFFORT-IS-SIZING). "Rows" counts collected test rows; "prod" counts production lines, `[R]` estimates.

| step | content | size | depends on |
|---|---|---|---|
| 1 | The readings: `Enclosure` and its quotient; `Printed`; `ReferenceReading` and `ProductionReading`; the `orpheus/reference/` package and its layer registration | small: about 130 prod lines, 2 new test files (90 rows), two one-line amendments (S5.9's union, the linter). **Written and red, §1.1** | — |
| 2 | `ExitCertificate` renamed `ExitReport`; the `Solution` slot `certificate` renamed `exit_report` | small, mechanical: 16 lines in 6 files name the class (`orpheus`, `tests`), 36 lines in 9 files read or pass the slot, 21 mentions on 8 theory pages and `CLAUDE.md` | — |
| 3 | The observables: `FluxIntegral`, `Ratio`, `Eigenvalue`, `PointValue`, the closed sum; content identity | small: about 120 prod lines, 2 test files | 1 |
| 4 | The answer representation: graded piecewise-Chebyshev panels on region-local charts | medium-large: about 350 prod lines, the most mathematics in P2 | — (parallel to 2, 3) |
| 5 | `PublishedSolution` (with the `Withdrawn` state) and `ReferenceCertificate` | medium | 1, 3 |
| 6 | `ReferenceSolution` and `read(observable)`; the first instance, the exact infinite medium | medium | 3, 4, 5 |
| 7a | `VerificationCertificate` and the comparison verb (agreement and order verdicts, the floors) | medium | 6 |
| 7b | The migration of the four `certify_agreement` files; `AgreementCertificate` retires | medium-large (it needs the trajectory-resolvent sphere and cylinder written as reference solutions, and the shape rows re-posed, §1.7b) | 7a |
| 8 | `Rate(cells, weight)`, the specification's resolving constructor | small-medium; RECOMMENDED to land with posing unit 3 (#526), §1.8 | 3 |

The order puts the laws that everything else reads (the readings, the observables) first and the migration last. Steps 2 and 4 are independent of 1 and 3 and can be taken in any position.

### 1.1 Step 1: the readings (gates written, run and mutated, `[M]` 2026-10-03)

**Artefacts** (untracked, `scratch/reference_architecture/p2/ta/`): the draft gate files `drafts/tests/gates/numerics/test_enclosure.py` and `drafts/tests/gates/reference/test_readings.py` (+ `__init__.py`), which land at the same paths; `s1_union_and_layer.diff`, the amendment of `tests/gates/numerics/test_content_identity.py` (the two rosters join S5.9's union, and `orpheus.reference` joins S5.9's package walk) and of `tests/gates/test_layer_imports.py` (the package registered; two cold entry points); `proto/`, a GATE-HARNESS PROTOTYPE written only so that the gates could be run, type-checked and mutated (it is not a design and not production code); `battery/` (the `-p` plugin `s1_battery.py`, its driver, `battery_summary.txt`, one log per arm); `run_absent.log`, `run_proto.log`, `run_proto_amended.log`; `probe_tight.py`. Every measurement was taken in a detached worktree of `317418da` with `.venv` linked, under `.venv/bin/python -O -m pytest … -p no:cacheprovider --color=no`, with a session-finish plugin printing `orpheus.__file__` (it resolves into the worktree).

**Measured first reds:**
- On `317418da` as it is, both files fail at collection: `ModuleNotFoundError: No module named 'orpheus.numerics.enclosure'` and `… 'orpheus.reference'` (`run_absent.log`). Every row is DEFINING; its tooth is its battery arm below.
- With the prototype present and the two amendments NOT applied: 89 of 90 rows pass and `test_r1_18_the_package_is_registered_with_the_linter` is red, "the linter does not know the package: its imports are never checked" (`run_proto.log`). This red is real on the landed tree, not a placeholder: `_check_source` returns `[]` for a package missing from `FORBIDDEN_EDGES`, so an unregistered package is never linted and its zero violations mean nothing.
- With the prototype present, `tests/gates/numerics/test_content_identity.py::test_s5_9_the_population_is_the_roster` is red, "content types with no roster entry: ['Enclosure']". It does NOT name `Printed`: S5.9 walks `_PACKAGES`, which omits the new package, so `Printed` would escape the roster silently. The amendment adds `"orpheus.reference"` to `_PACKAGES` as well as the two rosters.
- With the prototype and both amendments: the two new files, `test_content_identity.py` and `test_layer_imports.py` read 568 passed (`run_proto_amended.log`). `npx pyright` on the two drafts and the prototype: 0 errors (4 warnings: the unused expressions inside `pytest.raises` blocks).

**The gates** (every class is resolved on its module at run time, so a rebinding arm reaches every row, lessons `L100`):

| id | gate | kind | battery arm (tooth) | rests on |
|---|---|---|---|---|
| R1.1 | **Admission.** `Enclosure(value, bound)`: NaN or infinite value or bound `ValueError`, a negative bound `ValueError` ("the bound is a distance"), a complex, string, `None` or `bool` `TypeError`, each message naming the field (9 rows); positive legs: an exact value, integers stored as doubles, a negative value, `-0.0` bound stored `+0.0` | THEOREM (defining refusal) | A7, the negative-bound check removed → `[negative-bound]` red | S5.4 (the encoder), `parse_finite_real` |
| R1.2 | **Containment.** For every population row, every exact quotient `x/y` of members lies in the quotient's interval, checked in `Fraction` at the four corners (the extremes of `x/y` on a box excluding `y = 0`). Population: 10 hand rows (exact operands whose quotient is not a double, every sign pair with wide enclosures, a numerator straddling zero, a denominator 1e-15 from touching zero, `1e-300 / 1e300`) and 3000 draws, seed 20261003, magnitudes 1e-30 to 1e31, relative widths 0 to 0.99 (numerator) and 0 to 0.9 (denominator); `[M]` 3010 of 3010 form a quotient, 0 fail. Plus the named row: `Enclosure(1,0)/Enclosure(3,0)` has a positive bound and contains 1/3 | THEOREM | A1 (first-order bound) → both red; A2 (round to nearest, no outward step) → both red; A3 (two corners, right only for positive operands) → the population row red | R1.1 |
| R1.3 | **The centre** is `a.value / b.value` bit for bit, so a ratio observable reads the ratio of its readings | THEOREM (the design's convention) | A6, the hull's midpoint as centre → red | R1.2 |
| R1.4 | **Tightness.** The bound exceeds the exact half-width about the centre by at most 10.5 ulp of the largest exact corner. Derived: each operand endpoint is rounded once and stepped outward once (3u relative), a corner quotient is perturbed by 6u, rounded (u) and stepped (2u), so 9 ulp; the half-width's subtraction is rounded and stepped, 1.5 ulp more. `[M]` worst 7.84 over the population (`probe_tight.py`). A merely valid bound (1e300 encloses everything) would fail every verification floor, so the law has two sides | THEOREM (the outward-rounding contract) | A5, the bound doubled → red (also A6) | R1.3 |
| R1.5 | **The quotient's refusals.** A denominator interval holding zero (exactly zero, touching, straddling, straddling negative: 4 rows) raises `ZeroDivisionError` matching "contains zero"; `1e300 / 1e-300` raises `ValueError` matching "infinite" (never an enclosure with an infinite bound); a bare `float` or `int` divides neither way (`TypeError`), so an exact value is spelled `Enclosure(v, 0)` once (3 rows) | THEOREM (defining refusal) | A4, the zero check removed → 4 of 4 red; A8, a float accepted as an exact operand → `[float]`, `[int]` red (`[none]` stays green: `None` is no number under either) | R1.2 |
| R1.6 | **The fields** are exactly `(value, bound)`: the interval's ends are derived, never stored beside the centre (X4) | THEOREM (structure) | B5's sibling not run for `Enclosure` (owed post-landing) | — |
| R1.7 | **Content identity of `Enclosure`**: the roster entry (parts `value`, `bound`); each part moved by one ULP or a sign moves the digest, `==` and a set (4 legs); equal content is one value (int and float, `-0.0` bound, `-0.0` value: 3 pairs); the pickle round trip | THEOREM | P0, `content_parts` returns `()` → the population row and 4 of 4 legs red | S5.2, S5.3 (`test_content_identity.py`) |
| R1.8 | **The layer of `enclosure.py`.** (a) a fresh interpreter importing it and forming one quotient loads only the `orpheus` packages `{geometry, numerics}` (what a cold `import orpheus.numerics` loads, `[M]`) and no SymPy; (b) by AST, every `orpheus` import is under `orpheus.numerics`, activation: the import of `numerics.content` is seen | THEOREM (the layer contract) | C3 (the file on disk gains `import orpheus.data.citation`) → red; C6 (the script also imports `orpheus.transport`) → red | S6.21, S7.13 |
| R1.10 | **One claim, one value.** 7 pairs that print one claim (`1.0E0`/`1.0`, `+1.0`/`1.0`, `10E-1`/`1.0`, `1.234E-03`/`0.001234`, lowercase `e`, `-0.000`/`0.000`, surrounding spaces) are equal both ways, one hash, one set member, one digest; 4 pairs that are two claims (`1.0`/`1.00`, `10`/`1E+1`, a sign, a last digit) are two values | THEOREM | B1, the raw text stored → the 6 non-space pairs and the roster pair red (the arm strips spaces itself, so `[surrounding-space]` is not its target); P0 → the 4 two-claim rows red | R1.15 |
| R1.11 | **The value** is a `float` equal to `float(text)` (correctly rounded), 9 rows | THEOREM | not isolated by an arm (`value` is one call of the trusted parser) | — |
| R1.12 | **The half unit** is a `Decimal` equal, exactly, to half a unit in the last printed digit, against a HAND-WRITTEN table (9 rows, from `1.00000` → `0.000005` to `6.02214076E+23` → `5E+14`): the table is the independent route to the exponent arithmetic | REFERENCE (hand table) | B2, the half unit one decade off → 9 of 9 red | — |
| R1.13 | **Refusals of the text**: 11 non-decimal or non-finite texts (`""`, `1,0`, `abc`, `1.0.0`, `0x1p3`, `nan`, `NaN`, `sNaN`, `inf`, `-Infinity`, `1e400`) `ValueError` with three disjoint fragments ("not a printed decimal number", "not a finite printed number", "beyond the range of a double"); a `float`, a `Decimal`, `None` `TypeError` ("a printed value is its decimal text") | THEOREM (defining refusal) | B4, NaN admitted → `[nan]`, `[NaN]`, `[sNaN]` red | — |
| R1.14 | **The citation names a place**: a `Citation` with no locator `ValueError` ("locator"), a bare key `TypeError`; the positive leg keeps the citation | THEOREM (defining refusal) | B3, the locator check removed → red | `Citation`'s own gate |
| R1.15 | **One source.** `dataclasses.fields(Printed)` is exactly `(text, citation)`: no stored `value`, `digits` or `half_unit`. Content identity: roster (parts `text`, `citation`), 4 legs (a last digit, one more digit, another locator, another work), 1 pair, pickle | THEOREM | B5, a stored `digits` field → the fields row, the population row and R1.16 red; P0 → population and 4 of 4 legs red | S5.2, S5.3 |
| R1.16 | **The two sums are disjoint by role.** `ReferenceReading` names exactly `{Enclosure, Printed}`; `ProductionReading` names exactly `{Measured}`, and that class IS `orpheus.numerics.outcome.Measured`; no member of one is a subclass of a member of the other; an `Enclosure` and a `Printed` are not production readings, a `Measured` not a reference reading. An `Estimated` reference is therefore unspellable twice over: `Estimated` would join the production sum only | THEOREM (structure) | C1, `Measured` added to `ReferenceReading` → red | R1.6 |
| R1.17 | **`Estimated` is not built; `Measured` is defined once.** By AST over every `.py` under `orpheus/` (`[M]` 376 files parsed in the prototype worktree, which adds 3; more than 300 asserted): no class `Estimated` (the user's G3 ruling: it waits for a Monte Carlo consumer, and this row is edited ON PURPOSE when it lands); `Measured` defined exactly in `orpheus/numerics/outcome.py` (the positive control) | THEOREM (census) | C2, a class `Estimated` appended to `outcome.py` on disk → red | — |
| R1.18 | **The package's layer.** (a) The linter knows it: `FORBIDDEN_EDGES["reference"]` holds `mesh`, `transport`, `derivations` and the L3 packages; `numerics`, `geometry`, `data`, `mesh` and `specification` are forbidden to import it; `derivations` is not; `orpheus.reference` is a cold entry point. (b) A fresh interpreter importing `reference.reading` and building one `Printed` loads only `{reference, specification, numerics, data, geometry}` (`[M]` data, geometry, numerics, reference) | THEOREM (the layer contract) | (a) **red today** (above); C4, the registration popped → red. (b) C6 → red | the linter's own parametrised row |
| R1.19 | **Seed-stable digests.** `PYTHONHASHSEED` 1 and 2 print the same digests and hashes of an `Enclosure`, a quotient and a `Printed` (a `str` part is where a salted hash would leak); the control line, a `str` hash, differs | THEOREM | C5, a salted `hash` mixed into the subprocess's digest → red | S5.1 |

**The battery** (`battery/s1_battery.py`, selected by `S1_ARM`, installed at `pytest_configure` before collection; each arm checks that its mutant differs from the honest object on a probe taken before the patch and raises `Uninstallable` otherwise; textual arms copy the file aside and restore it, and the driver `diff -q`s after every arm: 0 restore failures). Scope: the two step-1 files, 90 rows. `[M]` 2026-10-03:

| arm | mutation | red set (target rows in bold) |
|---|---|---|
| none | — | 0 |
| **P0 (positive control)** | `content_parts` returns `()` for both types | 14: **2 population rows, 8 perturbation legs**, the 4 two-claim rows |
| A1 | first-order bound `(a.b·|b.v| + b.b·|a.v|)/b.v²` | 5: **R1.2 both**, R1.3, R1.4, the overflow row |
| A2 | round to nearest, no outward step | 2: **R1.2 both** |
| A3 | the two corners right only for positive operands | 1: **R1.2 population** |
| A4 | the zero check removed | 4: **R1.5 4 of 4** |
| A5 | the bound doubled | 1: **R1.4** |
| A6 | the hull's midpoint as centre | 2: **R1.3**, R1.4 |
| A7 | a negative bound admitted | 1: **R1.1 `[negative-bound]`** |
| A8 | a float divides as an exact enclosure | 2: **R1.5 `[float]`, `[int]`** |
| B1 | the raw text stored | 7: **6 one-claim pairs, the roster pair**. First run 0 red: the constructor's eager digest had read the canonical text and is cached by `id` (lessons `L101`); the arm drops the cache entry |
| B2 | half unit ×10 | 9: **R1.12 9 of 9** |
| B3 | locator check removed | 1: **R1.14** |
| B4 | NaN admitted | 3: **R1.13 `[nan]`, `[NaN]`, `[sNaN]`** |
| B5 | a stored `digits` field | 3: **R1.15 fields, population**, R1.16 |
| C1 | `Measured` in `ReferenceReading` | 1: **R1.16** |
| C2 | class `Estimated` minted | 1: **R1.17** |
| C3 | `enclosure.py` imports `orpheus.data` | 1: **R1.8** |
| C4 | the linter registration popped | 1: **R1.18 (a)** |
| C5 | salted digest in the subprocess | 1: **R1.19** |
| C6 | the cold scripts import `orpheus.transport` | 2: **R1.8, R1.18 (b)** |

21 arms, 0 blind after B1's repair. Owed after landing: the arms re-targeted on the production modules (the arithmetic arms as `inspect.getsource` transforms of the real `__truediv__`), and an R1.6 arm.

**What no step-1 gate can see.** Whether a production author's printed text is transcribed correctly from the paper: that is the `PublishedSolution` registry's review, and the per-value citation (R1.14) is what makes it checkable. Whether an enclosure a reference REPORTS is true: that is the `ReferenceCertificate`'s (step 5); `Enclosure` checks only its own arithmetic.

### 1.2 Step 2: `ExitCertificate` becomes `ExitReport`

The user's G4 ruling: two certificates (`ReferenceCertificate`, `VerificationCertificate`); production's exit self-report is not one. Proposed names: the class `ExitReport`, the `Solution` slot `exit_report` (the slot keeps no "certificate"; `record` already names the `IterationRecord`, so `report` alone would invite a confusion), the builder `_exit_report`. A pure rename: every gate that survives is classed by `retirement-audit` D (none is expected to DEMOTE or PROMOTE, since no comparison changes sides).

| id | gate | first red | tooth |
|---|---|---|---|
| R2.1 | AST census over `orpheus/` and `tests/` (input count printed): no identifier `ExitCertificate`, no attribute `.certificate` on a `Solution`-typed receiver, no keyword `certificate=` into `SolutionBase`; positive control: the census finds `ExitReport` in `numerics/outcome.py` | red on `317418da`: `[M]` 16 lines in 6 files name the class, 36 lines in 9 files read or pass the slot (`git grep -P '\.certificate\b|certificate='`) | re-insert one old spelling → red |
| R2.2 | The five-tuple's field names are `(problem, outcome, strategy, exit_report, record)`, read by `dataclasses.fields(SolutionBase)` | red today | — |
| R2.3 | The rename's prose: a grep of `ExitCertificate` over `docs/` (excluding `_build`), `CLAUDE.md` and `.claude/rules/` returns only past-tense history lines (`retirement-audit` E.19), and the `.rst` underline scan after a length-changing rename (F.23) is clean; `dead_references` after | 21 mentions on 8 pages + `CLAUDE.md` ("the frozen five-tuple problem, outcome, strategy, certificate, record") | — |

`ConvergenceCertificateError` (an exception: a within-group solve claimed convergence and the honest residual disagrees; 13 + 9 occurrences) also carries the word; whether it is renamed is a NEEDS item.

### 1.3 Step 3: the observables

`orpheus/numerics/observable.py`: `FluxIntegral(weight: MeshFreeFunction)`, `Ratio(numerator: Observable, denominator: Observable)`, `Eigenvalue()` (field-less, like `Fundamental`), `PointValue(position: float, group: int)`; `Observable` their closed sum. Physics-free, content identity, admitted eagerly. A weight's fit to a problem (its group count, its region count, the coordinates a `Symbolic` may depend on) is the specification's `_admit_datum`, reused at read time (one definition, never a second admission here).

| id | gate | kind | first red | tooth | rests on |
|---|---|---|---|---|---|
| R3.1 | **Closed.** `get_args(Observable) == {FluxIntegral, Ratio, Eigenvalue, PointValue}`; the content classes of the module are exactly those four; a test-side `match` with `assert_never` names each once. `Rate` is NOT a member and NOT defined in numerics (the user's Q2 ruling) | THEOREM (structure) | defining | a fifth content class minted in the module → red; `Rate` added to the sum → red | R1.16 |
| R3.2 | **Fields are the roles** (by name, resolved at run time): `FluxIntegral (weight)`, `Ratio (numerator, denominator)`, `Eigenvalue ()`, `PointValue (position, group)`; no `bool` field, no default on any datum | THEOREM | defining | a default weight → red | R3.1 |
| R3.3 | **Admission.** A weight that is not a `MeshFreeFunction` (an `ndarray`, a `float`, `None`) `TypeError` naming the weight; a `Ratio` operand that is not an `Observable` `TypeError`; a non-finite position `ValueError`; a negative or non-integer group `ValueError`/`TypeError` | THEOREM (defining refusal) | defining | each check removed → its rows red | S6.8, R1.1 |
| R3.4 | **Content identity**: roster rows for the four types (every part moved by one leg, equal-content pairs, pickle), the union amended; `Ratio(a, b) != Ratio(b, a)`; `Ratio` nests (`Ratio(Ratio(a, b), c)` is a value) | THEOREM | defining; S5.9 reds until the union is amended (as step 1) | P0 | S5.2, S5.3 |
| R3.5 | **Seed-stable digests and the RECORD bytes** of four canonical observables, pinned after landing ("pin after landing" until then, as S7.12) | THEOREM, RECORD | defining; RECORD red by design until pinned | a salted `str` part; a `__content_version__` bump | S5.1, S5.7 |
| R3.6 | **The layer**: as R1.8, `{geometry, numerics}`, no SymPy on import or on building a `FluxIntegral` of a `RegionwiseConstant` | THEOREM | defining | the module importing `orpheus.data` → red | R1.8 |

**Step 3's battery** (`[M]` 2026-10-03 on the production module at `9f8d1f7f`, before the review round of `cf9893e8` re-posed R3.3–R3.5). Artefacts: `scratch/reference_architecture/p2/ta/step3/`. Scope: 61 rows, baseline 1 red (the RECORD placeholder). Result: 14 arms, every target red, 0 blind. One gate was tightened: the weight-refusal rows matched only the bare `weight`, and the encoder's own refusal of a writeable array names the path `FluxIntegral.weight`. With the type check removed, the ndarray and list rows therefore stayed green. The rows now match `the weight is a mesh-free function`.

### 1.4 Step 4: the answer representation (re-specified 2026-10-03 against the investigator's prototype; gates written, run and mutated)

**WITHDRAWN from P2 (the user's ruling 2 of 2026-10-03, `reference_cache.md`, "The reading-bound prototype and the user's rulings").** The text, drafts and battery below are kept as P4's starting material, where graded panels return as each reference method's discretisation of its emission density (#566).

The 2026-10-03 first text of this section is superseded by this one; the R4 ids below are new.

**The value.** `orpheus/numerics/piecewise_chebyshev.py`, `ChebyshevPanels(breakpoints, coefficients, bounds)`, holds ONE region's answer on its chart `u ∈ [0, 1]`:
- the panel breakpoints, increasing strictly from 0 to 1;
- one `(terms, groups)` Chebyshev array per panel;
- one bound per panel and group, the claim of whoever built the value.

The bound is stored the way `Enclosure` stores its bound. The design had said "the bound is derived, the trailing zeros trimmed", and the test-architect found those two in conflict. The resolution rule needs the stored tail, while trimming removes exact zeros, so an exact polynomial stored trimmed would be refused. Two equal values (with and without appended zeros) would then be one accepted and one refused. Storing the bound removes the conflict:
- the rule lives in `fit`, which judges a SAMPLED function;
- a value built directly (an exact polynomial, or a reload) carries its builder's bound;
- trimming exact zeros is then harmless.

The methods:
- `fit(f, singular_ends, layers, slope, tau)` grades panels geometrically (ratio 0.15) toward the declared singular ENDS of the chart, at degree `max(2, round(slope·(j+1))) + 1` for the j-th panel from the singular end. These are the investigator's grading and degree law.
- `fit` samples `f` at the first-kind Chebyshev points. It refuses the whole region with `UnresolvedPanel(panel, interval)` unless, in every group, the coefficients fall below `tau` times the group's scale (the largest sampled magnitude) and stay there for the last three.
- `fit`'s bound per panel is twice the sum of that resolved tail.
- `evaluate(u, group)` reads a breakpoint from the panel on its right and `u = 1` from the last panel. It refuses a point outside `[0, 1]`, NaN included.
- `integral(interval, cell, power)` returns, per group, the integral of the field times `x**power` over a physical cell of the physical interval `[a, b]`. Powers 0, 1 and 2 serve the slab, the cylinder and the sphere; the `2π` and `4π` are the reader's.
- `to_arrays()` and `from_arrays()` give plain arrays for `.npz`.

A singular point is always a region edge (an interface or a non-mirror boundary), so declaring the singular set is declaring ends; the centre of a sphere is never one.

**Artefacts** (untracked, `scratch/reference_architecture/p2/ta/step4/`):
- `drafts/test_piecewise_chebyshev.py` and `drafts/_panel_fixtures.py`, which land in `tests/gates/numerics/`. The fixtures are F1, F1b, F2, F3 and F4, carried from the investigator's `common.py`; F5 needs the investigator's Nyström solver and is owed (§NEEDS).
- `proto/piecewise_chebyshev.py`, a GATE-HARNESS PROTOTYPE (the investigator's candidate (b) as a value; not a design).
- `battery/`, the logs, and the probes `probe_s4.py` and `probe_s4_slack.py`.
- Every production name is reached through the draft's adapter section, so a rename is one edit there.

**Measured first reds** (`[M]` 2026-10-03, worktree of `9f8d1f7f`, `-O`, `orpheus.__file__` verified):
- **Without the module:** collection fails with `ModuleNotFoundError: No module named 'orpheus.numerics.piecewise_chebyshev'`.
- **With the prototype:** 79 passed, the slow sphere row included (8.9 s). The first run had 1 red, R4.4 `[nan]`: the prototype's chart check `(u < 0) | (u > 1)` admits NaN. The row stays, as a requirement on the module.
- **S5.9's union:** S5.9 reds on `ChebyshevPanels` until the roster in the draft joins it.
- **pyright:** 0 errors.

| id | gate | kind | battery arm (tooth) | rests on |
|---|---|---|---|---|
| R4.1 | **An exact polynomial's cell read is its exact integral, under every measure, in group order.** The value is built directly with dyadic breakpoints and coefficients and bounds 0. The fixture has 3 panels of 7, 4 and 9 terms and 2 groups with different series. The parameter grid is 6 cells (the whole interval; on breakpoints; inside one panel; spanning three; the last panel to the edge; an empty cell) × 2 intervals (`[0, 1]` and the shell `[1, 2]`) × 3 powers. The reference is the exact rational integral (SymPy `Rational`). The tolerance is the derived rounding bound `γ_m · Σ_pieces length · Σ|c_j| · max|x|^power`, with `m = 4n + N + power + 8` (Clenshaw, Gauss sum, measure products, affine maps; Higham eq. 4.4 and Lemma 3.1). `[M]` the worst error is 0.10 of it over 72 group-rows (59 nonzero). Plus a cell outside the interval is refused | THEOREM (exact) | B1, Gauss one degree low → 30 of 36 rows red; B3, groups reversed → 30 red; B12, the cell check removed → the refusal row red | — |
| R4.2 | **The measure is activated.** These are R4.1's power-1 and power-2 rows. A slab-only fixture never reads the measure, so these rows are MANDATORY | THEOREM | B2, the measure dropped → 24 red, 0 of the 12 slab rows (the declared partial null holds) | R4.1 |
| R4.3 | **Groups**: R4.1's two-group fixture. A 1-group fixture cannot see a transposed group axis (lessons `L86d`) | THEOREM | B3 → 30 red | R4.1 |
| R4.4 | **Point evaluation** equals `chebval` of the panel's own series at 50 interior points per panel, per group, within 4 ulp of the sum of the series' magnitudes. A breakpoint is read from the right panel (activated: the two panels differ there); `u = 1` from the last. Points below 0, above 1 or NaN are refused | THEOREM | B4, a breakpoint read from the left → 2 red; B5, the local map reversed → 2 red (plus the 5 R4.6 rows); B11, the chart check removed → 3 red | — |
| R4.5 | **An undeclared singular point is refused, at the panel holding it.** F2, at 14 layers and slope 3. Declared (two regions split at the interface x = 1): accepted. Undeclared (one region `[0, 3]`), and misplaced by 1e-3 (`[0, 1.001]`, `[1.001, 3]`): `UnresolvedPanel`, and the refused panel's interval holds x = 1. `[M]` the prototype refuses panel 14, `u ∈ [0.075, 0.5]` (x ∈ [0.225, 1.5]), undeclared, and panel 18 next to 1.001, misplaced. The investigator's probe7 agrees: 1 of 30 and 1 of 60 panels refused, 0 of 60 declared | THEOREM (defining refusal) | B10, the grading ignored → the declared row is refused (red), plus 5 rows; B6, one trailing coefficient instead of three → the undeclared rows stay refused (the arm is aimed at R4.8) | — |
| R4.6 | **The declared bound holds on every panel of the closed-form fields.** Each panel's true maximum error, over 200 uniform points plus 240 log-clustered toward its ends, is at most its bound. Fields: F1 (2 mfp), F1b (40 mfp), F2 (both sides of the interface), F4 groups 0 and 1 (smooth diffusion), and F3 the sphere (`slow`; its centre not declared singular, its interface and surface declared). `[M]` on the prototype the bounds are about 1e-8 and the true errors about 1e-14 (bound, not estimate). The investigator measured true/bound 0.011–0.19 on the same rule | REFERENCE (closed forms) | B8, the bound taken as the last coefficient → 5 of 5 fast field rows red; B5 → 5 red | R4.5 |
| R4.7 | **A read is within the bounds times the measure.** For the smooth field `e^{-x} cos 3x` (fit at tau 1e-12), each read of the whole interval and of a partial cell, on `[0, 1]` and `[1, 2.5]`, under powers 0, 1 and 2, differs from mpmath (30 digits) by at most `Σ_k b_k ∫_piece |x|^power dx`, plus 1e-13 of the integral of the absolute field (the rounding allowance, `[R]`). This is the bound the reference solution's `FluxIntegral` reading reports (step 6). MANDATORY curvilinear rows | REFERENCE (mpmath) | not separately armed (B2 is R4.1's) — owed | R4.1 |
| R4.8 | **The rule, exactly.** One unsingular panel of 8 terms with series `[1, .5, .25, .1, 1e-3, 1e-12 ×3]` is resolved, with bound `2·3e-12` (to 6e-15). The same series with `c_5 = 1e-3` (two trailing) is refused | THEOREM | B6, one trailing → the refusal row red (and R4.9); B8 → the bound row red | — |
| R4.9 | **The rule is relative.** `fit(2^40 f)` has bounds exactly `2^40` times f's. A 1e-20 field resolves with bound ≤ 1e-30, and its two-trailing twin is refused | THEOREM | B7, an absolute tau → red | R4.8 |
| R4.10 | **Content identity and storage.** Population (breakpoints, coefficients, bounds); one-ULP legs on each part; appended trailing zeros are one value; bounds `-0.0` and `0.0` are one value; pickle. The `.npz` round trip (loaded with `allow_pickle=False`) is equal, digest-equal and `array_equal` everywhere. 8 malformed values are refused (NaN coefficient, non-increasing breakpoints, a chart not ending at 1, a negative or infinite bound, mis-shaped bounds, panels disagreeing on groups, one panel too few) | THEOREM | P0 → population and 4 of 4 legs red; B9, no trimming → the appended-zeros pair red | S5.2–S5.4 |
| R4.11 | **The layer.** A cold interpreter that fits, reads and evaluates loads only `{geometry, numerics}`, and no SymPy. By AST every `orpheus` import is under numerics, so no geometry TYPE reaches the module | THEOREM | B13, the file imports `orpheus.data` → red | S1.8 |

**The battery** (`battery/s4_battery.py`: textual arms built from `inspect.getsource` of the prototype and exec'd into the live module before collection; the on-disk arm restored, `diff -q` clean). Scope: the file without the `slow` row, 78 rows, baseline 0 red. Result: 14 arms (P0, B1–B13), every target row red, 0 blind (red sets in `battery/logs/`). Owed on the production module: the same arms re-targeted, and an R4.7 arm.

**What no step-4 gate sees.** Whether the singular set a reader DECLARES is the right one: the representation can only refuse the panel that cannot resolve. The reader's grading from the specification's geometry is step 6's gate. Whether the interpolated function is the reference's own natural extension is the derivation's claim (P3/P4).

### 1.5 Step 5: `PublishedSolution` and `ReferenceCertificate` (re-specified 2026-10-03; gates written)

This section was re-specified on the user's two rulings of 2026-10-03, "Step 5 ruled: no ladder certifies" and "The reading-bound prototype and the user's rulings", and on the orchestrator's API ruling of the same day. The first text, built on `ConvergedLadder`, is in git; the R5 ids below are new.

**The values** (all are `ContentIdentity`, with roster entries):
- **The one Withdrawal.** `Withdrawal(reason, issue)` moves to `orpheus/reference/withdrawal.py`. P0's lock (`withdrawn_generator`, `GeneratorWithdrawn`, the opt-in variable) stays in `derivations/common/withdrawal.py` and imports it. There is no re-export shim: the importers re-point.
- **The standing.** `Current()` and `Standing = Current | Withdrawal`: one sum for a publication and for a certificate.
- **The publication.** `PublishedSolution(specification, printed: Mapping[Observable, Printed], standing)`.
  - `read(observable)` returns the stored `Printed`.
  - It raises `NotPrinted(LookupError)` for an observable the publication does not print, and refuses every read of a withdrawn publication, naming the issue.
  - Construction refuses an observable the specification cannot pose, through the ONE function `admit_observable(observable, specification)` in `orpheus/specification/specification.py`, which the step-6 reader reuses.
- **The printed enclosure.** `Printed.enclosure()` is the printed claim as an outward-rounded `Enclosure` (review note C4).
- **The two establishments.** `Exact(expression, by)` holds the `srepr` text of an exact SymPy expression; a quadratic irrational qualifies. Its `.enclosure()` is the correctly rounded double and `|fl(x) − x|`, rounded up. `DerivedBound(enclosure, method)` is an a-posteriori bound with computable constants. `Establishment = Exact | DerivedBound`.
- **The claim.** `Claim(target, established)`, one per observable, so a target without an establishment is unspellable.
- **The corroboration.** `Corroboration(anchor, independence)`. The anchor is a CURRENT `PublishedSolution` in step 5; step 6 adds `ReferenceSolution`.
- **The refinement.** `Refinement(observable, members: ((parameter, Enclosure), ...))` is falsifying evidence only.
- **The certificate.** `ReferenceCertificate(claims, corroborations, refinements, standing)`. Its `.state` is DERIVED and never stored:
  - it is the `Withdrawal` whenever the standing is one;
  - else it is `Invalid(reasons)` when a claim's bound exceeds its target, when an anchor's printed enclosure is disjoint from a claim's enclosure, or when a refinement's enclosures share no common point;
  - else it is `Valid()`.

**The refinement law is band-free, by construction.** Correct bounds all contain the exact answer, so their enclosures intersect. A disjoint pair proves a bound wrong. The observed orders are returned (`observed_orders()`) and decide nothing. This replaces the withdrawn R5c.5, whose band and safety factor had no derivation.

**Artefacts:** `scratch/reference_architecture/p2/ta/step5/` holds:
- `drafts/_step5.py`: the adapter, every production name once;
- `drafts/test_published.py` and `drafts/test_reference_certificate.py`, which land in `tests/gates/reference/`;
- `s5_union.diff`, which adds the two rosters to S5.9's union.

**Measured first reds** (`[M]` 2026-10-03, worktree of `99862745`, `-O`): both files fail at collection, `ModuleNotFoundError: No module named 'orpheus.reference.published'` (and `.certificate`). pyright reports 0 errors. The battery runs on the production branch.

| id | gate | kind | intended tooth |
|---|---|---|---|
| R5w.1 | One class `Withdrawal` under `orpheus/` and `tests/`, in `orpheus/reference/withdrawal.py`. AST census; input count printed, more than 1000 files; control: `Printed` is found | THEOREM (X4) | a second class minted → red |
| R5w.2 | The modules that bind the name (`derivations.common.withdrawal`, `peierls_nystrom`, `tests._harness.registry`, `tests.conftest`) bind that object, by identity. The lock's machinery stays in derivations | THEOREM | a local copy in one importer → red |
| R5w.3 | The moved law: positive and negative legs (the full gate `tests/gates/test_withdrawal.py` re-points) | THEOREM | the issue check removed → red |
| R5e | `Printed.enclosure()`. Over 10 texts, the printed interval `[t − h, t + h]`, taken in exact decimals, lies inside the enclosure. Each end is within 2 ulp plus 2 ulp of the centre, and the centre is `Printed.value` | THEOREM | no outward step → `0.1` red; the half unit dropped → every row red |
| R5p.1 | `read` returns the stored `Printed`, and an equal-content observable reads the same. An unprinted flux weight or point raises `NotPrinted`, "does not print" | THEOREM | `None` returned instead → red |
| R5p.2 | A withdrawn publication refuses every read, naming the issue; `Current` reads | THEOREM | the standing check removed → red |
| R5p.3 | Admission: a non-specification, an empty mapping, a non-observable key, a non-`Printed` value (a float, an `Enclosure`) and a non-standing are each refused, naming the field | THEOREM (defining refusal) | each check removed → its row red |
| R5p.4 | Unposable observables are refused: an `Eigenvalue` of a source question; a `PointValue` on the infinite medium, outside the geometry, or with a group out of range; a weight with 3 regions or 1 group. Four positive legs construct | THEOREM | each arm of `admit_observable` removed → its row red |
| R5p.5 | ROUTE: `admit_observable` is rebound to a decoy, and a posable AND an unposable publication raise the decoy's error. Activation: the decoy ran at least twice | THEOREM (route, X4) | a second admission inside the publication → the unposable row raises its own error → red |
| R5p.6 | Content identity: `Current` and `PublishedSolution`, with legs on another specification, a printed digit, a locator, an observable and the standing; the mapping's insertion order is not content; pickle | THEOREM | P0 |
| R5c.1 | `Exact` over 5 expressions (1/3, 2/7, 3/4, the golden ratio, a 2-group-k-like quadratic irrational). Against mpmath at 120 digits: the value is the correctly rounded double, the interval contains x, the bound is at most one ulp, and the bound is 0 iff x is a double. `Rational(2, 6)` and `Rational(1, 3)` are one value. NOT claimed: every two spellings of one number are one value (identity is by SymPy's canonical text). Refused: non-`srepr` text, a free symbol, `oo`, `I`, a non-SymPy call, an empty `by`, a float | THEOREM | the bound set to 0 → the inexact rows red; the whitelist removed → the `__import__` row red |
| R5c.2 | `DerivedBound` is its enclosure and its method; a blank method and a non-`Enclosure` are refused | THEOREM | — |
| R5c.3 | A claim meets its target or the state is `Invalid`, naming the observable. Bound below and equal → `Valid` (closed); one ulp above and far above → `Invalid`. Admission: a target ≤ 0 or infinite; an establishment that is a bare `Enclosure` or a production `Measured`; no claims; a non-observable key | THEOREM | `<` for `≤` → the equal row red; the target check removed → the above rows red |
| R5c.4 | The anchor law, decided exactly. Inside → `Valid`; overlapping at the edge by 1e-12 → `Valid`; disjoint by one printed unit, and below → `Invalid`, naming the citation. An exactly touching row cannot be built: a printed half unit `5·10^-k` is never a double | THEOREM | the intersection test dropped → the disjoint rows red |
| R5c.5 | Admissible anchors: a withdrawn publication, a blank independence note, an `Enclosure`, a float, a `Printed` and a production `Measured` are each refused. An anchor sharing no observable with the claims is refused ("no observable"): a vacuous corroboration reads as evidence (X1) | THEOREM (structure) | each refusal removed → its row red |
| R5c.6 | Refinements falsify: model `J_h = 1 + h²` with bounds `2h²` (`Valid`, observed order 2 returned); a disjoint pair → `Invalid`; a perfect order never makes `Valid` a bound above target; a wild order with intersecting enclosures stays `Valid`. Admission: one member; non-decreasing parameters; a non-enclosure member; an unclaimed observable | THEOREM | the common-point test dropped → the disjoint row red; the order made to decide → the wild row red |
| R5c.7 | A `Withdrawal` standing IS the state (by identity) whatever the evidence | THEOREM | the standing read after the evidence → red |
| R5c.8 | `dataclasses.fields(ReferenceCertificate)` is exactly `(claims, corroborations, refinements, standing)`; no stored state | THEOREM (structure) | a stored `state` field → red |
| R5c.9 | No `ConvergedLadder`, `Extrapolated`, `RichardsonEstimate` or `LadderBound` class under `orpheus/` (AST; control: `DerivedBound` found); `Establishment` names exactly `{Exact, DerivedBound}` | THEOREM (the ruling) | a struck class minted → red |
| R5c.10 | Content identity of `Exact`, `DerivedBound`, `Claim`, `Corroboration`, `Refinement` and `ReferenceCertificate`: every part moved, equal pairs, pickle | THEOREM | P0 |
| R5.12 | A cold import of the three step-5 modules loads only `{reference, specification, numerics, data, geometry}`: the reference package never imports derivations | THEOREM (layer) | a derivations import → red |


**Step 5's battery** (`[M]` 2026-10-03, on the production modules at `1db08615`; artefacts in `scratch/reference_architecture/p2/ta/step5/battery/`). Arms are textual transforms of `certificate.py`, `published.py` and `reading.py`, exec'd into the live modules; three on-disk arms were restored and diff-checked. Scope: the two step-5 files, 131 rows, baseline 0 red. Result: P0 (27 red) and 24 arms, every target row red, 0 blind. The targets: target check removed (R5c.3 above-target rows), `>` made `>=` (the equal row), anchor law dropped (R5c.4's 2 disjoint rows), common point dropped (R5c.6 disjoint), standing read after the evidence (R5c.7), `Exact` bound zero (4 inexact R5c.1 rows), srepr whitelist bypassed (2 admission rows), the withdrawn, vacuous-anchor, unclaimed-refinement and unordered-parameter refusals (1 each), the observed order made decisive (the wild row), unprinted read returning `None` (R5p.1), withdrawn read allowed (R5p.2), admission skipped (6 R5p.4 rows plus R5p.5), a second local admission (R5p.5 alone), the empty publication admitted (R5p.3), half unit dropped (10 R5e rows plus R5c.4's edge row), no outward step (R5e `0.99996`, `0.6123`, `0.1`), a second `Withdrawal` class, a copied binding, `published.py` importing derivations, a `ConvergedLadder` minted (1 each). OWED: production's extra law "the claim's enclosure must meet the refinement's common part" (`certificate.py`, `ReferenceCertificate.state`) has no row of its own (the common-point arm disables it too).

### 1.6 Step 6: `ReferenceSolution` and `read(observable)` (re-specified 2026-10-03 against `a3ff64d0`; gates landed `4b724f04`)

**The values** (the orchestrator's API ruling of 2026-10-03):
- `Derivation`, a `runtime_checkable` Protocol, `establish(observable: Eigenvalue | Linear) -> Exact | DerivedBound`: the natural extension evaluated WITH its derived bound, on demand, for any admissible observable. It raises `NotCertified(LookupError)` when it can derive no bound.
- `ReferenceSolution(specification, derivation, certificate: ReferenceCertificate | None)`. All three fields are required. It is not a content value until P3 keys it by its derivation's trace.
- **Construction (qa F2 of step 5).** Every claimed observable passes `admit_observable`. Every anchor answers the holder's specification. A claimed `Ratio` is refused, because a ratio's reading is the quotient of its operands' readings: one definition.
- **`read(observable)`.** `admit_observable` runs first. A `Ratio` reads as `read(numerator) / read(denominator)`. Otherwise `read` returns `derivation.establish(observable).enclosure()`. When the certificate claims the observable, the claim's enclosure and the established one must share a point, or the read is refused ("disagrees"). A derived enclosure is a guarantee on its own (G3); the certificate adds targets and corroboration. `read` answers on an `Invalid` or `Withdrawn` certificate.
- **The first instance.** `ExactInfiniteMediumDerivation` and `exact_infinite_medium_reference(specification)`, in `derivations/common/exact_homogeneous.py`. The k question reads k∞ in the k chart. Flux integrals are read per unit volume under the gauge ⟨νΣf, φ⟩ = 100. Every reading is `Exact`.

**The gates** (`tests/gates/reference/_step6.py`, `test_reference_solution.py`):

| id | gate |
|---|---|
| R6.1 | The fields, all required; the test double is a `Derivation`. Refused: a non-specification, an object without `establish`, a non-certificate, a claimed `PointValue` on the medium (admission), an anchor of another specification, a claimed `Ratio` |
| R6.2 | An unposable observable (a point on the medium, a weight with 2 regions or 3 groups) is refused BEFORE `establish` runs. ROUTE: with `admit_observable` rebound to a decoy, a posable and an unposable read both raise the decoy's error |
| R6.3 | A read equals `establish(observable).enclosure()`, asked once. With a claim that agrees, the reading is the ESTABLISHED enclosure, not the claim's. An unclaimed observable reads, with and without a certificate. An exact `7/5 + 10⁻⁹` against an exact claim `7/5` is refused, "disagree". An `Invalid` or `Withdrawn` certificate still reads |
| R6.4 | ROUTE: a ratio is exactly one `Enclosure.__truediv__` call, equal to `read(a)/read(b)`, and contains the exact ratio |
| R6.5 | RE-POSED WHEN THE USER RULES ON THE UNCERTIFIED READING: an observable the derivation cannot bound raises `NotCertified`, "certif", through `read`, with and without a certificate |
| R6.6 | For 3 mixtures (fuel at 2 groups with upscatter and (n,2n), fission-only at 2 groups, fuel at 4 groups) the certificate is `Valid`. Every reading contains its exact value: k∞ (bound ≤ one ulp), each group flux, the group-0/group-1 ratio, and an unclaimed weighted sum (`(g+1)/4`). The factory refuses a slab, a fixed-source medium, and an eigen question along scattering emission ("k-eigenvalue": the row added after arm E3, see below) |
| R6.7 | A cold import of `reference.solution` loads no derivations |

**Declared limits.** R6.6 compares the reference against `exact_infinite_medium_of`, the function the factory itself calls. So it is a THREADING gate (gauge, group order, weighted sum, chart), not a value gate: the k∞ VALUE rests on `tests/gates/homogeneous/test_kinf_exact_reference.py` (certified against the defining equations). The factory's claims are established by the same derivation `read` consults, so the claim-versus-established check is vacuous for this family by construction (X4). R6.3's disagreement row exercises that check with a test double.

**The battery** (`[M]` 2026-10-03, production at `5af745a8`, byte-identical to `4b724f04`; `scratch/reference_architecture/p2/ta/step6/battery/`). Arms were built from the source of `solution.py` and `exact_homogeneous.py` and exec'd into the live modules; one arm edits a file on disk (restored). Scope: 27 rows, baseline 0 red. 16 arms; 15 turn their target rows red. One was blind and is now repaired: E3 (the factory's question check removed) left the fixed-source refusal green, because the claimed Eigenvalue's admission refuses that specification first (a twin guard). The new row, an eigen question along scattering emission, reddens under E3 alone (`test_reference_solution_with_e3_row.py`). Declared blindness: E2 (groups reversed) leaves the 4-group row green, because `fuel(4)` is group-symmetric (equal fluxes); the 2-group rows catch it.

### 1.7a Step 7a: `VerificationCertificate` and the comparison verb

`orpheus/reference/verification.py`: the certificate holds the observable, production's `ProductionReading`, the reference's `Enclosure`, the tolerance, the algebraic-error `Evidence` (an explicit argument, no default), the verdict kind (`Agreement` or `Order`) and the verdict; returned, never stored.

| id | gate | kind | first red | tooth |
|---|---|---|---|---|
| R7.1 | **A non-`Valid` reference is refused** (`Invalid`, `Withdrawn`: a typed error naming the state and its issue), BEFORE any reading is compared | THEOREM (defining refusal) | defining | the state check removed → red |
| R7.2 | **Role types.** A reference reading that is `Printed` is refused (agreement is against an enclosure); a production reading that is an `Enclosure` or a `Printed` is refused (`TypeError`); a static leg: pyright on a two-line snippet passing a `Measured` as the reference reading reports an error, and its well-typed twin reports 0 (an annotation has no runtime witness, lessons §1, `L59d`) | THEOREM (structure) | defining | the runtime check removed → the runtime row red; the annotation widened → the pyright row red |
| R7.3 | **The reference floor**: `b_ref ≤ tol / 10`, else the certificate reports the floor as the reason, and `require()` raises for the floor FIRST (the order `AgreementCertificate` established, kept as the migration's contract); the boundary row `b_ref == tol/10` holds | THEOREM | defining | the floor dropped → the floor rows read "agrees" → red |
| R7.4 | **The agreement verdict** is `|m − v_ref| + b_ref ≤ tol`, decided EXACTLY (`Fraction`), so a row at the boundary is decidable: one row exactly at `tol` agrees, one ULP above disagrees; one row with `|m − v_ref| ≤ tol < |m − v_ref| + b_ref` disagrees | THEOREM | defining | `+ b_ref` dropped → the third row agrees → red |
| R7.5 | **The algebraic error is recorded under agreement**: `NotYet(564, …)` and `NotApplicable` give a verdict and are carried on the certificate unchanged | THEOREM | defining | refusing `NotYet` under agreement → red |
| R7.6 | **The order verdict requires the algebraic floor**: each rung's `Measured` or `Certified` algebraic error at or below `tol / 10`, else refused as unestablished (`NotYet(564, …)` refused here); the observed orders are RETURNED, one per pair of rungs (a structure, not a bool), and the verdict compares them with the declared order and band | THEOREM | defining | the floor dropped → the `NotYet` order row gives a verdict → red |
| R7.7 | **Branch coverage** of the certificate as `test_certificate_reaches_every_branch` does today (no reference bound, a loose bound, disagreement, agreement), re-posed on the new type, so the migrated row keeps its branches | THEOREM | defining | — |

### 1.7b Step 7b: the migration of `certify_agreement` (adjusted 2026-10-03 to the user's ruling 3)

The trajectory-resolvent sphere and cylinder stay UNCERTIFIED in P2: no derived bound exists for either; that is #566 and #516's work, in P4. Under ruling 3 neither can anchor a `VerificationCertificate`. Their rows therefore re-pose as COMPARISONS against an uncertified reference, keeping their current tolerances: a RECORD-kind comparison that names the reference uncertified, never a verification claim. The cylinder's strict xfail (awaiting a bound) stays, re-keyed to the refusal R7.1 raises for a non-`Valid` reference. The shape rows' re-posing as per-cell ratio observables stands. `AgreementCertificate`, `certify_agreement` and `CYLINDER_3REG_REFERENCE_BOUND` still retire with the step; the sphere's ladder bound (`sphere_3reg_reference_bound`) becomes the tolerance's documented provenance, not a certificate.

### 1.8 Step 8: `Rate(cells, weight)`

`Rate(cells: CellCoefficient, weight)` lives in `orpheus/specification/` (it names the data layer's `CellCoefficient`, which `numerics` cannot import); `Specification.resolve(observable)` maps it to the canonical `FluxIntegral` (for a fission-emission cell, the weight contracted with the emission spectrum and multiplied by `νΣ_f` region by region; for a scattering-emission cell, contracted with the transfer matrix). The gates the user named: two spellings of one functional (a `Rate` and the `FluxIntegral` it resolves to; `every(fission emission)` and the explicit cells) digest equal after resolution; a cell on a material outside the specification's materials is refused (reusing `CellCoefficient.resolve`'s refusal: a route gate); ≥ 2 groups with a non-trivial emission spectrum, since `χ = (1, 0)` nulls the contraction.

**Recommendation: land it with posing unit 3 (#526), not in P2.** Today's cells are emission channels only (the P1 ruling); removal cells (capture, an absorber search, the `X_x` rates of the k identity) wait for #526's reaction grid. A resolver written now covers half the cell kinds and is re-opened by #526 for the other half. **`FluxIntegral` alone suffices for P2's first consumers**: the step-7b migration reads the eigenvalue and region- or cell-averaged group fluxes (R6.6 and §1.7b), and no P2 row reads a reaction rate.

## 2. The decisions the design left to the test-architect, with reasons

1. **What `Rate`'s `cells` ranges over:** one `CellCoefficient` (a set of reaction cells), never a region set. Spatial restriction lives in the weight (a `RegionwiseConstant` with zero rows), so an observable has one spelling; a region set beside the weight would be a second definition of the same restriction (X4). Superseded in part by the user's Q2 ruling: `Rate` resolves to `FluxIntegral` (§1.8).
2. **The homes** (`test_layer_imports.py`'s layers): `numerics/enclosure.py` (`Enclosure`); `numerics/outcome.py` (`ProductionReading` beside `Measured`; `ExitReport`); `numerics/observable.py` (the closed sum); `numerics/piecewise_chebyshev.py` (the panels); `specification/` (`Rate` and its resolution, readable by production and reference alike); `reference/` (`Printed`, `ReferenceReading`, `PublishedSolution`, `ReferenceCertificate`, `ReferenceSolution`, `VerificationCertificate`, and P0's `Withdrawal` moved down, R5c.6).
3. **The rename:** `ExitCertificate` → `ExitReport`, the slot `certificate` → `exit_report` (§1.2).
4. **The migration set** of `certify_agreement`: the 4 files of §1.7b, with the sphere `Valid` and the cylinder `Invalid`, and the shape rows re-posed as per-cell ratio certificates.

## 3. The ladder, as a graph

Foundations reused: S5.1–S5.4, S5.7, S5.9 (the encoder, the seed harness, the RECORD pattern, the roster), S6.8, S6.21 (the mesh-free functions, the layer probe), S7.13, S8's admission (`_admit_datum`, `CellCoefficient.resolve`), the linter's parametrised row. New rungs, bottom up: R1.1 → R1.2 → R1.3 → R1.4 (the enclosure); R1.15 → R1.10 (the printed text); R1.6 → R1.16 (the sums); R3.1 → R3.2–R3.6 (the observables, on R1.16 and S6.8); R4.1 → R4.2–R4.4, R4.5 → R4.6 → R4.7 (the panels); R5p.* and R5c.* on S1 and S3; S6.* on S3, S4, S5; S7.* on S6; §1.7b's migrated rows on S7. Each test declares its edges with `@pytest.mark.rests_on` (the step-1 drafts carry theirs). Improved tests: the four `certify_agreement` files (§1.7b); `test_content_identity.py`'s union and package walk; the linter's registration. No step-1 test duplicates an existing one: the census found no value-with-bound type in the tree.

## 4. Refuted or rejected candidates

- **`Printed(value, digits, citation)`:** rejected (Q1): ambiguous `digits`, two fields that can disagree, and a float that cannot carry a trailing zero.
- **A midpoint centre for the quotient:** rejected for R1.3: the centre would move with the operands' bounds, so a ratio would not read as the ratio of its readings.
- **A linearised bound for the quotient** (`(α|b| + β|a|)/b²`): refuted as an enclosure, `[M]` A1 reds the containment law on the population; it under-encloses whenever the denominator's relative width is not small.
- **A valid-only quotient law** (containment without tightness): rejected; a bound of 1e300 satisfies it and fails every floor (R1.4 is the second side).
- **Region sets as `Rate`'s cells:** rejected, a second spelling of the weight's restriction (§2.1).
- **`Rate` in `numerics` with an opaque key:** rejected after Q2: `Rate` resolves through the materials, so it lives with the specification.
- **The step-7b shape rows kept as a `max` over cells:** rejected: not a linear observable, so no enclosure propagates to it.

## NEEDS

1. **`ConvergenceCertificateError`**: rename with step 2 (proposed `ConvergenceClaimError`), or keep (it names the convergence claim an `IterationRecord` makes, a third subject)? The user asked for few certificates.
2. **R5c.5's two declared inputs**: the acceptance band on a ladder's observed order and the safety factor on its Richardson estimate. Roache's GCI factor 1.25 is from a paywalled source not read (the G6 literature note); a derivation, or the user's ruling, is owed before step 5's gates are final.
3. **Step 7b's size**: the migration needs the trajectory-resolvent sphere written into panels (its natural extension) and Garcia 2021 Table 5 entered as a `PublishedSolution`; confirm it stays in P2 rather than opening P4's first family.
4. **A reading of a non-`Valid` reference**: R6.1 lets `read` answer on an `Invalid` or `Withdrawn` reference (only the `VerificationCertificate` refuses). Confirm, or rule that `read` refuses too.
5. **The step-1 prototype's names are the test-architect's** (`Enclosure.__truediv__` with outward rounding via `math.nextafter`; `Printed.value`, `Printed.half_unit`; `ProductionReading` in `numerics/outcome.py`). The drafts name them, so a rename is one edit in each of the 2 draft files and in `battery/s1_battery.py`.

## Review notes carried into later steps (the elegance review of steps 1–2, 2026-10-03)

- **Step 7a (S4):** `ProductionReading = Measured` shares its class with four diagnostic producers (`solver.py` `_balance_evidence` arms, `solution.py:347`), so `exit_report.balance` narrowed to `Measured` would type-check as a production reading passed in the wrong slot. The `VerificationCertificate` therefore takes `(answer, observable)` and calls `answer.read(observable)` itself, so the reading verb is the only producer of a production reading; the algebraic-error `Evidence` remains its separate explicit argument.
- **Step 5 (C4):** `Printed.half_unit` (a `Decimal`) and `Enclosure.bound` (a `float`) both spell the distance to the exact value. At the first consumer that compares a `Printed` with anything (the `ReferenceCertificate`'s published anchor), give `Printed` one verb, `enclosure()`, rounded outward, rather than an `isinstance` branch in each consumer.
- **Renamed in step 2's review round:** the `Evidence` member `Certified` is `Asserted` (G4 reserves "certificate" for the two certificates; the member records a bound the convergence-claim check asserted). §0 item 3's floor therefore reads `Measured` or `Asserted`. The pre-existing weld in the admissibility member (the measured k stored in `Asserted.bound`) is #565.
- **Step 6 (the elegance review of step 3, finding 1):** an observable's fit to a problem is admitted ONCE, by the specification: `_admit_observable(observable, spec)` beside `_canonical_question`, an exhaustive `match` (a `FluxIntegral` through `_admit_datum`; a `Ratio` through both operands; a `PointValue`'s group against `n_groups` and its position against the geometry's extent, refused on the infinite medium; an `Eigenvalue` only for an `Eigen` question). Every answer's `read` calls it; gate R6.5 widens to cover each arm (`[M]` today `PointValue(1e9, 999)` is admitted, with no problem to fit).
- **Step 3, settled in its review round:** a `Ratio`'s operands are `Linear = FluxIntegral | PointValue` (a quotient of two linear functionals is the only one independent of an eigen answer's scale), so ratios do not nest; R3.3/R3.4/R3.5 re-posed accordingly. The four members are `@final`. The closed-sum refusal is one parser, `scalars.parse_member`, and the owner prefix is restored (`FluxIntegral: the weight is a mesh-free function ...`).

- **Step 6 (qa of step 5, F2):** a `ReferenceCertificate` holds no specification, so its corroborations cannot check "the same problem". The `ReferenceSolution` that holds the certificate enforces it: every anchor's specification equals the holder's (for `PublishedSolution` anchors as well as `ReferenceSolution` ones), and every claimed observable passes `admit_observable`. `[M]` today a certificate claiming k = 1/3 on a slab, anchored by an infinite-medium publication printing 0.33333, reads `Valid`.
- **Correction (qa of step 5, F4):** R5c's anchor law is decided on the OUTWARD-ROUNDED ends (`Enclosure.ends()`, `common_part`), not exactly: it never calls two intersecting intervals disjoint, and it can miss a disjointness of about one ulp. The review round replaced the three pairwise checks with one family law (every enclosure of one observable shares a point, Helly in one dimension), keeping the specific reasons.
