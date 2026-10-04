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
- `Derivation`, a `runtime_checkable` Protocol, `evaluate(observable: Eigenvalue | Linear) -> Evaluation`, with `Evaluation = Establishment | Uncertified` (`Establishment = Exact | DerivedBound`): the natural extension, evaluated on demand for any admissible observable, with its derived bound where the family has one. [REMEDIED 2026-10-03 by step 7b.1] The verb was `establish -> Exact | DerivedBound`, and a family without a bound raised `NotCertified(LookupError)`; step 7b.1 retired both (§1.7b.1).
- `ReferenceSolution(specification, derivation, certificate: ReferenceCertificate | None)`. All three fields are required. It is not a content value until P3 keys it by its derivation's trace.
- **Construction (qa F2 of step 5).** Every claimed observable passes `admit_observable`. Every anchor answers the holder's specification. A claimed `Ratio` is refused, because a ratio's reading is the quotient of its operands' readings: one definition.
- **`read(observable)`.** `admit_observable` runs first. A `Ratio` reads as `read(numerator) / read(denominator)`. Otherwise `read` evaluates the observable: an `Exact` or a `DerivedBound` reads as its enclosure, an `Uncertified` as itself (§1.7b.1). When the certificate claims the observable, the claimed observable is read once at construction, through the same route: it must not read `Uncertified`, and the claim's enclosure and the established one must share a point, or the reference is refused ("disagrees"). A derived enclosure is a guarantee on its own (G3); the certificate adds targets and corroboration. `read` answers on an `Invalid` or `Withdrawn` certificate.
- **The first instance.** `ExactInfiniteMediumDerivation` and `exact_infinite_medium_reference(specification)`, in `derivations/common/exact_homogeneous.py`. The k question reads k∞ in the k chart. Flux integrals are read per unit volume under the gauge ⟨νΣf, φ⟩ = 100. Every rational reading (k, a flux integral of a tabulated weight) is `Exact`; a symbolic weight whose value `Exact` refuses as `Uncertifiable` (cancellation) reads `Uncertified`, and every other refusal of `Exact` propagates (§1.7b.1). [REMEDIED 2026-10-03 by step 7b.1] This read "Every reading is `Exact`".

**The gates** (`tests/gates/reference/_step6.py`, `test_reference_solution.py`):

| id | gate |
|---|---|
| R6.1 | The fields, all required; the test double is a `Derivation`. Refused: a non-specification, an object without `evaluate` (and one spelling the retired `establish`), a non-certificate, a claimed `PointValue` on the medium (admission), an anchor of another specification, a claimed `Ratio` |
| R6.2 | An unposable observable (a point on the medium, a weight with 2 regions or 3 groups) is refused BEFORE `evaluate` runs. ROUTE: with `admit_observable` rebound to a decoy, a posable and an unposable read both raise the decoy's error |
| R6.3 | A read equals `evaluate(observable).enclosure()`, asked once. With a claim that agrees, the reading is the ESTABLISHED enclosure, not the claim's. An unclaimed observable reads, with and without a certificate. An exact `7/5 + 10⁻⁹` against an exact claim `7/5` is refused, "disagree". An `Invalid` or `Withdrawn` certificate still reads |
| R6.4 | ROUTE: a ratio is exactly one `Enclosure.__truediv__` call, equal to `read(a)/read(b)`, and contains the exact ratio |
| R6.5 | RE-POSED by step 7b.1 (the user's ruling of 2026-10-03; §1.7b.1): an observable the derivation evaluates `Uncertified` reads `Uncertified`, equal to what `evaluate` returned, with and without a certificate; a ratio with an uncertified operand, in either position or both, reads `Uncertified(a.value / b.value)` through one `Uncertified` dunder; a certificate claiming an uncertified observable is refused "uncertified", BEFORE the agreement check; what a derivation returns is admitted as an `Evaluation` (`TypeError`, "evaluation") |
| R6.6 | For 3 mixtures (fuel at 2 groups with upscatter and (n,2n), fission-only at 2 groups, fuel at 4 groups) the certificate is `Valid`. Every reading contains its exact value: k∞ (bound ≤ one ulp), each group flux, the group-0/group-1 ratio, and an unclaimed weighted sum (`(g+1)/4`). The factory refuses a slab, a fixed-source medium, and an eigen question along scattering emission ("k-eigenvalue": the row added after arm E3, see below) |
| R6.7 | A cold import of `reference.solution` loads no derivations |

**Declared limits.** R6.6 compares the reference against `exact_infinite_medium_of`, the function the factory itself calls. So it is a THREADING gate (gauge, group order, weighted sum, chart), not a value gate: the k∞ VALUE rests on `tests/gates/homogeneous/test_kinf_exact_reference.py` (certified against the defining equations). The factory's claims are established by the same derivation `read` consults, so the claim-versus-established check is vacuous for this family by construction (X4). R6.3's disagreement row exercises that check with a test double.

**The battery** (`[M]` 2026-10-03, production at `5af745a8`, byte-identical to `4b724f04`; `scratch/reference_architecture/p2/ta/step6/battery/`). Arms were built from the source of `solution.py` and `exact_homogeneous.py` and exec'd into the live modules; one arm edits a file on disk (restored). Scope: 27 rows, baseline 0 red. 16 arms; 15 turn their target rows red. One was blind and is now repaired: E3 (the factory's question check removed) left the fixed-source refusal green, because the claimed Eigenvalue's admission refuses that specification first (a twin guard). The new row, an eigen question along scattering emission, reddens under E3 alone (`test_reference_solution_with_e3_row.py`). Declared blindness: E2 (groups reversed) leaves the 4-group row green, because `fuel(4)` is group-symmetric (equal fluxes); the 2-group rows catch it.

### 1.7a Step 7a: `VerificationCertificate` and the comparison verbs (re-specified 2026-10-03; gates landed `9db19fd8`)

**The values** (`orpheus/reference/verification.py`, the orchestrator's API ruling of 2026-10-03):
- `ProductionAnswer`, a Protocol with `read(observable) -> ProductionReading`.
- `verify_agreement(answer, observable, reference, tolerance, algebraic_error) -> VerificationCertificate` and `verify_order(answers, observable, reference, tolerance, algebraic_errors, order, band) -> OrderVerification` (G4: two certificates only, so not "OrderCertificate"). Neither verb has a default.
- Both verbs hold the reference `Valid` BEFORE reading anything: `ReferenceNotValid` names a missing certificate, an `Invalid` state (its first reason) or a withdrawal (its issue).
- Both call `answer.read(observable)` themselves (review note S4); a reading that is not `Measured` is refused.
- `VerificationCertificate` is a frozen dataclass, not `ContentIdentity`. `floor_holds` (b_ref ≤ tol/10) and `agrees` (|m − v_ref| + b_ref ≤ tol) are decided exactly. `require()` raises for the floor first. The algebraic error is recorded under agreement, never required.
- `verify_order` requires every algebraic error to be `Measured` or `Asserted` ≤ tol/10 (`Unestablished` otherwise; `NotYet` and `NotApplicable` refused), the reference bound ≤ tol/10, and every error ≥ tol ("unresolved"). It returns the observed orders, and `holds` against the caller's declared order and band.
- **The production reading.** `HomogeneousResult.read` (`orpheus/homogeneous/solver.py`), returning `Measured`:
  - Eigenvalue reads k∞.
  - A one-region `FluxIntegral` reads per unit volume in the gauge νΣf·φ = 100; a `Symbolic` weight is read through `without`.
  - A `Ratio` reads as the quotient of its operands' readings.
  - `PointValue` is refused.
- **DECLARED LIMIT.** Nothing pairs the answer with `reference.specification`: production results hold no specification until P4's projection, so the caller pairs them.

**Gates** (`tests/gates/reference/_step7.py`, `test_verification.py`, R7.1–R7.10, 44 rows):
- R7.1: the Valid refusals, before any read.
- R7.2: the roles, the arguments, and no defaults.
- R7.3: the floor, inclusive, raised first.
- R7.4: agreement exactly at dyadic boundaries.
- R7.5: any `Evidence` is recorded.
- R7.6: the order verdict and its three floors.
- R7.7: returned, never stored.
- R7.8: `HomogeneousResult.read` and its refusals.
- R7.9: production homogeneous against the exact medium, three mixtures, at 1e-13 relative (`[M]` production within 1.65e-16 relative), and the X1 negative: νΣf scaled by 1 + 1e-11 disagrees.
- R7.10: the layer.

**Battery** (`[M]` 2026-10-03 at `9db19fd8`, also run at `b445b68a`; `scratch/reference_architecture/p2/ta/step7a/battery/`). 23 arms, baseline 0 red, every target red, 0 blind:
- Validity unchecked: 4 rows. Invalid accepted: 1. Withdrawn accepted: 1. Reading before the validity check: 4.
- The reading-type check dropped: 4. The algebraic-error check dropped: 2.
- The floor made strict: 1. The floor dropped: 1. The reference bound dropped from agreement: 2. Agreement made strict: 2. Agreement checked before the floor: 1.
- `NotYet` establishing: 2. The algebraic floor dropped: 2. Unresolved errors admitted: 1. The band ignored: 1. The order verb's reference floor dropped: 1.
- The reading ignored: 3, including the R7.9 negative.
- Ratio reversed: 3. Groups reversed: 3. A point read as 0: 1. The region check dropped: 1. 1/k returned: 4.
- `verification.py` importing homogeneous: 1.
- Declared: the 4-group R7.9 row is blind to a group reversal (`fuel(4)` has equal group fluxes).

### 1.7b.1 Step 7b part 1: the uncertified reading (specified, written, run and mutated 2026-10-03)

**The ruling** (the user, 2026-10-03, "Step 7b ruled" in `.claude/plans/reference_cache.md`). A reference whose family cannot derive a bound reads `Uncertified(value)`, a third member of `ReferenceReading`. The `VerificationCertificate` refuses it. A test compares against it through an explicit uncertified comparison at its current tolerance, so the weaker claim is visible in the code.

**The values** (the orchestrator's declared API, 2026-10-03; the landed code follows it):
- `reading.Uncertified(value)`: `@final`, frozen, `ContentIdentity`; the value is admitted by `parse_finite_real`. It has no `enclosure()`. `Uncertified / (Enclosure | Uncertified)` and `Enclosure / Uncertified` are `Uncertified(a.value / b.value)`, so a quotient with an uncertified operand drops the other operand's bound. `ReferenceReading = Enclosure | Printed | Uncertified`.
- `solution`: `NotCertified` is retired. The `Derivation` verb is `evaluate(observable) -> Evaluation`, with `Evaluation = Establishment | Uncertified` (review round). `read` returns `Enclosure | Uncertified`. Construction reads every claimed observable through `_read`, the one route from evaluation to reading (review round), and refuses one that reads `Uncertified`, before the claim-versus-establishment agreement. The reader admits what an open-protocol derivation returns as an `Evaluation`. `orpheus.reference` re-exports `Uncertified`.
- `certificate.Uncertifiable(ValueError)` (review round): the ONE refusal of `Exact` meaning "cannot certify" (strict `evalf` raising `PrecisionExhausted`). It carries `approximation`, the expression's non-strict `evalf` at the working precision `_EXACT_DIGITS` (60), as a float. Every other refusal of `Exact` (not a constant, not real, not provably finite) is a defect of the expression.
- `ExactInfiniteMediumDerivation.evaluate` catches `Uncertifiable` alone and returns `Uncertified(refusal.approximation)`; every other refusal propagates. Its unreachable `PointValue` arm raises `ValueError` ("position"). The factory matches on `evaluate`'s result (`certify` is retired). The uncertified value carries no guarantee: the rows pin the fallback's DEFINITION (the working precision), not an accuracy.
- `verification`: both verbs raise `ReadingUncertified(ValueError)` ("uncertified") when a `Valid` reference reads the observable `Uncertified`, before production is read; `ReferenceNotValid` names the standing only (no certificate, or one that is not `Valid`) (review round). `compare_uncertified(answer, observable, reference, tolerance) -> UncertifiedComparison` has no defaults. It refuses a certificate that exists and is not `Valid` before anything is read. It compares exactly what the verbs cannot anchor on: an `Uncertified` reading, or ANY reading of a reference with no certificate, an `Enclosure` included (its bound unused). It refuses an `Enclosure` read from a reference whose certificate is present, `ReadingCertified` (naming `verify_agreement`), so the two verbs cover every state of a reference in good standing and a test is forced up to a verification the day its family is certified (review round; the first version refused every certified reading). Production is read last, `Measured` only. `agrees` is `|m − v| ≤ tol`, decided exactly; `require()` raises `Disagreement` naming the "unverified reference value".

**The gates** (88 new rows and 1 strict xfail after the review round, 2 amended; files `tests/gates/reference/test_readings.py`, `test_reference_solution.py`, `test_verification.py`, adapters `_step6.py`, `_step7.py`):

| id | rows | gate | first red on `7d0258d1` |
|---|---|---|---|
| R1.16 (amended) | 1 | `ReferenceReading` names exactly `{Enclosure, Printed, Uncertified}`; an `Uncertified` is not a production reading | the alias names two members |
| R1.19 (amended) | 1 | the seed-stability script prints an `Uncertified` digest and hash too (four values, seeds 1 and 2) | the import fails |
| R1.15 (roster) | 7 | `Uncertified` joins the readings roster, and so S5.9's union: parts `(value,)`, perturbations (another value, one ulp, the sign), equal pairs (`-0.0`/`0.0`, `2`/`2.0`), pickle | collection error: the module has no `Uncertified` (it also stops `test_content_identity.py` from collecting) |
| R1.20 | 14 | admission: a finite real stored as a `float` (an int, a NumPy double, `-0.0` folded); NaN, ±inf, an int beyond 2⁵³, a bool, a string, `None`, a complex refused with the shared parser's keyed message, owned by `Uncertified`; fields `(value,)`, `__final__`, frozen, a `ContentIdentity`, no `enclosure`; not equal to an `Enclosure`, a `Printed` or a `Measured` of the same number (a set holds all four) | collection error |
| R1.21 | 16 | the quotient: six operand pairs (both orders, a non-dyadic `1/3`, a negative) give `Uncertified` with the float quotient of values; division by a zero value (an uncertified `0.0` or `-0.0`, an enclosure of value 0) raises `ZeroDivisionError`; a `Printed`, a bare float or a `Measured` operand is a `TypeError`; an overflowing quotient is refused, "infinite" | collection error |
| R1.22 | 1 | closure: the four quotients of `{Enclosure, Uncertified}` pairs lie in the three-member sum, with kinds `Enclosure, Uncertified, Uncertified, Uncertified` | collection error |
| R6.1 (+1) | 1 | a derivation spelling the retired verb `establish` is not a `Derivation` (`TypeError`) | the HEAD protocol accepts it |
| R6.5 (re-posed) | 10 | the reading is `Uncertified`, equal to `evaluate`'s, asked once, with and without a certificate; three ratio rows (route: one `Uncertified` dunder call); the claim refusal "uncertified", with a discrimination row whose claim would also disagree (message has no "disagree"); three non-`Evaluation` returns refused "evaluation", and `Evaluation` names exactly three members | the table derivation's verb is `evaluate`, so HEAD refuses it as a derivation |
| R6.8 | 10 + 1 xfail | the exact family reads the Machin zero `atan(1/2) + atan(1/3) − π/4` (constant, so admitted; `Exact` refuses it, asserted as the activation) `Uncertified`, equal to `evaluate`'s, within 1e-100 of 0 (`[M]` 1.94e-121), the certificate still `Valid`; control: k, an indicator weight and a non-dyadic weight still evaluate `Exact`, and a ratio over the Machin weight reads `Uncertified`; the `PointValue` arm raises `ValueError`, "position" Review round: the weight `(z + 10⁻¹³⁰, 0)` (qa F2) reads `Uncertified` within 1e-15 relative of the exact `10⁻¹³⁰ φ₀` (`[M]` 7.9e-17), closing the blind arm X2; the weight `(0, z + 10⁻ᵏ)` at k = 130 and 160 pins the working precision (`[M]` 1.7e-18 and 5.5e-17 at 60 digits; at `evalf`'s default 15 digits −1.21e-122 against +3.34e-128, and 3.6e35 times off); `1/z` and `acos(2)` raise a `ValueError` that is not `Uncertifiable` ("finite real constant"), never read uncertified; `Uncertifiable` is a `ValueError` raised on the Machin zero and not on `acos(2)`. STRICT XFAIL, the challenge (#568, P4): `(z + 10⁻²⁵⁰, 0)` within 1e-15 relative; `[M]` 2.8e59 times off at 60 digits, 3.7e-11 at k = 180. Ruled out of scope for 7b.1: no precision schedule earns a guarantee, and #568 certifies such values with an absolute radius from `evalf`'s accuracy record. | HEAD has no `evaluate`; the arm raised `NotCertified` |
| R7.11 | 3 | both verbs, and the certificate built directly, refuse an uncertified reading of a `Valid` reference, `ReadingUncertified` ("uncertified"; unrelated to `ReferenceNotValid` by subclassing, review round) with production unread (a counting double), the reference read (activation); a ratio over an uncertified flux is refused too; positive leg: the same reference verifies its certified observable | the table derivation is not a HEAD derivation |
| R7.12 | 4 | `Invalid` and `Withdrawn` certificates refused before any read (derivation and answer both uncalled); no certificate and a `Valid` one are compared, each side read once, the readings held (the answer's by identity) | HEAD has no `compare_uncertified` |
| R7.13 | 7 | an ENCLOSED reading of a reference whose certificate is present (exact k, a derived bound, the exact flux of a family whose other readings are uncertified) is refused, `ReadingCertified` naming `verify_agreement`, production unread; `ReadingCertified` and `ReferenceNotValid` are unrelated by subclassing; an uncertified ratio is compared. Review round: an enclosed reading of a reference with NO certificate (exact k, an exact flux, a derived bound 0.25 wider than the tolerance) is compared on its value alone: within tol/2 agrees, at 2·tol disagrees (these were refused before the review round) | as R7.12 |
| R7.14 | 6 | inclusive at a dyadic boundary on both sides, one ULP beyond refused; exactness: v = −0.1, m = 1e17, tol = 1e17, where the float distance rounds onto the tolerance (agrees in floats, disagrees exactly), with a one-ULP control | as R7.12 |
| R7.15 | 1 | `require()` raises `Disagreement` ("unverified reference value", review round) exactly when the comparison does not agree; no floor fires even at a tolerance of 2⁻¹⁰⁰⁰ | as R7.12 |
| R7.16 | 6 | production reading `Measured` only: an `Enclosure`, a `Printed`, an `Uncertified`, a bare float or an `Asserted` from `read` is a `TypeError` ("production reading"); a reading is not an answer | as R7.12 |
| R7.17 | 6 | tolerance zero, negative, infinite, a string; a non-reference refused; the verb and the class take exactly `(answer, observable, reference, tolerance)` with no defaults | as R7.12 |
| R7.18 | 1 | frozen, `@final`, not `ContentIdentity`, not a `VerificationCertificate`; fields the four arguments plus `reading` and `reference_reading` (`init=False`); `agrees` a property; no `floor_holds`; two calls, two objects | as R7.12 |

R7.10 (the layer: a cold import of `reference.verification` loads no method package) holds unchanged on the landed code. The first reds were measured `[M]` on a detached worktree of `7d0258d1` with the new test files copied in (`scratch/reference_architecture/p2/ta/step7b1/head_run.log`): 101 failed, 170 passed, and 1 collection error (`test_readings.py`). Every pre-existing row that builds the table derivation also reds there, because the adapter's verb is renamed; that red is the rename, not a defect. On the landed code: 369 of 369 rows pass (`.venv/bin/python -O -m pytest tests/gates/reference -p no:cacheprovider`), and pyright reports 0 errors and 0 warnings on `tests/gates/reference/`.

**The battery** (`scratch/reference_architecture/p2/ta/step7b1/battery/`, plugin `b7_battery.py`, driver `run_battery.sh`, `battery_summary.txt`, one log per arm; the first round's table is `battery_summary_round1.txt`). Module arms exec a transformed source into the live module. The `Uncertified` and `Exact` arms patch transformed methods onto the live class, because re-executing their module would mint a second class that the `Evaluation`, `Establishment` and `ReferenceReading` unions do not name (lessons `L104`). No file is written; the production and test sources were sha-identical before and after each run.

Round 1 (`[M]` 2026-10-03, before the review round): 32 arms; 31 reddened their targets; X2 (the fallback read as 0.0) was blind.

Review round (`[M]` 2026-10-03, uncommitted working tree after the review round; scope `tests/gates/reference`, 377 passed and 1 strict xfail, baseline 0 red): **36 arms, 36 redden their target rows, 0 blind.** The positive control S1 (an uncertified value read as an exact enclosure) reddens 18 rows (R6.5, R6.8, R7.11–R7.13); it reddened 27 before the review round, when R7.14–R7.16 ran on references that only an uncertified reading made comparable.
- `Uncertified`: admission dropped, 10 (R1.20, R1.21); infinity admitted, 3; the reverse quotient keeps the enclosure, 6; quotient reversed, 10; reverse quotient reversed, 5; any operand, 3; an `enclosure()` added, 3; not `@final`, 1; the sum without `Uncertified`, 2.
- `Exact` (new): the fallback at `evalf`'s default precision, 3 (the two working-precision rows, and the strict xfail turns XPASS, a red).
- `solution`: the claim check dropped, 2; after the agreement, 1 (the discrimination row); the evaluation unadmitted, 3.
- The exact family: the fallback removed, 5; the fallback read as 0.0 (X2), 3, the blindness closed by qa F2's row; every evaluation uncertified, 22; the catch widened back to `ValueError` (new), 2 (`1/z`, `acos(2)`); the point arm reads 0, 1.
- `verification`: the uncertified anchor admitted, 3 (R7.11); the certificate reads production first, 2; the order verb reads production first, 2; the comparison's standing unchecked, 2; checked after the read, 2; `ReadingCertified` unchecked, 3; `ReadingCertified` on a reference with no certificate (new), 3 (R7.13's no-certificate rows); the comparison reads production first, 5; a certificate required, 22; agreement strict, 3; agreement in floats, 1; `require` never raises, 1; the reading type unchecked, 6; `ReadingCertified` a subclass of `ReferenceNotValid`, 3; `ReadingUncertified` a subclass of `ReferenceNotValid` (new), 1; the tolerance unparsed, 4; the tolerance defaulted, 1.
- [REMEDIED 2026-10-03 by the review round] Round 1 declared X2 blind: the only uncertifiable witness then known was an exact zero. qa F2's nonzero witness `(z + 10⁻¹³⁰, 0)` closes it. The working-precision rows exposed that the 15-digit fallback could return a wrong-sign value; the review round moved the fallback to the working precision (carried by `Uncertifiable.approximation`), and the strict xfail records that any fixed precision fails for a deep enough offset (#568).

**Declared limits.** The table double, not a real family, carries every nonzero uncertified value: P2 ships no family without a bound, and the trajectory-resolvent sphere and cylinder are §1.7b's. Nothing pairs the answer with the reference's specification (P4).

### 1.7b.2 Step 7b part 2: the migration of `certify_agreement` (re-specified 2026-10-03 against the working tree of `feature/reference-uncertified-reading`)

Written against the landed tree (`main` `7d0258d1`) plus part 1's API as it stands uncommitted in the working tree (`Uncertified`, `Derivation.evaluate`, `ReferenceSolution.read -> Enclosure | Uncertified`, `compare_uncertified`, `ReadingCertified`; part 1's gates are §1.7b.1's, none are re-gated here). Probes, logs and the census are in `scratch/reference_architecture/p2/ta/step7b2/` (`census_ast.py`/`.txt`, `probe1`–`probe6`). Gate ids are `R7b2.<n>`, test functions `test_r7b2_<n>_…`.

**What part 2 is.** Under the user's ruling 3 (2026-10-03) no trajectory-resolvent reference has a derived bound (#566; the cylinder also #516), so none can anchor a `VerificationCertificate`, and under the step-5 ruling the sphere's ladder "bound" is no bound either: the sphere is uncertified too. Every row of the migration set either becomes an explicit `compare_uncertified` at its current tolerance, or stays a strict xfail whose expected failure is now the verbs' own refusal, or stays a RECORD. Nothing in part 2 is a verification claim.

#### The census (question 1)

`[M]` by AST (`census_ast.py`: every top-level `test_*` in the 4 files, the helper symbols it reads transitively through module-level helpers, its decorators). 26 test functions; 16 of 26 read a `_certified_agreement` or ladder symbol; 9 of 26 call `certify_agreement`; 5 carry `awaits_cylinder_bound` (7 collected cases); 4 call `assert_record`. Positive control and second instrument: Nexus `callers(certify_agreement)` returns the same 9 functions (the sphere k row is in both). The ladder module `tests/gates/derivations/_trajectory_resolvent_ladders.py` imports the helper in 2 functions (`_sphere_reference`, `_cylinder_reference`: `aba_xs_2g`, `ABA_RADII`) and owns `sphere_3reg_reference_bound` and `tolerance_for`.

The rows that change (11 collected cases in 9 functions, plus the 5 harness rows):

| row (file :: function) | compares today | tolerance and its provenance | `[M]` reading today | re-posed as |
|---|---|---|---|---|
| phase_c :: `test_sphere_3reg_k_against_trajectory_resolvent` | `abs(k_ref − k_sn)/k_sn`; SN Gauss–Legendre 32, 40 equal-volume cells, live; reference sphere (n_r, n_μ) = (36, 96) | 4e-3 relative = `tolerance_for(SN residual 1.5e-5, ladder estimate 3.2e-4)` | 6.51e-5 (k_ref 1.381169542313358, k_sn 1.381079639349644) | `compare_uncertified(sn, Eigenvalue(), ABA_SPHERE_REFERENCE, 4e-3 × 1.38)`; the ladder estimate is the tolerance's documented provenance, not a bound |
| phase_c :: `test_sphere_3reg_flux_shape_against_trajectory_resolvent` | `max_{g,i} |ŝ_SN − ŝ_ref| / max ŝ_ref`, ŝ the cell average over the SN cell gauged to unit total fission production; reference read as the NODAL φ through a per-region cubic spline | 2e-2 = `tolerance_for(7.3e-3, 1.4e-3)` | 4.361e-3 (nodal spline); 4.325e-3 through the natural extension | 80 observables `Ratio(FluxIntegral(cell_i, group g, weight 1/V_i), FluxIntegral(νΣf(r)))`, each `compare_uncertified` at τ = 2e-2 × 0.249 (0.249 = the largest reference ratio, 0.24937654699132977, rounded down) |
| phase_c :: `test_cylinder_3reg_k_against_trajectory_resolvent` | as the sphere row; SN folded 16×32 live; reference (24, 16, 32) | 8e-5 relative = `tolerance_for(3.0e-5, None)`; strict xfail, `raises=AssertionError` | 5.97e-4 relative (k_ref 1.231036749830859, k_sn 1.231772079284479): RED at its tolerance even with no floor | `verify_agreement(sn, Eigenvalue(), ABA_CYLINDER_REFERENCE, 8e-5 × 1.23, NotYet(564, …))`, strict xfail `raises=ReferenceNotValid` |
| phase_c :: `test_cylinder_3reg_flux_shape_against_trajectory_resolvent` | as the sphere shape row | 2e-2; strict xfail | 5.315e-3 (nodal spline); 4.881e-3 through the extension with the solver's own angular rule | `verify_agreement` over the 80 ratios, strict xfail `raises=ReferenceNotValid` (the first call refuses) |
| phase_c :: `test_cylinder_3reg_crosscheck_record` | RECORD of k_ref, k_sn, the k gap and the two shape gaps | band 2e-5 (relative for k, absolute for gaps) | the shape gaps move by 4.3e-4 under the extension | RECORD, every value read through `reference.read` / `answer.read`; `phase_c_shape_fast`/`_thermal` RE-BASELINED to the extension readings |
| standoff :: `test_cylinder_l1_sweep_vs_trajectory_resolvent` | `abs(k_sweep − k_ref)/k_ref`, folded 4×8, 40 cells, source iteration | 2e-3 = `tolerance_for(5.9e-4, None)`; strict xfail | 1.5e-5 | `verify_agreement`, strict xfail `raises=ReferenceNotValid` |
| standoff :: `test_cylinder_l1_refinement_against_reference[20, 40, 80]` | the same for sweep and Krylov at each nx | 2e-3; strict xfail | not measured (the rows never reach the comparison) | as above |
| standoff :: `test_cylinder_l1_reference_record` | RECORD k_ref, k_sweep(nx = 40) | 2e-5 | unchanged | RECORD; k_ref read as `reference.read(Eigenvalue()).value` |
| unified :: `test_unified_cylinder_l1_mr_2g_trajectory_resolvent` | `abs(k_unified − k_ref)/k_ref`, Krylov, folded 4×8, 40 cells | 2e-3; strict xfail | 1.5e-5 | `verify_agreement`, strict xfail `raises=ReferenceNotValid` |
| unified :: `…_trajectory_resolvent_record` | RECORD | 2e-5 | unchanged | RECORD; k_ref through `read` |
| harness :: `test_certificate_reaches_every_branch` | `AgreementCertificate`'s four branches | — | green | **DEAD**: its subject retires; deleted, never repaired. Its laws live on in R7.3/R7.4 (`verify_agreement`) and §1.7b.1 (`compare_uncertified`) |
| harness :: `test_tolerance_rule` | `tolerance_for` cases, then `certify_agreement(…).floor_holds` | — | green | IMPROVED: the floor property `tolerance_for(e, b) ≥ 10 b` asserted directly, exactly (`Fraction`) |
| harness :: `test_ladder_error_estimates`, `test_record_moves_red` | the ladder rules; `assert_record` | — | green | REUSED (the record helper moves with the ABA problem, below) |
| harness :: `test_cylinder_bound_rows_carry_the_shared_strict_xfail` | the shared mark: strict, `raises=AssertionError`; the 5 carriers | — | green | IMPROVED: `raises is ReferenceNotValid`; the carrier census unchanged (5 functions) |

**No relative tolerance in `compare_uncertified`**, so each relative tolerance τ_rel becomes τ_rel × K with K the recorded reference value TRUNCATED to three significant figures (1.38 sphere, 1.23 cylinder). Truncation only tightens: τ_rel × K ≤ τ_rel × k for both of today's divisors (k_sn and k_ref). Dividing the error by production's own reading in the test would re-spell the comparison verb test-side.

**Rows NOT re-posed, and why.** No row is given a new `compare_uncertified` that it does not hold today:
- The cylinder phase-C k row cannot be one: `[M]` probe6, its reading 7.35e-4 absolute is above 9.84e-5, so `compare_uncertified` at its tolerance is red, and any tolerance it passes widens the row.
- The cylinder 4×8 rows (standoff, unified) would pass at 2.46e-3 (`[M]` 1.5e-5), but that agreement is two errors cancelling: the 4×8 SN solve is 5.9e-4 from its own angular limit, and the reference is about as far from the converged SN answer (the 16×32 gap, 5.97e-4). A green made of cancelling errors is not evidence (X1), and today these rows hold no live comparison either (strict xfail plus RECORD).
- Gate 4.1, Gate 4.2 and the unified homogeneous row: not in the migration set (they never called `certify_agreement`); see the edge-row note below.

**The `verifies` markers** (orchestrator's ruling of 2026-10-03, relaying ruling 3). `verifies("sn-curvilinear-trajectory-resolvent-crosscheck")` is DROPPED from all 7 migrated physics functions: an uncertified comparison is never a verification claim, and a strict xfail verifies nothing (X3). The label's verifier count in `docs/theory/verification/matrix.rst` drops from 9 rows to 3 `[M]` (matrix line 802). The rows return when #566 (and #516 for the cylinder) give the family a derived bound in P4 and the strict xfails XPASS.

**The 3 surviving verifiers verify the homogeneous reduction only** (rider b). They are Gate 4.2's edge rows (`test_phase_d_trajectory_resolvent_crosscheck[sphere_2g_homogeneous_dd_n20, cyl_1g_homogeneous_folded_4x8_dd_n20, cyl_1g_homogeneous_folded_2x4_dd_n20]`): the SN snapshot's k against the trajectory resolvent's k on a uniform reflective medium, where both equal k∞ by the V_α1 / V_α1_cyl identities. A uniform medium's flux is flat, which nulls every spatial and angular redistribution term; `[M]` (the module docstring) the cylinder rows' flux moves 1.1e-10 under the closure mutation `tau := 0.7`, against 8.8e-2 on the 3-region cylinder; two of the three are 1-group (anti-pattern #3). So they verify the label's homogeneous reduction (k = k∞), not its claim. The label's own equation is stale besides: it states `‖φ_SN − φ_traj‖∞ ≤ 5e-4 on the 5 P0 curvilinear snapshots` (`curvilinear_one_group.rst:5278`), a flux-shape bound no current row asserts. **Proposed marking** (for the archivist, no new marker kind): keep the 3 edge rows' markers; beside the label on the theory page, state that its spatial and angular content is unverified until P4 certifies the family (#566, #516), that the 3 rows verify the homogeneous reduction only, and correct the stale equation (the 2026-09-26 rows compare fission-gauged cell averages and the eigenvalue at derived tolerances, now uncertified).

#### The trajectory-resolvent references as `ReferenceSolution`s (question 2)

**The specifications exist; no ruling is needed.** `[M]` probe1: `GeometrySpecification(Materials({0: A, 1: B}), StructuredGeometry.from_thicknesses(SPHERICAL | CYLINDRICAL, (0.5, 1.0, 0.5), (0, 1, 0), (BC.reflective,)), Eigen(CellCoefficient.every(Channel.FISSION_EMISSION)))` constructs for both coordinate systems (3 regions, 2 groups). It is the ONE definition of the A|B|A problem (today: `ABA_RADII` plus `aba_xs_2g` in the retiring helper, the snapshot generator's `_sphere_3region`/`_cylinder_3region`, and the standoff/unified files' own `_make_2g_mixture` spellings).

**The derivation** (proposed home `orpheus/derivations/continuous/trajectory_resolvent/reference.py`; `derivations` may import `reference`):
- `trajectory_resolvent_reference(specification, quadrature, *, max_iter, tol, initial_k) -> ReferenceSolution(specification, TrajectoryResolventDerivation(…), certificate=None)`. It routes through `Billiard(specification.geometry, materials, quadrature).solve_critical(...)`, which already maps a `StructuredGeometry` onto the multi-region sphere and cylinder solvers and forwards `max_iter`, `tol`, `initial_k` (P1 step 2b). It refuses an `InfiniteMediumSpecification`, a question that is not the k question, and every body `Billiard` refuses.
- **Lazy.** The solve runs on the first `evaluate`, once (cached on the derivation). Construction solves nothing, so the verbs' refusal (which precedes every read, `[M]` probe6: 0 evaluations under `verify_agreement` and `verify_order`) costs nothing.
- **The state is the emission density** q_g at the radial nodes, `(σ_sᵀ φ + χ νΣ_f·φ / k)/(4π)` per region, never the nodal φ (ruling 1: one reference, one answer).
- **`evaluate` returns `Uncertified` for every observable**: `Eigenvalue` reads k; `FluxIntegral(w)` reads Σ_g ∫ w_g φ_g^ext dV in the geometry's measure (4πr² dr, 2πr dr), by float quadrature per region SPLIT at the weight's own breakpoints (a `RegionwiseConstant`'s region edges; a `Symbolic` `Piecewise`'s condition boundaries, refused when they cannot be extracted); `PointValue` reads φ_g^ext(r). φ^ext is the transport integral of q along rays, with the angular integral split at every tangency of an interface or of an emission-density knot (the reading-bound prototype's split set; ruling 2: ordinary float quadrature).
- **One transport, not two.** The ray transport is the chord oracle's own body, carved so that the emission-density interpolant and the evaluation points are two arguments (today `r_nodes` is both the spline's knots and the evaluation points, `MultiRegionSphereChordOracle.apply_operator`). The solver calls it at the nodes, the reading anywhere. A second transport written for the reading would be a twin path (X4); R7b2.5 asserts there is one.

**`[M]` the reading against today's nodal reading** (probe3, probe5; sphere at (36, 96), cylinder at (24, 16, 32), 8-point Gauss–Legendre per SN cell in the geometry's measure). E0 is today's nodal spline, E1 the extension with the solver's own angular rule, E2 the transport integral (μ split at every tangency, 16 points per piece). Each value is the shape metric's normalisation, `max |Δ| / max ŝ_ref`.

| | sphere | cylinder |
|---|---|---|
| E1 − E0 | 3.2e-4 | 6.4e-3 |
| E2 − E0 | 3.4e-4 | not measured |
| shape gap, E0 → extension | 4.361e-3 → 4.325e-3 (E2) | 5.315e-3 → 4.881e-3 (E1) |
| the reading's own quadrature error: 8 → 16 points per cell | 6.0e-7 | — |
| the same: 16 → 32 per μ piece | 1.2e-8 | — |
| cost | 44 s for 320 points, 2 groups (E2) | 37 s (E1); E2 `[R]` 30 to 60 min (azimuthal splits at every tangency × axial nodes) |

The sphere's shape row passes at 2e-2 under both readings. The cylinder RECORD's shape keys move by 4.3e-4, twenty times the band, so they are re-baselined. The re-baseline is declared in the commit, with its cause (the reading changed from the nodal answer to the extension), never absorbed. **Declared cost:** the cylinder E2 reading enters the slow suite through its RECORD row. If it is above about 10 minutes, the RECORD drops its shape keys and keeps only the k keys; this is a sizing decision for the main agent, not a reason to read the nodal answer.

The per-cell ratio form reproduces today's metric exactly: `[M]` probe5, `max_{g,i} |m − v| / M` = 4.325303743809688e-3, equal as floats to `_shape_gap` on the same readings.

#### The production reading (question 3)

**Home: `Solution.read(observable) -> Measured` on the forward SN solution** (`orpheus/sn/solution.py`), the twin of `HomogeneousResult.read` (G1: the answer is the receiver of `read`; `ProductionAnswer` is the Protocol both satisfy):
- `Eigenvalue` reads `outcome.keff` on an `EigenOutcome`, and is refused on a `SourceOutcome`.
- `FluxIntegral(w)` reads Σ_g Σ_i φ_{g,i} ∫_{cell i} w_g dV. This is exact for the cell-average answer (piecewise constant), which is "today's mesh-bound `Functional` one level down" of G1: an `InnerProductFunctional` over the weight's cell integrals.
- `Ratio` reads the quotient of its operands' readings.
- `PointValue` is refused: a cell-average answer has no point evaluation.

**What the tree lacks:**
1. **The weight's cell co-vector**, ∫_{cell} w_g dV: a mesh fact (the cells and the coordinate measure), so its home is the mesh tier (`orpheus/mesh`, which may import `numerics`). For a `Symbolic` weight it is exact by SymPy integration of each group's expression in r times the measure, after `without(mu, phi)`. That reuses the one definition of "independent of the direction" and refuses a direction-dependent weight.
2. **A region map for a `RegionwiseConstant` weight.** `Mesh1D` holds `coord`, `edges`, `volumes`, `mat_ids` and the face laws, and no region index. Region is not material (the A|B|A regions 0 and 2 share material 0), so reading a region table through `mat_ids` silently merges regions. The `Mesher` partitions interval by interval and DROPS the region index, a lossy return (`coding-elegance` Pattern 4, the lossy-return corollary). Until the region index is kept or P4 projects the specification's geometry, the SN `read` refuses a `RegionwiseConstant` weight as a scope boundary naming that machinery. Part 2 needs only `Symbolic` weights: the cell indicators, and the gauge νΣ_f(r) piecewise in r. The gauge is the hand-resolved form of step 8's `Rate(fission production)`, declared as such.
3. **The answer–specification pairing**: P4's, the same declared limit as `HomogeneousResult.read`.

7b's needs fall short of P4: only item 1 is built. **Not a test-side adapter**: an adapter would be a second definition of the reading that every future SN verification re-spells (X4); the reading is the answer's verb.

**Review note (one definition).** "A ratio reads as the quotient of its operands' readings" would then be spelled three times: in `HomogeneousResult.read`, in `Solution.read` and in `ReferenceSolution._read`. Hoist it to one place at this step.

#### Admissibility of the shape observables (question 4)

`[M]` probe1, on both specifications: `Eigenvalue()`, the per-cell numerator `FluxIntegral(Symbolic.of(Piecewise((1/V_i, (r ≥ a_i) & (r < b_i)), (0, True)), 0))`, the denominator (a `RegionwiseConstant` or `Symbolic` νΣ_f) and their `Ratio` are all admitted. The weights: numerator, group g's indicator of SN cell i divided by V_i (in the geometry's measure, 4π/3 (b³ − a³) or π (b² − a²); `[M]` probe2 equal to `mesh.volumes` to 1 ulp), other group 0. Denominator, νΣ_{f,g}(r) piecewise on the breakpoints 0.5, 1.5, 2.0.

**A defect found** (`[M]` probe4, 4 of 4 shapes admitted on both coordinate systems): `admit_observable` admits a `FluxIntegral` whose weight depends on μ (and, on the cylinder, on the azimuth). A flux integral is a functional of the SCALAR flux, which has no direction, so such a weight has no meaning, and each reader would have to decide one. It is fixed in this step, in `_admit_datum`'s weight role (R7b2.7), not by each reader. Also admitted, and harmless: a weight supported off the geometry, a zero weight. A `Ratio` whose denominator weight is identically zero is admitted and then divides by zero at read time; it is noted, not gated.

#### The cylinder's strict xfail (question 5)

`[M]` probe6, against part 1's working-tree verbs: `verify_agreement` and `verify_order` on a reference with no certificate raise `ReferenceNotValid("verification: the reference has no certificate, so it cannot stand as a reference")` before any evaluation (0 calls). The re-keyed shared mark is `xfail(strict=True, raises=ReferenceNotValid, reason="#566/#516: the cylinder trajectory-resolvent family derives no bound, so its reference has no certificate and cannot anchor a verification")`. It is narrower than today's `raises=AssertionError`, which any bare `assert` before the floor would also have satisfied (mode 8(4)). Under `--runxfail` each of the 7 cases fails with that message. An XPASS needs the family's factory to return a `Valid` certificate (P4), not a test edit, so the flip is not a no-op (mode 8(10)).

#### The gates

Kinds: T = THEOREM, F = REFERENCE (structurally independent), R = RECORD; "route" = a path assertion.

| id | gate | kind | first red on today's tree | must redden under |
|---|---|---|---|---|
| R7b2.1 | The A|B|A specifications: built once; the SN fixtures' geometries (both 3-region snapshots, the standoff and unified meshes) have the specification's breakpoints, material map and boundary law; the three material spellings (`get_mixture`, `aba_xs_2g`, `_make_2g_mixture`) carry identical transport data (σ_t, σ_s, νΣ_f, χ: `array_equal`), so the rows pose one problem (X4; no gate asserts it today) | T | the helper does not exist | a σ_s transpose in one spelling; a breakpoint moved |
| R7b2.2 | The factory: a `ReferenceSolution` with `certificate=None`; refuses an infinite medium, a fixed-source question, an eigen question along scattering emission, a layered slab; constructing it solves nothing (spy on the two MR solvers: 0 calls) | T, route | `ImportError` | the solve made eager; a refusal arm dropped |
| R7b2.3 | Every evaluation is `Uncertified` (eigenvalue, a `Symbolic` and a `RegionwiseConstant` flux integral, a point value), and so is a ratio; `verify_agreement` and `verify_order` raise `ReferenceNotValid` "no certificate" with 0 solves; `compare_uncertified` accepts | T, route | `ImportError` | `evaluate` returning an `Exact`; validity checked after reading |
| R7b2.4 | The eigenvalue threads: `read(Eigenvalue()).value` is bit-identical to a direct `solve_greens_function_sphere_mr` / `…_cylinder_mr` call with the same parameters, at a cheap resolution (e.g. (8, 8), (8, 4, 8)) | T (threading) | `ImportError` | k read from another iterate; `initial_k` not forwarded |
| R7b2.5 | ONE transport: the derivation's angular flux at the radial nodes and the solver's angular rule equals the solver's final `psi_g` to 10 × its `tol` (the emission density is re-formed from the final φ, one power step later); its oracle call at the nodes is `array_equal` to the solver's own oracle output for the same density | T, route | `ImportError` | σ_s transposed in q; the 1/k dropped; the 1/(4π) dropped; the spline built on the evaluation points (the naive oracle reuse) |
| R7b2.6 | Quadrature splits at the weight's breakpoints: a cell indicator on (0.6, 0.7), inside region 1 and off every node, reads within 1e-9 relative of the same functional summed over explicit sub-intervals; a `Piecewise` whose breakpoints cannot be extracted is refused | T | `ImportError` | the breakpoints ignored (`[R]` an O(h) error, much larger than 1e-9) |
| R7b2.7 | Admission refuses a flux-integral weight that depends on μ, or on the azimuth on the cylinder; a point value and a ratio are unaffected | T | `[M]` admitted today (probe4) | the refusal dropped |
| R7b2.8 | The cell co-vector (mesh tier): Σ_i ∫ 1 dV equals Σ `volumes` (nulp at depth N) and equals the geometry's volume; the indicator of cell i gives V_i there and 0 elsewhere, exactly; a step inside a cell gives its partial volume in closed form; sphere and cylinder measures differ (r², r) | T | missing | the measures swapped; the edges shifted by one cell |
| R7b2.9 | `Solution.read`: k bit-identical to `outcome.keff`; refused on a source outcome; Σ_i of the per-cell indicator readings equals the whole-domain reading (nulp); the X1 negative, φ scaled by 1 + 1e-11 moves the reading; a ratio is one quotient; a point value and a `RegionwiseConstant` weight refused, the latter naming the region map | T, route | `AttributeError` | ratio reversed; groups reversed (2-group fixture); 1/k returned; the RegionwiseConstant read through `mat_ids` |
| R7b2.10 | The sphere rows, re-posed (table above), 1 + 80 comparisons, slow | F (uncertified) | the helper does not exist | today's documented arms re-run on the new rows: the reflective face realised as vacuum (k 2.6e1), a 1 % perturbation of the reference's moderator density (shape 7.2e-2), the one-spline reference (ERR-090) |
| R7b2.11 | The cylinder rows (5 functions, 7 cases): strict xfail on `ReferenceNotValid`; the shared-mark census of the harness row asserts `raises is ReferenceNotValid` and the 5 carriers | T (census), xfail | `raises is AssertionError` today | the mark made non-strict; `raises` widened to `Exception` |
| R7b2.12 | The RECORD rows read through `read` (one quantity in the record and the comparison); the k keys unchanged, the shape keys re-baselined to the extension | R | the shape keys move by 4.3e-4 (band 2e-5) | the documented arms of each RECORD (`tau := 0.7`, the vacuum-for-reflective law, the one-spline reference, a 1 % reference perturbation) |
| R7b2.13 | Layers: the new derivation module imports `reference` and no L3 package; `orpheus/reference` still imports no `derivations` (R6.7); `sn/solution.py` reads observables from `numerics` only | T | missing module | a `derivations` import inside `reference` |

**The battery, owed after the code lands**: one arm per "must redden" cell, in-process rebinding (a `-p` plugin), and the positive control, the reference's k scaled by 1 + 1e-3, which must redden R7b2.4 and R7b2.10. Rows rest on R7.* (the verbs), §1.7b.1 (`compare_uncertified`), R6.* (`ReferenceSolution`), S6.8 (the mesh-free functions), `test_trajectory_resolvent_regionwise_source.py::test_mr_oracle_first_leg_matches_the_line_integral` (the oracle's first leg against a line integral), and the Gate 4.2 edge rows.

#### What retires, and the three searches (question 6)

Retired in one commit, leaves first (`retirement-audit` G.24):
1. `CYLINDER_3REG_REFERENCE_BOUND`.
2. `AgreementCertificate` and `certify_agreement`.
3. The module `_certified_agreement.py`. Its ABA problem (the specifications, the cached references, the RECORD table, `assert_record`, the shared mark) moves to one helper; proposed `tests/gates/sn/verification/analytical/_aba_reference.py`.
4. `test_certified_agreement.py` renamed to the subject it keeps (the cross-check harness), with its one DEAD row deleted.

`sphere_3reg_reference_bound` is RENAMED (proposed `sphere_3reg_reference_ladder_estimate`), its docstring re-scoped from "bound" to "estimate, the tolerance's provenance". A name that says bound for something no ladder certifies is prose asserting (X3).

The three searches, `[M]` 2026-10-03:
- **Graph.** `callers(certify_agreement)` returns 9 test functions, all in the table.
- **Text** (`git grep -w`, docs build excluded; plus a `grep -r` over `.claude/` and `docs/` that includes untracked files):

  | symbol | tests | `.claude/plans` | other |
  |---|---|---|---|
  | `certify_agreement` | 18 lines in 5 files | 9 lines in 2 plans | — |
  | `AgreementCertificate` | 6 in the helper | 4 | 1 docstring, `tests/gates/reference/test_verification.py:182`, past tense (keeps) |
  | `CYLINDER_3REG_REFERENCE_BOUND` | 11 lines in 4 files | 1 | — |
  | `_certified_agreement` | 10 lines in 5 files, of which the ladder module has 2 | 1 | — |
  | `sphere_3reg_reference_bound` | 3 lines in 2 files | 1 | — |

  0 hits in theory prose for any of these symbols, and 0 in `.claude/agent-memory/` and `.claude/agents/`. The generated `docs/theory/verification/matrix.rst:55` (the harness file's row) and `docs/_generated/vv_audit.json` regenerate on the next build.
- **Direct constructors.** `AgreementCertificate(` appears at 1 site, inside `certify_agreement`.
- **The concept, not the symbol** (B.7). `docs/theory/methods/sn/curvilinear_numerics.rst:1312-1316` calls the sphere shape row "a plain L1 test" and the cylinder half a strict xfail "on #516". After the migration the first is an uncertified comparison and the second's expected failure is the refusal. The ERR-090 entry (`error_catalog.rst` near line 8592) names these rows. Both go to the archivist with the label note above.
- **After the delete.** Run `dead_references`; and `-n` as a set difference over the touched pages.

**Surviving gates** (`test-architect` §4), the ones the carve leaves untouched:
- The twin-path rows, slab rows, Gate 4.1 and Gate 4.2 keep their claims; they only re-import `aba_xs_2g`/`ABA_RADII` from the new helper.
- None is DEMOTED, PROMOTED or INVERTED. A strict xfail on the refusal is an expected refusal, not degradation pinned as the contract: it XPASSes on certification.
- The one DEAD row is listed above.

#### Refuted candidates, each with its structural reason

- **The sphere kept `Valid` on its ladder bound.** Refuted by the step-5 ruling: no ladder certifies. The sphere is uncertified like the cylinder.
- **The nodal φ through a per-region spline (E0) as the reference's reading.** Refuted by ruling 1: it is the second, nodal answer. FACT: it differs from the extension by 3.4e-4 (sphere) and 6.4e-3 (cylinder) in the shape metric.
- **The Nyström interpolant with the solver's own angular rule (E1) as the reading.** Refuted FOR the reading question: at the nodes it equals the nodal φ, so it is the nodal answer extended, and it carries #516's angular kink error. FACT: E1 is the transport integral's discretisation, the thing R7b2.5 uses to pin the threading.
- **`compare_uncertified` for the cylinder phase-C k.** Red at 8e-5 (`[M]` 7.35e-4 absolute); any passing tolerance widens it.
- **`compare_uncertified` for the cylinder 4×8 rows.** Green only by two errors of about 6e-4 cancelling. No live comparison exists there today.
- **A test-side SN adapter.** Refuted: the reading is the answer's verb (G1), and an adapter is a second definition (X4).
- **Reading a `RegionwiseConstant` weight through `mat_ids`.** Refuted: region is not material (regions 0 and 2 of A|B|A share material 0).
- **Gate 4.1 and the unified homogeneous row as `verify_agreement` against the exact infinite medium.** Refuted for 7b: the SN answer answers a `GeometrySpecification` and the exact reference an `InfiniteMediumSpecification`. Equating them is a reduction theorem (a closed homogeneous reflective body's k is k∞), a cross-specification pairing that nothing can check until P4's projection. FACT: those rows are exact-reference comparisons in all but the pairing, the first candidates when P4 lands.
- **The shape row kept as a `max` over cells.** Not a linear observable (§4). FACT: the per-cell ratios reproduce the max exactly (probe5).
- **A relative tolerance obtained by dividing by production's reading in the test.** Refuted: it re-spells the comparison verb test-side. The absolute τ_rel × K, with K truncated, is exact and only tightens.

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
