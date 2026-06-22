# Pyright cluster triage — making the checker signal, not noise (#226)

Measured 2026-06-18 (pyright 1.1.410, CLI, reading the committed `pyrightconfig.json`, over `orpheus/`).

## The two-population finding (the reframe)

Running **CLI** pyright on the exact files that flooded the session with
"Import … could not be resolved":

| Population | What | CLI count | Verdict |
|---|---|---|---|
| **A — the avalanche** | `reportMissingImports`, "X is not defined (unresolved)", cross-tree `SNMesh≠SNMesh` | **0** | pure IDE-langserver artifact |
| **B — real type debt** | every other rule | **691** | genuine findings |

**Population A does not exist in the type checker.** It is an artifact of the
pyright *langserver* (the IDE plugin): it reads config at startup only, cannot
receive config-change notifications, and roots at the main checkout even across
worktrees (see the `regen-pyrightconfig.sh` header + microsoft/pyright#10498,
anthropics/claude-code#27220). The streamed `<new-diagnostics>` blocks come from
that langserver; **the CLI is the oracle, the LSP stream is advisory.** #226 had
been conflating the artifact (A) with the debt (B) into one "734" number.

## Distribution of the 691 (CLI, `orpheus/`)

By rule: `reportArgumentType` 255 · `reportAttributeAccessIssue` 141 ·
`reportUndefinedVariable` 67 · Optional-family 67 (`Subscript` 39 +
`MemberAccess` 14 + `Call` 9 + `Operand` 5) · `reportIndexIssue` 28 ·
`reportCallIssue` 27 · `reportReturnType` 25 · `reportOperatorIssue` 23 ·
`reportAssignmentType` 14 · `reportPossiblyUnbound` 13 · tail …

By subsystem: **derivations 252** (SymPy reference solvers — billiard 77, peierls
30, mms 13) · **production 439** (sn 246, numerics 61, data 59, cp 22, geometry
17, …).

## The principle

Signal iff **(1) read from the correct oracle (CLI), and (2) every remaining
class is actionable.** Kill always-false classes at their root, in priority
order: **config/lifecycle > type stubs > code precision > narrow documented
suppression.** Never blanket-disable a rule that *sometimes* catches a real bug;
never line-`# type: ignore` what a root fix solves wholesale.

## The clusters — fix order = errors-cleared-per-edit, lowest-risk first

The top-3 rules (463 of 691) are **coupled**: one Unknown-typed symbol or one
too-loose annotation cascades into argtype at every call site + attribute errors
downstream. So fix the root, re-measure, and the count collapses in clusters —
this is far fewer than 691 independent edits.

### Cluster 1 — missing `TYPE_CHECKING` / forward-ref imports (MECHANICAL, do first)
~55 of 67 `reportUndefinedVariable`: `TimedFullField` (25), `AngularQuadrature`
(20), `AngularFlux` (3), `AngularSourceSink` (2), `BoundarySourceSink`/`SNMesh`/
`D1ClosedForm`/`FunctionSpace`/`Quadrature1D`/`Optional` (rest), plus `sp` (7) =
missing `import sympy as sp`.
- **Root**: string annotations (`"TimedFullField"`) whose `if TYPE_CHECKING:`
  import is missing/typo'd — REAL (no checker resolves them).
- **Blast radius**: each undefined annotation makes its symbol Unknown →
  cascades into `reportAttributeAccessIssue` (`.bulk`/`.values` on the Unknown
  value) and `reportArgumentType` at call sites. ~10 import-block edits clear
  the 55 *plus* an unmeasured chunk of the 141 + 255.
- **Fix**: add the missing imports. Zero behavior risk. **Re-measure after** —
  the total drops by more than 55.

### Cluster 2 — too-loose `ndarray` / typed-field annotations (SYSTEMATIC, high collapse)
`reportArgumentType` `ndarray`/`ndarray|None` (39) + `.bulk` on ndarray (5) +
`.real`/`.imag` on ndarray (19) + typed-protocol mismatches (`ScalarFluxFn` 7,
`LinearOperator` 6, `EigenvalueSolver` 4).
- **Root**: functions typed `-> np.ndarray` (or params `: np.ndarray`) that
  actually carry `TimedFullField` / typed fields; or arrays typed so loosely
  that `.real`/`.imag`/element ops are Unknown.
- **Blast radius**: one root annotation → every call site (argtype) + every
  attribute access. The Pattern-7-at-the-producer move, for types: fix the
  producer annotation, the consumers collapse.
- **Fix**: tighten the root annotations to the real types.

### Cluster 3 — `Mesh1D | Mesh2D` union not narrowed
`.centers` (10) + `.volumes` (6) + `.areas` (4) on `Mesh1D` (20) + `GeometryType`
arg (5).
- **Root**: dimension-specific attrs accessed on the union (or on `Mesh1D` where
  the value can be `Mesh2D`).
- **Fix**: narrow (`isinstance` / split d=1 vs d=2 handling) or hoist a shared
  Protocol/attr. Real correctness hygiene — the members genuinely differ.

### Cluster 4 — re-measure (Optional-family, Index, real/imag) AFTER 1+2
Optional-family (67) = None-safety gaps → guard at invariant points
(`assert x is not None` where an invariant guarantees it) or make the field
non-Optional. `reportIndexIssue` (28) and the ndarray `.real`/`.imag` likely
shrink once Clusters 1+2 remove the Unknown-typed cascades. **Re-bucket before
fixing** — much of this is downstream.

### Cluster 5 — derivations + SymPy (252, ISOLATE, lower priority)
billiard 77 / peierls 30 / mms 13: `ConvertibleToFloat`/`int`/`float`/`str`
argtype (SymPy `Float`/`Rational`/`Expr` where scalars expected), `.lhs`/`.rhs`
on `Expr`/`Boolean` (12), `sp` undefined (overlaps Cluster 1).
- **Root**: SymPy is poorly typed; pyright fights its expression algebra. This
  is the reference/MMS subsystem, not production.
- **Options**: (a) fix the genuinely-real ones (`import sympy as sp`, real
  scalar-vs-array mismatches); (b) for the irreducible SymPy noise, add a
  **per-directory `executionEnvironment` for `orpheus/derivations/**`** that
  downgrades the SymPy-fighting rules (`reportArgumentType` /
  `reportAttributeAccessIssue` → `"warning"` or off) — a SCOPED relaxation, NOT
  global, justified by "SymPy is untyped". Production stays full-strict.
- **Priority**: AFTER production (439 is the higher-value target).

### Cluster 6 — patterns pyright structurally can't model (NARROW suppression, last)
`reportCallIssue` `__init_subclass__(key=…)` registry, `reportOperatorIssue`
`__add__`/operator algebra (23), `.SO2`/`.OctahedralOh` on
`type[SubgroupOfO3]` enum/registry (12).
- **Root**: pyright can't infer `__init_subclass__` kwargs, dynamically-built
  dunders, or some enum/registry access.
- **Fix**: narrow `# pyright: ignore[ruleName]` on the exact line, WITH a comment
  naming why it is unmodellable. NEVER a file/global rule disable. These are the
  irreducible baseline.

## The gate (keep it signal forever)

1. **Read pyright from the CLI, never the LSP stream.** The `session-health`
   hook verifies the CLI is clean each session; the agent treats the IDE stream
   as advisory.
2. Drive **production (439)** to zero in cluster order 1→2→3→4, re-measuring
   between clusters (the cascades collapse faster than the raw counts suggest).
3. Then **derivations (252)**: fix-real + the scoped relaxed exec-env.
4. **Gate "no new CLI errors"** — a pre-commit / CI step running
   `pyright --outputjson`, failing if the error count exceeds a baseline. The
   baseline file holds only the Cluster-6 irreducibles, each line justified.
5. Stay in the current strictness until zero; tighten after.

## Recovery / re-measure command
`npx --no-install pyright --outputjson orpheus/` → bucket by `.rule` (and by
subsystem). Compare to this snapshot (691 / 0 imports) to track progress.
