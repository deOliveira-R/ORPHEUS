# #405 step 2, phase P3: the traced memo (verification specification)

Written by the test-architect on 2026-10-04 for workflow W3 (surgical carve), phase P3, against `docs/p3-traced-memo-plan` at `0e8740ef`. The design it specifies is `.claude/plans/reference_cache.md` from "P3 opened (2026-10-04)" to the end: the user's four rulings (architecture B, one fresh process per miss; the reading is an entry and its solve a child entry; the three clients; docstrings stripped) and the census findings (F1, the solve child is a memo at the solver functions, the signature binding). The form follows `.claude/plans/reference_p2_spec.md`: per step, a table of gates with their kind, their first red on the tree as it stands, the mutation that must redden each, and the rungs each rests on.

**Gate ids are `M<step>.<n>`** (M for the memo). `S<step>.<n>` ids are P1's (`reference_p1_spec.md`) and `R<step>.<n>` ids are P2's, and both are cited by later gates, so a third family avoids two meanings for one spelling. Test functions are named `test_m<step>_<n>_…`.

**Every gate is written, run and mutated** (§2), against a gate-harness PROTOTYPE written only so that the gates could run: `scratch/reference_architecture/p3/ta/wt/orpheus/numerics/traced_memo.py` (about 790 lines) and the client wiring in four files of the same worktree. The prototype is not a design and not production code: the main agent writes the module, and the gates reach it through one adapter file, `tests/gates/numerics/_traced_memo_api.py`, so a different spelling at the build is one edit there. Running the gates surfaced four defects of the prototype's own first draft (findings A, P, E, N) and three gaps in the ruled design (K, D, W), each now a gate (§5).

## 0. What this specification rests on

### 0.1 The design, as the gates read it

One primitive: a pure function of content-identified arguments, memoised on disk under `.cache/references/` (gitignored), decorated `@traced_memo`.

- **The lookup key** is the blake2b-256 of `encode((function id, bound arguments, argument types, platform tag))`. The function id is `module:qualname`. The arguments are bound by `inspect.signature(f).bind(...).apply_defaults()`, then each parameter the function DECLARES canonical is put in its canonical form (`traced_memo(canonical={"radii": as_float_array, ...})`), and the canonical arguments are what the child receives, not only what the key holds. The argument types are a tree of each top-level argument's type (an array's dtype; a container's element types), so that a value `content_digest` identifies with another, but the function can tell apart, is two entries (§5 finding K).
- **On a miss** the call runs in a FRESH interpreter (`sys.executable`, the parent's `-O` level, the parent's `sys.path`, the parent's working directory, the environment WITHOUT any `ORPHEUS_*` variable), traced by `sys.monitoring` `PY_START` from its first line. The arguments cross by a pickle that rebuilds every dataclass through its constructor (finding F1). The child calls the undecorated function, writes the entry atomically and exits; a nested memo call inside it is looked up like any other (a miss spawns a grandchild) and is recorded in the parent entry as a child reference.
- **The entry** is `entry.json` (the schema, the function id, the key, the manifest, the payload tree, the payload digest) and `payload.npz` (every array of the payload; written even when empty). The manifest is:
  - `functions`: `(path relative to its sys.path root, qualname of the OUTERMOST def, digest)` for every traced first-party code object; the digest is of the def's AST with docstrings stripped and line attributes excluded, so comments, layout and docstrings never move it;
  - `modules`: `(relative path, skeleton digest)` for every first-party module whose body ran; the skeleton is the module's AST with every function body replaced by `pass` and every docstring stripped, so it holds imports, constants, decorators, signatures (with their annotations and defaults), class attributes and dataclass fields;
  - `distributions`: `(name, version)` for every site-packages file that ran, mapped through `packages_distributions()` and, for files no top-level name maps (the editable install's import finder), through the distributions' `RECORD` files; a file of an EDITABLE distribution is hashed as source, since its version never moves;
  - `python`: `sys.version` and the cache tag; the platform tag is in the key;
  - `children`: `(function id, key, payload digest)` of every memo the generation called, hit or miss;
  - recommended, not yet ruled (§4 Q4): `data`, `(path, bytes digest)` of every file the generation opened for reading outside the stdlib and site-packages.
- **Validation** re-hashes the manifest against the files found on the CURRENT `sys.path`, imports and runs nothing, and recurses into the children (each must be a hit whose payload digest is the recorded one). A stale or corrupted entry is a miss, never served, and is overwritten by the regeneration.
- **The payload** is JSON with every float written as `float.hex`, numpy scalars with their dtype, arrays in the `.npz` by name, frozen dataclasses by their constructor fields (rebuilt through the constructor, so their laws re-run), and only a dataclass the function's return annotation names may be named by a payload; never pickle. Every load hands out fresh arrays (no memory shared with another load), read-only.
- **The bypass** `with bypass():` makes every memoised call run in the calling process, reading and writing nothing. A test that monkeypatches anything a generation runs reads under it.
- **The verdict** of a lookup is a closed sum `Hit | Absent | Stale(reasons) | Corrupt(reason)`, exposed by `memo.lookup(*args, **kwargs)` without generating. The gates read the verdict and the reasons, so "hit" and "miss" are observed, never inferred from timing (X3).

### 0.2 The API the gates assume (the adapter is the one edit)

`orpheus.numerics.traced_memo` (the recommended home, §4 Q1): `traced_memo`, `TracedMemo` (`__call__`, `lookup`, `key`, `__wrapped__`, `function_id`), `Hit`, `Absent`, `Stale`, `Corrupt`, `cache_root(path)`, `bypass()`, `default_root()`, `trace_call(f, *args, **kwargs) -> (value, Manifest)`, `Manifest` (the five fields above), `validate(manifest) -> tuple[str, ...]`, `function_digest(source, qualname)`, `skeleton_digest(source)`, the file classifier (adapter `classify`), `Unpinnable`, `Unencodable`, `encode_payload`, `decode_payload`, `platform_tag`. The clients (adapter `CLIENTS`): the two solvers at their current names; `trajectory_resolvent_reading` in `trajectory_resolvent/reference.py`; `exact_infinite_medium_reading` in `common/exact_homogeneous.py`.

### 0.3 The test-side instruments

- `_traced_memo_synthetic.Package`: a fresh generator package under the test's `tmp_path`, with a unique name, its root at the front of `sys.path` for that test only; every mutation witness edits this COPY. Each synthetic generator appends a line to a counter file per run, so a gate counts the generations that happened on the generator's side (X3), never from the memo's own report.
- `_traced_memo_api.SpawnCounter`: counts `subprocess.Popen` of `sys.executable` from the test process, the route of a generation (`vv-principles` #26), for the real clients, which carry no counter.
- The real-tree witnesses (M4.4) copy every `*.py` under `orpheus/` (383 files, 7.8 MB; the 810 MB of `orpheus/data` is not needed by the three clients) into `tmp_path` and drive a fresh interpreter whose `sys.path` starts at the copy (it asserts `orpheus.__file__` lies in the copy).

## 1. The step order, the sizes, and the gates

Sizes are sizing only, never a reason (`plan-authoring` EFFORT-IS-SIZING). Rows are collected test rows `[M]`; production lines are the prototype's, an `[R]` estimate of the real module.

| step | content | rows | depends on |
|---|---|---|---|
| 1 | The manifest: the normalised def digest, the module skeleton, the outermost-def map (with Python 3.14's generated code objects), the file classification, `trace_call` | 22 (`test_traced_memo_manifest.py`) | — |
| 2 | Validation against the checkout, and the payload codec | 18 (`test_traced_memo_validation.py`) | 1 |
| 3 | The key, the process boundary and the store: the memo end to end, the bypass, recursion, corruption, the environment scrub; and the data-file rider | 28 (`test_traced_memo_process.py`) + 1 (`test_traced_memo_data.py`, §4 Q4) | 1, 2 |
| 4 | The clients: 4a the two multi-region solvers memoised (the solve child); 4b the trajectory reading; 4c the exact reading; 4d the re-poses of the existing tests that patch a generator | 13 (`test_traced_memo_clients.py`) + the acceptance run (M4.7) + 5 re-posed rows | 3 |
| 5 | Ambient state: no environment read, no global mpmath precision write, the slab routing an explicit setting | 14 (`tests/gates/test_ambient_state.py`) | — (independent; RECOMMENDED before 4, §4 Q6) |

The order puts the leaves first: what a manifest holds (1), how it is checked (2), then the process that writes it (3), then the clients (4). Step 5 is independent of 1 to 4 and has real reds today; its env read is a hazard only for future clients (0 of 4 census workloads import `peierls_nystrom/cases.py`), so the brief's order holds, but nothing prevents landing it first.

**Measured first reds of steps 1–4** `[M]` 2026-10-04, `scratch/reference_architecture/p3/ta/first_red_pristine.log`: the five files in a pristine worktree of `0e8740ef` (`ta/wt0`) collect 82 rows and all 82 are red, 33 failed and 49 errors, every one `ModuleNotFoundError: No module named 'orpheus.numerics.traced_memo'` (78 direct, 4 inside the M4.4 driver). Every row of steps 1–4 is therefore DEFINING; its tooth is its battery arm (§2). Three rows also have a measured red on today's tree independent of the module's absence, named in their row (M3.4, M3.9, M4.3). Step 5's reds are real: 4 of 14 rows red on the pristine tree, 10 green (the census controls).

### 1.1 Step 1: the manifest

`[M]` 22 of 22 green on the prototype. `npx pyright` on all eight gate files (the five test files, the two helpers, `test_ambient_state.py`): 0 errors, 0 warnings.

| id | gate | kind | tooth (arm, §2) | rests on |
|---|---|---|---|---|
| M1.1 | **The def digest moves exactly with what can change an answer.** 15 edits of one source, each row asserting the def digest AND the skeleton digest: a docstring, a comment, blank lines above, a module docstring, a reformat across lines: neither moves; a body operator, a body constant: the def's moves, the skeleton's does not; a signature default, a decorator: both move; a neighbour's body, a method body: neither moves `target`'s | THEOREM (the normalisation) | A1 (docstrings kept) → `[docstring]`; A2 (line attributes kept) → the layout rows | — |
| M1.2 | **The skeleton** (the same 15 rows): a module constant, an import, a def added, a class attribute, a signature, a decorator move it; any function body does not | THEOREM | A3 (bodies kept) → every body row | M1.1 |
| M1.3a | A method is found by `Class.method`; a missing qualname is `None`, never "unchanged"; a second `def` of one name shadows the first, as Python binds it | THEOREM | P0 (indirect) | M1.1 |
| M1.3b | Nested code (an inner def, a lambda, a generator expression) is pinned by its OUTERMOST def; no `<locals>` qualname reaches the manifest | THEOREM | A5 (the code's own qualname) | M1.3a |
| M1.3c | **Python 3.14's generated code objects.** Evaluating a def's lazy annotations runs an `__annotate__` code object reported at the def's line; a PEP 695 generic runs `<generic parameters of …>`. Neither is refused, and neither pins the def's BODY: `annotated` (annotations read, body never run) is absent, while the skeleton, which keeps every signature, is present | THEOREM | A8 (the signature scope pinned by the body) → this row and M3.12a | M1.3b |
| M1.4 | **The manifest is what ran:** the called defs present, `unused` and `Box.other` absent (the negative leg a whole-module hash fails), the module skeleton present, numpy pinned at its installed version, no stdlib file and no `<string>` code among the functions | THEOREM | P0 → red | M1.3b |
| M1.5 | **Every shape of traced file is classified:** a numpy file is the `numpy` distribution, a stdlib file the interpreter, `<frozen …>` and `<string>` dropped, a source file hashed, a site-packages file of no distribution `Unpinnable` | THEOREM (defining refusal) | A7 (the orphan hashed as source) | — |
| M1.6 | Code compiled under a real file's name at a line where the file holds no such def is `Unpinnable` (validation could never re-hash it) | THEOREM (defining refusal) | not isolated (§2) | M1.3a |
| M1.7 | The manifest names the interpreter (`sys.version`, cache tag); the platform tag names the OS, the `-O` level and the BLAS | THEOREM | B4 (via M3.7) | — |

### 1.2 Step 2: validation, and the payload codec

`[M]` 18 of 18 green on the prototype.

| id | gate | kind | tooth | rests on |
|---|---|---|---|---|
| M2.1 | **The X1 witnesses at the manifest** (10 rows, each an edit of the synthetic copy, then `validate`): a traced body, a traced method, a lambda nested in a traced def, a constant the traced code reads, a constant it does not read (a DECLARED spurious miss: the skeleton is per module), each a reason naming what changed; an untraced def's body, an untraced method, a docstring, a comment added, a module never imported: valid | THEOREM | A4 (skeletons skipped) → the constant rows; A6 (whole-file digests) → the untraced rows; A3 → untraced rows | M1.1, M1.2, M1.4 |
| M2.2 | A traced def renamed away, and a module file deleted, are reasons | THEOREM | P0 | M2.1 |
| M2.3 | A manifest recording another numpy version, another interpreter, or a distribution no longer installed is stale, the reason naming it (the rows edit the manifest; the installed environment is not mutable from a gate) | THEOREM | not isolated (one comparison) | M1.7 |
| M2.4 | **Validation imports and runs nothing:** in a fresh interpreter, validating the manifest of a module whose import writes a marker file leaves the marker absent and the package unimported; the control imports it and the marker appears | THEOREM | A10 (validation imports each module) | M2.1 |
| M2.5 | An edit made after a validation in the same process is seen by the next one, and restoring the text restores validity (no per-path cache outlives an edit) | THEOREM | A9 (bytes cached per path) | M2.1 |
| M2.6 | **Floats cross the payload bit for bit:** 10 hand values (`-0.0`, the smallest subnormal, the largest double, `1 + ulp`, `1/3`, `0.1`, the smallest normal) and 20 000 random finite bit patterns (seed 20261004), compared by `struct.pack`; an `int` above 2⁵³, `bool`, `None`, `str` keep type and value | THEOREM | A11 (15 significant digits) | — |
| M2.7 | Arrays keep dtype, shape and bytes; each decode hands out a fresh, read-only array sharing no memory with the last; writing raises "read-only"; a numpy scalar keeps its numpy type | THEOREM | A12 (writeable) → also M3.14, M4.7; A13 (shared) | M2.6 |
| M2.8 | A dataclass payload is rebuilt through its constructor (a payload violating `Spec`'s law is refused with `Spec`'s own message), and only a type the function declares it returns may be named (another is refused, never imported) | THEOREM (defining refusal) | A14 (any type), A15 (no constructor) | M2.7 |
| M2.9 | An object array, a set, a list, a complex scalar, a function, a plain object: `Unencodable`, never a lossy write | THEOREM (defining refusal) | A16 (a list as a tuple) | M2.6 |

### 1.3 Step 3: the key, the process boundary and the store

`[M]` 28 of 28 green on the prototype, and M3.17 green once the prototype records the files it opened (§4 Q4).

| id | gate | kind | first red beyond absence | tooth | rests on |
|---|---|---|---|---|---|
| M3.1 | **A miss generates once, a hit never.** `Absent`; the first call: 1 generation, 1 entry; then `Hit`; the second call: 0 generations; both return the direct call's float bit for bit | THEOREM | — | B15 (a hit generates anyway) | M2.* |
| M3.2 | **The X1 witnesses through the memo** (4 rows): an edit to a traced def or to a constant it reads is `Stale`, and the next call generates and returns the NEW answer; an untraced def or a docstring edit is `Hit` and generates nothing | THEOREM | — | A4, A6, A3 | M2.1, M3.1 |
| M3.3 | **The signature binds one call to one key:** an omitted default, the default spelled, a positional by keyword, a keyword-only default: 1 key, 1 generation; another value, another key | THEOREM | the census `[M]` (`probe_legacy_key.py`): explicit default and omission digest differently without the binding | B1 | M3.1 |
| M3.4 | **The key separates what `content_digest` identifies and the function tells apart:** a float32 array and its float64 twin are 2 keys and 2 answers; `solve(8.0)` raises as the direct call does instead of being served `solve(8)` | THEOREM | `[M]` `probes/probe_key_hazards.log`: on `0e8740ef`, `content_digest` is EQUAL for the float64/float32 arrays, the float64/int64 arrays, `8`/`8.0` and `1`/`True`, while `f(x)=x/3` differs between the dtypes (`0.3333333333333333` against `0.33333334`) and `range(8.0)` raises | B2 (types dropped from the key) | M3.3 |
| M3.5 | A declared canonical form is both the key and what the child receives: a list and an array share 1 key and 1 generation, and the generator sees the array (its answer carries +0.5 only for an ndarray) | THEOREM | — | B3 | M3.3 |
| M3.6 | An argument with no content identity is refused naming its parameter (`ContentlessError`, `'x'`), before any generation | THEOREM (defining refusal) | — | not isolated (the encoder's own refusal, S5.x of P1) | M3.3 |
| M3.7 | Two functions given the same arguments are 2 keys; another platform tag is another key | THEOREM | — | B4 (platform dropped) | M3.3 |
| M3.8 | **The child is born clean and imports from the parent's `sys.path`:** the package root is on the test process's `sys.path` only (not `PYTHONPATH`, not the working directory) and the child still imports it (a worktree parent never spawns a main-tree child, lessons L99); the entry's modules hold the generator's skeleton AND `orpheus/numerics/content.py`'s (tracing began before the first first-party import) | THEOREM | — | B5 (the child keeps its default path) | M3.1 |
| M3.9 | **F1: arguments cross through their constructor:** a frozen dataclass argument's `__post_init__` is in the manifest, and an edit to it makes the entry stale | THEOREM | `[M]` census F1: in 2 of 2 traced trajectory workloads, 6 construction functions absent from the child's trace under the default pickle (7 of 7 present in the constructor-routed control) | B6 (default pickle) | M3.1 |
| M3.10 | **A patch in the parent never reaches an entry, and the bypass reads and writes nothing:** (a) with `helper` patched in the test process, the memo returns the HONEST value (the child imports unpatched code); (b) under `bypass()` the decoy is returned although a warm honest entry exists (nothing read), and no entry is written; (c) the control: outside the bypass the warm entry is served | THEOREM | — | B7 (bypass ignored) | M3.1 |
| M3.11 | A raising generator: the caller gets the child's exception type and message; no entry; the next call generates and raises again (a refusal is never cached) | THEOREM | — | B9 (the failure written as an entry) | M3.1 |
| M3.12a | **A parent pins its child by reference:** `outer` calls the memo `inner`; outer's manifest holds 1 child reference and NOT `inner`'s def | THEOREM | — | A8 → `inner` listed through `get_type_hints` (§5 finding A) | M1.3c, M3.1 |
| M3.12b | **Recursive validation** (2 rows): an edit only the child traced makes the parent `Stale` with a reason naming the child, when the child was cold inside the parent's generation and when it was WARM (a hit inside a generating process must record the reference) | THEOREM | — | B11 (children not validated) → both; B10 (references recorded only when generated) → `[warm-child]` only | M3.12a |
| M3.12c | A child regenerated to another answer makes the parent stale through the recorded payload digest (the child itself is a valid hit) | THEOREM | — | B12 (digest not compared) | M3.12b |
| M3.12d | Two calls of one child in one generation: 1 child generation, 1 reference | THEOREM | — | not isolated | M3.12a |
| M3.13 | **A corrupted entry is refused and regenerated, never served** (5 rows): a flipped byte, a truncated or deleted payload, a truncated `entry.json`, a payload value edited in place (the JSON parses; only the digest sees it); each `Corrupt`, then 1 more generation and the honest answer. (b) An edited manifest (a digest replaced) is `Stale` | THEOREM | — | B13 (digest not checked) → `[payload-value-edited]` | M3.1 |
| M3.14 | A served payload is fresh, read-only and bit-identical: miss, hit and direct call give the same bytes; the type is the declared dataclass; writing raises | THEOREM | — | A12 | M2.7, M3.1 |
| M3.15 | Two processes racing on one miss both return one answer and leave a `Hit` | THEOREM (DECLARED WEAK) | — | B14 (entry written before payload, in place): see §2 | M3.1 |
| M3.16 | **The child sees no `ORPHEUS_*` variable:** a generator reading one answers its default under the memo whatever the caller set; another variable passes; under `bypass()` the variable is seen (the control). The runtime backstop of M5.1, and the reason a withdrawn generator can never write an entry (§5 finding W) | THEOREM | — | B16 (the environment inherited) | M3.1 |
| M3.17 | **Data files the generation read** (the rider, §4 Q4): after an edit to `table.json`, which `lookup` read, the lookup is `Stale` and the call returns the new value | THEOREM | `[M]` on the prototype without the hook: `assert 'Hit' == 'Stale'`, the stale hit served | D1 (the record dropped) | M3.1 |

### 1.4 Step 4: the clients

`[M]` 13 of 13 green on the prototype with the wiring of §1.4a–c (89 s; the four M4.4 rows 12 s each, M4.5 26 s).

**4a, the solve child.** `solve_greens_function_sphere_mr` and `solve_greens_function_cylinder_mr` are decorated `@traced_memo(canonical=…)` with the six array parameters (`radii`, `sigma_t`, `sigma_s`, `nu_sigma_f`, `chi`, `initial_psi`) declared canonical as `np.asarray(·, dtype=float)`, which is each solver's own first reading of them (`[M]` `greens_function_cylinder.py`, the body's first statements), so canonicalising changes no answer. The return types (11 and 10 plain-data fields) are payloads as they are.

**4b, the trajectory reading.** `trajectory_resolvent_reading(specification, solver_quadrature, max_iter, tol, initial_k, quadrature, observable) -> Uncertified`, module level, memoised: it constructs `TrajectoryResolventDerivation` from those fields in the generating process and evaluates. `TrajectoryResolventDerivation.evaluate` calls it with its own init fields; the evaluation body moves to a method the memo calls (the prototype's `evaluate_here`). In the parent, `answer`, `scalar_flux` and the test API's `oracle(ref, g)` still solve on first access, through the memoised solver, so they read the solve child's entry rather than solving again.

**4c, the exact reading.** `exact_infinite_medium_reading(mixture, observable) -> Exact | Uncertified`: `ExactInfiniteMediumDerivation` is built from the MIXTURE (a `ContentIdentity`) and derives its `ExactInfiniteMedium`, so the reading is keyed on content; today it holds the medium, which has no content identity (`[M]` census Q2: a list field and `Fraction`s).

| id | gate | kind | first red beyond absence | tooth | rests on |
|---|---|---|---|---|---|
| M4.1 | (a) The four clients are traced memos where every caller reaches them (the module attribute that `api.SOLVERS` patches and `Billiard` calls). (b) **The route:** a sphere reference's first `read` writes a reading entry and a solve entry, starting interpreters (the activation leg: ≥ 1); a second reference over the same content reads with 0 interpreters started | THEOREM | (a) today the four are plain functions | C1 → (a); C5 (the reading not routed) → (b) | M3.1 |
| M4.2 | **One solve entry serves `Billiard` and a direct caller:** after the reference has read, the direct call with the census's spellings (a LIST of radii, the caller's own arrays) is a `Hit`, starts nothing, and returns the reference's k bit for bit | THEOREM | the census `[M]`: same k bit for bit (cylinder, `1.1703350703390714`), two keys without the binding | C3 (no canonical declaration) | M3.3, M3.5 |
| M4.3 | **F1 on the real client:** the reading entry's manifest holds `TrajectoryResolventDerivation.__post_init__`, `Billiard.__post_init__`, `_route`, `_layered_xs_payload`, `_read_isotropically`, `reference_body`, `_SphereRays.per_group`; not `solve_greens_function_sphere_mr`; exactly one child, the sphere solver | THEOREM | the census `[M]`: 0 of 6 construction functions in 2 of 2 traced children of a pickled derivation | C1 → the solver def appears; C5 → no reading entry | M3.9, M3.12a |
| M4.4 | **The X1 witnesses on the real tree** (4 rows, a copy of `orpheus/`, a neutral edit `_ = None` as a def's first statement): `_layered_xs_payload` → reading `Stale`, solve `Hit`; `solve_greens_function_sphere_mr` → solve `Stale` AND reading `Stale` (through its reference); `_CylinderRays.per_group` (never run by a sphere reading) → both `Hit`; a constant in `reference.py`'s skeleton (`_RAYS`) → reading `Stale` | THEOREM | — | B11 → `[the-solver]` reading stays `Hit`; A6 → `[rays…]`; A4 → `[…constant]`; C1 | M3.2, M3.12b, M4.3 |
| M4.5 | **A served reading is the fresh reading bit for bit:** the sphere's k, a point value and the flux integral of νΣf (`float.hex`), and the exact medium's enclosure and `Exact` expression string, memo against `bypass()` | THEOREM | — | A11 (floats rounded) | M2.6, M4.1 |
| M4.6 | **The re-posed spy rows** (4d): every row of `test_trajectory_resolvent_reference.py` that patches a generator reads under `bypass()`: the four `_SolveSpy` functions (`:136`, `:182`, `:218`, `:580`) and the `Symbolic.steps` counting spy (`:373`). With the memo on, the solve and the weight's steps run in a child, where no parent-side spy counts them | THEOREM (re-pose) | `[M]` acceptance run, cold, prototype wired: 3 rows INVERTED (red): `test_r7b2_2_construction_solves_nothing_and_the_solve_is_cached[spherical]`, `[cylindrical]` (`0 == 1`), `test_r7b2_6_the_reading_consults_symbolic_steps` ("never consulted"). The other three functions (`:182` refusals, `:218` verbs refuse before any solve, `:580` the constructor's refusals) assert `spy.calls == 0` and stay GREEN, DEMOTED silently: a regression that solved before refusing would go through the memo into a child and the spy would still read 0 | the bypass removed → the 3 inverted rows red; for the 3 demoted functions, an arm that reads the reference before the refusal (`api.reference(...).read(Eigenvalue())` prepended) must red under the bypass and is green without it (the demotion, measured at the build) | M3.10 |
| M4.7 | **The read-only acceptance run:** the 12 files holding the 21 consumer sites of the census (Q1c) and the A\|B\|A rows run green with the memo on, from a COLD cache root and again WARM (no consumer writes into a served array; a write would raise). The positive control is a row: writing into a served `phi_g` raises "read-only", and all six array fields of a served solve are read-only | THEOREM (acceptance) + its control | — | A12 → the control row | M2.7, M4.1 |
| M4.8 | An unconverged solve (`max_iter=2`) is a cached value (`converged` False, `iterations` 2; 3 consumer sites read them) and its reading a refusal never cached: `RuntimeError` "did not converge" crosses back; the second reading starts 1 interpreter (the reading), not 2; only the solve entry exists | THEOREM | — | B9 | M3.11 |
| M4.9 | The exact reference is built from 1 + G entries (G = 2: 3); a second construction starts nothing; a point value is refused with the derivation's own `ValueError` "no position", raised in the child | THEOREM | — | C4 (exact not routed) | M3.11 |
| M4.10 | The default root is `.cache/references` under the checkout and `git check-ignore` holds it; the control: `orpheus/__init__.py` is not ignored | THEOREM | green today (`.gitignore:11`, `.cache/`); it guards the next edit of `.gitignore` | `.cache/` removed from `.gitignore` | — |

**4d, the existing gates the carve touches** (`retirement-audit` D):
- The spy rows of M4.6: 3 INVERTED and 3 functions DEMOTED (green and blind) at the landing of 4b unless they enter `bypass()`; all re-posed in the same commit (`retirement-audit` D.14: the demotion reads green, which is why it is listed by name).
- `tests/gates/reference/test_verification.py:899` patches `Ratio.quotient`, which runs in `ReferenceSolution.read` in the parent, outside every entry (census Q6): untouched, and still a real gate.
- `tests/gates/reference/test_reference_solution.py:586, 630, 642, 657` call `ref.derivation.evaluate(...)` on the exact reference: they now read through the memo; `:657` expects the "no position" refusal, which M4.9 shows crosses the boundary with its type. Their reading is PROMOTED (they now also exercise the codec); their docstrings need no change.
- `api.SOLVERS` direct callers (`test_trajectory_resolvent_reference.py:244`, `:512`; `test_reference_body.py:288`) and the 16 legacy direct calls: now memoised; every one of them is in M4.7's acceptance run.

### 1.5 Step 5: ambient state

| id | gate | kind | first red `[M]` on `0e8740ef` (`ta/wt0`) | tooth (owed post-landing) | rests on |
|---|---|---|---|---|---|
| M5.1a | The environment census finds each spelling on its line (`os.environ.get`, `os.getenv`, `from os import environ`, `getenv as g`, `import os as _os`) and not `os.path` (6 rows) | THEOREM (the census's controls) | green | — | — |
| M5.1b | **No module under `orpheus/` reads the environment except the withdrawal switch and the memo** (which reads it only to hand it on without `ORPHEUS_*`, M3.16; found by this gate on the prototype, which first reddened it); the known reader (`withdrawal.py`) must be found; ≥ 300 files read | THEOREM (census) | RED: `1 of 383 modules read the environment: peierls_nystrom/cases.py: [88]` | re-introducing the read | M5.1a |
| M5.2a | The precision census finds `mp.mp.dps =`, `mpmath.mp.prec +=`, `setattr(mpmath.mp, 'dps', …)`, and not `with mpmath.workdps(30):` (4 rows) | THEOREM | green | — | — |
| M5.2b | **No module writes `mpmath.mp`'s precision globally** | THEOREM (census) | RED: 2 of 383, `cylinder_derivations.py:404` and `greens_function_slab.py:456` | re-introducing a write | M5.2a |
| M5.3 | With `mp.dps` set to 17, each of the two origins derivations leaves it 17 and still passes its own identity (fresh interpreter) | THEOREM (behaviour) | RED: `derive_bessel_wronskian_identity left mp.dps = 30` (`probes/probe_mp_dps.log`: both leave 30, 0.28 s and 0.10 s) | `workdps` removed | M5.2b |
| M5.4 | Two fresh interpreters, with and without `ORPHEUS_SLAB_VIA_E1=1`, see the same module-level values in `peierls_nystrom/cases.py` | THEOREM (behaviour) | RED: `{"_SLAB_VIA_UNIFIED": true} != {"_SLAB_VIA_UNIFIED": false}` | DECLARED WEAK: removing the flag makes both sides `{}`; M5.1b is the catcher, M5.4 the behavioural complement | M5.1b |

**The existing gates step 5 touches:** `test_peierls_multigroup.py`, `test_default_flag_is_unified` and `test_env_var_forces_native` (around `:565` and `:590`) pin the environment behaviour as the contract: INVERTED by the carve, re-posed onto the explicit setting in the same commit. And `peierls_nystrom/cases.py:231` documents `ORPHEUS_SLAB_VIA_UNIFIED=1`, a variable nothing reads (`[M]` the AST census finds one reader, of `ORPHEUS_SLAB_VIA_E1`): present-tense false today, fixed in the same commit.

## 2. The battery

`scratch/reference_architecture/p3/ta/battery/run_battery.py`: one TEXTUAL mutation per arm of the prototype or its wiring in the worktree (an in-process rebind cannot reach a generating child, which re-imports from disk: `vv-principles` #17's closing clause). Each arm asserts its anchor occurs the expected number of times (`Uninstallable` otherwise), runs the memo files (81 rows in pass 1; 82 in pass 2, with `test_traced_memo_data.py`) under `python -O -m pytest`, records the red set, restores every file from a pristine copy taken before the first arm, and byte-compares; the driver checks integrity first. Results: `battery/results.json`, one log per arm.

`[M]` 2026-10-04. Pass 1: 36 arms, on the prototype before the data-file hook, scope 81 rows (the four memo files), 2 to 4 minutes per arm. Pass 2: 6 arms, on the prototype with the hook, scope 82 rows (the data file added). Each `none` arm read 0 red (81 of 81 and 82 of 82 green). Integrity after both passes: every mutated file byte-equal to its pristine copy.

| arm | mutation | red | target rows reddened (in the arm's own red set) | other rows reddened, and why |
|---|---|---|---|---|
| **P0** (positive control) | the manifest records no def and no module | 22 | **M1.4** | 21 more across M1.3b/c, M2.1 stale rows, M2.2, M2.5, M3.2, M3.8, M3.9, M3.12b/c, M3.13b, M4.3, M4.4: every row that reads a manifest |
| A1 | def digests keep docstrings | 3 | **M1.1 `[docstring]`, M2.1 `[docstring]`, M3.2 `[docstring]`** | — |
| A2 | digests keep line attributes | 9 | **M1.1 `[blank-lines-above]`, `[reformat]`, M2.1 `[comment-added]`** | 6: every edit that shifts a later line (M1.1 body rows, M1.3a, two M4.4 rows): line-sensitivity, expected |
| A3 | the skeleton keeps function bodies | 9 | **M1.1 `[body-operator]`, M2.1 `[untraced-body]`, M3.2 `[untraced-body]`** | 6 more body rows (M1.1, M2.1 `[untraced-method]`, two M4.4) |
| A4 | validation skips the skeletons | 5 | **M2.1 `[constant-read]`, M3.2 `[constant-read]`, M4.4 `[reading-module-constant]`** | M2.1 `[constant-unread]`, M2.2 (the deleted file) |
| A5 | nested code pinned by its own `co_qualname` | 51 | **M1.3b** | CRASH-DOMINATED (`vv` #17(h)): a `<locals>` qualname is no key of the def table, so every generation raises; read as "the target is in the set", not as coverage of the 50 others |
| A6 | every def digest is the whole file's | 19 | **M2.1 `[untraced-body]`, `[untraced-method]`, M3.2 `[untraced-body]`, M4.4 `[rays-of-the-other-geometry]`** | 15 M1.1/M1.3a rows that compare a def's digest across edits elsewhere in the file |
| A7 | an orphan site-packages file hashed as source | 1 | **M1.5** | — |
| A8 | a def's signature scope pinned by its body | 3 | **M1.3c, M3.12a** | M4.3 (the reading's manifest lists the solver through `get_type_hints`: finding A on the real client) |
| A9 | validation caches bytes per path for the process | 13 | **M2.5** | 12 stale-witness rows that edit after a first validation in the same process: the hazard's real reach |
| A10 | validation imports each module it validates | 1 | **M2.4** | — |
| A11 | floats written at 15 significant digits | 3 | **M2.6, M4.5** | M3.4 (`divide`'s float32 answer) |
| A12 | served arrays writeable | 4 | **M2.7, M3.14, M4.7** | M2.8 (the rebuilt dataclass's array) |
| A13 | served arrays share memory between loads | 1 | **M2.7** | — |
| A14 | a payload may name any dataclass | 1 | **M2.8** | — |
| A15 | a dataclass rebuilt without its constructor | 1 | **M2.8** | — |
| A16 | a list encoded as a tuple | 1 | **M2.9** | — |
| A17 (pass 2) | code with no def at its line accepted as module code | 1 | **M1.6** | — |
| A18 (pass 2) | distribution versions not validated | 1 | **M2.3** | — |
| B1 | the key over the call as spelled (no binding at all) | 38 | **M3.3** | CRASH-DOMINATED: the child's argument assembly reads bound names; B1b is the clean arm |
| B1b (pass 2) | bound, defaults NOT applied (the census's two-keys case) | 7 | **M3.3, M4.2** | the four M4.4 rows and M4.8, whose direct-call lookups spell the defaults `Billiard` omits: the census's finding, measured on the real client |
| B2 | the key without the argument types | 1 | **M3.4** | — |
| B3 | the declared canonical form never applied | 7 | **M3.5** | M4.2, the four M4.4 rows, M4.8 (the list-of-radii direct call misses) |
| B4 | the key without the platform | 1 | **M3.7** | — |
| B5 | the child keeps its default `sys.path` | 30 | **M3.8** | every generation of a `tmp_path` package fails to import: the worktree hazard's real reach |
| B6 | arguments pickled by default | 1 | **M3.9** | — |
| B7 | the bypass ignored | 2 | **M3.10** | M3.16 (its bypass control leg) |
| B9 | a raising generation writes an entry | 4 | **M3.11** | M3.4, M4.8, M4.9 (each has a refusal row) |
| B10 | a child reference recorded only when generated | 1 | **M3.12b `[warm-child]`** | — (`[cold-child]` green, as designed: the arm's half) |
| B11 | children never validated | 4 | **M3.12b both, M3.12c, M4.4 `[the-solver]`** | — |
| B12 | a child's payload digest not compared | 1 | **M3.12c** | — |
| B13 | the payload digest not checked | 3 | **M3.13 `[payload-value-edited]`** | `[payload-byte-flipped]`, `[payload-truncated]` (the npz still parses after those) |
| B14 | the entry written in place, entry before payload, 0.3 s apart | 11 | **NONE: M3.15 stays green** | DECLARED BLIND, confirmed: the race row cannot see a non-atomic write (both racers recover through the `Corrupt` path). The 11 reds are the in-place write's own breakage (an existing directory, M3.13's regeneration). The atomic write is ungated (§6) |
| B15 | a hit generates anyway | 9 | **M3.1, M4.1b, M4.2** | 6 rows counting generations |
| B16 | the child inherits `ORPHEUS_*` | 1 | **M3.16** | — |
| B17 (pass 2) | child references not deduplicated | 1 | **M3.12d** | — |
| B18 (pass 2) | a contentless argument keyed by its `repr` | 1 | **M3.6** | — |
| C1 | the sphere solver not memoised | 10 | **M4.1, M4.3, M4.4 `[the-solver]`** | 7 more M4 rows |
| C3 | the solver's canonical arrays not declared | 6 | **M4.2** | the four M4.4 rows, M4.8 |
| C4 | the exact reading not routed through its memo | 1 | **M4.9** | — |
| C5 | the trajectory reading not routed through its memo | 7 | **M4.1b, M4.3** | the four M4.4 rows, M4.8 |
| D1 (pass 2) | the files a generation read not recorded | 1 | **M3.17** | — |

**The verdict:** 42 arms besides the two `none` runs; 41 of 42 redden their target row, and the one that does not (B14) was declared blind before it ran and is listed in §6. Two arms are crash-dominated (A5, B1); B1 has a clean twin (B1b), A5 does not, and M1.3b's only clean witness is therefore P0 (a declared gap: a mutant that pins nested code by the right digest under the wrong name was not constructed). Rows with no arm of their own: M1.3a (P0 only), M1.7 (B4 through M3.7), M2.2 (P0, A4), M4.10 (the `.gitignore` edit is not run: an arm on a tracked file of the main tree).

**Step 5's teeth are its first reds:** on the pristine tree the four behaviour and census rows are red (§1.5), and a prototype of the step (`with mp.workdps(30):` at both sites; the slab flag a constant) turns 14 of 14 green (`[M]` `wt`, 3.6 s). The same prototype reddens `test_peierls_multigroup.py::TestSlabViaUnifiedRoutingInfrastructure::test_env_var_forces_native`, the INVERTED row §1.5 names (`test_default_flag_is_unified` stays green: it asserts the default, which holds).

## 3. The measurement owed

**What the user ruled P3 is for: the development cadence.** The protocol, after the build, on the host `.venv`, serially, nothing else running:

1. **The A|B|A reference rows** (`test_trajectory_resolvent_reference.py`, the five `tests/gates/sn/verification/analytical/` files naming `aba_reference`, `test_unified_matvec_cylinder.py`) and **the multi-region solver tests** (`test_peierls_greens_function_mr.py`, `_cylinder_mr.py`, `_cylinder_mr_xverif.py`, `test_reference_body.py`, `test_peierls_greens_function_garcia2021.py`): each file under `python -O -m pytest --durations=0` three times: COLD (an empty `.cache/references`), WARM (the same root), WARM again. Report per file: wall time cold, warm, the ratio, and the row whose time moved most; the `slow` rows separately (`-m slow`), since the canonical run deselects none of them but CI may.
2. **The same files with the memo bypassed** (today's arithmetic) as the denominator, from the same commit.
3. **The cache size:** `du -sk .cache/references` after the cold run, the entry count per function id, and the largest payload.
4. **The fixed costs:** one interpreter start for a trivial synthetic generator (the child's import closure: `orpheus.numerics` pulls in about 120 first-party modules); one warm validation of the A|B|A reading entry (176 functions, 123 modules on the prototype).
5. Report every number with its command, the commit, and the repeat count; timings as the minimum of the repeats (`instrument-doctrine` X2, draw-stable statistics).

**Measured now, on the prototype** `[M]` 2026-10-04, `probes/`, one run each, concurrent with the baseline run below (so upper bounds):
- The sphere A|B|A at the gates' resolution (`n_r=8, n_mu=8, n_traj_quad=16`): cold `Eigenvalue` 8.98 s (two interpreters: the reading and the solve), cold `PointValue` 1.96 s (one: the reading; the solve was a hit); warm 0.02 s each. A direct call of the solver afterwards: `Hit`, 0.01 s, k `0x1.61b0767fbf7edp+0` equal to the reading's. Three entries, 97 582 bytes (`probes/probe_traj_client_sphere.log`).
- The exact medium (mixture A, 2 groups): cold construction and one read 6.30 s (three interpreters), warm 0.117 s and 0.113 s. **Today's cost of the same is 0.01 s** (census Q1a). The exact client is about 10 times SLOWER warm and about 600 times slower cold: §4 Q5.
- The reading entry's manifest: 176 functions, 123 modules, 6 distributions (`charset-normalizer`, `h5py`, `mpmath`, `numpy`, `scipy`, `setuptools`: what the import closure ran, not only what the answer used; each a spurious miss on an upgrade, never a stale hit).
- Validation cost: reading and hashing the 121 first-party files of the census's `traj_sphere` workload takes 8.5 ms per pass; parsing them once takes 268.5 ms (`probes/probe_key_hazards.log`). The prototype caches the parse by the file's BYTES, so the first validation in a process costs about 0.3 s and every later one about 10 ms.

**Today's baseline** (no memo), `ta/baseline/summary.txt`, `--durations=0` per file: `[M]` 2026-10-04, `0e8740ef`, one run per file, `--durations=0`, run CONCURRENTLY with the battery and the acceptance run (three serial streams on 10 cores): upper bounds, not the protocol's numbers.

| file | today (no memo) | memo, cold (empty root, files in this order) | memo, warm |
|---|---|---|---|
| `test_trajectory_resolvent_reference.py` (46) | 414.7 s | 424.3 s (3 red: M4.6) | 55.4 s (3 red) |
| `test_crosscheck_harness.py` (5) | 2.6 s | 2.8 s | 2.1 s |
| `test_l1_standoff_slab_cylinder.py` (14) | 1949.5 s | 1920.0 s | 1063.0 s |
| `test_aba_specification.py` (7) | 2.2 s | 1.8 s | 1.7 s |
| `test_phase_c_crosscheck.py` (9) | 1149.1 s | 511.5 s | 41.8 s |
| `test_unified_matvec_cylinder.py` (32) | 870.8 s | 55.9 s | 53.6 s |
| `test_peierls_greens_function_mr.py` (5) | 169.7 s | 161.1 s | 1.8 s |
| `test_peierls_greens_function_cylinder_mr.py` (10) | 559.2 s | 566.2 s | 2.1 s |
| `test_peierls_greens_function_cylinder_mr_xverif.py` (1) | 23.6 s | 25.0 s | 2.0 s |
| `test_reference_body.py` (59) | 20.5 s | 14.5 s | 2.3 s |
| `test_peierls_greens_function_garcia2021.py` (17) | 5.8 s (pytest) | 7.0 s | 6.8 s |
| `test_eigenvalue_finalize_reconstruction.py` (95) | 144.7 s (pytest) | 147.0 s | 53.3 s |
| `tests/gates/reference/` (386) | 7.6 s (pytest) | 60.5 s | 14.5 s |
| `tests/gates/sn/test_solution_read.py` (15) | 6.6 s (pytest) | 13.2 s | 8.1 s |

Read with care: "cold" is cold for the FIRST file only; later files of the cold pass meet entries the earlier ones wrote (the A|B|A cylinder solve is shared by the reference rows, `phase_c` and `unified_matvec`, which is why those two are already fast cold). The last four rows' "today" figures are pytest's own time on the pristine worktree (`baseline/extra_pristine.txt`), the others `/usr/bin/time`'s real time. The cache after both passes: 8 724 KB (`du -sk`). **The acceptance (M4.7):** 14 targets, 785 rows, cold and warm: every row green except the three M4.6 predicts, so 0 consumers wrote into a served array (a write raises). **The exact client** (`tests/gates/reference/`, `test_solution_read.py`) is SLOWER under the memo, cold and warm (§4 Q5).

## 4. Open questions for the user, each with a recommendation

1. **The home of the primitive.** Recommended: `orpheus/numerics/traced_memo.py`, beside `content.py`, whose encoder it extends (the key) and whose layer it shares (no physics; `derivations` may import `numerics`). The alternative is a new input-tier package (`orpheus/memo/`), which would need its registration in `tests/gates/test_layer_imports.py` and S5.9's roster walk (lessons L102), for one module.
2. **The key's type tree** (§5 finding K). The ruled design keys on the arguments' content digest; `content_digest` follows `==` by the user's ruling of 2026-10-02 (`8 == 8.0`, an int array equals its float twin), and a function can tell those apart, so a pure content key serves one caller another's answer. Recommended: the key adds the type tree of the top-level arguments (an array's dtype, a container's element types; nothing inside a `ContentIdentity` value, whose constructor already canonicalises), at the price of a spurious miss for `8` against `8.0`; and the six array parameters of the solvers are declared canonical, which is what makes the list-and-array callers share an entry (M4.2). The alternative, normalising every array-like globally, is unsound for a parameter where a list and a tuple differ.
3. **The platform tag.** Recommended: `sys.implementation.cache_tag`, `sys.platform`, `platform.machine()`, the `-O` level, and numpy's BLAS name (Accelerate against OpenBLAS is where #504's last-bit differences come from). The `-O` level is there because a bare `assert` in a generator is stripped under `-O`, so the two levels can answer differently. The thread count of the BLAS is NOT in it: whether it moves the last bits here is unmeasured, and M4.5 (memo against in-process, same machine) is the gate that would show it.
4. **Data files** (M3.17, §5 finding D). The trace records code, not the files a generation reads. Today 0 of the 3 clients read a data file (`[M]` AST census: 4 files under `orpheus/` read data, all in `micro_xs` and the Sood cache), so P3 is sound without it; P4's generators will read nuclear data. Recommended: land the audit hook on `open` in the generating process in P3, with M3.17, so the manifest pins every file a generation opened for reading by its bytes. Its blind spot must be stated where it is built: a C library that opens a file itself (HDF5 through `h5py`) raises no Python audit event, so a generator reading HDF5 must take the file's digest as an argument, a rule P4 inherits.
5. **The exact client.** Measured: warm about 10 times slower than computing (0.11 s against 0.01 s for a construction and a read), cold 6.3 s against 0.01 s; on the suites that use it, `tests/gates/reference/` takes 7.6 s today, 60.5 s cold and 14.5 s warm under the memo, and `test_solution_read.py` 6.6 s, 13.2 s and 8.1 s (§3; concurrent runs, upper bounds). The user ruled it a client for uniformity. Recommended: keep it a client only if uniformity is worth that cost (one primitive for every `ReferenceSolution` producer; the exact reading is then gated by the same M-rows); otherwise exclude it explicitly, with a sentence on its theory page that a reading cheaper than one interpreter start is not memoised, and keep M4.9 as the gate that it is not. My recommendation is to EXCLUDE it: memoising a computation cheaper than its own validation buys nothing the cache exists for, and the uniformity is already provided by the shared `Derivation` protocol.
6. **The order of step 5.** Recommended: land step 5 first, not last. It has the only reds that exist on today's tree, it is independent, and M3.16 (the environment scrub) and M5.1b are two halves of one claim, better landed together.
7. **The withdrawal switch and the memo** (§5 finding W). `ORPHEUS_RUN_WITHDRAWN` decides whether a withdrawn generator runs at all. A memo around a withdrawn generator would, without M3.16, write an entry while the switch is set and serve it to a later caller without it, bypassing the refusal. Recommended: M3.16's scrub (the child never sees the switch, so a withdrawn generator always refuses inside a child) plus a rule P4 inherits: a withdrawn generator is never memoised (a gate in P4: no `@traced_memo` on a function reachable from `withdrawn_generator`).
8. **`ExactInfiniteMediumDerivation`'s field** (4c). Recommended: it holds the mixture and derives the medium (a `cached_property` or a field `init=False`). The alternative keeps the medium and keys the reading on `specification.mixture` passed separately, which gives the derivation two sources for one quantity (X4).
9. **The bypass's spelling.** Recommended: a context manager, `with bypass():`, entered by the five re-posed spy rows (a `pytest` fixture wrapping it is test-side sugar). Not an environment variable: that would be the ambient state step 5 removes.

## 5. Findings of the prototype (each now a gate)

- **A, the signature scope pins a body.** On Python 3.14 (PEP 649) `typing.get_type_hints(f)` runs a generated `__annotate__` code object reported at `f`'s def line. A manifest that maps every code object to the def spanning its line therefore lists `f` whenever anything reads its annotations. The prototype's first run listed the CHILD's def `inner` in the PARENT's manifest, through the payload decoder's `get_type_hints`, so the recursive witness M3.12b passed through the parent's own manifest instead of through its child reference: green for the wrong reason, and caught only by M3.12a's negative leg. The signature scope is pinned by the skeleton, which keeps every signature (M1.3c).
- **P, PEP 695 scopes.** `<generic parameters of find_factor>` (`orpheus/numerics/space.py:1207`) is a code object with no def of its name at its line; the first prototype refused it as `Unpinnable` and every generation failed. 3.12+ generates such code from signatures and class bodies; it is pinned by the AST that holds the signature (M1.3c).
- **E, the editable finder.** Every child runs `.venv/.../site-packages/__editable___orpheus_0_1_0_finder.py`, which no top-level name maps to a distribution; the first prototype refused it and every generation failed. It belongs to the editable `orpheus` distribution through its `RECORD`, and an editable distribution's version never moves, so it is hashed as source. Building a full file-to-distribution index through `importlib.metadata` took 9.4 to 10.9 s per process (24 484 files); searching the raw `RECORD` texts for the one file takes 0.23 s.
- **N, `np.float64` is a `float`.** A codec testing `isinstance(v, float)` before `np.generic` writes a numpy scalar as a Python float, and `k_eff` comes back with another type (M2.7's last assertion; the first prototype failed it).
- **K, the key's identifications** (§4 Q2; M3.4).
- **D, data files** (§4 Q4; M3.17).
- **W, the withdrawal switch** (§4 Q7; M3.16).
- **O, an operating-system file.** With the data record on, every client entry (`[M]` 5 of 5 entries: the sphere reading, its solve, the three exact readings) records one data file, `/System/Library/CoreServices/SystemVersion.plist`, read during the import closure (macOS's version, read by a platform query). An OS update therefore misses every entry. That is sound (an OS update can change libm and the BLAS) and is declared here so it does not read as a defect.
- **The stale docstring** `cases.py:231` (`ORPHEUS_SLAB_VIA_UNIFIED`, read by nothing), §1.5.

## 6. What no gate here can see

- **Native state.** A C extension's global state (a BLAS thread pool, a C library's precision mode) is neither code the trace records nor an argument. M4.5 compares memo and in-process readings on one machine; it cannot see a difference that both sides share.
- **A file a C library opens** (HDF5): no audit event (§4 Q4).
- **The editable install's mapping**: if the finder maps `orpheus` to another checkout, the child imports from wherever the finder points only when `sys.path` does not already hold `orpheus`; M3.8 shows the path wins for a package on `sys.path`, and the main tree is always on it under `pytest`'s `pythonpath`. Not gated beyond that.
- **Correctness of the answer itself.** The memo makes a reading cheap; it certifies nothing about the reading. A wrong reference served from the cache is exactly as wrong as computed (the families' certificates are P4's, #566).
- **The atomic write.** M3.15 cannot see a non-atomic write (arm B14, declared and confirmed blind): two racers both recover through the `Corrupt` path. The write is atomic in the prototype by one `os.replace` of a staged directory; no gate pins that structure.
- **M1.3b's clean witness.** Its only arms are P0 and the crash-dominated A5; a mutant pinning nested code by the right digest under a wrong name was not constructed.
- **A spurious miss is not a defect** and is not gated except where declared (M2.1 `[constant-unread]`; the six distributions of §3): every rule here errs toward a miss.

## 7. Artefacts

All under `scratch/reference_architecture/p3/ta/` (untracked):
- `wt/`: a detached worktree of `0e8740ef` with `.venv` linked, holding the prototype (`orpheus/numerics/traced_memo.py`), the client wiring (`greens_function.py`, `greens_function_cylinder.py`, `trajectory_resolvent/reference.py`, `common/exact_homogeneous.py`), the step-5 prototype (`cylinder_derivations.py`, `greens_function_slab.py`, `peierls_nystrom/cases.py`), and the gate files, which land at the same paths: `tests/gates/numerics/_traced_memo_api.py`, `_traced_memo_synthetic.py`, `test_traced_memo_manifest.py`, `test_traced_memo_validation.py`, `test_traced_memo_process.py`, `test_traced_memo_data.py`, `test_traced_memo_clients.py`, and `tests/gates/test_ambient_state.py`.
- `wt0/`: a pristine worktree of `0e8740ef` with only the gate files copied in: the first reds (`first_red_pristine.log`).
- `battery/`: the driver, `results.json`, one log per arm, `pristine/`.
- `probes/`: `probe_key_hazards.py` (+ log), `probe_mp_dps.py` (+ log), `probe_traj_client.py` (+ sphere log).
- `baseline/`: today's per-file timings (`summary.txt`, one log per file, `extra_pristine.txt`).
- `wt2/` and `acceptance/`: the M4.7 acceptance run, the prototype as before the data hook, cold then warm (`summary.txt`, one log per file and pass).
- `final_green_run.log`: the six gate files on `wt`, 92 passed and the four step-5 reds (before the step-5 prototype); after it, `test_ambient_state.py` 14 of 14.
- `battery/results_first_pass.json` (pass 1), `battery/results.json` (both passes), `battery_driver.log`, `battery_driver_pass2.log`, `pristine_before_data_hook/` (pass 1's pristine copies).
- `plugins/where_plugin.py`: prints at session finish which tree `orpheus` and the memo module resolved from.

Run a gate file: `cd scratch/reference_architecture/p3/ta/wt && .venv/bin/python -O -m pytest -p no:cacheprovider tests/gates/numerics/test_traced_memo_process.py`.
