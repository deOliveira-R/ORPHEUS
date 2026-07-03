# Pyright signal cleanup (#226) — "only signal, not noise"

**Goal (user, 2026-06-21):** make pyright give the agent *only real, actionable
signal* — (1) **kill the noise**: the streamed `<new-diagnostics>` LSP avalanche
of phantom `reportMissingImports` / "not defined" that floods every turn must
stop or match the CLI; (2) **drive real errors toward zero** so a genuine pyright
hit is always a real bug, never background hum.

**Oracle:** `npx pyright` (CLI). The streamed LSP `<new-diagnostics>` is advisory
and currently NOISY (microsoft/pyright#10498 dual-identity). **Constraint (hard,
user):** NO `# type: ignore`, NO blanket suppression — fixes are principled
(declare dynamically-set attrs statically [runtime-inert]; fix the real type;
correctly parametrise a Protocol). Per-cluster commit + **re-measure** (counts are
coupled — one `Unknown`-typed root cascades into argtype+attribute at every call
site, so raw counts overstate the independent-edit count).

**Branch:** fresh off `main` @ `7aebab3` (e.g. `refactor/pyright-signal-2`; note
an older `refactor/pyright-signal` branch predates the #257 campaign — do NOT
reuse it without rebasing onto current main). ff-merge per cluster-batch.

---

## Progress log

### Session 2026-06-22 — orpheus/ **502 → 470**, all merged + pushed to `main` (`5856160 → 11fa279`)

Ratchet re-baselined 706 (stale) → taut, then per-cluster commits (each: tests green under `-O`, count drops, 0 net-new `# type: ignore`, ratchet re-tightened):

- **Ratchet re-baseline** (`af13734`): the #257 + foundation merges had burned orpheus/ 706→502 without re-tightening — closed a 204-error regression hole. Baseline now tracks the live per-module floor (`tests/_harness/pyright_baseline.json`).
- **B1 warm-up** (`fe90b77`, −16): the 11 `reportUndefinedVariable` were REAL missing-import annotations in `derivations/` (sympy stays lazy → a module-level `TYPE_CHECKING` block; `Quadrature1D` folded into the already-runtime `common.quadrature` import; `TimedFullField`) — fixed regardless of the B6 derivations exec-env decision. Cluster D = 6 `bool_`/`csr`/`csc` return narrowings (`bool(...)`, build `csr_matrix` from triplets directly, `np.asarray(spsolve(...))`).
- **B2 §3 MoC plotters** (`3ed5fec`, −6): `plot_moc_rays`/`plot_moc_mesh` were a retired-Cartesian relic (`geom.n_cells`/`delta`/`mat_map`) that `AttributeError`'d **live** (demo_moc.py calls them). Rewritten against the actual reverse-Wigner-Seitz ring model (`pitch = R_outer·√π` equal-area square, concentric `radii`, real `Track.entry/exit` chords). User reviewed the rendered figures. **#266 filed** for the third plotter `plot_moc_spatial_flux` + the broadly-stale `demo_moc.py` (`default_pwr()` doesn't exist either) — a bigger rewrite + keff re-validation; it is the **last orpheus/(root) error**.
- **B2 §2 AngularMeasure** (`aa93757`, −4): widened the deliberately-lean geometry `AngularMeasure` Protocol with the genuine curvilinear-contract members `eta`/`mu_z`/`level_indices` (types matching the concrete `Quadrature` properties); retired 2 pre-existing `# type: ignore`.
- **B2 §1 PoleAngularClosure retirement — DEFERRED to #236** (user decision, 2026-06-22). The explorer's premise was WRONG, verified against main: the ABC `PoleAngularClosureBase` is `__call__`-only; the strategy methods (`precompute_psi_state`/`cell_contribution`/`angular_adjoint`/`level_indices`/`psi_half_seed`) are `MorelMontryAngularSweep`-specific; `IdentityAngularClosure` (cartesian) can't carry them as abstractmethods; #236's ABC-contract work (`6fdb0fe`) is **unmerged** (merge-status trap). Retiring the orphan would require narrowing consumers to the concrete (= B3 work on loss_representation.py) and would be superseded once #236 hoists the methods to the ABC — so keep the orphan Protocol until #236 lands. **#248 updated.** The independent half (AngularMeasure) shipped in `aa93757`.
- **B3 full_field_space** (`11fa279`, −6): `_require_blocks()` now RETURNS the narrowed `(bulk_space, trace_space)` pair so the three callers bind non-None locals (it was `-> None`, which can't narrow `self.x`). Runtime-inert. The remaining `.bulk`/`.boundary`-on-`NDArray` + the `inner_product` override error in that file are one deeper mismatch — `FunctionSpace` is typed for `NDArray` while `FullFieldSpace` carries composite fields → **B4** generic-hierarchy. LESSON: matching an override param *name* to a base whose *type* is wrong only surfaces the type mismatch harder (the rename was reverted).

### Session 2026-06-22 (cont., post-compaction) — orpheus/ **470 → 440** (`b36ca11 → bff6135`, on branch `refactor/pyright-signal-2`)

Two ROOT-CAUSE structural fixes (both more than local narrowing — they make the type tell the truth), each verified `-O` + ratchet-retightened:

- **`_bases.py` TraceSpace.layout** (`382fe5f`, −3): `BoundaryField.from_face_arrays` read `layout = space.layout` after the `space is None` guard, but `TraceSpace.layout` is `Optional` (bare-constructor footgun). Added an explicit `layout is None` guard (parse-don't-validate, mirrors `_face_row` / `_require_blocks`). #236-independent.
- **`StreamingTerms` required fields** (`dae34e0`, −24 across diamond.py + cell_balance.py): the 10 curvature fields were typed `float | None = None`, contradicting the file's own **Step 2.5 invariant** (all populated for every geometry; slab neutral, curvilinear physical — to retire the `alpha_in is None` branch). Made all 10 required `float` (Pattern 4). **explorer audit** confirmed all 3 production factory arms populate every field by keyword + no production None-consumer. Riders: deleted the now-dead `abs_mu/volume is None` guard in `linear_discontinuous.py`; fixed the self-contradicting docstring (one para said the branch was retired, a later stale one said it was live); migrated 5 test constructors that omitted `mu_start`. diamond/cell_balance/reduced_operator/linear_discontinuous all reach 0.
- **`AngularMeasure` Protocol read-only properties** (`b578286`, −3): the Protocol declared `mu_x`/`weights`/`N`/`eta`/`mu_z`/`level_indices` as mutable attrs → invariant → a concrete read-only `@property` (Quadrature) fails to match ("property is not assignable to ndarray"). Declared them read-only `@property` (correct covariant-read contract; static-only, Protocols never instantiate). Cleared geometry.py:414/420/427 AND the test-suite flood. 509 geometry+sweep tests green.
- Ratchet: `9c6f037` (470→443), `bff6135` (443→440). Tests verified: `tests/sn/sweep/core/` (460), `tests/geometry/test_reduced_operator.py` + sweep/core (509).

**FINDING — #250 is real and pre-existing on main.** Two sphere bit-identity snapshot reds (`test_sphere_{1g,2g}_apply_bit_identical`, 100% mismatch, rel ~3.4) surfaced in the broad run. Verified via a **worktree at clean `b36ca11`** (PYTHONPATH override; editable install pins imports) — they fail IDENTICALLY on untouched main, so NOT caused by this work. Already tracked: **#250** ("re-baseline stale SPHERE snapshots — 5 curvilinear reds from the W1 clamp fix `b2d8a6d`"). Out of pyright scope; a V&V re-baseline (numerics-investigator confirms W1 values, then re-capture `pre_t4_snapshots.npz`). The slab arms already moved to `assert_regression`; the sphere arms still `assert_allclose` the stale npz.

**Two more clusters this session** (`114b679 → 105dbf8`, merged+pushed; ratchet 440→428):
- **registry SelectionLog typed closure** (`457d2e3`, −8): `select_quadrature` splatted a union-valued `dict(...)` (`**log_template`) into two `SelectionLog(...)` constructions (8 errors = 4 fields × 2 sites). Replaced with one typed closure `_selection_log(chosen_spec, chosen_parameters)` — every field passed with its precise type. registry.py 8→0; 82 quadrature tests green.
- **Axis1D Protocol read-only properties** (`77aaf7c`, −4): SAME variance fix as AngularMeasure — `edges`/`coord` were mutable attrs (invariant) rejecting the frozen-dataclass / `@property` implementers; declared read-only `@property` (matching the already-property `n`/`endpoints`/`bc`). axis.py 4→0 + cleared the `coord_system(...)`/`from_axes(...)` test flood. 97 axis+sweep tests green.
- **Protocol-variance sweep is now EXHAUSTED in production** — a precise grep for `is not read-only in protocol` / `invariant because it is mutable` (non-deriv) leaves only **1** straggler (solver.py:753, a different union shape). The two real instances (AngularMeasure, Axis1D) are caught.
- **geometry.py dag_walk narrowing** (`d6c99d9`, −8, ratchet 428→420): the two-entry-mode `dag_walk` (`ordinate_idx` XOR `direction_sign`) + `dag_walk_cell_indices` + `_representative_ordinate` couldn't narrow `mu_level_idx`/`ordinate_idx`. Added `_require_mu_level(mu_level_idx) -> int` (Pattern 2, single-sources the 3 cylindrical guards, replaces the 2 early compound guards) + an `ordinate_idx is None` post-dispatch narrowing. geometry.py 12→4. **Full SN sweep suite green (679 passed)** — hot path, all geometries verified. No test asserted the moved guard messages (checked first).

### Session 2026-07-02 (campaign RESTART, task #139) — scoped ratchet **145 → 132** via the #236 re-stage

Campaign restarted post-#226-merge on `refactor/pyright-burndown` (approved plan:
`~/.claude/plans/reactive-moseying-cake.md`). **User reframe: the burn is an
architecture-smell detector** (the 18 `cast(` sites are flags, 3 families mapped); **Level 1
scope = the minimal production architecture** (data/geometry/numerics/transport/sn/homogeneous);
thermal_hydraulics+kinetics NEVER lifted (possible split out of ORPHEUS). Measured un-silenced:
five ignore entries already clean (cp/geometry/homogeneous/moc/__init__), 23-error tail, derivations 232.

**Phase 0 — the #236 branch design audit** (user: evaluate before touching the gated 45):
verdict SOUND (main's own docstrings confess the τ-ownership mislabel; live c-from-τ twin;
lying Protocol contracts = the C5 cluster; post-fork trait-vocabulary convergence). ~90%
auto-merges (merge-tree probe). **USER: re-stage on main NOW** → executed on
`feature/sn-spatial-angular-restage` as gated commits (task #148): `82cd210` (Phase-1
traits+pairing+honest curvilinear selection), `cdb3cd1` (closure-owned τ/c + CellVisit stamps,
**sn 99→91**, baseline→137), `9b93db7` (geometry-τ retirement; StreamingTerms loses the M-M
trio; branch's 2 dead type-ignores dropped), `4f9e9b3` (separability gate, 6 legs), `81422f4`
(#248 Protocol retirement + ABC strategy-trio hoist + #249 L2 single-source, **sn 91→86**,
baseline→**132** = 21/86/25; closes #248+#249). Landing catches: the orphaned
`PoleAngularClosure` import in streaming.py (#238 leftover), the post-fork `from_material_mesh`
annotation, the post-fork LD-slope `_l2_2d` call. Docs merge landed as `607b548` (archivist;
+1078 theory lines, 12 dead Protocol xrefs repointed — retirement-audit leg 2; Sphinx `-E -W` 0);
elegance-enforcer PASS-WITH-NITS (folded). **MERGED TO MAIN @ `607b548`** (ff; both #236
branches deleted; #236 CLOSED; #248/#249 close at push; docs follow-up → #286; full
affected-tree 3345/0). `refactor/pyright-burndown` recreated off the new main — **resume at
U1 from ratchet 132** (the brief = task #139). Remaining in the two pole files = the
lazy-Optional cluster (__init__ 39 + pole_angular_closure 26) → C4/C5 proper.

**U1 DONE (2026-07-02, post-compaction): the Level-1 minimal-architecture scope is complete
— PUSHED.** `a6cdf85` un-silenced geometry+homogeneous (measured 0 on lift); `ffab593`
un-silenced data and fixed its 6 at the root: `compute_macro_xs`'s three accumulators
(SigP/SigS/Sig2) got their typed empty-sum identities (`np.zeros(NG)` /
`csr_matrix((NG,NG))` as `sum(..., start=)` — retires the `if fissile_indices else` branch;
FP-exact), and matpro's `k_He/k_Kr/k_Xe` got the honest scalar domain `float | np.ndarray`
(bodies already self-coerce). tests/data 209/0. Ratchet stays **132** (21/86/25) with
data/geometry/homogeneous guarded at 0 via the key-union. ff-merged to main @ `ffab593` and
**pushed** (`1729647..ffab593` — the 6 held #236 commits published; #248/#249 closed at
push; per-cluster push cadence RESUMED by the user). Next: C1 (numerics 16, task #142).

**C2 DONE (2026-07-02): ratchet 115 → 110 — MERGED @ `46d8f79`, PUSHED.** numerics → ZERO.
`FunctionSpace` became `Generic[Carrier]` (PEP-696 `default=Any` — bare operator slots hold
either realization; `FunctionSpace[Domain]` slot-specialization = declared follow-on); the
metric quartet split SURFACE (Carrier-typed, so FullFieldSpace's overrides are exactly the
base at Carrier=CompositeField — LSP-clean by construction) from REALIZATION (`_diagonal_*`
private bodies); `norm` genuinely carrier-generic. The duck-typed composite contract made
static as `CompositeField(Protocol)` (+`__dataclass_fields__` leaf so `replace()` checks).
Enforcer PASS-WITH-NITS; `_recombine → recombine` rename filed as **#287**. (Lesson: pyright
runs override checks against the base signature regardless of self-specializations — the
generic surface is the honest home, not `self: FunctionSpace[NDArray]` scoping.)

**C3 DONE (2026-07-02): ratchet 110 → 83 — MERGED @ `03a33a1`, PUSHED.** transport 24→1,
sn 86→82 free. Ledger −5 casts −2 ignores. **The F2 carve** (user-flagged smell): the SNMesh
narrowing lives on the FIELD declarations; erased at ONE seam (`FullField.bulk` cross-method
slot) → new `FullField.mesh` property reads the one mesh off the BOUNDARY leaf (the
`__post_init__` identity invariant made readable) → 4/5 casts retired; the 5th
(`as_per_ordinate`) = irreducible true boundary interlocked with #288. fission's composite
arm PARSES its real contract (AngularFlux-only). harmonic_frame: closed-grid typed
@overload + isinstance dispatch (open singledispatch registry retired); `basis` covariantly
narrowed to `SphericalHarmonicBasis`; `from_galerkin` = parse boundary. Root fixes:
`Frame.conjugate → OperatorProduct` (annotation lagged its own docstring) →
`full_scatter_kernel` typed, `apply_transpose` static; type-matches-input @overloads on
`LegendreMomentScattering.apply` + `add_iso/add_n2n`; `build_aniso_source` param narrowed
(D-I.2 runtime retirement, annotation lied since). TraceSpace `layout`+`omega_dot_n`
REQUIRED kw_only (bare-ctor guard retired). material_xs_field cell-index cache re-declared
VARIADIC `tuple[np.ndarray, ...]` (np.where arity = mesh ndim; a first-attempt (ix,iy)
unpack broke three 1-D tests — the runtime gate caught it). **THE #288 FINDING**: the #207
cross-class dunder (`Scalar+Angular→Angular`) is statically UNSPELLABLE under LSP (return
must be ≤ Self for the same-class arm) — dunders stay deliberately untyped, docstrings point
at #288 (3 resolutions, USER decision); the 1 residual transport error = the accepted
visible cost. Gates: 1133/0 + sn slice 276/0; enforcer PASS-WITH-NITS (conditions folded).
**NEXT: C4 (task #145, sn ~22 + F1 composition-leg casts) → C5 (task #146, pole-file
lazy-Optional ~60) from ratchet 83 = sn 82 + transport 1.**

**C4 DONE (2026-07-03): ratchet 83 → 61.** sn 82→60 (the residual 60 = EXACTLY the two C5
pole files: loss_representation/__init__ 34 + pole_angular_closure 26 — every other sn file
ZERO), transport 1 (the #288 cost, unchanged). Ledger −6 casts, 0 casts/ignores added.
**THE F1 CARVE (the keystone)**: OperatorSum/OperatorProduct/ScaledOperator generic over
their LEG types (`Generic[Domain, Codomain, SummandA, SummandB]` / FactorA·FactorB /
ScaledOperand) with PEP-696 defaults (`LinearOperator[Domain, Codomain]`) — legs COVARIANT
(a pinned composition must upcast to the defaulted `__add__` contract) ⟹ read-only
properties over `Final` storage (probe-validated on CLI before production: covariant+default
TypeVars, Final storage, the upcast, @overload overrides returning the pinned subclass —
all 0-error; subclass attribute redeclaration REJECTED, confirming the property design).
The sn compositions pin their legs (InvertibleOperator = OperatorSum[FF, FF, Streaming,
Multiplication]; Scheduled = [..., Invertible, ScaledOperator[FF, FF, SNMaskedBoundary]];
WindowedSweep = OperatorProduct[FF, TFF, BulkAnalysis, Sweep]) → ALL 6 F1 casts retired
incl. the :137 nested double-cast (now `self.b.op`). The fusion rules got typed @overloads
(`L+C→InvertibleOperator`, `(L+C)−B_lower→ScheduledInvertibleOperator`, `P@A⁻¹→WindowedSweep`).
**THE 17**: solver — dead `_zero_within_group_fission` RETIRED (zero callers; doc row
updated); dead `mesh/quadrature/materials` params retired from both fixed-source paths (the
#18-remnant — the bodies only `del`'d them); `_apply_default_bcs`'s per-class replace-kwargs
arm → new `Axis1D.with_uniform_bc(bc)` Protocol verb (each axis owns its endpoint structure);
the builder family (`_within_group_triple/si/krylov`, `_select_si_resolvent`, `_maybe_window`)
typed with honest operator types. sweep_graph — the slice/ndarray-union slab addressing
rebuilt per-branch (`shifted` lives with its representation); **`_CellSolve` split into
`_CellSolveAngular`/`_CellSolveMoment`** (kw_only frozen, REQUIRED buffers — the exactly-one
__post_init__ guard became unrepresentable; `_SweepEmit` gained transitional
`angular_buffers()/moment_buffers()` narrowing accessors, the C5 forward-split noted).
augmented_mesh — pole-closure family got an ABSTRACT `__init__(sn_mesh)` construction
contract; `streaming(axis)` guard rekeyed to the direct `_streaming_axes is None` fact.
**THE ROLE-ERASURE FAMILY → #289 FILED**: `FullField.bulk/.boundary` erase the leaf ROLE
(F2's sibling) — C4 resolved all demand sites with loud isinstance PARSES (evaluate_residual,
boundary._apply_faces, Solution.boundary_flux, the 4 driver exits reifying the
solve-echoes-template contract), inventoried in #289 (FullField leaf-generics carve, own
gated slot; notes the driver V-conflation: operator carrier FullField vs iterate carrier
TimedFullField). **FullField.__add__/__sub__ params honest-widened `T`→`FullField`** (the
runtime `_check_partner` accepts ANY flavor — load-bearing for the timed−timeless stencil;
FullField is its own root, no LSP constraint). Gates: full numerics+transport+sn 3070/0
(9:21) + solve/harness slice 68/0; Sphinx -E -W 0 (doc row touched); enforcer
PASS-WITH-NITS (5 nits folded: #289 cited, Krylov bulk local bound, `_within_group_si`
return typed, _SweepEmit forward-note, tombstone trimmed). Test migration: 4 angular-mode
`_CellSolve(` test ctors → `_CellSolveAngular(`. Local task #18 closes as folded (the
mesh-union remnant retired with the dead params). **NEXT: C5 (task #146) from ratchet 61 =
sn 60 (34 + 26, the two pole files ONLY) + transport 1.**

**C5 DONE (2026-07-03): ratchet 61 → 1 — LEVEL 1 COMPLETE.** sn 60→**ZERO** (both pole
files); the surviving 1 = the accepted #288 transport residual. **F-A — the M-M unbound
legacy mode RETIRED** (option (a), the audit was decisive: zero production callers of any
unbound surface; SNMesh always binds; the base's own guards branded it "legacy"; the
production `is None` fallback at the 1-D scan walk was statically DEAD and would have raised
2 lines later): `MorelMontryAngularSweep(sn_mesh)` REQUIRED, the `if sn_mesh is None` state
block + `_require_mesh_bound` + 4 base RuntimeError guards + the 7 base `| None` state
declarations all retired; dead `_sn_mesh`/`_V`/`_N` storage dropped; Identity's apologetic
`| None` annotations + 2 dead vestigial assigns dropped. The pure-algebra kernel moved BACK
to module level (`compute_psi_half_per_level` + `_psi_half_grid_single_level` — it was born
a free function; Phase 2.2 class-hosted it only to serve the unbound mode). **The SNMesh
closure override became a CLASS param** (`type[PoleAngularClosureBase]` — an instance could
only ever be unbound/foreign-bound since the mesh doesn't exist yet; Pattern 4); the
verbatim-vs-default branch folded to one `closure_cls(self)`. `PsiHalfAngleSeed.is_linear` →
`ClassVar` (the F-D Protocol fix). **F-B — `_SweepEmit` split into the closed type family**
(`_SweepEmitAngular`/`_SweepEmitMoment`, REQUIRED buffers, polymorphic `pure_z`; the
exactly-one guard + C4's transitional accessors retired); consumers: windowed interior
isinstance→`_CellSolve` subclass, oracle interior loud angular-only parse, scan-march binds
per-mode `emit_row`/`finish_octant` closures ONCE (mode-blind row loop, staging choreography
preserved); the sweep chain's return became the honest mode-keyed
`tuple[ndarray, ndarray | None]` across 8 surfaces (CumprodScan + _OneDimScanWalk keep the
narrow always-angular pair — union chosen over 16-signature overload ceremony; both
production callers discard slot 2). **The 1-D scan walk's two-arm preamble moved into its
body arms** (slab vs curvilinear bind their own face views; no cross-arm Optionals; the
curvilinear arm gained a loud M-M isinstance parse — the scan IS the M-M thread). Singles:
`moment_scan_closure` hoisted to the scheme ABC (raising-default capability idiom);
`A_downstream` honest-widened `ndarray | float`; 2 `list[slice | int]` edge indices;
`fi_bar` unconditional (the existing `outflow_inner_bar` discard idiom). **FOUR stale
`type: ignore` retired** ([attr-defined] on the base-property `level_indices`, [index]×2 on
`inflow_full`, [union-attr] on `reduced.coord` — the last three proven dead by the enforcer
with `reportUnnecessaryTypeIgnoreComment`). Retirement audit: graph+grep+ctor legs run;
derivations probes' dead `or MorelMontryAngularSweep()` fallbacks fixed; docs re-pointed
(`:meth:`→`:func:` ×5 incl. the source file); a C5 Development-history entry added.
**Test migration** (coverage-preserving, enforcer-verified): test_compute_psi_half_per_level
rewired to the module function (+ the kw-only pin now covers `psi_half_seed`);
`redistribution_via_live_path` dropped its closure param; bound instances via the new shared
`make_tiny_spherical_sn_mesh`; registry `create(..., sn_mesh=)`. **The phase_c flat-flux
xfail leg was INERT** (the factory built an unbound instance that CRASHED — swallowed by
non-strict xfail; the "cylinder telescopes cleanly" claim was never a live observation):
C5 made it live and the bound run shows BOTH geometries fail → `strict=True` with the live
observation recorded. Gates: full tests/sn 1937/0 (11:09) + touched-file re-run 81/0;
Sphinx `-E -W` 0; enforcer PASS-WITH-NITS (3 folded: the 2+1 dead ignores, the source-file
`:meth:` refs, the xfail tightening). Ledger this cluster: −4 ignores, 0 casts added.

**CAMPAIGN CLOSED (2026-07-03, same session): the user STOPPED the standalone Level-2
grind.** The pivot: **the #290 diffusion-integration campaign** (local task #149) — diffusion
is the genuine consumer of scalar flux + scalar composites that exercises the architecture on
both angular AND scalar carriers, opens DSA for SN (#2), and carries the pyright burn as a
BY-PRODUCT of correct integration (the diffusion ignore lifts with the work; baseline gains
diffusion:0). Level-2 tail parked in task #147 (cp/moc/__init__ zero-lifts, mc/fuel/plotting,
derivations→#19's fork; TH+kinetics NEVER). Milestone comment on #226; the truthful-model
commit-trailer ruling recorded ([[feedback-truthful-model-attribution]]). This log is the
campaign's terminal record.

**C1 DONE (2026-07-02): ratchet 132 → 115 — MERGED @ `95a2f8b`, PUSHED.** numerics 21→5
(residual 5 = the C2 carrier carve), transport 25→24 (side effect). All root-cause designs,
ledger −2 casts / −4 ignores / 0 added: gmres dead-arm retirement (+ the over-broad except);
`_CarrierMatvecOperator(spla.LinearOperator)` subclass at the sole scipy boundary;
`_DisplacementLeaf` Protocol (the L1↛L2 duck-face, consumer-minimal); KEigenvalue
boundary guard; **the `_check_partner` parse-don't-validate carve** (`_is_same_class →
TypeIs[Self]` + return-the-proven-partner — pyright rejects Self-vs-Self param overrides
on contravariance, the covariant RETURN is the honest home; 3 transport overrides chain
`partner = super()._check_partner(other)`); `dual() → FunctionSpace` (reflexive double-dual)
+ `DualSpace.primal` required kw_only; `_as_unit` pint-stub boundary (**closes #258**);
`__pow__` self-annotation `LinearOperator[Domain, Domain]` (endomorphic precondition in the
signature; 2 casts retired). F3 triage: :806/:1509/iteration-:611 = documented true
boundaries, KEPT. Gates: numerics+transport 1133/0, sn slice 270/0, enforcer
PASS-WITH-NITS (doc nit folded; `T`-vs-`Self` collapse deferred to next dunder edit;
`_flux_role.py` mixin ignores = different root cause, C3/C4 scope). Next: C2 (task #143).

### [historical] Resume point of the pre-#226 session — orpheus/ at **420**
The clean single-root-cause autonomous wins are now **spent**. geometry.py's 4 residuals are distinct one-offs (`Mesh2D.areas` L369; `pole_angular_closure` assign/positional L465 = #236; one `None` subscript L768). Remaining `orpheus/sn` (excl. #236-deferred `pole_angular_closure` 34 + `loss_representation` 43) is **heterogeneous** and falls into harder buckets that need user-steering or deep flow work:
- **B4 — `LinearOperator[V]` generic → an ARCHITECTURAL CARVE (DETOUR), evolved into the TYPED-CARRIER GRID carve.** Not a typing patch. Progress (2026-06-22):
  - **Foundation LANDED (`272caa9`):** RC1 `block_role` ClassVar→instance + RC3 `_as_boundary`/`realize` retyped to `LinearOperatorMixin` (+ `α·B` dunder) — **−12 reds (420→408), zero runtime**. (Phase 0 net = `ab75d02`.)
  - **ACTIVE plan: `.claude/plans/typed_carrier_grid_carve.md`** (approved, cold-pickup §0). Operators speak typed `Field` carriers by completing the (angular⊗moment)×(flux⊗source) grid; the domain/codomain two-param split lands here (Phase A) with a real consumer (non-endomorphic `HarmonicMomentProjection`/`Reconstruction`). SUPERSEDES Phase 1 of `operator_inverse_algebra_carve.md`.
  - **Follow-on: `operator_inverse_algebra_carve.md` Phases 2–5** (inverse-as-operator: `solve`→`inverse().apply`, retire `CAP_SOLVE`) builds on Phase A.
  - Net: `issue_226_b4_operator_generics_verification.md`; map: `agent-memory/explorer/issue_226_operator_generics_map.md`. SURGICAL / main-agent-direct + user-steered — do NOT dispatch method-implementer. When these merge, B4's 18 reds clear by architecture; resume at B5 below.
- **B5 — union dispatch**: scattering.py (13: `ndarray | ScalarSourceSink | None`, `.mesh`/`.values`/`.spatial_moments_per_axis` on raw ndarray), solver.py (`BulkField` vs `AngularField`). Needs source-build flow understanding.
- **units.py (5)** — pint-stub-fighting (`PlainUnit | Unknown` vs declared `Unit`); third-party-stub category like the SymPy backlog. Low value; defer with B6.
- Then **B6**+**Workstream C** (derivations ~228), **Workstream A** (user-gated hook).

Operating rules in force: `npx pyright` CLI is the oracle (the streamed `<new-diagnostics>` lags edits — never trust it over a CLI run); NO `# type: ignore`; re-tighten `pyright_baseline.json` after each cluster (`python -m tests._harness.pyright_ratchet --update`); per-cluster ff-merge + push (user authorized this cadence). **Corrected paths**: ratchet test = `tests/test_pyright_ratchet.py`; `test_reduced_operator.py` is under `tests/geometry/`.

---

## Current state (triage, 2026-06-21, `npx pyright` 1.1.410)

Full repo = **2353 errors / 19 warnings**, 774 files. By directory:

| dir | errors | nature |
|---|---|---|
| `tests/` | 1403 | test-stub / `**kwargs`-splat / loose-fixture idioms |
| `orpheus/` | **502** | production — the real target (260 non-derivation + 242 SymPy) |
| `derivations/` (repo-root) | 267 | SymPy-fighting |
| `scratch/` | 129 | throwaway |
| `student_resources/` | 37 | teaching |
| `examples/` | 12 | demos |
| `tools/` | 3 | |

**No regression** — the #257 branch *reduced* `orpheus/` 538→502 (−36, all `orpheus/sn`). The 2353 is full-repo scope; historical "538" was always `orpheus/`-only (`pyrightconfig.json` has no `include`/`exclude` → bare `npx pyright` walks the whole root).

Production (`orpheus/`) rule histogram: reportArgumentType 224 · reportAttributeAccessIssue 98 · reportOptionalSubscript 38 · reportReturnType 25 · reportOperatorIssue 19 · reportAssignmentType 15 · reportCallIssue 14 · reportOptionalMemberAccess 14 · reportUndefinedVariable 11 · …

---

## Workstream A — NOISE elimination (do FIRST; it's the "not noise" half)

The streamed `<new-diagnostics>` reportMissingImports avalanche ("Import 'orpheus.transport.full_field' could not be resolved", ".spatial.scheme", etc.) is **0-real in the CLI for `orpheus/`** — it's the langserver dual-identity artifact (the worktree/exec-env rooting + microsoft/pyright#10498). The fix lives in `.claude/hooks/regen-pyrightconfig.sh` (the per-worktree `executionEnvironment` design) — the two `.claude/worktrees/*` checkouts already contribute 0 errors via that mechanism.

**A1. Diagnose the residual LSP noise.** Why does the LSP still stream phantom import errors when the CLI is clean? Confirm whether the in-editor/agent langserver is using the regenerated `pyrightconfig.json`, whether stale worktree exec-envs or a missing root mapping cause the dual-identity, and whether microsoft/pyright#10498 has a config-level mitigation. Dispatch explorer + read the hook + reference_pyright_lsp_rooting memory.
**A2. Make the LSP agree with the CLI.** Extend/repair `regen-pyrightconfig.sh` (or the langserver invocation) so the streamed diagnostics match `npx pyright orpheus/` — phantom imports gone. This is the single biggest agent-UX win (every turn is currently polluted). **This file is in `.claude/hooks/` (commit-protected) — propose the change to the user, do not self-commit it.**
**A3. Pin the scope.** Decide + document the canonical analysis scope: should bare `npx pyright` analyze the whole repo, or should `pyrightconfig.json`/`pyproject.toml [tool.pyright]` set `include = ["orpheus"]` (production) with tests/derivations/scratch under separate, relaxed exec-envs? Pinning scope makes "the count" stable and meaningful (no more 2353-vs-538 confusion).

---

## Workstream B — real-error reduction (the "only signal" half)

Attack order = highest signal-cleared-per-unit-effort, re-measuring between big clusters. Each cluster: explorer/enumerate → fix → tests green + count drops + 0 net-new ignores → commit.

**B1. Cluster D + undefined-vars (~17, trivial, zero-risk) — warm-up.** `numpy.bool_`/`csr_array` vs annotated return → `bool(...)`/correct annotation (~6). The 11 `reportUndefinedVariable` = missing `import sympy as sp` (7×), `Quadrature1D`/`TimedFullField` TYPE_CHECKING forward-refs (4×). Mechanical, runtime-inert.

**B2. Cluster C — Protocol/ABC missing concrete attrs (~12, low-risk).** `PoleAngularClosure` missing `level_indices`/`precompute_psi_state`/`cell_contribution` (8×); `AngularMeasure`, `MOCMesh` (rest). Declare the members on the ABC/Protocol. ⚠ VERIFY against the #248 `PoleAngularClosure`-Protocol retirement first (some accesses may be on a stale name).

**B3. Cluster A — Optional/None not narrowed (~104, the biggest single bucket).** reportOptionalSubscript 38 + OptionalMemberAccess 14 + OptionalOperand 5 + the `... | None` argtype/operator errors. **Two files hold ~77:** `orpheus/sn/loss_representation.py` (43) + `orpheus/sn/spatial/pole_angular_closure.py` (34). Root: lazily-populated `Optional`/`= None` fields (`_alpha_per_level`, `_tau_per_level`, `_dAw_per_level`, `psi_half_seed`) subscripted after an invariant guarantees they're set; + `float | None` params used in arithmetic. Fix per-site (judgement): declare the real non-Optional type + init in `__post_init__`/builder, OR `assert x is not None` at the invariant-established point, OR make the param non-Optional. ⚠ Respect the lazy-init contract — don't force eager allocation where the None sentinel is load-bearing. **Re-measure after** — narrowing one lazy field clears cascading downstream argtype/operator errors.

**B4. Cluster B — `LinearOperator[V]` under-parametrised generic (~24+, architecturally coupled).** `LinearOperator[Unknown]` loses `.solve`/`.apply_transpose`; `block_role` "cannot assign"; `TensorProductOperator`/`IncomingSourceOperator` not assignable (7×). Concentrated in `boundary_realizer.py` (9), `operator.py` (3), `iteration.py` (3), `solver.py` (3). Fix = ONE coherent hierarchy edit per `.claude/agent-memory/explorer/issue_226_operator_generics_map.md` (unbounded `V=TypeVar("V")`): parametrise the `LinearOperator`/`Mixin`/`OperatorSum` family `Generic[V]`; declare `block_role` as a settable field on the composers; surface `.solve` via the resolvent type / a Protocol. Type-level, runtime-inert; re-measure (likely clears uncounted downstream argtype collapse).

**B5. Cluster E — stringly/union dispatch (~10, ergonomics).** `str|SubgroupOfO3|int|dict` arg in `numerics/quadrature/registry.py` (8×); `Mesh1D|Mesh2D|tuple[Axis1D,...]` mesh arg (2×); `Quadrature→AngularQuadrature` (3×). Tighten signatures / overload / narrow with `isinstance`.

**B6. Cluster F — `orpheus/derivations` SymPy backlog (242, ISOLATE, last).** Unchanged from main; SymPy `Expr`/`Float`/`Rational` fight pyright. Fix the genuinely-real ones, then a **SCOPED per-directory `executionEnvironment`** relaxation for `orpheus/derivations/**` (NOT global) — the triage doc's recommended handling. Lowest production value.

---

## Workstream C — the non-production trees (1403 tests + 267 root-derivations + 129 scratch + …)

Lower value, larger count. Decide policy rather than grind every one:
- **`tests/`**: many are test-stub idioms (`_StubQuad` not a `Quadrature`, `**kwargs`-splat into typed factories, `BC.vacuum` enum access, `float`→`int` arg looseness). Triage: fix the cheap structural ones (the `solve_*(geometry, **dict)` splat collapses pyright cleanly into a typed helper — see the S10a peierls-test refactor precedent); for the rest, consider a relaxed test-tree exec-env (tests are not shipped product). DO NOT let test-stub noise gate production signal.
- **`derivations/` (repo-root), `scratch/`, `student_resources/`, `examples/`**: scratch/teaching — strong candidates for scoped exec-env relaxation or `exclude`, NOT line-by-line fixes.

The principle: production (`orpheus/`, ex-derivations) → **zero real errors**; the SymPy/test/scratch trees → **scoped relaxation** so they don't drown the signal, with the genuinely-real ones fixed.

---

## Sequencing
1. **A1–A3 noise-rooting FIRST** (kill the streamed avalanche — makes all subsequent work legible; the hook change goes to the user since `.claude/hooks/` is protected).
2. **B1 → B2** (zero/low-risk warm-up, ~29 cleared).
3. **B3** (the 104 Optional/None bucket; re-measure).
4. **B4** (the LinearOperator-Generic coherent edit; re-measure).
5. **B5**, then **B6** (derivations scoped relaxation).
6. **Workstream C** policy (scope/relax the non-production trees).

## Verification (every batch)
- `npx pyright orpheus/` count strictly drops; report before/after.
- 0 net-new `# type: ignore` (grep the diff).
- The relevant test suites stay green under `-O` (route around the documented baseline reds #250/#232/#212).
- After A2: the streamed `<new-diagnostics>` for an `orpheus/` edit matches the CLI (no phantom imports).

## Open investigations (dispatch at execution)
- explorer: the LSP dual-identity root cause + whether `regen-pyrightconfig.sh` fully resolves it (A1).
- explorer: confirm the #248 `PoleAngularClosure` retirement status before B2.
- the `issue_226_operator_generics_map.md` memo is the B4 design input — re-read it.

## Relations
#226 (the pyright-signal effort), the prior `.claude/plans/pyright_cluster_triage.md` (superseded counts: 691→502 orpheus, undefined-var 67→11, missing-imports 0), `reference_pyright_lsp_rooting.md`, `issue_226_operator_generics_map.md`. #258 (units.py pint-stub debt) + #254/#255 feed here.
