---
name: nd-generic-carve-tuple-axis
description: The d-generic SN sweep carve (C3) replaces named per-axis pairs with positionally-ordered tuples; the axis↔index convention + the WavefrontFlux typing destination, and where each belongs (kernel vs walk).
metadata:
  type: project
---

The dimension-agnostic wavefront carve (C3, branch `worktree-sn-nd-layout`,
`SweepDependencyGraph.from_cartesian(shape)` + d-generic `DiamondDifference`
cell kernel) produces a new recurring elegance shape. Reviewed C3.0/C3.1/C3.2a
(`5b941eb`/`9b75374`/`9063487`).

**Why:** the d-generic API replaces the named `psi_in_x`/`psi_in_y` pair with
`psi_in: tuple[ndarray, ...]` (element `a` = face normal to spatial axis `a`).
This is correct for d-genericity but creates two standing review axes.

**How to apply (the durable rulings):**

1. **Two-family SSOT is NOT a twin — accept the split.** The Cartesian batch
   kernel (`diamond.py cell_kernel_batch`: `denom = Σ_t + Σ_a s_a`, multi-axis,
   no curvature/no angular-redistribution) and the curvilinear `cell_balance_for_streaming`
   (`denom = 2|μ|·A_down + (ΔA/w)·c_out + Σ_t·V`, single-axis + Morel–Montry)
   share only the SHAPE `denom·ψ̄ = q + numer`, not the algebra. Forcing them
   through one helper = Pattern-6 abstract-over-the-difference trap (Cartesian
   would carry phantom neutral curvature into the hot path). The "single source
   of the DD cell math" docstring claim is scoped to the CARTESIAN family; the
   curvilinear family lives in `cell_balance.py`. Verdict ACCEPT.

2. **The positional axis convention (psi_out[0]=x, [1]=y) belongs typed at the
   WALK, raw at the KERNEL.** `WavefrontFlux` (`transport/fields/wavefront_flux.py`)
   is the already-built typed interior 1-cochain with `face(axis)` views + an
   `axes` property — the destination for the axis↔index map. But do NOT type the
   KERNEL arg as `WavefrontFlux`: the kernel touches a transient `(N_oct,ng,n_diag)`
   anti-hyperplane SLICE, not the whole cochain (category error), and `WavefrontFlux`
   is explicitly storage-granularity-indifferent (cross-domain-attacker rejected
   per-face objects at this granularity). The win is at C3.2b's d-generic
   `SweepCellSlice`: route the per-axis gather/scatter through `WavefrontFlux.face(a)`
   so the positional tuple order becomes the typed `axes` property. Until then,
   pin the convention with a one-line docstring contract at the kernel definition
   site (Pattern 7).

3. **`| None` optional-optimization fields are the right model for ONE
   accelerator.** `window_slots`/`window_edges: tuple|None` (None for d≠2) with
   fail-loud guards in `apply_windowed`/`residual_windowed` is correct: the
   windowed walk is a pure optimization pinned bit-identically against the
   all-d full-field `apply`/`residual` oracle; a d≠2 graph is COMPLETE, just
   without that accelerator (method-unavailable, not object-malformed). When the
   d=1 cumprod accelerator lands (C3.4) there will be TWO `|None` optimization
   fields → that is the unify-after-two trigger to consider a
   `SweepOptimization|None` sum type. Do NOT pre-build (one instance now).

4. **Transitional partial-function shims with a PLAN-tracked retirement are
   ACCEPT.** `OctantLabel.sign_x`/`sign_y` = `signs[0]`/`signs[1]` (IndexError on
   d=1 `.sign_y`) + `streams_in_2d` alias are explicitly slated to retire in
   C3.5 (per-axis FaceLabel + twin retirement), tracked in
   `.claude/plans/nd_layout_foundation.md` C3 sub-staging — the skill's exact
   "temporary needs a removal trigger" requirement. Audit at C3.4/C3.5 pickup
   that no d=1 path reaches `.sign_y` before the shims retire.

**Recurring FIX-NOW tells in d-generic carves:**
- A kwarg shadowing a same-scope local of a DIFFERENT type: `s = slice_args`
  (the SweepCellSlice) adjacent to `self.cell_kernel_batch(..., s=(sx,sy))`
  (the streaming-coeff tuple). Rename the PARAM (`s`→`streaming`), not the
  established `s = slice_args` accessor idiom.
- Positionally-ordered tuple crossing a boundary with axis↔index identity
  unstated. The closure `out[a]=2·avg−in[a]` is axis-DIAGONAL → a transposed
  tuple passes any symmetric (Δx=Δy) test and only fails on asymmetric grids.
  Demand the positional contract at the definition site.

**C3.2b ruling — WHERE the typed axis↔buffer map lands (the ruling-#2 follow-through).**
Reviewed C3.2b (uncommitted working tree on `worktree-sn-nd-layout`). Ruling #2 above said
"route the per-axis gather/scatter through `WavefrontFlux.face(a)` so the positional tuple
order becomes the typed `axes` property." The implementer did NOT store a `WavefrontFlux`
inside `SweepCellSlice`; instead the typed map lands at the ORCHESTRATOR boundary
(`tuple(wavefront.face(a) for a in wavefront.axes)` in `sweep.py:_sweep_2d_full_field`,
`operator.py:_apply_2d..full_field`) and the slice/kernel carry the raw octant-projected
tuple in that axis order. **This is the CORRECT call, not a shortcut** — and ruling #2's
"at C3.2b's d-generic SweepCellSlice" phrasing was too prescriptive about the storage site.
The reason is load-bearing and verified: the buffers handed to the walk are OCTANT-RESTRICTED
copies `psi_faces[a][oct_idx].copy()` → shape `(N_oct, …)`, NOT the full mesh `(N, …)`.
`WavefrontFlux._check_partner` (`wavefront_flux.py:182`) is mesh-bound and the type's
`__post_init__` requires `values.shape == (space.layout.total_size,)` keyed on the FULL mesh
— so wrapping an `N_oct`-buffer would be a type LIE (illegal-state in the other direction,
Pattern 4 inverted). The axis↔buffer binding a type would enforce is established ONCE at the
producer (Pattern 7 normalise-at-definition: the hardcoded `psi_x=face(0);psi_y=face(1)` pair
is retired in favour of the cochain's own `.axes` map), and the octant projection + slice
preserve that order positionally. The `SweepCellSlice` "Why a raw per-axis tuple, not a
WavefrontFlux" docstring section states this correctly. **The right level for the type IS the
whole-mesh cochain; the slice is a transient sub-block where the type cannot validly exist —
do NOT demand a fictitious octant-submesh `WavefrontFlux`.** (One residual nit: `str_axes` is
re-spelled positionally `(str_x, str_y)` / `(streaming_x, streaming_y)` adjacent to the
`.axes`-derived `psi_faces` — a SECOND, untyped axis-order spelling that the `.axes` map does
not gate. Harmless at d=2, a latent transpose habitat when a 3rd accelerator/axis lands;
flagged CONCERN not VIOLATION since both are length-2 in the only live path.)

**C3.2b ruling — the gather/scatter factoring is a GENUINE Pattern-2 win.**
`_gather_cell_inputs`/`_scatter_outgoing_faces`/`_cell_face_selector` (`diamond.py`) collapse
a real twin: `update_batch` (solve) and `residual_batch` (apply) had byte-identical front-half
(incoming faces + streaming + Σt + Q gather) and back-half (outgoing scatter), differing ONLY
in the cell algebra (`cell_kernel_batch` vs `residual_kernel_batch`). Two concrete consumers
existed BEFORE the extraction → not premature (Pattern 6 satisfied). The `2 + axis` lattice
offset in `_cell_face_selector` is the SAME edge-locating arithmetic as
`WavefrontFlux._edge_slot` (`wavefront_flux.py:228`, "skips the (N,ng) prefix") — consistent,
not a new third spelling. Verdict ACCEPT/PASS.

**C3.3 ruling (`bf6b55c`) — `is_1d` becomes the ONE partition predicate; `not is_1d` is honest TODAY.**
C3.3 made the geometry octant build d-generic (`itertools.product((-1,+1), repeat=ndim)` +
`from_cartesian(spatial_shape, label=OctantLabel(signs))`) and flipped `is_1d` from `ny==1` to
`ndim==1`, then routed all 5 streaming-operator dispatch gates (`_compute_LpC`,
`_compute_decomposition`, `_compute_LpC_transpose`, `apply`, `apply_transpose`) from `ny>1` to
`not sn_mesh.is_1d`. **This RESOLVES the "two spellings of one partition" tell** (AGENT.md §7):
`ny>1` (operator) and `ny==1` (is_1d) were two spellings; now there is ONE predicate
(`is_1d ≡ ndim==1`) with `is_1d` getting its FIRST production consumers. STRENGTH, not a smell.
- **The `not is_1d` vs explicit `ndim==2` question (review axis #1): NOT a defect today.**
  `SNMesh` wraps ONLY `Mesh1D`/`Mesh2D` (no `Mesh3D` dataclass exists — geometry.py:745
  "no Mesh3D dataclass exists today"), so for any CONSTRUCTIBLE `SNMesh`, `not is_1d ≡ ndim==2`
  exactly. The d=3 trap (where `_apply_2d_cartesian` is 2-D-specific but `not is_1d` would also
  match d=3) is UNREACHABLE in production — a 3-D mesh cannot be built. Latent-but-gated:
  when `Mesh3D` lands, `apply`'s `not is_1d` gate would misroute d=3 into the 2-D FD kernel.
  CONCERN-not-VIOLATION because (a) unconstructible today, (b) `_apply_2d_cartesian` would
  shape-mismatch loudly on a (N,ng,nx,ny,nz) field, (c) the C3 campaign's own d=3 admission
  pins (`test_sweep_graph_nd_admission`) are the tracked trigger. The honest end-state is a
  3-way dispatch (chain-scan / window-walk / full-field-walk) keyed on `ndim`, not a binary
  `is_1d`; that lands with the d≥3 walk, not here. Bug-habitat: a future `Mesh3D` author who
  trusts `not is_1d` to mean "2-D". Mitigation already present: the RAISE sites correctly read
  "multi-D Cartesian … not yet wired" (not "2-D"), so only the live `apply`/`_apply_2d` DISPATCH
  carries the latent assumption, and it is documented as the 2-D path.
- **Eager d=1 graph build (review axis #2) — JUSTIFIED, not speculative.** `_setup_cartesian`
  now builds 2 genuine d=1 chain graphs for every 1-D slab even though production 1-D still uses
  the cumprod scan (consumer `_wavefront_1d_sweep` lands C3.4). This is NOT Pattern-6 speculation:
  (a) the graphs are a PURE FUNCTION of mesh topology, built once at construction, cost ~5 lines
  of arange, zero per-sweep cost; (b) the consumer is KNOWN and IMMINENT (next sub-step, plan-tracked);
  (c) the alternative — gating the build on dimensionality — would REINTRODUCE the per-d special-case
  the carve exists to remove. Building all `2^d` octants uniformly IS the d-genericity. The pre-C3.3
  state already eagerly built 4 (phantom) graphs for 1-D; C3.3 makes them 2 GENUINE ones. Net:
  less dead weight, not more. PASS.
- **`from_cartesian_2d` retirement — COMPLETE + CLEAN.** Zero live call sites remain
  (`grep from_cartesian_2d(` → empty); the method body is deleted; the 3 surviving textual
  references are docstrings/comments that correctly DESCRIBE it as retired (not stale calls);
  13 test sites migrated to `from_cartesian((nx,ny), label=)`; the `_needs_spine`
  `xfail(strict=False)` scaffold + "not yet landed" framing retired now that the builder is
  permanent. Test migration done (not delete-only) per `feedback_retirement_means_test_migration`.
- **`_legacy_2d_anti_diagonal` golden — GENUINELY STRUCTURALLY-INDEPENDENT.** Verified the two
  constructions differ: production `from_cartesian` does `np.indices(shape).reshape(d,-1)` →
  `local.sum(axis=0)` → boolean-mask `level_of == k` (full-lattice partition by index-sum); the
  golden does the explicit `local_i ∈ [max(0,k-ny+1), min(nx-1,k)]`, `local_j = k - local_i`
  per-level anti-diagonal recurrence. Different algorithms (set-partition-by-sum vs per-level
  index-range recurrence) that MUST agree — a real oracle, not a tautological re-call of the
  builder. The non-square `(5,7)`/`(4,3)` params catch an x↔y transpose the square grid hides.
  This is the textbook "re-derive the golden, don't call the retired method" pattern. PASS.
- **Reads-like-the-domain (review axis #4) — PASS.** `for signs in itertools.product((-1,+1),
  repeat=self.ndim)` reads aloud as exactly "the 2^d streaming sign-signatures over d axes". No
  stringly-typed dispatch, no boolean-flag parameter, no procedural transcription introduced.
  The 1-D test's `np.testing.assert_array_equal`/`pytest.fail` (not bare `assert`) is correct
  `-O`-safe hygiene (Mode-8), and the test docstring honestly scopes the surrounding bare-assert
  gap as pre-existing.

**Verdict: PASS.** One tracked CONCERN (the `not is_1d` d=3 latent gate, unconstructible today,
to revisit when `Mesh3D` + the d≥3 walk land — it is the SAME 3-way-dispatch destination the
campaign already targets, so no new issue needed; the plan's d≥3 admission pins are the trigger).

**Done well this carve (reinforce):** the explicit left-fold
(`for s_a,in_a in zip(s,psi_in): denom=denom+s_a; numer=numer+s_a*in_a`)
chosen over `sum()` to preserve `((sigt+s_0)+s_1)` IEEE-754 association for d=2
bit-identity — thoroughly documented at the point of the what (both kernel
docstrings + commit + 6920-ULP-unchanged gate). This is how a bit-identity
micro-decision should read. The C3.0→C3.2a scaffold-to-real transition retired
the `if dd is None: pytest.fail` import-keepalive hack and wired B5(d=3)/B6c(d=1)
to the real kernel against a structurally-independent hand oracle.
