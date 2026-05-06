# SN Reshape — Operator Algebra Architectural Campaign

**Date**: 2026-05-05 (updated 2026-05-06)
**Branch**: `refactor/sn-operator-algebra`
**Scope**: SN module + numerics primitives + geometry additions
**Pilot for**: cross-solver migration sequence SN → MoC → CP → MC
**Status**: **Phase 0 + Phase 1 + Phase 2 (DD-only) complete (Issues 1-7 + 8 + 9-DD); Wave C-extension (LD/EC/Step) and Wave D (Phase 3) pending**

---

## Progress (as of 2026-05-06)

Campaign branch `refactor/sn-operator-algebra` ahead of `main` by 13 commits.
Bit-identical regression contract held throughout: 11/11 frozen snapshots at
`tests/sn/regression/snapshots/` remain `np.array_equal`-bit-identical to
the pre-reshape baseline.

### Pre-reshape (already on `main` before campaign branch)

- Issue 16 (regression suite, GH [#165](https://github.com/deOliveira-R/ORPHEUS/issues/165)) — 11 frozen snapshots installed at `tests/sn/regression/snapshots/`
- Anisotropic curvilinear MMS infrastructure with 2 `xfail(strict=True)` ERR-026 tripwires at [tests/sn/l1_analytical/test_mms_curvilinear_aniso_dd_convergence.py](../../tests/sn/l1_analytical/test_mms_curvilinear_aniso_dd_convergence.py) — designed to flip to xpass when Issues 11/12 close ERR-026, forcing marker removal
- vv-principles SKILL extended with failure mode #7 (MMS simplification bias) at [.claude/skills/vv-principles/SKILL.md](../skills/vv-principles/SKILL.md)
- 15 L1 homogeneous + 12 anisotropic MMS foundation tests (~6.2s combined) per-commit gate
- GH [#149](https://github.com/deOliveira-R/ORPHEUS/issues/149) — divide-warning at `solver.py:192` filed for Issue 13/15 owner during Wave H
- 18 reshape issues filed at GH #150–167 on milestone "SN Reshape (Wave 1)"

### Wave A — Phase 0 numerics primitives (DONE)

| Plan # | GH # | Status | Commit | Summary |
|---|---|---|---|---|
| 1 | [#150](https://github.com/deOliveira-R/ORPHEUS/issues/150) | ✓ DONE | `60d3932` | `numerics/operator.py` — LinearOperator Protocol with capability tags; `LinearOperatorMixin` for dunder algebra; OperatorSum/Product/Scaled/Identity/Zero composers; `as_scipy_linop` adapter; **38 foundation tests** |
| 2 | [#151](https://github.com/deOliveira-R/ORPHEUS/issues/151) | ✓ DONE | `2f67853` | `numerics/measure.py` — DiscreteMeasure with tensor product / pushforward / restrict / direct sum; BundleMeasure for MoC (Wave 2); `gauss_legendre`/`gauss_chebyshev`/`equispaced` 1D rules; **22 L1 + 17 foundation tests** |

### Wave B — Phase 0 (Issues 3-5) + Phase 1 (Issues 6-7) (DONE)

| Plan # | GH # | Status | Commit | Summary |
|---|---|---|---|---|
| 3 | [#152](https://github.com/deOliveira-R/ORPHEUS/issues/152) | ✓ DONE | `724547d` | `numerics/symmetry.py` — `SubgroupOfO3` enum + static containment lattice + `is_invariant`; orbit-closure check with O_h (48 elements) and I_h (120 elements) generator sets via Rodrigues + BFS closure; **71 foundation tests** |
| 4 | [#153](https://github.com/deOliveira-R/ORPHEUS/issues/153) | ✓ DONE | `0a9aa94` | `numerics/quadrature/{rules_1d,rules_sphere,rules_product}.py` + thin-adapter refactor of `orpheus/sn/quadrature.py` (502 → 464 LOC). Hot-path attributes (`mu_x`, `mu_y`, `weights`, `level_indices`, `_ref_*`) cached as numpy views on construction; `reflection_index(axis)` reimplemented as cached pushforward. **First SN-touching change**, regression contract held: 11/11 bit-identical |
| 5 | [#154](https://github.com/deOliveira-R/ORPHEUS/issues/154) | ✓ DONE | `60f9fb2` | `numerics/quadrature/registry.py` — `QuadratureSpec` dataclass with structural flags + populated registry + `select_quadrature(geometry, target_degree, **structural)` precedence-chain selector with `SelectionLog` explainability; **37 foundation tests** |
| 6 | [#155](https://github.com/deOliveira-R/ORPHEUS/issues/155) | ✓ DONE | `e4276e1` | `geometry/reduced_operator.py` — Bailey 2009 connection-coefficient lift (α dome recursion, ΔA/w redistribution, M-M τ_mm clamp) into `ReducedStreamingOperator` + `StreamingTerms` + factories (slab/cylindrical/spherical). **Hash equality with `SNMesh._setup_*` outputs** verified via `np.array_equal`; 29 parametrized hash-equality tests across N ∈ {4, 8, 16, 32} (sphere) and (n_mu, n_phi) ∈ {(2,4), (4,4), (4,8)} (cylinder) |
| 7 | [#156](https://github.com/deOliveira-R/ORPHEUS/issues/156) | ✓ DONE | `e93fe47` | `geometry/boundary.py` — `ResolvedBC` Protocol with tensor-decomposition framing `R = Σ_α G_α ⊗ A_α`; concrete `VacuumBC`, `SpecularBC`, `WhiteBC`, `PeriodicBC`, `AlbedoBC`, `MixedBC`. `SNMesh.BC_REGISTRY` factories return `ResolvedBC`; sweep at `orpheus/sn/sweep.py` calls `apply_to_incoming(...)`. WhiteBC/PeriodicBC ship as primitives only (not yet wired into `solve_sn`). Backward-compat `__eq__` shim on VacuumBC/SpecularBC accepts string comparisons (transitional). **19 primitive tests, regression contract held: 11/11 bit-identical** |

Wave A + Wave B verification gates (end-of-session):

- 398 numerics + geometry tests on integrated campaign HEAD
- 11/11 regression snapshots bit-identical (629s)
- 27 + 2 xfail safety net (L1 homogeneous + anisotropic MMS foundation + ERR-026 tripwires)
- Sphinx clean (-W); audit clean (22 orphans, 36/38 ERR coverage — both unchanged from baseline)
- Full L1 sweep exit 0

### Wave C — Phase 2 cell-update strategy (Issue 8 + DD-portion of Issue 9) (DONE)

| Plan # | GH # | Status | Commit | Summary |
|---|---|---|---|---|
| 8 | [#157](https://github.com/deOliveira-R/ORPHEUS/issues/157) | ✓ DONE | `e999c15` | `orpheus/sn/spatial/cell_update.py` — `CellUpdate` `@runtime_checkable Protocol` (slot-style traits `is_linear`/`is_positivity_preserving` + `update(streaming_terms, total_xs, source, upstream_state) → CellResult`). `UpstreamState` and `CellResult` `@dataclass(frozen=True, slots=True)` carry per-cell `(ng,)`-shaped state. `StreamingTerms` extended additively with `volume` and `abs_mu` fields populated by all 3 factories (slab/sphere/cylinder). Slab vs curvilinear discriminated via `streaming_terms.alpha_in is None`; cylindrical pure-azimuthal `|η|<1e-15` degenerate case via `outgoing_spatial_flux=None`. **13 protocol tests + 16 vol/abs_mu tests = 29 new foundation tests; existing 29 hash-equality tests preserved (45 total geometry tests on `test_reduced_operator.py`)** |
| 9 (DD) | [#158](https://github.com/deOliveira-R/ORPHEUS/issues/158) | ✓ DONE (DD portion) | `3b1fc75` | `orpheus/sn/spatial/diamond.py` — `DiamondDifference` strategy as **bit-identical extraction** of existing inlined sweep math. Single geometry-polymorphic class with 3 branches: slab (`alpha_in is None`, `s = 2·source/denom` matching the cumprod path's `bQ = 2·source_coeff·Q_1d` semantics), curvilinear (`alpha_in is not None, abs_mu ≥ 1e-15`, mirrors `sweep.py:350-361` operation order verbatim), cyl-degenerate (`alpha_in is not None, abs_mu < 1e-15`, no spatial DD closure). `is_linear=True, is_positivity_preserving=False` (Lewis & Miller §5.3). **11 hand-calc bit-identical tests via `np.array_equal` against scalar formulas typed verbatim from sweep.py:117-123/350-361/533-546**. LD/EC/Step deferred to Wave C-extension |

Wave C verification gates (end-of-session):

- 69 spatial + geometry tests (13 protocol + 11 diamond + 45 geometry incl. 16 new vol/abs_mu)
- 11/11 regression snapshots bit-identical — sweep is provably untouched in Wave C, so the contract holds by construction; verified twice (post-R1 worktree run 617s + post-R2 agent run 599s)
- 27 + 2 xfail safety net intact (ERR-026 tripwires gated for Wave D, did NOT flip)
- Sphinx clean (-W exit 0); audit unchanged from Wave-B baseline (23 orphans of 234 testable labels — denominator grew from 231 to 234 because Wave C added 3 new equation labels [`dd-slab-scalar`, `dd-mm-closure-constants`, `dd-curvilinear-scalar`]; all 3 covered by `@pytest.mark.verifies` in test_diamond.py so orphan count unchanged)
- 36/38 ERR coverage preserved

### Operational notes (lessons from Wave A + Wave B + Wave C)

1. **Worktree-base bug** (discovered Wave B Round 1, mitigated for the rest): background-dispatched worktrees may come up at `main`'s HEAD (06b46f2) instead of the orchestrator's current branch HEAD. Mitigation: brief every method-implementer / general-purpose dispatch with explicit detection (`git status && git log --oneline -3 && ls orpheus/...`) and recovery (`git rebase refactor/sn-operator-algebra`) instructions. Issues 3, 4, 5, 7 all hit this and recovered cleanly via rebase.
2. **Cherry-pick conflict pattern** (consistent across all 3 rounds): conflicts on `docs/verification/matrix.rst` (Sphinx auto-regenerates with different test-count totals per branch). Resolution: `git checkout --ours` on the conflicting file; complete the cherry-pick; re-run `sphinx-build`; commit the regenerated matrix as a `chore(matrix)` commit.
3. **Bit-identical contract is the campaign's success criterion** — and held: every SN-touching commit (Issues 4 and 7) was gated by `pytest -m regression -q` (629–640s, 11/11 bit-identical) before merge. The agent verifications + post-merge orchestrator verifications agreed every time.
4. **Sub-agent assignment heuristic**: `general-purpose` for software-only primitives without published-math content (Issues 1, 3, 5); `method-implementer` for issues with bit-identical contracts or published-math grounding (Issues 2, 4, 6, 7). Dispatching parallel pairs (one general-purpose + one method-implementer) per round avoided sequential bottlenecks.
5. **`type:refactor` label does not exist** in the GitHub repo. Reshape issues marked `type:refactor` in the plan (#153, #156, #157, #159, #161, #162, #164) were filed with `type:improvement` instead. Acceptable substitution.
6. **Wave C duplicate-edit side effect** (Round 1 only): the `general-purpose` agent's worktree-isolated dispatch landed file edits in BOTH the agent's worktree AND the main repo (probably a tooling oversight when the agent's CWD-tracker was reset between Bash calls). Mitigation: orchestrator stashed the duplicated main-repo edits before cherry-picking the agent's worktree commit; the cherry-pick was clean since the stash removed the working-tree conflict surface. Round 2 (`method-implementer`) did NOT exhibit this. Worth surveilling on future `general-purpose` worktree dispatches.
7. **Wave C source-semantics confirmation** (Round 2 lesson): the `source` parameter passed to `CellUpdate.update(...)` follows the sweep's call-site convention — `source = Q · V · weight_norm` (already weight-normalized AND already volume-multiplied). The slab cumprod path bakes `2 · weight_norm · dx` into `source_coeff` so the per-cell scalar form becomes `s = 2 · source / denom`; the curvilinear path inserts `source` directly into the numerator as `QV[i]`. Both branches verified bit-identical via `np.array_equal` hand-calc tests typed verbatim from `sweep.py:117-123/350-361`.

### Remaining waves (next sessions)

| Wave | Phase | Issues | Description |
|---|---|---|---|
| C-ext | 2 | 9 (LD, EC, Step) | Concrete cell updates beyond DD: `LinearDiscontinuous`, `ExponentialCharacteristic`, `Step` strategies. Each with its own MMS spatial-convergence test. Most natural to ship sequentially (one round per strategy) AFTER Wave D's unified sweep dispatches via `cell_update`, so the strategies can be verified end-to-end as production cell updates rather than shipped in isolation |
| D | 3 | 10, 11, 12, 13 | SN core reshape: SNMesh refactor, SNStreamingOperator, unified sweep (consumes `DiamondDifference`), ScatteringOperator + FissionOperator. **Closes ERR-026** via Issues 11+12 (the 2 xfail-strict tripwires at `tests/sn/l1_analytical/test_mms_curvilinear_aniso_dd_convergence.py` flip xpass when curvilinear sweep produces O(h²) on anisotropic ansatz) |
| E | 4 | 14, 15 | Iteration as operator algebra (SourceIteration, KEigenvalue) + SNSolver migration. **Closes #96, #97 via Issue 15** |
| F | 5 | 17, 18 | Symmetry-preservation + reciprocity invariant tests; Sphinx documentation campaign |

ERR-026 closure is the implicit success criterion for Issues 11 + 12: the
2 anisotropic MMS xfail-strict tripwires at
`tests/sn/l1_analytical/test_mms_curvilinear_aniso_dd_convergence.py`
will flip from xfail to xpass when the curvilinear sweep produces
O(h²) convergence on the anisotropic ansatz, forcing marker removal.

---

## Campaign context

The codebase has a working multigroup, anisotropic-scattering, steady-state SN
solver with three coordinate systems (Cartesian slab, cylindrical, spherical).
The math is correct — Bailey et al. (2009) angular redistribution, Morel-Montry
weighted diamond closure for spherical, per-level azimuthal redistribution for
cylindrical, all referenced in code.

What's missing is **architectural unification across the deterministic solver
family**. The pieces are consistent inside SN but the abstractions don't span
SN/MoC/CP/diffusion. This campaign installs the cross-cutting primitives —
`LinearOperator`, `DiscreteMeasure`, `ReducedStreamingOperator`,
`CellUpdate` — using SN as the pilot. MoC, CP, and MC migrate sequentially
in their own campaigns once SN is proven.

The architectural narrative is: the transport equation has differential,
integral, and variational forms. SN/MoC discretize the differential form;
CP discretizes the integral form; PN/diffusion discretize spectral truncations.
At the operator-algebra level they all expose `(L, S, F)` triples consumed
uniformly by source iteration, Krylov, and eigenvalue solvers. That's the
end state.

## Confirmed decisions

1. **1D cumprod fast path lives inside the unified sweep algorithm**, not as
   a separate dispatched function. It's a DD-specific algebraic optimization;
   keep it but express it as such.
2. **Issues #96 / #97 close as a side effect of the reshape**: the BiCGSTAB
   path stops building a separate FD operator and instead does matrix-free
   Krylov on `SNStreamingOperator.apply`. Behavioral change in the
   BiCGSTAB convergence is acceptable.
3. **`numerics/eigenvalue.py:power_iteration` lingers** until CP, diffusion,
   and other consumers migrate to operator-algebra iteration. No hard
   removal date. **Note**: full migration of all solvers (SN → MoC → CP → MC)
   IS planned, sequentially. The new infrastructure must be designed
   forward-compatibly with all four.

## Cross-solver migration sequence

This SN reshape is **Wave 1**. Subsequent waves consume the same primitives:

| Wave | Solver | Notes |
|------|--------|-------|
| 1 | SN | this campaign — installs primitives |
| 2 | MoC | uses `BundleMeasure`, ray-quadrature DiscreteMeasure |
| 3 | CP | dense-matrix `LinearOperator`, scalar-flux operator |
| 4 | MC | uniform interface for verification, shares `Geometry` |

Primitives must be designed in this campaign without MoC/CP/MC blockers.
`BundleMeasure` ships in Phase 0 even though SN doesn't use it.

## Architectural commitments (invariants)

These are non-negotiable across all issues. Any deviation is a session
failure per Cardinal Rule 2 (architecture is critical).

1. **`LinearOperator` is the cross-cutting abstraction** for L (streaming +
   collision), S (scattering), F (fission). Capability-tagged
   (`apply` / `solve` / `apply_transpose`), composable, capability-checked
   at construction. Operators that lack a capability raise on attempted
   call, not on import.

2. **`DiscreteMeasure` is the quadrature primitive.** Composable via tensor
   product (`__mul__`), pushforward, restrict, direct sum (`__add__`).
   `AngularQuadrature` becomes a thin SN-domain adapter over it.

3. **Connection coefficients live on geometry**, exposed through
   `ReducedStreamingOperator`. Cell updates consume `StreamingTerms`,
   never raw α / ΔA / τ arrays.

4. **Cell updates are polymorphic** with a curvilinear-aware signature:
   `update(streaming_terms, total_xs, source, upstream_state) → CellResult`.
   `upstream_state.angular_upstream` is `None` for Cartesian, a per-cell
   per-group array for curvilinear.

5. **Production schemes are positivity-preserving by construction**: LD,
   exponential characteristic, step. DD kept only as comparison artifact,
   marked `is_positivity_preserving = False`.

6. **Existing patterns are extended, not replaced.** `BC_REGISTRY`,
   `SNMesh`, `Mesh1D.from_geometry`, the geometry-layer two-tier split
   (reference vs production solvers), all stay. The reshape adds layers
   above and refactors internals; it does not replace the public shape.

7. **Behavioral regression is the gating contract.** Existing SN tests
   produce bit-identical results when `cell_update == DiamondDifference`.
   Snapshot-frozen scalar-flux + k_eff outputs are captured BEFORE the
   reshape begins (Issue 16, executed first).

## Out of scope

- MoC, CP, MC reshapes (separate campaigns)
- Time-dependent / kinetics (steady-state stays; transient is a downstream
  campaign building on this primitives layer)
- New angular quadrature families (LS_N already implemented; Lebedev variants
  beyond what scipy provides are not expanded here)
- New geometry kinds (HSPH, ANN, 2D curvilinear) — these slot in cleanly
  once `ReducedStreamingOperator` exists
- Anisotropic-scattering Pℓ ordering changes (current implementation in
  `SNSolver._build_aniso_scattering` is preserved as-is, just lifted to
  `ScatteringOperator`)

---

## Issue plan

Each subsection is a self-contained spec for a GitHub issue. Format:

- Module label, type, phase, dependencies, complexity (S / M / L)
- Context (the WHY — Cardinal Rule 3)
- Acceptance criteria (checkbox list — gating for PR merge)
- Files affected
- Design notes (decisions, gotchas, references)

A sub-agent should translate each subsection 1:1 into a GitHub issue with
the appropriate `module:` and `type:` labels.

### Phase 0 — Numerics primitives

#### Issue 1: `numerics/operator.py` — `LinearOperator` protocol

- **Module**: `module:numerics`
- **Type**: `type:feature`
- **Phase**: 0
- **Depends on**: none
- **Complexity**: M

**Context**: Today `scipy.sparse.linalg.LinearOperator` is used only inside
the BiCGSTAB Krylov path (`orpheus/sn/operator.py`). The sweep, scattering
source, and fission source live as bare functions / methods. Lifting them
to a uniform `LinearOperator`-shaped interface lets eigenvalue and Krylov
code consume any solver method (SN/MoC/CP) without knowing which transport
discretization is below. This is the foundation for the cross-solver
migration sequence.

**Acceptance criteria**:

- [ ] `LinearOperator` Protocol in `orpheus/numerics/operator.py` with
      `apply(x)`, optional `solve(b)`, optional `apply_transpose(x)`,
      and `capabilities: frozenset[str]` property
- [ ] Composition primitives: `OperatorSum`, `OperatorProduct`,
      `ScaledOperator`, `IdentityOperator`, `ZeroOperator` — each computes
      its own capability set from constituents (sum requires `apply` from
      both; product requires composition of capabilities)
- [ ] Operator algebra: `__add__`, `__sub__`, `__mul__` (scalar),
      `__matmul__` (operator product) on `LinearOperator`
- [ ] Capability mismatches raise `MissingCapability` at composition time,
      not at call time
- [ ] Adapter `as_scipy_linop(op)` for Krylov consumption via scipy
- [ ] Unit tests on synthetic operators (numpy matrices wrapped) covering:
      composition, capability propagation, scipy interop

**Files**:

- new: `orpheus/numerics/operator.py`
- new: `tests/numerics/test_operator.py`

**Design notes**: Follow the runtime-checkable Protocol pattern already
established for `AngularQuadrature` in `orpheus/sn/quadrature.py`. The
capability set is the key idea — many operators have no efficient `solve`
(S, F), and forcing them to provide stubs is harmful. Don't abstract
shape/dtype — leave that to numpy duck-typing for now.

---

#### Issue 2: `numerics/measure.py` — `DiscreteMeasure` primitive

- **Module**: `module:numerics`
- **Type**: `type:feature`
- **Phase**: 0
- **Depends on**: none
- **Complexity**: M

**Context**: Quadratures are mathematically discrete measures. The natural
operations (tensor product → product quadrature, pushforward → coordinate
change, restriction → half-range) all have direct mathematical content.
Today these are implicit: `ProductQuadrature.create(n_mu, n_phi)` constructs
the product internally without exposing the tensor-product structure.
Promoting `DiscreteMeasure` to a primitive lets consumers compose 1D rules
into 2D rules into S² rules cleanly, and is **required for MoC's bundle
measures** (Wave 2).

**Acceptance criteria**:

- [ ] `DiscreteMeasure` in `orpheus/numerics/measure.py`: `nodes`,
      `weights`, `space` (str/enum tag), `integrate(f)`, `__mul__` (tensor
      product), `__add__` (direct sum), `pushforward(phi)`,
      `restrict(predicate)`, `n_points`
- [ ] `BundleMeasure` for disintegrated measures — base measure plus
      per-base-point fiber measure factory. Required for MoC; not used
      by SN in this campaign but lives here so the abstraction is correct
      from day one
- [ ] 1D primitives: `gauss_legendre(n)`, `gauss_chebyshev(n)`,
      `equispaced(a, b, n)` — each returns a `DiscreteMeasure`
- [ ] Optional metadata fields populated lazily: `invariance_group`
      (filled by Issue 3), `degree_of_exactness`
- [ ] Unit tests: integration of polynomials of known degree, tensor
      product equivalence to manual construction, pushforward correctness
      under invertible maps

**Files**:

- new: `orpheus/numerics/measure.py`
- new: `tests/numerics/test_measure.py`

**Design notes**: Don't try to enforce `Space` types via Python generics —
not expressive enough without runtime overhead. Use `space: str | enum`
as a runtime tag. Composition operations return new `DiscreteMeasure`s
with metadata combined sensibly. `BundleMeasure` is critical to ship now
so MoC migration doesn't have to revisit Phase 0.

---

#### Issue 3: `numerics/symmetry.py` — Subgroups of O(3) and invariance

- **Module**: `module:numerics`
- **Type**: `type:feature`
- **Phase**: 0
- **Depends on**: Issue 2
- **Complexity**: S

**Context**: Quadrature selection (Issue 5) needs subgroup containment
lattice logic: geometries declare their symmetry group, quadratures
declare theirs, selection requires `G_geom.is_subgroup_of(G_quad)`.
The lattice is finite for ORPHEUS-relevant cases — a static lookup
suffices; no need for generator-based machinery.

**Acceptance criteria**:

- [ ] `SubgroupOfO3` enum-backed in `orpheus/numerics/symmetry.py` with
      named entries: `Trivial`, `Z2`, `SO2`, `O2`, `OctahedralOh`,
      `IcosahedralIh`, `SO3`, `O3`, plus parameterized `Cn(n)`,
      `Dnh(n)` for hex/lattice future use
- [ ] `contains(self, other) -> bool` implementing subgroup containment
      via static lookup table
- [ ] `is_invariant(self, measure: DiscreteMeasure) -> bool` checks that
      measure nodes form a union of orbits with consistent weights
- [ ] Unit tests: containment lattice (Z2 ⊂ O_h ⊂ O3, etc.); invariance
      checks on existing quadratures (Lebedev → O_h, LS_N → O_h,
      Gauss-Legendre on μ → SO2)

**Files**:

- new: `orpheus/numerics/symmetry.py`
- new: `tests/numerics/test_symmetry.py`

**Design notes**: Don't overengineer. Static dict mapping
`(G_a, G_b) → bool` for containment is fine. Generator-based machinery
deferred until ORPHEUS needs novel discrete groups (hex / triangular
lattices: C_6v, D_6h — add when consumed).

---

#### Issue 4: Refactor `AngularQuadrature` to use `DiscreteMeasure`

- **Module**: `module:numerics`, `module:sn`
- **Type**: `type:refactor`
- **Phase**: 0
- **Depends on**: Issues 2, 3
- **Complexity**: M

**Context**: `orpheus/sn/quadrature.py` has 4 working quadrature classes
(GL1D, Lebedev, LS_N, Product) pre-dating the `DiscreteMeasure`
abstraction. Each exposes a similar but ad-hoc interface (`mu_x`, `mu_y`,
`weights`, `N`, `reflection_index`, `spherical_harmonics`). Refactoring
them to be `DiscreteMeasure`-backed — with the existing
`AngularQuadrature` Protocol kept as a domain-specific adapter — gives
composability without breaking SN's heavy attribute-access consumption.

**Acceptance criteria**:

- [ ] `orpheus/numerics/quadrature/` package: `rules_1d.py`,
      `rules_sphere.py`, `rules_product.py`
- [ ] Each rule function returns a `DiscreteMeasure` with
      `invariance_group` and `degree_of_exactness` populated
- [ ] `orpheus/sn/quadrature.py` refactored: existing classes (GL1D,
      Lebedev, LS_N, Product) become thin adapters wrapping
      `DiscreteMeasure` and caching SN-specific fields (`mu_x`, `mu_y`,
      `mu_z`, `level_indices`, `_ref_x` etc.) on construction
- [ ] `reflection_index(axis)` re-implemented as `pushforward` with
      coordinate negation; result cached
- [ ] `spherical_harmonics(L)` stays where it is (SN-specific)
- [ ] **All existing SN tests pass without modification.** This is the
      gating constraint — no behavior change

**Files**:

- new: `orpheus/numerics/quadrature/{__init__,rules_1d,rules_sphere,rules_product}.py`
- refactor: `orpheus/sn/quadrature.py`
- existing tests: unchanged

**Design notes**: The bridge issue between the new primitive and existing
SN consumption. Trick: SN consumers index into `quad.mu_x`, `quad.weights`
etc. heavily. Don't break that. The adapter caches array views from the
backing `DiscreteMeasure` on construction.

---

#### Issue 5: Quadrature registry with G + V + structural tags + selector

- **Module**: `module:numerics`
- **Type**: `type:feature`
- **Phase**: 0
- **Depends on**: Issues 3, 4
- **Complexity**: S

**Context**: With G-tagged quadratures and a subgroup lattice, automated
selection becomes the precedence chain: G compatibility → V compatibility
→ structural compatibility → minimum points. `solve_sn` keeps explicit
quadrature-passing as the canonical API; the registry adds
`select_quadrature(geometry, ...)` as a convenience and as a documentation
artifact. The structural tags are themselves teaching content.

**Acceptance criteria**:

- [ ] `QuadratureSpec` dataclass: `name`, `factory`, `invariance_group`,
      `degree_of_exactness`, structural flags (`positive_weights`,
      `axis_aligned`, `level_structured`, `half_range_clean`)
- [ ] `quadrature_registry` populated with all current rules
- [ ] `select_quadrature(geometry, target_degree, **structural)
      → DiscreteMeasure` implementing the precedence chain
- [ ] Selection log returned alongside the chosen quadrature so the choice
      is explainable / loggable
- [ ] Unit tests: slab → GL1D (SO(2)-invariant after isotropy reduction);
      sphere → GL1D on μ_r; cylinder → ProductQuadrature; Cartesian-like
      2D → LS_N or Lebedev with explainable preference

**Files**:

- new: `orpheus/numerics/quadrature/registry.py`
- new: `tests/numerics/test_registry.py`

**Design notes**: Don't make selection mandatory. Explicit
quadrature-passing in `solve_sn` stays canonical; the registry is
opt-in convenience. The structural tags double as Sphinx documentation
content.

---

### Phase 1 — Geometry layer additions

#### Issue 6: `ReducedStreamingOperator` per geometry

- **Module**: `module:geometry`
- **Type**: `type:feature`
- **Phase**: 1
- **Depends on**: Issue 1
- **Complexity**: L

**Context**: Connection coefficients (α, ΔA/w, τ_mm) currently live in
`SNMesh._setup_spherical` and `SNMesh._setup_cylindrical`. They are
SN-specific consumers of geometry-side math. MoC and CP will need the
same reduced-operator concept (different consumption pattern, same
underlying math). Lifting it to geometry now removes future duplication.
Per Cardinal Rule 2: shared concepts between solvers means the codebase
needs an architectural overhaul.

**Acceptance criteria**:

- [ ] `ReducedStreamingOperator` class in
      `orpheus/geometry/reduced_operator.py`
- [ ] Properties: `requires_upstream_angular_state: bool`,
      `angular_marching_axis: Literal["mu", None]`
- [ ] Method: `streaming_terms(cell_idx, direction_idx, mu_level_idx=None)
      → StreamingTerms`
- [ ] `StreamingTerms` carries everything the cell update needs in this
      geometry: chord lengths, face areas, ΔA factor, α coefficients,
      τ_mm — geometry-dependent shape
- [ ] Concrete factory functions (or subclasses): `slab_streaming(mesh)`,
      `cylindrical_streaming(mesh, angular_measure)`,
      `spherical_streaming(mesh, angular_measure)`
- [ ] Output is bit-identical to current `SNMesh._setup_*` (verified by
      hash-equality of arrays in tests)
- [ ] Unit tests against precomputed reference values matching current
      `SNMesh.alpha_half`, `SNMesh.redist_dAw`, `SNMesh.tau_mm`

**Files**:

- new: `orpheus/geometry/reduced_operator.py`
- new: `tests/geometry/test_reduced_operator.py`

**Design notes**: The current `SNMesh` couples geometry math (connection
coefficients) with SN math (quadrature). Split: connection-coefficient
*values* depend on quadrature points (they're integrals of α from the
angular discretization), so the factory takes both mesh and angular
measure. Output is coordinate-system-aware. Per-level structure for
cylindrical (consuming `angular_measure.level_indices` via the adapter
from Issue 4) is preserved.

---

#### Issue 7: Boundary conditions as tensor decompositions

- **Module**: `module:geometry`, `module:sn`
- **Type**: `type:refactor`
- **Phase**: 1
- **Depends on**: Issue 4
- **Complexity**: M

**Context**: BCs are currently `BC(kind, params)` dataclasses with
SN-specific resolution in `SNMesh.BC_REGISTRY` returning a string tag.
Math says BCs are tensor decompositions: R = Σ_α G_α ⊗ A_α where G_α is
geometric (permutation/index map) and A_α is response (albedo, scalar
amplitude). Lifting the resolved BC to a tensor-decomposed object
makes specular / white / albedo / periodic / mixed all uniform, and
lets multi-region interfaces reuse the same primitives.

**Acceptance criteria**:

- [ ] `ResolvedBC` Protocol/ABC in `orpheus/geometry/boundary.py`:
      `apply_to_incoming(angular_flux_outgoing, quadrature)
      → angular_flux_incoming`
- [ ] Concrete: `VacuumBC` (zero), `SpecularBC` (rank-1: permutation
      `pushforward` from `DiscreteMeasure` + albedo scalar), `WhiteBC`,
      `PeriodicBC`, `AlbedoBC`, `MixedBC` (rank-N sum)
- [ ] `SNMesh.BC_REGISTRY` updated: factories return `ResolvedBC`
      instances, not string tags. The factory pattern stays
- [ ] Sweep code in `orpheus/sn/sweep.py` updated to call
      `resolved_bc.apply_to_incoming(...)` instead of branching on
      string kind
- [ ] All existing SN tests pass with bit-identical outputs for `vacuum`
      and `reflective`
- [ ] Unit tests for `WhiteBC` and `PeriodicBC` (currently unsupported
      by SN — adding support is a downstream win that this issue enables)

**Files**:

- new: `orpheus/geometry/boundary.py`
- modify: `orpheus/sn/geometry.py`, `orpheus/sn/sweep.py`

**Design notes**: This is where BC tensor-decomposition framing pays off
architecturally. Use the existing `BC_REGISTRY` pattern — just change
what the factories return. White and periodic come essentially free
once specular is correct.

---

### Phase 2 — SN cell update strategies

#### Issue 8: `CellUpdate` ABC with curvilinear-aware signature

- **Module**: `module:sn`
- **Type**: `type:refactor`
- **Phase**: 2
- **Depends on**: Issue 6
- **Complexity**: M

**Context**: Current sweep dispatches by `sn_mesh.curvature` ("spherical"
/ "cylindrical" / None) and inlines the cell update equation per
geometry. To swap spatial schemes (LD, exponential, characteristic)
systematically, the per-cell update needs to be a strategy. The
signature must accommodate curvilinear coupling — the previous
μ-direction's flux feeds into the connection term.

**Acceptance criteria**:

- [ ] `CellUpdate` ABC in `orpheus/sn/spatial/cell_update.py`:
      `update(streaming_terms, total_xs, source, upstream_state)
      → CellResult`
- [ ] Properties: `is_linear: bool`, `is_positivity_preserving: bool`
- [ ] `UpstreamState` dataclass: `spatial_upstream` (always present),
      `angular_upstream` (None for Cartesian, per-cell per-group array
      for curvilinear)
- [ ] `CellResult` dataclass: `cell_average_flux`,
      `outgoing_spatial_flux`, `outgoing_angular_state` (None for
      Cartesian)
- [ ] No concrete strategies yet (Issue 9). Just the ABC and
      protocol-level tests confirming any concrete strategy can be
      substituted

**Files**:

- new: `orpheus/sn/spatial/__init__.py`,
       `orpheus/sn/spatial/cell_update.py`
- new: `tests/sn/spatial/test_cell_update_protocol.py`

**Design notes**: Extract the cell-update equation from existing
`_sweep_1d_spherical` (lines ~230–290 in `orpheus/sn/sweep.py`).
The numer/denom assembly per cell is the cell update; the loop
structure around it is the sweep. This issue establishes the contract;
Issue 9 implements concrete strategies; Issue 12 rebuilds the sweep
around the new ABC.

---

#### Issue 9: Concrete cell updates — DD, LD, ExponentialCharacteristic, Step

- **Module**: `module:sn`
- **Type**: `type:feature`
- **Phase**: 2
- **Depends on**: Issue 8
- **Complexity**: L

**Context**: Production schemes need positivity by construction. LD
(linear discontinuous) and exponential characteristic are workhorses;
step is the robust fallback. DD is kept as comparison artifact and
explicitly marked non-positive.

**Acceptance criteria**:

- [ ] `DiamondDifference` in `orpheus/sn/spatial/diamond.py` — implements
      the existing math from `_sweep_1d_spherical` /
      `_sweep_1d_cylindrical` / `_sweep_2d_wavefront`.
      `is_positivity_preserving = False`. Default-deselected for
      production; available for verification
- [ ] `LinearDiscontinuous` in `orpheus/sn/spatial/linear_discontinuous.py`
      — slope-limited LD as recommended default.
      `is_positivity_preserving = True`. Cite Lewis & Miller §5.3 in
      docstring
- [ ] `ExponentialCharacteristic` in `orpheus/sn/spatial/exponential.py`
      — analytical exponential within cells, optimal for thin/thick
      optical regimes. `is_positivity_preserving = True`
- [ ] `Step` in `orpheus/sn/spatial/step.py` — first-order step
      characteristic, robust last resort. `is_positivity_preserving = True`
- [ ] All four pass MMS verification — link to
      `derivations/continuous/mms/sn/`
- [ ] Spatial convergence tests per scheme on a benchmark
      slab/sphere/cylinder problem demonstrating expected order

**Files**:

- new: `orpheus/sn/spatial/{diamond,linear_discontinuous,exponential,step}.py`
- new: `tests/sn/spatial/test_*.py`

**Design notes**: Don't ship all four in one PR. Sequence:
DD first (move existing code, verify bit-identical results), then LD,
then exponential, then step. Each cell update works for slab AND
curvilinear via the connection-term inputs in
`streaming_terms.alpha_in / alpha_out / tau` and
`upstream_state.angular_upstream`. **The sequencing matters**: DD
landing first preserves the regression contract throughout the rest
of the campaign.

---

### Phase 3 — SN core reshape

#### Issue 10: Refactor `SNMesh` to consume `ReducedStreamingOperator`

- **Module**: `module:sn`
- **Type**: `type:refactor`
- **Phase**: 3
- **Depends on**: Issue 6
- **Complexity**: M

**Context**: `SNMesh._setup_spherical` and `_setup_cylindrical` compute
connection coefficients inline. Issue 6 lifted that math to
`geometry/reduced_operator.py`. Now `SNMesh` becomes a thinner
aggregator: hold mesh, hold quadrature (DiscreteMeasure-backed), hold
reduced operator, hold resolved BCs. No streaming math in `SNMesh`
itself.

**Acceptance criteria**:

- [ ] `SNMesh.__init__` calls
      `geometry.reduced_streaming_operator(mesh, quadrature)` instead
      of `_setup_*` methods
- [ ] All `SNMesh` consumers (sweep, BiCGSTAB operator) read connection
      coefficients via `sn_mesh.reduced.streaming_terms(...)` instead of
      `sn_mesh.alpha_half` etc.
- [ ] Backward-compatible attribute access via properties for one
      release: `sn_mesh.alpha_half` returns
      `sn_mesh.reduced.alpha_half` with a `DeprecationWarning`
- [ ] All existing SN tests pass

**Files**:

- refactor: `orpheus/sn/geometry.py`

**Design notes**: After this issue, `SNMesh` knows nothing about α
coefficients or τ weights — it just holds and routes. The
geometry-vs-SN split becomes visible in code.

---

#### Issue 11: `SNStreamingOperator` as `LinearOperator`

- **Module**: `module:sn`
- **Type**: `type:feature`
- **Phase**: 3
- **Depends on**: Issues 1, 8, 9, 10
- **Complexity**: L

**Context**: The sweep is the implementation of L⁻¹. Wrapping
sweep + streaming-residual + adjoint sweep in a `LinearOperator`
interface lets eigenvalue solvers and Krylov code consume L uniformly
across solver methods. **This is also where #96 / #97 close**: `apply`
and `solve` use the same cell update strategy, so the FD/sweep
inconsistency disappears by construction.

**Acceptance criteria**:

- [ ] `SNStreamingOperator` in `orpheus/sn/operator.py` (replacing the
      BiCGSTAB-only `transport_operator_matvec_*` functions) implements
      `LinearOperator`
- [ ] `solve(q)` invokes the unified sweep (Issue 12)
- [ ] `apply(psi)` computes forward streaming-collision (Ω·∇ + Σ_t)ψ
      using the same cell update as the sweep — needed for residuals
      and matrix-free Krylov
- [ ] `apply_transpose(psi)` is the adjoint sweep: reversed Ω,
      transposed cell update
- [ ] `capabilities = frozenset({"apply", "solve", "apply_transpose"})`
- [ ] Reciprocity test: ⟨Lψ, φ*⟩ = ⟨ψ, L*φ*⟩ to round-off
- [ ] Closes #96 and #97 (subject to Issue 15 landing)

**Files**:

- refactor: `orpheus/sn/operator.py`

**Design notes**: The current `operator.py` has a documented
inconsistency between the BiCGSTAB FD operator and the DD sweep
(issues #96 / #97). This issue makes them the same operator: matrix-free
`apply` is the algebraic dual of the sweep using the same cell update.
The behavioral change in BiCGSTAB convergence is acceptable per
campaign decision 2.

---

#### Issue 12: Unified sweep — single algorithm, parameterized

- **Module**: `module:sn`
- **Type**: `type:refactor`
- **Phase**: 3
- **Depends on**: Issues 6, 8, 10
- **Complexity**: L

**Context**: Currently `transport_sweep` dispatches into 4 separate
sweep functions (`_sweep_1d_cumprod`, `_sweep_1d_spherical`,
`_sweep_1d_cylindrical`, `_sweep_2d_wavefront`). Dispatch is by
`sn_mesh.curvature`, conflating "what coordinate system" with "what
discretization scheme." After this refactor: 2 sweep algorithms
(Cartesian flow-ordered + curvilinear μ-marching) parameterized by
cell update; choice between them comes from
`reduced_streaming_operator.requires_upstream_angular_state`.

**Acceptance criteria**:

- [ ] `orpheus/sn/sweep.py` rewritten with `_cartesian_sweep(L, q)` and
      `_curvilinear_sweep(L, q)`. Dispatch on
      `L.geometry.reduced.requires_upstream_angular_state`, not on
      string-comparing `curvature`
- [ ] Per-cell math delegated to `L.cell_update.update(...)` with
      appropriate `streaming_terms` and `upstream_state` assembled
- [ ] **1D cumprod fast path preserved as an optimization within the
      Cartesian sweep** when `cell_update == DiamondDifference` and
      quadrature is GL1D — algebraic identity holds for DD specifically;
      keep for performance, but as a fast path inside the unified
      algorithm, not a separate dispatched function
- [ ] 2D wavefront diagonal scheduling preserved
- [ ] All existing SN tests pass with bit-identical results when
      `cell_update == DiamondDifference`
- [ ] New tests with `cell_update == LinearDiscontinuous`,
      `ExponentialCharacteristic`, `Step` exercise the polymorphism

**Files**:

- rewrite: `orpheus/sn/sweep.py`

**Design notes**: Most behavior-preserving issue in the campaign.
Bit-identical output for DD case is the gating contract. The 1D
cumprod recurrence is an algebraic identity holding for DD only —
keep it as fast path inside the unified Cartesian sweep.

---

#### Issue 13: `ScatteringOperator` and `FissionOperator` as `LinearOperator`s

- **Module**: `module:sn`
- **Type**: `type:refactor`
- **Phase**: 3
- **Depends on**: Issue 1
- **Complexity**: M

**Context**: Scattering and fission sources are currently methods on
`SNSolver` (`_add_scattering_source`, `_build_aniso_scattering`,
`compute_fission_source`). Lifting them to `LinearOperator`s makes the
operator algebra `(L − S − F) ψ = q` explicit and lets eigenvalue
solvers consume the pieces uniformly across solver methods.

**Acceptance criteria**:

- [ ] `ScatteringOperator(LinearOperator)` in `orpheus/sn/scattering.py`:
      holds materials + quadrature; `apply(psi)` returns scattering
      source (P0 + Pℓ); `capabilities = frozenset({"apply"})`
- [ ] `FissionOperator(LinearOperator)` in `orpheus/sn/fission.py`:
      holds materials; `apply(psi)` returns χ ⊗ νΣ_f acting on flux
      moments. The rank-1-in-energy tensor structure χ ⊗ νΣ_f is
      reflected in the implementation
- [ ] (n,2n) folded into `ScatteringOperator` (consistent with existing
      `_add_n2n_source` placement in `SNSolver`)
- [ ] All existing SN tests pass

**Files**:

- new: `orpheus/sn/scattering.py`, `orpheus/sn/fission.py`
- modify: `orpheus/sn/solver.py` (remove source-construction methods)

**Design notes**: Pℓ anisotropic scattering is already implemented in
`SNSolver._build_aniso_scattering`. **Move** the implementation; don't
rewrite the math. The Galerkin projection on Y_ℓm is what's happening —
make it explicit in the docstring.

---

### Phase 4 — Iteration as operator algebra

#### Issue 14: New iteration primitives in `numerics/iteration.py`

- **Module**: `module:numerics`
- **Type**: `type:feature`
- **Phase**: 4
- **Depends on**: Issue 1
- **Complexity**: M

**Context**: Current `EigenvalueSolver` Protocol bundles transport
solving with scattering iteration in `solve_fixed_source`. This
conflates concerns and makes adjoint, Krylov-on-the-outer-iteration,
and contour-integral eigenvalue solvers awkward. New shape: iteration
consumes `LinearOperator` instances directly.

**Acceptance criteria**:

- [ ] `SourceIteration(L, S, F, preconditioner=None)` in
      `orpheus/numerics/iteration.py` — fixed-source solver
      `solve(q_ext) → psi`. `(I − L⁻¹·S)·ψ = L⁻¹·(F·ψ + q)` implemented
      as a fixed-point iteration
- [ ] `KEigenvalue(L, S, F, eigenvalue_method="power")` — eigenvalue
      solver consuming operators directly. Other methods to land
      separately: `"contour_integral"` (the FEAST-style extension you
      had identified)
- [ ] Existing `power_iteration` in `numerics/eigenvalue.py` deprecated
      via `DeprecationWarning` but kept functional for transitional
      CP/diffusion compatibility
- [ ] Both new primitives consume any `LinearOperator` triple — agnostic
      to whether L is SN, MoC, or CP
- [ ] Tests against synthetic operators (numpy matrices wrapped) plus
      tests using SN operators from Issues 11, 13

**Files**:

- new: `orpheus/numerics/iteration.py`
- modify: `orpheus/numerics/eigenvalue.py` (add deprecation, keep
  functional)

**Design notes**: Don't delete the old `EigenvalueSolver` protocol —
CP and diffusion still use it and migrate in their own waves.
`power_iteration` can linger until the cross-solver migration sequence
completes.

---

#### Issue 15: Migrate `SNSolver` to operator-algebra iteration

- **Module**: `module:sn`
- **Type**: `type:refactor`
- **Phase**: 4
- **Depends on**: Issues 11, 13, 14
- **Complexity**: L

**Context**: With L, S, F as operators and `KEigenvalue` consuming them,
`SNSolver` shrinks dramatically. Most of its current ~600 lines becomes
adapter code; physics moves to operators.

**Acceptance criteria**:

- [ ] `SNSolver` rewritten in `orpheus/sn/solver.py`: constructs
      `SNStreamingOperator` (L), `ScatteringOperator` (S),
      `FissionOperator` (F); delegates to `KEigenvalue(L, S, F)` for
      eigenvalue path; delegates to `SourceIteration(L, S, F)` for
      fixed-source path
- [ ] `solve_sn` and `solve_sn_fixed_source` public APIs preserved
      verbatim
- [ ] BiCGSTAB inner solver path becomes: GMRES around L using
      `L.apply` matrix-free, instead of building a separate FD operator.
      **Closes #96 and #97**
- [ ] All existing SN tests pass with bit-identical results in the
      `inner_solver="source_iteration"` path
- [ ] BiCGSTAB / GMRES path tests pass without the documented
      inconsistencies (confirmed via the reciprocity invariant test
      from Issue 11)

**Files**:

- rewrite: `orpheus/sn/solver.py`
- simplify: `orpheus/sn/operator.py`

**Design notes**: Issue that makes the architectural campaign visible
at the API surface. After it lands, `SNSolver` is a thin coordinator
and heavy lifting is in compositional operators. Closing of #96 / #97
is a major correctness win.

---

### Phase 5 — Tests + documentation

#### Issue 16: Behavioral regression suite

- **Module**: `module:tests`
- **Type**: `type:test`
- **Phase**: 5 (executed first)
- **Depends on**: none — this issue runs BEFORE the reshape begins
- **Complexity**: M

**Context**: Cardinal Rule 1 (correctness is critical) means the reshape
must not change physics. A frozen-reference regression suite is the
gating contract for every subsequent issue.

**Acceptance criteria**:

- [ ] Snapshot scalar flux + k_eff outputs from a representative
      pre-reshape test set (slab, sphere, cylinder; multigroup;
      Pℓ scattering at orders 0, 1, 3; vacuum + reflective BCs) into
      a frozen reference at `tests/_regression/sn_reference/`
- [ ] Pytest fixture compares post-reshape output to frozen reference
      at machine precision when `cell_update=DiamondDifference`
- [ ] Tagged `@pytest.mark.l1` and `@pytest.mark.regression`
- [ ] CI gates on regression for any PR touching `orpheus/sn/`,
      `orpheus/geometry/`, `orpheus/numerics/`

**Files**:

- new: `tests/_regression/sn_reference/` (snapshot data)
- new: `tests/sn/test_regression.py`

**Design notes**: **Generate snapshots ON main BEFORE branching for the
reshape work.** This is the very first action — Issue 16 lands on main,
then the campaign branches off. Otherwise the reference drifts with the
reshape.

---

#### Issue 17: Symmetry-preservation and reciprocity invariant tests

- **Module**: `module:tests`
- **Type**: `type:test`
- **Phase**: 5
- **Depends on**: Issues 11, 12
- **Complexity**: M

**Context**: New tests exercising structural properties the operator
algebra is supposed to guarantee — properties that weren't expressible
before because the abstractions weren't there.

**Acceptance criteria**:

- [ ] G-symmetry preservation: G-symmetric problem + G-invariant
      quadrature → discrete solution exactly G-invariant (to round-off).
      Concrete tests for slab(Z2), sphere(SO(3)), cylinder(SO(2)×R)
- [ ] Reciprocity: forward × adjoint detector responses match across
      source-detector swap. Uses `apply_transpose` on
      `SNStreamingOperator`
- [ ] Conservation: per-cell balance to round-off, verified against
      `derivations/discrete/sn/balance/`
- [ ] Capability honesty: every `LinearOperator` in the codebase
      correctly reports what it can do (auto-detected by attempting
      `solve` / `apply_transpose` and checking failure mode)
- [ ] Tagged `@pytest.mark.l1`, `@pytest.mark.verifies(...)` linked to
      Sphinx equation labels per V&V harness convention

**Files**:

- new: `tests/sn/test_invariants.py`

**Design notes**: Where the architectural rigor pays off — invariants
that were not previously expressible become testable.

---

#### Issue 18: Sphinx documentation for the operator algebra architecture

- **Module**: `module:docs`
- **Type**: `type:docs`
- **Phase**: 5
- **Depends on**: All implementation issues (lands progressively)
- **Complexity**: L

**Context**: Cardinal Rule 3 — Sphinx is the LLM's brain. The reshape
introduces a new architectural narrative (operator algebra unifying
SN/MoC/CP/diffusion) that needs a dedicated theory page, plus updates
to existing pages. Documentation IS the bridge to the AI-assisted
methodology presentation.

**Acceptance criteria**:

- [ ] New `docs/theory/operator_algebra.rst` — unifying narrative:
      differential / integral / variational forms; L/S/F as
      `LinearOperator`s; sweep as L⁻¹; how this scales to MoC/CP/diffusion
- [ ] New `docs/theory/discrete_measures.rst` — quadratures as discrete
      measures; tensor product / pushforward / restrict; symmetry-V
      framework; registry rationale with selection examples
- [ ] Update `docs/theory/discrete_ordinates.rst` — Key Facts, equations,
      how new architecture maps onto existing math; reference
      Bailey et al. (2009) where present
- [ ] Update `docs/theory/geometry.rst` — orbit-space classification,
      reduced operators, connection coefficients, where they live in
      code
- [ ] All cross-referenced via Nexus per Cardinal Rule 3
- [ ] Sphinx builds clean

**Files**:

- new: `docs/theory/operator_algebra.rst`,
       `docs/theory/discrete_measures.rst`
- modify: `docs/theory/discrete_ordinates.rst`,
          `docs/theory/geometry.rst`

**Design notes**: The math from the planning conversation thread
(orbit-space classification, fiber bundle for rays,
Galerkin = multigroup = PN = FE structure, quadratures as discrete
measures with G–V duality) is the substantive content. Use it.
Documentation campaign overlaps the last few implementation waves.

---

## Dependency waves (parallelization plan)

```
Wave A (parallel, no deps): Issue 1, Issue 2, Issue 16
Wave B (parallel after A): Issue 3, Issue 4 [needs 2,3], Issue 5 [needs 3,4]
                           Issue 6 [needs 1], Issue 7 [needs 4]
Wave C (after B):          Issue 8 [needs 6]
Wave D (parallel after C): Issue 9 [needs 8], Issue 10 [needs 6]
Wave E (after D):          Issue 11 [needs 1,8,9,10], Issue 13 [needs 1]
Wave F (after E):          Issue 12 [needs 6,8,10]
Wave G (parallel after E): Issue 14 [needs 1]
Wave H (after F+G):        Issue 15 [needs 11,13,14] — closes #96, #97
Continuous:                Issue 17 [after 11,12], Issue 18 [progressive]
```

**Issue 16 lands on `main` BEFORE the campaign branch is created.**

Realistic single-developer sequencing with Claude Code: ~10–14 sessions
for the reshape, plus documentation campaign overlapping the last few
waves.

## Conventions for sub-agent issue creation

When the sub-agent translates these specs to GitHub issues:

1. **Title format**: short imperative, e.g. `numerics: LinearOperator
   protocol with capability tags`
2. **Labels**: copy `module:` and `type:` labels from each spec
3. **Body**: copy Context + Acceptance Criteria + Files + Design Notes
   verbatim. Sub-agents and fresh Claude Code sessions read the body.
4. **Cross-references**: use `Depends on #NN` and `Closes #NN`
   (latter for issues 96 / 97, closed by Issue 15)
5. **Milestone**: `SN Reshape (Wave 1)` — create this milestone first
6. **Branch convention**: each issue gets a topic branch
   `<type>/<scope>/<short-description>` per the project's git workflow
   (e.g. `feat/numerics/linear-operator-protocol`)
7. **Commit convention**: Conventional Commits
   `<type>(<scope>): <summary>` per `CLAUDE.md` git workflow section
8. **Issue close**: `Closes #NN` trailer in commit body; PR uses
   `git merge --ff-only`

## Long-term migration sequence (after this campaign)

After SN reshape merges to main, subsequent waves consume the same
primitives. Each wave is a separate campaign with its own plan
document under `.claude/plans/`.

| Wave | Solver | Estimated complexity | Key consumer of new primitives |
|------|--------|---------------------|-------------------------------|
| 2 | MoC | L | `BundleMeasure`, ray-quadrature pushforward, `LinearOperator` |
| 3 | CP | L | dense `LinearOperator`, scalar-flux ops, BC tensor decomp |
| 4 | MC | M | `Geometry`, `BoundaryCondition`, sampling from `DiscreteMeasure` |
| 5 | Diffusion | S | `LinearOperator` (already FD-based), retire `EigenvalueSolver` Protocol |

After Wave 5: `numerics/eigenvalue.py:power_iteration` and the old
`EigenvalueSolver` Protocol can be deleted.

## Risk register

- **Behavioral regression slip**: catch by Issue 16 before branching
- **Curvilinear cell-update API insufficient**: surface during Issue 9
  implementation; if so, revise Issue 8 ABC and reconsider downstream
  issues
- **MoC/CP primitive blockers discovered late**: mitigated by Issue 2
  shipping `BundleMeasure` upfront; if other shape mismatches surface
  during MoC reshape, those become primitive revisions in Wave 2
- **Performance regression** in 1D cumprod: Issue 12 keeps the algebraic
  optimization; benchmark in regression suite to catch slowdowns
- **Pℓ anisotropic correctness drift**: behavioral regression suite
  includes Pℓ tests at orders 0, 1, 3 to gate this

## References

- `CLAUDE.md` — project-level Cardinal Rules and conventions
- `docs/theory/discrete_ordinates.rst` — current SN theory page
- `derivations/discrete/sn/balance/` — symbolic derivation of cell
  balance equations
- Bailey et al. (2009) — `_setup_spherical` / `_setup_cylindrical`
  references for angular redistribution and Morel-Montry closure
- Lewis & Miller, Computational Methods of Neutron Transport (1984) —
  §4.5 (M-M closure for RZ), §5.3 (LD discretization)
- GitHub issues #96, #97 — operator-sweep inconsistency, closed by
  Issue 15
