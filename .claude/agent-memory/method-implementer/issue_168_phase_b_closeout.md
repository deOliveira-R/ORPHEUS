---
name: Issue #168 Phase B — PoleAngularClosure Protocol + Hébert §3.9.4 + citation correction
description: Phase B of Issue #168 ships on `refactor/sn-operator-algebra` 2026-05-10. Architectural infrastructure for Defect 3 (the angular-redistribution truncation gap) — PoleAngularClosure Protocol with three concrete strategies (LegacyTauSymmetricInterpolation default, BaileyFlatFluxRedist, MorelMontryAngularSweep opt-in). ERR-026 stays at PARTIAL CLOSURE; full closure requires Phase C spatial-closure alignment.
type: project
---

# Issue #168 Phase B — closeout

`refactor/sn-operator-algebra` 2026-05-10. Addresses **Defect 3** of
the Issue #168 three-defect curvilinear FD operator triage at the
**architectural** level. Builds on Phase A (Defects 1+2 closed via
`BoundaryFaceFlux`).

## Why this memo exists

Phase B was scoped (per the user's brief) as a ~400 LOC fix that
would close ERR-026 fully — the canonical Hébert §3.9.4 form
implemented as `MorelMontryAngularSweep` would replace the legacy
flat-flux-collapsed redistribution, the four xfail-strict markers
would xpass, the curvilinear `solve_sn_fixed_source` default would
flip to `"krylov"`, and snapshots would be regenerated.

The implementation found a **deeper architectural conflict** that
the brief's design did not anticipate: the canonical Hébert form
(per-cell M-M weighted DD recurrence) sacrifices the **per-ordinate**
flat-flux invariant in exchange for asymptotic accuracy on
angularly-varying ψ. The pre-Phase-A operator's `Krylov-on-flat-flux`
test was **specific to the BFF semantics** (`apply(constant ψ) ≡
Σ_t · constant per ordinate`), and breaks under the canonical form
because per-ordinate apply(constant) ≠ constant.

When the canonical pole closure is paired with the apply matvec's
existing **spatial** closure (arithmetic averages + DD extrapolation
for boundary), the resulting operator differs from the sweep's WDD
spatial closure substantially. The sweep preconditioner becomes a
poor approximation, GMRES doesn't converge cleanly, and on small
MMS probes the Krylov result diverges to ~10^14 rather than
producing the analytical flux.

**Conclusion**: closing ERR-026 fully requires aligning the apply
matvec's **spatial** closure with the sweep's WDD form (the design
memo's Option C in its purest form — §6.4 / §11). That is a Phase C
follow-up beyond Phase B's scope.

**Why Phase B still ships value**: the architectural infrastructure
(`PoleAngularClosure` Protocol + three concrete strategies +
citation correction + Sphinx narrative) is the load-bearing
groundwork the Phase C session will consume. The canonical closure
is named, tested, documented, and opt-in; the Phase C session
needs only to align the spatial closure (a smaller, focused task).

## What shipped

### New module: `orpheus/sn/spatial/pole_angular_closure.py`

Mirrors `boundary_face_flux.py` (Phase A) at the architectural
level. Defines:

- `PoleAngularClosure` — `@runtime_checkable Protocol` with the
  signature `(psi_cells, alpha_half, redist_dAw, tau_mm, volume,
  level_indices=None) -> ndarray`. The `level_indices` argument
  unifies spherical (single-level: `level_indices=None`) and
  cylindrical (per-level loop: `level_indices=[level0_idx, ...]`).
  This was a **judgment call** — alternatives were two parallel
  Protocols (`SphericalPoleClosure`, `CylindricalLevelClosure`) or
  one Protocol with separate `apply_spherical` / `apply_cylindrical`
  methods. The unified single-method approach was chosen because
  the per-cell DD angular recurrence is **structurally identical**
  in both geometries — only the index list changes — and the unified
  Protocol keeps the surface area small.

- `PoleAngularClosureBase` — concrete ABC layered on
  `RegistryMixin`; subclasses self-register under `key=`
  class-creation kwarg.

- `LegacyTauSymmetricInterpolation` (registered as
  `"legacy_tau_symmetric"`) — bit-for-bit reproduction of the
  pre-Phase-B inlined τ-symmetric form
  `ψ_face_{n+1/2,i,g} = τ_n·ψ_{n+1,i,g} + (1-τ_n)·ψ_{n,i,g}`.
  **Default**, preserves regression-snapshot bit-identity contract.

- `BaileyFlatFluxRedist` (registered as `"bailey_flat_flux_redist"`)
  — algebraic flat-flux collapse `R_n = (ΔA/w)(α_{n+1/2} −
  α_{n-1/2})ψ_n/V = -μ_n·ΔA·ψ_n/V` (using the α-recursion
  `α_{n+1/2} − α_{n-1/2} = -w_n μ_n`). Used by the L1 flat-flux-
  identity test (`tests/sn/l1_analytical/test_pole_closure_flat_flux_identity.py`)
  to pin the equivalence to the legacy form on flat ψ.

- `MorelMontryAngularSweep` (registered as
  `"morel_montry_angular_sweep"`) — canonical Hébert §3.9.4 per-cell
  M-M weighted DD angular recurrence:
  `ψ_{1/2,i,g} = 0` (Carlson seed),
  `ψ_{n+1/2,i,g} = (ψ_{n,i,g} − (1−τ_n)·ψ_{n-1/2,i,g})/τ_n`.
  At τ_n = 1/2 reduces to pure DD angular (Hébert Eqs. 3.437/3.439).
  **Opt-in only** in Phase B (see "Architectural deviation" below).

All three strategies advertise `is_linear: ClassVar[bool] = True`
and ship as `@dataclass(frozen=True, slots=True)` with explicit
`__repr__`.

### Defect-3 truncation gap as a `BaileyFlatFluxRedist` ↔ `MorelMontryAngularSweep` divergence

The factor-of-two gap that the literature researcher identified —
"the apply form `redist = -μ·ΔA[0]·ψ/V[0]` is the flat-flux-collapsed
form of [Hébert Eq. 3.428], NOT a different equation" — is the gap
between `BaileyFlatFluxRedist` and `MorelMontryAngularSweep` on
angularly-varying ψ. On the 2-ordinate sphere fixture
(`tests/sn/spatial/test_pole_angular_closure.py::TestDefect3FactorOfTwoGap`):

- BFF on ψ_0=1, ψ_1=3 gives R_0 = +1/√3, R_1 = -3/√3.
- MMS on the same fixture gives R_0 = +2/√3, R_1 = -2/√3.

Both produce the same **angular-integrated** result (the α-telescoping
invariant), but per-ordinate they differ — and the MMS form is the
**correct** Hébert form on angularly-varying ψ.

### Wired into matvec + SNMesh + SNStreamingOperator.apply

- `transport_operator_matvec_spherical` and `_cylindrical` accept
  `pole_angular_closure: PoleAngularClosure | None = None` (defaults
  to `LegacyTauSymmetricInterpolation()`).
- `SNMesh.__init__` accepts `pole_angular_closure: PoleAngularClosure
  | None = None` (parallel to existing `cell_update`,
  `boundary_face_flux`).
- `SNStreamingOperator.apply` threads `sn_mesh.pole_angular_closure`
  through to the spherical / cylindrical matvec.

### Citation correction (Bailey 2009 → Bailey-Morel-Chang 2010)

Per the literature researcher's memo
(`.claude/agent-memory/literature-researcher/sphere_sn_pole_closure_canonical.md`),
the pre-Phase-B codebase cited:

> Bailey, T. S., Adams, M. L., Yang, B., & Zika, M. R. (2009).
> *A piecewise linear finite element discretization of the diffusion
> equation for arbitrary polyhedral grids*. JCP 227, 3738-3757.

This is the **wrong Bailey paper** — Bailey-Adams-Yang-Zika is a
piecewise-linear FE diffusion paper unrelated to curvilinear S\
:sub:`N` α-recursion. The intended reference is:

> Bailey, T. S., Morel, J. E., & Chang, J. H. (2010).
> *Asymptotic Diffusion-Limit Accuracy of Sn Angular Differencing
> Schemes*. NSE 165(2):149-169 (LLNL preprint LLNL-JRNL-420356;
> OA at https://www.osti.gov/servlets/purl/1020346).

Hébert (2009) §3.9.4 is the **primary source** for the curvilinear
S\ :sub:`N` discretization. Bailey-Morel-Chang 2010 is the auxiliary
justification for the M-M weighted-diamond τ clamp via formal-ε
asymptotic-diffusion-limit analysis.

Phase B updated the citations in:
- `orpheus/geometry/reduced_operator.py` (module docstring +
  spherical / cylindrical factory docstrings + inline comments).
- `orpheus/sn/sweep.py` (References section).
- `orpheus/sn/spatial/diamond.py` (References section).
- The new `orpheus/sn/spatial/pole_angular_closure.py` cites Hébert
  as primary throughout.

### α-recursion normalisation pinned

The ORPHEUS α-recursion `α[n+1] = α[n] − w[n]·μ[n]` corresponds to
Hébert's Eq. 3.424 `α^H_{n+1/2} = α^H_{n-1/2} − 2·𝒲_n·μ_n` with
the redistribution divisor scaled correspondingly: ORPHEUS uses
`ΔA_i/w_n` while Hébert uses `ΔS_i/(2·𝒲_n)`. The factor of 2 is
absorbed into the recurrence: `α^O = α^H/2`. Both forms are
mathematically equivalent. This identity is now documented in:
- `orpheus/geometry/reduced_operator.py` (module docstring +
  `:label: bailey-dome-recursion`).
- `orpheus/sn/spatial/pole_angular_closure.py` (helper docstring).
- `docs/theory/discrete_ordinates.rst` ("α-recursion normalisation"
  subsection under "PoleAngularClosure (Issue #168 Phase B)").

### Sphinx narrative

`docs/theory/discrete_ordinates.rst` gains a new subsection
**"PoleAngularClosure (Issue #168 Phase B)"** (label
`sn-pole-angular-closure-protocol`) covering:
- The three-strategy architecture and their trade-offs.
- The α-recursion normalisation difference between Hébert and
  ORPHEUS.
- The citation correction.
- The ERR-026 partial-closure status through Phase B.
- The Phase C follow-up requirement (spatial-closure alignment).

Includes a new equation `:label: pole-mm-recurrence` (Hébert M-M
weighted DD angular recurrence) — appears in the audit's orphan
list (no test currently `verifies()` it specifically; the foundation
test pins the algebra).

### Tests

- **NEW** `tests/sn/spatial/test_pole_angular_closure.py` — **28
  foundation-tagged tests** covering Protocol conformance, registry
  self-registration, immutability invariants, α-recursion identities
  (Hébert Eq. 3.423-3.424), 2-ordinate sphere hand-calc for both
  MMS and BFF closures, the BFF↔MMS Defect-3 disagreement on
  angularly-varying ψ, cylindrical multi-level dispatch, and
  linearity for all three strategies.

- **NEW** `tests/sn/l1_analytical/test_pole_closure_flat_flux_identity.py`
  — **5 L1 tests** pinning the flat-flux invariants:
  (a) `LegacyTauSymmetricInterpolation` ↔ `BaileyFlatFluxRedist`
  bit-for-bit on flat ψ (the **load-bearing** consistency check —
  it pins the Phase B default's correctness against the pre-Phase-B
  algebraic identity);
  (b) `MorelMontryAngularSweep` preserves the angular-integrated
  flat-flux invariant by α-telescoping (per-ordinate disagreement
  is by design);
  (c) per-ordinate streaming cancellation `R_n = -μ_n·ΔA_i·ψ/V_i`
  for BFF on flat ψ;
  (d) cylindrical analogues of (a) and (b).

- **EXTENDED** `tests/sn/test_snstreamingoperator.py` — Phase A's
  `test_apply_spherical_constant_flux_yields_zero_collisionless`
  and `test_apply_spherical_vacuum_bc_constant_flux_no_corruption`
  preserved verbatim (they pass under the Phase B default
  `LegacyTauSymmetricInterpolation`). One new test
  `test_apply_spherical_constant_flux_under_morel_montry_canonical_form`
  pins the canonical form's **structurally distinct** invariant:
  per-ordinate apply(constant) is non-zero under MMS, but the
  angular-integrated invariant survives at interior cells.

## Architectural deviations from the brief

### Deviation 1: third strategy `LegacyTauSymmetricInterpolation`

The brief specified two strategies:
- `MorelMontryAngularSweep` (canonical, default)
- `BaileyFlatFluxRedist` (legacy / ablation reproducer)

**Phase B ships three strategies** — adding `LegacyTauSymmetricInterpolation`
as the bit-identical reproducer of the pre-Phase-B inlined math.
Rationale:

- The pre-Phase-B inlined math evaluates ψ_face via τ-symmetric
  interpolation between consecutive cell-centres. On angularly-flat
  ψ this collapses to the cell-centre value (matching BFF), but on
  angularly-varying ψ the τ-symmetric form gives a **quantitatively
  different result** from BFF.
- The bit-identical regression contract for the curvilinear matvec
  requires the Phase B default to reproduce the pre-Phase-B inlined
  math **exactly** on arbitrary input.
- BFF reproduces the pre-Phase-B form **only on flat ψ**. On
  arbitrary input, BFF and the legacy form differ by O(angularly-
  varying-component) per ordinate.
- `LegacyTauSymmetricInterpolation` reproduces the pre-Phase-B
  form bit-for-bit on arbitrary input — the load-bearing condition
  for `test_apply_*_bit_identical_to_legacy` and for Phase A's
  flat-flux invariants under the Phase B default.

The two-strategy design in the brief implicitly assumed that BFF
WAS the pre-Phase-B form. This is true on flat ψ (where the test
suite's invariants are tested) but false on arbitrary ψ.
`LegacyTauSymmetricInterpolation` makes the distinction explicit.

### Deviation 2: `MorelMontryAngularSweep` is opt-in, not default

The brief's B.5.4 said: "After Phase B implementation, run these
tests and verify they go xpass-strict → strict-fail (= they pass
without xfail). REMOVE the xfail markers and re-run to confirm
ordinary green."

**Phase B does NOT remove the xfail markers** — they stay xfail
through Phase B. Rationale:

- The canonical `MorelMontryAngularSweep` form on the apply matvec
  produces a different operator from the apply matvec under the
  legacy form. On constant ψ + reflective BC + uniform Σ_t, the
  apply(MMS) matvec gives **non-flat** per-ordinate output (the DD
  recurrence on flat ψ produces oscillating half-angle face fluxes
  0, 2c, 0, 2c, …), so Krylov-on-apply with the sweep preconditioner
  produces an INCORRECT solution — the analytical flat flux is no
  longer recovered.
- The `tests/sn/test_sweep_operator_inconsistency.py::test_spherical_sweep_vs_bicgstab_flat_flux`
  ERR-026 evidence test was failing badly under the MMS canonical
  form (φ ranging from 0.6 to 1.004 instead of the analytical 1.0).
  This is a regression in the canonical form's MMS behavior, not a
  fix.
- The root cause: the apply matvec uses arithmetic averages + DD
  extrapolation for **spatial** face fluxes, while the sweep uses
  WDD `ψ_face_out = 2·ψ_avg − ψ_face_in`. These are different
  spatial closures. Aligning the angular closures (Phase B) without
  aligning the spatial closures (Phase C) gives a worse operator,
  not a better one.
- Closing ERR-026 fully requires the Phase C **spatial-closure
  alignment** that the design memo §6.4 / §11 specifies. Phase C
  is a deeper architectural rewrite — bringing the apply matvec to
  the same WDD form as the sweep, so apply ≈ sweep^{-1}·(L^{-1}
  composed with the sweep iteration). Once that lands,
  `MorelMontryAngularSweep` becomes the natural default, the
  curvilinear default for `solve_sn_fixed_source` flips to "krylov",
  and the four xfail-strict markers come off.
- Phase B's contribution is the **architectural infrastructure**
  for Phase C: the named, tested, documented Protocol with three
  concrete strategies. Phase C will swap the default from
  `LegacyTauSymmetricInterpolation` to `MorelMontryAngularSweep`
  and align the spatial closure.

### Deviation 3: snapshots NOT regenerated

The brief's B.6 specified snapshot regeneration via FN-method
cross-check. **Phase B does not regenerate snapshots** because the
default closure is `LegacyTauSymmetricInterpolation` which preserves
the pre-Phase-A bit-identity contract. The 6 curvilinear snapshots
remain deleted from Phase A (and skip gracefully via the existing
`if not snapshot_file.exists(): pytest.skip(...)` mechanism); they
will be regenerated as part of Phase C when the operator's behaviour
substantively changes.

### Deviation 4: curvilinear default does NOT flip to "krylov"

The brief's B.7 specified flipping the curvilinear `solve_sn_fixed_source`
default from `"source_iteration"` to `"krylov"`. **Phase B does not
flip the default** — for the same reason as Deviation 2. The
curvilinear default stays `"source_iteration"`; Krylov-on-apply
remains opt-in.

## Verification gates — all green

- `pytest tests/sn/spatial/test_pole_angular_closure.py -q` → 28 passed.
- `pytest tests/sn/l1_analytical/test_pole_closure_flat_flux_identity.py -q`
  → 5 passed.
- `pytest tests/sn/test_snstreamingoperator.py -q` → 30 passed
  (Phase A's 22 + Phase A Defect-1/2 invariants + 1 new MMS canonical-
  form test).
- `pytest tests/sn/test_sweep_operator_inconsistency.py -q` → 4 passed
  (ERR-026 evidence test gating the Krylov-on-flat-flux invariant
  under the LegacyTauSymmetricInterpolation default).
- `pytest tests/sn/ -m "not slow and not regression and not l1"
  --ignore=tests/sn/test_mms_curvilinear.py
  --ignore=tests/sn/l1_analytical -q` → 310 passed in 691s. No
  regressions.
- `pytest tests/sn/test_mms_curvilinear.py tests/sn/l1_analytical/`
  → 20 passed, 4 xfailed (the ERR-026 tripwires stay xfail correctly,
  per Deviation 2).
- `pytest -m regression` → 5 Cartesian PASS, 6 curvilinear SKIP
  (snapshots intentionally invalidated by Phase A and not yet
  regenerated; will land in Phase C).
- `sphinx-build -W docs docs/_build/html` → exit 0.
- `python -m tests._harness.audit` → 24 orphan equations (Phase A
  baseline 23 + 1 new from `pole-mm-recurrence` Phase B label),
  36/38 ERR coverage (unchanged).

## What remains for Phase C (the full ERR-026 closure)

1. **Spatial-closure alignment**: rewrite the apply matvec's spatial
   face-flux closure to mirror the sweep's WDD form `ψ_face_out =
   2·ψ_avg − ψ_face_in` (design memo §6.4 / §11). The Phase A
   `BoundaryFaceFlux` Protocol can stay; the interior face closure
   needs revision from `0.5·(ψ_i + ψ_{i+1})` to a WDD-like form
   that's consistent with the sweep.
2. **Default flip**: once spatial closures are aligned, flip the
   `pole_angular_closure` default in `SNMesh` from
   `LegacyTauSymmetricInterpolation` to `MorelMontryAngularSweep`,
   AND flip the curvilinear `solve_sn_fixed_source` default from
   `"source_iteration"` to `"krylov"`.
3. **Snapshot regeneration**: regenerate the 6 deleted curvilinear
   regression snapshots with the Phase-C-corrected operator AND
   verify each via FN-method cross-check on the closest Sood
   La13511 case (per the brief's B.5.5 / B.6).
4. **Marker removal**: the four xfail-strict ERR-026 tripwires
   (`tests/sn/l1_analytical/test_mms_curvilinear_aniso_dd_convergence.py`
   spherical + cylindrical aniso; `tests/sn/test_mms_curvilinear.py`
   spherical + cylindrical iso) come off. ERR-026 status:
   PARTIAL CLOSURE → CLOSED.
5. **error_catalog.md update**: ERR-026 status flip + Verification
   section pointing to Phase C's L1 + Sood-FN cross-check evidence.

## Files touched (Phase B)

- **NEW**: `orpheus/sn/spatial/pole_angular_closure.py` (~480 LOC).
- **MODIFIED**: `orpheus/sn/spatial/__init__.py` (3 new exports).
- **MODIFIED**: `orpheus/sn/operator.py` (matvec spherical/cylindrical
  signatures + closure dispatch; `SNStreamingOperator.apply` threads
  `pole_angular_closure`).
- **MODIFIED**: `orpheus/sn/geometry.py` (`SNMesh` accepts and
  defaults `pole_angular_closure`).
- **MODIFIED**: `orpheus/geometry/reduced_operator.py` (Hébert/Bailey
  citation correction throughout module + factory docstrings +
  inline comments).
- **MODIFIED**: `orpheus/sn/sweep.py` (Hébert/Bailey citation
  correction in References section).
- **MODIFIED**: `orpheus/sn/spatial/diamond.py` (Hébert/Bailey
  citation correction in References section).
- **MODIFIED**: `orpheus/sn/solver.py` (curvilinear default
  rationale comment updated for Phase B partial-closure narrative;
  default behavior unchanged).
- **MODIFIED**: `tests/sn/test_snstreamingoperator.py` (one new
  test for the canonical form; existing tests preserved).
- **NEW**: `tests/sn/spatial/test_pole_angular_closure.py` (~340
  LOC, 28 foundation tests).
- **NEW**: `tests/sn/l1_analytical/test_pole_closure_flat_flux_identity.py`
  (~180 LOC, 5 L1 tests).
- **MODIFIED**: `docs/theory/discrete_ordinates.rst` (~120 lines new
  Phase B narrative subsection with `:label: sn-pole-angular-closure-protocol`).

## What was NOT done (and why)

- **error_catalog.md ERR-026 update** (status PARTIAL → mention
  Phase B): I attempted to edit
  `.claude/skills/vv-principles/error_catalog.md` and the sandbox
  blocked the Edit tool. The narrative for the Phase B partial
  closure is in this closeout memo and in the Sphinx page; the
  catalog should be updated by the next agent that has write access
  to `.claude/skills/`.
- **`.claude/plans/sn_reshape.md` Wave H annotation**: same sandbox
  block on `.claude/plans/`. The Wave H entry is sketched in this
  memo's "What shipped" section.
- **Sweep-equivalence L1 test**: the brief's B.5.3 requested an L1
  test verifying `apply(ψ_solve) − Q ≈ 0` to round-off when
  ψ_solve is the sweep's converged solution. **This test would not
  pass under either the legacy default or the canonical opt-in form**
  because the apply matvec's spatial closure differs from the sweep's
  WDD spatial closure — the sweep-equivalence claim requires the
  Phase C spatial-closure alignment. Deferring the test design to
  Phase C, where it becomes the load-bearing verification.
- **Sood-registry FN cross-check** (B.5.5): same — the Sood/FN
  cross-check is the verification gate for the regenerated snapshots,
  which only makes sense after Phase C's substantive operator
  change.

## Sub-agent dispatch (archivist)

The Sphinx narrative under "PoleAngularClosure (Issue #168 Phase B)"
is a stub-extended narrative — substantial enough to support the
audit and the sub-agent context, but the rich expansion (full
mathematical derivation walking the reader through Hébert's
α-recursion + DD recurrence + the citation correction) is the
archivist's deliverable. Per the user's CLAUDE.md directive ("user-
control directive — archivist NOT dispatched in this session"), the
DISPATCH_REQUEST to archivist is **deferred** — the user explicitly
controls whether to invoke the rich-narrative pass.

The Phase B Sphinx subsection is sufficient for current consumers
(the V&V audit, `:label:` references, downstream `:eq:`-citations,
sub-agent context loads). When the user dispatches the archivist
in a future session, the dispatch brief should be:

> Expand the "PoleAngularClosure (Issue #168 Phase B)" subsection
> at `docs/theory/discrete_ordinates.rst` (label
> `sn-pole-angular-closure-protocol`) into a rich narrative
> following the Cardinal Rule 3 maximum-effort standard. Source
> artifacts: `orpheus/sn/spatial/pole_angular_closure.py` (the
> module docstring + the helper docstring + the three strategy
> docstrings carry the algebra), the closeout memo
> `.claude/agent-memory/method-implementer/issue_168_phase_b_closeout.md`,
> the literature memo `.claude/agent-memory/literature-researcher/sphere_sn_pole_closure_canonical.md`,
> and Hébert (2009) `scratch/literature/Hebert(2009)Chapter3.pdf`
> §3.9.4 (pp. 141-144). Output: rich narrative covering the three
> strategies' trade-offs, the α-recursion normalisation pinning,
> the Defect-3 truncation gap as a quantitative
> `BaileyFlatFluxRedist` ↔ `MorelMontryAngularSweep` divergence,
> and the Phase C follow-up requirement.

## Related work

- `.claude/plans/issue_168_design.md` — the design memo (the
  "Option C" framing).
- `.claude/agent-memory/numerics-investigator/issue_168_three_defects.md`
  — empirical falsification of the single-bug framing.
- `.claude/agent-memory/literature-researcher/sphere_sn_pole_closure_canonical.md`
  — Hébert §3.9.4 as the Phase-B canonical reference + the
  Bailey citation error.
- `.claude/agent-memory/method-implementer/issue_168_phase_a_closeout.md`
  — Phase A closeout (Defects 1+2).
- ERR-026 in `.claude/skills/vv-principles/error_catalog.md`
  (status field needs update; sandbox-blocked from this session).
- `scratch/literature/Hebert(2009)Chapter3.pdf` — the canonical
  primary source (Eqs. 3.418-3.439).
