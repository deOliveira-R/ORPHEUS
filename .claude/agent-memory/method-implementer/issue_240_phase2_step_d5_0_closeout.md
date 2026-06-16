---
name: issue-240-phase2-step-d5-0-closeout
description: #240 Phase 2 Step D5-0 routing honesty — mint the scheme trait transverse_coupling_is_facewise (DD True / LD False default) and narrow ScanMarch.supports d≥2 arm to read it, closing the CONFIRMED-LIVE 2-D Cartesian LD→ScanMarch misroute (inline DD silently dropping LD's bilinear slope). Bit-identical for every exercised path: the strict gate's pre-existing set stays 505/1/4 (only +7 new D5-0 tests). Honest interim state: 2-D LD now routes to MovingFrontierWindow whose LD cell_kernel_batch raises the d=1-only NotImplementedError (D5b closes it).
metadata:
  type: project
---

# #240 Phase 2 Step D5-0 — routing honesty (the 2-D LD misroute fix)

**Branch** `feature/sn-space-angle-tier2` 2026-06-16 (NOT committed — main
agent runs elegance-enforcer + qa on the diff, then commits). Host env,
`.venv/bin/python -O -m pytest`.

## What was wrong (CONFIRMED LIVE pre-fix, probed this session)

`ScanMarch.supports(2-D Cartesian LD mesh).ok == True` →
`default_for(2-D LD) == ScanMarch`. ScanMarch's row-march interior runs
**inline DiamondDifference with no scheme dispatch**, so a 2-D LD mesh
**silently computed DD, dropping LD's bilinear slope**. The gate conflated
`is_affine_scannable` (a SINGLE-axis prefix-scannability claim — LD: `True`)
with the CROSS-axis separability the `d≥2` row-march actually requires
(LD: `False`, because LD's multi-D closure is bilinear — a 1st-order slope
moment per axis, not a 0th-order face value).

## The architecture (SETTLED — cross-domain-attacker Frame 1, brief-mandated)

Minted a SCHEME-named trait (NOT strategy-named — `is_scan_march_compatible`
was the test-architect's SUPERSEDED proposal; naming a scheme property after
one consuming strategy is a frame-leak — see
`.claude/agent-memory/cross-domain-attacker/d5_trait_and_mms_frames.md`
Frame 1):

```python
transverse_coupling_is_facewise: ClassVar[bool] = False  # on DiscretizationSchemeBase
```

- **Default `False`** (conservative opt-in — a scheme must DECLARE separable/
  facewise transverse coupling; a new scheme that forgets it is safely
  excluded from the d≥2 scan-march and falls back to the wavefront).
- **`True`** on `DiamondDifference` (its transverse term is the 0th-order face
  value `s_y·ψ_{y,in}` folded into the scan source → tensor-product separable).
- **`False` (default)** on `LinearDiscontinuous` (one-line docstring note on its
  trait block: bilinear per-axis slope → d-D closure non-separable → rides the
  DAG wavefront in d≥2, #38/D5b).
- `Step` is a docstring stub only (`scheme.py`'s `class Step`); NOT created.
  Base docstring notes Step will be `True` when built (slopeless upwind shares
  DD's facewise structure).
- Source cited in the trait docstring: Maginot–Ragusa–Morel 2016 (UBLD) /
  Adams 2001 (LD's irreducible multi-D coupling); DD/box facewise separability
  = standard tensor-product central-difference (Lewis & Miller §§4.5, 8).

The trait was ALSO added to the `@runtime_checkable DiscretizationScheme`
Protocol attribute set (`is_linear`/`is_positivity_preserving`/
`is_affine_scannable`/`transverse_coupling_is_facewise`) — this is what made
the synthetic-strategy conformance tests (below) require the new declaration.

## The supports() narrowing

`ScanMarch.supports` (`loss_representation.py`) split into two arms (was one
conflated boolean):

- **1-D arm UNCHANGED**: `mesh.scheme.is_affine_scannable` (LD's 1-D scan IS
  valid on the scan-march; just not the default — `CumprodScan` precedes it).
- **d≥2 arm**: `mesh.is_cartesian and mesh.ndim == 2 and
  mesh.scheme.transverse_coupling_is_facewise` (kept `ndim==2` — 3-D row-march
  is future scope #227). Reason string is pedagogically load-bearing: "2-D
  scan-march requires a scheme whose transverse coupling is facewise
  (separable…) — …; Linear-Discontinuous's bilinear slope coupling needs the
  wavefront…".

**Verified honest interim state** (probed live, post-fix):
`default_for(2-D LD)` falls through `CumprodScan`(1-D→no) → `ScanMarch`(now→no)
→ `MovingFrontierWindow` (`_DAGWavefront.supports`: `is_cartesian and ndim==2`
→ yes), whose LD sweep then RAISES the existing
`NotImplementedError("…supports d=1 (slab/1-D) only; got d=2…")` from
`LinearDiscontinuous._kernel_terms` (`linear_discontinuous.py:431`, the brief
said `:430`; the message is the loud d=1-only signal). Silent-DD-wrong →
LOUD not-yet-implemented. This is the correct interim state; D5b closes the
raise by supplying the multi-D bilinear LD kernel.

## Files touched (give to elegance-enforcer + qa)

PRODUCTION:
- `orpheus/sn/spatial/scheme.py` — Protocol attr-block + docstring
  (`:351` is_affine_scannable note + NEW `transverse_coupling_is_facewise`
  attr ~`:363`); `DiscretizationScheme` declares the 4th attr ~`:392`; Base
  ClassVar + full docstring inserted after the `is_affine_scannable` ClassVar
  (~`:695`, before `_registry_base`).
- `orpheus/sn/spatial/diamond.py` — NEW `transverse_coupling_is_facewise:
  ClassVar[bool] = True` + docstring, inserted after the `is_affine_scannable`
  ClassVar block (~`:143`, before `def update`).
- `orpheus/sn/spatial/linear_discontinuous.py` — one-line WHY note appended to
  the `is_affine_scannable` docstring block (`:279-291`); LD LEAVES the new
  trait at the `False` default (no new ClassVar — opt-out by inheritance).
- `orpheus/sn/loss_representation.py` — `ScanMarch.supports` (`:1225-1233` →
  split into the two-arm predicate).

TESTS:
- `tests/sn/sweep/core/test_unified_sweep_dispatch.py` — imports
  (`LinearDiscontinuous`, `MovingFrontierWindow`, `FullFieldWavefront`); NEW
  `_2d_ld_sn_mesh(nx=4, ny=3)` fixture (non-square, LS-S4, LD scheme);
  `TestD3SupportsMatrix._fake` migrated (added `facewise` kwarg →
  `scheme.transverse_coupling_is_facewise`); 4 NEW D5-0 routing tests +
  3 NEW `TestSchemeTraitProbe` strategy-free trait probes.
- `tests/sn/sweep/core/test_discretization_scheme_protocol.py` — 2 synthetic
  strategies (`IdentityDiscretizationScheme`, `FakeCurvilinearStrategy`) gained
  `transverse_coupling_is_facewise: ClassVar[bool] = False` (Protocol grew a
  required attr; conforming fixtures must declare it). NOT a value change —
  a contract-completeness fix; the 2 conformance `isinstance` tests turn green.

SPHINX STUB:
- `docs/theory/loss_representations.rst` — NEW `.. _loss-rep-scanmarch-facewise:`
  subsection (label + `:meth:`/`:class:` cross-refs to
  `ScanMarch.supports` / `DiscretizationSchemeBase` / `DiamondDifference` /
  `LinearDiscontinuous` + the test-gate class names + closeout-memo path +
  TODO marker). Flags the STALE inline `ScanMarch.supports` code block in
  `loss-rep-selection` (still shows the pre-D5-0 conflated predicate) +
  the `default_for` table (needs the scheme-dependence note) for D6.

## The trait's final name + docstring (one-line)

`transverse_coupling_is_facewise: ClassVar[bool] = False` — "does the multi-D
cell closure couple a NON-swept axis only through a face value (0th-order
trace), so the d-D update factors into independent per-axis 1-D scans chained
by scalar traces?" DD `True`, LD `False` (default). Distinct from
`is_affine_scannable` (single-axis); forward-reused by the diffusion ADI /
line-SOR preconditioner (#240's confirmed next consumer) — the reason it is
named for the SCHEME property, not the `ScanMarch` strategy.

## The supports() predicate (final)

```python
if mesh.is_1d:
    return Compatibility(mesh.scheme.is_affine_scannable, "…1-D…")
return Compatibility(
    mesh.is_cartesian and mesh.ndim == 2
    and mesh.scheme.transverse_coupling_is_facewise,
    "2-D scan-march requires a scheme whose transverse coupling is facewise…",
)
```

## New tests (7)

`TestD3SupportsMatrix` (routing honesty):
- `test_scan_march_refuses_2d_non_facewise_scheme_fake` — fake `facewise=False`
  refused + fake `facewise=True` admitted (the anti-pattern-#11 pair on the
  fake).
- `test_scan_march_refuses_2d_ld_real_mesh` — `ScanMarch.supports(2-D LD).ok
  is False` on a real `SNMesh(Mesh2D, LS-S4, scheme=LinearDiscontinuous())`.
- `test_2d_ld_default_for_routes_to_wavefront` — `default_for(2-D LD)` is
  `MovingFrontierWindow`/`FullFieldWavefront`, never `ScanMarch`.
- `test_2d_ld_sweep_raises_not_silently_dd` — `solve_sn_fixed_source(2-D LD)`
  RAISES `NotImplementedError` (the headline honesty claim: no silent DD).

`TestSchemeTraitProbe` (strategy-free, the cross-domain-attacker's
second-consumer discriminator — proves the trait is NOT strategy-entangled):
- `test_dd_reports_facewise_transverse_coupling` — DD `True` standalone.
- `test_ld_reports_non_facewise_transverse_coupling` — LD `False` standalone.
- `test_facewise_distinct_from_affine_scannable` — LD diverges (affine_scannable
  `True`, facewise `False`); DD coincides (both `True`). The split is OBSERVABLE
  only on a scheme where the two answers differ — LD.

All Mode-8 safe (`pytest.fail`, never bare `assert` — `-O` canonical).

## VERIFICATION (verbatim)

Baseline (pre-change), strict bit-identity gate:
```
505 passed, 1 skipped, 4 xfailed, 2 warnings in 85.47s (0:01:25)
```

Strict bit-identity gate POST-change
(`tests/sn/sweep/core tests/sn/solve -W "error::tests.sn.regression._regression_assert.DriftWarning"`):
```
512 passed, 1 skipped, 4 xfailed, 2 warnings in 85.55s (0:01:25)
```
512 = 505 pre-existing + 7 new D5-0 tests (which live INSIDE `tests/sn/sweep/
core`). Skip/xfail counts UNCHANGED (1/4). Proof the pre-existing 505 are
untouched — deselecting the 7 new tests:
```
505 passed, 1 skipped, 7 deselected, 4 xfailed, 2 warnings in 86.12s (0:01:26)
```
So the trait + supports change altered NO computed value on any exercised path
(bit-identical); the only delta is the 7 added tests.

Dispatch + spatial:
```
58 passed, 1 warning in 47.52s        (tests/sn/sweep/core/test_unified_sweep_dispatch.py tests/sn/spatial)
```
New D5-0 + trait-probe tests (all 7) PASSED (verbose run); the 4 migrated-fake
`TestD3SupportsMatrix` pins stay green.

LD green floor (spec baseline `24 passed, 1 xfailed`):
```
24 passed, 1 xfailed, 1 warning in 4.94s
```

Broader route-around set (`operators spatial sweep/core sweep/cartesian_2d
solve` with the spec `-k`):
```
1100 passed, 6 skipped, 7 deselected, 5 xfailed, 4 warnings in 147.34s (0:02:27)
```

V&V audit:
```
audit exit: 0
```

Sphinx `-W` build:
```
build succeeded.
sphinx exit: 0
```
(No warning references the new stub label `loss-rep-scanmarch-facewise`, the
changed scheme files, or the trait. Pre-existing `ld-cartesian-1d` / `ld-slab`
orphan-label warnings are unrelated — the LD theory section is still a D6
`.. todo:` stub; those labels are scheduled for D6 minting.)

## Lessons / skill-evidence

1. **Adding an attribute to a `@runtime_checkable` Protocol is a structural
   contract change that breaks `isinstance` conformance for EVERY structural
   implementor that doesn't declare it** — including synthetic test fixtures.
   In Python 3.12+ `@runtime_checkable` `isinstance` checks DATA-member
   presence, not just methods (confirmed: the 2 conformance tests went red with
   `assert False` until the fixtures declared the trait). This is the correct,
   expected blast radius of a Protocol-level trait mint — NOT a test bug. Worth
   a `coding-elegance` Pattern-4 note: "minting a trait on a `@runtime_checkable`
   Protocol forces every structural implementor (incl. test mocks) to declare
   it — that propagation IS the make-illegal-states-unrepresentable guarantee
   working; the conformance-test red is the guarantee firing, not a regression."
2. **The brief's "strict gate MUST stay 505/1/4" had an implicit assumption the
   new tests land OUTSIDE the gate's path** — they didn't (the dispatch test is
   in `tests/sn/sweep/core`). The correct generalization of the gate's intent is
   "the PRE-EXISTING set stays 505/1/4 and skip/xfail are unchanged" — verified
   by deselecting the new tests. A new-test-count delta on a directory-scoped
   bit-identity gate is expected; the load-bearing invariant is that no
   pre-existing test's value moved (here: literally bit-identical, the change is
   a routing predicate that touches no computed flux on any exercised path).

## Owed downstream

- **D6 archivist DISPATCH** (rich narrative): expand the
  `loss-rep-scanmarch-facewise` stub; refresh the STALE inline
  `ScanMarch.supports` code block + `default_for` outcome table in
  `loss-rep-selection`; thread the facewise/slope-wise transverse-coupling
  narrative (Maginot–Ragusa–Morel 2016 / Adams 2001) into the LD theory
  `.. todo:` stub on `discrete_ordinates`. (DISPATCH_REQUEST emitted with
  `followup: false`.)
- **D5a** (next, per the test-architect sequence): the N-D DD scan-march fold
  onto the coefficient model — makes ScanMarch scheme-generic so the D5-0 EXCL
  (LD never reaches inline DD) is enforced by the ABSENCE of inline DD, not just
  by `supports()`.
- **D5b** (the new math): the multi-D bilinear LD kernel closes the
  `linear_discontinuous.py:431` raise that D5-0 turned into the honest signal.
