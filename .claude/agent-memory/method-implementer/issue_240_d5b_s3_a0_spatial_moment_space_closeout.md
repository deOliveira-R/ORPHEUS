---
name: issue-240-d5b-s3-a0-spatial-moment-space-closeout
description: #240/#38/#37 D5b-S3-A0 — mint first-class SpatialMomentSpace (within-cell tensor-Legendre DG moment factor, the SPATIAL sibling of SphericalHarmonicSpace) + the OPTIONAL spatial-moment factor on the field-space factories (AngularField/ScalarField/HarmonicMomentField), gated "append iff per_axis>1" so DD/Step/LD all byte-identical at default. Unblocks the closeout's typed-field space-shape BLOCKER. Construct-general-only (capability default-OFF; NO production field selects the axis yet). DONE+VERIFIED. NOT committed.
metadata:
  type: project
---

# #240 D5b-S3-A0 — first-class SpatialMomentSpace + optional field-space factor

**Branch** `feature/sn-space-angle-tier2`. **NOT committed** (main agent reviews + commits).
Host env, `.venv/bin/python -O`. Canonical `python -O -m pytest`; NEVER all `tests/sn` (#212).

This is the typed-field-space half the prior dispatch
(`issue_240_d5b_s3_a_inc_c_closeout.md`) surfaced as a HARD PREREQUISITE
(its §BLOCKER). The USER chose the design: a **first-class typed space factor**
(NOT a raw `int` trailing axis), following the `SphericalHarmonicSpace` /
`TensorProductSpace` precedent. The scattering-lift half (the 3 einsum `...`
subscripts in `material_xs_field.py`) was ALREADY landed (uncommitted) by that
prior dispatch — LEFT UNTOUCHED here.

## WHAT LANDED (all DONE + VERIFIED)

### 1. NEW `SpatialMomentSpace` (the type)
`orpheus/numerics/spaces/spatial_moment_space.py` (NEW). A peer of
`SphericalHarmonicSpace`, the **spatial** moment carrier (within-cell
tensor-Legendre DG basis `{1, P₁}` per axis) vs the **angular** one (SH
over direction). Orthogonal axes — naming each distinct dispels the
"moment" collision.
- `@dataclass(frozen=True) class SpatialMomentSpace(FunctionSpace)`, fields
  `per_axis: int = 1`, `ndim: int = 0`; `shape == (per_axis ** ndim,)`
  checked in `__post_init__` (Pattern 4). Euclidean (`inner_product_weights
  is None` — the cell-mass metric lives on the UBLD operator, not the space,
  per #207 "units do not live on the space").
- Factory `from_per_axis(per_axis, ndim)` — the `SphericalHarmonicSpace.from_L`
  analogue (size in `shape`; `(per_axis, ndim)` descriptive metadata).
- Props `n_moments` (= `shape[0]`), `average_moment_index` (single-sourced
  from `_ubld.AVERAGE_MOMENT`, NOT a re-spelled `0`).
- `__eq__`/`__hash__` delegate to `FunctionSpace` (size-identity by
  `(name, shape)`; the `(per_axis,ndim)` factorisation is NOT in identity —
  `(pa=4,d=1)` == `(pa=2,d=2)`, both `(4,)`; they never coexist on one mesh).
- Module-level `spatial_moment_tail(n_cell_moments)` — the CELL analogue of
  `_ubld.face_moment_tail`, **DELEGATES to it** (Pattern 7 single-source of the
  "append iff > 1" policy: `() if n==1 else (n,)`).
- **CIRCULAR-IMPORT FIX (surprise vs precedent):** `_ubld` is imported LAZILY
  (deferred inside `spatial_moment_tail` + `_average_moment_index`), NOT at
  module level. A top-level `from orpheus.sn.spatial._ubld import ...` triggers
  `orpheus.sn.__init__` → `orpheus.transport` → `harmonic_moment_field` → back
  to this module (partially-initialised) → ImportError. The numerics layer must
  not depend on the SN package at IMPORT time; `_ubld` is leaf (numpy-only) so
  the deferred call-time import is cheap and still single-sources the constants.
  `SphericalHarmonicSpace` did NOT hit this (it imports `numerics.basis`, not SN).

### 2. NEW `TensorProductSpace.find_factor(factor_type)` (closes #207 gap)
`orpheus/numerics/space.py`. The `space.find_factor(SphericalHarmonicSpace).L`
query the `HarmonicMomentField` docstrings DOCUMENTED (#207) but that was NEVER
IMPLEMENTED — referenced in 3 docstrings, `grep`/Nexus found no method. Minted
now: returns the first factor `isinstance(f, factor_type)`, raises `KeyError`
if absent (Pattern 4 — the query is a structural assertion, not a silent None).
Verified `hmf.space.find_factor(SphericalHarmonicSpace).L` now WORKS (was a
latent broken claim). Both moment factors (SH + Spatial) are queryable by type.

### 3. OPTIONAL `spatial_moments` factor on the field-space factories
The 3 carriers from the closeout's carrier analysis:
- `AngularField` (un-windowed `AngularFlux` iterate) — `_space_for_mesh` /
  `from_mesh` / `zeros_on` gained `*, spatial_moments: int = 1`.
- `ScalarField` (the `ScalarSourceSink` scattering accumulator) — same.
- `HarmonicMomentField` (windowed iterate) — `from_mesh_and_L` /
  `zeros_for_mesh_and_L` gained `*, spatial_moments: int = 1`; a `spatial_moments`
  dataclass FIELD (default 1) so `_phase_space_shape` can compute the widened shape.
- Shared single-source `BulkField._compose_spatial_moments(space, mesh,
  spatial_moments_per_axis)` — composes `space * SpatialMomentSpace.from_per_axis(...)`
  via `*` EXACTLY as `from_mesh_and_L` composes the SH factor, gated `if
  spatial_moment_tail(per_axis**ndim) == (): return space`.
- **`_phase_space_shape` cross-check fix (the ripple):** `BulkField.__post_init__`
  validates `space.shape == _phase_space_shape()`. For Angular/Scalar I added
  `BulkField._spatial_moment_tail_of(space)` (reads the `SpatialMomentSpace`
  factor OFF the field's own space via `find_factor` → the space is the single
  source of truth for the moment width); `_phase_space_shape` appends it to the
  base `(N/ng,*spatial)` prefix. `HarmonicMomentField`/`MomentDisplacement` use
  their stored `spatial_moments` field instead (they already carry `L`, breaking
  the uniform signature, so a field is natural).

### 4. `MomentDisplacement` widened (the ONLY surprise ripple)
`orpheus/transport/displacements/moment_displacement.py`. The
`FluxRole._mint_displacement` (`flux ⊖ flux`) COPIES EVERY init dataclass field
into `_DISPLACEMENT_CLS`. Adding `spatial_moments` to `HarmonicMomentField` →
the mint passed `spatial_moments` to `MomentDisplacement.__init__` → `TypeError`
(13 strict-gate reds). FIX: `MomentDisplacement` gained the same
`spatial_moments: int = 1` field + the widened `_phase_space_shape` (mirrors
its flux sibling). `AngularDisplacement`/`ScalarDisplacement` needed NO change —
they inherit `AngularField`/`ScalarField._phase_space_shape` (which reads the
tail off the shared space), and `AngularFlux`/`ScalarFlux` carry no
`spatial_moments` field for the mint to copy.

### 5. Sphinx anchor (scope item 4)
`docs/theory/discrete_ordinates.rst`: NEW `.. _spatial-moment-space:` section
+ `:label: spatial-moment-space-size` eq + `:mod:` cross-refs + a `.. todo::`
for the archivist (the typed-space-half narrative the prior scattering-lift
stub's TODO had flagged as owed). Build exit 0; the `spatial-moment-space`
anchor rendered (3 hits in HTML); NO new warning references my label/files (the
only warnings are PRE-EXISTING `SyntaxWarning` in `test_projection_operators.py`
/`test_fission_operator.py` + the pre-existing `ld-cartesian-1d`/`ld-slab`
verifies-without-matching-eq skips).

## CONSTRUCT-GENERAL-ONLY — NO production field carries the axis yet
The `spatial_moments` parameter defaults to `1` EVERYWHERE and is NOT auto-read
from `mesh.scheme.spatial_basis_per_axis`. This is the deliberate
construct-general / select-narrow discipline (scope item 3): the CAPABILITY
exists, default-OFF, so DD/Step AND LD field shapes are unchanged this step.
Auto-reading the scheme would silently widen LD fields and break LD byte-id
before the consumers that FILL the axis exist (a Pattern-4 violation — an axis
no producer fills is an illegal state). The iterate/cell-emit/source seams that
thread the scheme's `spatial_basis_per_axis` here (SELECTING the axis for LD)
are S3-A proper (the NEXT sub-step). **Mutation-verified:** making
`_compose_spatial_moments` auto-read the scheme turns the `[ld]` byte-id gates
RED — the gate has teeth on exactly that mistake.

## GATES (all GREEN)
- **NEW foundation tests, 34P -O:**
  `tests/numerics/test_spatial_moment_space.py` (22 — space layer: shape/metadata/
  `from_per_axis`/`find_factor` round-trip + raises-when-absent + composition
  shape + `per_axis==1` no-widening + `average_moment_index`==`_ubld.AVERAGE_MOMENT`
  + `n_moments`==`per_axis**ndim` + size-identity eq) +
  `tests/sn/spatial/test_spatial_moment_field_space.py` (12 — byte-id-at-default
  for DD AND LD on all 3 carriers [the negative control] + widened d=1/d=2 shapes
  + both-moment-factors-coexist + from_mesh roundtrip + wrong-shape raises).
  **Mode-8/L26 SAFE:** ALL asserts are `np.testing.*` / `pytest.fail` / a
  `_check()` helper / `pytest.raises` — NO bare `assert` (confirmed bare assert
  is a NO-OP under `-O`). Mutation-verified: break the `per_axis**ndim` size law
  → 9 reds under `-O`; auto-read the scheme → 2 `[ld]` byte-id reds.
- **DD/Step bit-identity strict gate (THE negative control):**
  `tests/sn/sweep/core tests/sn/solve -W error::DriftWarning` = **513P/1skip/4xf**,
  IDENTICAL to the documented S2 baseline pre==post (re-confirmed 513, NOT 562).
- `tests/numerics` **644P** (incl my 22 + the `find_factor` addition).
- `tests/transport` **204P** (the `MomentDisplacement` + `_bases` ripple clean).
- `tests/sn/operators tests/sn/spatial` **573P / 7F** — the 7F are the EXACT
  documented PRE-EXISTING reds (sphere 1-D matvec SPH ×3 + `Face 'ymin' mu_y` ×2
  + sphere curvilinear apply ×2), `git stash`-confirmed identical at clean tree.
  ZERO new failures from my changes.
- Space + harmonic-moment-field floor **76P**.

## FILES CHANGED (mine, this session)
- `orpheus/numerics/spaces/spatial_moment_space.py` — NEW (the type + tail policy).
- `orpheus/numerics/spaces/__init__.py` — export `SpatialMomentSpace` + `spatial_moment_tail`.
- `orpheus/numerics/space.py` — NEW `TensorProductSpace.find_factor` (closes #207).
- `orpheus/transport/fields/_bases.py` — `_compose_spatial_moments` +
  `_spatial_moment_tail_of` helpers; `spatial_moments` param on Angular/Scalar
  `_space_for_mesh`/`from_mesh`/`zeros_on`; widened `_phase_space_shape`.
- `orpheus/transport/fields/harmonic_moment_field.py` — `spatial_moments` field +
  widened `_phase_space_shape` + `from_mesh_and_L`/`zeros_for_mesh_and_L` param.
- `orpheus/transport/displacements/moment_displacement.py` — `spatial_moments`
  field + widened `_phase_space_shape` (the mint-copy ripple fix).
- `docs/theory/discrete_ordinates.rst` — the `spatial-moment-space` stub.
- `tests/numerics/test_spatial_moment_space.py` — NEW (22 foundation).
- `tests/sn/spatial/test_spatial_moment_field_space.py` — NEW (12 foundation).
NOT MINE (pre-existing uncommitted, LEFT UNTOUCHED): `material_xs_field.py` (the
scattering lift), the other `docs/`/`.claude/` working-tree edits.

## OWED (S3-A proper, the NEXT dispatch)
The iterate carriers SELECTING `spatial_moments=mesh.scheme.spatial_basis_per_axis`
(both `AngularFlux` un-windowed + `HarmonicMomentField` windowed) + cell-emit φ̂
accumulation (`_CellSolve.cell` sweep_graph:883-895) + 2 source seams (d≥2
`_ubld_system` genuine `(2^d,ng)` Q; d=1 scan `D1ClosedForm.kernel_rhs`/`schur_xV`
slope) + d=1 matvec slope iterate + GATES 1/2/5. The space + factory CAPABILITY
this dispatch built makes those straightforward (the factory call site just
passes the scheme's per-axis count; the validator already accepts the widened
space). archivist DISPATCH emitted (followup:false).

## LESSON (candidate for the crosswalk / coding-elegance Pattern 4/7)
Minting a first-class typed FACTOR (not a bare int axis) for a quantity that
travels between solver steps is the Pattern-4 win the closeout's BLOCKER asked
for: a φ̂-carrying field becomes a LEGAL, typed, `find_factor`-queryable state
instead of an opaque widened ndarray. TWO ripples to expect when adding a
dataclass FIELD (not just a space factor) to a Field leaf: (a) the
`FluxRole._mint_displacement` field-copy → the sibling Displacement needs the
SAME field (else `flux ⊖ flux` raises `TypeError`); (b) the `_phase_space_shape`
cross-check validator must learn the new width — prefer reading it OFF the space
(`find_factor`) when the leaf has no field, vs the stored field when it does
(`HarmonicMomentField`/`MomentDisplacement` already carry `L`). The
circular-import trap is the third: a numerics-layer space that single-sources a
constant from the SN package must DEFER the import (the numerics layer must not
import SN at module-load time).
