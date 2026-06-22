---
name: issue-245-moment-layout-relocation
description: #245 numerics/moment_layout.py relocation — PASS-WITH-NITS; the recurring "single-source claim names the re-export not the post-move home" staleness from AXIS_NAMES-style down-relocations
metadata:
  type: project
---

#245 (branch `refactor/sn-foundation-cleanup`, uncommitted at review) relocated the
physics-free moment-layout policy — `AVERAGE_MOMENT` (Kronecker slot-0 average index) +
`face_moment_tail(n)` ("append a trailing moment axis iff n>1") — from `orpheus/sn/spatial/_ubld.py`
UP-reach to a NEW leaf `orpheus/numerics/moment_layout.py`, killing the numerics→`orpheus.sn`
deferred-import band-aid. **Verdict: PASS-WITH-NITS.**

**Why:** band-aid / layering-inversion cleanup, byte-identical, foundation-cleanup pass.
**How to apply:** the durable lesson below recurs on EVERY AXIS_NAMES-style down-relocation.

## What was RIGHT (don't re-litigate)
- HOME = top-level `numerics/moment_layout.py` next to sibling `numerics/face_layout.py` (both
  physics-free layout descriptors). NOT `numerics/spaces/` (holds FunctionSpace classes, not layout
  primitives); NOT co-homed in `spatial_moment_space.py` (a primitive consumed by TWO siblings —
  the SpatialMomentSpace CAPABILITY half + the `_ubld` REALIZATION half — must sit BELOW both, not
  inside one, else mirror-image inversion). Taxonomy = layer × kind.
- RE-EXPORT (`_ubld` keeps a `# noqa: F401` down-import + `__all__` retention) is IDIOMATIC, NOT a
  Pattern-2 twin. Discriminator = OBJECT IDENTITY: `_ubld.X is moment_layout.X` (verified True for
  both symbols) → one definition, one alias, zero divergence surface. Precedent = `orpheus/sn/axis.py:72`
  re-exporting `AXIS_NAMES` down from `numerics.face_layout` (identical noqa+__all__ idiom). Keeping
  ~5 SN importers (`loss_representation`/`linear_discontinuous`/`sweep_graph`/`solver`) on `._ubld` is
  correct — symbols belong next to the UBLD Kronecker primitives they name at SN-side use; mass-rename
  = blast radius for zero structural gain.
- Cycle GENUINELY broken: importing `spatial_moment_space` pulls ZERO `orpheus.sn` modules (smoke-verified).
  Deferred import + `_average_moment_index()` wrapper FULLY RETIRED (grep-clean), not paralleled.
- New module docstring = exemplary Cardinal-Rule-3 (WHY/two-consumers/predecessor/AXIS_NAMES precedent/
  leaf-guarantee-no-recur).

## ⭐ STANDING LESSON (the recurring smell) — single-source CLAIM names the RE-EXPORT, not the post-move HOME
A down-relocation modeled on AXIS_NAMES leaves a STALE-CLAIM blast radius in the SN consumers that import
THROUGH the re-export: their docstrings still assert "**single-sourced from** `_ubld.X`" / "**single source
via** `_ubld.X`". After the move `_ubld.X` is a RE-EXPORT, not the source — the claim is strictly false.
- Discriminator: an IMPORT through `_ubld` is FINE (intended ergonomic). A SINGLE-SOURCE CLAIM ("the policy
  lives HERE") is an architectural assertion that must name the REAL home (`moment_layout`). A pure
  cross-ref to `_ubld` for something that GENUINELY lives there (`assemble_ubld`, the Kronecker primitive —
  `linear_discontinuous.py:310/461`) is correct, LEAVE IT.
- Bug habitat: the false pointer sends a future policy-editor to `_ubld` where there is no definition to
  edit → tempts a re-inlined local copy in `_ubld` → re-creates the exact twin the move eliminated.
- #245 instances flagged (Pattern-7 / anti-#11): `loss_representation.py:350` ("Single source via
  `_ubld.face_moment_tail`") + `transport/fields/_bases.py:180` ("single-sourced from `_ubld.face_moment_tail`
  via `spatial_moment_tail`"). The two CLAIMS in `spatial_moment_space.py` itself WERE correctly repointed
  to `moment_layout`; the SN-consumer claims were missed (they import via the re-export → invisible to the
  diff's own files).
- Cosmetic sibling: the foundation test was repointed in IMPORT + docstring BODY but its NAME still encodes
  the old source (`test_average_moment_index_matches_ubld_single_source` → should be `...matches_moment_layout
  _single_source`). A test name is a contract label; grep-for-`_ubld`-coverage and grep-for-`moment_layout`-
  coverage both mislead.
- GENERAL RULE for future AXIS_NAMES-style moves: after relocating, grep the WHOLE tree for "single-source"/
  "single source"/"single-sourced" + the moved symbol name; every such CLAIM must name the post-move home,
  even in files the move didn't otherwise touch (they reach it through the re-export, so the diff is blind
  to them). The false claims hide where the diff doesn't reach — same shape as the Phase2-StepC deleted-symbol-
  cite lesson ([[issue-236-phase2-stepc-tau-retirement]]), here for a RE-EXPORT instead of a deletion.
