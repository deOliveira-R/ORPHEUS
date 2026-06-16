---
name: issue-240-phase2-step-d5-0-scheme-trait
description: "#240 Ph2 D5-0 PASS-WITH-NITS — transverse_coupling_is_facewise trait kills silent 2-D LD→ScanMarch(inline-DD) misroute; the runtime_checkable Protocol data-member footgun ruling"
metadata:
  type: project
---

# #240 Phase 2 Step D5-0 — scheme-capability trait (routing-honesty fix)

Branch `feature/sn-space-angle-tier2`, reviewed PRE-commit. VERDICT **PASS-WITH-NITS** (1 blocking).

**The fix.** A 2-D Cartesian LinearDiscontinuous mesh SILENTLY routed to `ScanMarch`,
whose row-march interior runs INLINE DiamondDifference with no scheme dispatch → computed
DD, dropped LD's bilinear slope. Root cause = a DIMENSIONAL predicate (`is_cartesian and
ndim==2`) standing in for a SCHEME-CAPABILITY question. Fix mints `transverse_coupling_is_facewise:
ClassVar[bool]=False` (Base `scheme.py:721`), `True` on DD (`diamond.py:143`), LD inherits
default. `ScanMarch.supports` (`loss_representation.py:1227`) splits: 1-D arm reads
`is_affine_scannable` (unchanged — 1-D LD scans via CumprodScan), d≥2 arm reads the new trait.

**Why the architecture is RIGHT (review the execution, not the choice).**
- Pattern 4 textbook: illegal LD/scan-march pairing now unrepresentable. `False` default =
  correct conservative (mirrors sibling `is_affine_scannable: ClassVar=False`; both capability
  traits default exclusionary → forgotten opt-in is "slow but correct," never "fast but wrong").
- Honesty claim VERIFIED e2e: wavefront LD `cell_kernel_batch` RAISES real `NotImplementedError`
  (`linear_discontinuous.py:437`) on d≥2, so misroute → LOUD refusal not different-silent-wrong.
- Pattern 2 SSOT: value in 2 places (Base default + DD override = minimal), predicate read at
  ONE prod site (`:1246`). No second spelling in registry/default_for. Dispatcher+guard share
  the one `supports` predicate → cannot drift.
- Trait NAMED for scheme-property (separable transverse coupling) NOT strategy (`is_scan_march_compatible`
  rejected) → diffusion ADI/line-SOR preconditioner (#240's confirmed next consumer) reuses it,
  no rename. The forward-reuse rationale is the cross-domain-attacker Frame 1.
- Docstrings EXEMPLARY (Base carries numerical-PDE statement + is_affine_scannable contrast +
  per-scheme table + forward-reuse + Maginot-Ragusa-Morel 2016/Adams 2001/Lewis-Miller lit).
- Tests EXEMPLARY: negative+positive both-ways pin (anti-pattern-#11) + fake/real-mesh layering
  + strategy-free `TestSchemeTraitProbe` (reads trait off class, no ScanMarch/supports/mesh in
  scope → PROVES genuine scheme property, not frame leak). Non-square nx=4,ny=3 = x↔y axis-swap
  corner defense. 37/37 green.

**BLOCKING NIT — the runtime_checkable Protocol data-member RULING (user flagged).**
RULE AGAINST keeping the trait on `@runtime_checkable DiscretizationScheme` Protocol
(`scheme.py:402`); keep it Base-only (`DiscretizationSchemeBase`). Reason = `@runtime_checkable`
checks data-member PRESENCE not type/value → broke 2 fixtures (visible symptom) AND latently a
scheme declaring `transverse_coupling_is_facewise="yes"` passes `isinstance(x, DiscretizationScheme)`
then reads as truthy bool in supports → re-opens a NARROWER silent-misroute (the exact class the
trait was minted to close, Pattern-4 self-defeat). NO production `isinstance(x, DiscretizationScheme)`
exists (only the conformance test) → Protocol data members buy nothing in prod, cost the
presence-only false-gate. Fix = drop from Protocol body+docstring (keep 1-line xref to Base
ClassVar), revert the 2 fixture declarations (they conformed only because Protocol demanded
presence). ALT = keep 4 traits symmetric on Protocol BUT strengthen `TestProtocolConformance` to
assert each is a genuine `bool`. NOTE the 3 pre-existing traits (is_linear/is_positivity_preserving/
is_affine_scannable) ALSO sit on the Protocol with the same presence-only weakness — optional
follow-up (`module:sn`), do NOT bundle into D5-0.

**2 non-blocking nits.** (a) `.. todo::` in loss_representations.rst flagging the STALE pre-D5-0
supports code block = anti-pattern-#11 EXCEPTION (carries removal trigger #240 D6 + accurate
brief body) — confirm D6 is a tracked GitHub issue not just the closeout memo. (b) LD docstring
"slope-wise (1st-order)" prose adjacent to the 1-D-Schur-elimination note → one disambiguating
clause (eliminated per-axis in 1-D, cross-axis slope survives d≥2). Cosmetic.

**Lesson promoted to AGENT.md candidate: runtime_checkable + data members = presence-only gate.**
A capability ClassVar belongs on the concrete Base (real SSOT, typed-and-defaulted), NOT on a
runtime_checkable Protocol (presence-only, weaker contract masquerading as validation). When a
Pattern-4 trait exists to make a misroute unrepresentable, a presence-only Protocol gate re-opens
a narrower version of that same class. Check for production `isinstance(x, Protocol)` consumers —
if none, the Protocol data members are pure cost.
