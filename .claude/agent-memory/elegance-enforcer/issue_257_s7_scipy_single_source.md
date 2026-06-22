---
name: issue-257-s7-scipy-single-source
description: #257 S7 PASS — single ORPHEUS↔scipy Krylov boundary (_as_scipy_linop); retired the zero-caller as_scipy_linop twin
metadata:
  type: project
---

#257 **S7** — scipy single-source plumbing carve. PASS (`feature/field-typed-operator-algebra`,
HEAD `9ccfec7`, UNCOMMITTED, behavioral-NEUTRAL claim). The cleanest carve in the S-series:
a textbook Pattern-2 consolidation with a textbook Pattern-3 naming win folded in.

**Why:** Cardinal Rule 2 — the live ORPHEUS↔scipy Krylov boundary was THREE `spla.LinearOperator`
construction sites across `numerics/`: 2 inline closures inside `KrylovAcceleration.solve`
(`A_matvec`/`M_matvec`, each duplicating `_ravel`/`_unravel_like`) + a PUBLIC `as_scipy_linop` in
`operator.py` that was a FLAT-ONLY TWIN with ZERO production callers (only 5 tests + 3 doc refs).

**How to apply / what was verified:**
- ⭐ The factoring is GENUINE, not procedural. `_as_scipy_linop(carrier_matvec, template, n)` extracts
  the ravel-WRAP (the true shared primitive) and exposes the carrier matvec + template as the HONEST
  variation points: A=`loss_minus_gains`/`solution_template` (flux space), M=`self._preconditioner`/
  `q_ext` (source space). The A↔M template ASYMMETRY (flux vs source space) is load-bearing and
  correctly preserved — B.5.2 (iterate lives in solution space, b/precond-input in source space).
  The two closures' ONLY real difference (carrier matvec + template) became two args; everything they
  SHARED (the flat lift/ravel) became the one primitive. Verified: exactly ONE `spla.LinearOperator(`
  construction in all of `orpheus/` now (was 3) — `iteration.py:228`. Both callers route through it.
- ⭐ Pattern-3 win bundled in: the old anonymous `A_matvec` FLAT closure (`psi_flat`→ravel→compose→
  ravel) became the named CARRIER-space `loss_minus_gains(psi: V) -> V` reading as the within-group
  system matvec `(L+C − S − B)·ψ` — the math is now visible because the ravel plumbing moved OUT to the
  adapter. Master-standard improvement (reads like the domain).
- Pattern-5/6: dropping the retired twin's `rmatvec`/`dtype`-param was RIGHT (defer abstraction). Both
  live callers are float + matvec-only; the rmatvec support was TEST-ONLY (the 3 deleted rmatvec/
  transpose tests exercised a capability no production caller used). `dtype=float` hardcode honest.
- Retirement CLEAN: def + module-docstring ¶ + orphaned `spla` import + `__all__` (both operator.py and
  numerics/__init__.py) + 5 orphaned tests + 3 doc xrefs all removed. NO dangling. Capability symbols
  `_has`/`CAP_APPLY`/`CAP_APPLY_TRANSPOSE`/`CAP_SOLVE`/`MissingCapability` confirmed STILL heavily used
  in operator.py (≥40 sites) — NOT over-eagerly removed. Test fixtures `NoApplyOperator`/
  `matrix_apply_only`/`MissingCapability` retain other live consumers — surgical deletion.
- SI path (`SourceIteration.solve`, 329-575) confirmed legitimately scipy-FREE — every scipy ref
  (122/210/228/777+/802/821) is inside `KrylovAcceleration` (576+) or the adapter. NOT a missed target
  (SI is a fixed-point loop, no linear-operator handoff).
- Gates: test_operator.py 55✓ (post-deletion), Krylov-path 18✓. pyright: the 3 `reportCallIssue` at the
  new `:228` (`matvec`/`dtype`/`shape` overload-resolution noise) are a PRE-EXISTING scipy
  `LinearOperator.__init__` STUB artifact — IDENTICAL construction `spla.LinearOperator((n,n),
  matvec=...,dtype=float)` existed at HEAD on TWO sites; S7 reduced the stub-noise footprint 2→1.
  The operator.py `block_role` ClassVar errors are the known dataclass-field-trap (memory #5), outside
  the deleted region.

NITS (both OPTIONAL, non-blocking): (N1) `_as_scipy_linop` docstring says "system matvec `(L − Σ gᵢ)`"
while the `loss_minus_gains` comment more honestly says `(L+C − S − B)·ψ` — minor internal inconsistency
on the SAME operator (C folds into L's `.apply` for SN within-group; both spellings are correct but the
adapter docstring drops the `+C`). (N2) the `except TypeError` fallback comment at `:814` still names
"`A_matvec` / `M_matvec`" — the closures it references were renamed to `loss_minus_gains` /
`self._preconditioner` by this very diff → stale identifier in a load-bearing comment (the
docstring-asserts-a-name-the-code-no-longer-has tell, mild form).

STANDING TELLS reinforced: (1) a public utility with ZERO production callers + test-only extra
capability (rmatvec) is a genuine twin to RETIRE, not keep "for future migration" — the Issue-15
migration rationale in the old docstring was a Pattern-11 untriggered-temporary; (2) extracting the
SHARED plumbing (ravel-wrap) while NAMING the carrier matvec simultaneously satisfies P2+P3 — the
plumbing leaving the closure is WHAT makes the math visible; (3) when a carve renames inline closures,
grep the surrounding error-handling comments for the OLD closure names (stale-comment habitat).
