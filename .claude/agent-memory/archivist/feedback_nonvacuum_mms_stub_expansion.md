---
name: Non-vacuum prescribed-inflow MMS stub expansion (Phase 4 / O.2b 4.6)
description: Pattern for expanding a method-implementer's 4-:label: MMS stub into rich narrative when the novelty is a BC-trace property (non-vacuum inflow), not a new operator; the cross-domain-attacker frame memo carries the load-bearing rationale (Q1 linear-μ-fully-activates, Q2 native-P1, H1 pole regularity)
metadata:
  type: feedback
---

# Phase 4 / O.2b 4.6 non-vacuum MMS expansion (2026-06-06)

Expanded the `_sn-mms-nonvacuum` stub in `docs/theory/discrete_ordinates.rst`
(4 `:label:` blocks + `.. todo::`) into ~700 LoC rich narrative. Branch
`refactor/field-role-typing`. Sphinx `-W` clean (0→0 warnings, gate satisfied).

**Why this one was distinctive:** the verification novelty was NOT a new
operator or a new discretisation — it was a single BOUNDARY-TRACE PROPERTY
(`a0>0` makes the inflow trace `γ₋ψ≠0`). The angular form, interior
operator, redistribution term, and quadrature were ALL shared with the
already-documented Phase 3.6 anisotropic cases. So the narrative's spine is
"what is the single structural delta, and why does it light up an untested
path" — not "here is a new derivation". The expansion frames everything
around that delta.

## The source-reading order that worked

1. **Method-implementer closeout memo** (`phase4_o2b_46_nonvacuum_mms_closeout.md`)
   — the measured-numbers table + decisions A-E + the T3 stagnation values.
2. **cross-domain-attacker frame memo** (`phase4_o2b_4_6_mms_ansatz_frame.md`)
   — THE load-bearing rationale source for the hardest prose: Q1 (linear-in-μ
   FULLY activates redistribution because the closure is LINEAR → quadratic-μ
   adds NO new structural coverage; settles "do we need μ²?" = NO), Q2 (the
   ansatz IS the native P1 = P0⊕P1 Legendre element), H1 (the pole-regularity
   `B(0)=0` hard constraint, with the (r/R) prefactor mechanism). This memo
   is where the "WHY" lives that the SymPy docstrings only gesture at.
3. **SymPy source** (`mms/sn.py`) — the canonical algebra-of-record; the
   `derive_*` docstrings give the substitution chain verbatim.
4. **All three test files** — the measured evidence (orders, L2, reldiffs),
   the bypass-solve mechanism (`_within_group_triple` + `_select_si_resolvent`
   + flux `initial_guess`), the three-assertion value-vs-rate structure, the
   Mode-7 activates/nulls declaration (already written IN the test docstrings).
5. **Existing sphere-aniso section** (`sn-mms-curvilinear-aniso-verification`)
   for tone/depth template + the `sn-mms-spherical-aniso-qext` label to
   cross-ref (the "same closed form" Cardinal-Rule-2 claim).

## Narrative shape (12 sections under the `----`)

Key-Facts admonition (`:class: important`, 7 bullets) → 10-row gate
list-table at the head → vacuum-automatic motivation (with the γ₋ψ=0
derivation) → ansatz/P1/why-linear-μ-enough → slab source derivation
(full \begin{aligned} substitution) → spherical source / Cardinal-Rule-2
reuse → H1 pole regularity (`^^^^` sub-subsection) → non-vacuum lever
(a0>0) → affine-BC-to-RHS framing → bypass path + B.5.2 flux/source bridge
→ converged-value-not-rate (vv anti-pattern #5) → Mode-7 term map list-table
→ T3 xfail rationale + stagnation table + T3g green companion → T4 vv Mode 9
→ verification chain Branch1/Branch2/L1 + L11 + "does NOT verify eigenvalue".

## Cross-ref targets that resolved cleanly (verified BEFORE writing)

- `:ref:affine-bc-form` AND `:eq:affine-bc-form` — boundary_conditions.rst
  line 158/166. THE perfect target for the affine-BC-to-RHS framing (the q
  inhomogeneous term). Both the section anchor and the eq label exist.
- `:class:~orpheus.geometry.boundary.PrescribedInflow` — exists at
  `geometry/boundary/prescribed_inflow.py`. Used to say "this BC bridge
  exists but 4.6 does NOT touch it".
- `:eq:sn-mms-spherical-aniso-qext` / `-psi` — the Phase 3.6 sibling, for
  the "same closed form, only A,B differ" Cardinal-Rule-2 claim.
- `:ref:sn-pole-angular-closure-protocol` — the Carlson/Morel-Montry closure
  the redistribution rides on (line 2277).
- `:ref:sn-mms-curvilinear-isotropic-verification` /
  `-aniso-verification` — the vacuum-automatic predecessors.
- All `:func:`/`:class:` to `mms/sn.py` resolve by DEFINITION site (the
  module IS automodule'd somewhere in api/). Test `:func:` refs use the
  dotted `tests.derivations....` / `tests.sn.verification.analytical....`
  paths and resolve.

## Gotchas / disciplines confirmed this task

- **vv-status equation labels: -psi → documented, -qext → verified.** The
  ansatz `:label:`s (`-psi`/`-sph-psi`) are DEFINITIONAL (the angular flux is
  imposed, not solved) → annotate `documented`; they correctly land in the
  matrix.rst "Documented-only equations" list AND in the verified table (the
  tests' `verifies(...)` cite them too — dual presence is correct, NOT an
  orphan). The source `:label:`s (`-qext`/`-sph-qext`) are DERIVED closed
  forms verified by SymPy `simplify==0` + bit-equal cross-check → `verified`.
- **vv-status comment uses the `.. (vv-status rationale) <category>:` prefix
  THEN `.. vv-status: <label> <status>` on the next comment line** (matching
  the existing line 5372-5373 pattern in this file). Categories used:
  `definition` (the ansatz) and `derivation` (the source).

(Build-gate, title-marker-ladder, and memory-dir discipline are in
AGENT.md "Build-Gating & Cross-Ref Reality".)

## Quality self-assessment (Directive 3)

| Dimension | Score | Notes |
|-----------|-------|-------|
| Derivation depth | 5 | full slab substitution \begin{aligned}, γ₋ψ=0 vacuum-automatic derivation, φ=A under Σwμ=0, the linear-μ-fully-activates argument written out |
| Cross-references | 5 | every `:func:`/`:class:` to mms/sn.py + all 9 test funcs linked; `:eq:`/`:ref:` to affine-bc-form, the 3.6 sibling, the pole-closure protocol |
| Numerical evidence | 5 | 10-row gate table + T3 stagnation table (2.37/2.42/2.43e-2) + T1/T2 orders + T4 reldiffs 1.3e-13…5.6e-13 |
| Failed approaches | 5 | the dropped-q.boundary failure mode (converges O(h²) to WRONG limit), the slab-style B(0)≠0 wrong-on-sphere, the slab-only Mode-7 trap, the μ²-not-needed settling |
| Code traceability | 5 | every equation → its `derive_*` SymPy func + its Branch-2 factory + its test gate |
| Derivation source | 5 | every algebraic claim cites the originating `derive_nonvacuum_*` / `_spherical_anisotropic_symbolic` SymPy function; the frame memo Q1/Q2/H1 carries the design rationale |

5/5 across the board: this was an ideal case — closeout memo + frame memo +
SymPy + rich test docstrings gave near-complete fact density. The frame
memo's Q1 (linear-μ-fully-activates) was the single most valuable input;
without it the "do we need μ²" answer would have been hand-waved.
