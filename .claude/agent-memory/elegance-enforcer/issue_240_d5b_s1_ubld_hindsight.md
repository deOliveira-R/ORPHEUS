---
name: issue-240-d5b-s1-ubld-hindsight
description: HINDSIGHT placement+twin audit of the landed D5b-S1 UBLD deliverables (Branch-1 sympy + Branch-2 numpy + LD consumers) — PASS, two trivial nits, one escalation
metadata:
  type: project
---

HINDSIGHT architecture audit (not per-diff) of the LANDED D5b-S1 UBLD work on
`feature/sn-space-angle-tier2` (commits `cb84b7b` Branch-1, `69b19c9` Branch-2).
Three artifacts: `orpheus/derivations/discrete/sn/ld_ubld.py` (sympy algebra-of-record),
`orpheus/sn/spatial/_ubld.py` (numpy primitive + d=1 fast path), `linear_discontinuous.py`
(3 prod consumers of `d1_closed_form` @ 340/460/581). Verified: 5 oracles PASS,
16 conformance tests pass. Verdict **PASS**, two trivial follow-ups + one escalation.

**Why:** the main agent + author asked for a critical re-examination against confirmation
bias — did it land in the RIGHT PLACE, are there twins.

**How to apply (the durable rulings — reuse on the S2 carve):**

- **A1 placement of `d1_closed_form` in `_ubld.py` (KEEP, ESCALATED over the brief's framing).**
  The deciding axis is NOT "one-class consumer → maybe move into the class." It is
  **co-location with the reduction it is proven == to**. `d1_closed_form` is the analytic
  Schur of `per_cell_solve`'s d=1 2×2; its correctness contract is "== the dense primitive
  at d=1." Moving it onto `LinearDiscontinuous` splits a proven-equal pair across a module
  boundary → the day S2 generalizes the dense primitive, the d=1 fast path orphans and drifts.
  It is UBLD-algebra-LD-consumes (L13 scheme-consumes-primitive), not LD-private-behavior.
  Generalize: a fast path lives WITH its reference oracle, not with its single consumer.

- **A2 dense primitive (`assemble_ubld`/`per_cell_solve`/`assemble_inflow_axis`) in S1 with
  ZERO prod consumer is NOT Pattern-6 premature (KEEP, more permissive than the brief feared).**
  It is the reference-oracle retirement-exception (my memory note #3): a fuller/more-expensive
  view kept as the verification oracle is correct iff kernel shared + foundation-tier pin +
  corner probed. ALL THREE hold: shared Kronecker kernel; `@foundation` 1e-12 pins; and
  `test_d2_exact_on_bilinear` exercises the xy cross-coupling the d=1 paths are BLIND to and
  CATCHES ERR-060 (missing |μ_axis| on multi-axis inflow). The dense primitive is what makes
  the d=1 fast path's "single-source the math" a NUMERICALLY-verified claim, not symbolic-only.
  Shipping it in S2 instead would be WORSE (d=1 fast path pinned only to the symbolic oracle
  across the structural-independence line). "Test-only production-shaped numpy" is the right
  S1 state when it is a live oracle for code that ships today.

- **B1 sympy↔numpy twin = JUSTIFIED algebra-of-record (KEEP). Do NOT single-source the four
  factors across branches** — that defeats structural independence (two witnesses → one).
  Drift IS pinned both halves: Branch-1 self-pins via `simplify(diff)==0`; Branch-2 pinned to
  Branch-1 by `test_d1_assembled_matrices_match_symbolic`; d=2 numpy pinned independently by
  exact-on-bilinear physics. `balance.py` confirmed clean (DD/WDD scalar-balance, no Kronecker
  / no LD mass / no UBLD factors; only overlap is the generic `a·ψ_in+b` shape + `s=2μ/Δx`,
  which is DD's own recurrence). **The one CONCERN: the numpy-vs-symbolic MATRIX-equality pin
  fires only at d=1.** d≥2 numpy assembly is pinned only by physics (exact-on-bilinear), never
  entry-for-entry vs symbolic. Acceptable now; S2-plan bullet = add d2-A==d2-A when d≥2 ships
  (a symmetric xy-block error could survive bilinear-exactness in one branch).

- **B2 the ONE real duplicated expression (COLLAPSE, do-now trivial):** `schur_xV` and `scan_xV`
  both recompute `self.g * V` (`_ubld.py:346` `mu_Adown` vs `:371` `m`) AND `V * self.eff_denom`
  (`:348` vs `:372` verbatim). These ARE the two ×V-convention named intermediates (m=|μ|A_down,
  S=V·eff_denom). Extract `D1ClosedForm._xV(V) -> (m, S)`; both consume it. Latent today (byte-id);
  live the moment ÷V→×V convention changes (curvilinear `V_eff`, or the D6 `w(Σ)` Péclet blend
  touching the Schur diagonal). Same class as standing tell #2 but the collapse is trivial
  (one helper, same class) → no reason to defer.

- **B3 `derive_*` oracles hand-transcribe prod closed forms = NECESSARY L11 cross-check (KEEP)**
  — importing prod would make the proof vacuous (`x==x`). NIT: they cite prod by LINE NUMBER
  (`ld_ubld.py` comments @ 345/395/437) → use SYMBOL cites (recurring team tell, also in the
  Branch-1 note).

- **B4 no other twin.** DD's `2.0*s`/`_DD_W` = DD's own diamond w=½ (separate scheme, same
  `outgoing_face_from_average` primitive at different w). The shared blend/reconstruction lives
  ONCE on `DiscretizationSchemeBase` (predates carve). `D₂'=θ·V·d2` extra-θ crosswalk
  single-sourced as named `d2p`, oracle-pinned — a WIN (the historical extra-θ bug-habitat is
  now named+single-sited).

- **`assemble_inflow_axis` `axis∈{0,d-1}` else-raise = honest-per-phase capability boundary,
  NOT anti-#7.** Latent: S2 d≥3 interior-axis interleave is un-derived; the two boundary
  reshapes (`:242-248`) must collapse to ONE general insert then (do not pre-build — Pattern 6
  applies THERE).

No new issue (the 2 S2 carry-forwards belong in the existing #240/#38 S2 plan). Links:
[[issue_240_d5b_s1_ld_ubld_branch1]] [[issue_240_d5b_s1_ld_ubld_branch2]]
