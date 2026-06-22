---
name: w1-mm-tau-clamp-unclamp
description: W1 sphere M-M tau-clamp removal review — SSOT-at-producer PASS, doc-drift deferred-to-W5 is the legitimate anti-pattern-11 exception, named-intermediate inline is correct
metadata:
  type: project
---

W1 of the curvilinear-aniso-pole-clamp program (`.claude/plans/curvilinear_aniso_pole_clamp_program.md`),
branch `fix/curvilinear-aniso-pole-and-clamp`. Removed the `max(0.5,min(1.0,τ_raw))`
clamp from `spherical_streaming` M-M weight; retained it (with explanatory comment) in
`cylindrical_streaming`. Verdict: **PASS-WITH-NITS**.

**Why:** the clamp was an over-conservative positivity floor mis-cited to L&M §4.5;
Bailey-Morel-Chang 2010 Eq.43 IS the unclamped weight (unique exact-on-linear-in-μ).
100% spurious on physical fields. Correctness-of-design fix.

**How to apply (durable rulings for this file family):**

1. **SSOT is single-sourced AT THE PRODUCER.** `tau_mm` is built ONCE in
   `spherical_streaming` (`reduced_operator.py:654-657`) and inherited by EVERY consumer:
   `pole_angular_closure.py:635` (`self._tau_per_level = (reduced.tau_mm,)`),
   `sweep_cache.py:297`, `diamond.py:165/240`, `cell_balance.py:306`. NO consumer
   recomputes the clip. The SI-sweep path and the Krylov-matvec twin BOTH inherit the one
   producer value — Pattern-2/Pattern-7 satisfied by construction. When reviewing any τ /
   α-dome / ΔA-w change here, grep `\.tau_mm` to confirm no second computation site; the
   answer should always be "all reads, one write."

2. **Named-intermediate (Pattern 3): dropping `tau_raw` is CORRECT.** Once there is no
   clamped/raw distinction, the per-iteration expression has no second consumer and no
   distinct identity. `tau_mm` IS the named result (the M-M weight = fractional position of
   μ_n in its half-angle interval). Inlining the expression into `tau_mm[n]` removes a name
   that named nothing-distinct. Re-introducing `tau_raw =` would be ceremony.

3. **Sphere-unclamped / cylinder-clamped is a PRINCIPLED divergence, NOT a Pattern-2 twin.**
   The cylinder has a genuine STRUCTURAL singularity the sphere provably lacks: the
   most-inward azimuthal ordinate sits exactly on the level boundary (`eta[0]==eta_edge[0]
   ==-sinθ` bit-exactly at the producer → `τ_raw=0` → ÷0 in the recurrence). Sphere dome is
   non-singular (`dmu=w_n>0` on GL). Do NOT recommend unifying — premature
   (feedback_unify_after_two_instances). The clamp difference is a math-fact difference, not
   a copy-paste twin. Cylinder's real fix = 2-D (η,φ) closure, tracked #229.

4. **The doc-drift is the legitimate anti-pattern-11 exception (deferred, tracked).** The
   commit makes the sphere code contradict (a) its OWN module docstring
   `reduced_operator.py:58-69` (`:label: morel-montry-clamp` still shows clip), (b) the
   theory page `docs/theory/structured_geometry.rst:252-262` (same label, clipped), (c)
   `docs/verification/matrix.rst:741` lists `morel-montry-clamp` verified, (d) the
   consumer-side STALE COMMENT `cell_balance.py:310` ("clamped to (½,1] for curvilinear").
   The program EXPLICITLY assigns all doc surfaces to W5 (archivist), W1 step 2 scopes only
   the inline comment. So it is tracked, not orphaned. BUT: this is a transient
   same-file/same-equation-label contradiction; flag as a NIT to (i) tighten the
   `cell_balance.py:310` stale comment now (it is solver code, one line, in W1's blast
   radius), and (ii) ensure W5 retires the duplicated `morel-montry-clamp` label (it lives
   in BOTH the .py docstring AND the .rst — itself a pre-existing SSOT smell for an
   equation-label). Equation labels duplicated across .py-docstring + .rst is a latent
   drift habitat independent of W1.

5. **Linearity preserved = the deep architectural reason this static removal is right.** A
   ψ-dependent (dynamic) negative-flux fixup would make the operator nonlinear and break the
   linear-Krylov / SI≡Krylov twin identity. The clamp removal is config-time (static), so
   both twins stay linear and stay identical. This is WHY "drop the static clamp" beats "add
   a fixup" — record it as the load-bearing rationale.

Cross-ref: [[wave_o_operator_algebra]] (the SI≡Krylov twin discipline this preserves),
project_curvilinear_sn_cluster (the #229/#9 cluster).
