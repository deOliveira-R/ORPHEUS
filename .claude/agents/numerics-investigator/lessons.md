# Numerics Investigator Lessons

Read this at the START of every invocation.

---

## 2026-04-05 — Cylindrical DD divergence

- **Symptom**: heterogeneous keff diverges with refinement (1.15→0.90→0.52)
- **6 wrong hypotheses tested**: reverse sweep, step closure, starting
  direction, bidirectional sweep, scaled α, zero redistribution
- **Root cause**: wrong α recursion (cumsum(+w·ξ) → cumsum(−w·η)) +
  missing ΔA/w geometry factor. Both broke per-ordinate flat-flux
  consistency.
- **What worked**: fixed-source diagnostic (step 3) revealed flux spike
  at r=0. Component isolation (step 4) with α=0 proved spatial streaming
  was correct. Per-ordinate analysis (step 5) proved the balance equation
  itself was wrong.
- **Key lesson**: the cascade would have found this in steps 3-5 if
  followed systematically. The 6 wrong hypotheses came from skipping
  straight to "try fixes" without isolating the component first.
