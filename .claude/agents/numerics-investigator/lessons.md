# Numerics Investigator Lessons

Read this at the START of every invocation.

---

## L1: Never skip to "try fixes" — the cascade exists for a reason

Wasted 6 wrong hypotheses on cylindrical DD divergence by guessing at
fixes before isolating the broken component. Steps 3-5 (fixed-source,
component isolation, per-ordinate analysis) would have identified the
root cause directly. **Always run the cascade in order.**

## L2: Curvilinear redistribution is the first suspect

In cylindrical/spherical DD, the alpha recursion and geometry factors
(signs, delta-A/w scaling) are where bugs hide. A diverging-with-refinement
keff in curvilinear geometry means the balance equation itself is wrong
— go straight to per-ordinate flat-flux consistency (step 5) after
confirming spatial streaming works alone (step 4, alpha=0).
