# ERR-006 — Convergence is not correctness (cylindrical/spherical SN)

**The bug.** The curvilinear discrete-ordinates sweep had two
simultaneous defects — an α-recursion that cumulated `+w·ξ`
(azimuthal cosine) instead of `−w·η` (radial cosine), and a missing
`ΔA_i / w_m` geometry factor on the redistribution term — that
combined to produce a Richardson-stable solver which converged
smoothly to the wrong physics.

**Evidence that existed before it was caught.** Twenty tests passed:
the homogeneous eigenvalue was exact to machine precision; particle
balance closed by construction; conservation was exact; flux
non-negativity held; single sweeps ran to completion with finite
values; convergence under mesh refinement was monotone at the
expected order on the homogeneous problem. Every test the team had
written said the solver was correct.

**Why that evidence didn't catch it.** All twenty tests probed
regimes where the redistribution term either cancelled (flat flux
on a homogeneous medium) or was self-consistent with the same buggy
operator (Richardson extrapolation of the system under test). The
solver was internally consistent — it converged to *its own*
discretised equation — but the discretised equation was wrong. The
heterogeneous-refinement signal (k_eff: 1.15 → 0.90 → 0.52 → 0.25)
only emerged once mesh refinement was combined with a heterogeneous
problem and at least two energy groups, because only that
combination forces redistribution to do net work that the bug then
mis-routes.

**What evidence class would have caught it.** Convergence at the
expected order is necessary, **NEVER** sufficient — **instead**,
require either (a) a structurally independent analytical reference
(one whose derivation does NOT pass through the discretised
operator under test), or (b) a per-ordinate flat-flux consistency
check that audits the balance equation term-by-term against a
known closed form. The latter (L0-SN-003 in the catalog) would
have fired immediately: with the missing `ΔA/w` factor, the
per-ordinate streaming + redistribution residual is non-zero on
flat flux even though the *summed* balance still telescopes to
zero. Per-ordinate consistency is the fundamental L0 criterion for
curvilinear SN, and self-referential references (Richardson on the
system under test) cannot detect bugs that are themselves
convergent. See vv-principles reference.md §1 (structural
independence) and §4.4 (semi-analytical 2-step ladder, where the
"reduction correctness" step is exactly what fails here).

**References.** ERR-006 in `error_catalog.md`;
`docs/theory/discrete_ordinates.rst` §Investigation History (six
failed hypotheses before root cause); numerics-investigator agent
memory L1, L2 (diagnostic cascade order; curvilinear redistribution
as prime suspect); `@pytest.mark.catches("ERR-006")` regression
test in `tests/sn/test_cylindrical.py`.
