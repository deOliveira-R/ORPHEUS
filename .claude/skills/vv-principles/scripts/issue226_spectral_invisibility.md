# Issue #226 (step 5b) — the gate's functional is blind to its stabiliser

**The bug.** No shipped solver bug — a verification-PLAN bug, caught
pre-implementation: the step-5b teeth table claimed the k∞ value gates
would red a factor-swap mutation (spelling `F·A⁻¹` instead of `A⁻¹F`
in the homogeneous multiplication operator `K`) and a
resolvent-transpose mutation (`Aᵀ` for `A`). Had the carve shipped
with those gates credited as the catchers, an entire class of
composition-order and transpose bugs in the K assembly would have been
permanently invisible — green forever, at any tolerance.

**Evidence that existed before it was caught.** The full k-level gate
stack, all green and all tight: the cross-engine consistency gate
(`dominant_eigenpair` vs `direct_eigenvalue` at rtol=1e-12), the
structurally-independent SymPy closed-form k∞ anchor
(`test_kinf_exact`, 1e-12), and a planned mutation row asserting these
gates red the swap. By every conventional standard — independent
reference, tight tolerance, planned teeth — the k gates looked like
committed catchers.

**Why that evidence didn't catch it.** The measured quantity
k∞ = λ_max(K) is a SPECTRAL functional, and the mutation class sits
inside its invariance group: `A·(A⁻¹F)·A⁻¹ = F·A⁻¹` (the swapped
product is *similar* to the true one) and `eig(Mᵀ) = eig(M)` (the
transpose shares the spectrum). Measured: |Δk| = 0.0 **exactly** for
both mutations, while the matrix itself moves O(1)
(‖ΔK‖ ≈ 1.46 swap / 1.43 transpose). This is not Mode-10 sub-floor
magnitude — no tolerance, refinement, or regime change exposes it,
because the functional annihilates the error *algebraically*.
Tightness and structural independence of the reference are irrelevant
when the functional itself cannot distinguish the mutated object from
the true one. Note the epistemic twist: the plan HAD the right
mutation with the WRONG expected outcome — the empirical
run-the-mutation discipline alone would have reported "gate green"
and needed this diagnosis anyway; only the design-time algebra says
*why*, and says it before any code exists.

**What evidence class would have caught it.** A spectrum-level gate is
necessary, NEVER sufficient for mutations inside the spectrum's
stabiliser — **instead, pin the OBJECT**: a matrix-level intrinsic
gate `K.as_matrix() ≡ np.linalg.solve(A, F)` against an
independently-posed reference (committed as the step-5b catcher; both
mutations move it O(1)). Operationally, at gate-design time enumerate
the measured functional's invariance group (similarity + transpose for
spectra; per-term cancellation for balance/telescoping sums; global
scaling for normalised shapes) and intersect it with the threat model
— any mutation class inside the stabiliser needs an object-level or
out-of-stabiliser gate, then a red mutation (SKILL.md test-design
Mode 12; anti-patterns #3 and #8 are prior instances of the same
lens). Imminent application: the #276 A4 daggered eigenvalue —
`eig(Kᵀ) = eig(K)` by construction, so "k* matches k" carries zero
mutation coverage on the adjoint machinery; verify with
eigenvector/bilinear functionals (φ*-weighted reaction rates,
biorthogonality) instead.

**References.** Full similarity derivation + the T1–T5 mutation
measurements: `docs/theory/foundations/infinite_medium.rst`, the
`spectral-invisibility` section. The committed catcher:
`tests/homogeneous/test_homogeneous.py::test_K_operator_as_matrix_is_the_resolvent`.
Origin: #226 taxonomy step 5b, commit `394d8c0`
(`refactor/inverse-as-operator`); the corrected teeth row is recorded
in the step-5b verification plan §8 stamp and the commit message. No
ERR-NNN entry: the bug never shipped — this is a test-design failure
mode (SKILL.md Mode 12), not a solver bug.
