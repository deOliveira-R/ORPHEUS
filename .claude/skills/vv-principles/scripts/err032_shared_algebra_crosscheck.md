# ERR-032 — Procedural independence is not structural independence

**The bug.** The slab white-BC analytical reference shipped with a
factor-of-two algebra error: the partial-current integral was
reduced using the identity `∫E₂ = 1 − E₃` when the correct
antiderivative is `∫E₂ = ½ − E₃`. The algebraically correct closed
form collapses to `φ ≡ 1/Σ_t` (the Wigner-Seitz identity); the
buggy form did not.

**Evidence that existed before it was caught.** Two derivations of
the white-BC flux — an analytical closed form and a fixed-point
iteration of the partial-current balance — agreed to 1e-39. They
were written independently, in different code paths, by different
reasoning chains. By every conventional standard of "two-source
verification," the formula was confirmed.

**Why that evidence didn't catch it.** Both derivations reduced the
same class of volume integrals using the same antiderivative table
from the same human's working memory. They were *procedurally*
independent (different code, different control flow) but
*structurally* coupled (both applied `∫E₂ = 1 − E₃` at the same
algebraic step). Two derivations of the same closed form, no
matter how independently coded, share whatever upstream identity
reduces the integrand. Agreement to 1e-39 reflected the shared
identity, not the correctness of either branch.

**What evidence class would have caught it.** Cross-checks **MUST**
be *structurally* independent — different integrand or different
identity, not the same identity expressed in two control flows.
**NEVER** accept agreement between two derivations of the same
closed form as verification — **instead**, force the cross-check
to come from a different structural angle: pair a kernel-based
check (row-sum of `K_ij`, particle balance, conservation identity)
with a closed-form check (eigenvalue, asymptotic limit). The
row-sum identity that ultimately exposed ERR-032 — `Σ_j K_{ij}·1
= (1/(2 Σ_t))[2 − E₂(τ_i) − E₂(τ_i')]` — was derived from the
*kernel* itself, not from a volume integral, and so used the
correct `½ − E₃` identity by construction. The factor of ~2.2
between the row-sum and the buggy φ was the first-order
fingerprint. See vv-principles reference.md §1 (procedural vs
structural independence is the canonical worked example here).

**References.** ERR-032 in `error_catalog.md`; main agent
lessons.md L11; commit `2538cfe` (introduced); the row-sum L1
test `tests/derivations/test_peierls_reference.py::TestSlabKernelRowSum`
(`@pytest.mark.catches("ERR-032")`).
