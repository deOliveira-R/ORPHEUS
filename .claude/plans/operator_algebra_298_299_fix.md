# operator_algebra #298 + #299 — correctness fixes (NEXT phase, queued 2026-07-20)

**Status:** ✅ DONE (2026-07-20). **#298 fixed @ `639cec9e`** — the block matrix's (s,s)
entry corrected to `A_ss` + the trace-resolved refinement equation
(`bc-extraction-trace-blocks`: L_full unit-lower-triangular `[[I,0],[−T,I]]` vs B
strictly-upper `[[0,R],[0,0]]` on `V_inflow ⊕ V_outflow`; the identity diagonal =
triangularity = sweep-as-forward-substitution, verified against `streaming.py` O.4a.2 +
`boundary.py`); matrix.rst regenerated (321 documented). **#299 fixed @ `4425516c`** — the
audit widened the issue's 3 sites to **5** (+ the published `green_operator.py` /
`operator.py` docstrings) + the boundary-law census correction (reflective/albedo are
angular permutations rank N/2, NOT rank-1 — the old prose contradicted its own G_α table);
canonical statement at the `smw-low-rank-exception` anchor, compact pointers elsewhere.
Both commits gated `-E -W` green (0 warnings). **Phase H's X2/X3 blockers cleared** (X1 was
`018ecb7b`). #298/#299 close on merge via `Closes` trailers (branch unpushed); breadcrumb
left on #300.

> ⚠ **Content moved during the reframe.** Both issues cite the pre-reframe
> `docs/theory/operator_algebra.rst:<line>`; that page is now
> `docs/theory/foundations/operator_algebra.rst`, and some content SPLIT to deep-dive pages.
> The CURRENT locations below were verified 2026-07-20.

---

## #298 — `L_full`'s (s,s) block is written `0` but the table says identity

**Now lives on `docs/theory/foundations/boundary_conditions.rst`** (the augmented-operator / BC
block algebra moved off `operator_algebra` during the skeleton/reframe):

- **L423–431** — three block matrices: `[[A_bb, 0],[0,0]]`, `[[A_bb, A_bs],[A_sb, 0]]` labelled
  **`L_full (FULL)`**, and `[[0,0],[0,A_ss]]`. The `L_full` **(s,s) block reads `0`**.
- **L447–458** — the table: "`L_full` → `A_bb` (streaming) + `A_bs` (inflow seeds the …) +
  **outflow row keeps the self-consistency defect** + **inflow row carries the identity**
  `I·ψ.inflow`".
- **L549–556** — a "the in-flight plan prose said `L_full` should emit the *raw* streamed outflow —
  this is **wrong**" note (bears on the outflow-row structure; read it before touching the matrix).

**The defect + the CARE POINT (do NOT do a naive `0 → I`):** the block matrix's (s,s) `0`
contradicts the table's "inflow row carries identity `I`". But the (s,s) block is **not uniformly
`I` or `0`** — per the table it is **identity on the inflow sub-row, the streaming
self-consistency defect on the outflow sub-row**. The fix is to make the (s,s) block faithfully
represent that inflow-identity / outflow-defect structure (the identity is load-bearing: it is what
makes the trace rows diagonal ⟹ the augmented operator **triangular** ⟹ sweep = forward
substitution — exactly what the root page will derive). READ L420–560 carefully; determine the
correct (s,s) structure from the CODE — `numerics/coupled_system.py::build_within_group_system` +
the `CoupledOperator` block algebra (`coupled_block_operator.rst` records a measured `A_ss = 5.0`,
`A_bs = 7.5`, forward-substitution #284 certificate) — reconcile the matrix with the table, and
**check whether any downstream triangularity prose derived a conclusion from the `0`.**

## #299 — Sherman–Morrison–Woodbury dismissed blanket (true of C, FALSE of B)

**Now in THREE places** (the reframe's splits spread the claim — consider stating it ONCE
canonically and `:ref:`-ing from the others, per single-source-of-truth):

- **`operator_algebra.rst:762–763`** — "SMW applies only under low-rank structure, which the dense
  collision diagonal is not" (the original).
- **`operator_algebra.rst:1457`** — a second SMW mention (Composition-algebra region, near the
  P4-T1 fold).
- **`operator_inverse_family.rst:397`** — a third instance (split to the inverse-family deep-dive).

**The fix (per #299's suggested wording):** scope each dismissal to **`C`** (and `L`), and add that
**`B` IS low-rank** — rank-1 per face for white / albedo (one isotropic re-entry mode); an ordinate
permutation (rank N/2) for specular 1-D slab — and that **CP already exploits exactly this**
(`cp/solver.py:395`; see **#300**, the Woodbury-borrowing algorithmic win this sentence currently
blocks). Apply at all three sites (or one canonical + refs).

## Gate + close
- `-E -W` clean each commit; run the whitespace-flattened consuming-prose sweep if any label moves.
- Trailers: `Closes #298` / `Closes #299`.
- **On landing: Phase H unblocks** (X2/X3 clear; X1 already fixed). Phase H still also wants Phase I's
  literature survey as an input before it is authored.

## Provenance
Both surfaced by an adversarial cross-domain-attacker structural review (2026-07-14) while
verifying `documentation_corpus_architecture.md` §3.6 (the root-page unification). Companion:
#300 (the Woodbury borrowing #299 unblocks). Sibling already-fixed defects from the same review:
`eigenvalue.py` docstring @ `018ecb7b`; the corpus-wide stale `(L+C−S−F/k)` spelling @ `0ca0d378`.
