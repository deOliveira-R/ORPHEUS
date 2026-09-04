---
name: scalar-to-stack-retype-rulings
description: Reviewing a data-layer retype "one matrix → a Legendre/index stack" (#426 step 1) — the three recurring smells (pad spelled once per layer with one spelling unreachable; the retype makes the empty stack representable; a rationale comment naming a structural invariant that is false on one shipped input) and the probes that measure each.
metadata:
  type: project
---

# Retyping a scalar data field into a STACK — the three smells, and how to measure them

Fact: reviewed the uncommitted production diff of #426 step 1 (`fix/n2n-anisotropy`,
2026-09-03) — `Isotope.sig2` / `Mixture.Sig2` went `csr_matrix → list[csr_matrix]`
over ℓ, the ingest stopped truncating every scattering channel at P2, and the HDF5
store went to format 2. The findings generalise to ANY "one value → an indexed stack"
retype (a moment stack, a per-order table, a per-mode list).

**Why:** this shape recurs whenever a lossy data layer is made lossless, and all
three smells are invisible to a green test suite — the tree was 313-passed green in
the data tree with all three live.

**How to apply:** run the three probes below before grading anything else.

## 1. The PAD RULE gets one spelling per consumer layer — count the WITNESSES

"An order the stack lacks is zero" is one domain statement. A retype spells it once
per layer that reads the stack. In this diff there were **four**: the macroscopic sum
(`_macroscopic_stack`, skip-term), the per-cell gathers (`_gather_legendre`/
`_gather_n2n`, `else np.zeros`), and a NEW early-return guard in `interp_sig_s`.

⭐ The measurement that decides which are real: **instrument each spelling and count
how many times it fires on the raggedest SHIPPED input.** Here, on borated water
(H/B stacks length 1, O length 7): `_macroscopic_stack` padded **12** terms per
build; the `interp_sig_s` guard fired **0 of 28** calls and is structurally
unreachable, because the one production call site loops `range(len(iso.sigS))`.
⟹ a guard that lands with no case it catches (`plan-authoring` §6c) — and worse, its
docstring claimed to BE the mechanism that makes the mixture sum work, which the
other spelling actually does (`coding-elegance` anti-pattern #20).

⚠ The direction: an unreachable pad reads as *defensive care*, and silently
returning a zero for an out-of-range index is the `getattr(x, "k", default)` hazard —
it fails in the DEFAULT's direction. A list's natural `IndexError` is the honest
refusal. Delete the redundant spelling; keep the one with witnesses.

## 2. The retype MAKES an illegal state representable that the scalar could not spell

A `csr_matrix` field cannot be *absent* — "no channel" is the zero matrix. A
`list[csr_matrix]` can be `[]`, and `Sig2[0]` is read unconditionally by every
P0 consumer. `[M]` `Mixture(Sig2=[])` constructs fine and dies with a bare
`IndexError: list index out of range` inside `absorption_xs`, miles from the
fixture; ragged blocks (`[(2,2), (5,7)]`) are accepted too.

⟹ a stack retype owes its type a `__post_init__` law: **non-empty, square, uniform
`ng` across every block of every stack.** Before demanding it, run anti-pattern
#18's leg (ii) — spy on `__post_init__` and exercise the in-tree factories to prove
the law refuses NO legal value (`[M]` here: 15 constructions, 5 distinct
(ng, len SigS, len Sig2) shapes, 0 violations).

## 3. ⭐⭐ A rationale comment naming a STRUCTURAL invariant is a claim over the whole corpus — and the counterexample is usually the input that MOTIVATED the comment

The sharpest find. A new comment justified a `+1e-30` epsilon as *"what keeps every
order's sparsity pattern equal to P0's, which the sigma-zero interpolation
assumes"*. `[M]` over all 13 shipped isotopes: true on 12, **FALSE on NA023**
(L0 nnz 16809, L1..L6 13196) — because NA023's MT=91 stores NL = 1, so its
positions enter L0 only. And the named consumer does not assume that at all:
`interp_sig_s` reads one pattern across the **σ₀ columns of ONE order**, which the
epsilon does protect and which holds everywhere.

So the code is CORRECT and the prose names the wrong invariant — the inverse of the
usual tell (§inst-knowledge #7). NA023 is the exact isotope whose short section
motivated the comment and the exact one carrying the diff's single bit-identity
exception.

⟹ **any comment of the form "X keeps property P true" is a universal; run it over
every shipped input before accepting it.** One loop over the store. The bug habitat
is specific: a maintainer either writes a gate asserting P (red on NA023, "fixed" by
re-adding the retired pad — reopening the truncation the whole step retires), or
trims the epsilon on ℓ≥1 "since only P0's pattern matters", silently breaking the
invariant the code DOES depend on.

## 4. The twin the diff CREATES by parameterising an already-padded gather

`_gather_sig2(self)` (no order arg, no pad) became `_gather_n2n(self, order)` — a
byte-for-byte copy of the pre-existing `_gather_legendre` modulo `SigS`→`Sig2`
(verified by AST-unparse comparison, both the gathers and the `_n_legendre*`
counters). The retype is what turns a non-twin into a twin; look for it whenever a
second channel is promoted to the shape an existing channel already had. The file's
own `_gather_vector(attr: str)` was already the collapse target.

## 5. Reusable probes (all ran in seconds, all in-process)

- AST-unparse two method bodies with the channel name normalised → literal-twin proof.
- `ast.walk` the module for the guard's call sites and read the loop BOUND → reachability.
- Monkeypatch the producer + the guard, build the raggedest shipped mixture, count fires.
- Spy `__post_init__`, exercise the factories, print the distinct shape tuple set.
- Open every `*.h5` in the store and print the per-channel order depth — this is what
  proved F1(b) landed (sigS NL=7 on 13 of 13; sig2 NL=7 on 9, NL=1 on B_010/H_001/NA023).

Related: [[symmetry-realization-carve-rulings]] (the guard does not move with the
mathematics), [[coupled-block-boundary-unweld-rulings]] (twin delivery single-sourced
at the operator).
