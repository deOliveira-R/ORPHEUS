---
name: domain-op-l2-promotion-asymmetry-law
description: Documenting a domain operation born from an L2 promotion + data/behavior type split; the asymmetry-law contract pattern for a deferred sibling op
metadata:
  type: feedback
---

Recipe for documenting a **domain operation** (a transform on a
`Solution`/state object, NOT a solver step) that lands together with an
**L2 promotion** of supporting modules and a **data/behavior type
split**. Instance: `Solution.homogenize` (#267,
`refactor/operator-inverse-algebra`) — flux·volume-weighted spatial
homogenization → `MaterialMesh`; `axis`/`material_xs_field` promoted from
`orpheus.sn` to L2 `orpheus.transport.mesh`; `SNMesh(MaterialMesh)`.

**Rule:** a domain-op theory section is defined by *what it preserves*,
not by an averaging recipe. Lead with the preservation identity (give it
THE `:label:`), then DERIVE every weight from it. Pair it with an
**asymmetry-law `.. list-table::`** that contracts the deferred sibling.

**Why:** the page is the LLM's brain — a future session must be able to
(a) prove the op's claim, (b) understand WHY each weight choice is
forced not chosen, and (c) implement the deferred sibling against the
documented contract. See [[lessons-L004]] (vv-status on decomposition
labels), [[lessons-L010]] (V&V vocabulary).

**How to apply** (the 6-part section shape, all from the method
docstring as algebra-of-record — no SymPy module exists for a
domain-algebra op, the closed-form identity IS the source of truth):

1. **Preservation identity first, as THE verifies-target.** The defining
   property (here `Σ_{R,g}·Φ_{R,g} = Σ_{i∈R} V_iΣ_{i,g}φ_{i,g}`) gets a
   `:label:` the brief will bind a `@pytest.mark.verifies` to. Leave it
   UNTAGGED (orphan) — it is a real solver claim with a real L0 test,
   just not yet wired; tagging it `documented` would be wrong. Every
   OTHER eq-label in the section is a derivation-decomposition step →
   `.. vv-status: <label> documented` with a rationale comment naming the
   identity it serves (L-004). One verifies-target + N documented steps.
2. **Derive each weight from the identity — forced, not chosen.** Show
   the flux·volume weight `w=V·φ` falls OUT of solving the identity for
   Σ_eff; note any other weight breaks preservation at interfaces (the
   regime the op exists for). This pre-empts "why this average?".
3. **Flag the variable-swap trap as a `.. warning::`.** Matrix channels
   (`SigS[g_from,g_to]`) weight by the SOURCE group g_from — the
   out-scatter rate scales with the source population. Weighting by the
   sink group is vv Mode 2 (`SigS` vs `SigS^T`), invisible on
   flat-flux/single-material, caught only by a heterogeneous multi-group
   gate whose reference loop weights by g_from EXPLICITLY (L11 structural
   independence). Always name the Mode + why it hides + what catches it.
4. **Special-case the simplex channel.** χ is a probability simplex, NOT
   a rate — weight it by PRODUCTION rate (convex combo of simplices = a
   simplex). Cross-ref the existing multi-fissile χ_mix law as the
   analogue (`:eq:`emission-spectrum-chi-mix``), don't re-derive.
5. **One-line balance-preservation argument.** Every removal channel
   shares the one per-(R,g) weight → homogenized balance is the same
   convex avg of fine balances → zero, automatically. Note this is
   ANOTHER reason the weighting is forced.
6. **Asymmetry-law table for the deferred sibling.** When the op has a
   sibling deferred to a later slice (homogenize/space ↔ condense/energy),
   a `.. list-table::` contracting BOTH (collapses-what, mesh-coupling,
   return-type) IS the spec the sibling must satisfy when it lands. The
   structural reason for different return types (mesh-COUPLED →
   geometry+materials carrier; mesh-DECOUPLED → portable `dict[Mixture]`)
   is the load-bearing content — spell out WHY, not just WHAT.

**L2-promotion cross-ref hygiene:** a `git mv` from L3 `sn/` to L2
`transport/mesh/` makes EVERY `:class:`/`:mod:`/`:func:`/`:meth:` on the
old path stale. These were ALREADY plain-text (the modules were never
automodule'd — page convention), so `-W` is blind; fix on correctness
grounds via grep-gate across docs + code + TESTS (the retirement-audit
3-search blast radius — tests carry stale docstring paths too). New
domain-op refs (`Solution.homogenize`, `MaterialMesh`) ALSO render
plain-text by the same convention — match it, do NOT half-surface 1-2
leaves with a new automodule (L-002). Changelog row + a `:ref:` to the
new section closes the loop (create the labelled section in the SAME
edit or the intra-doc `:ref:` warns under `-W`).
