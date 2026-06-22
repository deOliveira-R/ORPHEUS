---
name: Plan-cluster triage with built-in close-out comments
description: When a cluster of plans (4 quadrature + 2 specular here) has already been close-out-commented on its umbrella issues, the triage decisions collapse to DELETE for almost everything; the rare exception is unshipped tooling that was never tracked, which needs a fresh issue
type: feedback
---

When 6 plans are clustered around a topic and the umbrella issues already have post-cleanup §22.9-style close-out comments, the **default decision is DELETE for all plans whose content the close-outs already capture**. Verify by reading the close-out comment top-to-bottom against each plan's acceptance criteria — when the audit table maps every plan-§ to a shipped commit, the plan is fully archived on GitHub already.

**Why:** Re-archiving in scratch/archive directories adds noise without preserving anything new. The §22.9 audit comment for #133 (and the per-issue rows mirrored on #134/#135/#136) was specifically designed during the post-#138 cleanup to be the durable record. Plans that match the audit are by-definition redundant.

**How to apply:**

1. **Read the umbrella close-out comment first** (here: #133's §22.9 audit comment). It has the per-primitive landing map + commit chain + acceptance-criterion table.
2. **For each plan, check** if every acceptance criterion has a "**DONE — commit X**" row. If yes → DELETE. The plan even sometimes embeds its own "Status: COMPLETE" header confirming this (see `visibility-cone-substitution-rollout.md`).
3. **For plans whose work was abandoned** (e.g., specular Phase 5 / continuous-µ / method-of-images), check if the close-out narrative comment captured the retreat. If yes → DELETE. The image-series approach for sphere is captured in `Round 1 Front A` of the #133 retreat comment; the cylinder image-series 2-D lattice is documented as never-attempted; the mode-space `R_specular = (1/2)M⁻¹` derivation shipped in §peierls-specular-bc.
4. **The rare exception is genuinely unshipped tooling** with no GitHub tracking. For `precision-floor-tool.md` (T3.4 / G-P2.1), no implementation existed in code, no issue had been filed. Cardinal Rule 3 says NEVER let an improvement opportunity pass undocumented — file a new issue with the plan content as the body, then delete the plan. Use the `Provenance` footer to point back at the original plan file path so future readers can trace the lineage.

**Observed pattern**: 4 of 4 quadrature plans + 2 of 2 specular plans were DELETE. The single exception (precision-floor) became a new issue (#145). The triage rubric here is more skewed than the 4-plan peierls cluster (which had several UPDATE EXISTING COMMENT actions) because the §22.9 audit was unusually thorough — every plan acceptance criterion was already row-mapped to a commit.

**Cross-check before deletion**: `grep -rn <plan-key-concept> orpheus/ docs/ tests/` — if nothing implements the concept and no issue exists, file the issue. If implementations exist, the plan is shippable-redundant and DELETE is safe.
