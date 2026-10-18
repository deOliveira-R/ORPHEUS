---
name: verification-graph-claim-trace-frames
description: W5 2026-10-18 attack on verification_tree_refinement.md — claim = statement+mentions (Lean), rests_on welds executes/reads edges (Shake), certificate needs argument nodes (GSN), memo ChildPin early cutoff lost
metadata:
  type: project
---

New problem class: a VERIFICATION GRAPH (tests, claims, verdict cache, certificate). First sighting of each frame
in this class; nothing proposed to the skill yet. Memo: `scratch/verification_tree/pass2/cross_domain_attack.md`.

- **Lean statement/proof split**: claim = (predicate, domain, band) with a statement digest; labels are the
  constants it MENTIONS (108 of 553 function-level `verifies` functions name ≥2 labels, [M] AST census). An edge
  pins the lower statement digest. A loosened band is an "assumption changed" state. `#print axioms` gives the
  certificate's sorry/trusted-base query.
- **Shake à la carte**: the memo is a constructive trace whose `ChildPin` keys on payload digest, but `lookup`
  does not bring the child up to date before judging the parent, so early cutoff is lost. Entries store the key
  digest, not the arguments. Early cutoff on a verdict bit is useless across executes-edges: cone(upper) ⊇
  cone(lower).
- **GSN / Assurance 2.0**: dependency (diagnostic) ≠ inference (supports); kill row = confidence argument, which
  goes stale separately from the verdict.
- **Kill matrix = formal context** (Galois connection): mutant subsumption reads columns, test dominance reads
  rows.

**Why:** second sightings (another verification or test-suite design brief) go to the skill.
**How to apply:** on any verdict-cache, DAG-of-tests or certificate brief, measure labels per test and the type
of each edge (executes vs reads) first. Related: [[reference-object-equation-tower-frames]].
