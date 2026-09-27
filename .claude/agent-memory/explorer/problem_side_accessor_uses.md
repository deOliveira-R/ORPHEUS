---
name: problem-side-accessor-uses
description: What the problem-side loss/production/fission/lhs/rhs reads are FOR (2026-09-26, HEAD 2b4d7424); the k estimator is the sides' Rayleigh quotient; the SN gauge straddles the sides
metadata:
  type: project
---

Census of 2026-09-26 (report `scratch/posing_sequence/accessor_census.md`, instruments beside it). Durable shape:
- Every Strategy use of `fission`/`production` needs only `−balance_derivative` (KEigenvalue already reads only `posing.pencil.rhs`).
- The method-tier SN k (fission / absorption+leakage−n2n) == `EigenPosing.rayleigh` on the sides to 2e-13 (Σ₂-carrying probe).
- The SN production gauge = IRR(νΣf) + n2n emission: it straddles the two sides. No side accessor serves it, and IRR alone does not. Absorption is not a layer-1 term, so a term selector alone does not serve it either. That makes two mechanisms (a channel-rate functional extended to emission channels, and the question's sides), not one.
- `_pair(1, ·)` on the SN carrier is a coordinate sum: ⟨1,Fψ⟩ = N_ord × fission rate.
- About half of the test reads of `problem.fission` are handles to the F TERM, not question reads.
**Why:** it feeds the #405 ruling "no physical accessor on a question".
**How to apply:** re-verify the line anchors before citing them; methods mean an AST census (Nexus callers returned 0, with 10 unresolved).
