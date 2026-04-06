# Literature Researcher Lessons

Read this at the START of every invocation.

---

## 2026-04-05 — Cylindrical SN formulation research

- **Most useful source**: Bailey, Morel & Chang (2009) — had the exact
  equations needed (Eq. 50 for α, Eq. 74 for M-M weights).
- **Notation pitfall**: Bailey uses μ for the radial cosine in
  cylindrical (our η = mu_x). Our code's mu_x/mu_y/mu_z map to
  different symbols in different references.
- **Key finding**: the ΔA/w geometry factor is present in the standard
  formulation but rarely explained in textbooks (Lewis & Miller omit
  the derivation). Bailey's asymptotic analysis is the clearest
  explanation of WHY it's needed.
- **Export control**: some production code manuals have the complete
  equations but cannot be referenced. Use only publicly available
  papers and textbooks.
