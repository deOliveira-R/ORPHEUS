# Literature Researcher Lessons

Read this at the START of every invocation.

---

## Cylindrical SN: Bailey is the primary reference

- Bailey, Morel & Chang (2009): **Eq. 50** = α recursion, **Eq. 74** = Morel-Montry weights. These are the exact equations ORPHEUS implements.
- Bailey's μ = radial cosine = ORPHEUS `mu_x` (our η). In slab/spherical contexts, μ means the *streaming* cosine instead. Always state which geometry when mapping μ.
- The ΔA/w geometry factor in the curvilinear balance equation is critical but poorly explained in textbooks. Bailey's asymptotic analysis (Section 4) is the clearest derivation of why it's needed. Lewis & Miller omit this.
