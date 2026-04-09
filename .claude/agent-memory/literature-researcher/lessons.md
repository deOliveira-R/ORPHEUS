# Literature Researcher Lessons

Behavioral corrections only. AGENT.md has source priorities,
notation mapping, and extraction procedure -- never duplicate here.

---

## L-001 -- Bailey is the canonical cylindrical SN reference

Bailey, Morel & Chang (2009) NSE: **Eq. 50** = alpha recursion,
**Eq. 74** = Morel-Montry weights. These are the exact equations
ORPHEUS implements for curvilinear discrete ordinates.

**Notation trap**: Bailey's mu = radial cosine = ORPHEUS `mu_x`.
But in slab/spherical contexts, mu means the *streaming* cosine.
Always state which geometry when mapping mu.

**Rule**: When asked about curvilinear SN balance equations or the
DeltaA/w geometry factor, start with Bailey Section 4, not Lewis
& Miller (who omit the asymptotic analysis that explains *why*
the geometry factor is needed).
