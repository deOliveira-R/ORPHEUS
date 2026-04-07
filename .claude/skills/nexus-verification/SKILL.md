---
name: nexus-verification
description: "Use when the user wants to check verification status, documentation coverage, or doc-code drift. Examples: \"What's verified?\", \"Which docs are stale?\", \"What equations have no tests?\", \"Documentation coverage report\""
---

# Verification & Documentation Quality with Nexus

## When to Use

- "Which equations are verified?"
- "What's the documentation coverage?"
- "Which docs are stale?"
- "What equations have no implementing code?"
- "What code has no tests?"
- V&V status assessment
- Documentation quality review

## Workflow

```
1. nexus verification_coverage()                    → Full equation → code → test map
2. nexus verification_coverage({status_filter: "implemented"})  → Gaps: code but no tests
3. nexus staleness()                                → Docs that drifted from code
4. nexus session_briefing()                         → Combined overview
```

## Checklist

```
- [ ] nexus verification_coverage() for full V&V status
- [ ] Review "implemented" entries — equations with code but no tests (verification gaps)
- [ ] Review "documented" entries — equations with no implementing code
- [ ] Review "orphan_code" entries — code with no equation (undocumented theory)
- [ ] nexus staleness() to find docs needing updates
- [ ] Create GitHub Issues for each gap found
```

## Coverage Status Values

| Status | Meaning | Action |
|--------|---------|--------|
| **verified** | Equation + code + test | Good — fully traced |
| **tested** | Code + test, no equation link | Add IMPLEMENTS documentation |
| **implemented** | Equation + code, no test | Write a test — verification gap |
| **documented** | Equation only, no code | Either implement or mark as future work |
| **orphan_code** | Code with no equation | Document the theory behind it |

## Tools

**verification_coverage** — the V&V map:
```
nexus verification_coverage({status_filter: "implemented"})
→ math:equation:alpha-cylindrical — has code (sweep_cylindrical) but no test
→ math:equation:surface-to-surface — has code (CPMesh._compute_slab_rcp) but no test
```

**staleness** — doc-code drift:
```
nexus staleness()
→ STALE: doc:theory/discrete_ordinates
→   sn_sweep.sweep_cylindrical modified 3 days ago, doc unchanged for 2 weeks
→   Affected symbols: sweep_cylindrical, compute_alpha_dome
```

**session_briefing** — combined overview:
```
nexus session_briefing()
→ 5 stale docs, 12 verification gaps, 452 unresolved references
→ Priority: update discrete_ordinates.rst
```
