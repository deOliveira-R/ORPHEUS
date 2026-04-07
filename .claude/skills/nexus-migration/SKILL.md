---
name: nexus-migration
description: "Use when the user wants to plan a dependency migration, find all uses of an external library, or plan a technology switch. Examples: \"Plan numpy to jax migration\", \"Find all scipy usage\", \"What depends on numpy.ndarray?\""
---

# Dependency Migration with Nexus

## When to Use

- "Plan the numpy → jax migration"
- "Find everything that uses scipy"
- "What would it take to replace matplotlib?"
- "Show me all numpy.ndarray usage"
- Planning technology switches or dependency upgrades

## Workflow

```
1. nexus query({text: "numpy", node_types: "external"})       → Find all external dependency nodes
2. nexus migration_plan({from_dep: "numpy", to_dep: "jax"})   → Phased migration plan
3. Review phases: leaf functions first, core last
4. Review doc_updates: documentation pages that reference the dependency
5. Execute phase by phase, running nexus retest() after each
```

## Checklist

```
- [ ] nexus migration_plan() for the full plan
- [ ] Review Phase 1 (leaf functions) — safe to change first, no downstream deps
- [ ] Review Phase 2 (mid-level) — moderate blast radius
- [ ] Review Phase 3 (core) — high blast radius, change last
- [ ] Note documentation pages that need updating
- [ ] Execute phase by phase
- [ ] After each phase: nexus retest() and run tests
- [ ] After all phases: nexus staleness() to catch stale docs
```

## Tools

**migration_plan** — the primary tool:
```
nexus migration_plan({from_dep: "scipy.special", to_dep: "custom_functions"})
→ Phase 1 (leaf): derivations._kernels.e3_kernel (0 upstream callers)
→ Phase 2 (mid): collision_probability.CPMesh._compute_slab_pij (3 callers)
→ Phase 3 (core): collision_probability.CPMesh.compute_pinf_group (14 callers)
→ Doc updates: theory/collision_probability.rst, api/collision_probability.rst
```

**query with type filter** — find all external nodes for a package:
```
nexus query({text: "numpy", node_types: "external"})
→ numpy.ndarray, numpy.array, numpy.linalg.solve, numpy.zeros, ...
```

**impact** — blast radius of a specific dependency:
```
nexus impact({target: "py:class:numpy.ndarray", direction: "upstream", edge_types: "type_uses"})
→ 298 functions use numpy.ndarray in their type annotations
```
