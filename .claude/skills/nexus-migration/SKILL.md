---
name: nexus-migration
description: "Use when the user wants to plan a dependency migration, find all uses of an external library, or plan a technology switch. Examples: \"Plan numpy to jax migration\", \"Find all scipy usage\", \"What depends on numpy.ndarray?\""
---

# Dependency Migration with Nexus

IMPORTANT: This skill is the dedicated tool for migration planning.
It replaces Grep for dependency analysis and usage tracking.

## Workflow

```
1. migration_plan({from_dep: "numpy", to_dep: "jax"})  → Phased plan
2. Review phases: leaf first (safe), core last (high risk)
3. Execute phase by phase
4. After each phase: retest() → run tests
5. After all phases: staleness() → catch stale docs
```

## Key Tool

**migration_plan** — phased dependency migration:
```
migration_plan({from_dep: "scipy.special", to_dep: "custom_functions"})
→ Phase 1 (leaf): derivations._kernels.e3_kernel (0 upstream callers)
→ Phase 2 (mid): collision_probability.CPMesh._compute_slab_pij (3 callers)
→ Phase 3 (core): collision_probability.CPMesh.compute_pinf_group (14 callers)
→ Doc updates: theory/collision_probability.rst
```

See [../nexus-exploring/reference.md](../nexus-exploring/reference.md) for full reference.
