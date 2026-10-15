---
name: c1-three-types-content-identity-census
description: PanelBasis/ConcentricPartition/RigidMotion before ContentIdentity (#405 P3 C1) — RigidMotion already value-equal and aliases caller arrays; 0 runtime eq/hash in production; S5.9 roster gate reds on any new ContentIdentity class
metadata:
  type: project
---

Tree `9e903a4d` (2026-10-14); a present-tense line is a claim to re-verify. Census
`scratch/characteristic_architecture/p3/c1/census.md`, probes beside it (`eqhash_spy.py`,
`mint_counter.py`, `run_spy.sh`).

- `RigidMotion` is `eq=False` but owns an exact `__eq__`/`__hash__` (`_exact_key`): NOT identity-compared.
  It stores the caller's array by reference (`.linear is Q`), so a freeze must copy first.
- Runtime spy, all trees not-slow: 0 production eq/hash/pickle of the three; 43 test-side motion `==`
  (deck holders' generated eq). Symmetry code uses tolerance (`close_group` rounded keys,
  `approximately_equals`), never motion `==`.
- Hidden gate: `test_content_identity.py` S5.9 walks `orpheus.geometry` and `...characteristic`; any new
  `ContentIdentity` class reds until it has a roster `Entry` and takes the mixin's dunders.
- Instrument lessons: `cls.__bases__ = (ContentIdentity,)` is refused by CPython (deallocator differs) — simulate
  with a stand-in subclass with the mixin first in the MRO. Measure "today" BEFORE any patch in the same process.
  `_ConstructorPickler` already sends every non-ContentIdentity dataclass through its constructor.
