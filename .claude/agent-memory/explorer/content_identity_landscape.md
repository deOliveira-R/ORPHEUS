---
name: content-identity-landscape
description: Which input-layer types have content eq/hash, which raise, the byte encoders that exist, and the probes that measure them (P1 step 5, #405)
metadata:
  type: project
---

Tree of `1dc31163` (2026-10-01); a present-tense line is a claim to re-verify. Census in
`scratch/reference_architecture/p1step5/census.md`; probes beside it (`probe_types.py`,
`eqhash_spy.py` + `run_shards2.sh`).

- Content identity today: `Mixture` (bytes; -0.0 distinct; NaN admitted; explicit stored
  CSR zero distinct), typed laws (dataclass float `==`), `StructuredGeometry` (eq ok, hash
  raises only through a `BC`), `CellsBy*`, `CellEdges`, `MaterialMesh` (tuple keys).
- Raise or identity: `BC` hash (params dict, mutable even on `BC.vacuum`), `FaceLaws`
  (Mapping: unhashable, equals a plain dict), `Mesh1D` (unhashable, pinned by
  `test_unhashable_until_step_5`), `Mesh2D` (`==` RAISES on twins, aliases caller arrays),
  `Materials` (`eq=False`, and UNPICKLABLE: mappingproxy) and so `MaterialMesh` unpicklable.
- `VacuumInflow`/`ReflectiveBoundary` equal their kind STRING but hash differently.
- `Partition` and `Region` do not exist; `Mesh1D` stores coord/edges/volumes/mat_ids/face_laws.
- Byte encoders: `Axis._structural_bytes` (tagged, length-prefixed, repr for scalars so
  1 ≠ 1.0) and three unprefixed trace-space twins; Sood `_hash_params` (JSON, no consumer).
- Instrument lessons: a regex Compare census misses helper-returned receivers (`_mesh()`,
  `_two()`), so it is a lower bound; the runtime eq/hash spy is the census of record.
  `pytest -n` fails collection here (param ids carry addresses): shard by directory.
  `MaterialMesh._contractibility_key` reads `Mixture._identity_key` directly, bypassing the
  dunders, so a dunder spy does not see it. See [[census-predicates-bound-reference-and-activation-traceback]].
