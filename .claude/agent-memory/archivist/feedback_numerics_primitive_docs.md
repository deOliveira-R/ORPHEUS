---
name: Numerics primitive docs (Wave 0 SN performance plan)
description: Patterns for documenting generic numerics primitives that lift from per-solver code into orpheus/numerics/ — naming-hierarchy rationale, cross-method consumer tables, tensor-product vs numpy-primitive distinction
type: feedback
---

When documenting a Wave 0-style lift of solver-internal helpers into a generic
`orpheus/numerics/` module, follow this template (galerkin_projection.rst,
spherical_harmonics.rst, operator_algebra.rst tensor-product section):

**Rule**: A theory page for a generic numerics primitive MUST contain (1) the
mathematical content (full derivation, convention, identity), (2) a
cross-method consumer table (every place the primitive is or will be consumed,
with status: live / pending), (3) the naming-hierarchy rationale when an ABC
chain encodes a discipline (e.g. `Galerkin` vs `Petrov-Galerkin`), and (4) a
"what was tried and discarded" history if the primitive replaces an earlier
implementation.

**Why**: The user's repeated steering — "galerkin projection is generally
useful — same machinery is used in cross-section energy condensation. Catalog
everything that needs to be implemented and where" — drove the cross-method
consumer table. Without it, future sessions implementing PN / energy-
condensation rebuild the primitive instead of reusing it. The naming-hierarchy
rationale was demanded explicitly: "why HarmonicMomentProjection(GalerkinProjection)
rather than just HarmonicMomentProjection(ProjectionOperator)" — the type
signal carries the discipline without docstring lookup.

**How to apply**:
- New theory page per primitive (or related primitive cluster): one for
  `galerkin_projection`, one for `spherical_harmonics`, extend existing
  `operator_algebra` for the tensor-product section. Don't conflate; each has
  its own audience.
- `tensorial-framing` label as section anchor in `operator_algebra.rst` is
  the canonical destination for "tensor product of operators" cross-refs from
  `discrete_measures.rst` partition_by, `discrete_ordinates.rst` Pℓ Galerkin
  note, etc.
- Numpy relationship: explicitly contrast np.einsum / np.kron / np.tensordot
  ("array primitives — implementation layer") vs TensorProductOperator
  ("operator-algebra type — algebra layer"). Use a 2-row list-table to make
  the layering visible.
- DO NOT add `automodule:: orpheus.numerics.X` directives to api/numerics.rst
  for modules that have rich `:label:` docstrings. The labels collide with
  the theory-page labels (Sphinx "duplicate label" warning). Cross-link via
  `:mod:`/`:class:` instead, plus a paragraph in api/numerics.rst pointing to
  the theory pages. (Retroactive lesson: pre-existing docstring-embedded
  `:pydata:` roles in operator.py also surface as ERROR when automodule
  renders them — work that is out of scope for a doc-extension session.)
- Title overline length: count display chars (`echo -n "..." | wc -c`) and
  match `=` count exactly. Off-by-one is silently failing or warning.
