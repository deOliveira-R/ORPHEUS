---
name: container-ownership-dof-criterion
description: The decidable criterion for "what should this container own?" — does removing it move the DOF SET or the GRAM? — plus the admissibility-seat clause that decides whether the container survives at all. Derived on SNMesh (campaign 1), method-agnostic.
metadata:
  type: project
---

# "Who owns X?" is decidable — the DOF-set/Gram criterion + the admissibility seat

**Fact.** A container-ownership question ("should `SNMesh` hold the scheme? the
frame? everything?") is not a taste question. Two clauses decide it, and they
select the same set:

> **(i)** The container owns exactly what determines the **DOF SET and its
> GRAM** — the discrete measures and the basis refinements that change what the
> unknowns *are*.
> **(ii)** Of what survives (i), the container's only irreducible BEHAVIOR is
> adjudicating **mutual admissibility** of the choices before any space or
> operator exists.
> Everything else is a **function ON the space** (→ the medium/data layer) or a
> **binding of an operator TO the space** (→ that operator's construction).

**Operational form (greppable, no taste argument):** *remove the constituent —
did the SHAPE of the state vector or the METRIC `G` change?* YES ⟹ carrier data.
NO ⟹ operator data.

**Why:** derived 2026-08-21 attacking `SNMesh` (`scratch/sn_posing_ontology_and_names.md`,
every claim `file:line`-grounded on `feature/cs1-energy-space`). It decides in one
stroke what four separate arguments had decided four different ways:

| datum | verdict | by the criterion |
|---|---|---|
| angular quadrature | CARRIER | ordinates ARE the DOFs |
| spatial mesh | CARRIER | cells + volumes |
| scheme's **element** (`spatial_basis_per_axis`, `moment_mass_diagonal`, `θ`) | CARRIER | appends the `2^d` tail to ψ AND enters `G_bulk` |
| scheme's **closure** (`update`/`residual`/scan coeffs) | OPERATOR (streaming) | moves neither shape nor metric |
| the SH **frame** | OPERATOR (scattering) | a morphism BETWEEN representations; owned by the operator whose eigenbasis it is (L-009) |
| stencil, pole closure, sweep schedule, realized BC laws, gauge | OPERATOR | ditto |
| materials / `mat_map` | MEDIUM | functions ON the space |

⭐ **The answer to the standing symmetry attack** ("if the frame is not
mesh-owned, why is the scheme?"): frame and scheme are not comparable objects.
**Measures are carrier data; frames are operator data.** The frame's spatial
counterpart is the *closure* (and it should leave, as the frame did); the
quadrature's spatial counterpart is the *element* (and it should stay, as the
quadrature does). The symmetry is real; restoring it removes the closure, not
the element.

⭐ **The scheme is a conjunction, but it does NOT fission.** Its two faces share
one parameter (`θ` is the element's Legendre mass AND the closure's Galerkin
coefficient — `[M]` `linear_discontinuous.py:294`/`:334`/`:622`). Two objects
sharing one datum is Smell #16 shape 2. The correct move is **mint-and-forget**:
one scheme class, the space reads the element face and then forgets the scheme,
the streaming binding reads the closure face. Same shape as
`frame.conjugate` → the bound operator keeps `(M, Λ, R)` and forgets the Frame.
⚠ Forgetting is free only when the retained arrows are **closed under the
operations you will later ask for** — a Petrov-Galerkin binding forgets its
frame at the price of not being able to produce its own adjoint (CS4a-R / L-018).

⭐ **Whether the container survives AT ALL is the dependency-preservation
question** (relational normalization, Frame 4 of the attack). Group each
attribute by its determinant; the decomposition is lossless-join by construction
(every group shares `axes × quadrature`). The ONLY thing that can force a
relation to survive is a constraint spanning two decomposed relations — here the
four cross-choice refusals (geometry×rule, geometry×scheme, quadrature×geometry,
scheme×closure). ⟹ **the container is worth ~40 lines and one invariant, or it
is worth nothing.** There is no middle version that is not the God object again.

**How to apply:** on any "what should this hold / should it hold everything"
dispatch — (1) enumerate the container's ACTUAL contents from the tree before
accepting the requester's framing (on `SNMesh` the premise named 3 measures; the
tree held 6 further objects: a bound strategy, a realized-operator table, a
projector, a stencil, a traversal order, a stringly-typed geometry twin);
(2) run the remove-it-does-`shape`-or-`G`-move test per constituent; (3) compute
the determinant table and find the constraints that span groups — those name the
survivor's reason to exist; (4) if no constraint spans, say plainly that the
honest outcome is NO CLASS, and say it must be RULED rather than defaulted into
(a later session re-derives "there is no such object" as an omission).
