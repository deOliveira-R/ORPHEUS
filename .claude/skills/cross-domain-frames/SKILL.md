---
name: cross-domain-frames
description: 'Proactively use whenever evaluating a first-pass implementation, derivation, proof, or method choice to detect whether it matches the problem''s native mathematical structure. Provides the trigger table mapping structural features (curvilinear coordinates, rotational symmetry, compositional boundaries, stochastic sampling, fixed-point iteration, integral operators, scale separation) to mathematical frames (differential geometry, group theory, tensor networks, Feynman-Kac, spectral theory, number theory, harmonic analysis, asymptotic analysis, etc.), the cross-method pollination map, the elegance detector, and validated precedents from ORPHEUS. Examples: "Is there a cleaner formulation", "Review this derivation for native structure", "What foreign frames apply", "Is this method missing borrowings from adjacent methods", "Which mathematical structure does this problem actually have"'
---

# Cross-Domain Frame Detection

IMPORTANT: Frame detection is NOT critique. This skill matches
structural features of a problem to entries in a fixed trigger
table. A match is a match or it is not. There is no "balanced
view" on whether a trigger fires.

## When to Use

- After a first-pass implementation — the minimal-effort path
  rarely matches native structure
- After a derivation — to check whether a cleaner frame exists
- When choosing between methods — to detect cross-method
  borrowings
- When multiple problem variants share structure that is not
  yet factored out
- When the elegance detector (Part C in reference.md) fires
  two or more smells

## When NOT to Use

- Numerical bug diagnosis → numerics-investigator
- V&V coverage auditing → qa or test-architect
- Literature lookup → literature-researcher
- Code exploration → explorer

## Workflow

1. Extract structural features → enumerated list
2. Trigger table lookup (Part A) → frame candidates
3. Cross-method pollination (Part B) → borrowing candidates
4. Elegance detector (Part C) → smell count
5. Reformulate per candidate → concrete sketch + first test
6. Assess elegance payoff → four-criteria score
7. Produce structured output → see "Output Shape" below
8. Update memory → grow or sharpen the library

## Four Elegance Criteria

Every reformulation candidate MUST be assessed against these
four. A reformulation that hits zero criteria is rejected output.

| Criterion                 | Fires when                                                                                                                                            |
| ------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------- |
| **Structure-exposing**    | The reformulation makes a previously-hidden mathematical structure explicit (e.g., redistribution term IS a connection coefficient, not a correction) |
| **Expressive**            | The reformulation states the same content more compactly, or composes naturally where the native form requires special cases                          |
| **Structurally-simpler**  | The reformulation reduces the number of distinct objects, special cases, or ad-hoc constants                                                          |
| **Algorithmic-advantage** | The reformulation changes asymptotic cost, enables parallelism, avoids matrix formation, or improves convergence rate                                 |

Elegance is strongest when multiple criteria hit simultaneously.
A reformulation that hits only one criterion is weak. Two is
standard. Three or four is a high-value match worth pursuing
immediately.

## Tables (in reference.md)

See [reference.md](reference.md):

- **Part A: Mathematical field → structural trigger table**
  - A.1 Geometry and topology
  - A.2 Algebra and representation
  - A.3 Analysis (functional, spectral, harmonic)
  - A.4 Probability, stochastics, information
  - A.5 Optimization, dynamics, control
  - A.6 Asymptotic and multiscale
  - A.7 Integral equations and special structure

- **Part B: Cross-method pollination map**
  - Pollination sources for CP, MOC, MC, SN, diffusion,
    eigenvalue iteration, sensitivity, coupled multiphysics,
    resonance self-shielding

- **Part C: Elegance detector**
  - Structural and semantic smells that indicate the native
    frame has not been found

## Precedents (in scripts/)

Validated matches (confirmed in ORPHEUS work):

- [scripts/validated_bc_tensor_network.md](scripts/validated_bc_tensor_network.md)
  — Boundary conditions as a tensor network; avoids dense matrix
- [scripts/validated_unified_geometry.md](scripts/validated_unified_geometry.md)
  — Slab / annulus / hollow sphere as one parameterized kernel

Candidate matches (high-prior, not yet validated in ORPHEUS):

- [scripts/candidate_cylindrical_connection.md](scripts/candidate_cylindrical_connection.md)
  — Cylindrical α-redistribution as a connection coefficient
- [scripts/candidate_feynman_kac.md](scripts/candidate_feynman_kac.md)
  — CP / MOC / MC as discretizations of one path integral

Template for new precedent files:
[scripts/\_template.md](scripts/_template.md)

## Output Shape

STRUCTURAL FEATURES

feature 1
feature 2
...

ELEGANCE DETECTOR HITS

smell N: one-line explanation
...
(or: no smells fired)

FRAME CANDIDATES
Frame: <name from reference.md Part A>
Trigger: <structural feature from enumerated list>
Reformulation: <concrete sketch — native objects, not gesture>
Elegance payoff:

Structure-exposing: <specific structure made explicit>
Expressive: <specific gain in expressiveness>
Structurally-simpler: <specific simplification>
Algorithmic-advantage: <specific asymptotic / computational gain>
(list only criteria that fire; omit ones that do not)
First test: <smallest experiment that discriminates from current>
Structural attack on current: <what the frame exposes that the
native formulation misses>
Precedent: <link to scripts/ if a similar match exists>

(repeat per frame)
CROSS-METHOD POLLINATION
Current method: <CP / MOC / MC / SN / diffusion / eigenvalue /
sensitivity / coupled / resonance>
Borrowings:

From <adjacent method>: <what it brings, trigger, first test>
...

UNEXPLORED

<frame> — <one-line reason no trigger matched>

...

A response missing any required section is a failure. Report
the gap explicitly rather than padding.

## Growth Protocol

The skill is not static. Each invocation should grow or sharpen
the library.

**When a frame-problem match produces concrete reformulation payoff:**

1. Create a new file in `scripts/validated_<name>.md` using
   the template
2. If the trigger that fired is new, add it to the relevant
   subsection of reference.md Part A
3. If the problem class is new, add it to reference.md Part B

**When a frame-problem match is speculated but produces no
concrete payoff:**

1. Add the rejection to agent memory under "checked and
   low-signal for <problem class>"
2. Do NOT add it to the trigger table — the table grows only
   from confirmed payoffs

**When an elegance smell is discovered that is not in Part C:**

1. Add it to Part C with a specific example that shows the
   smell pointing to a real reformulation
2. Cross-link from the reformulation precedent in scripts/

The library grows from evidence, not speculation. Speculative
entries degrade the signal of the table and lead to
noise-output from the cross-domain-attacker agent.
