---
status: validated | candidate
project: <project name>
modules: [<module 1>, <module 2>, ...]
date_validated: <YYYY or "pre-chat">
primary_criteria: [<elegance criterion>, ...]
secondary_criteria: [<elegance criterion>, ...]
requires_validation: <true | false>
---

# Precedent: <short name>

## Status (candidate files only)

What validation would look like and why this hasn't been done
yet. Omit this section for validated precedents.

## Problem (current / minimal formulation)

Describe the current implementation / derivation. Include the
symptoms (elegance smells from reference.md Part C) that
suggest the wrong frame.

## Structural trigger

Name the triggers from reference.md Part A that justify the
frame candidate. Trigger names must match reference.md exactly.

## Frame

The mathematical frame being applied. Name it and briefly
identify the core objects from the frame (manifold, group,
measure, tensor network, etc.).

## Reformulation (sketch)

Concrete restatement of the problem in the frame's native
objects. "Apply frame X" is rejected. Specific mathematical
content is required.

## Elegance payoff

Assess against the four criteria from SKILL.md. Omit criteria
that do not fire. For each that fires, state specifically how.

## Concrete first test

Smallest experiment that discriminates the reformulation from
the current formulation. An experiment that cannot fail is
rejected.

## Literature path (optional)

References for the frame and for any existing partial
applications to this problem class.

## Transferable pattern

How to recognize when this reformulation applies in other
contexts. This is the critical section — it makes the
precedent useful beyond its specific instance.
