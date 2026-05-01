"""Path 1 — discrete-solver derivations.

Each sub-package mirrors a production solver in :mod:`orpheus` and
holds the **symbolic discretisation** that solver commits to. These
are *not* references — they are the explicit form of the discrete
equations that the production code must satisfy.

The references they are verified against live in
:mod:`orpheus.derivations.continuous`.
"""
