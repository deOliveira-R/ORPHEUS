"""The reference specification: materials, geometry and a question, as one keyed value (#405 P1 step 8).

An input-tier package: it composes :mod:`orpheus.data`, :mod:`orpheus.geometry`
and :mod:`orpheus.numerics`, and imports nothing above them (no mesh, no
transport, no method package, no derivations), so the derivations build a
specification and production consumes it.
"""

from orpheus.specification.specification import Specification

__all__ = ["Specification"]
