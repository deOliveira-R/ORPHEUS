"""The reference solutions: their readings, and later the solutions and certificates (#405 P2).

An input-tier package above :mod:`orpheus.specification`: it may import the
specification, :mod:`orpheus.numerics`, :mod:`orpheus.data` and
:mod:`orpheus.geometry`, and nothing above them (no mesh, no transport, no
method package, no derivations). The derivations write reference solutions,
so they import this package; production reads them through it.
"""

from orpheus.reference.reading import Printed, ReferenceReading

__all__ = ["Printed", "ReferenceReading"]
