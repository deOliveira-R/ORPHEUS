r"""Shared primitives for the singular-eigenfunction expansion family.

Currently empty — the cylinder is the only inhabitant of the
:mod:`...singular_eigenfunction` package family. This subpackage exists as
the reservation point for primitives that surface when slab or sphere
variants (Mitsis 1963, Smith 1969) are added.

The dispersion-function root primitive
:func:`...fn_method.core.dispersion.case_nu0` is shared from the F_N
package; this is acceptable cross-package reuse because the dispersion
function :math:`\Lambda(\nu) = 1 - c\nu\,\mathrm{atanh}(1/\nu) = 0`
is a **medium property** common to every singular-eigenfunction
treatment of the 1G isotropic transport equation, not an
algorithmic-pillar primitive of either package.
"""
from __future__ import annotations
