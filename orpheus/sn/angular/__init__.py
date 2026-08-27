r"""The angular half of the discrete-ordinates discretization.

The curvilinear redistribution operator **factors** —

.. math::

   \mathcal{R} \;=\; R_{\rm spatial} \;\otimes\; A_{\rm angular}(\tau,\ \alpha,\ w)

— a spatial pairing against an angular operator, and this package is the
second factor.  It exists because the factorization *is* the ownership
map: the measure comes from the chart, the basis from the discretization
scheme, and the angular partner from the quadrature weights.  A module
that closes the angular axis therefore belongs beside its siblings on
that axis, not inside the package that walks the spatial one.

* :mod:`~orpheus.sn.angular.closure` — the angular-redistribution
  closure family (Morel--Montry weighted-DD, the neutral identity), and
  the :math:`\tau` / half-angle-grid machinery it is built from.

History: this package was created 2026-08-26 (the un-weld arc, P2).
``closure`` lived at ``orpheus.sn.angular.closure`` until then
— a path that carried two falsehoods.  **"pole"** named the special case
(the pole cell) for a family that closes the *whole* angular axis, and
:class:`~orpheus.sn.angular.closure.IdentityAngularClosure` never sees a
pole at all; **"sweep"** is the traversal package
(:mod:`orpheus.sn.sweep` — ``scan``, ``cache``, ``pairing``), and a
closure is a *discretization* object the traversal consumes, not a
traversal artefact.  The move is a pure relocation: no symbol was
renamed and no value changed.

⚠ The angular factor's OTHER half — the :math:`\alpha` dome and the
starting direction (``AngularRedistribution``, ``angular_redistribution``,
``alpha_dome``) — still lives in :mod:`orpheus.geometry.reduced_operator`
and joins this package at P4.  It could not move here in P2: its three
callers, the ``*_streaming`` factories, stay in ``geometry/``, so moving
it would force a module-scope ``geometry -> sn`` import and a measured
circular-import failure of ``import orpheus.geometry``.  P4 dissolves
those callers first, which is why the two halves land in that order.
"""

from .closure import (
    IdentityAngularClosure,
    MorelMontryAngularSweep,
    PoleAngularClosureBase,
    default_angular_closure_class,
)

__all__ = [
    "IdentityAngularClosure",
    "MorelMontryAngularSweep",
    "PoleAngularClosureBase",
    "default_angular_closure_class",
]
