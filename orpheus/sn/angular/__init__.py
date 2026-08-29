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
* :mod:`~orpheus.sn.angular.redistribution` — the **member-independent**
  half: the :math:`\alpha` dome with its admission contract, the
  :class:`~orpheus.sn.angular.redistribution.AngularRedistribution` factor and
  its single producer.  Every closure member reads it; none of them owns it.

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

✅ **The angular factor's OTHER half arrived 2026-08-28 (P4.2)** — see
:mod:`~orpheus.sn.angular.redistribution` in the list above.  It could not
move in P2, and the reason is worth keeping: its three callers, the
``*_streaming`` factories, were still in ``geometry/``, so moving it forced a
module-scope ``geometry -> sn`` import — a declared ``FORBIDDEN_EDGES``
violation *and* a measured circular-import failure of ``import
orpheus.geometry``.  ⛔ The same refutation landed a second time against the
order ``P4.2 -> P4.4``: the callers are *intra-file* references, so the split
itself creates the edge and no import grep can see it coming.  P4.4 moved the
factories to :mod:`orpheus.sn.mesh.reduced_operator` first, which turned the
call into ``sn -> sn`` and made this landing trivial.
"""

from .closure import (
    IdentityAngularClosure,
    MorelMontryAngularSweep,
    AngularClosureBase,
    default_angular_closure_class,
)

__all__ = [
    "IdentityAngularClosure",
    "MorelMontryAngularSweep",
    "AngularClosureBase",
    "default_angular_closure_class",
]
