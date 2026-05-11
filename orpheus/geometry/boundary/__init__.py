r"""Tensor-decomposed boundary conditions for transport solvers.

A boundary condition on a transport interface is, in full generality,
a linear operator :math:`R` mapping the *outgoing* angular flux at the
boundary to the *incoming* angular flux:

.. math::

    \psi_{\text{in}}(\Omega)
    = (R \psi_{\text{out}})(\Omega)
    = \sum_\alpha \bigl(G_\alpha \, \psi_{\text{out}}\bigr)(\Omega)
      \cdot A_\alpha,
    :label: bc-tensor-decomposition

where each term :math:`G_\alpha \otimes A_\alpha` is the tensor product
of a **geometric operator** :math:`G_\alpha` (a permutation, a
pushforward, an angular average, a spatial wrap-around, …) with a
**scalar amplitude** :math:`A_\alpha` (typically an albedo
:math:`\in [0, 1]`).

Most boundary conditions of practical interest are **rank-1**:

* :class:`VacuumBoundaryOperator` -- :math:`R = 0` (the empty sum,
  rank 0; algebraically the trivial case of the decomposition).
* :class:`SpecularBoundaryOperator` -- :math:`R = G_{\text{refl}} \cdot \alpha`
  where :math:`G_{\text{refl}}` is the angular-permutation operator
  that maps ordinate :math:`\Omega_n` to its reflected partner
  :math:`(\Omega_n \cdot \hat{n})` and :math:`\alpha \in [0, 1]` is the
  specular albedo. Equivalent to a
  :meth:`~orpheus.numerics.measure.DiscreteMeasure.pushforward` under
  the reflection map.
* :class:`WhiteBoundaryOperator` -- :math:`R = G_{\text{diff}} \cdot \alpha`
  where :math:`G_{\text{diff}}` is the cosine-weighted angular average
  over the outgoing hemisphere, broadcast isotropically over the
  incoming hemisphere (Lambertian reflection). Rank-1 in *angle* even
  though the geometric operator is an integral, not a permutation.
* :class:`PeriodicBoundaryOperator` -- :math:`R` is a *spatial*
  pushforward (wrap-around to the opposite face) with :math:`\alpha = 1`.
  Rank-1 in space; the angular structure is identity.
* :class:`AlbedoBoundaryOperator` -- :math:`R = I \cdot \alpha` where
  :math:`I` is the angular identity. Rank-1; the geometric operator is
  trivial.

Mixed and partial-current boundaries are **rank-N** sums of the above
primitives:

* :class:`MixedBoundaryOperator` -- a list of
  ``(weight, BoundaryOperator)`` pairs whose ``apply`` is the linear
  combination of the components.
  ``MixedBoundaryOperator([(0.3, SpecularBoundaryOperator), (0.7, WhiteBoundaryOperator)])``
  realises the standard Marshak mixed boundary
  (Bell & Glasstone 1970, §1.5).

The abstract base :class:`BoundaryOperator` is what production solvers
consume today. Issue 9.6 lifted it from a duck-typed ``Protocol`` into
a concrete ``ABC`` inheriting
:class:`~orpheus.numerics.operator.LinearOperatorMixin` and
:class:`~orpheus.numerics.registry.RegistryMixin`: each concrete
subtype self-registers under its ``key=`` class-creation kwarg, gains
the operator-algebra dunders (``+`` / ``-`` / ``*`` / ``@``), and can
be looked up by string key via :meth:`BoundaryOperator.create`. The
``apply`` method (renamed from the directional ``apply_to_incoming``
it replaced) is the canonical entry point; sweeps and operators call
``bc.apply(outgoing, quadrature)`` with no string-kind branching at
the call site.

Specular reflection ships a load-bearing :meth:`apply_transpose` --
the index-permutation it realises is its own inverse for clean axis
reflections, so the transpose acts as the same permutation scaled by
``albedo``. The other concrete subtypes ship ``apply`` only; their
transposes can be added when a consumer demands them.

Package layout (Wave 4 split)
=============================

The boundary module was split into a package in Wave 4 of the
trace-law refactor (plan: ``.claude/plans/transient-giggling-cake.md``):

* :mod:`_base` -- :class:`BoundaryOperator` (legacy ABC) +
  :class:`BoundaryTraceLaw` (new method-agnostic ABC, Wave 3).
* :mod:`_errors` -- :class:`BoundaryError` + 8 typed subclasses
  (named-error catalog, Wave 3 / ERR-040..ERR-047).
* :mod:`_source` -- :class:`BoundarySource` Protocol +
  :class:`NoSource` + :class:`ConstantInflowSource` (Wave 3).
* :mod:`vacuum` -- :class:`VacuumBoundaryOperator`.
* :mod:`reflective` -- :class:`SpecularBoundaryOperator`.
* :mod:`white` -- :class:`WhiteBoundaryOperator`.
* :mod:`periodic` -- :class:`PeriodicBoundaryOperator`.
* :mod:`albedo` -- :class:`AlbedoBoundaryOperator`.
* :mod:`mixed` -- :class:`MixedBoundaryOperator` (scheduled for
  removal in Wave 11 -- replaced by Wave-0 ``OperatorSum`` algebra).

This ``__init__.py`` re-exports every public name so existing
``from orpheus.geometry.boundary import X`` import sites keep working
unchanged.

References
----------

* Lewis, E. E. & Miller, W. F. *Computational Methods of Neutron
  Transport*, §3.4 (boundary conditions in transport).
* Bell, G. I. & Glasstone, S. *Nuclear Reactor Theory*, §1.5
  (albedo, white, and Marshak boundary conditions).
"""

from __future__ import annotations

# ---------------------------------------------------------------------------
# Abstract bases.
# ---------------------------------------------------------------------------

from ._base import BoundaryOperator, BoundaryTraceLaw

# ---------------------------------------------------------------------------
# Wave 3 -- typed error catalog (ERR-040..ERR-047).
# ---------------------------------------------------------------------------

from ._errors import (
    BoundaryError,
    BoundaryGeometryMapNotMeasurePreservingError,
    BoundaryResponseNotPositiveError,
    BoundarySourceNotOnIncomingTraceError,
    IncomingOutgoingTraceClassificationError,
    ReflectionDidNotMapInflowToOutflowError,
    ReflectionNotInvolutiveError,
    SubmarkovViolationError,
    VacuumAppliedToOutgoingTraceError,
)

# ---------------------------------------------------------------------------
# Wave 3 -- prescribed-inflow source.
# ---------------------------------------------------------------------------

from ._source import BoundarySource, ConstantInflowSource, NoSource

# ---------------------------------------------------------------------------
# Wave 5 -- BoundaryRealizer Protocol + per-method registry.
# ---------------------------------------------------------------------------

from ._realizer import (
    BoundaryRealizer,
    BoundaryRealizerRegistry,
    BoundaryRealizerRegistryError,
)

# ---------------------------------------------------------------------------
# Legacy concrete BCs -- split into per-BC submodules in Wave 4.
# ---------------------------------------------------------------------------

from .albedo import AlbedoBoundaryOperator
from .mixed import MixedBoundaryOperator
from .periodic import PeriodicBoundaryOperator
from .reflective import SpecularBoundaryOperator
from .vacuum import VacuumBoundaryOperator
from .white import WhiteBoundaryOperator


__all__ = [
    # Abstract bases
    "BoundaryOperator",
    "BoundaryTraceLaw",
    # Errors
    "BoundaryError",
    "BoundaryGeometryMapNotMeasurePreservingError",
    "BoundaryResponseNotPositiveError",
    "BoundarySourceNotOnIncomingTraceError",
    "IncomingOutgoingTraceClassificationError",
    "ReflectionDidNotMapInflowToOutflowError",
    "ReflectionNotInvolutiveError",
    "SubmarkovViolationError",
    "VacuumAppliedToOutgoingTraceError",
    # Source
    "BoundarySource",
    "ConstantInflowSource",
    "NoSource",
    # Wave 5 -- realizer Protocol + registry
    "BoundaryRealizer",
    "BoundaryRealizerRegistry",
    "BoundaryRealizerRegistryError",
    # Legacy concretes
    "AlbedoBoundaryOperator",
    "MixedBoundaryOperator",
    "PeriodicBoundaryOperator",
    "SpecularBoundaryOperator",
    "VacuumBoundaryOperator",
    "WhiteBoundaryOperator",
]
