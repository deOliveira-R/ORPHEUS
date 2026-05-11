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

* :class:`VacuumInflow` (deprecated alias: ``VacuumBoundaryOperator``)
  -- :math:`R = 0`, :math:`q = 0` (the empty sum, rank 0;
  algebraically the trivial case of the decomposition).
* :class:`ReflectiveBoundary` (deprecated alias:
  ``SpecularBoundaryOperator``) -- :math:`R = G_{\text{refl}} \cdot
  \alpha` where :math:`G_{\text{refl}}` is the angular-permutation
  operator that maps ordinate :math:`\Omega_n` to its reflected
  partner :math:`(\Omega_n \cdot \hat{n})` and :math:`\alpha \in [0,
  1]` is the specular albedo. Equivalent to a
  :meth:`~orpheus.numerics.measure.DiscreteMeasure.pushforward` under
  the reflection map.
* :class:`WhiteBoundary` (deprecated alias: ``WhiteBoundaryOperator``)
  -- :math:`R = G_{\text{diff}} \cdot \alpha` where
  :math:`G_{\text{diff}}` is the cosine-weighted angular average over
  the outgoing hemisphere, broadcast isotropically over the incoming
  hemisphere (Lambertian reflection). Rank-1 in *angle* even though
  the geometric operator is an integral, not a permutation.
* :class:`PeriodicBoundary` (deprecated alias:
  ``PeriodicBoundaryOperator``) -- :math:`R` is a *spatial*
  pushforward (wrap-around to the opposite face) with
  :math:`\alpha = 1`. Rank-1 in space; the angular structure is
  identity.
* :class:`AlbedoBoundary` (deprecated alias:
  ``AlbedoBoundaryOperator``) -- :math:`R = I \cdot \alpha` where
  :math:`I` is the angular identity. Rank-1; the geometric operator
  is trivial.

Mixed and partial-current boundaries are **rank-N** sums of the above
primitives. Wave 11 removed the specialised ``MixedBoundaryOperator``
composer in favour of the Wave-0 operator algebra: the standard
Marshak mixed boundary is now expressed by realising each leaf
:class:`BoundaryTraceLaw` against the method space and composing the
realised primitives with the algebra dunders

.. code-block:: python

    spec_realized = realizer.realize(ReflectiveBoundary(axis="x"), ms)
    white_realized = realizer.realize(WhiteBoundary(axis="x", outward_sign=+1), ms)
    bc = 0.3 * spec_realized + 0.7 * white_realized

producing a Wave-0 :class:`~orpheus.numerics.operator.OperatorSum` of
:class:`~orpheus.numerics.operator.ScaledOperator` wrappers around the
realised leaves. Tree-walking realisation of a
``BoundaryTraceLaw``-rooted expression is provided by
:func:`orpheus.sn.boundary_realize.realize_recursively` (Wave 11 /
C11.2). See Bell & Glasstone 1970 §1.5 for the physics.

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
* :mod:`vacuum` -- :class:`VacuumInflow` (Wave-7 rename;
  ``VacuumBoundaryOperator`` retained as deprecated alias).
* :mod:`reflective` -- :class:`ReflectiveBoundary` (Wave-7 rename;
  ``SpecularBoundaryOperator`` retained as deprecated alias).
* :mod:`white` -- :class:`WhiteBoundary` (Wave-7 rename;
  ``WhiteBoundaryOperator`` retained as deprecated alias).
* :mod:`periodic` -- :class:`PeriodicBoundary` (Wave-7 rename;
  ``PeriodicBoundaryOperator`` retained as deprecated alias).
* :mod:`albedo` -- :class:`AlbedoBoundary` (Wave-7 rename;
  ``AlbedoBoundaryOperator`` retained as deprecated alias).
* :mod:`prescribed_inflow` -- :class:`PrescribedInflow` (Wave 7
  addition; the rank-0 source-only affine BC).

The ``mixed`` submodule was removed in Wave 11: the rank-N
composition is now Wave-0 ``OperatorSum``-algebra over realised
leaves. See :func:`orpheus.sn.boundary_realize.realize_recursively`
for tree-walking realisation of a ``BoundaryTraceLaw``-rooted
expression.

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
# Abstract base — Wave 7 merged the legacy ``BoundaryOperator`` ABC into
# :class:`BoundaryTraceLaw`. The legacy symbol is kept as a deprecated
# alias for backward compatibility with the 8 production import sites
# (``orpheus.sn.geometry``, ``orpheus.sn.operator``, ...) and the tests
# that import it. Remove in a future cleanup wave.
# ---------------------------------------------------------------------------

from ._base import BoundaryTraceLaw

# Deprecated alias for backward compat; remove in a future cleanup wave.
BoundaryOperator = BoundaryTraceLaw

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
# Concrete BCs -- split into per-BC submodules in Wave 4, renamed to
# Grand Report v3 vocabulary in Wave 7. The legacy
# ``*BoundaryOperator`` names are kept as deprecated aliases below for
# backward compatibility with the production import sites and tests
# that ``from orpheus.geometry.boundary import VacuumBoundaryOperator``.
# Remove the aliases in a future cleanup wave.
# ---------------------------------------------------------------------------

from .albedo import AlbedoBoundary
from .periodic import PeriodicBoundary
from .prescribed_inflow import PrescribedInflow
from .reflective import ReflectiveBoundary
from .vacuum import VacuumInflow
from .white import WhiteBoundary

# Deprecated aliases for backward compat (Wave 7 → cleanup wave).
AlbedoBoundaryOperator = AlbedoBoundary
PeriodicBoundaryOperator = PeriodicBoundary
SpecularBoundaryOperator = ReflectiveBoundary
VacuumBoundaryOperator = VacuumInflow
WhiteBoundaryOperator = WhiteBoundary


__all__ = [
    # Abstract bases
    "BoundaryOperator",  # deprecated alias for BoundaryTraceLaw
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
    # Concrete BCs (Wave 7 canonical names)
    "AlbedoBoundary",
    "PeriodicBoundary",
    "PrescribedInflow",
    "ReflectiveBoundary",
    "VacuumInflow",
    "WhiteBoundary",
    # Deprecated aliases (Wave 7 → cleanup wave)
    "AlbedoBoundaryOperator",
    "PeriodicBoundaryOperator",
    "SpecularBoundaryOperator",
    "VacuumBoundaryOperator",
    "WhiteBoundaryOperator",
]
