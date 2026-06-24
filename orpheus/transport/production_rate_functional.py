r"""The §5.6 production-rate functional :math:`\phi \mapsto p(\vec r)`.

The grand-report §5.6 *suffix law* (see
:mod:`orpheus.numerics.functional`) names three kinds of map by what
they return; this module ships the canonical transport-level
**Functional**: the per-cell fission emission density

.. math::

    p(\vec r) \;=\; \sum_{g'} \nu\Sigma_{f,g'}(\vec r)\,\phi_{g'}(\vec r) ,

the rate at which fission emits new neutrons per unit volume at
:math:`\vec r`, summed over the source groups :math:`g'`. It consumes a
flux field :math:`\phi` (group-resolved, ``(ng, *spatial)``) and returns
a **scalar-field** (the group axis collapsed, one number per cell): a
Functional applied *fiberwise over space*, in the §5.6 sense.

Where this sits in fission
==========================

Frame 3 of the coefficient-promotion analysis (the cross-domain-attacker
memo ``coefficient_field_promotion_frames.md``) decomposes the fission
operator as a composition

.. math::

    F \;=\; M_\chi \;\circ\; \mathrm{ProductionRate} \;\circ\; M_{\nu\Sigma_f} ,

where :math:`M_{\nu\Sigma_f}` multiplies the flux by the production cross
section, :math:`\mathrm{ProductionRate}` (this functional) contracts the
group axis to the per-cell scalar density :math:`p(\vec r)`, and
:math:`M_\chi` re-broadcasts that density across groups weighted by the
emission spectrum :math:`\chi`. This class is the **middle** of that
chain. In #257 S5 the functional is created and cross-checked against the
existing fission realisation, ready for S6 to wire the composition; the
:class:`~orpheus.numerics.operator.RankOneOperator`-backed fission kernel
is **not** rewired here (that principled, non-bit-identical step is S6).

Bit-identity with the legacy ``inner`` contraction
===================================================

The group contraction :math:`\sum_{g'} \nu\Sigma_{f,g'}\,\phi_{g'}` is
*exactly* the anonymous ``inner`` reduction inside
:meth:`orpheus.numerics.operator.RankOneOperator.apply`:

.. code-block:: python

    inner = (self.right * x_arr).sum(axis=self.axis, keepdims=True)

with ``right = νΣf`` and ``axis = 0`` (the group axis). :meth:`evaluate`
reproduces that reduction **byte-for-byte** — the same numpy primitive,
the same axis, the same ``keepdims=True`` — so the S6 composition will be
bit-identical to the legacy fission rank-1 by construction. Naming the
contraction turns the most physically central diagnostic in criticality
(the per-cell fission rate) from an anonymous einsum into a typed,
inspectable Functional (``coding-elegance`` Pattern 3).

The density carries NO volume measure
=====================================

:meth:`evaluate` returns the per-cell *density* :math:`p(\vec r)`, **not**
the volume integral :math:`\int p\,dV`. No cell volume is folded in, and
no :math:`\chi` is re-broadcast (that is the :math:`M_\chi` factor, a
separate S6 concern). The volume-weighted production *rate* (the integral)
is a distinct functional that lives in the eigenvalue drivers
(:func:`orpheus.numerics.iteration._default_production_estimator`); this
leaf stays the pure density so the suffix law reads cleanly:
"the production rate is :math:`\nu\Sigma_f` contracted against the flux
over groups."

References
----------

* Grand Report v3 §5.5–5.7 — the field hierarchy, the suffix law, and
  the fission composition :math:`F = M_\chi \circ \mathrm{ProductionRate}
  \circ M_{\nu\Sigma_f}`.
* ``.claude/agent-memory/cross-domain-attacker/coefficient_field_promotion_frames.md``
  — Frame 3 (the production rate as a Functional; the locality criterion
  separating this from the scattering / fission integral kernels).
* ``.claude/plans/issue_257_coefficient_field_promotion.md`` — S5.
* Lewis, E.E. & Miller, W.F. (1993). *Computational Methods of Neutron
  Transport*. ANS. §1.1 — the fission source.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from orpheus.numerics.field import Field
from orpheus.transport.fields.cross_section_field import CrossSectionField

__all__ = ["ProductionRateFunctional"]


@dataclass(frozen=True)
class ProductionRateFunctional:
    r"""The §5.6 functional :math:`\phi \mapsto \sum_{g'} \nu\Sigma_{f,g'}\,\phi_{g'}`.

    A :class:`~orpheus.numerics.functional.Functional`: it exposes only
    :meth:`evaluate` (field → scalar-field) and deliberately none of the
    :class:`~orpheus.numerics.operator.LinearOperator` surface (no
    ``apply``, no ``capabilities``). That disjointness makes "the
    production rate is a Functional, not an operator" a structural fact.

    Parameters
    ----------
    nu_sigma_f : CrossSectionField
        The production cross section :math:`\nu\Sigma_f` (units ``1/cm``),
        a per-group per-cell coefficient field. Obtain it from the typed
        accessor
        :meth:`orpheus.transport.mesh.material_xs_field.MaterialXSField.fission_production_field`.
        Its ``.values`` is the ``(ng, *spatial)`` array contracted against
        the flux over the leading (group) axis.
    """

    nu_sigma_f: CrossSectionField

    def evaluate(self, phi: Field | NDArray, /) -> np.ndarray:
        r"""Return the per-cell fission emission density :math:`p(\vec r)`.

        .. math::

            p(\vec r) \;=\; \sum_{g'} \nu\Sigma_{f,g'}(\vec r)\,
                \phi_{g'}(\vec r) .

        The group axis (axis 0) is contracted; ``keepdims=True`` leaves a
        length-1 leading axis so the output ``(1, *spatial)`` is
        byte-identical to the ``inner`` reduction inside
        :meth:`orpheus.numerics.operator.RankOneOperator.apply`. No cell
        volume measure is folded in (this is the density, not the
        integral) and no :math:`\chi` is re-broadcast (that is the S6
        :math:`M_\chi` composition).

        Parameters
        ----------
        phi : Field | NDArray
            The group-resolved flux carrier, ``(ng, *spatial)``. A typed
            flux :class:`~orpheus.numerics.field.Field` (carrying
            ``.values``) or a raw ``np.ndarray`` is accepted; the
            contraction reads the underlying array.

        Returns
        -------
        np.ndarray
            The per-cell production density, shape ``(1, *spatial)``.
        """
        phi_values = phi.values if isinstance(phi, Field) else np.asarray(phi)
        return (self.nu_sigma_f.values * phi_values).sum(axis=0, keepdims=True)
