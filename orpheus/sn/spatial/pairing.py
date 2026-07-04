r"""Pairing-validity predicates for the (spatial ⊗ angular) discretization product.

Issue #236 realizes the SN discretization as a tensor product of two
independently-selectable axes — a SPATIAL closure
(:class:`~orpheus.transport.spatial.scheme.DiscretizationScheme`) and an ANGULAR
redistribution closure
(:class:`~orpheus.sn.spatial.pole_angular_closure.PoleAngularClosureBase`).
Some properties of the discretization are properties of the *pair*, not of
either axis alone.  This module is the home for those pairing-validity
predicates: each reads only the class-level traits the two axes declare, so a
frontend can query a pairing's validity without constructing a solve.

The diffusion-limit predicate
=============================

Whether a discretization recovers the thick-diffusion limit FACTORIZES into a
SPATIAL condition and an ANGULAR condition — the two were established by two
separate papers, the strongest evidence they live on independent axes:

* SPATIAL (Larsen–Morel–Miller 1987, *J. Comput. Phys.* **69**): the scheme's
  thick limit must be a consistent diffusion discretization for the
  leading-order scalar flux.  Carried by
  ``DiscretizationScheme.diffusion_limit_consistent`` (Diamond Difference and
  full Linear-Discontinuous ``True``; Step ``False``).
* ANGULAR (Bailey–Morel–Chang 2010, *Nucl. Sci. Eng.* **165**): the
  redistribution closure's first-order functional :math:`\beta` (BMC Eq. (41))
  must vanish.  Carried by ``PoleAngularClosureBase.beta_first_order_consistent``
  (Morel–Montry ``True`` by BMC Eq. (42); the Cartesian identity closure
  vacuously ``True``).

Both must hold simultaneously — independence of *selection* is not independence
of *consequence*: a good Morel–Montry closure paired with bare Step still breaks
the limit (spatial ``False``), and full LD paired with a step-in-angle closure
on a curvilinear mesh still dips (angular ``False``).

The Cartesian collapse is **automatic, not a branch here**: the identity closure
declares ``beta_first_order_consistent = True`` (no redistribution ⇒
:math:`\alpha \equiv 0` ⇒ :math:`\beta \equiv 0`), so on a Cartesian mesh the
conjunction reduces to the spatial condition alone — the collapse lives in the
closure's trait, where it belongs, not in a geometry test in this predicate.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from orpheus.transport.spatial.scheme import DiscretizationScheme

    from .pole_angular_closure import PoleAngularClosureBase


def pair_diffusion_limit_consistent(
    scheme: "DiscretizationScheme",
    closure: "PoleAngularClosureBase",
) -> bool:
    r"""Whether the (spatial scheme, angular closure) PAIR recovers the thick-diffusion limit.

    The conjunction of the two per-axis conditions:

    .. math::

        \text{pair valid} \;=\;
        \underbrace{\texttt{scheme.diffusion\_limit\_consistent}}
            _{\text{Larsen–Morel–Miller 1987, spatial}}
        \;\wedge\;
        \underbrace{\texttt{closure.beta\_first\_order\_consistent}}
            _{\text{Bailey–Morel–Chang 2010, angular}}.

    See the module docstring for the factorization, the joint requirement, and
    the Cartesian collapse (encoded by the identity closure's
    ``beta_first_order_consistent = True``, not by a branch here).

    Both arguments are read for a class-level trait only, so either the class or
    an instance may be passed; production calls pass the instances
    ``mesh.scheme`` and ``mesh.pole_angular_closure``.

    Parameters
    ----------
    scheme :
        The spatial discretization scheme (e.g. ``DiamondDifference()``,
        ``LinearDiscontinuous()``).  Read for ``diffusion_limit_consistent``.
    closure :
        The angular pole-redistribution closure (``MorelMontryAngularSweep`` for
        curvilinear, ``IdentityAngularClosure`` for Cartesian).  Read for
        ``beta_first_order_consistent``.

    Returns
    -------
    bool
        ``True`` iff BOTH the spatial scheme and the angular closure are
        diffusion-limit-consistent.
    """
    return (
        scheme.diffusion_limit_consistent
        and closure.beta_first_order_consistent
    )


__all__ = ["pair_diffusion_limit_consistent"]
