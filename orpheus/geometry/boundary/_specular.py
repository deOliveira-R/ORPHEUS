r"""Certification of a **specular pairing** — three independent invariants.

A specular pairing is the ordinate permutation realizing
:math:`\Omega \mapsto \Omega - 2(\Omega\cdot\hat n)\hat n` on a quadrature —
since **G6.3 step 7d** derived HERE from the mirror MOTION
(:func:`_specular_pairing` →
:meth:`~orpheus.numerics.quadrature.Quadrature.ordinate_permutation`), the
same single source the realization's deck kernel reads. Until 7d these
checks read the precomputed ``reflection_index`` table while realization
derived the permutation from the law's motion — two objects for one
pairing, bit-identical by gate but structurally twinned; certifying the
very array the kernel narrows is what closes the split.
Two laws carry one, in **different tiers** of the affine form:

* :class:`~orpheus.geometry.boundary.ReflectiveBoundary` carries it as
  :math:`G` (:meth:`~orpheus.geometry.boundary.SelfPairedDeck.mirror`) — a
  symmetry of the domain;
* :class:`~orpheus.geometry.boundary.AlbedoBoundary` with a
  :class:`~orpheus.geometry.boundary.SpecularReturn` closure carries it as
  :math:`R` (:class:`~orpheus.geometry.boundary.SpecularReemission`) — a
  polished surface's constitutive return.

They realize to the same matrix and assert different physics (see
:class:`~orpheus.geometry.boundary.SpecularReemission`), but the *pairing*
they both stand on is one object with one set of invariants. This module is
that one set.

Why the checks moved here (campaign phase B3.4b)
================================================

They were three methods on :class:`~orpheus.geometry.boundary.ReflectiveBoundary`
— which was correct while reflective was the only law with a specular pairing.
The moment the user's 2026-08-01 ruling put a specular closure on
:class:`~orpheus.geometry.boundary.AlbedoBoundary`, keeping them there would
have meant **one route certified and the other not**: a wrong
``reflection_index`` table caught on ``ReflectiveBoundary(axis)`` and silently
realized on ``AlbedoBoundary(α, SpecularReturn(axis))``. The invariants belong
to the pairing, not to the law that historically owned it.

The three are INDEPENDENT — this is the catalog's lesson
=========================================================

ERR-045's entry states it verbatim: *"the inflow partition, the involution
property, and the inflow → outflow image are three independent invariants. All
three must hold; checking only one or two leaves a hole."* Concretely, each
check passes on a table the other two reject:

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - broken table
     - who catches it, and who does NOT
   * - pairs ordinates from different weight classes
     - **measure** (ERR-042) catches it; the involution does not (it can still
       be its own inverse) and the sign check does not (partners can still be
       opposite-signed)
   * - maps an inflow ordinate to **itself**
     - **inflow→outflow** (ERR-045) catches it; the involution does not (a
       self-map is trivially involutive) and the measure does not (a self-map
       trivially preserves :math:`m_n`)
   * - a non-involutive cycle of length 3
     - **involution** (ERR-044) catches it; the other two can both pass

So the aggregate :func:`assert_specular_pairing_valid` fires all three, and a
caller that needs only one (the polymorphic
``assert_geometry_map_measure_preserving`` hook the base template fires) reaches
the individual function rather than a re-implementation.

.. note::

   The base template's hook is named
   :meth:`~orpheus.geometry.boundary.BoundaryTraceLaw.assert_geometry_map_measure_preserving`,
   which reads the invariant as belonging to :math:`G`. Since B3.4b a pairing in
   :math:`R` needs the same check, so the hook's name is narrower than the
   concept. It is not *wrong* — it does check exactly what it says — and
   renaming a polymorphic hook with live overrides is deferred to **B4**, where
   the realizer starts dispatching on the factors and the tier question is
   answered structurally rather than by a method name.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from orpheus.numerics.face_layout import AXIS_NAMES
from orpheus.numerics.spaces.angular_trace_space import TANGENTIAL_EPS

from ._factors import _mirror_motion

if TYPE_CHECKING:
    from orpheus.numerics.quadrature import Quadrature


__all__ = [
    "assert_specular_pairing_involutive",
    "assert_specular_pairing_maps_inflow_to_outflow",
    "assert_specular_pairing_measure_preserving",
    "assert_specular_pairing_valid",
]


def _specular_pairing(
    quadrature: "Quadrature", axis: str, *, law_key: str,
) -> "np.ndarray":
    r"""The pairing under certification — derived from the mirror MOTION,
    the same single source realization reads (G6.3 step 7d).

    Each invariant below asserts ITS OWN claim on this array — the
    derivation's built-in certificate (bijection + bare-weight equality,
    from :meth:`~orpheus.geometry.transformation.RigidMotion.preserves`)
    is deliberately NOT treated as subsuming them: ERR-044's involution
    has no carrier in ``preserves`` (a bijection needn't be its own
    inverse), and ERR-042's measure is the TRACE one,
    :math:`w\,|\mu_a|`, not the bare :math:`w` the matcher compares.
    The checks keep their teeth; only the array's source moved.

    Raises
    ------
    BoundaryError
        When the rule is not closed under the mirror — no bijective,
        weight-preserving match of the ordinates onto their mirror
        images exists (e.g. an odd-:math:`n_\varphi` product rule has
        no x-mirror closure). This is the refusal that used to surface
        as ``reflection_index``'s missing-axis ``ValueError``; a pairing
        that does not exist can satisfy no invariant, so it is refused
        before any of the three is asked.
    """
    pi = quadrature.ordinate_permutation(_mirror_motion(axis, 3))
    if pi is None:
        from ._errors import BoundaryError

        raise BoundaryError(
            f"no specular pairing about axis {axis!r} exists on this "
            f"quadrature: the mirror does not permute its weighted "
            f"ordinate set (no bijective, weight-preserving match — "
            f"e.g. an odd-n_phi product rule has no x-mirror closure), "
            f"so a specular law cannot be certified, let alone "
            f"realized.",
            law=law_key,
        )
    return pi.indices


#: Relative tolerance for the ERR-042 measure-preservation check.
#: Reflection partners are constructed from the same symmetric node/weight
#: sets, so a correct table agrees to a few ULP; a wrong table (mispaired
#: weight classes) violates at O(1). Twelve orders of headroom on both sides.
_MEASURE_RTOL: float = 1e-12


def _axis_cosines(quadrature: "Quadrature", axis: str) -> "np.ndarray":
    r"""The ordinates' orbit-mean direction cosines along ``axis``.

    ``AXIS_NAMES.index`` rather than a local ``{"x": 0, ...}`` literal: the
    name→column correspondence has ONE home
    (:data:`~orpheus.numerics.face_layout.AXIS_NAMES`), and a second spelling
    of it is a twin waiting to disagree.

    ⭐ ``mean_axis_cosine``, not ``axis_cosines`` (phase 0.2): this feeds the
    cosine measure :math:`w\lvert\mu_a\rvert` on a face, which is a **flux**
    question — so on an axis a 1-D rule has suppressed, zero is the answer
    (nothing flows along it), not a missing one. `[M]` the full-suite census
    records this site reading a suppressed axis 5 times, so the distinction is
    live here rather than theoretical.
    """
    return quadrature.mean_axis_cosine(AXIS_NAMES.index(axis))


def assert_specular_pairing_measure_preserving(
    quadrature: "Quadrature", axis: str, *, law_key: str,
) -> None:
    r"""The pairing preserves the direction-cosine measure
    :math:`w(\Omega)\,|\Omega\cdot\hat n|` (ERR-042).

    Specular reflection is an isometry, so the pushforward of the discrete face
    measure under the index map :math:`\pi` must be the measure itself: with
    :math:`m_n = w_n\,|\mu_{a,n}|` on axis :math:`a`,

    .. math::

        m_{\pi(n)} \;=\; m_n \qquad \forall\, n .

    Violation means the reflected current at the face differs from the incident
    current — a persistent phantom net leakage at a purely reflective boundary,
    amplified by heterogeneous :math:`\Sigma_t` contrast.

    Checked DIRECTLY, independent of the involution property: an involutive
    table that pairs ordinates from different weight classes passes the
    involution check while breaking the measure — exactly the hole the pre-#52
    delegation (*"weight equality is implied by construction"*) left open.

    Raises
    ------
    BoundaryGeometryMapNotMeasurePreservingError
        When any ordinate and its partner carry different direction-cosine
        measures (relative tolerance :data:`_MEASURE_RTOL`; tangential
        ordinates carry zero measure, so a tangential↔non-tangential
        mispairing is caught by the same comparison).
    """
    perm = _specular_pairing(quadrature, axis, law_key=law_key)
    cosine_measure = quadrature.weights * np.abs(_axis_cosines(quadrature, axis))
    partner_measure = cosine_measure[perm]
    if not np.allclose(
        partner_measure, cosine_measure, rtol=_MEASURE_RTOL, atol=0.0
    ):
        from ._errors import BoundaryGeometryMapNotMeasurePreservingError

        worst = int(np.argmax(np.abs(partner_measure - cosine_measure)))
        raise BoundaryGeometryMapNotMeasurePreservingError(
            f"the specular pairing about axis {axis!r} does not preserve "
            f"the direction-cosine measure w·|μ_{axis}|: ordinate "
            f"{worst} carries m={cosine_measure[worst]:.6e} but its "
            f"partner {int(perm[worst])} carries "
            f"m={partner_measure[worst]:.6e} (ERR-042: wrong "
            f"pairing, or a quadrature whose weights "
            f"are inconsistent with its nodes).",
            law=law_key,
        )


def assert_specular_pairing_involutive(
    quadrature: "Quadrature", axis: str, *, law_key: str,
) -> None:
    r"""The pairing is an involution, :math:`\pi \circ \pi = \mathrm{id}`
    (ERR-044).

    An axis reflection applied twice returns every ordinate to itself. A table
    with a longer cycle is not a reflection at all, and the operator it
    realizes has no honest transpose — the ``.H`` that the mirror deck element
    (:meth:`~orpheus.geometry.boundary.SelfPairedDeck.mirror`) and
    :class:`~orpheus.geometry.boundary.SpecularReemission` both advertise as a
    theorem rests on exactly this property.

    Raises
    ------
    ReflectionNotInvolutiveError
        When the derived specular pairing is not its own inverse.
    """
    ref = _specular_pairing(quadrature, axis, law_key=law_key)
    if not np.array_equal(ref[ref], np.arange(quadrature.N)):
        from ._errors import ReflectionNotInvolutiveError

        raise ReflectionNotInvolutiveError(
            f"the specular pairing about axis {axis!r} is not an involution",
            law=law_key,
        )


def assert_specular_pairing_maps_inflow_to_outflow(
    quadrature: "Quadrature", axis: str, *, law_key: str,
) -> None:
    r"""Every non-tangential ordinate pairs with the OPPOSITE sign class
    (ERR-045).

    The pairing's geometric contract at a face with normal along axis
    :math:`a`: every inflow ordinate has an outflow partner with
    :math:`\mu_{a,\pi(n)} = -\mu_{a,n}` (and vice versa — the sign classes
    swap, so the check is face-side-free). A table that maps an inflow ordinate
    to itself, or to any same-sign partner, plants a length-1 self-loop in the
    sweep dependency graph that degenerates the sweep convergence theory — while
    passing BOTH the involution check (a self-map is trivially involutive) and
    the measure check (a self-map trivially preserves :math:`m_n`).

    Tangential ordinates (:math:`|\mu_a| \le`
    :data:`~orpheus.numerics.spaces.angular_trace_space.TANGENTIAL_EPS`) are
    exempt: they belong to neither half-trace, and a tangential self-map is the
    geometrically correct image.

    Raises
    ------
    ReflectionDidNotMapInflowToOutflowError
        When any non-tangential ordinate's partner is tangential or lies in the
        same sign class.
    """
    perm = _specular_pairing(quadrature, axis, law_key=law_key)
    mu_axis = _axis_cosines(quadrature, axis)
    partner_mu = mu_axis[perm]
    active = np.abs(mu_axis) > TANGENTIAL_EPS
    partner_opposite = (np.abs(partner_mu) > TANGENTIAL_EPS) & (
        np.sign(partner_mu) != np.sign(mu_axis)
    )
    bad = active & ~partner_opposite
    if np.any(bad):
        from ._errors import ReflectionDidNotMapInflowToOutflowError

        worst = int(np.flatnonzero(bad)[0])
        raise ReflectionDidNotMapInflowToOutflowError(
            f"the specular pairing about axis {axis!r} maps ordinate "
            f"{worst} "
            f"(μ_{axis}={mu_axis[worst]:+.6e}) to ordinate "
            f"{int(perm[worst])} "
            f"(μ_{axis}={partner_mu[worst]:+.6e}) — same sign "
            f"class instead of the outflow partner (ERR-045: wrong "
            f"pairing, or a non-axis-aligned "
            f"reflection that needs a different BC type).",
            law=law_key,
        )


def assert_specular_pairing_valid(
    quadrature: "Quadrature", axis: str, *, law_key: str,
) -> None:
    r"""All three invariants — the full certification of a specular pairing.

    Call this from a law's ``assert_realizable`` when the law realizes a
    specular pairing (the deck kernel derives it from the same mirror
    motion these checks certify), whichever tier the pairing sits in.
    Firing only a subset leaves one of the three holes tabulated in this
    module's docstring.
    """
    assert_specular_pairing_measure_preserving(quadrature, axis, law_key=law_key)
    assert_specular_pairing_involutive(quadrature, axis, law_key=law_key)
    assert_specular_pairing_maps_inflow_to_outflow(quadrature, axis, law_key=law_key)
