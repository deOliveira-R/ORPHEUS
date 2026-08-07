r"""Pure albedo boundary condition.

See :class:`AlbedoBoundary` for the algebraic definition. This class
was previously named ``AlbedoBoundaryOperator``; the legacy alias was
retired in Wave O step O.4a.1. ``AlbedoBoundary`` (the Grand Report
v3 vocabulary) is the sole live name.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Optional

from ._base import BoundaryTraceLaw
from ._factors import (
    BoundaryResponseKernel,
    SelfPairedDeck,
    ReemissionClosure,
    ScalarResponse,
    SpecularReemission,
)
from ._specular import assert_specular_pairing_valid

if TYPE_CHECKING:
    import numpy as np

    from orpheus.numerics.quadrature import Quadrature


__all__ = ["AlbedoBoundary"]


@dataclass(frozen=True)
class AlbedoBoundary(BoundaryTraceLaw, key="albedo"):
    r"""An absorbing **surface**: a fraction :math:`\alpha` of the arriving
    flux is returned, in an angular shape the caller chooses.

    .. math::

        \gamma_-\psi \;=\; \alpha\,C\,\gamma_+\psi,

    where :math:`C` is the **re-emission closure** — the pairing that says
    which outgoing direction feeds which incoming one. The law's two factors
    are :math:`G = \mathrm{id}` (a surface is not a quotient of the domain, so
    it fixes no geometry) and :math:`R = \alpha C` (all of it constitutive).

    The two degrees of freedom are independent
    ------------------------------------------

    *How much* comes back is :attr:`albedo`; *in what shape* is
    :attr:`reemission`. Every shape admits every amplitude, so they are
    separate parameters rather than a menu of their product:

    .. code-block:: python

        AlbedoBoundary(0.7, SpecularReturn(axis="x"))          # polished wall
        AlbedoBoundary(0.7, IsotropicReturn(axis="x", outward_sign=+1))  # matte
        AlbedoBoundary(0.7)                                    # shape UNSTATED

    Why the shape may be left unstated
    -----------------------------------

    ``reemission=None`` is not an omission — it is the honest spelling of
    :math:`R = \alpha\,I`, which is **complete on a scalar trace** and
    **under-determined on an angular one**:

    * **Diffusion** resolves one boundary degree of freedom, :math:`J^- =
      \alpha J^+`. There is no angular distribution to fix, so the bare
      amplitude IS the whole law and the realizer reads
      :attr:`~orpheus.geometry.boundary.ScalarResponse.amplitude` and stops.
      ``BC("albedo", albedo=…)`` means exactly what it always meant.
    * **SN** resolves the full hemisphere. There :math:`\alpha\,I` is a map
      :math:`\Gamma_+ \to \Gamma_+`, and :math:`G = \mathrm{id}` supplies no
      crossing to :math:`\Gamma_-`. Composing it anyway pairs incoming
      ordinate :math:`j` with outgoing ordinate :math:`j` — an artefact of
      **array position**, not geometry. The SN realizer therefore refuses it
      and names the two completions rather than picking one.

    The positional pairing is not merely *meaningless* — it is a
    **configuration-dependent accident**, which is worse. Measured on the
    ``xmax`` face, sorted index sets:

    .. code-block:: text

        product(2, 4)        positional pairing == the specular one   True
        level_symmetric(6)   positional pairing == the specular one   True
        gauss_legendre(4/8)  positional pairing == the specular one   False
        lebedev(17)          positional pairing == the specular one   False

    So before this closure existed, a bare albedo silently behaved like a
    mirror on two of the tree's quadratures and like nothing in particular on
    the others — with no signal at the seam. A user who validated on
    ``level_symmetric`` and ran on ``lebedev`` got different physics from the
    same law. That is the strongest argument for the refusal: the old answer
    was not a defensible default, it was a coincidence of index order.

    This is the method-realizability taxonomy's **angular-resolution** axis
    (``bc-method-realizability``) biting for the first time. It is a refusal of
    an *incomplete spelling*, not of the law: both completions are fully built
    and realized, and they route through the same bodies that realize
    :class:`~orpheus.geometry.boundary.ReflectiveBoundary` and
    :class:`~orpheus.geometry.boundary.WhiteBoundary`.

    Equivalences, by construction
    -----------------------------

    ``AlbedoBoundary(α, SpecularReturn(a)) ≡ ReflectiveBoundary(a, α)`` and
    ``AlbedoBoundary(α, IsotropicReturn(a, s)) ≡ WhiteBoundary(a, s, α)`` as
    *matrices* — they execute the same realizer body — while asserting
    different things: a surface's constitutive return versus a symmetry of the
    domain. See :class:`~orpheus.geometry.boundary.SpecularReemission` for why
    that distinction is load-bearing rather than pedantic.

    This is a **pure descriptor** (Issue #186 / B3 + β2) — it carries no
    ``apply`` method. Realise via
    :class:`~orpheus.sn.boundary.realizer.SNBoundaryRealizer` (angular) or
    :class:`~orpheus.diffusion.boundary_realizer.DiffusionBoundaryRealizer`
    (scalar).

    Rename history
    --------------
    Previously named ``AlbedoBoundaryOperator``. The legacy alias was
    retired in Wave O step O.4a.1 — ``AlbedoBoundary`` is the sole
    importable name.

    Parameters
    ----------
    albedo : float
        Fraction returned. ``0`` is vacuum, ``1`` returns everything.
    reemission : ReemissionClosure, optional
        The angular shape of the return —
        :class:`~orpheus.geometry.boundary.SpecularReturn` or
        :class:`~orpheus.geometry.boundary.IsotropicReturn`. ``None`` (the
        default) leaves it unstated, which is complete for a scalar method and
        refused by an angular one.
    """

    albedo: float = 0.0
    reemission: "Optional[ReemissionClosure]" = None

    # ── The affine form's two factors (B1) ──────────────────────────────
    @property
    def geometry_map(self) -> "SelfPairedDeck":
        r""":math:`G = I`, **unconditionally** — a surface is not a quotient.

        Whatever the closure, the deck transformation stays the group's
        identity: the domain does not continue on the other side of a wall, so
        nothing is identified with anything. In particular
        ``reemission=SpecularReturn(axis)`` does **not** put a mirror here —
        the specular pairing it asks for is a statement about the *wall*, and
        it lives in :attr:`response_kernel`. That is the user's 2026-08-01
        ruling and the discipline it enforces: precisely what is Geometric,
        precisely what is Response.

        It is also why the diffusion arm needs no geometric factor at all — on
        a scalar trace ``G`` is forced to the identity by dimension anyway.
        """
        return SelfPairedDeck.identity()

    @property
    def response_kernel(self) -> "BoundaryResponseKernel":
        r""":math:`R = \alpha\,C` — the whole law.

        With no closure this is :class:`~orpheus.geometry.boundary.ScalarResponse`
        — :math:`\alpha` alone, the SAME number the diffusion realizer's first
        stage returns (``_partial_current_albedo(law) -> float``). With one, it
        is that closure instantiated at this law's :attr:`albedo`, so
        :math:`\alpha` still has exactly one home.
        """
        if self.reemission is None:
            return ScalarResponse(self.albedo)
        return self.reemission.kernel(self.albedo)

    # ------------------------------------------------------------------
    # §16A.12 universal invariants — Wave 7 / C7.6 overrides.
    # ------------------------------------------------------------------

    def assert_response_positive_if_declared(self) -> None:
        r"""Albedo coefficient must be non-negative.

        Raises
        ------
        BoundaryResponseNotPositiveError
            When ``self.albedo < 0``.
        """
        if self.albedo < 0.0:
            from ._errors import BoundaryResponseNotPositiveError
            raise BoundaryResponseNotPositiveError(
                f"Albedo BC albedo={self.albedo} < 0",
                law="albedo",
            )

    def assert_submarkov(self) -> None:
        r"""Albedo satisfies the sub-Markov bound :math:`\alpha \le 1`.

        Raises
        ------
        SubmarkovViolationError
            When ``self.albedo > 1``.
        """
        if self.albedo > 1.0:
            from ._errors import SubmarkovViolationError
            raise SubmarkovViolationError(
                f"Albedo BC albedo={self.albedo} > 1",
                law="albedo",
            )

    def assert_realizable(
        self,
        quadrature: "Quadrature",
        *,
        inflow_indices: "Optional[np.ndarray]" = None,
    ) -> None:
        r"""Universal invariants + the sub-Markov bound (ERR-046) +, when the
        closure is specular, the pairing's three table checks.

        ``assert_submarkov`` is albedo-specific (not one of the
        §16A.12 universal five), so it joins the certification here —
        the same extension pattern as :class:`ReflectiveBoundary`'s
        table checks.

        **The specular clause is B3.4b's, and it closes a hole the phase would
        otherwise have opened.** ``AlbedoBoundary(α, SpecularReturn(axis))``
        realizes through the mirror MOTION the closure names
        (:attr:`~orpheus.geometry.boundary.SpecularReemission.motion` →
        the deck kernel, G6.3 step 7 — the same element
        :class:`ReflectiveBoundary`'s :math:`G` carries, read from one
        spelling), whose ERR-042/044/045 invariants were that law's methods.
        Adding a specular closure without adding this
        call would mean a wrong pairing is **caught on one route and silently
        realized on the other**, which is precisely the twin-path failure this
        campaign exists to remove. The certification now belongs to the pairing
        (:mod:`~orpheus.geometry.boundary._specular`) and both carriers fire it.

        Note the tier: this law's :math:`G` **is** the identity and **does**
        preserve the measure, so the base's
        :meth:`~orpheus.geometry.boundary.BoundaryTraceLaw.assert_geometry_map_measure_preserving`
        hook is correctly a no-op here. The invariant is fired as a
        law-specific extension because the pairing sits in :math:`R`. Reading
        that hook as "the pairing's check" is the conflation the user's ruling
        forbids.
        """
        super().assert_realizable(quadrature, inflow_indices=inflow_indices)
        self.assert_submarkov()
        kernel = self.response_kernel
        if isinstance(kernel, SpecularReemission):
            assert_specular_pairing_valid(
                quadrature, kernel.axis, law_key="albedo",
            )
