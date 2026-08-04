r"""Zero-flux (homogeneous Dirichlet) boundary law.

See :class:`ZeroFluxBoundary` for the algebraic definition. Minted in
issue #290 P3 (2026-07-03, user ruling 3) to give the textbook
Dirichlet idealization :math:`\phi_\Gamma = 0` its own honestly-named
law — the legacy diffusion island registered exactly this condition
under the name ``"vacuum"``, which it is NOT: vacuum means *zero
incoming current* (:math:`J^- = 0`, the Marshak realization of
:class:`~orpheus.geometry.boundary.vacuum.VacuumInflow`), while
zero-flux pins the scalar flux itself to zero AT the surface.
"""

from __future__ import annotations

from dataclasses import dataclass

from ._base import BoundaryTraceLaw
from ._factors import SelfPairedDeck, ScalarResponse


__all__ = ["ZeroFluxBoundary"]


@dataclass(frozen=True)
class ZeroFluxBoundary(BoundaryTraceLaw, key="zero_flux"):
    r"""Zero-flux Dirichlet boundary: :math:`\phi_\Gamma = 0`.

    The **mathematical idealization** of a bare surface — the classic
    textbook condition under which the bare-slab criticality problem
    has the closed-form sine mode vanishing exactly at the physical
    surface (no extrapolation length). It is an excellent TEST /
    reference condition (every analytic sine anchor in
    ``tests/diffusion`` is derived under it) but it is **not a
    physical transport boundary**: no material or void outside the
    domain produces it.

    In the partial-current basis (the scalar trace's ``(J⁺, J⁻)``
    pair, ``.claude/plans/diffusion_crosswalk.md``), the P1 closure
    dictionary :math:`\phi_\Gamma = 2(J^+ + J^-)` makes the law a
    member of the albedo family :math:`J^- = \mathcal{A}\,J^+` with

    .. math::

        \phi_\Gamma = 0
        \;\Longleftrightarrow\;
        J^- = -\,J^+
        \;\Longleftrightarrow\;
        \mathcal{A} = -1 .

    **Deliberately outside the sub-Markov family.** The physical
    albedo laws satisfy :math:`\mathcal{A} \in [0, 1]` and preserve
    partial-current positivity (:math:`J^\pm \ge 0`); zero-flux
    requires a NEGATIVE incoming partial current, which no
    non-negative angular flux can produce. Positivity is therefore a
    property of the *physical* laws, not a type invariant of the
    family — this class does not override
    :meth:`~BoundaryTraceLaw.assert_response_positive_if_declared`
    to raise, because the negative response IS its definition.

    Realization
    -----------
    * **Diffusion** — :math:`\mathcal{A} = -1`, i.e.
      ``ScaledOperator(-1, IdentityOperator)`` on the outflow trace
      (:class:`~orpheus.diffusion.boundary_realizer.DiffusionBoundaryRealizer`).
    * **SN** — REFUSED
      (:class:`~orpheus.sn.boundary.realizer.SNBoundaryRealizer`
      raises :class:`~orpheus.geometry.boundary.BoundaryError`): a
      transport boundary condition prescribes the inflow trace from
      the outflow trace, and a negative angular inflow is not
      representable in a non-negative angular flux. Use
      :class:`~orpheus.geometry.boundary.vacuum.VacuumInflow` for the
      physical zero-incoming law.

    Factors, and the tier caveat
    ----------------------------
    The two factors below describe the **realization**, which genuinely is
    affine: :math:`G = I`, :math:`R = -1` on the partial-current trace, exactly
    what the diffusion realizer builds. The **posing** is a tier above —
    :math:`\phi_\Gamma = 0` is the relation :math:`A_-\gamma_- + A_+\gamma_+ =
    0`, not a map from outflow to inflow, which is why SN cannot realize it and
    must refuse by hand rather than by type. Grand Report v3 names that tier
    ``BoundaryRelation``; issue **#177** builds it. Until then these factors are
    honest about the realization without claiming the law is affine.

    Naming ruling (#290, 2026-07-03)
    --------------------------------
    The legacy ``orpheus.diffusion`` island's ``BC_REGISTRY`` mapped
    the string ``"vacuum"`` to this Dirichlet condition
    (``_diff_bc_vacuum``: "Zero-flux (Dirichlet φ=0)") — stale,
    unfaithful naming. Vacuum MEANS :math:`J^- = 0`. The analytic
    references derived under :math:`\phi_\Gamma = 0` survive
    mathematically intact and are re-attributed to THIS law; they are
    never renamed away.

    References
    ----------
    * Lamarsh, *Introduction to Nuclear Reactor Theory* (1966) §9
      (bare-slab one-group criticality under :math:`\phi(\pm a)=0`).
    * Bell & Glasstone, *Nuclear Reactor Theory* (1970) §3.4
      (diffusion boundary conditions; the Marshak vacuum condition
      this idealization is often conflated with).
    """

    # No parameters: the law is stateless (like VacuumInflow). The
    # frozen dataclass gives value equality (all instances equal) and
    # registry construction via ``BoundaryTraceLaw.create("zero_flux")``.

    # ── The affine form's two factors (B1) ──────────────────────────────
    # These describe the REALIZATION (see "Factors, and the tier caveat"
    # above), not a claim that the law is posed affinely.
    @property
    def geometry_map(self) -> "SelfPairedDeck":
        r""":math:`G = I` — no relabeling; the sign lives entirely in `R`."""
        return SelfPairedDeck.identity()

    @property
    def response_kernel(self) -> "ScalarResponse":
        r""":math:`R = -1` — deliberately outside the sub-Markov :math:`[0,1]`.

        :math:`J^- = -J^+` is what makes :math:`\phi_\Gamma = 0`, and it is why
        SN refuses this law: a negative angular inflow is unrepresentable in a
        non-negative :math:`\psi`. The value is the diffusion realizer's
        ``ScaledOperator(-1.0, IdentityOperator)`` scalar, named.
        """
        return ScalarResponse(-1.0)
