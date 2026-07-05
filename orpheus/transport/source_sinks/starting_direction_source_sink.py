r"""The q½ starting-direction source block — the *source* role leaf of
:class:`~orpheus.transport.fields._bases.StartingDirectionField`.

The typed per-level source/sink block driving the starting-direction
equation (Hébert §3.9.4, Eq. 3.432: closed :math:`\mu = \mu_{\rm start}`
streaming :math:`+\;\sigma\psi_{1/2} = q_{1/2}`) — the third block of an
augmented SOURCE composite (a ``q_ext``, an operator ``.apply`` output)
on seed-carrying meshes (#282 route (a); presence per R12a).

Cells vs corner semantics (units note)
======================================

* The **cells** blocks ``(ng, nx)`` carry the within-group emission
  :math:`\bar q_{1/2}` folded to the starting direction (ruling R14:
  the full :math:`(-1)^\ell` Legendre fold) — a volumetric rate
  density, :data:`~orpheus.numerics.units.ANGULAR_RATE_UNITS`
  (``1/(cm³·s·sr)``), exactly like the bulk
  :class:`~orpheus.transport.source_sinks.angular_source_sink.AngularSourceSink`.
* The **corner** slots ``(ng,)`` are trace-like: the inflow corner
  carries the prescribed r = R inflow *datum* (an angular-flux value —
  the affine-BC ``q`` of the seed's identity row, mirroring the
  boundary trace's "on the trace a source does NOT pick up the
  volumetric cm⁻³" convention), and the outflow corner is the defect
  row's source slot. The declared :attr:`UNITS` is the dominant
  cells-block identity; the corner deviation is this documented note
  (the same per-slot exception the angular trace documents). UNITS is
  metadata — the arithmetic gate is class identity.

All storage, validation, algebra, ``cells``/``corner`` views, and the
``zeros_on`` / ``from_mesh`` factories are inherited from
:class:`StartingDirectionField`. Plain vector-space algebra (source
sums are closed) — no role mixin, like every SourceSink leaf.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, ClassVar

import numpy as np

from orpheus.numerics.units import ANGULAR_RATE_UNITS, Unit
from orpheus.transport.fields._bases import StartingDirectionField

if TYPE_CHECKING:  # pragma: no cover
    from numpy.typing import NDArray

    from orpheus.sn.mesh.augmented_mesh import SNMesh

__all__ = ["StartingDirectionSourceSink"]


@dataclass(frozen=True, eq=False, kw_only=True, repr=False)
class StartingDirectionSourceSink(StartingDirectionField):
    r"""L2 starting-direction source/sink — the q½ block of an augmented source composite.

    Parameters
    ----------
    values : NDArray
        Flat backing buffer, shape ``(space.shape[0],)``.
    space : StartingDirectionSpace
        The R12a-keyed space (canonically
        ``mesh.starting_direction_space``).
    mesh : SNMesh
        The SN phase-space carrier (the cross-mesh-arithmetic guard).

    Notes
    -----
    A thin role leaf — the class identity is what keeps seed source,
    flux, and displacement arithmetic from silently mixing (all three
    share the SAME ``mesh.starting_direction_space``, so it is the
    class gate, not the space gate, that rejects cross-role sums).
    """

    #: Dimensional identity (View-G): the cells blocks carry the folded
    #: volumetric emission rate ``1/(cm³·s·sr)`` —
    #: :data:`~orpheus.numerics.units.ANGULAR_RATE_UNITS`, shared with
    #: ``AngularSourceSink``. The corner slots deviate (trace-like flux
    #: datum) — see the module docstring's units note. Metadata, not the
    #: arithmetic gate.
    UNITS: ClassVar[Unit] = ANGULAR_RATE_UNITS

    @classmethod
    def from_angular_source(
        cls, angular_source_values: "NDArray", mesh: "SNMesh",
    ) -> "StartingDirectionSourceSink | None":
        r"""Fold a per-ordinate volumetric source into its q½ seed block.

        The ONE source-side birth factory of #282 route (a) (Pattern 2 —
        the solver cold-starts, the fixed-source rhs, and the operator-
        free :func:`~orpheus.sn.loss_representation.transport_sweep` all
        route through here): ``None`` on a non-carrying mesh; on a
        carrying mesh (1-D curvilinear, R12a) the starting-direction legs
        receive the value of the source at the starting direction
        :math:`\mu = \pm 1`, reconstructed from ALL its Legendre moments
        (Hébert Eq. 3.432, the R14 full :math:`(-1)^\ell` fold):

        .. math::

           \bar q_{1/2}(\mu = \pm 1)
             \;=\; \sum_\ell \tfrac{2\ell+1}{2}\,q_\ell\,(\pm 1)^\ell,
           \qquad
           q_\ell(r) \;=\; \sum_n w_n\,P_\ell(\mu_n)\,q_n(r),

        via the R14 helper
        (:func:`~orpheus.numerics.spaces.starting_direction_space.fold_moments_to_starting_direction`).
        The full fold is REQUIRED for an anisotropic source: even an
        isotropic trial flux :math:`\psi = A(r)` streams to a
        :math:`\mu`-linear source :math:`q = \mu A'(r) + \sigma_t A(r)`,
        whose value at :math:`\mu = -1` is :math:`\sigma_t A - A'` — the
        :math:`\ell = 1` term carries the :math:`-A'` that an
        :math:`\ell = 0`-only fold drops (which floored the anisotropic
        curvilinear MMS; #282 route (a)).  For an isotropic source the
        higher moments vanish and the fold collapses to
        :math:`\tfrac12 q_0` bit-exactly (so the isotropic eigenvalue /
        fixed-source paths are unchanged).  The corner slots stay zero:
        the inflow-corner datum is the BOUNDARY's job (vacuum ⇒ 0;
        reflective ⇒ the ``B`` corner arm into the SI rhs).

        Parameters
        ----------
        angular_source_values : NDArray
            The per-ordinate source in principled 1-D ``(N, ng, nx)``
            layout (carrying meshes are 1-D curvilinear).
        mesh : SNMesh
            The phase-space carrier (its ``starting_direction_space``
            is the R12a presence predicate + the flat layout; its
            ``pole_angular_closure.level_indices`` give each level's
            ordinate bundle for the per-level moment integration).
        """
        from orpheus.numerics.spaces.starting_direction_space import (
            fold_moments_to_starting_direction,
        )

        space = mesh.starting_direction_space
        if space is None:
            return None
        vals = np.asarray(angular_source_values)
        if vals.ndim != 3:
            raise ValueError(
                "StartingDirectionSourceSink.from_angular_source expects the "
                f"principled 1-D (N, ng, nx) per-ordinate layout; got shape "
                f"{vals.shape} (carrying meshes are 1-D curvilinear, R12a)."
            )
        mu = mesh.quad.mu_x
        weights = mesh.quad.weights
        level_indices = mesh.pole_angular_closure.level_indices
        seed = cls.zeros_on(mesh)
        for p in space.levels:
            ords = np.asarray(level_indices[p])
            mu_p = mu[ords]
            w_p = weights[ords]
            q_p = vals[ords]                                  # (M_p, ng, nx)
            # Legendre moments of the source over the level's μ-nodes:
            # q_ℓ = Σ_n w_n P_ℓ(μ_n) q_n, ℓ = 0 … M_p−1 (the full angular
            # content the level resolves; the fold reconstructs q(μ=±1)).
            legendre = np.polynomial.legendre.legvander(mu_p, ords.size - 1)
            moments = np.einsum("n,nl,ngx->lgx", w_p, legendre, q_p)
            for sign in (-1, +1):
                space.cells_view(seed.values, p, sign)[...] = (
                    fold_moments_to_starting_direction(moments, sign)
                )
        return seed
