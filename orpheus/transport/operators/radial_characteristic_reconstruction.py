r"""``A_BA``'s fold — reconstruct a bulk angular source at the closed μ=±1 rays.

This module hosts :class:`RadialCharacteristicReconstruction`, the shared
factor of the ``A_BA`` bulk→ray coupling block of the augmented SN system
(campaign :mod:`coupled block operator <orpheus.sn>`)::

    [ A_AA   A_AB ] [ transport ]      A_BA = Reconstruction ∘ Emission
    [ A_BA   A_BB ] [ ray       ]

**Home — a step-2 TRANSIENT in transport; it moves to sn at step 4.** This
operator's *native* home is beside its self-block sibling
:class:`~orpheus.sn.operators.radial_characteristic.RadialCharacteristicOperator`
(``A_BB``) in ``orpheus/sn/operators/radial_characteristic.py`` — both are
blocks of the one ψ½ coupled operator. It sits in ``transport/`` at step 2
ONLY because the model-generic
:class:`~orpheus.transport.operators.scattering.ScatteringOperator` and
:class:`~orpheus.transport.operators.fission.FissionOperator` still *consume*
it in their ``(ray, bulk)`` seed arms: a runtime ``transport → sn`` import is
forbidden (sn depends on transport, never the reverse), so a fold that S/F
consume must sit at/below transport.

That consumption is the **monolithic posing this campaign un-welds** — the ψ½
starting-direction ray is a *curvilinear-SN augmentation*, and a model-generic
gain has no business owning its projection. **Campaign step 4 reverses the
dependency**: S/F become pure bulk, ``A_BA = Fold ∘ K_iso`` is composed by the
coupled driver (the Wave-O #208 lagged-gain pattern), and this operator moves
to sn beside ``A_BB``. The ψ½ *data* types legitimately stay at the
transport/numerics data layer (this operator's
:class:`~orpheus.transport.source_sinks.RadialCharacteristicSourceSink` codomain,
the carrier space — generic transport ops like ``M[σ]`` and the residual act on
them), so only the *operator* migrates. Until then it needs SN solely for the
carrier-mesh *type* (``TYPE_CHECKING``), and the runtime ``transport → sn`` ban
holds. (User ruling 2026-07-08 — see the campaign plan step-4 lift deliverable.)
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Optional

import numpy as np

from orpheus.numerics.operator import LinearOperator
from orpheus.numerics.spaces.radial_characteristic_space import (
    fold_moments_to_radial_characteristic,
    fold_moments_to_radial_characteristic_transpose,
)

if TYPE_CHECKING:
    from orpheus.numerics.space import FunctionSpace
    from orpheus.numerics.spaces.radial_characteristic_space import (
        RadialCharacteristicSpace,
    )
    from orpheus.sn.mesh.augmented_mesh import SNMesh
    from orpheus.transport.fields._bases import RadialCharacteristicField
    from orpheus.transport.source_sinks import RadialCharacteristicSourceSink


__all__ = ["RadialCharacteristicReconstruction"]


class RadialCharacteristicReconstruction(LinearOperator):
    r"""Reconstruct a bulk angular source at the closed μ=±1 rays — the ``A_BA`` fold.

    The shared factor of the ``A_BA`` bulk→ray coupling block
    (``A_BA = RadialCharacteristicReconstruction ∘ Emission``, assembled in
    campaign step 4): given a bulk within-group angular source in its
    Legendre-moment representation, it evaluates that source at the two
    closed radial-characteristic rays :math:`\mu = \pm 1` and writes the
    result as the q½ ray source on every carried μ-level.

    **What it is (operator algebra).** The 1-D angular reconstruction
    :math:`\mathcal R` (Hébert Eq. 3.432) SAMPLED at the closed rays:

    .. math::

        \bar q_{1/2}(\mu = \pm 1)
          \;=\; \sum_\ell \frac{2\ell + 1}{2}\,q_\ell\,(\pm 1)^\ell ,

    the same :math:`(2\ell+1)/2\,(\pm 1)^\ell` weights the angular frame
    reconstructs with (:math:`P_\ell(\pm 1) = (\pm 1)^\ell`), evaluated at the
    rays rather than at the quadrature nodes. Its Euclidean transpose is the
    injection of a ray cotangent back into moment space
    (:func:`~orpheus.numerics.spaces.radial_characteristic_space.fold_moments_to_radial_characteristic_transpose`)
    — so it advertises ``is_adjointable = True`` (the shape of the shared
    :class:`~orpheus.transport.operators.isotropic_scattering.IsotropicScattering`
    kernel).

    **The single source of the S/F seed fold.** The scattering and fission
    composite operators each emit a q½ seed on a seed-carrying input — the
    ``(ray, bulk)`` block of their lagged gain. That block factors as this
    reconstruction after the operator-specific isotropic emission
    (:math:`K_{\rm iso}\varphi_0` for S, :math:`\chi\nu\Sigma_f\varphi/k`
    for F). Routing both forward arms AND the scattering seed adjoint
    through here makes the fold — and its transpose — single-sourced
    (Cardinal Rule 2), retiring the S-adjoint's hand-rolled :math:`\ell = 0`
    factor :math:`\tfrac12`. The forward reads a Legendre-moment source
    ``(n_moments, ng, nx)``; the production emitters feed the isotropic
    :math:`\ell = 0` emission (``n_moments = 1``), and the fold accepts any
    order so the anisotropic reach is testable before it is needed.

    **Broadcast across levels.** The same moment source is folded onto every
    carried level (an angularly-uniform-across-levels source — exact for the
    isotropic emission). Carrying meshes have EXACTLY one level (R12a), so
    "broadcast across levels" is "on the one ray"; the transpose therefore
    sums the per-level, per-sign cotangents. Corners stay zero: the fold
    writes only the cells legs (the inflow-corner datum is the boundary
    block ``B_b``'s job; scattering/fission are volumetric).

    Parameters
    ----------
    sn_mesh : SNMesh
        The seed-carrying geometry (1-D curvilinear, R12a). Supplies the ray
        carrier ``radial_characteristic_space`` (the codomain) and the
        cells-leg layout. A seedless mesh has no bulk→ray coupling → rejected.
    n_moments : int, default 1
        The operator's domain dimension — the number of angular Legendre
        moments it reconstructs from. ``1`` is the isotropic production reach
        (:math:`\ell = 0`); a larger order is the manufactured anisotropic
        case. Fixed at construction so the transpose has a well-defined
        codomain ``(n_moments, ng, nx)``.
    """

    def __init__(self, sn_mesh: "SNMesh", n_moments: int = 1) -> None:
        space = sn_mesh.radial_characteristic_space
        if space is None:
            raise ValueError(
                "RadialCharacteristicReconstruction: the mesh carries no "
                "radial-characteristic ray (radial_characteristic_space is "
                "None) — a seedless mesh (Cartesian / cylinder, R12a) has no "
                "bulk→ray coupling to fold."
            )
        if n_moments < 1:
            raise ValueError(
                f"RadialCharacteristicReconstruction: n_moments must be ≥ 1 "
                f"(at least the ℓ = 0 moment); got {n_moments!r}."
            )
        #: The augmented geometry (ray carrier + cells-leg layout).
        self.sn_mesh = sn_mesh
        #: The angular Legendre-moment order of the domain (``1`` = isotropic,
        #: the production reach; larger is the manufactured anisotropic case).
        self.n_moments = n_moments
        #: The ψ½ carrier — the codomain (non-None by the ctor guard).
        self._ray_space: "RadialCharacteristicSpace" = space

    # ── Predicates / spaces ───────────────────────────────────────────

    @property
    def is_adjointable(self) -> bool:
        # Both the forward fold and its Euclidean transpose are realized
        # (:meth:`apply` / :meth:`apply_transpose`) — the same shape as the
        # shared IsotropicScattering kernel. is_invertible inherits False (a
        # reconstruction from n_moments to the two rays is not square).
        return True

    @property
    def domain(self) -> Optional["FunctionSpace"]:
        # The bulk angular-moment source — an untyped ``(n_moments, ng, nx)``
        # intermediate (like K_iso's ndarray domain); no minted moment space.
        return None

    @property
    def codomain(self) -> Optional["FunctionSpace"]:
        return self._ray_space

    # ── Forward fold — reconstruct at μ=±1 ────────────────────────────

    def apply(self, moments: np.ndarray, /) -> "RadialCharacteristicSourceSink":
        r"""Reconstruct the bulk moment source at μ=±1 → the q½ ray source.

        Folds the Legendre-moment source ``moments`` (shape
        ``(n_moments, ng, nx)``) onto every carried level's cells legs at both
        closed rays via
        :func:`~orpheus.numerics.spaces.radial_characteristic_space.fold_moments_to_radial_characteristic`
        (the single-source fold). The corners stay zero — the fold is
        volumetric (the inflow-corner datum is ``B_b``'s job).

        Parameters
        ----------
        moments : np.ndarray
            The bulk within-group source in ℓ-leading Legendre moments,
            shape ``(n_moments, ng, nx)``. S and F pass the isotropic
            emission with a unit ℓ axis (``emission[None]``, ``n_moments=1``).

        Returns
        -------
        RadialCharacteristicSourceSink
            The q½ ray source — cells legs folded, corners zero.
        """
        from orpheus.transport.source_sinks import RadialCharacteristicSourceSink

        arr = np.asarray(moments)
        expected = (self.n_moments, self._ray_space.ng, self._ray_space.nx)
        if arr.shape != expected:
            raise ValueError(
                f"RadialCharacteristicReconstruction.apply expects a bulk "
                f"Legendre-moment source of shape (n_moments, ng, nx) = "
                f"{expected}; got {arr.shape}."
            )
        space = self._ray_space
        seed = RadialCharacteristicSourceSink.zeros_on(self.sn_mesh)
        for level in space.levels:
            for sign in (-1, +1):
                space.cells_view(seed.values, level, sign)[...] = (
                    fold_moments_to_radial_characteristic(arr, sign)
                )
        return seed

    # ── Euclidean transpose — inject a ray cotangent into moment space ─

    def apply_transpose(
        self, cotangent: "RadialCharacteristicField", /,
    ) -> np.ndarray:
        r"""Euclidean transpose — a ray cotangent → the bulk moment cotangent.

        The adjoint of :meth:`apply` with respect to the moments: sum the
        per-level, per-sign ray-cells cotangents expanded back onto the
        moments via
        :func:`~orpheus.numerics.spaces.radial_characteristic_space.fold_moments_to_radial_characteristic_transpose`
        (the SAME reconstruction weights the forward contracts — the sign is
        spelled once). This is the single source of the scattering seed
        adjoint (``∂S/∂ψ½`` cotangent → the ℓ = 0 bulk moment, which
        ``K_isoᵀ`` then scatters into the bulk).

        Parameters
        ----------
        cotangent : RadialCharacteristicField
            A cotangent on the ray source (role-erased). Must share this
            operator's mesh. Only the cells legs are read (the forward writes
            only cells); the corner cotangents are ignored.

        Returns
        -------
        np.ndarray
            The bulk moment cotangent, shape ``(n_moments, ng, nx)``.
        """
        if cotangent.mesh is not self.sn_mesh:
            raise ValueError(
                f"RadialCharacteristicReconstruction.apply_transpose: the "
                f"cotangent and the operator must share the same SNMesh "
                f"(mesh-identity invariant); got field mesh "
                f"{cotangent.mesh!r} vs operator mesh {self.sn_mesh!r}."
            )
        space = self._ray_space
        cot_vals = cotangent.values
        moment_bar = np.zeros(
            (self.n_moments, space.ng, space.nx), dtype=float,
        )
        for level in space.levels:
            for sign in (-1, +1):
                cells_bar = space.cells_view(cot_vals, level, sign)
                moment_bar += fold_moments_to_radial_characteristic_transpose(
                    cells_bar, sign, self.n_moments,
                )
        return moment_bar

    def __repr__(self) -> str:
        return (
            f"RadialCharacteristicReconstruction("
            f"n_moments={self.n_moments}, levels={self._ray_space.levels}, "
            f"ng={self._ray_space.ng}, nx={self._ray_space.nx})"
        )
