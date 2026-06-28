r"""The model-independent isotropic energy operators on the scalar flux.

Every transport model — discrete ordinates (SN), collision probability (CP),
method of characteristics (MoC), diffusion, the homogeneous 0-D solver, Monte
Carlo — builds the **isotropic** (:math:`\ell=0`) in-scatter source by the SAME
per-cell energy operation on the **scalar flux** :math:`\phi`:

.. math::

    Q^{\rm iso}_g(\vec r) \;=\; \sum_{g'} \Sigma_{s,0}(g'\to g)\,\phi_{g'}(\vec r)
        \;+\; 2\,\sum_{g'} \Sigma_{2n}(g'\to g)\,\phi_{g'}(\vec r) .

What differs per model is ONLY the *angular wrapper* around it — how :math:`\phi`
is obtained from the model's flux (SN: :math:`\phi=\sum_n w_n\psi_n`; diffusion /
homogeneous: :math:`\phi` is the native unknown; CP / MoC: per region / per FSR)
and how the scalar source is embedded back (SN: :math:`/W` isotropic broadcast;
diffusion / homogeneous: fold :math:`-\Sigma_{s,0}^{T}` into the loss matrix
:math:`A`; CP / MoC: through the transport kernel). The **energy operation is
model-independent** (cross-domain-attacker + explorer, 2026-06-28: the same
``einsum("fg,f->g")`` is re-implemented 6× across the solver families — a
Cardinal-Rule-2 single-concept-many-places situation).

This module owns that shared energy operation as two
:class:`~orpheus.numerics.operator.LinearOperator`\ s on the scalar flux:

* :class:`IsotropicScattering` — :math:`\Sigma_{s,0}` (P0 in-scatter);
* :class:`IsotropicN2N` — :math:`2\,\Sigma_{2n}` (the (n,2n) doubling, a DISTINCT
  *multiplication* channel — physics-faithful, and it also feeds the keff
  *production* numerator, not the loss — so it stays its own operator, summed with
  :class:`IsotropicScattering` for the in-scatter source).

Both are the **scalar (:math:`\ell=0`) realization** of the moment-space operators
:class:`~orpheus.transport.operators.scattering.LegendreMomentScattering` (at
:math:`\ell=0`) and
:class:`~orpheus.transport.operators.scattering.N2NMomentOperator`: they route
through the SAME per-material :class:`~orpheus.transport.mesh.material_xs_field.MaterialXSField`
verbs (``apply_p0_in_scatter`` / ``apply_n2n`` + the ``…_transpose`` twins), so the
cross-section DATA and the per-material dispatch live ONCE (``coding-elegance``
Pattern 2). The harmonic-moment
:attr:`~orpheus.transport.operators.scattering.ScatteringOperator.full_scatter_kernel`
is the permanent verification oracle for this scalar form.

Layout & carriers
=================

Bare ``np.ndarray`` in / out — the model-portable contract (CP / MoC / diffusion
feed raw scalar-flux arrays; SN passes ``phi.values`` and re-wraps). The action is
**spatial-moment-axis-agnostic** (the trailing ``…`` rides an LD ``2^d`` spectator
axis through, #240 D5b-S3), so a φ̂-carrying scalar flux scatters every spatial
moment by the same :math:`\Sigma_s`.

Capabilities
============

``{CAP_APPLY, CAP_APPLY_TRANSPOSE}``. No ``solve`` — the per-cell group-transfer
matrix is generally singular (a thermal group with no up-scatter source has a zero
row); it is *applied*, never inverted. The transpose IS advertised (campaign #276):
:math:`\Sigma_{s,0}^{T}` / :math:`2\Sigma_{2n}^{T}` (the group-axis flip) is the
group-asymmetric factor of the adjoint isotropic scattering source.

References
----------

* explorer (2026-06-28): cross-model grounding — the 6× / 5× duplication inventory.
* cross-domain-attacker (2026-06-28):
  ``iso_source_frame_conjugation_unification.md`` — relocate the shared energy
  operator; do NOT mint a ``ConstantBasis`` iso-frame (it forks ``R∘M``).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING

import numpy as np

from orpheus.numerics.operator import (
    BlockRole,
    CAP_APPLY,
    CAP_APPLY_TRANSPOSE,
    LinearOperator,
)

if TYPE_CHECKING:
    from orpheus.numerics.space import FunctionSpace
    from orpheus.transport.mesh.material_xs_field import MaterialXSField


__all__ = ["IsotropicScattering", "IsotropicN2N"]


def _values_of(phi: "np.ndarray | object") -> np.ndarray:
    """Read the bare ``(ng, *spatial[, 2^d])`` array off a flux carrier or ndarray."""
    return np.asarray(getattr(phi, "values", phi))


@dataclass(frozen=True)
class IsotropicScattering(LinearOperator):
    r"""The P0 isotropic in-scatter energy operator :math:`\Sigma_{s,0}` on the scalar flux.

    Per cell of material :math:`m`, :math:`(\Sigma_{s,0}^{T}\phi)_g =
    \sum_{g'}\Sigma_{s,0}^m(g'\to g)\,\phi_{g'}` — the per-material group-transfer
    matmul, vectorised over cells. The model-independent half of the isotropic
    in-scatter source (the other half is :class:`IsotropicN2N`); both route through
    :class:`~orpheus.transport.mesh.material_xs_field.MaterialXSField` so the
    per-material dispatch lives once.

    Parameters
    ----------
    mat_xs : MaterialXSField
        The macroscopic XS field — the single source of the per-material
        :math:`\Sigma_{s,0}` matrices and the cell-to-material map.
    space : FunctionSpace, optional
        Optional scalar-flux :class:`~orpheus.numerics.space.FunctionSpace` for the
        :class:`~orpheus.numerics.operator.OperatorSum` composition guard; ``None``
        (the default) leaves the operator space-anonymous (the guard skips it — the
        model-portable bare-ndarray contract).
    """

    mat_xs: "MaterialXSField"
    space: "FunctionSpace | None" = None
    capabilities: frozenset[str] = field(
        default_factory=lambda: frozenset({CAP_APPLY, CAP_APPLY_TRANSPOSE})
    )
    # A BULK energy operator (the scalar flux is the bulk block); no boundary action.
    # Class-level constant (unannotated so the dataclass does not treat it as a field).
    block_role = BlockRole.BULK

    def apply(self, phi: "np.ndarray | object") -> np.ndarray:
        r""":math:`\Sigma_{s,0}^{T}\phi` — the per-cell P0 in-scatter source (bare ndarray)."""
        arr = _values_of(phi)
        out = np.zeros_like(arr)
        self.mat_xs.apply_p0_in_scatter(out, arr)
        return out

    def apply_transpose(self, chi: "np.ndarray | object") -> np.ndarray:
        r""":math:`\Sigma_{s,0}\chi` — the group-flip transpose (the bare Euclidean :math:`A^{T}`)."""
        arr = _values_of(chi)
        out = np.zeros_like(arr)
        self.mat_xs.apply_p0_in_scatter_transpose(out, arr)
        return out

    def dense_per_material(self) -> dict[int, np.ndarray]:
        r"""Per-material operator matrix :math:`\Sigma_{s,0}^{T}` (``[g_to, g_from]``).

        The ``as_dense`` consumption mode for the LHS-fold solvers (diffusion /
        homogeneous build :math:`A = \mathrm{diag}(\Sigma_t) - \Sigma_{s,0}^{T} -
        2\Sigma_{2n}^{T}`): returns ``{mid: M}`` with ``M @ φ_cell == apply(φ)_cell``,
        i.e. ``M = sig_s0.T`` (``sig_s0`` is stored ``[g_from, g_to]``). Each entry
        is a fresh copy.
        """
        return {
            mid: np.ascontiguousarray(self.mat_xs.sig_s_legendre(mid)[0].T)
            for mid in self.mat_xs.materials
        }

    @property
    def domain(self) -> "FunctionSpace | None":
        return self.space

    @property
    def codomain(self) -> "FunctionSpace | None":
        return self.space


@dataclass(frozen=True)
class IsotropicN2N(LinearOperator):
    r"""The (n,2n) isotropic energy operator :math:`2\,\Sigma_{2n}` on the scalar flux.

    Per cell, :math:`(2\Sigma_{2n}^{T}\phi)_g = 2\sum_{g'}\Sigma_{2n}^m(g'\to g)\,
    \phi_{g'}`. A DISTINCT *multiplication* channel (each event emits two neutrons),
    NOT scattering — kept its own operator (physics-faithful; it also feeds the keff
    *production* numerator) and summed with :class:`IsotropicScattering` for the
    isotropic in-scatter source. Routes through
    :meth:`~orpheus.transport.mesh.material_xs_field.MaterialXSField.apply_n2n`
    (Pattern 2).

    Parameters
    ----------
    mat_xs : MaterialXSField
        The macroscopic XS field (the per-material :math:`\Sigma_{2n}` matrices).
    space : FunctionSpace, optional
        See :class:`IsotropicScattering`.
    """

    mat_xs: "MaterialXSField"
    space: "FunctionSpace | None" = None
    capabilities: frozenset[str] = field(
        default_factory=lambda: frozenset({CAP_APPLY, CAP_APPLY_TRANSPOSE})
    )
    block_role = BlockRole.BULK

    def apply(self, phi: "np.ndarray | object") -> np.ndarray:
        r""":math:`2\Sigma_{2n}^{T}\phi` — the per-cell (n,2n) source (bare ndarray)."""
        arr = _values_of(phi)
        out = np.zeros_like(arr)
        self.mat_xs.apply_n2n(out, arr)
        return out

    def apply_transpose(self, chi: "np.ndarray | object") -> np.ndarray:
        r""":math:`2\Sigma_{2n}\chi` — the group-flip transpose."""
        arr = _values_of(chi)
        out = np.zeros_like(arr)
        self.mat_xs.apply_n2n_transpose(out, arr)
        return out

    def dense_per_material(self) -> dict[int, np.ndarray]:
        r"""Per-material operator matrix :math:`2\Sigma_{2n}^{T}` — ``M @ φ == apply(φ)``."""
        return {
            mid: np.ascontiguousarray(2.0 * self.mat_xs.n2n_matrix(mid).T)
            for mid in self.mat_xs.materials
        }

    @property
    def domain(self) -> "FunctionSpace | None":
        return self.space

    @property
    def codomain(self) -> "FunctionSpace | None":
        return self.space
