r"""Back-compat shim: re-exports the real-SH evaluator from its new home.

The canonical home of the real-:math:`Y_\ell^m` evaluator is now
:class:`orpheus.numerics.basis.SphericalHarmonicBasis` (Wave 0 P1.1 of
the moment-space + layering plan). This module ships
:func:`evaluate_real_sh` as a back-compat alias so the per-quadrature
``spherical_harmonics(L)`` methods on
:class:`orpheus.numerics.quadrature.Quadrature` and the projection
factory :meth:`orpheus.numerics.projection.MomentProjection.from_measure`
continue to import the same function name they always have.

This shim deletes in P3.2 when those delegators are rewired to consume
:meth:`SphericalHarmonicBasis.evaluate` directly.

See :mod:`orpheus.numerics.basis.spherical_harmonic_basis` for the
convention, the addition-theorem label
:eq:`real-sh-addition-theorem`, the Gram-matrix identity
:eq:`sh-mass-matrix-diagonal`, and the cross-method consumer list.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

from orpheus.numerics.basis.spherical_harmonic_basis import (
    SphericalHarmonicBasis,
    _evaluate_real_sh,
)


__all__ = ["evaluate_real_sh", "SphericalHarmonicBasis"]


def evaluate_real_sh(
    L: int,
    mu_x: NDArray,
    mu_y: NDArray,
    mu_z: NDArray,
) -> np.ndarray:
    r"""Evaluate real spherical harmonics on a discrete direction set.

    Back-compat alias for the algorithm body of
    :meth:`SphericalHarmonicBasis.evaluate` /
    :meth:`SphericalHarmonicBasis.evaluate_from_components`. Preserves
    the legacy ``(L, mu_x, mu_y, mu_z)`` per-component signature.

    Parameters
    ----------
    L : int
        Maximum spherical-harmonic order. ``L == 0`` returns the
        :math:`P_0` table only; ``L < 0`` returns an all-zero array of
        zero shape.
    mu_x, mu_y, mu_z : NDArray, shape ``(N,)``
        Direction cosines of the ordinates :math:`\hat\Omega_n =
        (\mu_x, \mu_y, \mu_z)`. The polar axis is :math:`\mu_x` (so
        :math:`\mu_x = \cos\theta`); azimuth is measured in the
        :math:`(\mu_y, \mu_z)` plane.

    Returns
    -------
    np.ndarray, shape ``(N, L+1, 2L+1)``
        ``Y[n, l, l+m]`` is :math:`Y_\ell^m(\hat\Omega_n)` under the
        no-:math:`4\pi/(2\ell+1)`-prefactor convention. See
        :mod:`orpheus.numerics.basis.spherical_harmonic_basis` for the
        full convention.

    Notes
    -----
    Implementation chains through ``_evaluate_real_sh`` in the basis
    module — bit-identical to the legacy free-function shipped before
    P1.1.
    """
    return _evaluate_real_sh(L, mu_x, mu_y, mu_z)
