"""Angular quadrature for the Method of Characteristics.

Product quadrature: uniform azimuthal angles x Tabuchi-Yamamoto (TY)
polar angles.  The TY quadrature is optimised for the Bickley-function
integrals that arise when the 2-D MOC flat-source solution is integrated
over polar angle (Yamamoto et al., J. Nucl. Sci. Technol. 44(2), 2007).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from orpheus.numerics.roots_of_unity import roots_of_unity

# ── Tabuchi-Yamamoto tables ─────────────────────────────────────────
# Each entry: (sin_theta, weight) for one polar angle per half-space.
#
# Source, VERIFIED 2026-08-11 against the paper (now in the library, with an
# OCR sidecar): Yamamoto, Tabuchi, Sugimura, Ushio & Mori, "Derivation of
# Optimum Polar Angle Quadrature Set for the Method of Characteristics Based
# on Approximation Error for the Bickley Function", J. Nucl. Sci. Technol.
# 44(2), 129-136 (2007), doi:10.3327/jnst.44.129 — **Table 1**, "Tabuchi and
# Yamamoto (TY)" column, printed p. 133 (PDF p. 6).
#   ⛔ This comment said "Table 2" until 2026-08-11. Table 2 is a C5G7
#      k-effective comparison and carries no quadrature. The AUTHOR
#      attribution was right; only the table number was wrong.
#
# CONVENTION: the published weights sum to 1.000000 (the mass of
# `sin_theta d_theta` on (0, pi/2)), and the values below are those weights
# HALVED, so they sum to 0.5 per hemisphere and 1.0 over the full sphere.
# The published `sin theta` values are stored verbatim, unhalved.
#
# HOW THE PAPER DERIVED THEM (`[M]` printed p. 133, and §II): minimize the
# MAXIMUM absolute error of Ki_3 -- the paper's `E_{3,max}` -- by steepest
# descent from an initial guess; its Figs. 2-4 plot that residual over
# x in [0, 5], one decade smaller per added division. Our own residual against
# `orpheus.derivations.common.kernels.ki_n` reproduces exactly that: `[M]`
# max|r| = 1.279e-02 (L=1), 4.123e-04 (L=2), 3.361e-05 (L=3), and at L=3 it
# EQUIOSCILLATES with six alternations at ~3.35e-05 -- the Chebyshev signature
# of a minimax fit. ⟹ these 12 literals are DERIVABLE, and #355 is the issue
# to replace them with the derivation plus its own machine-checkable claim.
# (The 1- and 3-division Leonard columns of the same Table 1 were re-derived
# by these authors against `E_{2,max}`, i.e. Ki_2 -- a second target family.)
#
# See also Knott & Yamamoto, "Lattice Physics Computations", Handbook of
# Nuclear Engineering ch. 9 (2010), doi:10.1007/978-0-387-98149-9_9, which
# re-tabulates the set; it is what OpenMOC cites for its `TYPolarQuad`.

_TY_TABLES: dict[int, tuple[np.ndarray, np.ndarray]] = {
    1: (
        np.array([0.798184]),
        np.array([0.500000]),
    ),
    2: (
        np.array([0.363900, 0.899900]),
        np.array([0.212854 / 2, 0.787146 / 2]),
    ),
    3: (
        np.array([0.166648, 0.537707, 0.932954]),
        np.array([0.046233 / 2, 0.283619 / 2, 0.670148 / 2]),
    ),
}


@dataclass(frozen=True)
class MOCQuadrature:
    r"""Product quadrature for 2-D Method of Characteristics.

    Azimuthal: ``n_azi`` uniform midpoint angles in [0, pi) — the grid
    :math:`\varphi_k = \pi(2k+1)/(2n)`, i.e. the upper half of the
    STAGGERED periodic trapezoid on :math:`2n` points.  Each azimuthal
    angle defines a family of parallel tracks; the supplementary angle
    (phi + pi) is the same physical track traversed in the opposite
    direction.

    The PRIMARY azimuthal representation is the point pair
    ``(cos_phi, sin_phi)``, generated from the group action
    (:func:`~orpheus.numerics.roots_of_unity.roots_of_unity` at the
    rational :math:`(2k+1)/(4n)` of a full turn) rather than by
    evaluating trig on stored angles — so the dihedral symmetry is
    EXACT on floats: the azimuthal mirror :math:`\varphi \to \pi -
    \varphi` is the index map :math:`k \mapsto n-1-k` with
    ``cos_phi[n-1-k] == -cos_phi[k]`` and ``sin_phi[n-1-k] ==
    sin_phi[k]`` bit-for-bit, and odd ``n_azi``'s vertical member
    (:math:`\varphi = \pi/2`) carries ``cos_phi`` exactly ``0.0``
    (issue #325).  ``phi`` is the derived angle CHART kept for
    human-facing consumers; geometry consumes the points.

    Polar: Tabuchi-Yamamoto quadrature with ``n_polar`` angles per
    half-space.  The polar angle enters the MOC only through the
    effective 3-D optical thickness  tau / sin(theta_p).

    Attributes
    ----------
    n_azi : int
        Number of azimuthal angles in [0, pi).
    n_polar : int
        Number of polar angles per half-space (1, 2, or 3).
    phi : ndarray, shape (n_azi,)
        Azimuthal angles in [0, pi) (radians) — the chart of the exact
        points, spelled in the rational form ``pi*(2k+1)/(2n)``.
    cos_phi, sin_phi : ndarray, shape (n_azi,)
        The exact azimuthal direction components (the points on the
        half-circle; ``sin_phi > 0`` everywhere — the staggered grid
        never touches the axis).
    omega_azi : ndarray, shape (n_azi,)
        Azimuthal weights, each = 1 / n_azi, summing to 1.
    sin_polar : ndarray, shape (n_polar,)
        sin(theta_p) for each TY polar angle.
    omega_polar : ndarray, shape (n_polar,)
        TY polar weights (sum = 0.5 for one hemisphere).
    """

    n_azi: int
    n_polar: int
    phi: np.ndarray
    cos_phi: np.ndarray
    sin_phi: np.ndarray
    omega_azi: np.ndarray
    sin_polar: np.ndarray
    omega_polar: np.ndarray

    @classmethod
    def create(cls, n_azi: int = 16, n_polar: int = 3) -> MOCQuadrature:
        """Build an MOC product quadrature.

        Parameters
        ----------
        n_azi : int
            Number of azimuthal angles in [0, pi).  Must be >= 2.
        n_polar : int
            Number of TY polar angles per half-space (1, 2, or 3).
        """
        if n_azi < 2:
            raise ValueError(f"n_azi must be >= 2, got {n_azi}")
        if n_polar not in _TY_TABLES:
            raise ValueError(
                f"n_polar must be one of {sorted(_TY_TABLES)}, got {n_polar}"
            )

        # The midpoint grid phi_k = pi*(2k+1)/(2n) as EXACT points: the
        # rational (2k+1)/(4n) of a full turn, generated by the group
        # action so mirrors/axis members are bit-exact (#325).  The
        # angle chart is derived in the same rational spelling.
        odd_k = 2 * np.arange(n_azi, dtype=np.int64) + 1
        cos_phi, sin_phi = roots_of_unity(odd_k, 4 * n_azi)
        phi = np.pi * odd_k / (2.0 * n_azi)
        omega_azi = np.full(n_azi, 1.0 / n_azi)
        sin_polar, omega_polar = _TY_TABLES[n_polar]

        return cls(
            n_azi=n_azi,
            n_polar=n_polar,
            phi=phi,
            cos_phi=cos_phi,
            sin_phi=sin_phi,
            omega_azi=omega_azi,
            sin_polar=sin_polar.copy(),
            omega_polar=omega_polar.copy(),
        )
