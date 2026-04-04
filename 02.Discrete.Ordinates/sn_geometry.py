r"""Augmented geometry for S\ :sub:`N` discrete ordinates transport.

:class:`SNMesh` wraps a :class:`~geometry.mesh.Mesh1D` or
:class:`~geometry.mesh.Mesh2D` and precomputes the coordinate-specific
streaming stencil used by the transport sweep.

Currently only Cartesian coordinates are implemented.  Cylindrical and
spherical geometries will add curvature/angular-redistribution terms to
the stencil.
"""

from __future__ import annotations

import numpy as np

from geometry import CoordSystem, Mesh1D, Mesh2D
from sn_quadrature import AngularQuadrature


class SNMesh:
    """Augmented geometry for the discrete ordinates method.

    Wraps a :class:`~geometry.mesh.Mesh1D` or :class:`~geometry.mesh.Mesh2D`
    and precomputes the streaming stencil (diamond-difference coefficients
    that depend only on geometry and angular quadrature, not on cross
    sections).

    For Cartesian geometry the stencil stores:

    * ``streaming_x[n, i]`` = :math:`2|\\mu_{x,n}| / \\Delta x_i`
    * ``streaming_y[n, j]`` = :math:`2|\\mu_{y,n}| / \\Delta y_j`

    For future curvilinear geometries, additional curvature terms
    (:math:`\\alpha_n / r_i`) will be stored in ``self.curvature``.

    Parameters
    ----------
    mesh : Mesh1D or Mesh2D
        Base geometry.
    quadrature : AngularQuadrature
        Angular quadrature (Gauss–Legendre, Lebedev, etc.).
    """

    def __init__(
        self,
        mesh: Mesh1D | Mesh2D,
        quadrature: AngularQuadrature,
    ) -> None:
        self.mesh = mesh
        self.quad = quadrature

        # Normalise to (nx, ny) shaped arrays for both 1-D and 2-D
        if isinstance(mesh, Mesh1D):
            self.nx: int = mesh.N
            self.ny: int = 1
            self.dx: np.ndarray = mesh.widths
            self.dy: np.ndarray = np.array([1.0])
            self.mat_map: np.ndarray = mesh.mat_ids.reshape(mesh.N, 1)
            self._volumes: np.ndarray = mesh.volumes.reshape(mesh.N, 1)
        else:
            self.nx = mesh.nx
            self.ny = mesh.ny
            self.dx = mesh.dx
            self.dy = mesh.dy
            self.mat_map = mesh.mat_map
            self._volumes = mesh.volumes

        # Dispatch stencil setup by coordinate system
        match mesh.coord:
            case CoordSystem.CARTESIAN:
                self._setup_cartesian()
            case CoordSystem.CYLINDRICAL:
                raise NotImplementedError(
                    "Cylindrical SN transport not yet implemented. "
                    "Requires angular redistribution terms."
                )
            case CoordSystem.SPHERICAL:
                self._setup_spherical()

    # ── Properties ────────────────────────────────────────────────────

    @property
    def volumes(self) -> np.ndarray:
        """Cell volumes, shape (nx, ny)."""
        return self._volumes

    @property
    def is_1d(self) -> bool:
        """True if this is a 1-D mesh (ny == 1)."""
        return self.ny == 1

    # ── Stencil setup ─────────────────────────────────────────────────

    def _setup_cartesian(self) -> None:
        """Precompute Cartesian diamond-difference streaming coefficients.

        These are the purely geometric parts of the DD denominator:

        .. math::

            \\text{denom} = \\Sigma_t + \\frac{2|\\mu_x|}{\\Delta x}
                            + \\frac{2|\\mu_y|}{\\Delta y}

        Precomputing avoids per-ordinate per-cell divisions in the
        inner sweep loop.
        """
        mu_x = self.quad.mu_x
        mu_y = self.quad.mu_y

        # streaming_x[n, i] = 2|μ_x[n]| / dx[i] — shape (N_ord, nx)
        self.streaming_x: np.ndarray = (
            2.0 * np.abs(mu_x)[:, None] / self.dx[None, :]
        )
        # streaming_y[n, j] = 2|μ_y[n]| / dy[j] — shape (N_ord, ny)
        self.streaming_y: np.ndarray = (
            2.0 * np.abs(mu_y)[:, None] / self.dy[None, :]
        )

        # Curvature terms (None for Cartesian — placeholder for curvilinear)
        self.curvature = None

    def _setup_spherical(self) -> None:
        r"""Precompute spherical streaming stencil and angular redistribution.

        The 1-D spherical transport equation for ordinate *n*, cell *i* is:

        .. math::

            \frac{\mu_n}{V_i}
            \bigl[A_{i+\frac12}\psi_{i+\frac12}
                - A_{i-\frac12}\psi_{i-\frac12}\bigr]
            + \frac{1}{V_i}
            \bigl[\alpha_{n+\frac12}\psi_{n+\frac12}
                - \alpha_{n-\frac12}\psi_{n-\frac12}\bigr]
            + \Sigma_t \psi_{n,i} = Q_{n,i}

        Precomputed quantities:

        * ``face_areas`` — :math:`A_{i+1/2} = 4\pi r_{i+1/2}^2`
        * ``alpha_half`` — :math:`\alpha_{n+1/2} = \sum_{m=0}^{n} w_m \mu_m`

        The :math:`\alpha` coefficients satisfy :math:`\alpha_{1/2} = 0`
        and :math:`\alpha_{N+1/2} = 0` (by Gauss–Legendre antisymmetry).
        """
        mu = self.quad.mu_x
        w = self.quad.weights
        N = self.quad.N

        # Cell face areas: A_{i+1/2} = 4πr² at each edge
        self.face_areas: np.ndarray = self.mesh.surfaces  # (nx+1,)

        # Angular redistribution coefficients
        # α_{n+1/2} = Σ_{m=0}^{n} w_m μ_m
        alpha = np.zeros(N + 1)
        for n in range(N):
            alpha[n + 1] = alpha[n] + w[n] * mu[n]
        self.alpha_half: np.ndarray = alpha  # (N+1,)

        # Verify GL antisymmetry: α_{N+1/2} ≈ 0
        assert abs(alpha[N]) < 1e-12, (
            f"GL antisymmetry violated: α_{{N+1/2}} = {alpha[N]:.2e}"
        )

        # Cartesian stencil not used for spherical
        self.streaming_x = None
        self.streaming_y = None
        self.curvature = "spherical"
