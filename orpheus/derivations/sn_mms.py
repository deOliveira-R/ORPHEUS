r"""Method of Manufactured Solutions (MMS) cases for SN verification.

MMS is a **closed-form** construction of an L1 spatial-convergence test:
we pick a smooth angular flux :math:`\psi_n(x)` that satisfies the
vacuum boundary conditions, substitute it into the transport operator,
and algebraically compute the external source :math:`Q^{\text{ext}}`
that would have produced it. The solver is then run on this source;
any deviation of the numerical flux from :math:`\psi_n` is pure
spatial-discretisation error and must decay at the design order of
the scheme (:math:`\mathcal{O}(h^{2})` for diamond difference).

**1D slab ansatz** (vacuum BCs on :math:`[0, L]`, 1 group):

.. math::

    \psi_n(x) = \frac{1}{W}\,A(x),
    \qquad A(x) = \sin\!\left(\tfrac{\pi x}{L}\right)

where :math:`W = \sum_n w_n = 2` for Gauss–Legendre. The flux is
isotropic in angle, so the scalar flux recovered by any quadrature
order equals :math:`\phi(x) = A(x)` exactly — the test isolates
spatial error from angular quadrature error.

**Manufactured source**. Substituting into

.. math::

    \mu_n\,\psi'_n + \Sigma_t\,\psi_n
    = \frac{1}{W}\!\left(\Sigma_s\,\phi + Q^{\text{ext}}_n\right)

and solving for :math:`Q^{\text{ext}}_n`:

.. math::

    Q^{\text{ext}}_n(x)
    = \mu_n\,A'(x) + \bigl(\Sigma_t - \Sigma_s\bigr)\,A(x)
    = \mu_n\,\frac{\pi}{L}\cos\!\left(\tfrac{\pi x}{L}\right)
      + \bigl(\Sigma_t - \Sigma_s\bigr)\sin\!\left(\tfrac{\pi x}{L}\right).

The :math:`W` factor cancels because the ansatz is already divided
by :math:`W`; the solver divides the isotropic and anisotropic source
slots by :math:`W` internally, so what we hand it is already the
full residual.

The BCs :math:`A(0)=A(L)=0` imply :math:`\psi_n=0` on both faces
for every ordinate — vacuum BCs are satisfied automatically, so no
inflow-flux bookkeeping is required by the caller.

.. seealso::

   - :doc:`/theory/discrete_ordinates` — MMS verification section
     with the full derivation and convergence-rate argument.
   - :func:`orpheus.sn.solve_sn_fixed_source` — consumer of the
     external source produced here.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.sparse import csr_matrix

from orpheus.data.macro_xs.mixture import Mixture
from orpheus.geometry import Mesh1D
from orpheus.sn.quadrature import GaussLegendre1D


@dataclass(frozen=True)
class SNSlabMMSCase:
    r"""Closed-form MMS fixed-source problem for 1D slab SN verification.

    Attributes
    ----------
    name : str
        Unique identifier, e.g. ``"sn_mms_slab_sin"``.
    sigma_t, sigma_s : float
        Total and isotropic scattering macroscopic cross sections
        (1-group, cm\ :sup:`-1`). The absorption ratio
        :math:`c = \\Sigma_s/\\Sigma_t` controls source-iteration
        convergence; :math:`c<1` is required.
    slab_length : float
        Physical length :math:`L` of the slab in cm.
    materials : dict[int, Mixture]
        Material map consumable by :class:`orpheus.sn.SNSolver`.
    mat_id : int
        Material ID assigned to every mesh cell.
    quadrature : GaussLegendre1D
        Angular quadrature (shared across mesh refinements so the
        convergence study isolates spatial error).
    tolerance : str
        Expected convergence order, e.g. ``"O(h^2)"``.
    equation_labels : tuple[str, ...]
        Sphinx labels exercised by tests built from this case.
    """

    name: str
    sigma_t: float
    sigma_s: float
    slab_length: float
    materials: dict[int, "Mixture"]
    mat_id: int
    quadrature: GaussLegendre1D
    tolerance: str = "O(h^2)"
    equation_labels: tuple[str, ...] = (
        "transport-cartesian",
        "dd-cartesian-1d",
        "dd-slab",
    )

    # ── Manufactured solution ─────────────────────────────────────────

    def phi_exact(self, x: np.ndarray) -> np.ndarray:
        r"""Scalar flux :math:`\phi(x) = \sin(\pi x/L)`."""
        return np.sin(np.pi * np.asarray(x) / self.slab_length)

    def dphi_exact(self, x: np.ndarray) -> np.ndarray:
        r"""Derivative :math:`A'(x) = (\pi/L)\cos(\pi x/L)`."""
        L = self.slab_length
        return (np.pi / L) * np.cos(np.pi * np.asarray(x) / L)

    # ── Mesh + source construction ────────────────────────────────────

    def build_mesh(self, n_cells: int) -> Mesh1D:
        """Uniform Cartesian slab mesh with ``n_cells`` equal cells."""
        edges = np.linspace(0.0, self.slab_length, n_cells + 1)
        mat_ids = np.full(n_cells, self.mat_id, dtype=int)
        return Mesh1D(edges=edges, mat_ids=mat_ids)

    def external_source(self, mesh: Mesh1D) -> np.ndarray:
        r"""Per-ordinate external source :math:`Q^{\text{ext}}_n` on ``mesh``.

        Evaluated at cell centres to match the diamond-difference
        cell-average convention. Returned shape is
        ``(N, nx, 1, 1)`` — per ordinate, per cell, one energy group.
        """
        x = mesh.centers                          # (nx,)
        A = self.phi_exact(x)                     # (nx,)
        Ap = self.dphi_exact(x)                   # (nx,)
        mu = self.quadrature.mu_x                 # (N,)
        N = len(mu)
        nx = len(x)

        streaming = mu[:, None] * Ap[None, :]     # (N, nx)
        removal = (self.sigma_t - self.sigma_s) * A[None, :]  # (1, nx)
        Q = streaming + removal                   # (N, nx)
        return Q[:, :, None, None]                # (N, nx, 1, 1)


# ═══════════════════════════════════════════════════════════════════════
# Case builders
# ═══════════════════════════════════════════════════════════════════════

def _make_1g_mixture(sigma_t: float, sigma_s: float) -> Mixture:
    """Build a minimal 1-group non-fissile mixture with capture = Σ_t − Σ_s.

    The solver builds sig_a internally from absorption_xs =
    SigC + SigL + SigF + Sig2_out. With no fission / (n,2n) / (n,α),
    setting SigC = Σ_t − Σ_s gives absorption = Σ_t − Σ_s (exactly
    the pure-absorber fraction that completes the Σ_t balance).
    """
    if sigma_s >= sigma_t:
        raise ValueError(
            f"Need Σ_s < Σ_t for a physical mixture (got "
            f"Σ_t={sigma_t}, Σ_s={sigma_s})."
        )

    ng = 1
    SigS0 = csr_matrix(np.array([[sigma_s]], dtype=float))
    Sig2 = csr_matrix(np.zeros((ng, ng)))
    return Mixture(
        SigC=np.array([sigma_t - sigma_s]),
        SigL=np.zeros(ng),
        SigF=np.zeros(ng),
        SigP=np.zeros(ng),
        SigT=np.array([sigma_t]),
        SigS=[SigS0],
        Sig2=Sig2,
        chi=np.zeros(ng),
        eg=np.array([1e-5, 2e7]),
    )


def build_1d_slab_mms_case(
    sigma_t: float = 1.0,
    sigma_s: float = 0.5,
    slab_length: float = 5.0,
    n_ordinates: int = 16,
    mat_id: int = 1,
    name: str = "sn_mms_slab_sin",
) -> SNSlabMMSCase:
    r"""Build the canonical 1D slab MMS case.

    Default parameters give :math:`c = \Sigma_s/\Sigma_t = 0.5`
    (source iteration converges in ~40 sweeps to 1e-12) and a slab
    about 5 mean free paths thick, which fits several wavelengths
    of the :math:`\sin(\pi x/L)` ansatz without being so optically
    thick that the manufactured source amplitude is uninteresting.
    """
    materials = {mat_id: _make_1g_mixture(sigma_t, sigma_s)}
    quadrature = GaussLegendre1D.create(n_ordinates=n_ordinates)
    return SNSlabMMSCase(
        name=name,
        sigma_t=sigma_t,
        sigma_s=sigma_s,
        slab_length=slab_length,
        materials=materials,
        mat_id=mat_id,
        quadrature=quadrature,
    )


def all_cases() -> list[SNSlabMMSCase]:
    """Return every registered MMS case (currently just the default)."""
    return [build_1d_slab_mms_case()]
