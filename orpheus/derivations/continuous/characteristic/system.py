r"""The multigroup Galerkin system: every group's transport block, the emission, and the questions they answer.

The emission density of group :math:`g`,

.. math::

    q_g = \sum_{g'} S_{g \leftarrow g'} \phi_{g'}
        + \frac{1}{k} \sum_{g'} F_{g \leftarrow g'} \phi_{g'} + q^{\mathrm{ext}}_g,
    \qquad \phi_g = \mathcal K_g q_g,

lives on the group's emission support (the regions where something is
emitted into :math:`g`), and the flux everywhere. The unknown is the
emission, :math:`q_g = \sum_{j \in \mathrm{supp}_g} q_{g,j} u_j`: testing
:math:`q = (S + F/k)\,\mathcal K q` against the support's basis functions
gives the Galerkin equation

.. math::

    W_s q = (S + F/k)\,K q + W_s q^{\mathrm{ext}},

with :math:`K = \operatorname{blockdiag}(K_g)` the transport blocks
(:class:`~.assembly.GroupTransport`, shape ``(G N, sum M_g)``),
:math:`W_s = \operatorname{blockdiag}(W_{\mathrm{supp}_g})` the mass matrix
restricted to each support, and :math:`S`, :math:`F` the emission matrices
applied node by node (shape ``(sum M_g, G N)``): the cross sections are
constant on a region and no panel crosses one, so the emission of a basis
function is that basis function times a constant, exactly. The flux is
read from the emission, :math:`W_G \phi = K q` with :math:`W_G = I_G
\otimes W` (:meth:`GalerkinSystem.flux`).

**Why the emission and not the flux** (the user's ruling of 2026-10-07,
P1 step (b) fourth rung, after a measurement). Both forms have the same
eigenvalues and the same flux. Their transposes differ: the adjoint of the
flux form :math:`\phi = \mathcal K E \phi` is :math:`E^\ast \phi^\dagger`
(in an infinite medium :math:`\Sigma_t \phi^\dagger`; measured on a closed
two-group sphere, a group ratio of 2.147 against the adjoint flux's
1.073), while the adjoint of the emission form :math:`q = E \mathcal K q`
is the adjoint flux :math:`\phi^\dagger = \mathcal K E^\ast \phi^\dagger`
itself, the importance of a source, on the supports where a source can be
posed. Each :math:`K_g` restricted to its support is symmetric
(reciprocity) and its rows off the support are the transposes of the
columns the support never assembles, so :math:`(E K)^{\mathsf T}` is the
Galerkin matrix of :math:`\mathcal K E^\ast` on the same basis, and the
transpose needs no metric.

**Two splittings of one operator.** The k question splits
:math:`W_s - (S + F/k) K` with the fission alone scaled: the pencil
:math:`(W_s - S K, F K)`, :attr:`GalerkinSystem.pencil`, whose fundamental
mode, higher modes and, by transposition, their adjoints are the eigen
answers. The source questions split it at :math:`k = 1` with every
secondary emission in the gain, :math:`(W_s, (S + F) K)`,
:attr:`GalerkinSystem.source_pencil`, whose least solution is the sum of
the Neumann series over collisions; its subcriticality check then refuses
a body made supercritical by its (n,2n) emission with no fission, which
the k pencil's least solution would solve directly into a negative flux.

**The emission space** (:class:`EmissionSpace`) has one coefficient per
support node of each group, group-major. Its restriction :math:`R`, of
shape ``(sum M_g, G N)``, is the one spelling of that layout: :math:`W_s =
R W_G R^{\mathsf T}`, and the emission matrices are :math:`R` composed with
the emission per node, :math:`S = R\,E_{\mathrm{node}}(\Sigma_s +
2\Sigma_2)`.

**The system is posed, never handed its blocks.** Its fields are the
problem (the basis, the walls, the cross sections, where sources are
posed) and the resolution; the transport blocks are derived from them, so
a block assembled for another group or another cross section cannot be
put in.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from functools import cached_property

import numpy as np
from scipy.linalg import block_diag

from orpheus.derivations.common.dense_pencil import DensePencil

from .assembly import GroupTransport, LineRule, TransportResolution
from .basis import PanelBasis
from .cross_sections import RegionCrossSections
from .reading import PointRule, angular_flux
from .walls import Walls


@dataclass(frozen=True, eq=False)
class EmissionSpace:
    r"""The emission space: each group's coefficients on its support nodes, group-major, ``(sum M_g,)``.

    Attributes
    ----------
    supports:
        Each group's support, the basis indices of its nodes, ``(M_g,)`` each, increasing.
    region:
        The region of each basis node, ``(N,)``: what a refusal names.
    """

    supports: tuple[np.ndarray, ...]
    region: np.ndarray

    @property
    def n_nodes(self) -> int:
        """The basis size :math:`N`."""
        return self.region.size

    @property
    def n_groups(self) -> int:
        """The number :math:`G` of groups."""
        return len(self.supports)

    @property
    def size(self) -> int:
        r"""The dimension :math:`\sum_g M_g`."""
        return sum(support.size for support in self.supports)

    @cached_property
    def restriction(self) -> np.ndarray:
        r"""The restriction :math:`R` from a field ``(G N,)`` to the emission space, ``(sum M_g, G N)``: a 0/1 selection."""
        r = np.zeros((self.size, self.n_groups * self.n_nodes))
        on = np.concatenate([g * self.n_nodes + support for g, support in enumerate(self.supports)])
        r[np.arange(self.size), on] = 1.0
        return r

    def restrict(self, field: np.ndarray) -> np.ndarray:
        r"""A field's nodal coefficients ``(G, N)`` on the emission space, ``(sum M_g,)``.

        Refuses a field with a non-zero coefficient off its group's support,
        where the emission space has no coefficient: it would be dropped.
        """
        field = np.asarray(field, dtype=float)
        if field.shape != (self.n_groups, self.n_nodes):
            raise ValueError(
                f"a field has one coefficient per group and node, ({self.n_groups}, {self.n_nodes}); got {field.shape}"
            )
        dropped = field.ravel() - self.restriction.T @ (self.restriction @ field.ravel())
        if np.any(dropped != 0.0):
            off = {
                g: sorted(set(self.region[np.flatnonzero(row)].tolist()))
                for g, row in enumerate(dropped.reshape(field.shape)) if row.any()
            }
            raise ValueError(
                f"the field is non-zero off the emission support, in the regions {off} (by group): "
                "pose the system with those source regions"
            )
        return self.restriction @ field.ravel()

    def split(self, coefficients: np.ndarray) -> tuple[np.ndarray, ...]:
        """Coefficients on the emission space, one array per group on its support nodes."""
        ends = np.cumsum([support.size for support in self.supports])[:-1]
        return tuple(np.split(np.asarray(coefficients), ends))


@dataclass(frozen=True, eq=False)
class GalerkinSystem:
    r"""The Galerkin form :math:`W_s q = (S + F/k) K q + W_s q^{\mathrm{ext}}` of a multigroup transport problem.

    Attributes
    ----------
    basis:
        The panel basis of every group's flux and emission.
    walls:
        The body's walls.
    cross_sections:
        The regions' total cross sections and emission matrices.
    resolution:
        The transport blocks' resolution.
    source_regions:
        Where a source is posed, ``(n, G)``, region-major as the cross
        sections: each group's support is its emission support
        (:meth:`~.cross_sections.RegionCrossSections.emission_support`)
        widened by these regions. By default none.
    """

    basis: PanelBasis
    walls: Walls
    cross_sections: RegionCrossSections
    resolution: TransportResolution
    source_regions: np.ndarray | None = field(default=None)

    def __post_init__(self) -> None:
        shape = self.cross_sections.total.shape
        if shape[0] != self.basis.regions.n_regions:
            raise ValueError(f"one cross section row per region of the basis, {self.basis.regions.n_regions}; got {shape[0]}")
        regions = np.zeros(shape, dtype=bool) if self.source_regions is None else np.array(self.source_regions, dtype=bool)
        if regions.shape != shape:
            raise ValueError(f"the source regions are a region-major {shape} mask, as the cross sections; got {regions.shape}")
        regions.flags.writeable = False
        object.__setattr__(self, "source_regions", regions)

    @property
    def n_groups(self) -> int:
        """The number :math:`G` of groups."""
        return self.cross_sections.n_groups

    @cached_property
    def groups(self) -> tuple[GroupTransport, ...]:
        """Each group's transport block, on its own line rule, its columns on the group's support."""
        support = self.cross_sections.emission_support() | self.source_regions
        r = self.resolution
        return tuple(
            LineRule.of(self.basis, self.walls, self.cross_sections.total[:, g], r.line_points).transport(
                r.points, r.inner_points, support[:, g]
            )
            for g in range(self.n_groups)
        )

    # ── the spaces ───────────────────────────────────────────────────────

    @cached_property
    def emission(self) -> EmissionSpace:
        """The emission space: the block columns of every group."""
        return EmissionSpace(tuple(group.support for group in self.groups), self.basis.region)

    @cached_property
    def mass(self) -> np.ndarray:
        r""":math:`W_G = I_G \otimes W`, the flux's metric, ``(G N, G N)``."""
        return np.kron(np.eye(self.n_groups), self.basis.mass)

    @cached_property
    def emission_mass(self) -> np.ndarray:
        r""":math:`W_s = R W_G R^{\mathsf T}`, the emission's metric, ``(sum M_g, sum M_g)``."""
        r = self.emission.restriction
        return r @ self.mass @ r.T

    # ── the matrices ─────────────────────────────────────────────────────

    @cached_property
    def transport(self) -> np.ndarray:
        r""":math:`K = \operatorname{blockdiag}(K_g)`, ``(G N, sum M_g)``: the emission space's columns."""
        return block_diag(*(group.block for group in self.groups))

    def _per_node(self, emission: np.ndarray) -> np.ndarray:
        r""":math:`E_{\mathrm{node}}`: a per-region emission table ``(n, G, G)`` ``[to, from]`` at each node, ``(G N, G N)``."""
        n_nodes = self.basis.size
        per_node = np.zeros((self.n_groups, n_nodes, self.n_groups, n_nodes))
        node = np.arange(n_nodes)
        per_node[:, node, :, node] = emission[self.basis.region]                    # (N, G, G) [node, to, from]
        return per_node.reshape(self.n_groups * n_nodes, self.n_groups * n_nodes)

    @cached_property
    def scattering(self) -> np.ndarray:
        r""":math:`S = R\,E_{\mathrm{node}}(\Sigma_s + 2\Sigma_2)`, ``(sum M_g, G N)``."""
        return self.emission.restriction @ self._per_node(self.cross_sections.scattering)

    @cached_property
    def fission(self) -> np.ndarray:
        r""":math:`F = R\,E_{\mathrm{node}}(\chi \otimes \nu\Sigma_f)`, ``(sum M_g, G N)``."""
        return self.emission.restriction @ self._per_node(self.cross_sections.fission)

    # ── the questions ────────────────────────────────────────────────────

    @cached_property
    def pencil(self) -> DensePencil:
        r"""The k pencil :math:`(W_s - S K, F K)` on the emission: :math:`F`'s emission is the one scaled by :math:`1/k`."""
        return DensePencil(self.emission_mass - self.scattering @ self.transport, self.fission @ self.transport)

    @cached_property
    def source_pencil(self) -> DensePencil:
        r"""The source questions' pencil :math:`(W_s, (S + F) K)`: every secondary emission is gain."""
        return DensePencil(self.emission_mass, (self.scattering + self.fission) @ self.transport)

    def flux(self, emission: np.ndarray) -> np.ndarray:
        r"""The flux of coefficients on the emission space, ``(G, N)``: :math:`W_G \phi = K q`."""
        return np.linalg.solve(self.mass, self.transport @ emission).reshape(self.n_groups, self.basis.size)

    def source_emission(self, source: np.ndarray) -> np.ndarray:
        r"""The emission of a source's nodal coefficients ``(G, N)``, on the emission space ``(sum M_g,)``.

        The least solution of :math:`W_s q = (S + F) K q + W_s q^{\mathrm{ext}}`:
        the source and every collision's emission.
        """
        return self.source_pencil.least_solution(self.emission_mass @ self.emission.restrict(source))

    def fixed_source(self, source: np.ndarray) -> np.ndarray:
        r"""The flux of a source's nodal coefficients ``(G, N)``, ``(G, N)``: :meth:`flux` of its :meth:`source_emission`."""
        return self.flux(self.source_emission(source))

    # ── the readings at a point ──────────────────────────────────────────

    def point_flux(self, position: float, emission: np.ndarray) -> np.ndarray:
        r"""The scalar flux at the point of orbit coordinate ``position``, ``(G,)``: :math:`(\mathcal K_g q_g)(x)`.

        ``emission`` is on the emission space; each group's emission is
        transported to the point once more (:class:`~.reading.PointRule`),
        on that group's own lines and with its walls' currents.
        """
        r = self.resolution
        return np.array(
            [
                PointRule.of(self.basis, self.walls, self.cross_sections.total[:, g], position, r.line_points).row(group, r) @ q
                for g, (group, q) in enumerate(zip(self.groups, self.emission.split(emission), strict=True))
            ]
        )

    def angular_flux(self, position: float, directions: np.ndarray, emission: np.ndarray) -> np.ndarray:
        r"""The angular flux per steradian at the point, ``(..., G)``, in the directions at box coordinates ``directions``.

        ``directions`` are coordinates of :meth:`~orpheus.geometry.chart.Chart.directions_at`'s box at the point,
        ``(..., k)`` (:func:`~.reading.angular_flux`); ``emission`` is on the emission space.
        """
        r = self.resolution
        per_group = [
            angular_flux(self.basis, self.walls, self.cross_sections.total[:, g], group, position, directions, r) @ q
            for g, (group, q) in enumerate(zip(self.groups, self.emission.split(emission), strict=True))
        ]
        return np.stack(per_group, axis=-1)

    def response(self, detector: np.ndarray) -> np.ndarray:
        r"""The adjoint flux of a detector's nodal coefficients ``(G, N)``, on the emission space ``(sum M_g,)``.

        The least solution of
        :math:`(W_s - K^{\mathsf T} (S + F)^{\mathsf T}) \phi^\dagger = K^{\mathsf T} r`,
        the Galerkin form of :math:`\phi^\dagger = \mathcal K (r + (S + F)^\ast \phi^\dagger)`:
        the importance of a source. The detector's reading of a source
        :math:`q` is :math:`\langle r, \phi(q) \rangle = \phi^{\dagger\mathsf T} W_s q`
        on the emission space.
        """
        detector = np.asarray(detector, dtype=float)
        if detector.shape != (self.n_groups, self.basis.size):
            raise ValueError(
                f"a detector has one coefficient per group and node, ({self.n_groups}, {self.basis.size}); got {detector.shape}"
            )
        load = self.transport.T @ detector.ravel()                  # <u_i, K r>, by reciprocity
        return self.source_pencil.adjoint().least_solution(load)


__all__ = ["EmissionSpace", "GalerkinSystem"]
