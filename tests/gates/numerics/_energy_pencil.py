r"""A multigroup operator pencil built from numerics objects alone — the fixture
the numerics law tests (pencil, posing, gauge, outcome) pose their laws on.

Layer DAG Phase 2 (``.claude/plans/layer_dag.md`` ruling 23): these laws are
numerics', so their fixture is built from numerics (and values/geometry)
objects only; it imports nothing from ``data``, ``mesh``, ``transport``,
``homogeneous`` or ``sn``.

The pencil :math:`(A, F)` is the infinite-medium multigroup pencil spelled in
the operator algebra without the transport vocabulary:

* the space is the energy axis alone, ``FunctionSpace.of_axes(EnergyAxis.synthetic(ng))``
  (shape ``(ng,)``), or a deliberately WEIGHTED toy axis where a law must see
  the metric. The homogeneous pose's ``Energy ⊗ point`` (shape ``(ng, 1)``) is
  NOT reproduced: since U1b every action parses its operand at its domain, and
  ``SparseAssembledOperator`` on the two-axis space returns ``(ng,)`` for an
  ``(ng, 1)`` operand, which its own ``.H`` then refuses (``[M]`` 2026-11-01,
  ``specs_2r/proto/probe_fixture.py``). The point factor carries nothing the
  pencil, posing, gauge and outcome laws read;
* :math:`A = C - S`, the collision diagonal minus a lower-triangular
  (down-scatter) transfer matrix — a NON-symmetric, non-normal loss (the
  regime the homogeneous fixture carried);
* :math:`F = \sum_m |\chi_m\rangle\langle\nu\Sigma_{f,m}|`, rank 0, 1 or 2, one operator per dyad,
  summed in the algebra;
* each operator is a :class:`~orpheus.numerics.assembled_operator.SparseAssembledOperator`
  BOUND to the space on both ends, so ``.H`` is the metric adjoint and the
  ends law compares real spaces. A dyad is a bound dense leaf
  :math:`\chi_m \otimes \nu\Sigma_{f,m}` and not the production dyad
  (``elementary_operators.outer`` / ``RankOneOperator``), because that carries no bound ends (its binding is
  deferred, #330), so its ``.H`` is unposed and the ``.H`` laws could not be read on it.

The data are refused at the mint (:class:`GroupData`) unless they are a legal infinite medium whose loss is
lower-triangular with a positive diagonal, and the pencil is refused on a space that is not the data's group space.

The carrier is the ``(ng,)`` vector, the exact shape of the space (an element
takes the exact shape since U1b). The rate functionals are co-vectors,
``space.flat(w)`` (:class:`~orpheus.numerics.space.Covector`, the member of
:math:`V^*` minted by the Riesz lowering), and the gauge reads their
``evaluate``: the metric enters once, at the mint.

The numbers are chosen, not read from a data library: what the laws need is
the STRUCTURE (≥ 2 groups, non-symmetric loss, a supercritical fissile case,
a rank-0 and a rank-2 production), and a library value would re-import the
data layer for nothing.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import cached_property
from typing import Callable

import numpy as np
import scipy.sparse as sp

from orpheus.numerics.assembled_operator import SparseAssembledOperator
from orpheus.numerics.axis import Axis, BasisKind, Declared, EnergyAxis
from orpheus.numerics.fredholm import PredictedKernel
from orpheus.numerics.operator import LinearOperator
from orpheus.numerics.pencil import OperatorPencil
from orpheus.numerics.posing import K_MAP, EigenPosing
from orpheus.numerics.space import Covector, FunctionSpace

@dataclass(frozen=True)
class GroupData:
    """The multigroup coefficients of one infinite medium (chosen, not library)."""

    sigma_t: np.ndarray            # (ng,)  [1/cm]  total
    transfer: np.ndarray           # (ng, ng) [1/cm] S[g_to, g_from], lower-triangular (down-scatter)
    spectra: tuple[np.ndarray, ...]       # each (ng,), the emission spectra χ_m
    productions: tuple[np.ndarray, ...]   # each (ng,) [1/cm], the νΣ_f,m rows

    def __post_init__(self) -> None:
        """Refuse what is not a legal infinite medium with a lower-triangular, nonsingular loss (each a ValueError)."""
        ng = self.sigma_t.shape[0]
        if self.transfer.shape != (ng, ng) or len(self.spectra) != len(self.productions):
            raise ValueError("GroupData: transfer must be (ng, ng) and every spectrum needs its production row")
        if np.any(self.transfer < 0.0):
            raise ValueError("GroupData: a transfer cross section is negative")
        if np.any(np.triu(self.transfer, 1) != 0.0):
            raise ValueError("GroupData: the transfer is not lower-triangular (down-scatter only)")
        if np.any(self.sigma_t - np.diag(self.transfer) <= 0.0):
            raise ValueError("GroupData: the loss diagonal sigma_t - S_gg is not positive")
        if np.any(self.sigma_a < 0.0):
            raise ValueError("GroupData: the absorption sigma_t - out-scatter is negative")
        for chi in self.spectra:
            if abs(float(np.sum(chi)) - 1.0) > 1e-12 or np.any(chi < 0.0):
                raise ValueError("GroupData: an emission spectrum chi_m is not a probability (non-negative, summing to 1)")

    @property
    def ng(self) -> int:
        return int(self.sigma_t.shape[0])

    @property
    def sigma_a(self) -> np.ndarray:
        """Absorption = total − out-scatter (column sums of the transfer matrix)."""
        return self.sigma_t - self.transfer.sum(axis=0)


#: The fissile cases (rank-1 production). ``2g`` is supercritical (k∞ > 1) so
#: the α map's fundamental α is positive.
FISSILE: dict[int, GroupData] = {
    1: GroupData(np.array([1.0]), np.array([[0.5]]), (np.array([1.0]),), (np.array([0.75]),)),
    2: GroupData(
        np.array([0.30, 0.80]),
        np.array([[0.18, 0.00], [0.07, 0.70]]),
        (np.array([1.0, 0.0]),),
        (np.array([0.02, 0.21]),),
    ),
    4: GroupData(
        np.array([0.30, 0.60, 0.90, 1.20]),
        np.array([
            [0.12, 0.00, 0.00, 0.00],
            [0.10, 0.30, 0.00, 0.00],
            [0.03, 0.15, 0.55, 0.00],
            [0.01, 0.05, 0.25, 1.05],
        ]),
        (np.array([0.60, 0.35, 0.05, 0.0]),),
        (np.array([0.01, 0.02, 0.08, 0.15]),),
    ),
}


def non_fissile(ng: int) -> GroupData:
    """The same loss with NO production (rank 0)."""
    d = FISSILE[ng]
    return GroupData(d.sigma_t, d.transfer, (), ())


def two_spectra(ng: int) -> GroupData:
    """Two independent emission dyads (rank 2; ``ng ≥ 2``) — a rank the
    one-spectrum families cannot reach."""
    d = FISSILE[ng]
    chi2 = np.zeros(ng); chi2[-1] = 1.0
    nu2 = np.linspace(0.05, 0.01, ng)
    return GroupData(d.sigma_t, d.transfer, (d.spectra[0], chi2), (d.productions[0], nu2))


def energy_space(ng: int) -> FunctionSpace:
    """The energy axis under the counting metric (a synthetic group structure)."""
    return FunctionSpace.of_axes(EnergyAxis.synthetic(ng))


def weighted_space(weights: np.ndarray) -> FunctionSpace:
    """A deliberately non-physical toy: a WEIGHTED group axis, so ``.H`` must
    read a non-scalar Gram (the 0-D counting pose cannot see a metric error)."""
    w = np.asarray(weights, dtype=float)
    return FunctionSpace.of_axes(Axis("weighted-group", (w.shape[0],), w, origin=Declared(BasisKind.NODAL)))


def bound(matrix: np.ndarray, space: FunctionSpace) -> SparseAssembledOperator:
    """A matrix as an operator BOUND to ``space`` on both ends."""
    return SparseAssembledOperator(sp.csr_array(np.asarray(matrix, dtype=float)), domain=space, codomain=space)


@dataclass(frozen=True, eq=False)
class EnergyPencil:
    """The pencil, its posing, its rate functionals and its exact eigenpair, on one space."""

    data: GroupData
    space: FunctionSpace

    def __post_init__(self) -> None:
        if self.space.shape != (self.data.ng,):
            raise ValueError(f"EnergyPencil: the space {self.space!r} is not the {self.data.ng}-group space of the data")

    @classmethod
    def of(cls, data: GroupData) -> "EnergyPencil":
        """The counting pose: the data on the synthetic energy axis. A weighted pose is ``EnergyPencil(data, space)``."""
        return cls(data, energy_space(data.ng))

    # ── the operators ────────────────────────────────────────────────
    @cached_property
    def collision(self) -> SparseAssembledOperator:
        return bound(np.diag(self.data.sigma_t), self.space)

    @cached_property
    def transfer(self) -> SparseAssembledOperator:
        return bound(self.data.transfer, self.space)

    @cached_property
    def loss(self) -> LinearOperator:
        """:math:`A = C - S`, an ``OperatorSum`` of two bound leaves."""
        return self.collision - self.transfer

    @cached_property
    def production(self) -> LinearOperator:
        r""":math:`F = \sum_m \chi_m \otimes \nu\Sigma_{f,m}`: one bound leaf per dyad, summed in the algebra; the
        bound zero when there is no dyad (rank 0). Built from the coefficients, never from :attr:`matrices`."""
        dyads = [bound(np.outer(chi, nu), self.space) for chi, nu in zip(self.data.spectra, self.data.productions)]
        if not dyads:
            return bound(np.zeros((self.data.ng, self.data.ng)), self.space)
        total: LinearOperator = dyads[0]
        for dyad in dyads[1:]:
            total = total + dyad
        return total

    @cached_property
    def predicted_kernel(self) -> PredictedKernel:
        r"""The loss's kernel is trivial: :math:`A = C - S` with :math:`S` lower-triangular is lower-triangular
        with diagonal :math:`\Sigma_{t,g} - S_{gg} > 0`, so :math:`\det A \neq 0` (refused otherwise by
        :class:`GroupData`)."""
        return PredictedKernel.trivial(self.space, "the loss is lower-triangular with a positive diagonal")

    @cached_property
    def pencil(self) -> OperatorPencil:
        return OperatorPencil(self.loss, self.production, kernel=self.predicted_kernel)

    @cached_property
    def eigen_posing(self) -> EigenPosing:
        return EigenPosing(self.pencil, K_MAP)

    # ── the rate functionals: co-vectors, minted by the Riesz lowering ───
    @cached_property
    def production_covector(self) -> Covector:
        r""":math:`\flat\,\nu\Sigma_f`, the member of :math:`V^*` pairing to :math:`\langle\nu\Sigma_f, \cdot\rangle_V`."""
        return self.space.flat(sum(self.data.productions, np.zeros(self.data.ng)))

    @cached_property
    def absorption_covector(self) -> Covector:
        r""":math:`\flat\,\Sigma_a`."""
        return self.space.flat(self.data.sigma_a)

    @property
    def production_rate(self) -> Callable[[np.ndarray], float]:
        r"""The production co-vector's ``evaluate``: a bound method of ONE cached co-vector, so two gauges built
        from it compare equal (``[M]`` ``cov.evaluate == cov.evaluate`` is True)."""
        return self.production_covector.evaluate

    @property
    def absorption_rate(self) -> Callable[[np.ndarray], float]:
        return self.absorption_covector.evaluate

    # ── the exact pair (dense eig of A⁻¹F, the trusted primitive) ────────
    @cached_property
    def matrices(self) -> tuple[np.ndarray, np.ndarray]:
        """(A, F) assembled from the COEFFICIENTS by numpy — the reference route. Independent of the operators on both
        counts: :attr:`loss` and :attr:`production` are built from the same coefficients by another route (bound
        leaves summed in the algebra), and this route never reads ``as_matrix`` (the old 0-D fixture read both sides
        of its affine rows through it)."""
        ng = self.data.ng
        A = np.diag(self.data.sigma_t) - self.data.transfer
        F = sum((np.outer(c, n) for c, n in zip(self.data.spectra, self.data.productions)), np.zeros((ng, ng)))
        return A, F

    @cached_property
    def exact_pair(self) -> tuple[float, np.ndarray]:
        """(k, φ): the dominant eigenpair of A⁻¹F, φ summing to +1."""
        A, F = self.matrices
        w, V = np.linalg.eig(np.linalg.solve(A, F))
        i = int(np.argmax(np.real(w)))
        phi = np.real(V[:, i])
        return float(np.real(w[i])), phi / phi.sum()
