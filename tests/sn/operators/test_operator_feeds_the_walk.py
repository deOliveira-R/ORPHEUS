r"""The POSED operator feeds the walk — hub mutations after posing are inert.

P4.9b step-2's keystone (§6c witness) plus the read-set contract gate.
Spec: ``scratch/p4_9b_verification_plan.md`` §2.1/§2.2; design:
``scratch/p4_9b_design.md`` §§7–9.

⛔ RED BEFORE P4.9b step 2, by construction — and that is this gate's whole
point.  Pre-carve the representation reads ``self.mesh.scheme`` /
``self.mesh.pole_angular_closure`` at APPLY time (43 class-(ii) reads,
``scratch/p4_9b_row45_recount.md``), so an operator posed BEFORE the swap
still marches with the mutant.  [M] 2026-08-28 at 10314dfa, rel deviation
5.000e-02 (slab / scheme), 4.596e-02 (cyl fp(4,6) / closure), 5.313e-02
(cyl fp(4,8)), 1.196e-01 (sphere).  After the carve the operator holds the
pre-swap objects and every row must read ``np.array_equal``.

⚠ ``cell_contribution`` alone is an INSUFFICIENT mutation: the ``.solve``
route consumes ``advance_psi_half`` plus the minted scan constants (P4.9a
Q1), so a gate mutating only the matvec's per-cell arm reads
``array_equal=True`` and certifies nothing.  Both closure surfaces are
mutated below (the march AND a minted constant).

⚠ The mutants are SUBCLASSES: the walk's residual family dispatch
(``isinstance(closure, MorelMontryAngularSweep)``,
``loss_representation/__init__.py:4213``) refuses a proxy or duck-typed
stand-in (verification plan F6).

⚠ The drive builds ``(L + C).solve`` in this file — never
``tests/sn/_test_helpers.sweep_once``, which re-poses internally and
therefore cannot witness a route claim about an ALREADY-POSED operator.

These are software route/invariant claims (no theory ``:label:``), so the
file carries ``foundation`` and no ``verifies()``.
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.geometry import BC, CoordSystem, Mesh1D
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.angular.closure import MorelMontryAngularSweep
from orpheus.sn.mesh.augmented_mesh import SNMesh
from orpheus.sn.operators.streaming import StreamingOperator
from orpheus.transport.fields.angular_boundary_flux import AngularBoundaryFlux
from orpheus.transport.fields.angular_flux import AngularFlux
from orpheus.transport.operators.multiplication_operator import (
    MultiplicationOperator,
)
from orpheus.transport.spatial.diamond import DiamondDifference
from orpheus.transport.timed_full_field import TimedFullField
from tests.sn._test_helpers import placeholder_materials

pytestmark = [pytest.mark.foundation]

_NG = 2
_NX = 8
_MEMO_SLOTS = ("_geom_cache", "_coll_cache", "_pole_mirror_cache")


# ═══════════════════════════════════════════════════════════════════════
# Fixtures — the 4-config grid of verification plan §2.1 (F7: the scheme
# row activates on the slab, the closure rows on the curvilinear charts).
# ═══════════════════════════════════════════════════════════════════════


def _slab_mesh() -> SNMesh:
    mesh = Mesh1D(
        edges=np.linspace(0.0, 2.0, _NX + 1),
        mat_ids=np.zeros(_NX, dtype=int),
        coord=CoordSystem.CARTESIAN,
        bc_left=BC("vacuum"),
        bc_right=BC("vacuum"),
    )
    return SNMesh(mesh, Quadrature.gauss_legendre(8), placeholder_materials(ng=_NG))


def _cylinder_mesh(n_phi: int) -> SNMesh:
    mesh = Mesh1D(
        edges=np.linspace(0.01, 2.0, _NX + 1),
        mat_ids=np.zeros(_NX, dtype=int),
        coord=CoordSystem.CYLINDRICAL,
        bc_left=BC("reflective"),
        bc_right=BC("vacuum"),
    )
    quad = Quadrature.folded_product(n_mu=4, n_phi=n_phi)
    return SNMesh(mesh, quad, placeholder_materials(ng=_NG))


def _sphere_mesh() -> SNMesh:
    mesh = Mesh1D(
        edges=np.linspace(0.0, 2.0, _NX + 1),
        mat_ids=np.zeros(_NX, dtype=int),
        coord=CoordSystem.SPHERICAL,
        bc_left=BC("reflective"),
        bc_right=BC("vacuum"),
    )
    return SNMesh(mesh, Quadrature.gauss_legendre(8), placeholder_materials(ng=_NG))


def _het_sigma(sn: SNMesh) -> np.ndarray:
    rng = np.random.default_rng(20260828)
    return rng.uniform(0.4, 2.5, size=(sn.ng, *sn.spatial_shape))


def _het_rhs(sn: SNMesh) -> TimedFullField:
    rhs = TimedFullField.zeros(
        interior=AngularFlux,
        boundary=AngularBoundaryFlux,
        space=sn.full_field_space,
    )
    rng = np.random.default_rng(20260829)
    rhs.interior.values[...] = rng.uniform(0.1, 1.0, size=rhs.interior.values.shape)
    return rhs


def _drive(sn: SNMesh, L: StreamingOperator) -> np.ndarray:
    """``(L + C).solve(rhs)`` from THIS operator — no internal re-posing."""
    C = MultiplicationOperator.from_mesh(_het_sigma(sn), sn)
    psi = (L + C).solve(_het_rhs(sn))
    return np.asarray(psi.interior.values).copy()


def _drop_memos(sn: SNMesh) -> None:
    # The cached tables would mask a hub swap pre-carve; post-carve the
    # memo has moved off the mesh, but the drop stays so the gate's green
    # can never be "the cache answered" (verification plan §2.1).
    for slot in _MEMO_SLOTS:
        sn.__dict__.pop(slot, None)


# ═══════════════════════════════════════════════════════════════════════
# Mutants — SUBCLASSES (F6), both surfaces of each family (F1's traps).
# ═══════════════════════════════════════════════════════════════════════


class _MutantDD(DiamondDifference):
    """DD with two mutated surfaces: emission and the residual kernel."""

    def source_emission(self, *args, **kwargs):
        return super().source_emission(*args, **kwargs) * 1.05

    def residual_kernel_batch(self, *args, **kwargs):
        # The base returns (values, aux_tuple); mutate the value member.
        values, aux = super().residual_kernel_batch(*args, **kwargs)
        return values * 1.05, aux


class _MutantMM(MorelMontryAngularSweep):
    """M-M with two mutated surfaces: the march AND a minted constant."""

    def advance_psi_half(self, *args, **kwargs):
        return super().advance_psi_half(*args, **kwargs) * 1.05

    @property
    def c_out_per_ordinate(self) -> np.ndarray:
        return MorelMontryAngularSweep.c_out_per_ordinate.fget(self) * 1.05  # type: ignore[attr-defined]


def _mutant_closure(sn: SNMesh) -> _MutantMM:
    assert sn.reduced is not None
    return _MutantMM(sn.reduced.angular, sn.reduced.redistribution_pairing)


_ROWS = [
    pytest.param("scheme_slab", _slab_mesh, "scheme", id="scheme_slab"),
    pytest.param("closure_cyl_deg", lambda: _cylinder_mesh(6), "closure", id="closure_cyl_deg"),
    pytest.param("closure_cyl", lambda: _cylinder_mesh(8), "closure", id="closure_cyl"),
    pytest.param("closure_sphere", _sphere_mesh, "closure", id="closure_sphere"),
]


def _swap(sn: SNMesh, slot: str) -> None:
    if slot == "scheme":
        sn.scheme = _MutantDD()
    else:
        sn.pole_angular_closure = _mutant_closure(sn)


@pytest.mark.xfail(
    strict=True,
    reason="P4.9b step-2 keystone, RED until the walk consumes the "
    "OPERATOR's objects — pre-carve the representation reads the hub at "
    "apply time, so the post-pose swap moves the answer (rel ~5e-2).",
)
@pytest.mark.parametrize("row, factory, slot", _ROWS)
def test_hub_mutation_after_posing_is_inert(row, factory, slot):
    """[foundation] A posed operator's answer ignores later hub mutations.

    Activation leg FIRST: a freshly posed operator over the MUTANT hub
    must MOVE — proving the mutated surface is consulted on this fixture
    at all (an inert mutant would make the main leg pass vacuously).
    """
    # ── ACTIVATION: fresh pose over a mutant hub moves the answer ──
    sn_ref = factory()
    base = _drive(sn_ref, StreamingOperator.pose(sn_ref))
    sn_mut = factory()
    _swap(sn_mut, slot)
    _drop_memos(sn_mut)
    moved = _drive(sn_mut, StreamingOperator.pose(sn_mut))
    assert not np.array_equal(moved, base), (
        f"row {row}: the mutant surface is not consulted on this fixture "
        "— the main leg below would be vacuous"
    )

    # ── THE CLAIM: mutate the hub AFTER posing — the answer must not move ──
    sn = factory()
    L = StreamingOperator.pose(sn)
    _ = L.loss_representation  # realize the selection now
    before = _drive(sn, L)
    assert np.array_equal(before, base)
    _swap(sn, slot)
    _drop_memos(sn)
    after = _drive(sn, L)
    assert np.array_equal(after, before), (
        f"row {row}: a hub mutation AFTER posing changed the answer — the "
        "walk is still reading the mesh, not the operator"
    )


# ═══════════════════════════════════════════════════════════════════════
# The read-set gate — F5's partition as an executable allowlist (§2.2).
# ═══════════════════════════════════════════════════════════════════════

# Space/layout facts the ruling keeps HUB-side (F5); everything method-
# flavored must arrive through the operator.  Q4 (ruled): the strategy-
# selection predicates are operator-side, so they must LEAVE this set.
_ALLOWED_HUB_SCHEME_READS = {"spatial_basis_per_axis", "is_multi_moment"}
_ALLOWED_HUB_CLOSURE_READS: set[str] = set()


def _recording_subclass(base_cls, reads: list[str]):
    class _Recorder(base_cls):
        def __getattribute__(self, name):
            if not name.startswith("__"):
                reads.append(name)
            return super().__getattribute__(name)

    _Recorder.__name__ = f"_Recording{base_cls.__name__}"
    return _Recorder


@pytest.mark.xfail(
    strict=True,
    reason="P4.9b step-2 read-set gate, RED until the walk consumes the "
    "OPERATOR's objects — pre-carve the hub route carries the per-cell "
    "kernels (source_emission, residual_kernel_batch, precompute_psi_state, "
    "level_indices, cell_contribution...).",
)
@pytest.mark.parametrize(
    "factory, do_matvec",
    [
        pytest.param(_slab_mesh, True, id="slab"),
        pytest.param(lambda: _cylinder_mesh(6), False, id="cyl_deg"),
    ],
)
def test_hub_route_reads_only_space_facts(factory, do_matvec):
    """[foundation] Post-carve the HUB route carries only the F5 space facts.

    The instrument: after posing (the operator captures the REAL
    objects), the hub's slots are re-bound to recording SUBCLASSES
    (delegating, behavior-identical, F6-safe).  Reads arriving via
    ``mesh.scheme.X`` / ``mesh.pole_angular_closure.X`` hit the
    recorders; reads via the operator's own fields hit the real objects
    and are invisible — exactly the partition the gate asserts.

    ⚠ The matvec leg is asserted to RUN (result shape) — the architect's
    own probe carried a vacuous curvilinear-matvec arm and recorded the
    lesson here rather than the silence (vv #17).
    """
    sn = factory()
    L = StreamingOperator.pose(sn)
    _ = L.loss_representation
    scheme_reads: list[str] = []
    closure_reads: list[str] = []
    sn.scheme = _recording_subclass(type(sn.scheme), scheme_reads)()
    if sn.reduced is not None:
        rec_cls = _recording_subclass(type(sn.pole_angular_closure), closure_reads)
        sn.pole_angular_closure = rec_cls(
            sn.reduced.angular, sn.reduced.redistribution_pairing,
        )
    _drop_memos(sn)

    # Instrument canary — the recorder records (activation of the gate).
    _ = sn.scheme.spatial_basis_per_axis
    assert scheme_reads == ["spatial_basis_per_axis"]
    scheme_reads.clear()

    out = _drive(sn, L)
    assert out.shape[0] == sn.quad.N  # the drive ran
    if do_matvec:
        C = MultiplicationOperator.from_mesh(_het_sigma(sn), sn)
        state = _het_rhs(sn)
        applied = (L + C).apply(state)
        assert applied.interior.values.shape == state.interior.values.shape

    bad_scheme = set(scheme_reads) - _ALLOWED_HUB_SCHEME_READS
    bad_closure = set(closure_reads) - _ALLOWED_HUB_CLOSURE_READS
    assert not bad_scheme and not bad_closure, (
        "the hub route carried method-flavored reads (must arrive through "
        f"the operator): scheme={sorted(bad_scheme)} "
        f"closure={sorted(bad_closure)}"
    )
