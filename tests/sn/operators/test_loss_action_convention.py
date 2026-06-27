r"""S6.3 — the ``loss_action`` returns ``(L+C)ψ`` convention pin (#222).

After S6.3 the within-group loss action ``(L+C)ψ`` is the representation's walk
(:meth:`LossRepresentation.loss_action`); the operator's :meth:`apply` applies
the ONLY algebra glue, the Resolution-A collision subtraction
``L = (L+C) − C`` (``C = σ_t⊙``).  This file pins that convention as a TESTED
invariant (``sn_sweep_strategy.md`` §S6.1 locked-decision 2: the representation
returns the FULL within-group loss; the operator subtracts the collision
diagonal — returning ``L·ψ`` would re-couple ``C`` into every walk).

Two checks, complementary:

* **structural (NON-tautological).**  On a FULLY REFLECTIVE box a FLAT angular
  flux is the streaming fundamental mode ``L·ψ_flat = 0`` (the diamond-difference
  closure annihilates it — the SAME fact the analytical ``k_inf = 1.875`` anchor
  rests on).  So ``loss_action`` MUST return ``(L+C)ψ_flat = C·ψ_flat =
  σ_t·ψ_flat`` (collision only), and ``apply`` MUST return ``≈ 0``.  A
  representation that returned ``L·ψ`` (or a sign-flipped / double-counted ``C``)
  FAILS this — it is the check that ``loss_action`` is the FULL loss, not the bare
  streaming.  The ``−C`` glue checks alone are tautological (``apply`` is DEFINED
  as ``loss_action − σ_t·ψ``); this anchors the convention to an independent
  analytical fact.

* **the −C glue, cross-checked against an INDEPENDENT collision operator.**
  ``apply(ψ).bulk == loss_action(σ_t, ψ).bulk − C.apply(ψ).bulk`` on a ≥2G
  heterogeneous NON-FLAT config, where ``C`` is a SEPARATELY-constructed
  :class:`~orpheus.sn.operators.streaming.CollisionOperator` — so the subtrahend is verified
  to be exactly ``σ_t·ψ`` (the same σ_t), applied ONCE, not the production
  self-check.

``-O``-safe (vv Mode 8): ``np.testing`` / ``pytest.fail`` only.  ``foundation`` —
a software/algebra invariant (no theory ``:label:``).  ≥2G heterogeneous (vv
§H1/§H2: a 1-group uniform box makes ``σ_t·ψ`` flat and the ``(L+C)``-vs-``L``
distinction degenerate).
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.geometry import BC, CoordSystem, Mesh1D, Mesh2D
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.mesh.augmented_mesh import SNMesh
from orpheus.sn.operators.streaming import StreamingOperator
from orpheus.transport.operators.multiplication_operator import MultiplicationOperator
from orpheus.transport.fields.angular_flux import AngularFlux
from orpheus.transport.fields.boundary_flux import BoundaryFlux
from orpheus.transport.timed_full_field import TimedFullField
from tests.sn._test_helpers import placeholder_materials

pytestmark = [
    pytest.mark.foundation,
    # S5.5: the theory page labels for the math this file pins —
    # loss_action returns (L+C)psi (loss-rep-LpC) and the operator's
    # only glue is the Resolution-A subtraction (loss-rep-resolution-a).
    pytest.mark.verifies("loss-rep-LpC", "loss-rep-resolution-a"),
]


def _slab_2g(nx: int = 5) -> SNMesh:
    """1-D slab, reflective, 2G (CumprodScan → M_spatial._compute_LpC)."""
    mesh = Mesh1D(
        edges=np.linspace(0.0, 2.0, nx + 1),
        mat_ids=np.zeros(nx, dtype=int),
        coord=CoordSystem.CARTESIAN,
        bc_left=BC("reflective"), bc_right=BC("reflective"),
    )
    return SNMesh(mesh, Quadrature.gauss_legendre(4), placeholder_materials(ng=2))


def _cart2d_2g(nx: int = 4, ny: int = 5) -> SNMesh:
    """2-D Cartesian, reflective, 2G, NON-SQUARE (the d=2 representation
    walk — ScanMarch default since S6.9)."""
    mesh = Mesh2D(
        edges_x=np.linspace(0.0, 2.0, nx + 1),
        edges_y=np.linspace(0.0, 2.0, ny + 1),
        mat_map=np.zeros((nx, ny), dtype=int),
        coord=CoordSystem.CARTESIAN,
        bc_xmin=BC("reflective"), bc_xmax=BC("reflective"),
        bc_ymin=BC("reflective"), bc_ymax=BC("reflective"),
    )
    return SNMesh(mesh, Quadrature.level_symmetric(4), placeholder_materials(ng=2))


_CASES = {"slab_2g": _slab_2g, "cart2d_2g": _cart2d_2g}


def _zeros_state(sn: SNMesh) -> TimedFullField:
    return TimedFullField.zeros(bulk=AngularFlux, boundary=BoundaryFlux, mesh=sn)


@pytest.mark.parametrize("case", list(_CASES))
def test_loss_action_is_full_loss_LpC_flat_reflective(case):
    r"""[L0 structural] loss_action returns (L+C)ψ: flat-reflective → σ_t·ψ; apply → 0.

    The NON-tautological anchor.  A flat ψ on a fully-reflective box is the
    streaming fundamental (``L·ψ_flat = 0``), so ``loss_action`` MUST equal
    ``C·ψ_flat = σ_t·ψ_flat`` and ``apply.bulk`` MUST be ``≈ 0``.
    """
    sn = _CASES[case]()
    # Uniform σ_t so the flat-reflective fundamental is exact (L·ψ_flat = 0).
    sig_t = np.full((sn.ng, *sn.spatial_shape), 0.7)
    L = StreamingOperator(sn)        # pure σ-free streaming (#257 S8b)
    rep = L.loss_representation

    psi = _zeros_state(sn)
    psi.bulk.values[...] = 1.0   # FLAT in space + angle (the reflective fundamental)
    for face in psi.boundary.layout.faces:
        psi.boundary.face_view(face)[...] = 1.0   # flat trace (in = out = avg)

    lpc = rep.loss_action(sig_t, psi)
    lpsi = L.apply(psi)

    np.testing.assert_allclose(
        lpc.bulk.values, sig_t[None] * psi.bulk.values, rtol=0.0, atol=1e-12,
        err_msg=(
            f"[{case}] loss_action(flat) != σ_t·ψ — loss_action is NOT the FULL "
            "(L+C) loss (it must be C·ψ when streaming is annihilated)."
        ),
    )
    np.testing.assert_allclose(
        lpsi.bulk.values, 0.0, rtol=0.0, atol=1e-11,
        err_msg=(
            f"[{case}] apply(flat_reflective).bulk != 0 — the streaming "
            "fundamental L·ψ_flat = 0 is broken (k_inf annihilation)."
        ),
    )


@pytest.mark.parametrize("case", list(_CASES))
def test_pure_L_plus_C_recovers_loss_action_het(case):
    r"""[L0] the +C glue: (L_pure + C).apply == loss_action(σ_t), C INDEPENDENT (≥2G het).

    #257 S8b — the ``−C`` glue MOVED from inside ``L.apply`` to the
    ``OperatorSum`` (pure ``L`` + ``C``).  Two assertions pin the new
    composition:

    * ``(L + C).apply == loss_action(σ_t)`` BYTE-EXACT — the composite
      ``InvertibleOperator.apply`` IS ``loss_action(σ_t)`` (the full loss).
    * ``L.apply ≈ loss_action(σ_t) − C·ψ`` to FP ULP — the affine relation
      ``streaming_action(ψ) = M(σ_t)ψ − σ_t⊙ψ``, with ``C`` a SEPARATELY
      built :class:`CollisionOperator` so the subtrahend is the verified
      ``σ_t·ψ`` (same σ_t), applied ONCE.  pure ``L`` re-associates the FP
      reduction tree vs subtracting ``C·ψ`` off the σ-bearing matvec, so
      the agreement is to a few ULP, not bit-exact.
    """
    sn = _CASES[case]()
    rng = np.random.default_rng([sum(map(ord, case)), 2026])
    sig_t = rng.uniform(0.3, 3.0, size=(sn.ng, *sn.spatial_shape))   # heterogeneous
    L = StreamingOperator(sn)               # pure σ-free streaming (#257 S8b)
    C = MultiplicationOperator.from_mesh(sig_t, sn)        # INDEPENDENT collision operator
    rep = L.loss_representation

    psi = _zeros_state(sn)
    psi.bulk.values[...] = rng.standard_normal(psi.bulk.values.shape)   # NON-flat
    for face in psi.boundary.layout.faces:
        fv = psi.boundary.face_view(face)
        fv[...] = rng.standard_normal(fv.shape)

    lpc = rep.loss_action(sig_t, psi)
    lpsi = L.apply(psi)
    cpsi = C.apply(psi)
    composite = (L + C).apply(psi)          # InvertibleOperator.apply

    # non-vacuous guard: C·ψ must be non-trivial (else the pin is degenerate).
    if float(np.max(np.abs(cpsi.bulk.values))) < 1e-6:
        pytest.fail(f"[{case}] C·ψ ≈ 0 — the het/non-flat config degenerated.")

    # (1) the +C glue recovers the full loss BYTE-EXACT (the composite IS
    #     loss_action(σ_t)).
    np.testing.assert_array_equal(
        composite.bulk.values, lpc.bulk.values,
        err_msg=(
            f"[{case}] (L+C).apply != loss_action(σ_t) — the composite must "
            "recover the full within-group loss bit-exactly."
        ),
    )
    # (2) pure L.apply == loss_action(σ_t) − C·ψ to FP ULP (the affine
    #     relation streaming_action = M(σ_t) − σ_t⊙ψ).
    np.testing.assert_allclose(
        lpsi.bulk.values, lpc.bulk.values - cpsi.bulk.values,
        rtol=1e-12, atol=1e-12,
        err_msg=(
            f"[{case}] pure L.apply != loss_action − C·ψ — the affine "
            "streaming_action relation drifted beyond FP."
        ),
    )
