r"""**G1.1 / G1.3 / G1.4 / G1.5** — stage 1, *physics realization*: every
production leaf :math:`\{L, C, S, F, B\}` is **monomorphic** — it states one
arrow, refuses every other carrier the same way, and carries a metric-correct
Hilbert adjoint that cannot silently degrade.

Campaign: ``.claude/plans/operator_strategy_realization_campaign.md``.
Normative gate spec: ``.claude/plans/campaign_verification_plan.md`` §2.1.
Phase **P0** — gates only, NO production change.

Scope, and what is deliberately absent
======================================

* **G1.1** declared spaces — ``domain``/``codomain`` are concrete
  ``FunctionSpace``\ s, in VALUE and in ANNOTATION.
* **G1.3** uniform refusal — a wrong carrier is refused by a *typed*
  ``TypeError`` naming the operator and the carrier it wanted, never by a raw
  ``AttributeError`` leaking out of an unguarded attribute read.
* **G1.4** ``.H`` is the **G-adjoint**, not the Euclidean transpose:
  :math:`\langle A\psi,\varphi\rangle_G = \langle\psi, A^\dagger\varphi\rangle_G`.
* **G1.5** an **anonymous** leaf (built with no space) is unrepresentable, so
  ``.H`` can never see a ``None`` space.

**G1.2 (one arrow) and G1.6 (realizer purity) are OUT of scope for this
file.** G1.2 needs the carrier registry P2 builds — "exactly one carrier is
accepted" is only assertable once the set of carriers is enumerable; today
``F`` stands for 7 arrows and ``S`` for 6, and the *contract* (not the
dispatch spelling — #261 parks that ruling) can only be gated against a
registry that does not yet exist. G1.6 needs the realizer, which does not
exist at all.

**T = 1/v is absent from the leaf set** — MEASURED: ``git grep`` for
``inverse_velocity``/``InverseVelocity`` returns zero hits. The α-posing that
consumes it is campaign P4; when ``T`` lands it joins ``_LEAVES`` here.

The strict-xfail policy — these gates assert the POST-P1 world
==============================================================

Three of the four gates are **RED against ``main`` by design**. That is the
whole of ordering constraints O-2 and O-3: the leaf-arrow collapse RETIRES
arms and P1 makes the ``None``-space degradation *unrepresentable*, and a
degradation fixed before its catcher exists is a degradation nobody can prove
was ever possible. So the catcher is committed FIRST, red, as
``xfail(strict=True)``. Strict means an XPASS is a hard failure: the fix
cannot land silently, and the reader is told — in the marker's own ``reason``
— which campaign phase flipped it and to delete the marker.

**The strict-xfail set IS this file's todo list.** Rows that are already green
ship WITHOUT a marker, as a regression floor.

Measured state of the tree (``main`` @ ``b0a003b4``, host ``.venv``)
====================================================================

===========  ==============================  ===================================
gate         today                           measurement
===========  ==============================  ===================================
G1.1 value   **GREEN** on the SN ladder       all 5 leaves × 4 geometries declare
             **RED** off it (**R1**)          a concrete ``FullFieldSpace``; the
                                              model-generic construction
                                              (``homogeneous/solver.py:193``)
                                              reports ``domain is None``
G1.1 ann.    **RED** (R1's static face)        every leaf annotates
                                              ``FunctionSpace | None`` /
                                              ``Optional[FunctionSpace]``
G1.3         **RED on ``B`` only** (**R6**)    ``L``/``C``/``S``/``F`` raise a
                                              ``TypeError`` naming themselves;
                                              ``SNBoundaryOperator`` raises
                                              ``AttributeError: 'CoupledField'
                                              object has no attribute
                                              'interior'`` (``boundary.py:343``)
G1.4         **GREEN**                         residual ≤ 4.5e-14 across every
                                              leaf × geometry (rtol 1e-12 ⟹ 22×
                                              headroom)
G1.5         **RED** for ``C``/``S``/``F``      all three construct happily with
             (**R1**/**R2**)                   no space, and the ``.H`` they then
                                              build is **bit-identical** to the
                                              bare Euclidean transpose
===========  ==============================  ===================================

The mutation register (§7) — every gate's teeth
===============================================

**M-8** (G1.1) and **M-9** (G1.3): *today's state IS the mutation.* The
register says so explicitly, and it is why no monkeypatch is written for
them — the xfail rows below ARE the mutation proof, executed on every run.
When P1/P2 flips them green, the mutation "restore ``domain -> FunctionSpace
| None`` and build anonymously" / "restore the raw ``AttributeError``" is
literally a revert of that phase, and the strict XPASS at flip time is the
evidence the gate was live the whole way.

**M-10** (G1.4) is a real in-process mutation and is written out in full
below, in **both halves**:

* :func:`test_reciprocity_metric_is_load_bearing` — drop the metric inside
  ``_AdjointOperator.apply`` and the residual jumps to **O(1)** on the
  non-uniform curvilinear legs (MEASURED sphere ``L`` 1.40 / ``S`` 1.01 /
  ``F`` 1.46; cylinder 1.11 / 0.947 / 5.83e-2), against a clean ≤ 4.5e-14.
* :func:`test_a_globally_constant_metric_makes_reciprocity_blind` — the same
  mutation on a deliberately **flat-metric** slab is silent for EVERY leaf
  (MEASURED ≤ 4.3e-16). That is the config-blindness proof, and it is not
  optional: without it the gate could be passing for a reason that has
  nothing to do with the metric.

Two REFUTATIONS of the plan's stated G1.4 config (measured, not assumed)
=======================================================================

The gate spec asks for "non-uniform ``h``, curvilinear so ``V_cell`` spans an
order of magnitude", with the control leg phrased as "the uniform-``h`` leg
stays green". Both halves needed correcting against the tree:

1. **"Uniform ``h``" is NOT the blind config.** A uniform-``h`` slab under
   ``gauss_legendre(4)`` still REDs under M-10 at 1.3e-1 (``L``), 4.0e-1
   (``S``), 2.7e-1 (``F``) — because the *quadrature weights* vary, and the
   metric is :math:`G = V_{\rm cell}\,w_n`. The blindness is not about ``h``;
   it is about the metric being a **global constant**, since
   :math:`G = c\,\mathbb{1} \Rightarrow G^{-1}A^{\mathsf T}G = A^{\mathsf T}`
   identically. The control leg here is therefore built to BE that constant
   (:func:`_flat_metric_slab`), and it asserts the constancy as a
   precondition so it cannot silently stop being the blind config.
2. **M-10 constrains only ``{L, S, F}``.** ``C`` and ``B`` are
   *exercised-but-unconstrained* (``vv-principles`` **Mode 10**) — and
   necessarily so, not by a fixture accident:

   * ``C`` is a ``MultiplicationOperator``: **diagonal**, so
     :math:`C^{\mathsf T} = C` and :math:`G^{-1}CG = C` for ANY diagonal
     metric. No configuration exists in which dropping the metric moves
     ``C``'s reciprocity.
   * ``B`` with a specular law is a signed **permutation** of the trace that
     maps :math:`\mu \to -\mu`, preserving both :math:`|\Omega\cdot n|` and
     :math:`w_n` — so it commutes with the trace metric. ``SNMesh`` accepts
     only ``reflective``/``vacuum`` face laws (MEASURED: ``white`` raises
     ``ValueError``), and the vacuum law is the zero map, which commutes with
     everything. So there is no *reachable* config either.

   Rather than claim teeth these rows do not have, the honest closure is
   :func:`test_reciprocity_row_is_non_vacuous`: a SECOND mutation (scale the
   adjoint by 2) reddens every leaf including ``C`` and ``B`` at exactly
   0.5 relative, proving no G1.4 row is a dead gate — it pins
   ``apply_transpose``'s structure even where it cannot see the metric.

   ``B``'s metric weighting is not left unpinned: it is covered by the
   existing L11 wrong-metric control in
   ``tests/sn/operators/test_g_adjoint_reciprocity.py``, which drops
   :math:`|\Omega\cdot n|` from the *reference* metric on the composite
   ``(L + C - B)``. Not duplicated here.

Why this file builds its own fixtures rather than reusing the neighbours
========================================================================

``tests/sn/architecture/_config.py`` is the campaign's shared config home and
is reused wherever it serves — but its ``anisotropic_mixture`` is
deliberately **non-fissile** (``SigP = 0``, ``chi = 0``), because the
within-group record it was built for carries no fission piece. The LEAF
ladder includes ``F``, and on a non-fissile mixture ``F.apply(ψ) ≡ 0``, which
would make every ``F`` row of G1.4 the tautology ``0 == 0``. So the fissile
row lives here. **When the P0 files consolidate, ``_fissile_anisotropic_2g``
belongs beside ``anisotropic_mixture`` in ``_config.py``** (Pattern 2 — one
home for the mandatory config).

``tests/sn/operators/test_removal_form_matvec_sweep.py``'s ``_slab`` /
``_sphere`` / ``_cyl`` / ``_cart2d`` do not serve either, for the same
reason plus one more: they carry ``placeholder_materials``, whose ``SigS``,
``chi`` and ``SigP`` are all identically zero — so BOTH ``S`` and ``F`` would
be the zero operator and both G1.4 rows would be vacuous. They are also
uniform-``h``. The structurally-independent metric oracle
:func:`~tests.sn._test_helpers.g_inner` IS reused (it is the designated
shared home, and it is built from raw ``omega_dot_n`` / ``quad.weights`` /
``volumes`` rather than from the production metric under test — anti-R1).

Mandatory configuration (§8), checked per gate
==============================================

≥2G everywhere · **asymmetric** ``SigS`` in both Legendre orders (a
group-index transpose is invisible on a symmetric operator — Mode 12) ·
``Sig2 ≠ 0`` and ``SigL ≠ 0`` via a **direct ``Mixture``** construction
(``make_mixture`` hardcodes ``SigL = 0`` and defaults ``Sig2 = 0`` — lessons
L1) · **heterogeneous**, 2 regions · **P1 anisotropic** · **non-uniform**
radial/axial ``h``, **non-square** ``nx != ny`` in 2-D · **mixed**
vacuum + reflective BC with a **non-zero** trace on the reflecting face ·
fixed-seed **random, non-flat** probe states filled in BULK *and* TRACE ·
the geometry ladder slab **and** sphere **and** cylinder **and** 2-D
Cartesian (slab alone is the degenerate curvilinear case).

Runtime mode
============

Canonical ``python -O -m pytest``. Every check is ``pytest.fail`` /
``np.testing.assert_*`` — a function call, which fires under ``-O`` (a bare
``assert`` in a *helper* module would be stripped to a no-op; ``vv-principles``
Mode 8). Mutations are applied with ``monkeypatch``, IN PROCESS — never by
touching a file on disk.

Marks
=====

``foundation`` — software/architecture invariants of the operator algebra.
No theory ``:label:``, no ``verifies(...)`` (the verifies⊥level doctrine).
"""
from __future__ import annotations

import inspect
from typing import TYPE_CHECKING

import numpy as np
import pytest
from scipy.sparse import csr_matrix

from orpheus.data.macro_xs.mixture import Mixture
from orpheus.derivations.common.xs_library import get_mixture
from orpheus.geometry import BC, CoordSystem, Mesh1D
from orpheus.geometry.mesh import Mesh2D
from orpheus.numerics import operator as _operator_module
from orpheus.numerics.coupled_system import CoupledField
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.mesh.augmented_mesh import SNMesh
from orpheus.sn.operators.boundary import SNBoundaryOperator
from orpheus.sn.operators.streaming import StreamingOperator
from orpheus.transport.fields.angular_boundary_flux import AngularBoundaryFlux
from orpheus.transport.fields.angular_flux import AngularFlux
from orpheus.transport.mesh.material_mesh import MaterialMesh
from orpheus.transport.operators.fission import FissionOperator
from orpheus.transport.operators.multiplication_operator import (
    MultiplicationOperator,
)
from orpheus.transport.operators.scattering import ScatteringOperator
from orpheus.transport.timed_full_field import TimedFullField
from tests.sn.architecture._config import anisotropic_mixture
from tests.sn._test_helpers import g_inner

if TYPE_CHECKING:
    from orpheus.numerics.operator import LinearOperator

pytestmark = pytest.mark.foundation


# ═════════════════════════════════════════════════════════════════════════
# Contracts, tolerances and the xfail reasons — named, never inlined
# ═════════════════════════════════════════════════════════════════════════

#: G1.4's contract (gate spec §2.1). MEASURED worst clean residual across the
#: whole ladder: 4.50e-14 (``C`` on 2-D Cartesian) — 22× headroom.
_RECIPROCITY_RTOL = 1e-12

#: A mutation must move the residual by ORDERS, not by a factor. MEASURED
#: smallest constrained jump: 5.83e-02 (``F``, cylinder) — 58× headroom over
#: this floor, and 12 orders above the clean residual.
_MUTATION_FLOOR = 1e-3

_SEED_X = 20260729
_SEED_Y = 20260730

_R1_XFAIL = pytest.mark.xfail(
    strict=True,
    reason=(
        "R1 (annotation face) — a leaf's `domain`/`codomain` are still "
        "ANNOTATED `FunctionSpace | None`; the Optional stays legal until "
        "campaign 1's CS4 flips it to mandatory. (The VALUE rows this "
        "constant used to guard were deleted at CS1 step 3b: the "
        "model-generic path now threads a real space — the successor floor "
        "is tests/homogeneous/test_operator_spaces.py.) WHEN THIS XPASSES: "
        "CS4 has landed — delete this marker."
    ),
)

_R2_XFAIL = pytest.mark.xfail(
    strict=True,
    reason=(
        "R1/R2 — an ANONYMOUS leaf (no space) constructs happily, and the "
        "`.H` it then builds silently degrades to a bare Euclidean transpose: "
        "`_AdjointOperator.apply` (numerics/operator.py:1221-1227) applies the "
        "metric only when the inner space is non-None. MEASURED bit-identical "
        "to `apply_transpose`. Flipped by campaign P1 (construction without a "
        "space RAISES, so `.H` can never see a None space). WHEN THIS "
        "XPASSES: P1 has landed — delete this marker."
    ),
)

_R6_XFAIL = pytest.mark.xfail(
    strict=True,
    reason=(
        "R6 — the carrier guard is NON-UNIFORM across leaves: L/C/S/F raise a "
        "typed TypeError naming themselves and the carrier they wanted, while "
        "SNBoundaryOperator leaks a raw `AttributeError: 'X' object has no "
        "attribute 'interior'` from the unguarded read at "
        "sn/operators/boundary.py:343. Flipped by campaign P1/P2 (the leaf "
        "declares its arrow and refuses uniformly). WHEN THIS XPASSES: the "
        "guard has been unified — delete this marker."
    ),
)


# ═════════════════════════════════════════════════════════════════════════
# Materials — a DIRECT Mixture build (lessons L1: `make_mixture` nulls
# `SigL` AND `Sig2`, and offers no P1 channel)
# ═════════════════════════════════════════════════════════════════════════

def _fissile_anisotropic_2g(
    *,
    sig_t: "list[float]",
    sig_s0: "list[list[float]]",
    sig_s1: "list[list[float]]",
    sig_f: "list[float]",
    chi: "list[float]",
    nu: float = 2.6,
) -> Mixture:
    r"""This file's fissile row of the campaign's ONE mixture builder.

    Delegates to :func:`~tests.sn.architecture._config.anisotropic_mixture`
    (Pattern 2 — one home for the mandatory config; a second ``Mixture``
    builder inside the package whose whole purpose is "one home" would be
    the twin it exists to prevent).  What this row adds over a within-group
    fixture is the channels G1.4 needs ACTIVATED and ``make_mixture`` nulls:

    * ``sig_f``/``chi``/``nu`` — producing, so ``F.apply(psi) != 0`` and the
      ``F`` reciprocity rows are not the tautology ``0 == 0``;
    * ``SigL != 0`` — the (n,alpha) channel;
    * an **asymmetric** ``Sig2 != 0`` — ``S``'s (n,2n) term, a separate
      ell=0 block outside the Legendre fold.

    The non-fissile default in ``_config`` is correct for the *within-group*
    gates (fission enters as ``q_ext`` at the outer), which is exactly why
    this row must state its extra channels rather than inherit them.
    """
    return anisotropic_mixture(
        sig_t, sig_s0, sig_s1,
        sig_f=sig_f, chi=chi, nu=nu,
        sig_l=[0.004, 0.011],                      # (n,alpha) — non-zero
        sig_2=[[0.0, 0.03], [0.01, 0.0]],          # asymmetric (n,2n)
    )


def _two_region_fissile() -> "dict[int, Mixture]":
    """Two materials — heterogeneous nulls no spatial-distribution bug."""
    return {
        0: _fissile_anisotropic_2g(
            sig_t=[1.1, 2.3],
            sig_s0=[[0.38, 0.10], [0.05, 0.90]],
            sig_s1=[[0.02, 0.01], [0.00, 0.04]],
            sig_f=[0.02, 0.31], chi=[0.95, 0.05],
        ),
        1: _fissile_anisotropic_2g(
            sig_t=[0.9, 2.6],
            sig_s0=[[0.22, 0.03], [0.12, 1.10]],
            sig_s1=[[0.05, 0.02], [0.01, 0.03]],
            sig_f=[0.05, 0.12], chi=[0.80, 0.20],
        ),
    }


# ═════════════════════════════════════════════════════════════════════════
# The geometry ladder — NON-UNIFORM h everywhere (§8)
# ═════════════════════════════════════════════════════════════════════════

#: Deliberately non-uniform, and geometrically graded: on a sphere this makes
#: ``V_cell`` span 3.36e3 (MEASURED), so the bulk metric ``V·w`` is very far
#: from the constant that would make G1.4 a Euclidean-transpose check.
_NONUNIFORM_EDGES = np.array([0.0, 0.12, 0.35, 0.80, 1.30, 2.00])
_TWO_REGION_IDS = np.array([0, 0, 1, 1, 1])


def _slab() -> SNMesh:
    """1-D Cartesian, non-uniform ``h``, mixed reflective/vacuum, GL S4."""
    mesh = Mesh1D(
        edges=_NONUNIFORM_EDGES, mat_ids=_TWO_REGION_IDS,
        coord=CoordSystem.CARTESIAN,
        bc_left=BC("reflective"), bc_right=BC("vacuum"),
    )
    return SNMesh(
        mesh, Quadrature.gauss_legendre(n_ordinates=4), _two_region_fissile(),
    )


def _sphere() -> SNMesh:
    r"""1-D spherical — **the load-bearing G1.4 leg**.

    ``V_cell`` spans 3.36e3 (MEASURED) and ``gauss_legendre(4)`` weights are
    non-constant, so the metric :math:`G = V_{\rm cell} w_n` varies along
    BOTH axes. This is the config in which M-10 reddens ``L``, ``S`` and
    ``F`` at O(1). A ``level_symmetric`` rule would NOT do: its weights are
    all equal (MEASURED — one unique value), which leaves ``S`` and ``F``
    metric-blind because a purely spatial metric commutes with a
    space-diagonal, angle-mixing operator.

    Reflective at ``r = R`` so ``B`` acts on a live, non-zero trace.
    """
    mesh = Mesh1D(
        edges=_NONUNIFORM_EDGES, mat_ids=_TWO_REGION_IDS,
        coord=CoordSystem.SPHERICAL, bc_right=BC("reflective"),
    )
    return SNMesh(
        mesh, Quadrature.gauss_legendre(n_ordinates=4), _two_region_fissile(),
    )


def _cylinder() -> SNMesh:
    r"""1-D cylindrical, ``product(n_mu=4, n_phi=8)``.

    The second curvilinear leg, on a DIFFERENT quadrature family, so the
    M-10 result is not an accident of one rule's weight distribution.
    ``V_cell`` spans 160 (MEASURED); the product rule's ``mu`` weights vary.
    """
    mesh = Mesh1D(
        edges=_NONUNIFORM_EDGES, mat_ids=_TWO_REGION_IDS,
        coord=CoordSystem.CYLINDRICAL, bc_right=BC("reflective"),
    )
    return SNMesh(
        mesh, Quadrature.folded_product(n_mu=4, n_phi=8), _two_region_fissile(),
    )


def _cart2d() -> SNMesh:
    """2-D Cartesian, NON-SQUARE and non-uniform on both axes, mixed BC.

    ``nx != ny`` with unequal spacings is the ``x↔y``-swap catcher;
    ``level_symmetric`` avoids the ``mu_y == 0`` GL rank mismatch and the
    ERR-056 axis-aligned degeneracy.
    """
    mat_map = np.zeros((4, 3), dtype=int)
    mat_map[2:, :] = 1
    mesh = Mesh2D(
        edges_x=np.array([0.0, 0.1, 0.35, 0.9, 1.6]),
        edges_y=np.array([0.0, 0.2, 0.5, 1.4]),
        mat_map=mat_map, coord=CoordSystem.CARTESIAN,
        bc_xmin=BC("reflective"), bc_xmax=BC("vacuum"),
        bc_ymin=BC("reflective"), bc_ymax=BC("vacuum"),
    )
    return SNMesh(
        mesh, Quadrature.level_symmetric(sn_order=4), _two_region_fissile(),
    )


def _flat_metric_slab() -> SNMesh:
    r"""The CONFIG-BLINDNESS control: a **globally constant** metric.

    Reciprocity is blind to the metric exactly when :math:`G = c\,\mathbb{1}`,
    since then :math:`G^{-1}A^{\mathsf T}G = A^{\mathsf T}` for every ``A``.
    Three constants must coincide for that to hold on the bulk ⊕ trace
    composite:

    * ``V_cell``  — uniform ``h`` on a Cartesian slab;
    * ``w_n``     — ``gauss_legendre(2)``, whose two weights are BOTH exactly
      1 (the only SN-usable rule in the tree with equal weights);
    * ``|Ω·n| w`` — the 2-point GL nodes are :math:`\pm 1/\sqrt3`, so the
      trace weight is :math:`1/\sqrt3` on every ordinate and face.

    Choosing ``h = 1/√3`` makes the bulk constant equal the trace constant,
    so the metric is ONE number across the whole composite (MEASURED: a
    single unique value, 0.5773502691896258).

    This is not a contrivance — it is the exact statement of what the
    non-uniform curvilinear leg is buying. :func:`_assert_metric_is_constant`
    pins the property so the leg cannot silently stop being blind.
    """
    h = 1.0 / np.sqrt(3.0)
    mesh = Mesh1D(
        edges=np.arange(6) * h, mat_ids=_TWO_REGION_IDS,
        coord=CoordSystem.CARTESIAN,
        bc_left=BC("reflective"), bc_right=BC("reflective"),
    )
    return SNMesh(
        mesh, Quadrature.gauss_legendre(n_ordinates=2), _two_region_fissile(),
    )


_GEOMETRIES = {
    "slab": _slab,
    "sphere": _sphere,
    "cylinder": _cylinder,
    "cart2d": _cart2d,
}

#: The curvilinear legs M-10 is measured on — the only ones whose metric
#: varies along BOTH the spatial and the angular axis.
_CURVILINEAR = ("sphere", "cylinder")


# ═════════════════════════════════════════════════════════════════════════
# The leaf set and the probe states
# ═════════════════════════════════════════════════════════════════════════

#: Stage-1 production leaves (gate spec §2.1). ``T = 1/v`` is absent from the
#: tree — see the module docstring.
_LEAVES = ("L", "C", "S", "F", "B")

#: The leaves whose constructor admits an ANONYMOUS build (an optional space
#: parameter defaulting to ``None``). ``L`` and ``B`` derive their space from
#: the ``SNMesh`` they are handed and cannot be anonymous — pinned by
#: :func:`test_mesh_derived_leaves_carry_no_anonymous_construction_surface`.
_ANONYMOUS_CAPABLE = ("C", "S", "F")

#: MEASURED (M-10, sphere + cylinder): the leaves whose reciprocity the
#: metric actually constrains. ``C`` (diagonal) and ``B`` (a metric-preserving
#: specular permutation) commute with the metric ALGEBRAICALLY, so no config
#: makes them see it — see the module docstring's refutation #2.
_METRIC_CONSTRAINED = ("L", "S", "F")


def _leaf_set(sn_mesh: SNMesh) -> "dict[str, LinearOperator]":
    r"""The five production leaves, built exactly as the SN solver builds them.

    ``L``/``C`` mirror ``build_streaming_collision``
    (``sn/coupled_system.py:376``), ``S``/``F`` mirror ``SNSolver.__init__``
    (``sn/solver.py:1035`` / ``:1041``) and ``B`` mirrors
    ``build_within_group_system`` (``sn/coupled_system.py:464``) — each
    threading ``sn_mesh.full_field_space``, which is precisely why the G1.1
    VALUE rows are green on this ladder while the model-generic construction
    is not.
    """
    mat_xs = sn_mesh.material_xs_field()
    return {
        "L": StreamingOperator(sn_mesh),
        "C": MultiplicationOperator(
            coefficient=mat_xs.total_cross_section_field,
            space=sn_mesh.full_field_space,
        ),
        "S": ScatteringOperator.from_solver_data(
            mat_xs=mat_xs, quadrature=sn_mesh.quad, scattering_order=1,
            space=sn_mesh.full_field_space,
        ),
        "F": FissionOperator.from_solver_data(
            mat_xs=mat_xs, space=sn_mesh.full_field_space,
        ),
        "B": SNBoundaryOperator(sn_mesh),
    }


def _random_composite(sn_mesh: SNMesh, *, seed: int) -> TimedFullField:
    """Fixed-seed random state, filled in BULK **and** TRACE.

    A flat ψ nulls the streaming coupling and a zero trace nulls ``B``; both
    would make the reciprocity rows measure a smaller operator than the one
    they name (§8, probe-state lever).
    """
    rng = np.random.default_rng([seed, 7])
    state = TimedFullField.zeros(
        interior=AngularFlux, boundary=AngularBoundaryFlux, mesh=sn_mesh,
    )
    state.interior.values[...] = rng.standard_normal(state.interior.values.shape)
    for face in state.boundary.layout.faces:
        view = state.boundary.face_view(face)
        view[...] = rng.standard_normal(view.shape)
    return state


class _AlienCarrier:
    """A carrier no leaf's arrow accepts, and that is all it is."""


def _wrong_carrier(kind: str, sn_mesh: SNMesh) -> object:
    """The two wrong carriers G1.3 probes.

    ``alien`` is a bare object — the minimal statement of "not my domain".
    ``coupled_field`` is the R6 evidence carrier and the one that actually
    arises: the campaign's coupled grid hands ``CoupledField``\\ s around, and
    ``SNBoundaryOperator``'s unguarded ``psi.interior`` read on one is
    literally the measured ``AttributeError``.
    """
    if kind == "alien":
        return _AlienCarrier()
    return CoupledField(systems=(_random_composite(sn_mesh, seed=_SEED_X),))


# ═════════════════════════════════════════════════════════════════════════
# Shared reciprocity machinery — ONE spelling (Pattern 2)
# ═════════════════════════════════════════════════════════════════════════

def _reciprocity_residual(
    op: "LinearOperator", sn_mesh: SNMesh, x: TimedFullField, y: TimedFullField,
) -> "tuple[float, float]":
    r"""``(relative defect, |⟨Ax,y⟩_G|)`` for the G-adjoint identity.

    Returns the magnitude alongside the defect so every caller can check the
    row is not vacuously satisfied by both sides being zero.

    The inner product is :func:`~tests.sn._test_helpers.g_inner`, built from
    raw ``volumes`` / ``quad.weights`` / ``omega_dot_n`` — **structurally
    independent** of the production metric under test. Evaluating with the
    production metric would be a false green by construction: a wrong
    internal metric ``G'`` satisfies its own reciprocity trivially.
    """
    lhs = g_inner(op.apply(x), y, sn_mesh)
    rhs = g_inner(x, op.H.apply(y), sn_mesh)
    scale = max(abs(lhs), abs(rhs), 1e-300)
    return abs(lhs - rhs) / scale, abs(lhs)


def _drop_the_metric(self, y):
    """**M-10** — ``_AdjointOperator.apply`` with the metric removed.

    The honest body is
    :math:`(A^{*}y)_V = G_V^{+}\\odot \\mathrm{apply\\_transpose}(G_W \\odot y)`
    (``numerics/operator.py:1221-1227``); this returns the bare Euclidean
    transpose instead — which is EXACTLY what the operator already does today
    whenever its space is ``None`` (R2).
    """
    return self.inner.apply_transpose(y)


def _double_the_adjoint(self, y):
    """The NON-VACUITY mutation — a correct adjoint, scaled by two.

    Metric-agnostic on purpose: it moves the identity for every leaf,
    including the ones whose reciprocity commutes with the metric, so it can
    prove no G1.4 row is a dead gate.
    """
    return _TRUE_ADJOINT_APPLY(self, y) * 2.0


_TRUE_ADJOINT_APPLY = _operator_module._AdjointOperator.apply


def _assert_metric_is_constant(sn_mesh: SNMesh) -> float:
    """Precondition of the blindness leg: the metric IS one number.

    Reddens if the control fixture ever stops being flat — which would make
    :func:`test_a_globally_constant_metric_makes_reciprocity_blind` pass for
    a reason unrelated to the claim it makes.
    """
    weights = np.asarray(sn_mesh.quad.weights, dtype=float)
    volumes = np.asarray(sn_mesh.volumes, dtype=float)
    cosines = np.abs(np.asarray(sn_mesh.angular_trace.omega_dot_n, dtype=float))
    entries = np.concatenate([
        np.outer(weights, volumes).ravel(),
        (cosines * weights[None, :]).ravel(),
    ])
    unique = np.unique(np.round(entries, 12))
    if unique.size != 1:
        pytest.fail(
            f"the blindness control leg's metric is NOT a global constant: "
            f"{unique.size} distinct entries {unique[:6]} — the leg no longer "
            f"demonstrates config-blindness, because a varying metric is "
            f"exactly what the load-bearing leg uses."
        )
    return float(unique[0])


# ═════════════════════════════════════════════════════════════════════════
# G1.1 — declared spaces
# ═════════════════════════════════════════════════════════════════════════

@pytest.mark.parametrize("geometry", list(_GEOMETRIES), ids=list(_GEOMETRIES))
@pytest.mark.parametrize("leaf", _LEAVES)
def test_leaf_declares_both_function_spaces(leaf, geometry):
    r"""**G1.1**, value leg — a production leaf states BOTH ends of its arrow.

    GREEN today on the SN ladder (MEASURED: all 20 rows), so this ships as a
    **regression floor**, not a red gate: every solver-side construction
    threads ``sn_mesh.full_field_space``. Its job is to keep it that way
    while P1 moves the spaces from optional to mandatory — a phase that
    touches every leaf's constructor is exactly when a thread can be dropped.

    The red half of G1.1 lives in
    :func:`test_model_generic_leaf_declares_a_space` (**R1**, the
    model-generic construction) and
    :func:`test_leaf_space_annotation_is_not_optional` (R1's static face).
    """
    sn_mesh = _GEOMETRIES[geometry]()
    op = _leaf_set(sn_mesh)[leaf]
    for role, space in (("domain", op.domain), ("codomain", op.codomain)):
        if space is None:
            pytest.fail(
                f"{type(op).__name__}.{role} is None on the {geometry} "
                f"production construction — the leaf does not state its "
                f"arrow, so `.H` degrades to a Euclidean transpose (R2) and "
                f"`as_matrix` cannot derive a basis shape."
            )


# CS1 step 3b (2026-08-21): ``test_model_generic_leaf_declares_a_space``
# (4 strict-xfail rows, C/F × 2g/4g) was DELETED here. The C rows XPASSed
# the moment the ``from_mesh`` mesh-default chain gained the ``bulk_space``
# arm (strict ⟹ forced); the F rows could never XPASS (the body's bare
# ``from_solver_data(mat_xs=…)`` stays space-less under the ruled
# no-default-derivation) and were deleted in the same commit on the
# retired-mirror warrant: the production line they mirrored now threads
# ``space=``. The successor gate for all four rows is the positive floor
# ``tests/homogeneous/test_operator_spaces.py::``
# ``test_every_homogeneous_operator_reports_the_same_space``.


def _domain_annotation(leaf_cls: type) -> "tuple[str, str]":
    """``(owning class, raw return annotation)`` for the ``domain`` property.

    Read off the OWNING class in the MRO rather than the leaf, so the failure
    message points at the declaration a fix must edit. Annotations are
    strings under ``from __future__ import annotations``; they are compared
    textually on purpose — resolving the forward reference would require
    importing the ``TYPE_CHECKING``-only ``FunctionSpace`` at runtime and
    would tell us nothing extra.
    """
    for klass in leaf_cls.__mro__:
        prop = vars(klass).get("domain")
        if prop is None:
            continue
        fget = getattr(prop, "fget", prop)
        annotation = getattr(fget, "__annotations__", {}).get("return")
        return klass.__name__, str(annotation)
    return "<not found>", "<not found>"


@pytest.mark.parametrize("leaf", _LEAVES)
@_R1_XFAIL
def test_leaf_space_annotation_is_not_optional(leaf):
    r"""**G1.1 / R1**, static face — ``domain`` is not typed ``Optional``.

    RED today for all five leaves. MEASURED annotations::

        L  StreamingOperator       Optional['FunctionSpace']
        C  MultiplicationOperator  'FunctionSpace | None'
        S  ScatteringOperator      'FunctionSpace | None'
        F  FissionOperator         'FunctionSpace | None'
        B  SNBoundaryOperator      Optional['FunctionSpace']

    The value leg and the annotation leg are genuinely different claims and
    both are needed. A value that *happens* to be set is a runtime accident:
    the ``| None`` is what makes the composability guard *skippable*
    (``OperatorSum`` skips the domain check when either side is ``None`` —
    ``operator.py:582``) and what makes the ``.H`` metric branch
    *conditional*. P1 removes the Optional, which turns "the space is set"
    from a fact about one call site into a fact about the type.
    """
    sn_mesh = _sphere()
    owner, annotation = _domain_annotation(type(_leaf_set(sn_mesh)[leaf]))
    if "None" in annotation or "Optional" in annotation:
        pytest.fail(
            f"{owner}.domain is annotated {annotation} — an OPTIONAL space. "
            f"While it is optional the composability guard is skippable and "
            f"the `.H` metric application is conditional (R1/R2)."
        )


# ═════════════════════════════════════════════════════════════════════════
# G1.3 — the refusal is uniform, typed, and names the operator
# ═════════════════════════════════════════════════════════════════════════

_G13_ROWS = [
    pytest.param(
        leaf, geometry, carrier,
        marks=[_R6_XFAIL] if leaf == "B" else [],
        id=f"{leaf}-{geometry}-{carrier}",
    )
    for leaf in _LEAVES
    for geometry in _GEOMETRIES
    for carrier in ("alien", "coupled_field")
]


@pytest.mark.parametrize(("leaf", "geometry", "carrier"), _G13_ROWS)
def test_wrong_carrier_refusal_is_typed_and_names_the_operator(
    leaf, geometry, carrier,
):
    r"""**G1.3** — every leaf refuses a foreign carrier the SAME way.

    Three requirements, each measurable:

    1. the refusal is a ``TypeError`` — a *decision*, not an accident of the
       first attribute that happened to be missing;
    2. the message names the OPERATOR, so the reader learns which leaf in a
       composition rejected the input;
    3. the message names the expected carrier (``FullField``), so the reader
       learns what to hand it instead.

    **RED on ``B``, GREEN on the other four** — and shipping the sixteen
    green rows is the point, not padding: they are the **control leg** that
    makes the ``B`` red *attributable*. Without them a reviewer cannot tell
    "one leaf is non-conforming" from "the contract is aspirational".

    ``streaming.py:153`` is the shape the other four already have and the one
    ``B`` must reach — it names the method, the expected carrier, the carrier
    received, AND the remediation. The gate is satisfied by any dispatch
    spelling (``singledispatchmethod``, ``match``, a registry); #261 parks
    that ruling and no assertion here reads the mechanism.

    **Scope note (a deliberate narrowing of the written spec).** §2.1 asks
    for the message to name "the expected SPACE". MEASURED: not one leaf's
    message names a ``FunctionSpace`` — they all name CARRIERS. Requiring the
    space today would redden all five rows and destroy the control leg,
    turning an attributable single-leaf defect into an unattributable
    blanket one. The stronger wording belongs with P1, when the leaf's arrow
    is a declared, non-optional thing to name.

    **M-9**: today's state IS the mutation — this row's first run is its own
    mutation proof (§7).
    """
    sn_mesh = _GEOMETRIES[geometry]()
    op = _leaf_set(sn_mesh)[leaf]
    probe = _wrong_carrier(carrier, sn_mesh)

    try:
        op.apply(probe)  # type: ignore[arg-type]
    except TypeError as exc:
        message = str(exc)
        if type(op).__name__ not in message:
            pytest.fail(
                f"{type(op).__name__} refused a {carrier} carrier with a "
                f"TypeError that does not name the operator: {message!r} — a "
                f"reader cannot tell which leaf of a composition rejected the "
                f"input."
            )
        if "FullField" not in message:
            pytest.fail(
                f"{type(op).__name__}'s refusal does not name the expected "
                f"carrier: {message!r} — the reader is told what is wrong but "
                f"not what to hand it instead (contrast streaming.py:153)."
            )
    except Exception as exc:
        # Deliberately broad: the exception CLASS is the measurement. Naming
        # the classes we expect would hide exactly the ones we must catch.
        pytest.fail(
            f"{type(op).__name__} refused a {carrier} carrier with a raw "
            f"{type(exc).__name__}: {str(exc)[:120]!r}. The refusal must be a "
            f"typed TypeError naming the operator and the expected carrier "
            f"(R6 — boundary.py:343 reads `psi.interior` unguarded)."
        )
    else:
        pytest.fail(
            f"{type(op).__name__} ACCEPTED a {carrier} carrier — the leaf is "
            f"not monomorphic; it has no arrow at all."
        )


# ═════════════════════════════════════════════════════════════════════════
# G1.4 — `.H` is the G-adjoint
# ═════════════════════════════════════════════════════════════════════════

@pytest.mark.parametrize("geometry", list(_GEOMETRIES), ids=list(_GEOMETRIES))
@pytest.mark.parametrize("leaf", _LEAVES)
def test_hilbert_adjoint_reciprocity(leaf, geometry):
    r"""**G1.4** — :math:`\langle A\psi,\varphi\rangle_G = \langle\psi,A^\dagger
    \varphi\rangle_G` on a genuinely non-trivial metric.

    GREEN today (MEASURED worst residual 4.50e-14, contract 1e-12), so this
    is a **regression floor**. Its whole value is in its teeth, which are the
    two mutation gates below — and in the config, which is load-bearing:
    the bulk metric is :math:`V_{\rm cell} w_n` and the trace metric
    :math:`|\Omega\cdot\hat n| w_n`, so a *constant* metric would cancel from
    both sides and reduce this identity to a bare Euclidean-transpose check
    — blind to precisely the bug it exists to catch
    (``vv-principles`` Mode 12; ERR-067 is the catalogued instance of a
    degenerate metric putting an error class inside the measured functional's
    stabiliser).

    The non-vacuity guard is not decoration: on a non-fissile mixture ``F``
    is the zero operator and this identity degenerates to ``0 == 0``, which
    is why this file builds a fissile mixture (module docstring).
    """
    sn_mesh = _GEOMETRIES[geometry]()
    op = _leaf_set(sn_mesh)[leaf]
    x = _random_composite(sn_mesh, seed=_SEED_X)
    y = _random_composite(sn_mesh, seed=_SEED_Y)

    residual, magnitude = _reciprocity_residual(op, sn_mesh, x, y)
    if magnitude < 1e-8:
        pytest.fail(
            f"{type(op).__name__} on {geometry}: |<Ax,y>_G| = {magnitude:.3e} "
            f"— the operator is (near-)zero on this fixture, so reciprocity "
            f"holds vacuously. A channel the mixture should activate is off."
        )
    if residual > _RECIPROCITY_RTOL:
        pytest.fail(
            f"{type(op).__name__} on {geometry}: G-adjoint reciprocity "
            f"defect {residual:.3e} > {_RECIPROCITY_RTOL:.0e}. `.H` is not "
            f"the Hilbert adjoint of `apply` under the phase-space measure "
            f"(V_cell*w_n on the bulk, |Omega.n|*w_n on the trace)."
        )


@pytest.mark.parametrize("geometry", _CURVILINEAR)
def test_reciprocity_metric_is_load_bearing(geometry, monkeypatch):
    r"""**M-10**, first half — drop the metric and G1.4 REDs by O(1).

    Mutates ``_AdjointOperator.apply`` to return ``inner.apply_transpose(y)``
    bare, IN PROCESS via ``monkeypatch`` — never by editing a file, because
    this working tree carries uncommitted-by-policy state that a
    ``git checkout`` would destroy.

    MEASURED jumps against a clean ≤ 4.5e-14:

    ==========  ======  ======  ========
    geometry    L       S       F
    ==========  ======  ======  ========
    sphere      1.40    1.01    1.46
    cylinder    1.11    0.947   5.83e-2
    ==========  ======  ======  ========

    Two geometries on two different quadrature families (GL and product), so
    the result is a property of the metric rather than of one rule's weight
    distribution.

    The **control leg is inside this test**: the unmutated residual is
    re-measured first and must be under contract. Without it a
    still-broken baseline — off-tolerance both before and after — would
    *mimic* "caught", which is the ERR-067 closure trap in one line.

    Only ``{L, S, F}`` are asserted. ``C`` and ``B`` are metric-INVARIANT by
    algebra, not by fixture accident (module docstring, refutation #2);
    asserting a jump for them would be asserting a falsehood, and asserting
    silence would calcify. Their liveness is closed by
    :func:`test_reciprocity_row_is_non_vacuous`.
    """
    sn_mesh = _GEOMETRIES[geometry]()
    leaves = _leaf_set(sn_mesh)
    x = _random_composite(sn_mesh, seed=_SEED_X)
    y = _random_composite(sn_mesh, seed=_SEED_Y)

    clean = {
        name: _reciprocity_residual(leaves[name], sn_mesh, x, y)[0]
        for name in _METRIC_CONSTRAINED
    }
    off_contract = {n: r for n, r in clean.items() if r > _RECIPROCITY_RTOL}
    if off_contract:
        pytest.fail(
            f"CONTROL LEG BROKEN on {geometry}: the UNMUTATED reciprocity is "
            f"already off contract for {off_contract} — a mutation cannot be "
            f"credited with catching anything against a broken baseline "
            f"(the ERR-067 closure trap)."
        )

    monkeypatch.setattr(
        _operator_module._AdjointOperator, "apply", _drop_the_metric,
    )
    mutated = {
        name: _reciprocity_residual(leaves[name], sn_mesh, x, y)[0]
        for name in _METRIC_CONSTRAINED
    }
    silent = {n: r for n, r in mutated.items() if r < _MUTATION_FLOOR}
    if silent:
        pytest.fail(
            f"M-10 is SILENT on {geometry} for {silent} (floor "
            f"{_MUTATION_FLOOR:.0e}) — dropping the metric from "
            f"`_AdjointOperator.apply` does not move reciprocity, so this "
            f"config cannot distinguish the G-adjoint from the Euclidean "
            f"transpose. Check the metric still varies along BOTH the spatial "
            f"and the angular axis (a level_symmetric rule has constant "
            f"weights and leaves S and F blind)."
        )


def test_a_globally_constant_metric_makes_reciprocity_blind(monkeypatch):
    r"""**M-10**, second half — the CONFIG-BLINDNESS proof.

    On a metric that is one global constant, :math:`G = c\,\mathbb{1}` gives
    :math:`G^{-1}A^{\mathsf T}G = A^{\mathsf T}` identically, so dropping the
    metric changes **nothing** — MEASURED: every leaf stays under 4.3e-16,
    against 1.40 / 1.01 / 1.46 on the sphere.

    This is what makes the sphere/cylinder configuration *load-bearing*
    rather than decorative. Had G1.4 been written on a flat-metric fixture,
    it would have been a Euclidean-transpose check wearing a G-adjoint's
    name: green forever, catching nothing.

    It also **refutes the gate spec's phrasing**, which names the blind leg
    "the uniform-``h`` leg". MEASURED: a uniform-``h`` slab under
    ``gauss_legendre(4)`` is NOT blind — M-10 moves it 1.3e-1 (``L``),
    4.0e-1 (``S``), 2.7e-1 (``F``), because the *quadrature weights* vary
    even when ``V_cell`` does not. Blindness needs the metric constant along
    BOTH axes, which is why :func:`_flat_metric_slab` pins ``h = 1/√3``
    against a 2-point Gauss-Legendre rule and
    :func:`_assert_metric_is_constant` guards the property.
    """
    sn_mesh = _flat_metric_slab()
    constant = _assert_metric_is_constant(sn_mesh)
    leaves = _leaf_set(sn_mesh)
    x = _random_composite(sn_mesh, seed=_SEED_X)
    y = _random_composite(sn_mesh, seed=_SEED_Y)

    monkeypatch.setattr(
        _operator_module._AdjointOperator, "apply", _drop_the_metric,
    )
    for name in _LEAVES:
        residual, _ = _reciprocity_residual(leaves[name], sn_mesh, x, y)
        if residual > _RECIPROCITY_RTOL:
            pytest.fail(
                f"{name} on the flat-metric slab (G == {constant:.17g} "
                f"everywhere) moved {residual:.3e} when the metric was "
                f"dropped — but a constant metric cancels from both sides of "
                f"reciprocity exactly. Either the fixture is no longer flat "
                f"or `.H` does something beyond G^-1 A^T G."
            )


@pytest.mark.parametrize("leaf", _LEAVES)
def test_reciprocity_row_is_non_vacuous(leaf, monkeypatch):
    r"""Every G1.4 row can RED — including the metric-invariant ones.

    M-10 constrains only ``{L, S, F}``; ``C`` and ``B`` commute with the
    metric by algebra (module docstring, refutation #2). Left there, two of
    the twenty G1.4 rows would be *exercised-but-unconstrained*
    (``vv-principles`` Mode 10) with nothing said about it — a gate that runs
    and can never fail.

    So this ships a second, metric-agnostic mutation: scale the adjoint by
    two. It reddens EVERY leaf at exactly 0.5 relative (MEASURED), which
    proves each row still pins ``apply_transpose``'s structure — the ``C``
    and ``B`` rows are honest gates on the transpose, they simply cannot see
    the metric, and now that is a measured statement rather than an
    assumption.

    ``B``'s metric weighting itself is pinned elsewhere: the L11
    wrong-metric control in
    ``tests/sn/operators/test_g_adjoint_reciprocity.py`` drops
    :math:`|\Omega\cdot\hat n|` from the reference metric on the composite
    ``(L + C - B)``. Not duplicated here.
    """
    sn_mesh = _sphere()
    op = _leaf_set(sn_mesh)[leaf]
    x = _random_composite(sn_mesh, seed=_SEED_X)
    y = _random_composite(sn_mesh, seed=_SEED_Y)

    monkeypatch.setattr(
        _operator_module._AdjointOperator, "apply", _double_the_adjoint,
    )
    residual, _ = _reciprocity_residual(op, sn_mesh, x, y)
    if residual < _MUTATION_FLOOR:
        pytest.fail(
            f"{type(op).__name__}: doubling `.H` moved reciprocity by only "
            f"{residual:.3e} (floor {_MUTATION_FLOOR:.0e}) — this G1.4 row "
            f"cannot red, so it is a dead gate whatever it reports."
        )


# ═════════════════════════════════════════════════════════════════════════
# G1.5 — an anonymous leaf is unrepresentable
# ═════════════════════════════════════════════════════════════════════════

@pytest.mark.parametrize("leaf", _ANONYMOUS_CAPABLE)
@_R2_XFAIL
def test_leaf_without_a_space_refuses_construction(leaf):
    r"""**G1.5 / R1 / R2** — building a leaf with NO space must RAISE.

    RED today for all three of ``C``, ``S``, ``F``: each has an optional
    space parameter defaulting to ``None`` (MEASURED signatures) and each
    constructs happily without one. The consequence is R2 and it is silent:
    ``_AdjointOperator.apply`` (``numerics/operator.py:1221-1227``) applies
    the metric only ``if inner_codomain is not None``, so an anonymous leaf's
    ``.H`` is a **bare Euclidean transpose** wearing the Hilbert adjoint's
    name. MEASURED on the homogeneous ``F``: ``F.H.apply(φ)`` is
    ``array_equal`` to ``F.apply_transpose(φ)``.

    Guarding a caller cannot fix this — the degradation happens where no
    caller looks. P1 item 3 makes it **unrepresentable** instead: with the
    space mandatory, ``.H`` can never see a ``None``.

    Why the body is shaped this way (the xfail-strict false-positive, lesson
    L4). A ``strict=True`` xfail is satisfied by failing for ANY reason, so a
    test that fails on a stale fixture is a FALSE xfail — green suite, wrong
    reason. **This test was written that way first and caught itself**: the
    ``S`` row failed on a ``ValueError`` out of ``np.einsum`` (an anonymous
    ``ScatteringOperator``'s ``.H`` will not take the meshless ``(ng, 1)``
    probe ``C``/``F`` accept), so it "xfailed" while asserting nothing about
    R2. The degradation demonstration is therefore BEST-EFFORT and its
    outcome — including a failure to run — is reported as *evidence text*,
    never as the test's verdict. The ONLY route to passing is the constructor
    raising; the ONLY route to failing is the ``pytest.fail`` below. When P1
    lands, the ``except`` returns, this XPASSes, and strict mode makes that a
    hard failure that forces the marker's removal.
    """
    mixture = get_mixture("A", "2g")
    mat_xs = MaterialMesh.from_materials({0: mixture}).material_xs_field()
    builders = {
        "C": lambda: MultiplicationOperator(
            coefficient=mat_xs.total_cross_section_field,
        ),
        "S": lambda: ScatteringOperator.from_solver_data(
            mat_xs=mat_xs, quadrature=Quadrature.gauss_legendre(n_ordinates=4),
            scattering_order=0,
        ),
        "F": lambda: FissionOperator.from_solver_data(mat_xs=mat_xs),
    }
    try:
        op = builders[leaf]()
    except (TypeError, ValueError):
        return  # P1 has landed: an anonymous leaf is unrepresentable.

    probe = np.asarray(mixture.chi, dtype=float).reshape(mixture.ng, 1) + 0.3
    degraded = "not adjointable"
    if op.is_adjointable:
        try:
            adjoint_image = np.asarray(op.H.apply(probe))
            transpose_image = np.asarray(op.apply_transpose(probe))
            degraded = str(bool(np.array_equal(adjoint_image, transpose_image)))
        except Exception as exc:
            # Deliberately broad: this is EVIDENCE, never the verdict. Any
            # failure here must not become the xfail's reason (lesson L4).
            degraded = (
                f"undemonstrated on this probe ({type(exc).__name__}); the "
                f"claim under test is the CONSTRUCTION, not the probe"
            )
    pytest.fail(
        f"{type(op).__name__} CONSTRUCTED with no function space "
        f"(domain={op.domain!r}); its `.H` is bit-identical to the bare "
        f"Euclidean transpose: {degraded}. The metric-less adjoint is "
        f"representable (R1/R2)."
    )


@pytest.mark.parametrize("leaf", ["L", "B"])
def test_mesh_derived_leaves_carry_no_anonymous_construction_surface(leaf):
    r"""``L`` and ``B`` cannot be anonymous — which is why they are excluded.

    GREEN today, and it is the **scope control** for
    :func:`test_leaf_without_a_space_refuses_construction`. Without it, that
    test's ``{C, S, F}`` parametrisation reads as an arbitrary list; with it,
    the list is *complete* — the other two take an ``SNMesh`` and derive
    their space from it, so there is no anonymous surface to refuse.

    It has real teeth: adding an optional ``space=None`` parameter to
    ``StreamingOperator``/``SNBoundaryOperator`` (an easy thing to do while
    "generalising" a constructor during P1/P2) opens the R2 degradation on
    two more leaves and REDs this row immediately.
    """
    constructor = {"L": StreamingOperator, "B": SNBoundaryOperator}[leaf]
    parameters = inspect.signature(constructor.__init__).parameters
    optional_space = [
        name for name, param in parameters.items()
        if param.default is None
        and any(token in name for token in ("space", "domain", "codomain"))
    ]
    if optional_space:
        pytest.fail(
            f"{constructor.__name__}.__init__ grew an optional space "
            f"parameter {optional_space} — a leaf that derives its space from "
            f"its SNMesh now admits an ANONYMOUS construction, which is the "
            f"R2 `.H`-degrades-silently surface (see "
            f"test_leaf_without_a_space_refuses_construction)."
        )
