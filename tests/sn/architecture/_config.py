r"""The campaign's MANDATORY CONFIGURATION — one home, per gate reuse.

``.claude/plans/campaign_verification_plan.md`` §8 is a *table of levers*,
each of which breaks one specific blindness.  Re-spelling that table inside
every gate file would be a Pattern-2 twin **inside the gate itself** — the
exact failure O-11 names for the AC-c fingerprint helper.  So the fixtures
live here, once, with the blindness each lever breaks stated on the line
that pulls it.

Why each lever (§8, and the anti-pattern it defeats):

=========================  ==========================================================
lever                      the blindness it breaks
=========================  ==========================================================
**≥2G**                    ``k = νΣf/Σa`` is flux-shape-independent, so a 1-group
                           value claim cannot see a spatial/angular/scattering error
                           (``vv-principles`` anti-pattern #3).
**asymmetric** ``SigS``    a group-index or block-index transpose is *invisible* on a
                           symmetric operator — Mode 12 (the error class sits inside
                           the measured functional's invariance group).
**non-uniform / non-**     an ``x↔y`` swap; and a *constant* metric cancels from both
**square mesh**            sides of a reciprocity identity, degrading it to a bare
                           Euclidean-transpose check.
**heterogeneous**          flat flux nulls redistribution and every
(≥2 regions)               spatial-distribution bug.
**P1 anisotropic,**        ``make_mixture`` hardcodes ``SigL = 0`` and defaults
``Mixture`` built          ``Sig2 = 0``; an all-isotropic fixture nulls ``(𝕀 − P_iso)``
DIRECTLY                   — the exact term the σ_r gate measures (lessons L1).
**mixed vacuum +**         a vacuum-only probe nulls the reflection-partner trace rows
**reflective BC**          that ``B`` acts on.
``level_symmetric``        an axis-aligned ``product`` quadrature gives each face one
**quadrature**             octant — the ERR-056 degeneracy that hides a shared-face
                           Gauss-Seidel bug.
**fixed-seed random,**     flat ψ nulls the streaming coupling; a zero trace nulls
**non-flat state,**        ``B``.
**bulk AND trace**
=========================  ==========================================================

Helper modules are NOT collected by pytest, so a bare ``assert`` here would
be stripped by the canonical ``python -O`` invocation and assert *nothing*
(``vv-principles`` Mode 8).  Every check below is therefore ``pytest.fail``
or ``np.testing.assert_*``.
"""
from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import pytest
from scipy.sparse import csr_matrix

from orpheus.data.macro_xs.mixture import Mixture
from orpheus.geometry import BC, CoordSystem, Mesh1D
from orpheus.geometry.mesh import Mesh2D
from orpheus.numerics.coupled_system import CoupledField, CoupledOperator
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.coupled_system import build_within_group_system
from orpheus.sn.mesh.augmented_mesh import SNMesh
from orpheus.transport.full_field import FullField
from orpheus.transport.operators.multiplication_operator import (
    MultiplicationOperator,
)

if TYPE_CHECKING:
    from numpy.typing import NDArray

    from orpheus.numerics.operator import LinearOperator
    from orpheus.sn.coupled_system import WithinGroupSystem
    from orpheus.sn.operators.streaming import StreamingCollisionOperator
    from orpheus.transport.operators.scattering import ScatteringOperator


# ── materials ────────────────────────────────────────────────────────────

def anisotropic_mixture(
    sig_t: "list[float]",
    sig_s0: "list[list[float]]",
    sig_s1: "list[list[float]] | None" = None,
    *,
    sig_f: "list[float] | None" = None,
    chi: "list[float] | None" = None,
    nu: float = 2.6,
    sig_l: "list[float] | None" = None,
    sig_2: "list[list[float]] | None" = None,
) -> Mixture:
    r"""A ``Mixture`` built DIRECTLY — never through ``make_mixture``.

    **The campaign's one mixture builder.** ``make_mixture``
    (``xs_library.py:56``) hardcodes ``SigL = 0`` and defaults ``Sig2`` to an
    all-zero sparse matrix — and the SHIPPED tables pass no ``sig_s1``, so
    every channel the shipped pairs null goes *vacuously green* in a
    mutation catcher (lessons L1). (Wording corrected at CS4a-R QA-F7:
    ``make_mixture`` itself DOES accept a :math:`P_1` channel via
    ``sig_s1=`` — the nulling is the shipped DATA's, not the signature's.) Each optional
    argument here exists to un-null one of them, and each is off by default
    so a gate activates only what it actually constrains.

    What each argument activates
    ----------------------------
    ``sig_s0`` / ``sig_s1``
        Keep them **asymmetric** in the group index
        (``s0[0,1] != s0[1,0]``) so a group-axis transpose inside
        :math:`S` / :math:`S^\dagger` is observable rather than annihilated
        (Mode 12), and pass ``sig_s1`` so the :math:`P_1` moment blocks are
        built at all — an all-isotropic fixture nulls
        :math:`(\mathbb{1} - P_{\text{iso}})`, the single term the σ_r
        gates exist to measure.
    ``sig_f`` / ``chi`` / ``nu``
        Make the mixture **producing** (``SigP = ν·Σ_f > 0``), so
        ``F.apply(ψ) != 0`` and an ``F`` reciprocity row is not the
        tautology ``0 == 0``. ``chi`` must be a probability simplex for a
        producing mixture (``Mixture.__post_init__`` enforces it).
        **Omit both** for a within-group fixture: the within-group system
        carries no fission piece (fission enters as ``q_ext`` at the outer),
        so the default is non-fissile.
    ``sig_2``
        Non-zero and **asymmetric** exercises :math:`S`'s (n,2n) channel —
        a separate :math:`\ell = 0` term with a factor 2, outside the
        Legendre fold.
    ``sig_l``
        The (n,α) channel, non-zero on purpose where a gate needs it.
    """
    ng = len(sig_t)
    zero = np.zeros(ng)
    s0 = np.asarray(sig_s0, dtype=float)
    legendre = [csr_matrix(s0)]
    if sig_s1 is not None:
        legendre.append(csr_matrix(np.asarray(sig_s1, dtype=float)))
    fission = zero.copy() if sig_f is None else np.asarray(sig_f, dtype=float)
    n_alpha = zero.copy() if sig_l is None else np.asarray(sig_l, dtype=float)
    n2n = (
        csr_matrix((ng, ng)) if sig_2 is None
        else csr_matrix(np.asarray(sig_2, dtype=float))
    )
    if (sig_f is None) != (chi is None):
        pytest.fail(
            "sig_f and chi go together: a producing mixture needs an "
            "emission spectrum, and a non-producing one must have a null "
            "chi (Mixture.__post_init__ enforces the simplex/null law)."
        )
    return Mixture(
        SigC=np.asarray(sig_t, dtype=float) - s0.sum(axis=1) - fission - n_alpha,
        SigL=n_alpha,
        SigF=fission,
        SigP=nu * fission,
        SigT=np.asarray(sig_t, dtype=float),
        SigS=legendre, Sig2=n2n,
        chi=zero.copy() if chi is None else np.asarray(chi, dtype=float),
    )


#: Two heterogeneous 2G materials with ASYMMETRIC ``SigS`` in BOTH Legendre
#: orders (``s0[0,1] = 0.25 ≠ 0.02 = s0[1,0]``) — so a group-index transpose
#: is observable (Mode 12), and P1 present so ``(𝕀 − P_iso) ≠ 0``.
_REGION_A = ([1.2, 2.0], [[0.30, 0.25], [0.02, 0.90]], [[0.08, 0.05], [0.0, 0.20]])
_REGION_B = ([0.9, 2.6], [[0.22, 0.40], [0.01, 1.10]], [[0.05, 0.09], [0.0, 0.30]])


def _two_region_materials() -> dict[int, Mixture]:
    return {
        0: anisotropic_mixture(*_REGION_A),
        1: anisotropic_mixture(*_REGION_B),
    }


# ── meshes ───────────────────────────────────────────────────────────────

def cart2d_seedless() -> SNMesh:
    r"""2-D Cartesian, NON-SQUARE, heterogeneous 2G, P1, MIXED BC, S4 LS.

    **This is the only geometry on which R7 is observable.**  The boundary
    Gauss-Seidel arm of :func:`~orpheus.sn.solver._select_si_splitting` fires
    only on ``is_cartesian and not is_1d``; on a slab it falls back to Jacobi
    and the twin is *invisible*.  (The campaign's first probe of R7 used a
    slab and got the wrong answer — see
    :func:`test_a_slab_hides_r7_the_documented_trap`.)

    ``nx != ny`` is the ``x↔y``-swap catcher; ``level_symmetric`` avoids
    both the ``mu_y == 0`` GL rank mismatch and the ERR-056 axis-aligned
    degeneracy.
    """
    nx, ny = 4, 5
    mat_map = np.zeros((nx, ny), dtype=int)
    mat_map[nx // 2:, :] = 1                      # heterogeneous: 2 regions
    mesh = Mesh2D(
        edges_x=np.linspace(0.0, 2.0, nx + 1),
        edges_y=np.linspace(0.0, 3.0, ny + 1),
        mat_map=mat_map,
        coord=CoordSystem.CARTESIAN,
        bc_xmin=BC("reflective"), bc_xmax=BC("vacuum"),
        bc_ymin=BC("reflective"), bc_ymax=BC("vacuum"),
    )
    return SNMesh(
        mesh, Quadrature.level_symmetric(sn_order=4), _two_region_materials(),
    )


def sphere_carrying() -> SNMesh:
    r"""1-D spherical, heterogeneous 2G, P1, mixed BC — the CARRYING arm.

    A seed-carrying mesh (System B exists), so ``build_within_group_system``
    returns the 2×2 coupled grid on all three of ``loss`` /
    ``implicit_operator`` / ``explicit_gains``.  This is the arm whose shape
    the campaign's P3 makes the seedless arm match.
    """
    n = 6
    mesh = Mesh1D(
        edges=np.linspace(0.0, 3.0, n + 1),
        mat_ids=np.array([0, 0, 0, 1, 1, 1]),
        coord=CoordSystem.SPHERICAL,
        bc_left=BC("reflective"), bc_right=BC("vacuum"),
    )
    return SNMesh(
        mesh, Quadrature.gauss_legendre(n_ordinates=8), _two_region_materials(),
    )


def slab_seedless() -> SNMesh:
    """1-D Cartesian slab — the R7 **trap** control (G-S falls back to Jacobi)."""
    n = 8
    mesh = Mesh1D(
        edges=np.linspace(0.0, 4.0, n + 1),
        mat_ids=np.array([0, 0, 0, 0, 1, 1, 1, 1]),
        coord=CoordSystem.CARTESIAN,
        bc_left=BC("reflective"), bc_right=BC("vacuum"),
    )
    return SNMesh(
        mesh, Quadrature.gauss_legendre(n_ordinates=8), _two_region_materials(),
    )


def isotropic_slab(*, c: float = 0.9, sig_t: float = 1.0, n: int = 40) -> SNMesh:
    r"""Homogeneous 2G slab, **P0 self-scatter ONLY**, VACUUM both faces.

    The σ_r-fold fixture.  Every channel except within-group :math:`P_0` is
    zero, so the fold's error is exactly :math:`-\Sigma_{s0}(\mathbb{1} -
    P_{\text{iso}})` with nothing else contaminating it — and vacuum BC
    makes ``B ψ ≡ 0`` (MEASURED: ``|Bψ|∞ = 0.0``), so the Mode-9 control
    leg isolates the projection mechanism and nothing else.

    ``c = 0.9`` is the #215 divergence regime (``ρ = c/(1−c) = 9`` in the
    infinite medium; ≈ 6.91 measured with leakage).
    """
    mesh = Mesh1D(
        edges=np.linspace(0.0, 4.0, n + 1),
        mat_ids=np.zeros(n, dtype=int),
        coord=CoordSystem.CARTESIAN,
        bc_left=BC("vacuum"), bc_right=BC("vacuum"),
    )
    scatter = [[c * sig_t, 0.0], [0.0, c * sig_t]]
    mats = {0: anisotropic_mixture([sig_t, sig_t], scatter)}
    return SNMesh(mesh, Quadrature.gauss_legendre(n_ordinates=8), mats)


# ── the posed record ─────────────────────────────────────────────────────

def record_for(sn_mesh: SNMesh, *, scattering_order: int = 1) -> "WithinGroupSystem":
    """The posed within-group record — the ONE construction site."""
    return build_within_group_system(
        sn_mesh, sn_mesh.material_xs_field(),
        scattering_order=scattering_order,
    )


# ── probe states ─────────────────────────────────────────────────────────

def random_state(
    record: "WithinGroupSystem", *, seed: int, angularly_flat: bool = False,
) -> CoupledField:
    r"""A fixed-seed random state on the record's OWN carrier.

    Carrier-generic by construction: built from ``record.space.zeros()`` and
    filled through ``from_flat``, so the seedless (1×1) and carrying (2×2)
    arms take the same body and no gate has to know which arm it is on.

    ``angularly_flat=True`` broadcasts one ordinate's spatial pattern across
    every ordinate, giving :math:`P_{\text{iso}}\psi = \psi` — the **Mode-9
    degenerate subspace**, shipped only as a control leg (§3.2).
    """
    template = record.space.zeros()
    rng = np.random.default_rng([seed, 7])
    state = CoupledField.from_flat(
        rng.standard_normal(template.to_flat().size), template,
    )
    if not angularly_flat:
        return state
    values = system_a(state).interior.values       # (N, ng, *spatial)
    values[...] = np.broadcast_to(values[:1], values.shape)
    return state


def seedless_implicit(
    record: "WithinGroupSystem",
) -> "StreamingCollisionOperator":
    r"""``record.implicit_operator``, narrowed to the seedless arm's type.

    **This helper exists because of a defect, and it retires with that
    defect.**  ``implicit_operator`` is typed ``CoupledOperator |
    StreamingCollisionOperator`` — the union IS the shape asymmetry P3
    removes (the seedless arm is a bare operator where the carrying arm is a
    1-block grid), and it is why ``record.implicit_operator - anything``
    does not type-check.  When P3 makes both arms a ``CoupledOperator``, the
    union collapses and every call here becomes unnecessary.
    """
    from orpheus.sn.operators.streaming import StreamingCollisionOperator

    implicit = record.implicit_operator
    if not isinstance(implicit, StreamingCollisionOperator):
        pytest.fail(
            f"expected the SEEDLESS arm's bare (L+C); got "
            f"{type(implicit).__name__} — this fixture is carrying, or the "
            f"seedless arm's shape changed (P3 landing?)."
        )
    return implicit


def scattering_gain(record: "WithinGroupSystem") -> "ScatteringOperator":
    r"""``record.explicit_gains[0]``, narrowed to the scattering operator.

    **Also a defect marker.**  ``explicit_gains: tuple[LinearOperator, ...]``
    erases every role: the ``(S, N2N, B_a)`` convention with ``B_a`` LAST is a
    *positional* claim the type system cannot state, which is why the driver
    re-asserts it at runtime (``sn/solver.py:856`` raises a ``TypeError`` if
    gain 0 is not a ``ScatteringOperator``).  A tuple of roles is the thing
    P3 replaces with a gain GRID, after which the role is the slot.
    """
    from orpheus.transport.operators.scattering import ScatteringOperator

    gain = record.explicit_gains[0]
    if not isinstance(gain, ScatteringOperator):
        pytest.fail(
            f"the seedless record's FIRST gain must be the ScatteringOperator "
            f"(the builder's positional (S, N2N, B_a) convention); got "
            f"{type(gain).__name__} — the gain order changed and every "
            f"positional consumer of explicit_gains is now wrong."
        )
    return gain


def system_a(state: CoupledField) -> FullField:
    """The System-A member of a coupled state, narrowed to its real type.

    ``CoupledField.systems`` is typed by the ``SystemField`` Protocol — the
    member's ``interior``/``boundary`` composite structure is a fact about
    System A specifically, not about every system.  Narrow once, here, so no
    gate reaches through an un-narrowed member.
    """
    member = state.systems[0]
    if not isinstance(member, FullField):
        pytest.fail(
            f"System A must be a FullField (bulk ⊕ trace); got "
            f"{type(member).__name__} — the coupled space's member order "
            f"changed and every gate reading .systems[0] is now wrong."
        )
    return member


# ── the splitting law ────────────────────────────────────────────────────

def split_image(
    record: "WithinGroupSystem",
    state: CoupledField,
    *,
    implicit: "LinearOperator | None" = None,
    gains: "tuple[LinearOperator, ...] | None" = None,
) -> "NDArray":
    r"""``(M − Σ Nᵢ) x``, flattened and aligned with ``A x``.

    The carrier bridge lives here and nowhere else: on the carrying arm
    ``M`` and every gain are ``CoupledOperator``\ s that consume the coupled
    state directly; on the seedless arm they are System-A operators, so the
    state's single member is unwrapped.  **That asymmetry IS the shape
    defect P3 removes** — when it goes, this branch goes with it.
    """
    implicit = record.implicit_operator if implicit is None else implicit
    gains = record.explicit_gains if gains is None else gains
    if isinstance(implicit, CoupledOperator):
        image = implicit.apply(state).to_flat()
        for gain in gains:
            image = image - gain.apply(state).to_flat()
        return image
    member = system_a(state)
    image = implicit.apply(member).to_flat()
    for gain in gains:
        image = image - gain.apply(member).to_flat()
    return image


def reconstruction_residual(
    record: "WithinGroupSystem",
    state: CoupledField,
    *,
    implicit: "LinearOperator | None" = None,
    gains: "tuple[LinearOperator, ...] | None" = None,
) -> float:
    r"""``‖A x − (M − ΣN) x‖∞ / ‖A x‖∞`` — the splitting law's defect."""
    loss_image = record.loss.apply(state).to_flat()
    defect = loss_image - split_image(
        record, state, implicit=implicit, gains=gains,
    )
    denominator = float(np.max(np.abs(loss_image)))
    if denominator == 0.0:
        pytest.fail("‖A x‖∞ == 0 — the probe state does not excite A")
    return float(np.max(np.abs(defect))) / denominator


def sigma_s0_times_identity(
    sn_mesh: SNMesh, scattering: "ScatteringOperator",
) -> "MultiplicationOperator":
    r"""The σ_r fold's WRONG operator: :math:`\Sigma_{s0}^{g\to g}\,\mathbb{1}`.

    Diagonal in angle.  The right one is
    :math:`\Sigma_{s0}\,P_{\text{iso}}` (what ``S.foldable_part()`` realizes).
    The two agree **exactly** on an angularly-flat flux and nowhere else —
    which is the whole of #215 (46–56 % silent flux error) in one sentence.
    """
    per_material = scattering.foldable_sigma()          # {mid: (ng,)}
    mat_ids = np.asarray(sn_mesh.mat_map).reshape(sn_mesh.spatial_shape)
    field = np.zeros((sn_mesh.ng, *sn_mesh.spatial_shape))
    for mid, sigma in per_material.items():
        mask = mat_ids == mid
        for group in range(sn_mesh.ng):
            field[group][mask] = sigma[group]
    return MultiplicationOperator.from_mesh(field, sn_mesh)
