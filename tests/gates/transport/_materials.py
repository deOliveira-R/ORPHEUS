r"""Material helpers shared by the transport gates and every method's gates (layer DAG Phase 4, class Fx).

Re-homed from ``tests/gates/sn/_test_helpers.py`` (``placeholder_materials``, ``material_xs_from_raw``,
``CURVED_REFLECTIVE``, ``reflective_role``), ``tests/gates/sn/architecture/_config.py`` (``anisotropic_mixture``) and
``tests/gates/cp/test_verification.py`` (``_make_mixture_with_n2n``, renamed ``mixture_with_n2n``): their bodies read
only data, mesh, geometry and transport, so their home is the lowest gate layer whose vocabulary suffices. The old
modules import them back (re-export), so their method-side consumers do not change.
"""
from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import pytest
from scipy.sparse import csr_matrix

from orpheus.data.macro_xs.mixture import Mixture
from orpheus.data.macro_xs.synthetic import get_xs
from orpheus.geometry.deck import MirrorDeck
from orpheus.transport.boundary_conditions import WallCondition

if TYPE_CHECKING:
    from orpheus.transport.mesh.material_xs_field import MaterialXSField

def placeholder_materials(
    ng: int = 1, mat_ids: tuple[int, ...] = (0,),
) -> dict:
    """Build a placeholder ``{mat_id: Mixture}`` dict for SNProblem tests.

    Parameters
    ----------
    ng
        Number of energy groups.  All Mixtures will report this value
        via :attr:`Mixture.ng`.
    mat_ids
        Material ids to include in the dict.  Default ``(0,)`` covers
        the common case where the mesh's ``mat_ids`` / ``mat_map`` is
        zeros.

    Returns
    -------
    dict[int, Mixture]
        Each entry has ``SigT = ones(ng)`` and all other cross sections
        zero.  Suitable for SNProblem tests that don't compute physical
        quantities from the materials.
    """
    from orpheus.data.macro_xs.mixture import Mixture
    z = np.zeros(ng)
    z_mat = csr_matrix(np.zeros((ng, ng)))
    return {
        int(mid): Mixture(
            SigC=z.copy(),
            SigL=z.copy(),
            SigF=z.copy(),
            SigP=z.copy(),
            SigT=np.ones(ng),
            SigS=[z_mat],
            Sig2=[z_mat],
            chi=z.copy(),
        )
        for mid in mat_ids
    }


def _as_stack(channel) -> "list":
    """A single matrix is the P0 stack; a list/tuple is the stack as given."""
    return list(channel) if isinstance(channel, (list, tuple)) else [channel]


def material_xs_from_raw(
    *,
    sig_s: "dict[int, list[np.ndarray]]",
    sig2: "dict[int, np.ndarray | list[np.ndarray]] | None" = None,
    cells_by_mat: "dict[int, tuple[np.ndarray, np.ndarray]]",
    ng: int,
    nx: int,
    ny: int = 1,
) -> "MaterialXSField":
    """Build a REAL ``MaterialXSField`` from raw per-material scattering data.

    The production-path replacement for the retired
    ``MaterialXSField._synthetic_for_tests`` (campaign 1 CS4b, Q7 ruling —
    option 3): real :class:`Mixture` objects carry the given Legendre
    lists and (n,2n) matrices — or (n,2n) Legendre STACKS, a list per
    material, since #426 step 2 — a real ``Mesh2D`` paints ``cells_by_mat``
    into its ``mat_map``, and the field is the production
    hub's own ``MaterialMesh.mat_xs`` — so the per-material dispatch tests
    exercise the true admission + lazy dense-cache path (including the
    EE-4 read-only cache freeze) instead of pre-populated caches on a
    duck-typed mesh stub, which the space-carrying accessors of CS4b
    cannot serve.

    The grid is always rank-2 ``(nx, ny)`` (the stub's convention — the
    consuming fixtures shape their arrays accordingly), unit-pitch
    Cartesian. Non-scattering channels are zero and ``SigT = 1``: these
    fixtures assert scattering dispatch, never a physical balance.
    Every cell must be painted — a real carrier has no unassigned cells.
    """
    from orpheus.mesh import Mesh2D
    from orpheus.transport.mesh.material_mesh import MaterialMesh
    from orpheus.data.macro_xs.mixture import Mixture

    mat_map = np.full((nx, ny), -1, dtype=int)
    for mid, (ix, iy) in cells_by_mat.items():
        mat_map[np.asarray(ix), np.asarray(iy)] = mid
    if (mat_map < 0).any():
        raise ValueError(
            "material_xs_from_raw: cells_by_mat must paint every cell of "
            f"the ({nx}, {ny}) grid — a real carrier has no unassigned "
            "cells."
        )

    z = np.zeros(ng)
    materials = {
        int(mid): Mixture(
            SigC=z.copy(),
            SigL=z.copy(),
            SigF=z.copy(),
            SigP=z.copy(),
            SigT=np.ones(ng),
            SigS=[csr_matrix(np.asarray(s)) for s in mats],
            Sig2=[
                csr_matrix(np.asarray(block))
                for block in _as_stack(
                    sig2[mid] if sig2 is not None and mid in sig2
                    else np.zeros((ng, ng))
                )
            ],
            chi=z.copy(),
        )
        for mid, mats in sig_s.items()
    }
    mesh = Mesh2D(
        edges_x=np.arange(nx + 1, dtype=float),
        edges_y=np.arange(ny + 1, dtype=float),
        mat_map=mat_map,
        face_roles={face: MirrorDeck() for face in ("xmin", "xmax", "ymin", "ymax")},
    )
    return MaterialMesh(mesh, materials).mat_xs


#: The unit specular wall a CURVED body's face takes where a Cartesian body takes a mirror deck (rule R2: a curved face carries a
#: tag, never a deck). Tag ``0`` is the default role of the shared curvilinear builders below, and this is its table: a test that
#: builds a curved body with the default role declares ``boundary_conditions=CURVED_REFLECTIVE`` at its solve.
CURVED_REFLECTIVE = {0: WallCondition.specular(1.0)}


def reflective_role(coord, tag: int = 0):
    """A mirror deck on a Cartesian body, the tag ``tag`` (default ``0``, see :data:`CURVED_REFLECTIVE`) on a curved one.

    A test whose own table maps another tag to the unit specular wall passes that tag.
    """
    from orpheus.geometry import CoordSystem

    return MirrorDeck() if coord is CoordSystem.CARTESIAN else tag


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
    (``macro_xs/synthetic.py:56``) hardcodes ``SigL = 0`` and defaults ``Sig2`` to an
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
        SigS=legendre, Sig2=[n2n],
        chi=zero.copy() if chi is None else np.asarray(chi, dtype=float),
    )


def mixture_with_n2n(ng_key: str = "2g") -> Mixture:
    """Region A with nonzero (n,2n) cross section.

    Creates a material where group 0 has a (n,2n) reaction that
    produces 2 neutrons into the same group.  This tests that the
    solver correctly accounts for the extra neutron production.

    The Sig2 matrix encodes the transfer: sig2[g'->g].
    """
    xs = get_xs("A", ng_key)
    ng = len(xs["sig_t"])
    sig_s = xs["sig_s"].copy()

    # (n,2n) from group 0 to group 0 (fast self-multiplication)
    sig2 = np.zeros((ng, ng))
    sig2[0, 0] = 0.01  # small but nonzero

    # sig_t must include the (n,2n) reaction cross section
    # The (n,2n) XS out of group g is sum over g' of sig2[g, g']
    sig2_out = sig2.sum(axis=1)
    sig_t = xs["sig_c"] + xs["sig_f"] + sig_s.sum(axis=1) + sig2_out

    # Synthetic XS for the (n,2n) consistency test — no physical grid (Phase E).
    sig_s_list = [csr_matrix(sig_s)]
    return Mixture(
        SigC=xs["sig_c"].copy(),
        SigL=np.zeros(ng),
        SigF=xs["sig_f"].copy(),
        SigP=(xs["nu"] * xs["sig_f"]).copy(),
        SigT=sig_t,
        SigS=sig_s_list,
        Sig2=[csr_matrix(sig2)],
        chi=xs["chi"].copy(),
    )
