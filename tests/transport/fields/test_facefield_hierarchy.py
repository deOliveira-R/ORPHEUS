r"""L0 — the codim-1 :class:`FaceField` hierarchy (C2b of the FaceField carve).

C2b introduced :class:`~orpheus.transport.fields._bases.FaceField` — the shared
codim-1 parent that owns the flat-buffer / :class:`FaceLayout` discipline ONCE —
and re-parented the two codim-1 loci under it as **siblings**:

* :class:`BoundaryField` (``FaceField[str]``) — the SPATIAL faces;
* :class:`StartingDirectionField` (``FaceField[tuple[int, int, str]]``) — the ψ½
  ANGULAR-edge pole seed.

This module pins the two load-bearing structural facts the carve must preserve:

1. **Sibling, NOT child.** The :class:`~orpheus.transport.full_field.FullField`
   composite discriminates its boundary slot by ``isinstance(·, BoundaryField)``;
   the ψ½ pole MUST FAIL that test. If a refactor made
   :class:`StartingDirectionField` a *child* of :class:`BoundaryField` (or
   collapsed :class:`BoundaryField` into :class:`FaceField`), the composite would
   silently accept a pole field in its boundary slot — a wrong-slot bug invisible
   to every value gate.
2. **The metric descends PER LEAF, never on the ABC** (ERR-067). :class:`FaceField`
   unifies *structure* only; the Hilbert metric lives on each leaf's space (pole
   ``V_cell``; trace ``|Ω·n̂|·w``). Unifying a single metric onto the shared ABC
   would re-introduce the ghost-metric class of bug.

All assertions are ``np.testing.assert_*`` (function calls — they fire under the
canonical ``python -O``; vv Mode 8).
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.transport.fields._bases import (
    BulkField,
    Field,
    FaceField,
    BoundaryField,
    AngularBoundaryField,
    ScalarBoundaryField,
    StartingDirectionField,
)
from orpheus.transport.fields.angular_boundary_flux import AngularBoundaryFlux
from orpheus.transport.fields.starting_direction_flux import StartingDirectionFlux
from tests.sn._test_helpers import make_tiny_spherical_sn_mesh


pytestmark = [pytest.mark.foundation]


def _sphere_mesh():
    r"""A sphere-GL SNMesh — carries ONE seed-carrying μ-level (a ψ½ block)
    AND an angular boundary trace, so both codim-1 loci are constructible."""
    return make_tiny_spherical_sn_mesh()


# ── A. The class-level hierarchy ─────────────────────────────────────────


class TestFaceFieldIsTheCodim1Parent:
    def test_every_face_family_descends_from_facefield(self):
        for family in (
            BoundaryField,
            StartingDirectionField,
            AngularBoundaryField,
            ScalarBoundaryField,
        ):
            np.testing.assert_(
                issubclass(family, FaceField),
                f"{family.__name__} is not a FaceField",
            )
            np.testing.assert_(issubclass(family, Field))

    def test_facefield_is_codim1_disjoint_from_the_bulk_locus(self):
        # FaceField is codim-1; BulkField is codim-0. Neither is under the other.
        np.testing.assert_(not issubclass(FaceField, BulkField))
        np.testing.assert_(not issubclass(BulkField, FaceField))


class TestPoleIsSiblingNotChild:
    r"""THE load-bearing guard for the FullField boundary-slot discriminator."""

    def test_pole_is_not_a_boundaryfield_subclass(self):
        # If this flips, isinstance(pole, BoundaryField) in FullField.__post_init__
        # starts returning True and the composite accepts the pole in its
        # boundary slot — silently wrong.
        np.testing.assert_(not issubclass(StartingDirectionField, BoundaryField))
        np.testing.assert_(not issubclass(BoundaryField, StartingDirectionField))

    def test_pole_and_boundary_share_the_facefield_parent(self):
        np.testing.assert_(issubclass(StartingDirectionField, FaceField))
        np.testing.assert_(issubclass(BoundaryField, FaceField))


# ── B. Constructed-instance discrimination ───────────────────────────────


class TestConstructedInstancesDiscriminate:
    def test_constructed_pole_fails_the_boundary_isinstance(self):
        mesh = _sphere_mesh()
        pole = StartingDirectionFlux.zeros_on(mesh)
        boundary = AngularBoundaryFlux.zeros_on(mesh)
        # The exact test FullField.__post_init__ runs on its boundary slot.
        np.testing.assert_(not isinstance(pole, BoundaryField))  # the guarded fact
        np.testing.assert_(isinstance(pole, FaceField))
        np.testing.assert_(isinstance(boundary, BoundaryField))
        np.testing.assert_(isinstance(boundary, FaceField))


# ── C. __post_init__ fires exactly once (no double validation via MRO) ───


class TestPostInitFiresOnce:
    def test_facefield_post_init_runs_exactly_once_per_construction(self, monkeypatch):
        r"""The re-parent must not double-run the layout/shape validation.

        Each leaf's ``__post_init__`` calls ``super().__post_init__()`` exactly
        once up the (single-inheritance) chain to :meth:`FaceField.__post_init__`.
        A spy on the ABC's validator confirms one call per constructed field —
        a regression that re-declared ``__post_init__`` on the intermediate
        :class:`BoundaryField` (double validation) would count 2.
        """
        mesh = _sphere_mesh()
        calls = {"n": 0}
        original = FaceField.__post_init__

        def counting_post_init(self) -> None:
            calls["n"] += 1
            original(self)

        monkeypatch.setattr(FaceField, "__post_init__", counting_post_init)

        StartingDirectionFlux.zeros_on(mesh)
        np.testing.assert_equal(calls["n"], 1, err_msg="pole: post_init != once")

        calls["n"] = 0
        AngularBoundaryFlux.zeros_on(mesh)
        np.testing.assert_equal(calls["n"], 1, err_msg="boundary: post_init != once")


# ── D. ERR-067 — the metric descends per leaf, never on the ABC ──────────


class TestMetricIsPerLeaf:
    @pytest.mark.catches("ERR-067")
    def test_metric_descends_per_leaf_not_on_the_facefield_abc(self):
        r"""FaceField unifies STRUCTURE only — the metric stays per-leaf.

        ERR-067 tripwire: the pole's state metric is SPD ``V_cell`` (never the
        retired ghost ``0``), and the two codim-1 leaves carry STRUCTURALLY
        DISTINCT metrics (pole ``V_cell`` vs trace ``|Ω·n̂|·w``). A refactor that
        pulled a single uniform metric onto the shared :class:`FaceField` ABC
        would collapse them to equal — this reddens.
        """
        mesh = _sphere_mesh()
        pole_w = StartingDirectionFlux.zeros_on(mesh).space.inner_product_weights
        trace_w = AngularBoundaryFlux.zeros_on(mesh).space.inner_product_weights

        # The pole's state metric exists and is SPD (nonzero) — NOT the ghost 0.
        np.testing.assert_(pole_w is not None, "pole carries no state metric")
        np.testing.assert_(
            bool(np.all(np.asarray(pole_w) > 0.0)),
            "pole metric is not strictly positive (ghost G_sd=0 regression?)",
        )
        np.testing.assert_(trace_w is not None, "trace carries no metric")

        # Per-leaf ⇒ the two metrics are NOT the same array (would be equal only
        # if a single uniform metric were unified onto the shared ABC).
        pw, tw = np.asarray(pole_w), np.asarray(trace_w)
        unified = pw.shape == tw.shape and np.array_equal(pw, tw)
        np.testing.assert_(
            not unified,
            "pole and trace metrics are identical — the per-leaf metric may "
            "have been unified onto the FaceField ABC (ERR-067 regression).",
        )
