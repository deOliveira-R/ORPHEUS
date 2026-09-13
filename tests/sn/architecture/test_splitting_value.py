r"""The Strategy VALUE's defining laws — :class:`orpheus.sn.splitting.Splitting`.

The consumers campaign's step 2, C1 (2026-09-13): the splitting ``A = M − N``
of the within-group loss stopped being a field of the Problem's posed record
and became a VALUE minted from the record's factors
(:meth:`~orpheus.sn.splitting.Splitting.from_schedule`, the one labelling
site).  Every math-bearing type ships a test of its DEFINING laws, not just
its usage (``feedback_test_intrinsic_properties``); these are the value's.

What each row pins, and the measurement it was written against
(``scratch/_consumers/planning/test_architect_step2_delta.md`` §B1-bis, the
``_config`` fixtures, ``[M]`` 2026-09-13):

* **the law** — ``(M − ΣN − A)·x`` is EXACTLY ``0.0`` on the seedless arm for
  both schedules (one flat operator sum; the Gauss-Seidel split writes
  disjoint rows so no addition is reordered) and ``≤ 2.842e-14`` absolute /
  ``2.087e-16`` relative to ``‖A·x‖∞`` on the CARRYING arm over 40 draws (the
  2×2 grid re-associates the sums) — so the carrying contract is a
  RELATIVE band, never ``array_equal`` (a bit-identity claim there is a
  false red, ``vv-principles`` §bit-identity);
* **the partition** — a term carries exactly one label.  Double-labelling a
  piece and DROPPING it read the SAME law defect (``2.778702e+00`` both on the
  2-D box), so the law alone cannot tell them apart; the partition is
  therefore refused AT CONSTRUCTION;
* **transfer invariance** — moving the absorbed ``B_lower`` from the implicit
  to the explicit side leaves the law at ``0.0`` (the sign travels WITH the
  term as data, ``LossTerm.sign``; had it been positional the transfer would
  read ``2·|piece·x|``) — the property P6's ``transfer`` rides on;
* **teeth** — dropping a lagged piece moves the law by ``|piece·x|``
  (``6.746694e-01`` on the 2-D box for ``S``); flipping a term's sign by
  ``2·|term·x|`` (``1.349339e+00``); a corrupt sign (``±2``) is unspellable;
* **identity** — the Jacobi value's members ARE the record's factors (``is``),
  the Gauss-Seidel value's ``M`` is the scheduled composite; a carrying
  record refuses a sequenced schedule; the entry string is resolved ONCE and
  a 1-D mesh falls back to Jacobi under ``"gauss_seidel"``.
"""
from __future__ import annotations

import dataclasses

import numpy as np
import pytest

from orpheus.sn.loss_representation.sweep_schedule import (
    SweepSchedule,
    reflective_faces,
)
from orpheus.sn.operators.boundary import SNMaskedBoundaryOperator
from orpheus.sn.operators.scheduled_invertible import ScheduledInvertibleOperator
from orpheus.sn.splitting import LossTerm, Splitting, resolve_schedule
from tests.sn.architecture._config import (
    cart2d_seedless,
    random_state,
    record_for,
    slab_seedless,
    sphere_carrying,
)

pytestmark = pytest.mark.foundation

_SEED = 20260913
_SCHEDULES = ("jacobi", "gauss_seidel")


def _value(build_mesh, schedule: str) -> Splitting:
    sn_mesh = build_mesh()
    record = record_for(sn_mesh)
    return Splitting.from_schedule(record, resolve_schedule(sn_mesh, schedule))


# ═══════════════════════════════════════════════════════════════════════
# The law  A = M − N, per value
# ═══════════════════════════════════════════════════════════════════════


class TestLawTheSplittingLaw:
    @pytest.mark.parametrize("schedule", _SCHEDULES)
    def test_seedless_law_is_exactly_zero(self, schedule: str) -> None:
        """One flat operator sum over the record's own terms — 0 ULP."""
        value = _value(cart2d_seedless, schedule)
        state = random_state(value.system, seed=_SEED)
        np.testing.assert_array_equal(
            value.law_residual(state), 0.0,
            err_msg=f"M − N != A on the seedless arm under {schedule!r}",
        )

    def test_carrying_law_within_the_grid_reassociation_band(self) -> None:
        """The 2×2 grid re-associates the sums: a RELATIVE band, measured
        ``2.087e-16`` over 40 draws (``≈ 0.94 ε``); gated at 8 ε."""
        value = _value(sphere_carrying, "jacobi")
        state = random_state(value.system, seed=_SEED)
        scale = float(np.max(np.abs(value.system.loss.apply(state).to_flat())))
        residual = float(np.max(np.abs(value.law_residual(state))))
        assert scale > 0.0, "the probe state does not excite A"
        assert residual <= 8.0 * np.finfo(float).eps * scale, (
            f"carrying-arm law residual {residual:.3e} exceeds 8 ε · ‖A·x‖∞ "
            f"= {8.0 * np.finfo(float).eps * scale:.3e}"
        )

    def test_transfer_across_the_boundary_is_law_invariant(self) -> None:
        """Re-labelling the absorbed ``B_lower`` as explicit keeps the law at
        exactly 0 — the sign is DATA of the term, so a transfer needs no
        arithmetic (P6's ``transfer`` in embryo)."""
        value = _value(cart2d_seedless, "gauss_seidel")
        lower = [t for t in value.implicit_pieces if isinstance(t.operator, SNMaskedBoundaryOperator)]
        assert len(lower) == 1, "the G-S value must absorb exactly one boundary piece"
        moved = dataclasses.replace(
            value,
            implicit_pieces=tuple(t for t in value.implicit_pieces if t is not lower[0]),
            explicit_pieces=value.explicit_pieces + (lower[0],),
        )
        state = random_state(value.system, seed=_SEED)
        np.testing.assert_array_equal(
            moved.law_residual(state), 0.0,
            err_msg="a label transfer moved the law — the sign did not travel with the term",
        )


# ═══════════════════════════════════════════════════════════════════════
# The partition, and the teeth
# ═══════════════════════════════════════════════════════════════════════


class TestLawThePartition:
    def test_a_term_labelled_twice_is_refused_at_construction(self) -> None:
        """Double-label and drop read the SAME law defect (measured 2.78 on the
        2-D box), so the partition cannot be diagnosed from the law — it is
        refused where the value is made."""
        value = _value(cart2d_seedless, "gauss_seidel")
        lower = [t for t in value.implicit_pieces if isinstance(t.operator, SNMaskedBoundaryOperator)][0]
        with pytest.raises(ValueError, match="exactly one label"):
            dataclasses.replace(value, explicit_pieces=value.explicit_pieces + (lower,))

    @pytest.mark.parametrize("schedule", _SCHEDULES)
    def test_dropping_a_lagged_piece_reds_the_law(self, schedule: str) -> None:
        """Teeth: the law moves by ``|S·x|`` (``6.746694e-01`` measured)."""
        value = _value(cart2d_seedless, schedule)
        state = random_state(value.system, seed=_SEED)
        dropped = dataclasses.replace(value, explicit_pieces=value.explicit_pieces[1:])
        defect = float(np.max(np.abs(dropped.law_residual(state))))
        assert defect > 1e-3, f"dropping the first lagged piece moved the law by only {defect:.3e}"

    def test_flipping_a_terms_sign_reds_the_law(self) -> None:
        """Teeth: ``2·|S·x|`` (``1.349339e+00`` measured) — the in-class
        mutation the signed-term design creates."""
        value = _value(cart2d_seedless, "jacobi")
        state = random_state(value.system, seed=_SEED)
        first = value.explicit_pieces[0]
        flipped = dataclasses.replace(
            value,
            explicit_pieces=(LossTerm(first.operator, -first.sign),) + value.explicit_pieces[1:],
        )
        defect = float(np.max(np.abs(flipped.law_residual(state))))
        assert defect > 1e-3, f"flipping a term's sign moved the law by only {defect:.3e}"

    def test_a_corrupt_sign_is_unspellable(self) -> None:
        value = _value(cart2d_seedless, "jacobi")
        with pytest.raises(ValueError, match=r"\+1 or -1"):
            LossTerm(value.explicit_pieces[0].operator, 2)


# ═══════════════════════════════════════════════════════════════════════
# Identity and the labelling site
# ═══════════════════════════════════════════════════════════════════════


class TestLawTheLabelling:
    def test_jacobi_members_are_the_records_factors_by_identity(self) -> None:
        value = _value(cart2d_seedless, "jacobi")
        f = value.system.factors
        assert value.implicit is f.streaming_collision
        assert tuple(map(id, value.explicit)) == (id(f.scattering), id(f.n2n), id(f.boundary))
        assert not value.schedule.is_sequenced

    def test_gauss_seidel_absorbs_the_lower_boundary_half(self) -> None:
        value = _value(cart2d_seedless, "gauss_seidel")
        assert value.schedule.is_sequenced
        assert isinstance(value.implicit, ScheduledInvertibleOperator)
        assert value.implicit.invertible is value.system.factors.streaming_collision
        assert isinstance(value.explicit[-1], SNMaskedBoundaryOperator)
        # the collision gains ride both labellings unchanged
        assert value.explicit[0] is value.system.factors.scattering
        assert value.explicit[1] is value.system.factors.n2n

    def test_a_carrying_record_refuses_a_sequenced_schedule(self) -> None:
        sn_mesh = sphere_carrying()
        record = record_for(sn_mesh)
        sequenced = SweepSchedule.gauss_seidel(
            sn_mesh.ndim, sn_mesh.quad.octants, reflective_faces(sn_mesh),
        )
        with pytest.raises(ValueError, match="carrying mesh"):
            Splitting.from_schedule(record, sequenced)

    def test_the_entry_string_resolves_once_and_1d_falls_back_to_jacobi(self) -> None:
        slab = slab_seedless()
        assert not resolve_schedule(slab, "gauss_seidel").is_sequenced
        assert resolve_schedule(cart2d_seedless(), "gauss_seidel").is_sequenced
        with pytest.raises(ValueError, match="Unknown inner_schedule"):
            resolve_schedule(slab, "sor")
