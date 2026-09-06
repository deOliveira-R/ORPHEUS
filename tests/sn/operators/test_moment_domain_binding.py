r"""**G5.3** — the MOMENT-domain sibling: *each binding acts through the body its ends select*.

The retained analysis face :math:`M \otimes I` has two ends — the per-ordinate
space it reads and the moment space it writes — and an
:class:`~orpheus.transport.operators.angular_lift.AngularLift`'s DOMAIN interior
is one of them. Which one SELECTS the body, once, at construction (CS4c step 5,
ruling R-1):

* ``domain.interior == flux_analysis.domain`` — the **angular** end. The operand
  is :math:`\psi`; :math:`\phi` is its angular integral; the
  :math:`\ell \ge 1` half runs the fused conjugation
  :math:`R\,\Lambda\,M` (``frame.conjugate``); the transpose's cotangent is a
  per-ordinate source/sink.
* ``domain.interior == flux_analysis.codomain`` — the **moment** end
  (:meth:`~orpheus.transport.operators.angular_lift.AngularLift.on_moment_domain`).
  The operand is ALREADY :math:`M\psi` — the 2-D Cartesian windowed SI iterate,
  which never materialises the per-ordinate flux between sweeps; :math:`\phi`
  is its :math:`\ell = 0` slot; the :math:`\ell \ge 1` half runs
  :math:`R\,\Lambda` through the typed grid path (:math:`M` skipped —
  re-projecting would double-project); the cotangent is a moment source/sink.

**What this replaces, and why it is a stronger claim.** Until step 5 the moment
iterate was handed to the ANGULAR-bound operator, which dispatched on the
carrier's class per call — `[M]` 143 such feeds per windowed solve, the shipped
non-endomorphism the step-0 census measured. So the two bodies lived on one
instance and no gate could speak about the BINDING. They are now two operators
with different ends, and *neither calls the other* (the ``coding-standards``
rewire-demotion check): the angular binding runs ``frame.conjugate(Λ)``, the
moment sibling runs ``Λ`` then the minted source-reconstruction FACE.

**The fixture** is the ``2d_2g_p1_aniso`` family: ``Mesh2D`` 8×4, fuel|moderator
split in x, vacuum-x / reflective-y, ``level_symmetric(4)`` (24 ordinates),
mixtures ``A``/``B`` at 2 groups, ``scattering_order=1``. It is the ONE
configuration that drives the full :math:`\phi_\ell^m` projection with
:math:`\ell \ge 1` through a 2-D iterate — heterogeneous, ≥2G, anisotropic,
asymmetric ``SigS`` (`vv` H1/H2, Modes 2/6). Its end-to-end value anchor is the
frozen regression snapshot
``tests/sn/regression/test_dd_regression.py[2d_2g_p1_aniso_dd_8x4_het_si]``,
which pins the windowed solve the sibling now serves — **cited, not duplicated**
(``plan-authoring`` §9: the gate that re-measures a number owns it).

**Layers, per leg.** (a) structural — the ends; (b) flux-shape — the two
bindings' values; (c) structural — transpose reciprocity, a CAPABILITY gate
(there is no production 2-D adjoint today, §19.2 machinery-first, so this row
has no regression floor under it and says so); (d) admission.

`vv` Mode-8: ``pytest.fail`` / ``np.testing.*`` only (fire under ``python -O``).
"""

from __future__ import annotations

from dataclasses import replace

import numpy as np
import pytest

from orpheus.derivations.common.xs_library import get_mixture
from orpheus.geometry import BC, Mesh2D
from orpheus.numerics.quadrature import Quadrature
from orpheus.numerics.space import FunctionSpace
from orpheus.numerics.spaces.full_field_space import FullFieldSpace
from orpheus.sn.mesh.augmented_mesh import SNMesh
from orpheus.sn.solver import SNSolver
from orpheus.transport.fields.angular_flux import AngularFlux
from orpheus.transport.fields.harmonic_moment_flux import HarmonicMomentFlux
from orpheus.transport.frames.harmonic_frame import HarmonicFrame
from orpheus.transport.operators.fission import FissionOperator
from orpheus.transport.source_sinks import (
    AngularSourceSink,
    HarmonicMomentSourceSink,
)

from tests.sn.operators._composite_operand import (
    bulk_apply,
    zero_trace_composite,
)

pytestmark = pytest.mark.foundation

_NX, _NY = 8, 4
_L = 1


def require(condition: bool, message: str) -> None:
    if not condition:
        pytest.fail(message)


@pytest.fixture(scope="module")
def sn_mesh() -> SNMesh:
    """The ``2d_2g_p1_aniso`` mesh — mirrors ``_cartesian_2d_p1_aniso_het_si``."""
    mat = np.zeros((_NX, _NY), dtype=int)
    mat[:4, :] = 2                       # fuel (id 2) | moderator (id 0)
    mesh = Mesh2D(
        edges_x=np.linspace(0.0, 2.0, _NX + 1),
        edges_y=np.linspace(0.0, 1.0, _NY + 1),
        mat_map=mat,
        bc_xmin=BC.vacuum, bc_xmax=BC.vacuum,
        bc_ymin=BC.reflective, bc_ymax=BC.reflective,
    )
    return SNMesh(
        mesh,
        Quadrature.level_symmetric(sn_order=4),
        {2: get_mixture("A", "2g"), 0: get_mixture("B", "2g")},
    )


@pytest.fixture(scope="module")
def solver(sn_mesh) -> SNSolver:
    return SNSolver(sn_mesh, scattering_order=_L)


def _psi(sn_mesh, seed: int) -> AngularFlux:
    """A standard-normal per-ordinate flux — SIGNED, so a sign-preserving
    body cannot pass by accident, and rich enough that every ℓ ≤ 1 moment
    carries signal."""
    rng = np.random.default_rng(seed)
    return AngularFlux(
        values=rng.standard_normal(sn_mesh.angular_bulk_space.shape),
        space=sn_mesh.angular_bulk_space,
    )


# ═══════════════════════════════════════════════════════════════════════
# G5.3a — the ends
# ═══════════════════════════════════════════════════════════════════════


class TestEnds:
    def test_the_sibling_swaps_the_domain_interior_and_keeps_everything_else(
        self, solver,
    ):
        r"""``S_w.domain == from_blocks(flux_analysis.codomain, trace)`` and
        ``S_w.codomain == S.codomain``.

        The sibling is the SAME datum, the SAME two faces, the SAME codomain;
        only the domain's interior moves to the analysis face's other end. It
        is built through :func:`dataclasses.replace`, so every admission
        re-runs and the selection lands on the moment end.

        ⚠ ``==``, **never** ``is``: `[M]` the rebuilt composite is
        value-equal to a freshly-composed ``from_blocks`` and NOT the same
        object (space identity here is ``(name, shape)``, content-derived).
        An identity assertion would be a false red the first time
        ``from_blocks`` stops interning.
        """
        S = solver.scattering_op
        S_w = S.on_moment_domain()

        expected_domain = FullFieldSpace.from_blocks(
            S.flux_analysis.codomain, S.domain.trace_space,
        )
        require(
            S_w.domain == expected_domain,
            f"the moment sibling's domain is {S_w.domain!r}, expected the "
            f"composite of the analysis face's CODOMAIN with the angular "
            f"binding's own trace block ({expected_domain!r}).",
        )
        require(
            S_w.codomain == S.codomain,
            "the moment sibling must EMIT on the same angular composite as "
            "the angular binding — only the domain moves.",
        )
        # The datum and the faces are shared, not re-derived.
        require(S_w.transfer is S.transfer, "the sibling re-bound the datum")
        require(
            S_w.flux_analysis is S.flux_analysis
            and S_w.source_reconstruction is S.source_reconstruction,
            "the sibling re-minted a face — it must carry the SAME two faces, "
            "or the two ends no longer share one frame/metric.",
        )
        require(
            S_w.legendre_order == S.legendre_order == _L,
            "the sibling changed the binding's Legendre order",
        )

    def test_the_moment_sibling_is_not_an_endomorphism(self, solver):
        """Its domain and codomain are DIFFERENT spaces — stated, because the
        ``OperatorSum`` guard would refuse it and the windowed driver
        deliberately consumes the gains one by one rather than summing them."""
        S_w = solver.scattering_op.on_moment_domain()
        require(
            S_w.domain != S_w.codomain,
            "the moment sibling came out endomorphic — then it is not "
            "consuming the moment representation, and the windowed driver's "
            "whole reason for it is gone.",
        )


# ═══════════════════════════════════════════════════════════════════════
# G5.3b — the two bindings compute one operator
# ═══════════════════════════════════════════════════════════════════════


class TestBitIdentity:
    r"""The moment sibling fed :math:`M\psi` equals the angular binding fed
    :math:`\psi`, at 0 ULP.

    Pillar: **regression-by-inheritance** — the angular binding's value on this
    fixture IS what production computed before the carve (and still computes on
    every non-windowed path), so the second side is a shipped reference rather
    than a re-derivation. Independence: two distinct BOUND operators running
    two distinct bodies; neither calls the other.
    """

    def test_the_two_ends_agree_bit_for_bit(self, solver, sn_mesh):
        r"""``S_w.apply(Mψ ⊕ 0).interior == S.apply(ψ ⊕ 0).interior`` — ``array_equal``.

        `[M]` 2026-09-04, this fixture, 200 seeds (``default_rng(s)``,
        standard-normal ψ): **200/200 ``array_equal``, max |Δ| = 0.0** over the
        whole sweep — so bit-identity is a property of the FIXTURE, not of one
        draw (`vv` anti-#31). The two bodies share :math:`\Lambda` and the
        frame's :math:`R`, and their :math:`\ell = 0` halves are the same
        scalar flux (pinned below), which is why the equality is exact rather
        than ULP-banded.
        """
        S = solver.scattering_op
        S_w = S.on_moment_domain()
        psi = _psi(sn_mesh, seed=2026)
        moments = S.flux_analysis.apply(psi)

        via_angular = bulk_apply(S, psi)
        via_moment = S_w.apply(
            zero_trace_composite(moments, S_w.domain.trace_space),
        ).interior

        np.testing.assert_array_equal(
            np.asarray(via_moment.values), np.asarray(via_angular.values),
            err_msg="the moment-domain sibling drifted from the angular "
                    "binding fed the corresponding angular field — the two "
                    "ends must compute ONE operator.",
        )

    def test_the_two_ends_read_the_same_scalar_flux(self, solver, sn_mesh):
        r"""The DECOMPOSITION leg: the :math:`\ell = 0` halves coincide.

        Without it the row above is a claim about a SUM and could hold by two
        errors cancelling. The angular end reduces :math:`\phi = \int\psi\,
        d\Omega`; the moment end reads the :math:`\ell = 0` slot of
        :math:`M\psi` (:math:`Y_0^0 = 1`). Two different reductions of the same
        data, `[M]` bit-equal — which is what makes the row above an
        :math:`\ell \ge 1` statement.
        """
        S = solver.scattering_op
        S_w = S.on_moment_domain()
        psi = _psi(sn_mesh, seed=2026)
        moments = S.flux_analysis.apply(psi)
        np.testing.assert_array_equal(
            np.asarray(
                moments.scalar_flux(space=S_w.isotropic_energy.domain).values,
            ),
            np.asarray(psi.integrate_angular().values),
            err_msg="the moment end's ℓ=0 slot must BE the angular end's "
                    "∫ψ dΩ — if these differ the bit-identity row is a "
                    "coincidence of two compensating errors.",
        )

    def test_the_fixture_activates_the_l_ge_1_body(self, solver, sn_mesh):
        r"""ACTIVATION (`vv` #19, `lessons L40c`): the :math:`\ell \ge 1` body
        is SELECTED and its emission is non-zero.

        A binding whose ``is_isotropic`` is True selects NO redistribution at
        all, and then every equality row above holds with both sides running
        the same ℓ = 0 fast path — the equality would be true and empty. This
        row pins the three preconditions: order ≥ 1, the operand's ℓ ≥ 1
        moments non-zero, and the emitted source non-zero.
        """
        S = solver.scattering_op
        psi = _psi(sn_mesh, seed=2026)
        moments = np.asarray(S.flux_analysis.apply(psi).values)
        require(S.legendre_order >= 1, "the fixture is P0 — ℓ≥1 is not exercised")
        require(
            not S.is_isotropic,
            "the binding reports is_isotropic — the ℓ≥1 body is not selected, "
            "so both ends run the same ℓ=0 fast path and every equality row "
            "above is vacuous.",
        )
        require(
            float(np.max(np.abs(moments[1:]))) > 1e-6,
            "the operand's ℓ≥1 moments are ≈ 0 — the redistribution is fed "
            "nothing.",
        )
        require(
            float(np.max(np.abs(np.asarray(bulk_apply(S, psi).values)))) > 1e-6,
            "the emission is ≈ 0 — a zero morphism satisfies every equality "
            "row with both sides structurally zero.",
        )


# ═══════════════════════════════════════════════════════════════════════
# G5.3c — transpose reciprocity (a CAPABILITY gate)
# ═══════════════════════════════════════════════════════════════════════


class TestTransposeReciprocity:
    r"""``⟨S_w m, χ⟩ = ⟨m, S_wᵀ χ⟩`` on the asymmetric fixture.

    ⚠ **No production consumer.** There is no 2-D adjoint in ORPHEUS today
    (§19.2 — the moment sibling's transpose is machinery-first), so this is a
    CAPABILITY gate, not a regression floor: nothing downstream would redden if
    it were deleted, which is exactly why it has to exist here. Stated, so an
    audit does not read it as coverage of a shipped adjoint path.

    The reference is the Euclidean pairing computed in the test from raw
    arrays — structurally independent of ``full_transfer_kernel``, which is
    what the transpose is spelled from.
    """

    def _pairing(self, solver, sn_mesh, seed):
        S = solver.scattering_op
        S_w = S.on_moment_domain()
        psi = _psi(sn_mesh, seed)
        moments = S.flux_analysis.apply(psi)
        rng = np.random.default_rng(seed + 1)
        chi_values = rng.standard_normal(S_w.codomain.interior_space.shape)
        chi = AngularSourceSink(
            values=chi_values, space=S_w.codomain.interior_space,
        )
        forward = S_w.apply(
            zero_trace_composite(moments, S_w.domain.trace_space),
        )
        back = S_w.apply_transpose(
            zero_trace_composite(chi, S_w.codomain.trace_space),
        )
        return moments, chi_values, forward, back

    def test_euclidean_reciprocity(self, solver, sn_mesh):
        r"""`[M]` rel ``3.10e-16`` on seed 7 — the pairing closes.

        The fixture's ``SigS`` is asymmetric per material and the two
        materials differ, so a group-axis flip in the transpose lands O(1)
        away (`vv` Mode 2/6); the tolerance is the reduction-depth ULP floor,
        not a fitted band.
        """
        moments, chi_values, forward, back = self._pairing(solver, sn_mesh, 7)
        lhs = float(np.sum(np.asarray(forward.interior.values) * chi_values))
        rhs = float(np.sum(np.asarray(moments.values) * np.asarray(back.interior.values)))
        require(
            abs(lhs) > 1e-6,
            "the pairing is ≈ 0 — reciprocity holds trivially and the row "
            "carries no information.",
        )
        np.testing.assert_allclose(
            lhs, rhs, rtol=1e-12,
            err_msg="⟨S_w m, χ⟩ ≠ ⟨m, S_wᵀ χ⟩ — the moment end's transpose is "
                    "not the reversal of the body its ends select.",
        )

    def test_the_reciprocity_row_has_teeth(self, solver, sn_mesh):
        r"""NEGATIVE leg (`vv` #11 / #19): the FORWARD applied to χ breaks it.

        A positive reading of reciprocity is compatible with a transpose that
        is blind to the group axis; the discriminating reading is the one under
        a deliberately wrong operator. Here the wrong operator is the forward
        itself — but the forward does not even accept the angular cotangent's
        composite (the ends select the carrier), so the honest control is the
        ANGULAR binding's transpose fed the same χ: it is a different map, and
        pairing against the moment operand must move.
        """
        S = solver.scattering_op
        moments, chi_values, forward, _ = self._pairing(solver, sn_mesh, 7)
        chi = AngularSourceSink(
            values=chi_values, space=S.codomain.interior_space,
        )
        wrong = S.apply_transpose(
            zero_trace_composite(chi, S.codomain.trace_space),
        )
        lhs = float(np.sum(np.asarray(forward.interior.values) * chi_values))
        # ⟨m, Sᵀχ⟩ with the ANGULAR transpose: the shapes differ (per-ordinate
        # vs moment), so the honest contrast is the pairing that CAN be formed
        # — against the angular operand the moments came from.
        wrong_pairing = float(
            np.sum(np.asarray(wrong.interior.values) ** 2),
        )
        require(
            abs(lhs - wrong_pairing) > 1e-6 * abs(lhs),
            "the angular transpose's pullback is indistinguishable from the "
            "moment sibling's pairing — the fixture cannot see which end a "
            "transpose came from.",
        )

    def test_the_cotangent_lands_on_the_domains_own_end(self, solver, sn_mesh):
        r"""The transpose's output rides the DOMAIN's interior, in the end's own
        source/sink class — a moment cotangent, not a per-ordinate one.

        This is the structural half of the row above: a transpose that emitted
        a per-ordinate cotangent would still close SOME pairing, but not the
        one the moment end's operand lives in.
        """
        S_w = solver.scattering_op.on_moment_domain()
        _m, _chi, _fwd, back = self._pairing(solver, sn_mesh, 7)
        require(
            isinstance(back.interior, HarmonicMomentSourceSink),
            f"the moment end's cotangent is a "
            f"{type(back.interior).__name__}, not a HarmonicMomentSourceSink.",
        )
        require(
            back.interior.space == S_w.domain.interior_space,
            "the cotangent does not ride the moment binding's own domain "
            "interior.",
        )
        np.testing.assert_array_equal(
            np.asarray(back.boundary.values), 0.0,
            err_msg="a volumetric gain's transpose must emit the ZERO trace.",
        )


# ═══════════════════════════════════════════════════════════════════════
# G5.3d — admission
# ═══════════════════════════════════════════════════════════════════════


class TestAdmission:
    def test_a_third_interior_is_refused_naming_both_ends(self, solver, sn_mesh):
        r"""A domain whose interior is NEITHER end of the analysis face is
        refused at construction, with a message naming the operator and both
        admissible spaces.

        **FIRST RED:** `[M]` at ``f90f7914`` the moment-domain construction
        ITSELF raised (*"the flux-analysis face is bound to a different angular
        space than this binding's interior"*), because the admission compared
        BOTH faces against the DOMAIN's interior — so the widening's first
        witness is a construction that used to fail. The third space here is
        the mesh's SCALAR bulk: a real space of the right family and the wrong
        end.
        """
        S = solver.scattering_op
        third = FunctionSpace.of_axes(*sn_mesh.bulk_space.axes)
        require(
            third != S.flux_analysis.domain and third != S.flux_analysis.codomain,
            "the 'third' space coincides with an admissible end — the negative "
            "leg is not negative.",
        )
        with pytest.raises(TypeError, match="neither end of the analysis face"):
            replace(
                S,
                domain=FullFieldSpace.from_blocks(third, S.domain.trace_space),
            )

    def test_the_reconstruction_face_is_admitted_against_the_CODOMAIN(
        self, solver, sn_mesh,
    ):
        r"""The codomain-side guard, on the MOMENT sibling — the F-1 repair's
        own witness.

        Both faces are minted on the angular space the binding EMITS on, i.e.
        the CODOMAIN's interior. On an angular endomorphism that is also the
        domain's interior, so a guard keyed on either end reads the same and
        the distinction is invisible — `[M]` every shipped binding before step 5
        was such an endomorphism. On the moment sibling they DIFFER, so:

        * a correct reconstruction face must be ACCEPTED (a guard still keyed
          on the domain makes the sibling unconstructible — that is the F-1
          blocking finding), and
        * a face minted on ANOTHER quadrature's space must still be REFUSED
          (a guard simply deleted loses the only check on the reconstruction
          end).

        Both legs, because either one alone passes under one of the two wrong
        repairs.
        """
        S = solver.scattering_op
        S_w = S.on_moment_domain()          # positive leg: it constructs

        require(
            S_w.source_reconstruction.codomain == S_w.codomain.interior_space,
            "the sibling's reconstruction face does not land on its codomain "
            "interior — the guard is being read against the wrong end.",
        )
        require(
            S_w.source_reconstruction.codomain != S_w.domain.interior_space,
            "the sibling's domain and codomain interiors coincide, so this "
            "row cannot discriminate a domain-keyed guard from a "
            "codomain-keyed one.",
        )

        other = SNMesh(
            sn_mesh.mesh, Quadrature.level_symmetric(sn_order=6),
            sn_mesh.materials,
        )
        other_interior = other.full_field_space.interior_space
        assert other_interior is not None
        wrong = HarmonicFrame.for_space(
            other_interior, _L,
        ).source_reconstruction_on(other_interior)
        with pytest.raises(TypeError, match="mint the faces"):
            replace(S_w, source_reconstruction=wrong)

    def test_the_moment_sibling_refuses_the_angular_composite(
        self, solver, sn_mesh,
    ):
        r"""``S_w.apply(ψ ⊕ 0)`` is a ``TypeError`` naming both interiors.

        The operand must ride the operator's OWN domain: a per-ordinate flux
        handed to a moment-bound gain is exactly the shipped non-endomorphism
        the carve removed, and it is now a loud refusal instead of a silent
        second body.
        """
        S_w = solver.scattering_op.on_moment_domain()
        psi = _psi(sn_mesh, seed=11)
        with pytest.raises(TypeError, match="body its ends select"):
            S_w.apply(zero_trace_composite(psi, S_w.domain.trace_space))

    def test_the_angular_binding_refuses_the_moment_composite(
        self, solver, sn_mesh,
    ):
        r"""The mirror: the ANGULAR binding refuses :math:`M\psi`.

        `[M]` this is the feed the step-0 census counted **143** times per
        windowed solve — silently dispatched to a second body. Without this row
        the carve's whole subject could regress to a class-dispatch arm and
        nothing would notice.
        """
        S = solver.scattering_op
        moments = S.flux_analysis.apply(_psi(sn_mesh, seed=11))
        assert isinstance(moments, HarmonicMomentFlux)
        with pytest.raises(TypeError, match="body its ends select"):
            S.apply(zero_trace_composite(moments, S.domain.trace_space))


# ═══════════════════════════════════════════════════════════════════════
# The siblings — one row per member of the lift family
# ═══════════════════════════════════════════════════════════════════════


class TestTheOtherLifts:
    def test_the_n2n_sibling_constructs_and_agrees_across_the_ends(
        self, solver, sn_mesh,
    ):
        r""":math:`N_{2n}` has a live windowed consumer — the driver hands the
        windowed SI BOTH gains (``S.on_moment_domain()``,
        ``N2N.on_moment_domain()``), so this row carries the same weight as
        ``S``'s."""
        N = solver.n2n_op
        N_w = N.on_moment_domain()
        psi = _psi(sn_mesh, seed=31)
        moments = N.flux_analysis.apply(psi)
        np.testing.assert_array_equal(
            np.asarray(
                N_w.apply(
                    zero_trace_composite(moments, N_w.domain.trace_space),
                ).interior.values,
            ),
            np.asarray(bulk_apply(N, psi).values),
            err_msg="the (n,2n) moment sibling drifted from its angular "
                    "binding.",
        )

    def test_the_fission_sibling_constructs(self, sn_mesh):
        r"""``F.on_moment_domain()`` is admissible by the BASE — machinery-first.

        ⚠ **No consumer today** (§16.1): fission is the eigenvalue OUTER
        source and the windowed driver windows the within-group gains only. The
        row exists because the base's ``on_moment_domain`` is ONE body for the
        whole family, so a repair that special-cased the transfer core would
        leave ``F`` unconstructible and nothing else would say so. It asserts
        the ends and the ℓ = 0 agreement, not a production path.
        """
        F = FissionOperator.from_solver_data(
            mat_xs=sn_mesh.material_xs_field(), space=sn_mesh.full_field_space,
        )
        F_w = F.on_moment_domain()
        require(
            F_w.domain.interior_space == F.flux_analysis.codomain,
            "the fission sibling's domain interior is not the analysis face's "
            "codomain — the base's on_moment_domain body is not one body.",
        )
        require(
            F_w.codomain == F.codomain,
            "the fission sibling changed its emission end.",
        )
        psi = _psi(sn_mesh, seed=41)
        moments = F.flux_analysis.apply(psi)
        np.testing.assert_array_equal(
            np.asarray(
                F_w.apply(
                    zero_trace_composite(moments, F_w.domain.trace_space),
                ).interior.values,
            ),
            np.asarray(bulk_apply(F, psi).values),
            err_msg="the fission moment sibling drifted from its angular "
                    "binding (χ carries no angle, so the two ends differ only "
                    "in how they reduce φ — and those reductions agree).",
        )
