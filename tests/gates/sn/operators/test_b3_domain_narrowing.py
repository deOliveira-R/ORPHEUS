r"""Campaign phase **B3.2** — what the domain narrowing itself claims.

B3.2 retyped the realized SN boundary law from ``full-face → full-face`` to
:math:`\Gamma_+ \to \Gamma_-`. The re-posed C-1 gates live next door in
``test_sn_boundary_operator.py`` (RG-1 … RG-5); this module carries the two
claims that are ABOUT the migration rather than about the composite:

1. **Bit-identity** (:class:`TestBitIdentityAgainstTheRetiredExpression`) —
   the composite ``SNBoundaryOperator.apply`` / ``.apply_transpose`` output is
   ``np.array_equal`` to what the PRE-B3.2 code produced, for
   ``{vacuum, reflective}`` on slab, sphere, cylinder and a 2-D Cartesian
   mesh. The reference is **materialised from the retired expression** — a
   full-face ``np.take(x, perm, 0)[inflow]`` for reflective, zeros for vacuum
   — never by calling the new code twice. Old-vs-new is necessary, not
   sufficient, so it is paired with the closed-form anchors already
   downstream (``keff == k_inf``, ``φ == Q/Σ_t``) which stay green.

2. **The domain gate** (:class:`TestEveryLawsDomain`) — a
   parametrised statement that a realized operator's domain is :math:`\Gamma_+`,
   over EVERY declaration the realizer table admits.

Wave 2 S2b re-pose: a law is now a DECLARATION (a deck, or a ``WallCondition`` on a tag), posed through a real slab, and the
realized operator is the ``problem.bc`` entry. The retired expression below is keyed on the DECLARED kind (a mirror deck or a unit
specular tag against a vacuum tag), read off the mesh's roles and the file's condition table, never off a realized type, so the
old-vs-new comparison stays structurally independent of the code under test. The bare-albedo and zero-flux refusal rows are
re-posed (see the notes at their sites).

Why the reference must not be the new code
------------------------------------------

A regression contract that says "the new code agrees with the new code" is a
tautology wearing a diff. The retired expression is transcribed here in four
lines of numpy against the face's DECLARED kind (the role on the mesh), so it shares
nothing with the narrowed realizer above the ``np.take`` line: a wrong
``local_perm``, a wrong write target, or a dropped ``α`` all move one side
only. ``np.array_equal`` and not ``nulp``: the narrowing removes ROWS from a
gather, it changes no reduction tree, so any drift at all is a real change
(``vv-principles`` §bit-identity, criterion 3 — zero drift is the predicted
bound, so the contract stays exact).
"""
from __future__ import annotations

from dataclasses import replace

import numpy as np
import pytest

from orpheus.geometry import CoordSystem, StructuredGeometry
from orpheus.mesh import CellsByCount, Mesh2D, Mesher
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.boundary.realizer import SN_REALIZERS
from orpheus.sn.problem import SNProblem
from orpheus.sn.operators.boundary import SNBoundaryOperator
from orpheus.transport.fields.angular_boundary_flux import AngularBoundaryFlux
from orpheus.transport.fields.angular_flux import AngularFlux
from orpheus.transport.timed_full_field import TimedFullField
from tests._harness.references import mirror_partner_indices
from tests.gates.sn._test_helpers import placeholder_materials
from orpheus.geometry.deck import MirrorDeck, TranslationDeck
from orpheus.transport.boundary_conditions import WallCondition
from tests.gates._exterior import conditions_for

#: Exterior conditions of the tags this file's meshes carry (wave 2: a face's role is a deck or an int tag; the condition a tag
#: stands for is declared at the solve, keyword ``boundary_conditions``).
_CONDITIONS = {0: WallCondition.vacuum(), 1: WallCondition.specular(1.0)}
_V, _R, _RC = 0, MirrorDeck(), 1  # vacuum tag; mirror deck; the reflecting tag of a CURVED body (a response, rule R2)


pytestmark = [pytest.mark.foundation]


# ── Fixtures: every geometry the SN boundary machinery reaches ────────


def _sn_1d(coord: CoordSystem, bcs: tuple, nx: int = 4, ng: int = 2) -> SNProblem:
    geom = StructuredGeometry(
        coord=coord,
        breakpoints=(0.0, 2.0),
        mat_ids=(0,),
        face_roles=bcs,
    )
    mesh = Mesher(geom).partition(CellsByCount.uniform_volume(nx)).mesh
    quad = (
        Quadrature.folded_product(n_mu=2, n_phi=4)
        if coord is CoordSystem.CYLINDRICAL
        else Quadrature.gauss_legendre(n_ordinates=4)
    )
    return SNProblem(mesh, quad, placeholder_materials(ng=ng), boundary_conditions=conditions_for(mesh, _CONDITIONS))


def _sn_2d(ng: int = 2) -> SNProblem:
    mesh = Mesh2D(
        edges_x=np.linspace(0.0, 1.0, 5), edges_y=np.linspace(0.0, 1.0, 4),
        mat_map=np.zeros((4, 3), dtype=int),   # (nx, ny) — deliberately nx != ny
        face_roles={
            "xmin": 0,
            "xmax": MirrorDeck(),
            "ymin": MirrorDeck(),
            "ymax": MirrorDeck(),
        },
    )
    return SNProblem(mesh, Quadrature.level_symmetric(4), placeholder_materials(ng=ng), boundary_conditions=conditions_for(mesh, _CONDITIONS))


#: The bit-identity fixture set. Both laws SN admits × every geometry, and the
#: 2-D case is deliberately MIXED (vacuum on xmin, reflective elsewhere) and
#: RECTANGULAR (nx ≠ ny) so a face↔face swap or an axis transpose cannot hide.
_FIXTURES = {
    "slab_vacuum_reflective": lambda: _sn_1d(CoordSystem.CARTESIAN, (_V, _R)),
    "slab_reflective_reflective": lambda: _sn_1d(CoordSystem.CARTESIAN, (_R, _R)),
    "sphere_reflective": lambda: _sn_1d(CoordSystem.SPHERICAL, (_RC,)),
    # The only fixture with TANGENTIAL ordinates (4 of 8 per face).
    "cyl_reflective": lambda: _sn_1d(CoordSystem.CYLINDRICAL, (_RC,)),
    "cart2d_mixed": _sn_2d,
}


def _random_state(sn: SNProblem, seed: int = 11) -> TimedFullField:
    rng = np.random.default_rng(seed)
    z = TimedFullField.zeros(
        interior=AngularFlux, boundary=AngularBoundaryFlux, space=sn.full_field_space,
    )
    return replace(
        z,
        boundary=replace(
            z.boundary,
            values=rng.uniform(0.5, 2.0, size=z.boundary.values.shape),
        ),
    )


# ── The retired expression, transcribed ───────────────────────────────


def _declared_axis(sn: SNProblem, face: str) -> "str | None":
    """The reflecting axis a face DECLARES, ``None`` for a vacuum face: a mirror deck, or the tag declared the unit specular wall.

    Read off the mesh's roles and this file's condition table, never off a realized operator (the independence the old reference had
    by reading the law descriptor)."""
    role = sn.mesh.face_roles[face]
    if isinstance(role, MirrorDeck):
        return face[0]
    return face[0] if _CONDITIONS[role] == WallCondition.specular(1.0) else None


def _pre_b32_face_action(axis: "str | None", face_in, inflow, quad, transpose: bool):
    r"""The PRE-B3.2 per-face action, in four lines of numpy.

    Verbatim semantics of ``_reflect_trace`` before the narrowing::

        full = law.apply(face_in)                  # full face -> full face
        out[sel] = full[sel]                       # outflow rows DISCARDED

        masked = zeros_like(face_in); masked[sel] = face_in[sel]
        out[...] = law.apply_transpose(masked)     # FULL image written

    with the realized laws' bodies inlined:

    * **vacuum** was ``IncomingOrdinateMaskTensor(inflow)`` — zero the inflow
      rows, keep the rest. Diagonal and symmetric, so its transpose is itself.
    * **reflective(α)** was ``α · PermutationOperator(reflection_index(axis))``
      — ``apply`` is ``np.take(x, perm, 0)``, ``apply_transpose`` the same
      gather through ``inverse_perm``. Transcribed here at α = 1, the only
      amplitude the mirror has had since the reflective cleanup (the
      transcription read ``float(law.albedo)`` until then).

    Nothing here imports a B3.2 symbol: the input is the law DESCRIPTOR, the
    arithmetic is numpy, and the partner map comes from the independent
    geometric reference (:func:`~tests._harness.references.mirror_partner_indices`
    — until G6.3 §7d.3 it was the shared ``reflection_index`` datum, which
    B3.2 did not touch and the §7d retirement replaced).
    """
    out = np.zeros_like(face_in)
    if transpose:
        masked = np.zeros_like(face_in)
        masked[inflow] = face_in[inflow]
        if axis is None:
            full_t = masked.copy()
            full_t[inflow] = 0.0            # the mask is its own transpose
        else:
            perm = mirror_partner_indices(quad, axis)
            inverse_perm = np.argsort(perm)  # NOT production's cached table
            full_t = np.take(masked, inverse_perm, axis=0)
        out[...] = full_t
        return out
    if axis is None:
        full = face_in.copy()
        full[inflow] = 0.0
    else:
        perm = mirror_partner_indices(quad, axis)
        full = np.take(face_in, perm, axis=0)
    out[inflow] = full[inflow]
    return out


class TestBitIdentityAgainstTheRetiredExpression:
    r"""``B.apply`` / ``B.apply_transpose`` are BYTE-identical across B3.2."""

    @pytest.mark.parametrize("case_id", list(_FIXTURES))
    @pytest.mark.parametrize("method", ["apply", "apply_transpose"])
    def test_composite_is_bit_identical_to_pre_b32(
        self, case_id: str, method: str,
    ) -> None:
        """The phase's headline claim, per geometry and per direction."""
        sn = _FIXTURES[case_id]()
        psi = _random_state(sn)
        got = getattr(SNBoundaryOperator(sn), method)(psi)
        transpose = method == "apply_transpose"
        moved = 0
        for face in sn.angular_trace.layout.faces:
            expected = _pre_b32_face_action(
                _declared_axis(sn, face),
                psi.boundary.face_view(face),
                sn.angular_trace.inflow_indices_for_face(face),
                sn.quad,
                transpose=transpose,
            )
            np.testing.assert_array_equal(
                got.boundary.face_view(face), expected,
                err_msg=(
                    f"{case_id} face {face!r} {method}: the narrowed "
                    f"composite does NOT reproduce the pre-B3.2 value. The "
                    f"reference is the RETIRED full-face expression, so a "
                    f"mismatch here is a genuine value regression, not a "
                    f"stale expectation."
                ),
            )
            moved += int(np.count_nonzero(expected))
        # NON-VACUITY: an all-zero comparison would pass for any broken
        # implementation. Vacuum contributes zero by design in BOTH
        # directions, so the sum is taken over faces and every fixture here
        # carries at least one reflective face.
        assert moved > 0, (
            f"{case_id} {method}: the pre-B3.2 reference is identically zero "
            f"on every face, so this row compares 0 == 0 and proves nothing."
        )

    @pytest.mark.parametrize("case_id", list(_FIXTURES))
    def test_the_retired_reference_is_falsifiable(self, case_id: str) -> None:
        """POSITIVE CONTROL for the gate above — prove the reference can
        disagree, so its agreement is evidence rather than a coincidence of
        two zero arrays or a reference that echoes the SUT.

        Perturbs the input on the OUTFLOW rows only and requires the composite
        to move: those rows are the narrowed law's entire domain, so a
        composite that ignored them (the "I dropped γ₊ and read the wrong
        half" family) would be caught here and nowhere else in this module.
        """
        sn = _FIXTURES[case_id]()
        psi = _random_state(sn)
        B = SNBoundaryOperator(sn)
        base = B.apply(psi).boundary.values.copy()
        bumped = replace(
            psi,
            boundary=replace(psi.boundary, values=psi.boundary.values.copy()),
        )
        touched = False
        for face in sn.angular_trace.layout.faces:
            if _declared_axis(sn, face) is None:
                continue  # vacuum's image is zero for ANY input, by design
            outflow = sn.angular_trace.outflow_indices_for_face(face)
            bumped.boundary.face_view(face)[outflow] += 1.0
            touched = True
        assert touched, f"{case_id}: no reflective face to perturb"
        assert not np.array_equal(base, B.apply(bumped).boundary.values), (
            f"{case_id}: perturbing Γ₊ left the composite unchanged — the "
            f"boundary action does not read its own domain."
        )


# ── The seven-law domain gate: B3.4's todo list ───────────────────────


# ⏹ The strict-xfail todo mechanism this registry carried is RETIRED
# (G6.3 step 7, 2026-08-07) — its set is EMPTY. The narrowing campaign it
# tracked: B3.2 narrowed {vacuum, reflective}; B3.4a white and
# prescribed_inflow; B3.4b albedo-by-closure (the bare spelling REFUSED,
# pinned below); and finally step 7 narrowed PERIODIC — the last row, whose
# ``_B34_XFAIL`` marker measured "silently accepts a Γ₊ input and echoes it
# back, Γ₊ → Γ₊, invisible to a shape check because |Γ₊| == |Γ₋| (vv
# Mode 12)". The realized periodic law is now the wrap motion's bound
# PermutationOperator Γ₊(f') → Γ₋(f), so it passes BOTH legs genuinely; the
# marker XPASSed(strict) the moment the step landed — the flip mechanism
# working exactly as designed — and retired with its last row.

#: Every declaration the realizer table admits on a Cartesian face, with a representative amplitude per kind: a ``WallCondition``
#: (carried on a tag) or a deck. Wave 2 S2b: the rows are DECLARATIONS posed through a real slab (their realized operator is the
#: ``problem.bc`` entry), where the earlier rows were law objects handed to the realizer. The bare albedo and zero flux of the old
#: inventory are not spellable here: see the two notes below.
_LAWS = {
    "vacuum": WallCondition.vacuum(),
    "reflective_a1": MirrorDeck(),
    "white_a1": WallCondition.white(1.0),
    "white_a03": WallCondition.white(0.3),
    "prescribed": WallCondition.uniform(1.0),
    # The specular wall at three amplitudes (the zero map, a scaled permutation, the unit one) and the diffuse wall of amplitude
    # one half: the rows that replaced the ``albedo_*`` strict xfails at B3.4b, kept as declarations.
    "albedo_specular_0": WallCondition.specular(0.0),
    "albedo_specular_1": WallCondition.specular(1.0),
    "albedo_specular_05": WallCondition.specular(0.5),
    "albedo_white_05": WallCondition.white(0.5),
    # Narrowed at G6.3 step 7 — the LAST law to narrow. The fixture face is "xmax", so the realized arrow reads the PARTNER's
    # Γ₊("xmin") — Leg A's probe is |Γ₊|-shaped either way because |Γ₊(xmin)| == |Γ₊(xmax)| on the symmetric GL rule.
    "periodic": TranslationDeck(),
}

_QUAD = Quadrature.gauss_legendre(8)


def _posed_xmax(declared) -> SNProblem:
    """A slab whose ``xmax`` declares ``declared`` (a deck, or a condition on the tag 0), ``xmin`` its partner (the same deck for a
    translation, a vacuum wall otherwise)."""
    conditions: dict[int, WallCondition] = {}
    if isinstance(declared, TranslationDeck):
        left = right = declared
    elif isinstance(declared, MirrorDeck):
        left, right = 0, declared
        conditions[0] = WallCondition.vacuum()
    else:
        left, right = 1, 0
        conditions = {0: declared, 1: WallCondition.vacuum()}
    mesh = Mesher(StructuredGeometry.slab((0.0, 2.0), (0,), left=left, right=right)).partition(CellsByCount.uniform_width(4)).mesh
    return SNProblem(mesh, _QUAD, placeholder_materials(ng=1), boundary_conditions=conditions)


class TestEveryLawsDomain:
    r"""A realized operator's domain is :math:`\Gamma_+` — for every admitted declaration."""

    @pytest.mark.parametrize("law_id", list(_LAWS))
    def test_realized_law_has_gamma_plus_domain(self, law_id: str) -> None:
        r"""Two legs, and the second is why this gate is not Mode-12 blind.

        **Leg A — it CONSUMES** :math:`\Gamma_+`: ``apply`` on a
        :math:`|\Gamma_+|`-row probe succeeds and returns
        :math:`|\Gamma_-|` rows.

        **Leg B — it is NOT an ENDOMORPHISM of the face**: fed a full-face
        probe it must not emit ``N`` rows. Leg A alone cannot tell
        :math:`\Gamma_+ \to \Gamma_-` from :math:`\Gamma_+ \to \Gamma_+`,
        because ``|Γ₊| == |Γ₋|`` on EVERY quadrature × face pair in the tree
        (measured: gauss_legendre 4/5/8, product 2×4/3×4/4×8, lebedev 9/17,
        level_symmetric 4/6). The error class sits inside the shape
        functional's invariance group — vv Mode 12 — so the discriminator has
        to leave that functional. Leg B is exactly what caught ``albedo`` and
        ``periodic`` while they were deferred: both passed Leg A by accident
        (a full-face echo has the right shape), and their xfail rows flipped
        strictly when B3.4b / G6.3 step 7 narrowed them.

        **Two outcomes pass Leg B, and they are not equal in strength.** A
        narrowed operator either

        * **REFUSES** the full-face input (it validates its own domain) — the
          strictly stronger answer, and the one the white arm gives, since
          ``AngularAverageOperator.apply`` checks ``psi.shape[0]``; or
        * **returns** :math:`|\Gamma_-|` rows, which is enough to prove it is
          not an endomorphism even though it did not notice the wrong input.
          The zero map answers this way: it ignores its input entirely, so it
          has nothing to validate against.

        Both refute "endomorphism", which is Leg B's claim. Whether an operator also
        *validates* is the separate RG-3b gate in
        ``test_sn_boundary_operator.py`` — kept separate on purpose, so that
        this row cannot quietly credit the weaker property as the stronger one.
        """
        sn = _posed_xmax(_LAWS[law_id])
        trace = sn.angular_trace
        n_in = trace.inflow_indices_for_face("xmax").size
        n_out = trace.outflow_indices_for_face("xmax").size
        # Fixture non-vacuity, OUTSIDE the claim: if the face were degenerate
        # (|Γ₋| == N) Leg B could never fail and the row would be theatre.
        assert 0 < n_in < _QUAD.n_ordinates and 0 < n_out < _QUAD.n_ordinates

        op = sn.bc["xmax"]
        image = op.apply(np.ones((n_out, 1)))
        assert image.shape == (n_in, 1), (
            f"{law_id}: realized operator maps ({n_out}, 1) -> {image.shape}; the "
            f"narrowed contract is Γ₊ -> Γ₋, i.e. ({n_in}, 1) (the trailing axis is the group count of the problem: a bound zero "
            f"map emits its DECLARED codomain, not an echo of the probe's trailing shape)."
        )
        # Narrow catches only: a refusal of the WRONG SHAPE is the outcome
        # under test, so anything else must surface as a real error rather
        # than be scored as a pass (vv Mode 8, class 6 — a broad ``except``
        # turns a dead gate green).
        try:
            full_rows: "int | None" = int(
                np.asarray(op.apply(np.ones((_QUAD.n_ordinates, 1)))).shape[0]
            )
        except (ValueError, IndexError):
            full_rows = None  # refused the input — the stronger outcome
        assert full_rows != _QUAD.n_ordinates, (
            f"{law_id}: realized operator emits {_QUAD.n_ordinates} rows for a full-face "
            f"input — it is an ENDOMORPHISM of the whole face slot, not a "
            f"Γ₊ -> Γ₋ map. (|Γ₊| == |Γ₋| here, so the shape leg above cannot "
            f"see this — vv Mode 12.)"
        )

    def test_a_negative_return_is_refused_on_the_state_cone_axis(self) -> None:
        """Zero flux is not a deferred narrowing: SN has no transport realization for a negative angular inflow, and says so.

        Wave 2 S2b re-pose of ``test_zero_flux_is_refused_outright`` (which asserted ``BoundaryError.law == "zero_flux"``, the law
        key of the retired string registry and its ordering tripwire against the bare-albedo refusal: neither exists). The refusal
        is now the realizer table's UNREGISTERED-UNIT door: a condition with ``return_fraction < 0`` on a table whose state cone is
        ``psi >= 0`` is refused on the STATE_CONE axis, and SN does so without importing diffusion's ``ZeroFluxCondition``, which
        the probe stands for (``return_fraction == -1.0``, the value zero flux has).
        """
        from orpheus.transport.boundary_conditions import BoundaryRefused, RefusalAxis
        from tests.gates.transport._boundary_conditions_helpers import local

        mesh = Mesher(StructuredGeometry.slab((0.0, 2.0), (0,), left=0, right=1)).partition(CellsByCount.uniform_width(4)).mesh
        with pytest.raises(BoundaryRefused) as refused:
            SNProblem(
                mesh, _QUAD, placeholder_materials(ng=1),
                boundary_conditions={0: WallCondition.vacuum(), 1: local().NegativeProbe()},
            )
        assert refused.value.admission.axis is RefusalAxis.STATE_CONE

    # ``test_bare_albedo_is_refused_outright`` (B3.4b) is RETIRED with its subject, not re-posed. A bare albedo was refused because
    # :math:`R = \\alpha I` with :math:`G = \\mathrm{id}` names no pairing between the two half-traces. ``WallCondition`` is shaped
    # by its response network (``alpha * SpecularChain``, ``alpha * DiffuseChain``), so a response without a return closure has no
    # spelling and there is nothing to refuse: a test asserting its absence would be a gate no input can redden (X1). The three
    # amplitudes the old row parametrised over (the three pre-B3.4b branches) are the three ``specular`` rows of ``_LAWS``.

    def test_the_law_inventory_is_complete(self) -> None:
        """The gate covers EVERY unit type the SN realizer table registers.

        Without this, a unit registered tomorrow would silently escape the domain gate — the failure mode where a matrix looks
        exhaustive because nothing checks that it is. (Re-posed from the ``BoundaryTraceLaw.registry`` census: the table's own
        ``catalogue`` is the inventory now.)
        """
        from orpheus.numerics.operator import OperatorSum  # noqa: F401  (documentation of the sum case, not a row)
        from orpheus.transport.boundary_conditions import DiffuseChain, InflowTerms, SpecularChain, ZeroNetwork

        covered = {
            ZeroNetwork: WallCondition.vacuum(),
            SpecularChain: WallCondition.specular(0.5),
            DiffuseChain: WallCondition.white(0.5),
            InflowTerms: WallCondition.uniform(1.0),
        }
        registered = set(SN_REALIZERS.catalogue())
        missing = registered - set(covered)
        assert not missing, (
            f"unit type(s) {sorted(m.__name__ for m in missing)} are registered in SN_REALIZERS but absent from the B3.2 domain "
            f"gate — add a row to _LAWS and to this map."
        )
        for unit, declared in covered.items():
            assert declared in _LAWS.values(), f"{unit.__name__}: its representative {declared!r} is not a row of _LAWS"
