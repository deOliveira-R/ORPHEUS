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

2. **The seven-law domain gate** (:class:`TestEveryLawsDomain`) — a
   parametrised statement that a realized law's domain is :math:`\Gamma_+`,
   over EVERY law in the registry. The laws B3.4 has yet to narrow ship as
   ``xfail(strict=True)`` naming it: the marker set IS the todo list, so B3.4
   cannot land silently without deleting them.

Why the reference must not be the new code
------------------------------------------

A regression contract that says "the new code agrees with the new code" is a
tautology wearing a diff. The retired expression is transcribed here in four
lines of numpy against the law DESCRIPTOR (``sn.bc[face].law``), so it shares
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

from orpheus.geometry import (
    BC, Mesh1D, Mesh2D, Region, RegionMesh, StructuredGeometry,
)
from orpheus.geometry.boundary import (
    AlbedoBoundary,
    ConstantInflowSource,
    PeriodicBoundary,
    PrescribedInflow,
    ReflectiveBoundary,
    VacuumInflow,
    WhiteBoundary,
    ZeroFluxBoundary,
)
from orpheus.geometry.boundary._errors import BoundaryError
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.boundary.realizer import SNBoundaryRealizer
from orpheus.sn.mesh.augmented_mesh import SNMesh
from orpheus.sn.operators.boundary import SNBoundaryOperator
from orpheus.transport.fields.angular_boundary_flux import AngularBoundaryFlux
from orpheus.transport.fields.angular_flux import AngularFlux
from orpheus.transport.timed_full_field import TimedFullField
from tests.sn._test_helpers import face_method_space, placeholder_materials

pytestmark = [pytest.mark.foundation]


# ── Fixtures: every geometry the SN boundary machinery reaches ────────


def _sn_1d(geometry: str, bcs: tuple, nx: int = 4, ng: int = 2) -> SNMesh:
    geom = StructuredGeometry(
        geometry=geometry,
        regions=(Region(mat_id=0, outer_thickness_cm=2.0),),
        bcs=bcs,
    )
    mesh = Mesh1D.from_geometry(geom, region_meshes=(RegionMesh(n_cells=nx),))
    quad = (
        Quadrature.product(n_mu=2, n_phi=4)
        if geometry == "CYL"
        else Quadrature.gauss_legendre(n_ordinates=4)
    )
    return SNMesh(mesh, quad, placeholder_materials(ng=ng))


def _sn_2d(ng: int = 2) -> SNMesh:
    mesh = Mesh2D(
        edges_x=np.linspace(0.0, 1.0, 5), edges_y=np.linspace(0.0, 1.0, 4),
        mat_map=np.zeros((4, 3), dtype=int),   # (nx, ny) — deliberately nx != ny
        bc_xmin=BC.vacuum, bc_xmax=BC.reflective,
        bc_ymin=BC.reflective, bc_ymax=BC.reflective,
    )
    return SNMesh(mesh, Quadrature.level_symmetric(4), placeholder_materials(ng=ng))


#: The bit-identity fixture set. Both laws SN admits × every geometry, and the
#: 2-D case is deliberately MIXED (vacuum on xmin, reflective elsewhere) and
#: RECTANGULAR (nx ≠ ny) so a face↔face swap or an axis transpose cannot hide.
_FIXTURES = {
    "slab_vacuum_reflective": lambda: _sn_1d("SLB", (BC.vacuum, BC.reflective)),
    "slab_reflective_reflective": lambda: _sn_1d("SLB", (BC.reflective, BC.reflective)),
    "sphere_reflective": lambda: _sn_1d("SPH", (BC.reflective,)),
    # The only fixture with TANGENTIAL ordinates (4 of 8 per face).
    "cyl_reflective": lambda: _sn_1d("CYL", (BC.reflective,)),
    "cart2d_mixed": _sn_2d,
}


def _random_state(sn: SNMesh, seed: int = 11) -> TimedFullField:
    rng = np.random.default_rng(seed)
    z = TimedFullField.zeros(
        interior=AngularFlux, boundary=AngularBoundaryFlux, mesh=sn,
    )
    return replace(
        z,
        boundary=replace(
            z.boundary,
            values=rng.uniform(0.5, 2.0, size=z.boundary.values.shape),
        ),
    )


# ── The retired expression, transcribed ───────────────────────────────


def _pre_b32_face_action(law, face_in, inflow, quad, transpose: bool):
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
      gather through ``inverse_perm``.

    Nothing here imports a B3.2 symbol: the input is the law DESCRIPTOR, the
    arithmetic is numpy, and the only ORPHEUS datum is ``reflection_index``
    (which B3.2 did not touch).
    """
    out = np.zeros_like(face_in)
    if transpose:
        masked = np.zeros_like(face_in)
        masked[inflow] = face_in[inflow]
        if isinstance(law, VacuumInflow):
            full_t = masked.copy()
            full_t[inflow] = 0.0            # the mask is its own transpose
        elif isinstance(law, ReflectiveBoundary):
            perm = quad.reflection_index(law.axis)
            inverse_perm = np.argsort(perm)  # NOT production's cached table
            full_t = float(law.albedo) * np.take(masked, inverse_perm, axis=0)
        else:  # pragma: no cover - SN admits only the two
            raise AssertionError(f"unreachable law {type(law).__name__}")
        out[...] = full_t
        return out
    if isinstance(law, VacuumInflow):
        full = face_in.copy()
        full[inflow] = 0.0
    elif isinstance(law, ReflectiveBoundary):
        perm = quad.reflection_index(law.axis)
        full = float(law.albedo) * np.take(face_in, perm, axis=0)
    else:  # pragma: no cover
        raise AssertionError(f"unreachable law {type(law).__name__}")
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
                sn.bc[face].law,
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
            if isinstance(sn.bc[face].law, VacuumInflow):
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


_B34_XFAIL = (
    "B3.4b/c — this law still realizes FULL-FACE. B3.2 narrowed the two laws "
    "SN reaches from a mesh ({vacuum, reflective}); B3.4a narrowed white and "
    "prescribed_inflow. What remains is albedo — at EVERY α, since the α=0 "
    "ZeroOperator and α=1 IdentityOperator fast paths are endomorphisms too, "
    "so this is NOT limited to α ∉ {0,1} — and periodic. MEASURED: both "
    "silently accept a Γ₊ input and echo it back, i.e. Γ₊ → Γ₊, the wrong "
    "codomain, invisible to a shape check because |Γ₊| == |Γ₋| on every "
    "fixture (vv Mode 12). Each is blocked on a DESIGN ruling, not on "
    "plumbing: albedo is under-determined on an angular trace (R = α·I is a "
    "Γ₊ → Γ₊ endomorphism and its G = IdentityMap supplies no crossing), so "
    "B3.4b gives it an explicit re-emission closure carried in R; periodic's "
    "G reads the PARTNER face's Γ₊, which B3.4c builds (#183, #189). Delete "
    "this marker when the law narrows."
)

#: Every law in the registry, with a representative amplitude per kind.
#: ``zero_flux`` is absent by construction — SN refuses it outright, which is
#: pinned as its own negative below rather than as an xfail (a structural
#: refusal is not a deferred narrowing).
_LAWS = {
    "vacuum": (VacuumInflow(), False),
    "reflective_a1": (ReflectiveBoundary(axis="x", albedo=1.0), False),
    "reflective_a07": (ReflectiveBoundary(axis="x", albedo=0.7), False),
    # B3.4a narrowed these two. The law's axis / outward_sign must match the
    # fixture face ("xmax" ⇔ x, +1) or the realizer's orientation cross-check
    # fires — which is the point of that guard, not an inconvenience.
    "white_a1": (WhiteBoundary(axis="x", outward_sign=+1, albedo=1.0), False),
    "white_a03": (WhiteBoundary(axis="x", outward_sign=+1, albedo=0.3), False),
    "prescribed": (PrescribedInflow(source=ConstantInflowSource(value=1.0)), False),
    "albedo_0": (AlbedoBoundary(albedo=0.0), True),
    "albedo_1": (AlbedoBoundary(albedo=1.0), True),
    "albedo_05": (AlbedoBoundary(albedo=0.5), True),
    "periodic": (PeriodicBoundary(), True),
}


class TestEveryLawsDomain:
    r"""A realized law's domain is :math:`\Gamma_+` — for every law."""

    @pytest.mark.parametrize("law_id", list(_LAWS))
    def test_realized_law_has_gamma_plus_domain(
        self, law_id: str, request: pytest.FixtureRequest,
    ) -> None:
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
        to leave that functional. Leg B is exactly what catches ``albedo`` and
        ``periodic``, which pass Leg A by accident.

        **Two outcomes pass Leg B, and they are not equal in strength.** A
        narrowed law either

        * **REFUSES** the full-face input (it validates its own domain) — the
          strictly stronger answer, and the one B3.4a's white arm gives, since
          ``AngularAverageOperator.apply`` checks ``psi.shape[0]``; or
        * **returns** :math:`|\Gamma_-|` rows, which is enough to prove it is
          not an endomorphism even though it did not notice the wrong input.
          B3.4a's prescribed-inflow arm answers this way: it ignores its input
          entirely, so it has nothing to validate against.

        Both refute "endomorphism", which is Leg B's claim. Whether a law also
        *validates* is the separate RG-3b gate in
        ``test_sn_boundary_operator.py`` — kept separate on purpose, so that
        this row cannot quietly credit the weaker property as the stronger one.
        """
        law, deferred = _LAWS[law_id]
        if deferred:
            request.node.add_marker(pytest.mark.xfail(strict=True, reason=_B34_XFAIL))
        quad = Quadrature.gauss_legendre(8)
        space = face_method_space(quad, face="xmax")
        n_in, n_out = space.inflow_indices.size, space.outflow_indices.size
        # Fixture non-vacuity, OUTSIDE the claim: if the face were degenerate
        # (|Γ₋| == N) Leg B could never fail and the row would be theatre.
        assert 0 < n_in < quad.N and 0 < n_out < quad.N

        op = SNBoundaryRealizer().realize(law, space)
        image = op.apply(np.ones((n_out, 3)))
        assert image.shape == (n_in, 3), (
            f"{law_id}: realized law maps ({n_out}, 3) -> {image.shape}; the "
            f"narrowed contract is Γ₊ -> Γ₋, i.e. ({n_in}, 3)."
        )
        # Narrow catches only: a refusal of the WRONG SHAPE is the outcome
        # under test, so anything else must surface as a real error rather
        # than be scored as a pass (vv Mode 8, class 6 — a broad ``except``
        # turns a dead gate green).
        try:
            full_rows: "int | None" = int(
                np.asarray(op.apply(np.ones((quad.N, 3)))).shape[0]
            )
        except (ValueError, IndexError):
            full_rows = None  # refused the input — the stronger outcome
        assert full_rows != quad.N, (
            f"{law_id}: realized law emits {quad.N} rows for a full-face "
            f"input — it is an ENDOMORPHISM of the whole face slot, not a "
            f"Γ₊ -> Γ₋ map. (|Γ₊| == |Γ₋| here, so the shape leg above cannot "
            f"see this — vv Mode 12.)"
        )

    def test_zero_flux_is_refused_outright(self) -> None:
        """``zero_flux`` is not a deferred narrowing — SN has no transport
        realization for a negative angular inflow, and says so.

        Pinned here so the parametrised gate's law inventory is provably
        COMPLETE: registry = the ten rows above ⊎ ``zero_flux``.
        """
        quad = Quadrature.gauss_legendre(8)
        with pytest.raises(BoundaryError) as exc:
            SNBoundaryRealizer().realize(
                ZeroFluxBoundary(), face_method_space(quad, face="xmax"),
            )
        assert exc.value.law == "zero_flux"

    def test_the_law_inventory_is_complete(self) -> None:
        """The gate covers EVERY registered law kind.

        Without this, a law added to the registry tomorrow would silently
        escape the domain gate — the failure mode where a matrix looks
        exhaustive because nothing checks that it is.
        """
        from orpheus.geometry.boundary import BoundaryTraceLaw

        covered = {type(law).__name__ for law, _ in _LAWS.values()}
        covered.add("ZeroFluxBoundary")
        registered = {
            type(BoundaryTraceLaw.create(key)).__name__
            for key in BoundaryTraceLaw.registry
        }
        missing = registered - covered
        assert not missing, (
            f"law kind(s) {sorted(missing)} are registered but absent from "
            f"the B3.2 domain gate — add them to _LAWS (narrowed) or to the "
            f"zero_flux-style refusal negatives."
        )
