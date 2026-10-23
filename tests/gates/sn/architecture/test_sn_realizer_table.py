r"""Gates of the SN realizer table ``SN_REALIZERS`` of ``orpheus/sn/boundary/realizer.py`` (wave 2, slice S2b).

Spec: ``.claude/plans/wave2_boundary_carve.md`` section 2.6 ("The SN realizer table"), section 3 (laws L-A, L-D, L-E, L-L)
and section 6.2 (S2b gates 4, 5, 7 (the table half), 9 and 13). The table is read through the same door a problem uses,
``SN_REALIZERS.admit`` and ``SN_REALIZERS.realize(declared, site, space)``, on the method space of a face, so these gates
need no ``SNProblem`` and pose the realizer on the S16 slab face the laws name.

First reds, at the base tip of slice S2b: ``ModuleNotFoundError``/``ImportError`` from the first read of
``orpheus.sn.boundary.realizer.SN_REALIZERS`` (the module still exports ``SNBoundaryRealizer``; ``import orpheus.sn``
itself raises because ``orpheus.transport.method`` imports the retired ``BC``).

The ladder::

    test_the_table_registers_...                          foundation: the five realizers, sealed
    test_the_unregistered_door_names_the_state_cone_...   foundation: the door of ZeroFluxCondition (STATE_CONE), 659
        ^
    test_a_certification_failure_names_its_axis_at_the_table   edge: ordinate_pairing's refusal is an admission
    test_return_fraction_is_the_declared_...               edge: the realized fact IS the declared wall's
        ^
    test_mirror_naturality                                 interior: L-E on S16, every shared unit
    test_the_realizer_is_a_homomorphism                    interior: L-L, exact on parameters, ulp on matrices
    test_the_left_eigenvector_rider_through_the_realizer   interior: L-A (iv) through the real realizer

Mutations for the reviewers (each on a COPY of ``realizer.py``/``realized_operators.py``, the verdict per arm):

* M-report-one: ``SNSpecularChain`` builds ``alpha = 1.0`` whatever the weight: reddens
  ``test_return_fraction_is_the_declared_...`` (declared 0.7, built 1.0), the rider, and the ERR-094 gates
  (``test_dsa_low_order.py``, ``test_keff_estimator_gate.py``).
* M-orientation-blind (law L-E): the hemisphere of the Lambertian taken from the FIRST face the realizer was built
  for: accepted by the realizer, reddens ``test_mirror_naturality`` on the second face.
* M-rescale (law L-L): the diffuse image of a wall rescaled by the total ``rho = s + d``: inert on ``(0.3, 0.7)``
  (``rho = 1``), reddens on ``(0.3, 0.5)``.
* M-diffuse-admits-radial: ``SNDiffuseChain.admit`` returns ``Realized`` on a radial outer site: reddens the SN rows of
  ``test_admission_matrix.py`` and the construction refusal of ``test_sn_boundary_binding.py``.
* M-door-citing-654: the table's ``unregistered_issue`` is 654: reddens ``test_the_unregistered_door_...``.
"""

from __future__ import annotations

import math
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

from tests.gates.sn.architecture import _realized_boundary_helpers as H
from tests.gates.transport._boundary_conditions_helpers import local

pytestmark = pytest.mark.l1

_THIS = "tests/gates/sn/architecture/test_sn_realizer_table.py"
_FACES = ("xmin", "xmax")

#: The walls the shared family realizes on a Cartesian face, as ``(label, builder)``: the three units and a mix.
_WALLS = {
    "vacuum": lambda w: w.vacuum(),
    "specular 0.7": lambda w: w.specular(0.7),
    "specular 1.0": lambda w: w.specular(1.0),
    "white 0.3": lambda w: w.white(0.3),
    "white 1.0": lambda w: w.white(1.0),
    "both 0.3 0.5": lambda w: w.from_amplitudes(0.3, 0.5),
}


def _realize(declared, quad_name: str, face: str, trace=None):
    """``SN_REALIZERS.realize`` of ``declared`` on the Cartesian ``face`` of the slab trace of ``quad_name``."""
    quad = H.QUADRATURES[quad_name]()
    trace = trace if trace is not None else H.trace_of(quad)
    return H.sn_table().realize(declared, H.site(face), H.method_space(quad, face, trace)), quad, trace


# ═════════════════════════════════════════════════════════════════════════════
# The table itself
# ═════════════════════════════════════════════════════════════════════════════


@pytest.mark.foundation
def test_the_table_registers_the_five_realizers_and_is_sealed() -> None:
    r"""``SN_REALIZERS`` is a sealed ``RealizerTable`` of method ``"sn"`` on the non-negative-flux cone whose
    unregistered-unit door cites #659; it is filled by ``SNZero`` (``ZeroNetwork``), ``SNSpecularChain``,
    ``SNDiffuseChain``, ``SNSource`` (a ``SourceRealizer`` with NO ``realize``: a source is an array, never an operator)
    and ``SNDecks`` (the one ``DeckRealizer``). The catalogue lists the three units at ``Realized(True)``.

    Red input: the base tip exports ``SNBoundaryRealizer`` and no table (the import of ``SN_REALIZERS`` fails).
    Mutation: a table sealed before ``SNDecks`` registers; ``SNSource`` given a ``realize``.
    """
    import orpheus.sn.boundary.realizer as realizer

    B = H.bc()
    table = H.sn_table()
    assert table.method == "sn" and table.state_cone is B.StateCone.NONNEGATIVE_FLUX and table.unregistered_issue == 659
    for name, base in (("SNZero", B.ResponseRealizer), ("SNSpecularChain", B.ResponseRealizer), ("SNDiffuseChain", B.ResponseRealizer), ("SNSource", B.SourceRealizer), ("SNDecks", B.DeckRealizer)):
        assert issubclass(getattr(realizer, name), base), f"{name} is not a {base.__name__}"
    assert not hasattr(realizer.SNSource, "realize"), "a source realizer has no realize (decision 53)"
    assert not hasattr(realizer, "SNBoundaryRealizer"), "the 400-line isinstance ladder must be gone"
    catalogue = table.catalogue()
    for unit in (B.SpecularChain, B.DiffuseChain, B.ZeroNetwork):
        assert catalogue[unit] == B.Realized(True), f"{unit.__name__}: {catalogue.get(unit)!r}"


@pytest.mark.foundation
@pytest.mark.rests_on(f"{_THIS}::test_the_table_registers_the_five_realizers_and_is_sealed")
def test_the_unregistered_door_names_the_state_cone_without_importing_diffusion() -> None:
    r"""SN refuses a condition whose return fraction is negative (``ZeroFluxCondition``'s ``-1``) on the STATE_CONE
    axis, through the table's unregistered-unit door; a non-negative method-only condition with no entry is
    ``Refused(PLUMBING, #659)`` (the table's own issue, not CP's #654). A test-local condition with
    ``return_fraction == -1`` stands for the zero-flux type (it may not be imported: ``diffusion`` is S3), and a FRESH
    interpreter that imports the SN realizer has not imported ``orpheus.diffusion`` (the state cone keeps its axis
    without naming the diffusion type).

    Red input: the base tip has no table. Mutation M-door-citing-654.
    """
    B = H.bc()
    probes = local()
    table = H.sn_table()
    where = H.site("xmax")
    negative = table.admit(probes.NegativeProbe(), where)
    assert isinstance(negative, B.Refused) and negative.axis is B.RefusalAxis.STATE_CONE
    plain = table.admit(probes.Probe(0.5), where)
    assert isinstance(plain, B.Refused) and plain.axis is B.RefusalAxis.PLUMBING and plain.issue == 659
    with pytest.raises(B.BoundaryRefused) as info:
        table.require(probes.NegativeProbe(), where)
    assert info.value.admission.axis is B.RefusalAxis.STATE_CONE

    root = Path(os.path.abspath(__import__("orpheus").__file__)).parent.parent  # not resolve(): a symlinked tree must stay itself
    code = "import sys, orpheus.sn.boundary.realizer; sys.exit(1 if any(m == 'orpheus.diffusion' or m.startswith('orpheus.diffusion.') for m in sys.modules) else 0)"
    done = subprocess.run([sys.executable, "-O", "-c", code], cwd="/", env={"PYTHONPATH": str(root), "PATH": "/usr/bin:/bin"}, capture_output=True, text=True)
    assert done.returncode == 0, f"importing the SN realizer imports diffusion (or fails): {done.stderr[-300:]}"


@pytest.mark.foundation
@pytest.mark.rests_on(f"{_THIS}::test_the_table_registers_the_five_realizers_and_is_sealed")
def test_a_certification_failure_names_its_axis_at_the_table() -> None:
    r"""With a space the table asks ``ordinate_pairing``: an odd-``n_phi`` product rule has no x-mirror closure, so
    ``MirrorDeck`` and ``specular`` on an x-face are ``Refused(ANGULAR_RESOLUTION)`` and ``realize`` raises
    ``BoundaryRefused`` carrying that admission; without a space (the front end) the same declarations are
    ``Realized(True)`` (the pin of ``test_admission_matrix.py`` is site-level), and an even-``n_phi`` rule certifies
    (the control). The diffuse chain asks no pairing and is admitted on the odd rule.

    Red input: the base tip has no table. Mutation: an ``admit`` that skips the certification (the refusal moves to a
    ``realize`` that raises a bare ``OrdinatePairingError``, not a ``BoundaryRefused``).
    """
    from orpheus.geometry.deck import MirrorDeck

    B = H.bc()
    table = H.sn_table()
    odd = H.quadrature_class().product(n_mu=2, n_phi=3)
    even = H.quadrature_class().product(n_mu=2, n_phi=4)
    where = H.site("xmin", dimension=2)
    for quad, certifies in ((odd, False), (even, True)):
        space = H.method_space(quad, "xmin", H.trace_of(quad, ("xmin", "xmax", "ymin", "ymax")))
        for declared in (MirrorDeck(), B.WallCondition.specular(0.5)):
            admission = table.admit(declared, where, space)
            if certifies:
                assert admission == B.Realized(True), f"{declared!r} on the even rule: {admission!r}"
            else:
                assert isinstance(admission, B.Refused) and admission.axis is B.RefusalAxis.ANGULAR_RESOLUTION, f"{declared!r}: {admission!r}"
                with pytest.raises(B.BoundaryRefused) as info:
                    table.realize(declared, where, space)
                assert info.value.admission.axis is B.RefusalAxis.ANGULAR_RESOLUTION
            assert table.admit(declared, where) == B.Realized(True), "the front end (no space) is site-level"
    odd_space = H.method_space(odd, "xmin", H.trace_of(odd, ("xmin", "xmax", "ymin", "ymax")))
    assert table.admit(B.WallCondition.white(0.5), where, odd_space) == B.Realized(True)


# ═════════════════════════════════════════════════════════════════════════════
# Gate 9: the realized fact is the declared wall's, measured
# ═════════════════════════════════════════════════════════════════════════════


@pytest.mark.l1
@pytest.mark.catches("ERR-094")
@pytest.mark.rests_on(f"{_THIS}::test_the_table_registers_the_five_realizers_and_is_sealed")
@pytest.mark.parametrize("label", list(_WALLS))
@pytest.mark.parametrize("quad_name", ("gl8", "product(2,4)", "ls6", "lebedev17"))
@pytest.mark.parametrize("face", _FACES)
def test_return_fraction_is_the_declared_and_the_measured(label: str, quad_name: str, face: str) -> None:
    r"""G14 through the REAL realizer (ERR-094's collapse trigger made structural): for each wall of the shared family
    the realized face's ``return_fraction`` equals the declared wall's ``return_fraction`` (``math.fsum`` of its
    weights), and the measured ``J^-/J^+`` of the realized operator on a positive outflow equals both to ``1e-12``.
    The partial specular reflector is the row that matters: it permutes ordinates (so a mirror shortcut admits it) and
    returns ``0.7`` of what arrives.

    Red input: the base tip has no table. Mutation M-report-one (``SNSpecularChain`` builds ``alpha = 1.0``): reddens
    ``specular 0.7``, ``both 0.3 0.5`` and the ERR-094 gates of the DSA and the keff estimator.
    """
    B = H.bc()
    declared = _WALLS[label](B.WallCondition)
    op, quad, trace = _realize(declared, quad_name, face)
    assert op.return_fraction == declared.return_fraction
    m_out, m_in = H.trace_weights(trace, face)
    psi = H.positive_outflow(m_out.size)
    j_plus = float(np.sum(m_out * psi))
    j_minus = float(np.sum(m_in * np.asarray(op.apply(psi))))
    assert abs(j_minus / j_plus - declared.return_fraction) <= 1e-12, f"{label} on {quad_name}/{face}: measured {j_minus / j_plus!r}, declared {declared.return_fraction!r}"
    assert op.is_vacuum is (declared.return_fraction == 0.0)


@pytest.mark.l1
@pytest.mark.rests_on(f"{_THIS}::test_return_fraction_is_the_declared_and_the_measured")
@pytest.mark.parametrize("quad_name", ("gl8", "product(2,4)", "ls6", "lebedev17"))
@pytest.mark.parametrize("face", _FACES)
def test_the_left_eigenvector_rider_through_the_realizer(quad_name: str, face: str) -> None:
    r"""Law L-A (iv) through the real realizer, as the owed rider gate: ``op.apply_transpose(w|mu| on Gamma_-) ==
    rho * (w|mu| on Gamma_+)``, ``array_equal`` for the specular wall and ``ceil(N/3)`` ulp for the Lambertian wall,
    ``rho`` the DECLARED return fraction. The array identity is stronger than the one-vector ratio: an image scaled by
    ``alpha**2`` (M-alpha-squared) fails it where a ratio taken at ``alpha = 1`` cannot.

    Red input: the base tip has no table.
    """
    B = H.bc()
    for declared, exact in ((B.WallCondition.specular(0.7), True), (B.WallCondition.white(0.3), False)):
        op, _, trace = _realize(declared, quad_name, face)
        m_out, m_in = H.trace_weights(trace, face)
        got = np.asarray(op.apply_transpose(m_in))
        want = declared.return_fraction * m_out
        if exact:
            np.testing.assert_array_equal(got, want)
        else:
            ulps = np.max(np.abs(got - want) / np.spacing(np.abs(want)))
            assert ulps <= math.ceil(m_out.size / 3), f"{declared!r}: rider off by {ulps} ulp"


# ═════════════════════════════════════════════════════════════════════════════
# Gate 4: law L-E, the mirror is a natural isomorphism
# ═════════════════════════════════════════════════════════════════════════════


def _sigma(quad, trace) -> tuple[np.ndarray, np.ndarray]:
    r"""The slab mirror :math:`\sigma` between the two faces as permutation matrices ``(S_+, S_-)``: ``S_+`` carries
    ``Gamma_+(xmin)`` onto ``Gamma_+(xmax)`` and ``S_-`` carries ``Gamma_-(xmin)`` onto ``Gamma_-(xmax)``, built by
    matching direction cosines (``mu -> -mu``), not from any realizer."""
    mu = np.asarray(quad.mean_axis_cosine(0), dtype=float)

    def carry(src: np.ndarray, dst: np.ndarray) -> np.ndarray:
        p = np.zeros((dst.size, src.size))
        for j, g in enumerate(src):
            k = int(np.argmin(np.abs(mu[dst] + mu[g])))
            assert abs(mu[dst[k]] + mu[g]) < 1e-14
            p[k, j] = 1.0
        return p

    return (
        carry(np.asarray(trace.outflow_indices_for_face("xmin")), np.asarray(trace.outflow_indices_for_face("xmax"))),
        carry(np.asarray(trace.inflow_indices_for_face("xmin")), np.asarray(trace.inflow_indices_for_face("xmax"))),
    )


@pytest.mark.l1
@pytest.mark.rests_on(f"{_THIS}::test_return_fraction_is_the_declared_and_the_measured")
@pytest.mark.parametrize("label", [k for k in _WALLS if k != "vacuum"])
def test_mirror_naturality(label: str) -> None:
    r"""Law L-E: for the slab mirror :math:`\sigma`, :math:`S_-\,\mathrm{real}(\beta, {\rm xmin})\,S_+^{-1} =
    \mathrm{real}(\beta, {\rm xmax})` (``array_equal``, S16), two realizations compared through a THIRD object (the
    reversal built from direction cosines), for every shared unit and the mix, and for the mirror DECK. It catches the
    wrong-hemisphere class without reading ``face_normal``: the retired B3.4a tests' replacement witness.

    Red input: the base tip has no table. The law's first red is a MUTATION only (the tree refuses the wrong-hemisphere
    input at ``_checked_angular_average``): M-orientation-blind, a realizer taking the hemisphere from the first face
    it was built for, accepted by the realizer and reddening this row on ``xmax``.
    """
    B = H.bc()
    quad = H.QUADRATURES["gl16"]()
    trace = H.trace_of(quad)
    s_plus, s_minus = _sigma(quad, trace)
    declared = _WALLS[label](B.WallCondition)
    mats = {}
    for face in _FACES:
        op, _, _ = _realize(declared, "gl16", face, trace)
        mats[face] = H.dense(op, H.n_outflow(trace, face))
    np.testing.assert_array_equal(s_minus @ mats["xmin"] @ s_plus.T, mats["xmax"])


@pytest.mark.l1
@pytest.mark.rests_on(f"{_THIS}::test_mirror_naturality")
def test_mirror_naturality_holds_for_the_mirror_deck() -> None:
    r"""The same law for the deck: ``MirrorDeck`` realizes through ``SNDecks`` to the unit specular image on each face
    and the naturality holds. Red input: the base tip has no table."""
    from orpheus.geometry.deck import MirrorDeck

    quad = H.QUADRATURES["gl16"]()
    trace = H.trace_of(quad)
    s_plus, s_minus = _sigma(quad, trace)
    mats = {face: H.dense(_realize(MirrorDeck(), "gl16", face, trace)[0], H.n_outflow(trace, face)) for face in _FACES}
    np.testing.assert_array_equal(s_minus @ mats["xmin"] @ s_plus.T, mats["xmax"])


# ═════════════════════════════════════════════════════════════════════════════
# Gate 13: law L-L, the realizer is a homomorphism of the algebra
# ═════════════════════════════════════════════════════════════════════════════

#: ``(s, d)`` points: the carry witness ``mixed_30spec_70white`` (``s + d = 1``: inert to a rho-rescale) and the
#: DISCRIMINATING fixture ``0.3 S + 0.5 D`` (``s + d = 0.8``), plus a dyadic and a non-dyadic point.
_POINTS = ((0.3, 0.7), (0.3, 0.5), (0.25, 0.5), (0.1, 0.3), (0.7, 0.2))


def _ulps(a: np.ndarray, b: np.ndarray) -> float:
    spacing = np.spacing(np.maximum(np.abs(a), np.abs(b)))
    return float(np.max(np.where(spacing > 0, np.abs(a - b) / np.where(spacing > 0, spacing, 1.0), np.abs(a - b))))


@pytest.mark.l1
@pytest.mark.rests_on(f"{_THIS}::test_return_fraction_is_the_declared_and_the_measured")
@pytest.mark.parametrize("s,d", _POINTS)
def test_the_realizer_is_a_homomorphism(s: float, d: float) -> None:
    r"""Law L-L at two levels (decision 52).

    (a) EXACT on the parameters: ``real(s S + d D)`` is the ``RealizedSum`` of ``ScaledOrdinatePermutation(alpha=s)``
    and ``RankOneReemission(alpha=d)`` by parameter equality; ``real(c x)`` has the parameters of ``c`` times those of
    ``real(x)``; the image of the empty network is the ``ZeroMap``; ``real(x).return_fraction`` equals the wall's
    by ``math.fsum`` on both sides.

    (b) On dense matrices of the S16 slab face: ``array_equal`` for the specular image at every coefficient (its
    entries are 1.0) and for every image at dyadic coefficients; fewer than 5 ulp per entry otherwise (``[M]`` over
    N = 20000 draws at most 2 ulp for the scalar leg, below 5 for the additive leg).

    The mutation M-rescale (the diffuse image rescaled by the wall's total ``rho``) is inert on ``(0.3, 0.7)``, which
    is why that point is the carry witness and not the law's only fixture, and reddens on ``(0.3, 0.5)``.

    Red input: the base tip has no table.
    """
    B = H.bc()
    R = H.ro()
    face = "xmax"
    quad = H.QUADRATURES["gl16"]()
    trace = H.trace_of(quad)
    n = H.n_outflow(trace, face)
    wall = B.WallCondition.from_amplitudes(s, d)
    op, _, _ = _realize(wall, "gl16", face, trace)
    mirror = H.mirror_pairing(quad, face)
    spec = R.ScaledOrdinatePermutation(alpha=s, permutation=mirror, face=face, partner=face, trace=trace)
    diff = R.RankOneReemission(alpha=d, face=face, trace=trace)
    assert op == R.RealizedSum(terms=(spec, diff)), "the image of a sum is the sum of the images, in the canonical chain order"
    assert op.return_fraction == wall.return_fraction == math.fsum([s, d])

    half = 0.5 * B.WallCondition.specular(s)
    scaled, _, _ = _realize(half, "gl16", face, trace)
    assert scaled == R.ScaledOrdinatePermutation(alpha=0.5 * s, permutation=mirror, face=face, partner=face, trace=trace)
    zero, _, _ = _realize(B.WallCondition.vacuum(), "gl16", face, trace)
    assert isinstance(zero, R.ZeroMap) and zero.return_fraction == 0.0

    unit_spec = H.dense(_realize(B.WallCondition.specular(1.0), "gl16", face, trace)[0], n)
    unit_diff = H.dense(_realize(B.WallCondition.white(1.0), "gl16", face, trace)[0], n)
    got = H.dense(op, n)
    want = s * unit_spec + d * unit_diff
    dyadic = all(float(x) * 8 == int(x * 8) for x in (s, d))
    if dyadic:
        np.testing.assert_array_equal(got, want)
    else:
        assert _ulps(got, want) < 5, f"(s, d) = ({s}, {d}): {_ulps(got, want)} ulp"
    if d == 0.0 or s == 0.0:
        np.testing.assert_array_equal(got, want)
    np.testing.assert_array_equal(H.dense(_realize(B.WallCondition.specular(s), "gl16", face, trace)[0], n), s * unit_spec)
