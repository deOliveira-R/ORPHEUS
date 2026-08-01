r"""Freeze an INDEPENDENT reference image for each BC × quadrature case.

This module computes, for seven boundary-law cases, the incoming half-trace
:math:`\gamma_-\psi` that the *mathematics* of the law demands — from the
affine form :eq:`affine-bc-form`, the mirror isometry, and the Lambertian
partial-current balance — and freezes it under
``tests/geometry/snapshots/bc_equivalence_<case>.npz``. The companion harness
(``tests/geometry/test_bc_equivalence_snapshot.py``) re-runs each case's
realized SN operator on every pytest invocation and asserts it reproduces the
frozen reference at the per-case tolerance.

Run::

    python -m tests.geometry._generate_bc_equivalence_snapshots
    python -m tests.geometry._generate_bc_equivalence_snapshots --case vacuum_lebedev17

.. warning::

   **This module MUST NOT import** :class:`~orpheus.sn.boundary.realizer.SNBoundaryRealizer`
   **or anything else on the SN realization path, and the harness MUST NOT
   import the ``_*_image`` expression helpers below.** Both halves of that
   sentence are load-bearing; see :ref:`the inversion <bc-snapshot-inversion>`.
   (Importing the Layer-2 *law* descriptors from
   :mod:`orpheus.geometry.boundary` is fine and deliberate — they are
   declarative statements of which BC each case is about, held in the registry
   for the harness. The reference is derived from the law's EQUATION, never
   from the law object, and no ``_*_image`` function touches one.)

.. _bc-snapshot-inversion:

The 2026-08-01 inversion — the generator is no longer a recorder
================================================================

Until this rewrite ``_build_payload`` called ``SNBoundaryRealizer().realize(…)``
and froze ``op.apply(psi_out)``. The artefact was therefore **a recording of
production output**, and two consequences followed mechanically:

1. **The assertion was self-referential.** ``production == a recording of
   production`` is a regression LOCK: it detects *change*, and is worth exactly
   what the recorded code was right. It cannot say the value is correct
   (``vv-principles`` §bit-identity criterion 2).
2. **It broke on every signature change.** It did, twice — at B3.2 (the law
   narrowed to :math:`\Gamma_+ \to \Gamma_-`) and at B3.4a (white and
   prescribed inflow followed) — because a recorder is coupled to the very
   thing it records.

The frozen artefacts had one genuine virtue: they predated the narrowing, so
old-vs-new was a statement about two different code states rather than about
one. But that virtue is only **procedural** independence
(``vv-principles`` L11 / §structural independence): it certified *"the new
path agrees with the old path"*, and this campaign's premise is precisely that
the old path read the wrong half-trace.

The user's ruling (2026-08-01) was therefore to **migrate the schema and
re-anchor every case against an expression written from the mathematics**. All
seven artefacts were regenerated; the reference is now derived, not recorded:

* ``psi_in`` is computed by the ``_*_image`` functions below from the law's
  equation. Nothing in this module knows how the SN realizer is implemented.
* The generator cannot be broken by a realizer change again — it has no
  realizer coupling to break.
* The artefact **IS** the reference rather than a recording, so the harness's
  assertion is a genuine statement about correctness (an L1-shaped claim on a
  software invariant) instead of a drift lock.

What the frozen file still buys, now that the reference is computable
--------------------------------------------------------------------

It is the guard against **generator and production drifting together**. The
harness reads ``psi_in`` off disk and never calls the expressions below; so
changing an expression here does *not* silently move the gate — it moves it
only through a regeneration, which is a reviewable diff of committed binary
artefacts. That is why the two modules must not share the expression code
(and why the harness importing ``_specular_image`` would quietly undo the whole
inversion).

Schema
======

Uniform across every case — each artefact is self-describing:

===========================  ==================================================
key                          content
===========================  ==================================================
``psi_out``                   :math:`\gamma_+\psi`, shape ``(|Γ₊|, 5, 3)`` — the
                              law's INPUT under the narrowed
                              :math:`\Gamma_+ \to \Gamma_-` contract.
``psi_in``                    the independently-derived reference image, shape
                              ``(|Γ₋|, 5, 3)``.
``face``                      the face the law is installed on.
``inflow_indices``            :math:`\Gamma_-` — sorted global ordinate indices.
``outflow_indices``           :math:`\Gamma_+` — sorted global ordinate indices.
``case_id`` / ``description`` traceability.
===========================  ==================================================

``periodic_lebedev17`` additionally carries ``partner_face``,
``partner_outflow_indices`` and ``psi_out_partner``: a periodic law is a
two-face coupling, and its reference is a function of the PARTNER face's
outflow. Those keys are consumed, not decoration (this module's schema rule is
*store exactly what is persisted AND compared*).

``psi_out`` is drawn full-face with the same deterministic RNG the pre-B3.2
artefacts used and then RESTRICTED to :math:`\Gamma_+`. That is free
traceability: the stored rows are literally a subset of the numbers the retired
artefacts probed with, so the migration changed the schema and the reference —
never the probe.

Per-case tolerance
==================

Derived, never measured-and-rounded-up. The harness states each argument at its
own gate; in summary, with :math:`u = \varepsilon_{\text{mach}}`:

* **vacuum** — ``assert_array_equal`` against ``0``. The zero map performs no
  arithmetic.
* **specular** (both) — ``assert_array_equal``. The image is a gather, times a
  scalar; reduction depth 0, so the predicted drift is EXACTLY zero and a
  tolerance would admit the failure it is meant to catch.
* **white** (both) — ``rtol = |Γ₊| · u``. The reference contracts with
  :func:`numpy.tensordot` where production runs a broadcast-multiply-then-
  ``sum``; same value, different reduction ORDER. Every summand is positive
  (``ψ ∼ U(0,2)``, ``w > 0``, ``|Ω·n̂| > 0``) so the sum's condition number is
  1 and the classical ``(n−1)·u/2`` bound applies without amplification —
  ``|Γ₊| · u`` is that bound with room for the division.
* **mixed** — ``rtol = (|Γ₊| + 2) · u``: the white term's bound plus one
  rounding for its scaling and one for the sum. No cancellation (both weights
  and both terms positive), so again κ = 1.
* **periodic** — ``assert_array_equal`` (a gather), inside a strict xfail; see
  the harness.

The snapshot directory is committed — the snapshots ARE the verification
artefact. Regeneration is intentional only when the reference EXPRESSION or the
probe changes, and must carry an inline justification.
"""
from __future__ import annotations

import argparse
import hashlib
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Callable

import numpy as np

from orpheus.geometry.boundary import (
    PeriodicBoundary,
    ReflectiveBoundary,
    VacuumInflow,
    WhiteBoundary,
)
from orpheus.numerics.quadrature import Quadrature

if TYPE_CHECKING:
    from orpheus.geometry.boundary import BoundaryTraceLaw
    from orpheus.numerics.operator import LinearOperator

    #: What the harness hands :attr:`BCEquivalenceCase.compose`: its own
    #: realization step, ``law -> operator``. Typing the case's composition as
    #: a function OF the realizer is what lets a rank-N case state its own
    #: operator algebra (``0.3 * a + 0.7 * b``) in the registry while this
    #: module stays free of any realization import.
    Realize = Callable[["BoundaryTraceLaw"], "LinearOperator"]


SNAPSHOT_DIR = Path(__file__).parent / "snapshots"

# Input shape used for every snapshot: ``(N_ordinates, n_cells, n_groups)``.
# The trailing ``(5, 3)`` is arbitrary but NOT square — a group ↔ spatial axis
# transpose cannot hide behind it — and it is unchanged from the pre-B3.2
# artefacts, which is what makes the restricted probe traceable to them.
_PSI_TRAILING_SHAPE: tuple[int, int] = (5, 3)

#: Tangential band for the :math:`\Omega\cdot\hat n` classification, stated
#: here from its derivation rather than imported: a ×4 safety factor over the
#: IEEE-754 round-off bound of a 3-component unit-vector dot product
#: (:math:`d\,u` with :math:`d \le 3`, :math:`u = \varepsilon/2`). It equals
#: :data:`~orpheus.numerics.spaces.angular_trace_space.TANGENTIAL_EPS` today,
#: and the harness cross-checks the resulting index SETS against the live trace
#: space — so a divergence reds a gate instead of silently re-anchoring one.
#: (The B3.4a ``> 0.0`` classification twin is exactly this kind of drift.)
_TANGENTIAL_EPS: float = 4.0 * float(np.finfo(np.float64).eps)


# ─── seed derivation ─────────────────────────────────────────────────


def _seed_from_case_id(case_id: str) -> int:
    """Deterministic small seed from a case ID.

    Using ``hashlib.sha256`` (not ``hash()``) so the seed is stable
    across Python interpreter invocations (PYTHONHASHSEED-immune).
    Truncated to 16 bits — small enough for ``np.random.default_rng``
    to print cleanly in error messages, large enough that distinct case
    IDs collide rarely.
    """
    digest = hashlib.sha256(case_id.encode("utf-8")).digest()
    # Two-byte truncation; matches the "small integer derived from the
    # case ID" requirement in the plan.
    return int.from_bytes(digest[:2], byteorder="big")


def _generate_psi(quad: Quadrature, case_id: str) -> np.ndarray:
    """Deterministic FULL-FACE input draw for one snapshot case.

    ``np.random.default_rng(<seed>).uniform(0, 2, size=(N, 5, 3))`` — unchanged
    from the pre-B3.2 generator, so the :math:`\\Gamma_+` restriction stored in
    the migrated artefact is a subset of the numbers the retired artefacts
    probed with. The ``(0, 2)`` range is symmetric around 1.0, keeps inputs
    O(1), and — load-bearing for the white/mixed tolerance argument — is
    **non-negative**, which makes the cosine-weighted sum perfectly conditioned.
    """
    seed = _seed_from_case_id(case_id)
    rng = np.random.default_rng(seed)
    shape = (int(quad.N),) + _PSI_TRAILING_SHAPE
    return rng.uniform(0.0, 2.0, size=shape)


# ─── Layer-1 geometry: the half-traces, from the face convention ─────


@dataclass(frozen=True)
class HalfTraces:
    r"""The two directional half-traces of ONE boundary face.

    Everything here is derived from the documented geometric convention
    (:mod:`orpheus.numerics.spaces.angular_trace_space`): a face named
    ``<axis><min|max>`` has outward unit normal :math:`\mp\hat e_{\rm axis}`,
    the signed projection is :math:`\Omega\cdot\hat n = \mathrm{sign}\cdot
    \mu_{\rm axis}`, and an ordinate is INFLOW when that projection is
    :math:`< -\epsilon`, OUTFLOW when :math:`> +\epsilon`, and tangential
    (in NEITHER half) in between.
    """

    quadrature: Quadrature
    face: str
    axis: str
    outward_sign: int
    inflow: np.ndarray
    outflow: np.ndarray

    @property
    def omega_dot_n(self) -> np.ndarray:
        r""":math:`\Omega_n\cdot\hat n_f` for every ordinate."""
        mu = np.asarray(getattr(self.quadrature, f"mu_{self.axis}"), dtype=float)
        return float(self.outward_sign) * mu

    @property
    def cosine_weight(self) -> np.ndarray:
        r"""The boundary measure :math:`w_n\,|\Omega_n\cdot\hat n_f|`.

        This is the discrete partial-current weight — the same measure the
        trace space installs as its inner product — and the only weight the
        Lambertian law's normalisation may use.
        """
        return np.asarray(self.quadrature.weights, dtype=float) * np.abs(
            self.omega_dot_n,
        )


def _half_traces(quad: Quadrature, face: str) -> HalfTraces:
    """Classify a quadrature's ordinates against one face's outward normal."""
    axis, suffix = face[0], face[1:]
    if suffix not in ("min", "max") or axis not in ("x", "y", "z"):
        raise ValueError(f"unparseable face name {face!r}")
    sign = +1 if suffix == "max" else -1
    mu = np.asarray(getattr(quad, f"mu_{axis}"), dtype=float)
    omega_dot_n = float(sign) * mu
    return HalfTraces(
        quadrature=quad,
        face=face,
        axis=axis,
        outward_sign=sign,
        inflow=np.flatnonzero(omega_dot_n < -_TANGENTIAL_EPS).astype(np.intp),
        outflow=np.flatnonzero(omega_dot_n > +_TANGENTIAL_EPS).astype(np.intp),
    )


def _mirror_partner(quad: Quadrature, axis: str) -> np.ndarray:
    r"""The specular pairing, from the isometry — NOT from the quadrature's table.

    The mirror fixing a face with normal :math:`\hat n` is the unique ambient
    isometry :math:`\Omega \mapsto \Omega - 2(\Omega\cdot\hat n)\hat n`, i.e.
    "negate the axis cosine, keep the other two". This function realizes that
    map on the ordinate set by matching the transformed direction cosines back
    against the quadrature's own, requiring an EXACT and UNIQUE hit.

    Exactness is not an optimism: MEASURED 2026-08-01 on ``gauss_legendre(4)``
    / ``(8)``, ``lebedev(17)``, ``level_symmetric(4)`` / ``(6)`` and
    ``product(2,4)``, every reflected direction is bit-for-bit present in the
    ordinate table. Demanding it (rather than falling back to a nearest
    neighbour) is what keeps this an independent derivation: a tolerant match
    would quietly agree with a WRONG table, which is the failure mode the
    harness's specular rows exist to catch.
    """
    if axis not in ("x", "y", "z"):
        raise ValueError(f"unknown mirror axis {axis!r}")
    n = int(quad.N)
    cosines = np.stack(
        [
            np.asarray(getattr(quad, f"mu_{a}", np.zeros(n)), dtype=float)
            for a in ("x", "y", "z")
        ],
        axis=1,
    )
    reflected = cosines.copy()
    reflected[:, "xyz".index(axis)] *= -1.0
    partner = np.empty(n, dtype=np.intp)
    for k in range(n):
        hit = np.flatnonzero(np.all(cosines == reflected[k], axis=1))
        if hit.size != 1:
            raise AssertionError(
                f"the mirror about {axis!r} sends ordinate {k} "
                f"{tuple(cosines[k])} to {tuple(reflected[k])}, which matches "
                f"{hit.size} ordinates exactly (expected 1). This quadrature "
                f"is not closed under the reflection, so the specular pairing "
                f"cannot be derived from the isometry — and a tolerant match "
                f"would fabricate one."
            )
        partner[k] = hit[0]
    if not np.array_equal(partner[partner], np.arange(n)):
        raise AssertionError(
            f"the derived pairing about {axis!r} is not an involution; a "
            f"mirror is its own inverse, so the derivation is wrong."
        )
    return partner


# ─── Layer-2 laws: the independent reference expressions ─────────────
#
# Each function below writes ONE boundary law's image from its equation. None
# of them imports, reads, or mirrors ``orpheus.sn.boundary`` — that is the
# whole point of the inversion, and a transcription of the realizer would be
# procedural independence wearing a new coat.


@dataclass(frozen=True)
class ReferenceInputs:
    r"""Everything an independent expression is allowed to read."""

    law_face: HalfTraces
    psi_out: np.ndarray
    partner_face: HalfTraces | None = None
    psi_out_partner: np.ndarray | None = None


def _zero_image(ri: ReferenceInputs) -> np.ndarray:
    r"""**Vacuum.** :math:`R = 0` and :math:`q = 0` in :eq:`affine-bc-form`, so
    :math:`\gamma_-\psi = 0`.

    The content is the CODOMAIN as much as the value: a vacuum boundary emits
    the zero of :math:`\Gamma_-`, which has :math:`|\Gamma_-|` rows — not an
    echo of the :math:`|\Gamma_+|` rows it consumed. The two happen to be
    equal-sized on every quadrature in the tree, which is why the shape must be
    written from :math:`\Gamma_-` deliberately rather than from the input.
    """
    return np.zeros(
        (ri.law_face.inflow.size,) + ri.psi_out.shape[1:], dtype=float,
    )


def _specular_image(ri: ReferenceInputs, *, albedo: float) -> np.ndarray:
    r"""**Specular reflection with albedo** :math:`\alpha`.

    A mirror boundary returns the flux arriving from the mirrored direction:

    .. math::

        \psi^-(\Omega) \;=\; \alpha\,\psi^+(\Omega'), \qquad
        \Omega' \;=\; \Omega - 2(\Omega\cdot\hat n)\hat n .

    So row :math:`j` of the image — the inflow ordinate ``inflow[j]`` — reads
    the outflow ordinate ``partner[inflow[j]]``, scaled by :math:`\alpha`.
    ``psi_out`` carries only the :math:`\Gamma_+` rows, so the read is
    expressed by scattering it back to a full-face buffer and gathering the
    mirrored indices; that keeps the expression in GLOBAL ordinate numbering,
    where the mirror is defined, instead of re-deriving a local permutation
    (which is the realizer's own problem and must not be this file's).

    ERR-045 is asserted rather than assumed: the mirror of an inflow ordinate
    must be an OUTFLOW ordinate. If it is not, the law is not a boundary law at
    this face and no reference exists.
    """
    ht = ri.law_face
    partner = _mirror_partner(ht.quadrature, ht.axis)
    sources = partner[ht.inflow]
    if not np.all(np.isin(sources, ht.outflow)):
        raise AssertionError(
            f"ERR-045 violated at face {ht.face!r} about axis {ht.axis!r}: the "
            f"mirror sends inflow ordinate(s) "
            f"{sorted(set(sources) - set(ht.outflow.tolist()))} outside Γ₊. A "
            f"reflection must exchange the half-traces."
        )
    full = np.zeros((int(ht.quadrature.N),) + ri.psi_out.shape[1:], dtype=float)
    full[ht.outflow] = ri.psi_out
    return albedo * full[sources]


def _lambertian_image(ri: ReferenceInputs, *, albedo: float) -> np.ndarray:
    r"""**White (Lambertian) re-emission with albedo** :math:`\alpha`.

    Two statements define it, and neither mentions an implementation:

    1. **Isotropy.** The re-emitted flux is uniform over :math:`\Gamma_-`:
       :math:`\psi^-(\Omega) \equiv C` for every incoming ordinate.
    2. **Partial-current balance.** The returned current is :math:`\alpha`
       times the incident one, in the boundary measure
       :math:`|\Omega\cdot\hat n|\,\mathrm{d}\Omega`:

    .. math::

        J^- \;=\; \alpha\,J^+, \qquad
        J^\pm \;=\; \textstyle\sum_{n\in\Gamma_\pm}
        w_n\,|\Omega_n\cdot\hat n|\;\psi_n .

    Substituting (1) into (2) gives the whole law:

    .. math::

        C \;=\; \alpha\,
        \frac{\sum_{n\in\Gamma_+} w_n |\Omega_n\cdot\hat n|\,\psi^+_n}
             {\sum_{n\in\Gamma_-} w_n |\Omega_n\cdot\hat n|} .

    Note the normalisation runs over :math:`\Gamma_-` — it comes from the
    *returned* current, which is where the physics puts it. Production
    normalises over :math:`\Gamma_+` (equivalently: it re-emits *the
    cosine-weighted average of the outflow*). The two spellings agree because
    the mirror is an exact measure-preserving bijection
    :math:`\Gamma_+ \leftrightarrow \Gamma_-`, so the two sums run over the
    same multiset — :func:`_assert_the_normalisation_choice_is_immaterial`
    checks that at generation time rather than trusting it, because a
    quadrature for which it failed would make this reference the correct one
    and production's the wrong one.
    """
    ht = ri.law_face
    cos_w = ht.cosine_weight
    incident_current = np.tensordot(cos_w[ht.outflow], ri.psi_out, axes=(0, 0))
    _assert_the_normalisation_choice_is_immaterial(ht)
    emitted = albedo * incident_current / cos_w[ht.inflow].sum()
    return np.broadcast_to(
        emitted, (ht.inflow.size,) + incident_current.shape,
    ).copy()


def _assert_the_normalisation_choice_is_immaterial(ht: HalfTraces) -> None:
    r"""Both hemispheres carry the same total boundary measure.

    An explicit ``raise``, NOT a bare ``assert``: this module is a generator,
    never COLLECTED by pytest, so its asserts are not AST-rewritten and
    ``python -O`` strips them to a no-op (``vv-principles`` Mode 8 — the
    hazard's genuine domain is exactly non-collected code like this).
    """
    cos_w = ht.cosine_weight
    plus, minus = cos_w[ht.outflow].sum(), cos_w[ht.inflow].sum()
    bound = max(ht.outflow.size, ht.inflow.size) * float(np.finfo(np.float64).eps)
    if abs(plus - minus) > bound * max(abs(plus), abs(minus)):
        raise AssertionError(
            f"face {ht.face!r}: Σ_Γ₊ w|Ω·n̂| = {plus!r} but Σ_Γ₋ w|Ω·n̂| = "
            f"{minus!r} (relative gap exceeds the {bound:.2e} reduction-order "
            f"bound). The quadrature is not symmetric about this face, so "
            f"'normalise by the emitted current' and 'normalise by the "
            f"incident current' are DIFFERENT laws — the first is the "
            f"physical one, and production implements the second."
        )


def _mixed_image(ri: ReferenceInputs) -> np.ndarray:
    r"""**Marshak-style mix** :math:`0.3\,\text{specular} + 0.7\,\text{white}`.

    The rank-N law is a convex combination of response kernels acting on the
    SAME outflow trace, so the image is the pointwise combination of the two
    images — computed here from the two expressions above, never from realized
    leaves. (The pre-migration generator composed
    ``0.3 * spec_realized + 0.7 * white_realized`` out of realizer output,
    which made this the most self-referential row in the file.)

    Both leaves are at :math:`\alpha = 1`; the 0.3/0.7 weights are the mixing
    coefficients of the law, not albedos.
    """
    return 0.3 * _specular_image(ri, albedo=1.0) + 0.7 * _lambertian_image(
        ri, albedo=1.0,
    )


def _periodic_image(ri: ReferenceInputs) -> np.ndarray:
    r"""**Periodic wrap.** The partner face's outflow IS this face's inflow.

    Periodicity identifies the two faces by a translation :math:`T`, so the
    deck transformation carries the trace across the domain WITHOUT touching
    the direction:

    .. math::

        \psi^-\!\left(\mathbf{r}_f, \Omega\right) \;=\;
        \psi^+\!\left(T\mathbf{r}_f, \Omega\right).

    A particle leaving through ``xmax`` in direction :math:`\Omega` re-enters
    at ``xmin`` in the SAME direction :math:`\Omega` — which is why
    :math:`\Gamma_-(\text{xmin})` and :math:`\Gamma_+(\text{xmax})` are the
    same ordinate SET (asserted below, not assumed), and why the map is the
    identity in the angular index while being a genuine motion in space.

    The load-bearing consequence: the image is a function of the PARTNER
    face's outflow and of nothing on this face. A law realized as a per-face
    endomorphism cannot express it, however it permutes.

    That is exactly what production did until campaign phase **B3.4c**, and it
    is why this case shipped as a strict xfail: the composition fed every law
    its own face's :math:`\Gamma_+`, so the row reddened at 98 % relative. Since
    B3.4c the geometry factor names the partner
    (:meth:`SpatialWrap.domain_face`) and the composition supplies it, the
    marker is deleted, and the harness row is LIVE — it asks the law which face
    to read and feeds the probe that answer names, so a regression selects the
    wrong probe and reds again.
    """
    if ri.partner_face is None or ri.psi_out_partner is None:
        raise AssertionError(
            "the periodic reference needs the partner face's outflow; a "
            "single-face probe cannot state the law."
        )
    if not np.array_equal(ri.law_face.inflow, ri.partner_face.outflow):
        raise AssertionError(
            f"Γ₋({ri.law_face.face}) = {ri.law_face.inflow.tolist()} but "
            f"Γ₊({ri.partner_face.face}) = "
            f"{ri.partner_face.outflow.tolist()}. A translation does not "
            f"change the direction, so a periodic pair must agree on the "
            f"ordinate set — these faces are not a periodic pair."
        )
    return np.asarray(ri.psi_out_partner, dtype=float).copy()


# ─── case registry ───────────────────────────────────────────────────


@dataclass(frozen=True)
class BCEquivalenceCase:
    """One snapshot case: a law, where it is installed, and its reference.

    The registry is the SINGLE source of truth for a case's identity: the
    harness imports these entries and builds its system-under-test from
    ``compose`` / ``build_quadrature`` / ``face`` / ``faces`` rather than
    re-spelling them, so the two sides cannot drift into testing different
    problems (``coding-elegance`` Pattern 2 — before this migration the law
    and quadrature of every case were written out twice, once here and once in
    the harness).

    ``compose`` is declared here and **never called by this module**: it is a
    function of the harness's own realization step, so a rank-1 case is
    ``lambda realize: realize(VacuumInflow())`` and the rank-N case states its
    operator algebra directly. The reference is derived from the law's
    EQUATION, not from any law object.
    """

    case_id: str
    description: str
    build_quadrature: Callable[[], Quadrature]
    #: ``realize -> operator``. See the class docstring.
    compose: "Callable[[Realize], LinearOperator]"
    face: str
    faces: tuple[str, ...]
    reference: Callable[[ReferenceInputs], np.ndarray]
    #: Periodic only: the face whose outflow feeds this face's inflow.
    partner_face: str | None = None

    @property
    def snapshot_path(self) -> Path:
        return SNAPSHOT_DIR / f"bc_equivalence_{self.case_id}.npz"


_FOUR_FACES = ("xmin", "xmax", "ymin", "ymax")


CASES: tuple[BCEquivalenceCase, ...] = (
    BCEquivalenceCase(
        case_id="vacuum_lebedev17",
        description=(
            "VacuumInflow + LebedevSphere(17) at xmin. R = q = 0, so the "
            "reference is the zero of Γ₋ — sized from the CODOMAIN, not "
            "echoed from the input."
        ),
        build_quadrature=lambda: Quadrature.lebedev(17),
        compose=lambda realize: realize(VacuumInflow()),
        face="xmin",
        faces=("xmin", "xmax"),
        reference=_zero_image,
    ),
    BCEquivalenceCase(
        case_id="specular_x_lebedev17",
        description=(
            "ReflectiveBoundary(axis='x', albedo=1.0) + LebedevSphere(17) at "
            "xmax. Reference: the mirror gather Ω ↦ Ω − 2(Ω·n̂)n̂, α = 1."
        ),
        build_quadrature=lambda: Quadrature.lebedev(17),
        compose=lambda realize: realize(
            ReflectiveBoundary(axis="x", albedo=1.0),
        ),
        face="xmax",
        faces=_FOUR_FACES,
        reference=lambda ri: _specular_image(ri, albedo=1.0),
    ),
    BCEquivalenceCase(
        case_id="specular_y_partial_07_LS6",
        description=(
            "ReflectiveBoundary(axis='y', albedo=0.7) + LevelSymmetricSN(6) at "
            "ymax. Reference: the mirror gather about y, scaled by 0.7 — the "
            "α-fold row, on a non-x axis and a non-max-x face."
        ),
        build_quadrature=lambda: Quadrature.level_symmetric(sn_order=6),
        compose=lambda realize: realize(
            ReflectiveBoundary(axis="y", albedo=0.7),
        ),
        face="ymax",
        faces=_FOUR_FACES,
        reference=lambda ri: _specular_image(ri, albedo=0.7),
    ),
    BCEquivalenceCase(
        case_id="white_xmax_LS4",
        description=(
            "WhiteBoundary(axis='x', outward_sign=+1, albedo=1.0) + "
            "LevelSymmetricSN(4) at xmax. Reference: isotropic re-emission "
            "with J⁻ = J⁺ in the |Ω·n̂| measure."
        ),
        build_quadrature=lambda: Quadrature.level_symmetric(sn_order=4),
        compose=lambda realize: realize(
            WhiteBoundary(axis="x", outward_sign=+1, albedo=1.0),
        ),
        face="xmax",
        faces=_FOUR_FACES,
        reference=lambda ri: _lambertian_image(ri, albedo=1.0),
    ),
    BCEquivalenceCase(
        case_id="white_xmin_partial_03_GL",
        description=(
            "WhiteBoundary(axis='x', outward_sign=-1, albedo=0.3) + "
            "GaussLegendre1D(8) at xmin. Same law on a 1-D quadrature whose "
            "Σw is 2, not 4π — the canary against a hard-coded normalisation "
            "(numerical-bug-signatures Signature 4)."
        ),
        build_quadrature=lambda: Quadrature.gauss_legendre(n_ordinates=8),
        compose=lambda realize: realize(
            WhiteBoundary(axis="x", outward_sign=-1, albedo=0.3),
        ),
        face="xmin",
        faces=("xmin", "xmax"),
        reference=lambda ri: _lambertian_image(ri, albedo=0.3),
    ),
    BCEquivalenceCase(
        case_id="mixed_30spec_70white_LS4",
        description=(
            "0.3 · ReflectiveBoundary(axis='x', albedo=1.0) + 0.7 · "
            "WhiteBoundary(axis='x', outward_sign=+1, albedo=1.0) + "
            "LevelSymmetricSN(4) at xmax. Reference: the pointwise convex "
            "combination of the two independent images."
        ),
        build_quadrature=lambda: Quadrature.level_symmetric(sn_order=4),
        # Rank-N: the composition IS the Wave-0 operator algebra over realized
        # leaves (``MixedBoundaryOperator`` was removed in Wave 11), so the
        # case states it directly rather than describing it.
        compose=lambda realize: (
            0.3 * realize(ReflectiveBoundary(axis="x", albedo=1.0))
            + 0.7 * realize(
                WhiteBoundary(axis="x", outward_sign=+1, albedo=1.0),
            )
        ),
        face="xmax",
        faces=_FOUR_FACES,
        reference=_mixed_image,
    ),
    BCEquivalenceCase(
        case_id="periodic_lebedev17",
        description=(
            "PeriodicBoundary + LebedevSphere(17), xmin fed by xmax. "
            "Reference: the PARTNER face's outflow at the same ordinates — a "
            "two-face coupling no per-face endomorphism can express (B3.4c)."
        ),
        build_quadrature=lambda: Quadrature.lebedev(17),
        compose=lambda realize: realize(PeriodicBoundary()),
        face="xmin",
        faces=("xmin", "xmax"),
        reference=_periodic_image,
        partner_face="xmax",
    ),
)


def case_by_id(case_id: str) -> BCEquivalenceCase:
    """The registry entry for ``case_id`` — the harness's lookup."""
    for case in CASES:
        if case.case_id == case_id:
            return case
    raise KeyError(
        f"unknown BC-equivalence case {case_id!r}; known: "
        f"{[c.case_id for c in CASES]}"
    )


# ─── snapshot writer ─────────────────────────────────────────────────


def _build_payload(case: BCEquivalenceCase) -> dict:
    """Compute one case's payload: the probe, the derived reference, the metadata.

    Nothing here realizes a boundary law. ``psi_out`` is the deterministic
    full-face draw restricted to :math:`\\Gamma_+`; ``psi_in`` is
    ``case.reference(...)``, written from the law's equation.
    """
    quad = case.build_quadrature()
    law_face = _half_traces(quad, case.face)
    full = _generate_psi(quad, case.case_id)
    psi_out = full[law_face.outflow]

    partner_face = psi_out_partner = None
    if case.partner_face is not None:
        partner_face = _half_traces(quad, case.partner_face)
        # A SECOND, independent draw. Reusing ``full`` would put identical data
        # on both faces, and a periodic law's whole content is that it reads
        # the OTHER face — a shared draw would let a per-face endomorphism look
        # correct on the rows that happen to coincide.
        partner_full = _generate_psi(quad, f"{case.case_id}:partner")
        psi_out_partner = partner_full[partner_face.outflow]

    psi_in = case.reference(
        ReferenceInputs(
            law_face=law_face,
            psi_out=psi_out,
            partner_face=partner_face,
            psi_out_partner=psi_out_partner,
        )
    )

    payload: dict = dict(
        psi_out=np.asarray(psi_out, dtype=np.float64),
        psi_in=np.asarray(psi_in, dtype=np.float64),
        face=np.array(case.face),
        inflow_indices=law_face.inflow.astype(np.int64),
        outflow_indices=law_face.outflow.astype(np.int64),
        case_id=np.array(case.case_id),
        description=np.array(case.description),
    )
    if partner_face is not None and psi_out_partner is not None:
        payload["partner_face"] = np.array(partner_face.face)
        payload["partner_outflow_indices"] = partner_face.outflow.astype(np.int64)
        payload["psi_out_partner"] = np.asarray(psi_out_partner, dtype=np.float64)
    return payload


def _array_equal_after_load(snapshot_path: Path, payload: dict) -> bool:
    """Check whether the on-disk snapshot already matches the payload."""
    if not snapshot_path.exists():
        return False
    try:
        existing = np.load(snapshot_path)
    except Exception:
        return False
    if set(existing.files) != set(payload):
        return False
    for key, expected in payload.items():
        if not np.array_equal(existing[key], expected):
            return False
    return True


def generate_one(case: BCEquivalenceCase) -> tuple[Path, str]:
    """Write the snapshot for one case; return (path, status).

    Status is ``"UNCHANGED"`` when an existing snapshot already matches
    the freshly-built payload (bit-for-bit on every array). Status is
    ``"REGENERATED"`` when the snapshot is new or differs from disk.
    """
    payload = _build_payload(case)
    SNAPSHOT_DIR.mkdir(parents=True, exist_ok=True)
    out = case.snapshot_path
    if _array_equal_after_load(out, payload):
        status = "UNCHANGED"
    else:
        np.savez_compressed(out, **payload)
        status = "REGENERATED"
    return out, status


def generate_all(case_ids: list[str] | None = None) -> list[tuple[Path, str]]:
    """Generate every (or the named subset of) snapshots."""
    written: list[tuple[Path, str]] = []
    for case in CASES:
        if case_ids is not None and case.case_id not in case_ids:
            continue
        path, status = generate_one(case)
        written.append((path, status))
    return written


# ─── CLI ─────────────────────────────────────────────────────────────


def _print_summary(written: list[tuple[Path, str]]) -> None:
    """Pretty-print the per-snapshot status table."""
    if not written:
        print("(no snapshots generated — case-id filter matched nothing)")
        return
    name_w = max(len(p.name) for p, _ in written)
    print(f"{'Snapshot':<{name_w}}  Status        Size")
    print(f"{'-' * name_w}  ------------  ---------")
    for path, status in written:
        size_kb = path.stat().st_size / 1024.0 if path.exists() else 0.0
        print(f"{path.name:<{name_w}}  {status:<12}  {size_kb:6.1f} KiB")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Freeze the independently-derived BC reference images "
            "(re-anchored 2026-08-01)."
        ),
    )
    parser.add_argument(
        "--case", action="append", default=None,
        help="Restrict to a specific case_id (repeatable).",
    )
    parser.add_argument(
        "--list", action="store_true",
        help="List available cases and exit.",
    )
    args = parser.parse_args(argv)
    if args.list:
        name_w = max(len(c.case_id) for c in CASES)
        for case in CASES:
            print(f"  {case.case_id:<{name_w}}  {case.description}")
        return 0
    written = generate_all(case_ids=args.case)
    _print_summary(written)
    return 0


if __name__ == "__main__":
    sys.exit(main())
