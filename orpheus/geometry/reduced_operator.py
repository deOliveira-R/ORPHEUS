r"""The per-(cell, direction) sweep packet — the last tenant of a dissolving module.

:class:`StreamingTerms` is the frozen packet a curvilinear S\ :sub:`N` cell
balance is handed: the chord length, the two face areas, the radial cosine and
its magnitude, the cell volume, and :math:`\Delta A / w` — the connection
factor already divided by the ordinate's quadrature weight.

⚠ **This module holds nothing else, and it is not staying.**  It is the residue
of the un-weld arc: P4.4 (2026-08-28) moved the connection-coefficient operator
and its three chart factories to :mod:`orpheus.sn.mesh.reduced_operator`, and
P4.2 moved the :math:`\alpha` dome, its admission contract, the angular factor
and the measure adapter to :mod:`orpheus.sn.angular.redistribution`.  What is
left is one frozen dataclass of floats.

⟹ **P4.3 moves it to** ``orpheus/transport/spatial/`` **and deletes this file.**
That is where both of its runtime importers already live
(:mod:`~orpheus.transport.spatial.cell_balance`,
:mod:`~orpheus.transport.spatial.scheme`), and the packet is an L2 object: a
spatial scheme consumes it, and `[M]` this module now imports **nothing** —
not ``numpy``, not :class:`~orpheus.geometry.coord.CoordSystem`, not
:class:`~orpheus.geometry.mesh.Mesh1D` — beyond ``dataclasses``.  A file in
``geometry/`` that names no geometry is the whole finding the un-weld arc was
built on.

⛔ **The docstring's own name is a leftover.**  ``reduced_operator`` described
an object that left at P4.4; the module keeps the filename only because
deleting it is P4.3's job, one step away.

.. warning::

   ``StreamingTerms`` is **not** "purely geometric", whatever an older
   revision of this file claimed.  ``mu``/``abs_mu`` are per-*direction*, and
   ``delta_A_over_w`` divides a geometric increment by a **quadrature
   weight** — a per-direction quantity over a quadrature weight is posing, not
   geometry.  That observation is what moved this module's other tenants out.
   ⚠ ``delta_A_over_w`` is itself scheduled to retire at P4.7: `[M]` it is
   built per ``(cell, direction)`` *inside* the sweep loop, so unlike
   :attr:`~orpheus.sn.sweep.cache.StreamingCoefficientCache.dA_w`
   (an ``(N, nx)`` array built once at solver init) the fusion buys no
   performance — and neither of its two sides owns it (Pattern 2).

See also
========

* :mod:`orpheus.sn.mesh.reduced_operator` — the chart connection whose
  ``streaming_terms()`` is this packet's only producer.
* :mod:`orpheus.transport.spatial.cell_balance` and
  :mod:`orpheus.transport.spatial.scheme` — its two consumers, and its home
  after P4.3.
* :mod:`orpheus.sn.angular.redistribution` — the angular factor that left here
  at P4.2, carrying the :math:`\alpha` mathematics and the literature
  references that used to head this file.
"""

from __future__ import annotations

from dataclasses import dataclass



# ═══════════════════════════════════════════════════════════════════════
# StreamingTerms — geometry-dependent shape per (cell, direction)
# ═══════════════════════════════════════════════════════════════════════

@dataclass(frozen=True, slots=True)
class StreamingTerms:
    """Per-(cell, direction) **purely geometric** inputs for a sweep cell update.

    All curvature fields are populated for **every** geometry
    (Issue #196 Phase G Step 2.5):

    * **Slab**: neutral-curvature values — ``face_area_inner =
      face_area_outer = 1.0``, ``delta_A_over_w = 0.0``.  The
      Morel–Montry angular weight is NOT carried here (Issue #236
      Step C — it is closure-owned); the neutral slab closure
      (IdentityAngularClosure) supplies τ = 1, α = 0 and stamps the
      derived c on the CellVisit, but ``upstream_state.angular_upstream``
      is ``None`` for slab so the M-M contribution never engages.  Plus
      the always-populated ``chord_length``, ``mu``, ``volume``,
      ``abs_mu``.
    * **Sphere / cylinder**: physically-populated curvature fields
      from the dome recursion (the M-M angular weight is closure-owned,
      not carried here — Issue #236 Step C).

    Cylindrical ``mu`` / ``abs_mu`` are read from the global ordinate
    ``mu_x[level_indices[mu_level_idx][direction_idx]]`` because
    cylindrical ``direction_idx`` is the within-level azimuthal index
    :math:`m \\in [0, M)`, not the global ordinate.

    Before Step 2.5, slab left the curvature fields as ``None`` and
    cell-update strategies branched on ``alpha_in is None`` to
    discriminate slab from curvilinear.  Step 2.5 retires that
    branch — slab carries neutral curvature so the unified
    cell-balance helper consumes the same packet regardless of
    geometry.  Geometry is data, not control-flow (Cardinal Rule 2
    + Pattern 4: make illegal states unrepresentable).

    Geometric, not direction-resolved
    =================================

    Per Cardinal Rule 2 (architecture), this primitive carries
    **purely geometric** labels — labels that are independent of the
    sweep's marching direction.  The two face-area fields are named
    by their geometric position relative to :math:`r=0`:

    * :attr:`face_area_inner` is :math:`A_{i-1/2}` — the face closer
      to the centre of the geometry (smaller :math:`r`).
    * :attr:`face_area_outer` is :math:`A_{i+1/2}` — the face farther
      from the centre (larger :math:`r`).

    These labels do **not** depend on which way the sweep is
    marching.  For an outward sweep (centre → boundary,
    :math:`\\mu > 0`) the inner face is the upstream / incoming face;
    for an inward sweep (boundary → centre, :math:`\\mu < 0`) the
    outer face is the upstream / incoming face.  The
    sweep-direction resolution — *which* of the two faces is
    "downstream" for a given visit — lives in the SN module
    (:class:`~orpheus.transport.spatial.scheme.CellVisit` packs the
    geometric :class:`StreamingTerms` together with the
    sweep-resolved :attr:`face_area_downstream`).  This module is
    geometry-only precisely because the sweep-direction concept is
    SN's alone — MoC, CP, diffusion and MC have no sweep.  (⛔ That
    separation is the true claim; the clause "and is reusable by future
    MoC / CP / diffusion" stood here until 2026-08-27 and is withdrawn.
    A geometry-only module is not thereby a shared one — see the ⛔
    block in the module docstring.)

    The trailing ``volume`` and ``abs_mu`` fields are populated by
    **all three factories** so that a downstream
    :class:`~orpheus.transport.spatial.scheme.DiscretizationScheme` strategy
    receives a self-contained per-cell, per-direction packet and need
    not reach back into ``SNMesh`` or the ``Quadrature``.  Every
    surviving curvature field is populated for **every** geometry
    (Step 2.5; slab carries neutral values) as a required ``float``
    (Pattern 4), so cell-update strategies consume the same packet
    regardless of chart — geometry is **data, not control-flow**.  Slab
    vs curvilinear is discriminated downstream by
    ``upstream_state.angular_upstream is None`` (slab has no half-angle
    state), NOT by any field on this geometry packet — Issue #236 Step C
    retired the former M-M ``alpha_in`` / ``alpha_out`` / ``tau_mm``
    fields (the angular closure now owns that data).
    """

    chord_length: float
    """Cell radial width (slab/sphere/cylinder all use ``mesh.widths[i]``)."""

    mu: float
    """Signed primary direction cosine for this ordinate.

    :math:`\\mu` for slab and sphere (axial); :math:`\\eta` for
    cylindrical 1-D radial sweeps (the radial direction cosine,
    with the global ordinate index resolved through
    :attr:`AngularMeasure.level_indices`).  Signed.

    Slab DD's flow direction reads off the sign of this field.
    Cell-update strategies that need only the magnitude use
    ``abs_mu`` instead.
    """

    face_area_inner: float
    """:math:`A_{i-1/2}` — area of the **inner** radial face
    (closer to :math:`r=0`).

    Geometric label, independent of sweep direction.  See class
    docstring "Geometric, not direction-resolved".
    """

    face_area_outer: float
    """:math:`A_{i+1/2}` — area of the **outer** radial face
    (farther from :math:`r=0`).

    Geometric label, independent of sweep direction.  See class
    docstring "Geometric, not direction-resolved".
    """

    delta_A_over_w: float
    """:math:`\\Delta A_i / w_n` — the geometry-redistribution factor."""

    # Issue #236 Step C: the Morel–Montry ``alpha_in`` / ``alpha_out`` /
    # ``tau_mm`` are NO LONGER carried on the per-cell geometry packet.  The
    # angular weight :math:`\\tau` is an angular-scheme property owned by the
    # MorelMontryAngularSweep angular closure; the derived constants
    # ``c_in`` / ``c_out`` and ``tau`` are stamped on
    # :class:`~orpheus.transport.spatial.scheme.CellVisit`.  The α-dome itself
    # survives on :class:`ReducedStreamingOperator` (``alpha_half`` /
    # ``alpha_per_level``).

    # ``mu_start`` retired 2026-08-26.  It was the middle link of a
    # three-link dead chain: AngularRedistribution.mu_start_per_level ->
    # StreamingTerms.mu_start -> StreamingCoefficientCache.mu_start -> nothing.
    # [M] the terminal had ZERO readers of any kind (no attribute access,
    # no getattr by name), so this field's only production consumer was
    # the WRITE into it -- while its own docstring claimed
    # ``MorelMontryAngularSweep`` consumed it.  The closure does consume a
    # starting direction; it reads the OWNER
    # (``AngularRedistribution.mu_start_per_level``), which is why the
    # claim read as true.  The chain was ERR-058's threading; the un-weld
    # gave the datum one owner and the thread became dead weight.

    volume: float
    """Cell volume :math:`V_i`.

    Populated by all three factories from
    :attr:`~orpheus.geometry.mesh.Mesh1D.volumes`: slab uses
    ``mesh.volumes[i]`` (which equals ``widths[i]`` for unit
    cross-section in 1-D Cartesian); sphere uses
    :math:`\\tfrac{4}{3}\\pi(r_{i+1}^3 - r_i^3)`; cylinder uses
    :math:`\\pi(r_{i+1}^2 - r_i^2)` (per unit axial length).  Carried
    on the streaming-terms packet so that the
    :class:`~orpheus.transport.spatial.scheme.DiscretizationScheme` cell-update
    contract receives :math:`V_i` directly without needing access to
    the underlying ``SNMesh``.
    """

    abs_mu: float
    """Absolute primary direction cosine.

    :math:`|\\mu|` for slab and sphere (axial direction cosine);
    :math:`|\\eta|` for cylindrical 1-D radial sweeps (the radial
    direction cosine).  Slab and sphere factories compute this as
    ``abs(quadrature.mu_x[direction_idx])`` — for those geometries
    ``direction_idx`` is the global ordinate index.  The cylindrical
    factory resolves the global ordinate through
    :attr:`AngularMeasure.level_indices` because cylindrical
    ``direction_idx`` is the within-level azimuthal index
    :math:`m \\in [0, M)`, not the global ordinate index.  By
    ORPHEUS convention, ``mu_x`` carries :math:`\\mu` for sphere and
    :math:`\\eta` for cylinder, and both are real-valued, so the
    absolute value is always well-defined.

    The cylindrical pure-azimuthal degenerate case (``abs_mu <
    1e-15``) is signalled by cell-update strategies via
    :attr:`~orpheus.transport.spatial.scheme.CellResult.outgoing_spatial_flux`
    set to ``None`` (no radial face flow on the cell).
    """


__all__ = ["StreamingTerms"]
