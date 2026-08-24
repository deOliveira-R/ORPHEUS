r"""B3.2 mutation harness — an in-process pytest plugin proving gate teeth.

Usage (SERIAL, canonical ``-O``)::

    ORPHEUS_B32=N1 .venv/bin/python -O -m pytest tests/sn/operators \
        -p no:randomly -p tests._harness.mutation_batteries.b3_2_boundary -q

``ORPHEUS_B32`` unset ⇒ the CONTROL leg: nothing is patched, everything must be
GREEN. Any other value installs the named mutation by monkeypatching
``SNBoundaryOperator._reflect_trace`` (or the realizer) IN PROCESS — no file on
disk is touched, so this can never destroy uncommitted work
(``process-discipline``: never ``git checkout`` a file you hold edits in).

Why this file is IN THE REPO
============================

The 31-mutation harness the B3 plan's acceptance gate names ("the sweep still
catches 30/31") lived only in a job-scratch directory and evaporated with the
session — so the gate was unrunnable by the next reader. Every mutation that
justifies a B3.2 gate's teeth lives here instead, next to the memo that reports
its measured colour (``scratch/b3_2_migration.md`` §2).

Moved out of ``scratch/`` at **B3.5** (2026-08-14): tracked-but-in-a-holding-pen
is not discoverable, and `[M]` the evaporation failure had by then happened a
SECOND time — ``tests/sn/operators/__pycache__/`` still carries
``conftest_mutate_kernel.*.pyo`` and ``conftest_mutate_pr.*.pyo`` whose sources
were never tracked at all. See ``README.md`` in this package for the mechanism's
boundary against ``tests/_mutation/`` (cosmic-ray), which is a different tool
answering a different question and whose revert step is a ``git checkout``.

Each mutation is the plausible transcription of a real B3.2 hazard, NOT an
arbitrary perturbation:

======  =====================================================================
id      the wrong code it simulates
======  =====================================================================
N1      write the law's image to the OUTFLOW rows (wrong write target)
N2      additionally pass ``γ₊`` through onto the output's outflow rows
N3      hand the law the WRONG half — ``γ₋`` where ``γ₊`` is meant. Spellable
        and silent because ``|Γ₋| == |Γ₊|`` on every reachable fixture
N4      hand the law the FULL face (the pre-B3.2 domain)
N5      leak onto the TANGENTIAL rows — the class the pre-B3.2 ``got[outflow]``
        leg was structurally blind to
N7      the transpose writes onto ``Γ₋`` instead of ``Γ₊``
N8      the split remap uses ``arange(sel.size)`` instead of ``to_local(sel)``
        — exactly correct in 1-D, wrong in 2-D
N9      the REFLECTIVE remap uses ``arange`` instead of ``γ₊.to_local`` — the
        SIBLING trap, and its discriminating fixture is the opposite one: a
        slab, where the mirror reverses order
M1      the reflective permutation collapses to the identity
M2      the reflective permutation is rolled by one
======  =====================================================================

``N8`` and ``N9`` together are why the crosswalk insists BOTH gates ship: the
same naive ``arange`` appears at two sites, and **neither site's fixture
covers the other's** (N8 is vacuous in 1-D, N9 is vacuous wherever
``perm[inflow]`` happens to be ``outflow`` in order — e.g. the cylinder).

A mutation is only evidence if it BIT (vv Mode-8 method warning). Every run
prints a one-line fingerprint of the composite's output so a no-op patch is
visible rather than being read as "the gate is blind".
"""
from __future__ import annotations

import hashlib
import os

import numpy as np

_MUT = os.environ.get("ORPHEUS_B32", "").strip()


def _mutated_reflect_trace(mut: str):
    """Return a ``_reflect_trace`` carrying mutation ``mut``.

    Transcribed from ``orpheus/sn/operators/boundary.py::_reflect_trace`` at
    B3.2 so each mutation is a ONE-LINE difference from production and the
    diff is readable. Kept deliberately close to the original body; the
    mutation points are marked ``# MUT``.
    """

    def _reflect_trace(self, boundary, method, faces=None, rows=None):
        from orpheus.transport.source_sinks import AngularBoundarySourceSink

        mesh = self.sn_mesh
        trace = mesh.angular_trace
        out_boundary = AngularBoundarySourceSink.zeros_on(mesh)
        face_laws = self._face_laws
        if faces is not None:
            unknown = set(faces) - set(face_laws)
            if unknown:
                raise ValueError(
                    f"_reflect_trace: face(s) {sorted(unknown)} are not "
                    f"boundary faces of this mesh; available faces: "
                    f"{sorted(face_laws)}."
                )
            face_laws = {f: face_laws[f] for f in faces}
        if rows is not None:
            face_laws = {f: law for f, law in face_laws.items() if f in rows}
        for face, law in face_laws.items():
            face_in = boundary.face_view(face)
            gamma_out = trace.outflow_restriction(face)
            gamma_in = trace.inflow_restriction(face)
            if method == "apply":
                if mut == "N3":
                    argument = gamma_in.apply(face_in)          # MUT wrong half
                elif mut == "N4":
                    argument = face_in                          # MUT full face
                else:
                    argument = gamma_out.apply(face_in)
                image = law.apply(argument)
                if rows is None:
                    if mut == "N1":
                        # MUT: scatter over Γ₊ instead of Γ₋.
                        out_boundary.face_view(face)[...] = (
                            gamma_out.apply_transpose(image)
                        )
                    else:
                        out_boundary.face_view(face)[...] = (
                            gamma_in.apply_transpose(image)
                        )
                    if mut == "N2":
                        # MUT: additionally pass γ₊ through to the outflow rows.
                        out_boundary.face_view(face)[gamma_out.indices] = (
                            gamma_out.apply(face_in)
                        )
                    if mut == "N5":
                        # MUT: leak onto the tangential rows.
                        n_face = face_in.shape[0]
                        tangent = np.setdiff1d(
                            np.arange(n_face),
                            np.union1d(gamma_in.indices, gamma_out.indices),
                        )
                        if tangent.size:
                            out_boundary.face_view(face)[tangent] = (
                                face_in[tangent]
                            )
                else:
                    sel = rows[face]
                    if mut == "N8":
                        pos = np.arange(np.asarray(sel).size)   # MUT naive remap
                    else:
                        pos = gamma_in.to_local(sel)
                    out_boundary.face_view(face)[sel] = image[pos]
            else:
                from orpheus.numerics.operator import MissingAdjoint

                if not getattr(law, "is_adjointable", False):
                    raise MissingAdjoint(
                        f"boundary face {face!r}: the realized law "
                        f"{type(law).__name__} has no Euclidean transpose"
                    )
                if rows is None:
                    restricted = gamma_in.apply(face_in)
                else:
                    sel = rows[face]
                    masked = np.zeros_like(face_in)
                    masked[sel] = face_in[sel]
                    restricted = gamma_in.apply(masked)
                image = law.apply_transpose(restricted)
                if mut == "N7":
                    # MUT: scatter the transpose image over Γ₋, not Γ₊.
                    out_boundary.face_view(face)[...] = (
                        gamma_in.apply_transpose(image)
                    )
                else:
                    out_boundary.face_view(face)[...] = (
                        gamma_out.apply_transpose(image)
                    )
        return out_boundary

    return _reflect_trace


def _mutate_realizer(mut: str) -> None:
    """Patch the reflective branch's reduced permutation."""
    from orpheus.numerics import operator as _op

    original = _op.PermutationOperator.__init__

    def patched(self, perm, axis=0, **kw):
        perm = np.asarray(perm)
        if mut == "M1":
            perm = np.arange(perm.size)
        elif mut == "M2":
            perm = np.roll(perm, 1)
        original(self, perm, axis=axis, **kw)

    _op.PermutationOperator.__init__ = patched


def _mutate_to_local() -> None:
    """N9 — the reflective narrowing reaches for ``arange`` instead of the
    half-trace space's own ``to_local``.

    Patched on the OWNER of the remap rather than in the realizer body, so the
    mutation is EXACTLY the code a call site would have written by hand if the
    remap had not been owned (crosswalk §9: owning it is what makes this
    failure mode unspellable in production — this harness re-spells it to
    measure what the gates would have caught).

    ⛔ **REPOINTED 2026-08-14 (B3.5), and it had been silently INERT.** This
    patched ``TraceRestrictionOperator.to_local`` until now. G6.5 (``0d99140c``)
    moved the remap onto :class:`AngularFaceTraceSpace` — "the operator is the
    *arrow*, but which global row sits where is a fact about the SUBSPACE" —
    and a bare class-attribute assignment **creates** a missing attribute
    rather than failing, so this leg installed a method nobody calls.
    `[M]` before the repoint: ``ORPHEUS_B32=N9`` gave *1167 passed*, byte-
    identical to CONTROL. A reader would have taken that green as "the gates
    are blind to the reflective-remap hazard" — a harness lying in the
    safe-looking direction (``vv-principles`` #17).

    The existence check below is what stops that recurring, and it is a
    ``raise`` rather than an ``assert`` **on purpose**: the canonical runner is
    ``python -O``, which strips ``assert`` at compile time, so an asserted
    guard here would be exactly as inert as the bug it guards against
    (``coding-standards``: a bare assert is not a contract).
    """
    from orpheus.numerics.spaces.angular_trace_space import (
        AngularFaceTraceSpace,
    )

    if not hasattr(AngularFaceTraceSpace, "to_local"):
        raise RuntimeError(
            "N9 cannot install: AngularFaceTraceSpace has no `to_local`. The "
            "remap has moved again — find its new owner and repoint this "
            "mutation. Do NOT leave this leg patching a name that is not "
            "there: assigning it would CREATE a dead attribute and the leg "
            "would run green, which reads as 'the gates are blind' rather "
            "than 'the mutation never happened'."
        )

    AngularFaceTraceSpace.to_local = (
        lambda self, global_rows: np.arange(np.asarray(global_rows).size)
    )


def pytest_configure(config) -> None:  # noqa: D401 - pytest hook
    if not _MUT:
        print("\n[b3_2_mutations] CONTROL leg — no mutation installed.")
        return
    from orpheus.sn.operators.boundary import SNBoundaryOperator

    if _MUT in {"M1", "M2"}:
        _mutate_realizer(_MUT)
    elif _MUT == "N9":
        _mutate_to_local()
    else:
        SNBoundaryOperator._reflect_trace = _mutated_reflect_trace(_MUT)
    print(f"\n[b3_2_mutations] variant={_MUT} installed; bite check: "
          f"{_bite_fingerprint()}")


def _bite_fingerprint() -> str:
    """Positive control — prove the patch CHANGED production numbers.

    A mutation that silently fails to install manufactures a false "this gate
    is blind" verdict (the vv Mode-8 METHOD WARNING). Printing a hash of the
    composite's forward and transpose output on a fixed fixture makes a no-op
    patch visible in the run's first line.
    """
    try:
        from dataclasses import replace

        from orpheus.geometry import BC, Mesh1D, Region, RegionMesh, StructuredGeometry
        from orpheus.numerics.quadrature import Quadrature
        from orpheus.sn.mesh.augmented_mesh import SNMesh
        from orpheus.sn.operators.boundary import SNBoundaryOperator
        from orpheus.transport.fields.angular_boundary_flux import AngularBoundaryFlux
        from orpheus.transport.fields.angular_flux import AngularFlux
        from orpheus.transport.timed_full_field import TimedFullField
        from tests.sn._test_helpers import placeholder_materials

        geom = StructuredGeometry(
            geometry="SLB",
            regions=(Region(mat_id=0, outer_thickness_cm=2.0),),
            bcs=(BC.vacuum, BC.reflective),
        )
        mesh = Mesh1D.from_geometry(geom, region_meshes=(RegionMesh(n_cells=4),))
        sn = SNMesh(mesh, Quadrature.gauss_legendre(4), placeholder_materials(ng=1))
        z = TimedFullField.zeros(
            interior=AngularFlux, boundary=AngularBoundaryFlux, space=sn.full_field_space,
        )
        rng = np.random.default_rng(7)
        psi = replace(
            z,
            boundary=replace(
                z.boundary,
                values=rng.uniform(0.5, 2.0, size=z.boundary.values.shape),
            ),
        )
        B = SNBoundaryOperator(sn)
        fwd = B.apply(psi).boundary.values
        tra = B.apply_transpose(psi).boundary.values
        return (
            f"FWD sha={hashlib.sha256(fwd.tobytes()).hexdigest()[:16]} "
            f"TRA sha={hashlib.sha256(tra.tobytes()).hexdigest()[:16]}"
        )
    except Exception as exc:  # pragma: no cover - diagnostic only
        return f"<bite check unavailable: {type(exc).__name__}: {exc}>"
