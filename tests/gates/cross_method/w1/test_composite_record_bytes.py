r"""G05m: the composite record keeps the base's bytes (N05, spec 7.6 row G05m, decision ON05-11, additions ON05-22).

The record ``composite_probes/composite_probes.npz`` (1195 arrays) was taken on the UNCHANGED base ``ccecdc365`` by ``_composite_probe_population.write_record``,
before any production edit of the unit, over the tracked population builder, which names neither class N05 retires. Its population: the 33 composites of the 15 harness
fixtures (21 full fields, 12 coupled spaces, 6 of them System B ray instances reached through ``radial_characteristic_field_space``), each in both charts, primal and dual
(``apply_metric``, ``apply_inverse_metric`` flat, through the typed carrier, and with a two-column payload), ``rank``, ``kernel_basis``, ``metric_diagonal``,
``inner_product``, ``norm``, the member round trip; and the additions of ON05-22: a lifted full field (the moment bulk of ``angular_lift.py``, the codomain of
``BulkAnalysisOperator``) on the slab and the box, a ``TimedFullField`` with one frame in its history, the kernels of the loss-kernel gauge through its hand slices of the
bulk size (the harness's 12 SN fixtures all predict an empty kernel, 0 summands, so two all-specular 2-D boxes are built by the population builder on which both predictors answer: 112 arrays), the
diffusion solver's predicted kernel (0 summands on the 3 fixtures: the count rows only) and the group of each coordinate of the diffusion solver's full field (taken then by the solver's own method, now by ``full_field_space.axis_index_of_coordinates("energy")``). It is re-taken by no one before N10.

Host-keyed: the arrays pass through LAPACK, so the stored bytes are compared on the stamped host and the gate FAILS on a foreign stamp with both stamps in the message
(a skip would be a silent green); it is not in the CI step. GREEN at the base tip by construction (the record is the base's, two captures byte-equal); the input that
reddens it is in the unit's tree: ``reverse the summand order of the sum's pack`` moves the metric rows of the typed chart and of the flat chart together (the record
holds the 42 arrays of the harness's (d) group, which the sum's metric must reach in slot order, and 585 more).

Three answers change by ruling: ``unit``, ``integral`` and ``has_coordinate_cone`` of a sum (ON05-10, ON05-25). The record holds no key for them; the key arm parses the
manifest, with a planted key as its positive control. The harness shares 42 arrays with the record (the ``d`` group of the composite spaces); the control row reads 42 of 42.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from tests.gates.cross_method.w1 import _composite_probe_population as pop
from tests.gates.cross_method.w1 import _snapshot as snap

pytestmark = [pytest.mark.foundation, pytest.mark.rests_on("tests/gates/cross_method/w1/test_snapshot_harness.py")]

RECORD_DIR = Path(__file__).resolve().parent / "composite_probes"
HARNESS_DIR = Path(__file__).resolve().parent / "snapshots"
RECORD_COUNT = 1195


def _require(cond: object, msg: str) -> None:
    if not cond:
        pytest.fail(msg)


def _stamp_must_match() -> None:
    stamp = tuple((RECORD_DIR / "STAMP.txt").read_text().splitlines()[0].split("\t"))
    live = snap.host_stamp()
    _require(stamp == live, f"the composite record cannot be read on this host: it was taken on {stamp} and this host is {live}. Re-take it on this host at the last green commit before comparing.")


def test_the_composite_record_has_the_bases_bytes() -> None:
    _stamp_must_match()
    record = np.load(RECORD_DIR / "composite_probes.npz")
    live = {f"{f}::{n}": a for (f, n), a in pop.capture_composite_probes().items()}
    _require(len(record.files) == RECORD_COUNT, f"the record holds {len(record.files)} arrays, expected {RECORD_COUNT}")
    _require(sorted(live) == sorted(record.files), f"the population's keys differ from the record's: {sorted(set(live) ^ set(record.files))[:6]}")
    bad = [
        k
        for k in sorted(live)
        if np.ascontiguousarray(live[k]).dtype != record[k].dtype or live[k].shape != record[k].shape or np.ascontiguousarray(live[k]).tobytes() != record[k].tobytes()
    ]
    _require(not bad, f"{len(bad)} of {len(live)} composite arrays differ from the base's bytes: {bad[:6]}")


def test_the_record_holds_no_key_for_the_three_answers_that_change_by_ruling() -> None:
    manifest_keys = [tuple(line.split("\t")[:2]) for line in (RECORD_DIR / "MANIFEST.txt").read_text().splitlines()]
    npz_keys = [tuple(k.split("::")) for k in np.load(RECORD_DIR / "composite_probes.npz").files]
    _require(len(manifest_keys) == RECORD_COUNT and sorted(manifest_keys) == sorted(npz_keys), "the manifest and the arrays name different keys")
    found = pop.record_keys_are_admissible(manifest_keys)
    _require(not found, f"the record holds keys for answers that change by ruling (ON05-10, ON05-25): {found[:4]}")


def test_the_key_predicate_reads_a_planted_key() -> None:
    # Positive control of the arm above: the predicate finds each forbidden verb in a planted manifest, and nothing in a clean one.
    for verb in pop.FORBIDDEN_VERBS:
        planted = [("sn_slab_ng1.full_field_space", f"primal.{verb}")]
        _require(pop.record_keys_are_admissible(planted) == planted, f"the key predicate does not read a planted {verb!r} key")
    _require(pop.record_keys_are_admissible([("f", "primal.apply_metric.flat")]) == [], "the key predicate reads a clean key as forbidden")


def test_the_record_shares_42_arrays_with_the_harness() -> None:
    # The harness's (d) group on a composite space is the record's flat apply_metric / apply_inverse_metric of the same seeded vector: 42 of 42 equal.
    _stamp_must_match()
    _stamp, lines = snap.read_manifest(HARNESS_DIR / "MANIFEST.txt")
    harness = {(f, n): sha for f, n, _shape, _dtype, sha in (line.split("\t") for line in lines)}
    record = np.load(RECORD_DIR / "composite_probes.npz")
    shared = equal = 0
    for key in record.files:
        fixture_and_space, quantity = key.split("::")
        fixture, _, space = fixture_and_space.partition(".")
        harness_space = {"full_field_space": "full_field_space", "radial_field_space": "radial_characteristic_field_space"}.get(space)
        for verb in ("apply_metric", "apply_inverse_metric"):
            if harness_space is not None and quantity == f"primal.{verb}.flat":
                shared += 1
                equal += harness.get((fixture, f"d.{harness_space}.{verb}")) == snap.sha1(record[key])
    _require(shared == 42 and equal == 42, f"the record shares {shared} arrays with the harness and {equal} are byte-equal; expected 42 of 42")
