r"""#709 (W3): a moment field's ``L`` and ``spatial_moments`` are the ones its space states.

``MomentField`` stores ``L`` and ``spatial_moments`` beside the ``space`` that already fixes both
(``head_layout.L`` and ``spatial_moments_per_axis``). Today any pair constructs. The refusal is ONE
``raise ValueError`` in ``MomentField.__post_init__`` tagged ``ELEGANCE-DEBT[guard] #709 (W3)``;
W3 retires it by deriving both from the space. Specification: ``.claude/plans/ontology_w0_spec.md``
5.1 to 5.3 (unit ``defects``).

Claim layer: construction law of a type (L0, a refusal). Kind: THEOREM for the admitted legs (for every
admitted input ``L`` and ``spatial_moments`` are the space's), RECORD of nothing.

Activation: the space is the slab moment space at order 1 (``head_layout.L == 1``,
``spatial_moments_per_axis == 1``), so ``L=5`` and ``spatial_moments=3`` each disagree with it on exactly
one of the two predicates; the two legs are disjoint (a build that drops predicate B still refuses ``L=5``).

Legs (5.3): (1) each refused input in a cold ``python -O`` child; (1b) the same in process under
``pytest.raises(ValueError, match=...)``; (2) the admitted twin; (3) the ``TypeError`` of a space that is not a
tensor product. Mutations (each first reads GREEN on the unmutated tree, so a reading that is red before the
mutation proves nothing): the refusal replaced by ``assert False, (`` in a copy (leg 1 reddens: the input
constructs under ``-O``); predicate B dropped in process (only the ``spatial_moments`` leg reddens).

Every assertion is a function call or an ``if ... pytest.fail`` (``python -O`` safe).
"""

from __future__ import annotations

import dataclasses
import pathlib

import numpy as np
import pytest

from orpheus.numerics.space import FunctionSpace
from orpheus.transport.fields import HarmonicMomentSourceSink
from orpheus.transport.fields._bases import MomentField
from orpheus.transport.fields.harmonic_moment_flux import HarmonicMomentFlux
from tests.gates import _refusal_gates as rg
from tests.gates.transport._phase_space_cases import moment_space, posing_quadrature, slab_mesh
from orpheus.numerics.factor import Kind, Variance, VarianceChart
from orpheus.numerics.space import TensorSpace

pytestmark = pytest.mark.foundation

FRAGMENT = "disagrees with its space"

_CHILD_IMPORTS = (
    "import numpy as np\n"
    "from orpheus.transport.fields.harmonic_moment_flux import HarmonicMomentFlux\n"
    "from tests.gates.transport._phase_space_cases import moment_space, posing_quadrature, slab_mesh\n"
    "mesh = slab_mesh()\n"
    "sp = moment_space(mesh, posing_quadrature(mesh), 1)\n"
)
_CHILD_INPUTS = {
    "L": "HarmonicMomentFlux(values=np.zeros(sp.shape), space=sp, L=5)",
    "spatial_moments": "HarmonicMomentFlux(values=np.zeros(sp.shape), space=sp, L=1, spatial_moments=3)",
}
_ADMITTED_TWIN = "HarmonicMomentFlux(values=np.zeros(sp.shape), space=sp, L=1)"
_BASES_FILE = "orpheus/transport/fields/_bases.py"

LEAVES = [HarmonicMomentFlux, HarmonicMomentSourceSink]


@pytest.fixture(scope="module")
def space():
    mesh = slab_mesh()
    return moment_space(mesh, posing_quadrature(mesh), 1)


# ----------------------------------------------------------------- (2) admitted legs


@pytest.mark.parametrize("leaf", LEAVES, ids=lambda c: c.__name__)
def test_the_admitted_twin_constructs_and_reads_the_spaces_values(leaf, space):
    """CONTROL: the space states ``L = 1`` and one spatial moment; the field built with those is admitted."""
    field = leaf(values=np.zeros(space.shape), space=space, L=1)
    if (field.head_layout.L, field.spatial_moments_per_axis) != (1, 1):
        pytest.fail("CONTROL INVALID: the fixture space does not state L = 1 and one spatial moment")
    if (field.L, field.spatial_moments) != (1, 1):
        pytest.fail("the admitted twin lost its own L or spatial_moments")


# ----------------------------------------------------------------- (1b) in process


@pytest.mark.parametrize("leaf", LEAVES, ids=lambda c: c.__name__)
def test_a_field_whose_L_disagrees_with_its_space_is_refused(leaf, space):
    """First red: ``L=5`` against a space that says 1 constructs today. The message names the class and both pairs."""
    with pytest.raises(ValueError, match=FRAGMENT) as info:
        leaf(values=np.zeros(space.shape), space=space, L=5)
    message = str(info.value)
    for named in (leaf.__name__, "head_layout.L", "spatial_moments_per_axis"):
        if named not in message:
            pytest.fail(f"the refusal does not name {named!r}: {message}")


@pytest.mark.parametrize("leaf", LEAVES, ids=lambda c: c.__name__)
def test_a_field_whose_spatial_moments_disagrees_with_its_space_is_refused(leaf, space):
    """First red: ``spatial_moments=3`` against a space that says 1 constructs today."""
    with pytest.raises(ValueError, match=FRAGMENT):
        leaf(values=np.zeros(space.shape), space=space, L=1, spatial_moments=3)


def test_replace_re_runs_the_check(space):
    """``dataclasses.replace`` re-runs ``__post_init__``, so a changed ``L`` on an admitted field is refused."""
    field = HarmonicMomentFlux(values=np.zeros(space.shape), space=space, L=1)
    with pytest.raises(ValueError, match=FRAGMENT):
        dataclasses.replace(field, L=5)


# ----------------------------------------------------------------- (3) the TypeError leg


def test_a_space_that_is_not_a_tensor_product_is_now_refused_at_construction_with_TypeError():
    """The #709 refusal reads ``head_layout``, so a space with no head is refused at construction (``TypeError`` from
    ``MomentField.head``); before #709 such a field constructed. W3's retirement of the tag keeps or re-homes this."""
    odd = TensorSpace.declared("moments", (2, 3, 2, 4), Variance(Kind.FUNCTION, VarianceChart.DENSITY))
    with pytest.raises(TypeError):
        HarmonicMomentFlux(values=np.zeros((2, 3, 2, 4)), space=odd, L=1)


# ----------------------------------------------------------------- (1) cold child, -O


def test_the_unmutated_child_constructs_the_admitted_twin():
    """CONTROL of the child legs: the twin imports and constructs under ``python -O`` (the harness runs)."""
    child = rg.run_child(_CHILD_IMPORTS, _ADMITTED_TWIN)
    if not rg.imported_from(child, rg.import_root()):
        pytest.fail(f"the child did not import orpheus from {rg.import_root()}:\n{child.stdout}\n{child.stderr}")
    if not rg.constructed(child, FRAGMENT):
        pytest.fail(f"the admitted twin was refused in the child:\n{child.stderr}")


@pytest.mark.parametrize("which", sorted(_CHILD_INPUTS))
def test_the_refusal_survives_python_O(which):
    """The refusal is a ``raise``: under ``-O`` the input is refused with the fragment on stderr."""
    child = rg.run_child(_CHILD_IMPORTS, _CHILD_INPUTS[which])
    if not rg.imported_from(child, rg.import_root()):
        pytest.fail(f"the child did not import orpheus from {rg.import_root()}:\n{child.stdout}\n{child.stderr}")
    if not rg.refused(child, FRAGMENT):
        pytest.fail(f"the {which!r} input was not refused under -O (exit {child.returncode}):\n{child.stderr}")


# ----------------------------------------------------------------- mutations


@pytest.mark.parametrize("which", sorted(_CHILD_INPUTS))
def test_mutation_the_refusal_replaced_by_an_assert_constructs_under_python_O(which, tmp_path):
    """Mutation of leg (1): in a copy, the ``raise`` after the #709 tag becomes ``assert False, (``; under ``-O``
    the refused input constructs again and the leg's reading flips (exit 0, fragment absent).

    The honest reading is taken first on the unmutated tree, so a leg that was red before the mutation cannot
    pass for a catch."""
    honest = rg.run_child(_CHILD_IMPORTS, _CHILD_INPUTS[which])
    if not rg.refused(honest, FRAGMENT):
        pytest.fail(f"the unmutated tree does not refuse {which!r}, so the mutation proves nothing:\n{honest.stderr}")
    root = rg.copy_orpheus_py(rg.import_root(), tmp_path)
    rg.remove_refusal(root, _BASES_FILE, 709)
    mutant = rg.run_child(_CHILD_IMPORTS, _CHILD_INPUTS[which], root=root)
    if not rg.imported_from(mutant, root):
        pytest.fail(f"the mutant did not import from the copy (never red evidence):\n{mutant.stdout}\n{mutant.stderr}")
    if not rg.constructed(mutant, FRAGMENT):
        pytest.fail(f"the mutant still refused {which!r} (exit {mutant.returncode}):\n{mutant.stderr}")


def test_mutation_dropping_predicate_B_reddens_only_the_spatial_moments_leg(space, monkeypatch):
    """In process: ``spatial_moments_per_axis`` made to read the field itself drops predicate B; the
    ``spatial_moments=3`` input then constructs, the ``L=5`` input is still refused (the legs are disjoint)."""
    with pytest.raises(ValueError, match=FRAGMENT):
        HarmonicMomentFlux(values=np.zeros(space.shape), space=space, L=1, spatial_moments=3)
    monkeypatch.setattr(MomentField, "spatial_moments_per_axis", property(lambda self: self.spatial_moments))
    HarmonicMomentFlux(values=np.zeros(space.shape), space=space, L=1, spatial_moments=3)
    with pytest.raises(ValueError, match=FRAGMENT):
        HarmonicMomentFlux(values=np.zeros(space.shape), space=space, L=5)


def test_this_file_runs_on_the_tree_it_names():
    """The unmutated legs read the tree pytest imported ``orpheus`` from; a copy run is a different ``import_root``."""
    if not (rg.import_root() / "orpheus" / "transport" / "fields" / "_bases.py").is_file():
        pytest.fail(f"import_root {rg.import_root()} holds no orpheus/transport/fields/_bases.py")
    if not pathlib.Path(rg.__file__).is_file():
        pytest.fail("the shared helper is missing")
