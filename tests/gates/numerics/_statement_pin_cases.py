"""The fixture table of the per-statement import-time pins (the cone store's increment B2, issue #689).

Spec: ``scratch/verification_tree/cone_store/b2/gate_spec.md``. Each :class:`Case` is a small package ``fx`` written
under a temporary root, an entry ``fxrun.run`` whose body is the case's ``run``, and one edit. The gate
(``test_traced_memo_statement_pins.py``) traces ``fxrun.run`` with :func:`~orpheus.numerics.traced_memo.trace_call` in
one fresh process, applies the edit, and validates the manifest in a second fresh process (the memo's source reads are
kept per process, so an in-process validation after an on-disk edit could read the bytes from before it). ``run``
imports ``fx`` INSIDE the traced call, so the modules' import-time code runs under the recorder, as it does in a
generating process.

Every expectation is written from the sources below by hand (no recorder decides it). Every name a case binds is
spelled so that it cannot collide with a name the import machinery, the standard library or numpy reads (the rule's
names are not qualified by module, so a common name such as ``path`` would be read by the import system itself).

Kinds of row, by the side of the rule each pins:

* ``recall``: the edit changes what the run would return or do, and only an import-time statement pin can see it (no
  def that ran changes). The manifest must go stale. Green at B1 (the module skeleton moves) and after B2.
* ``precision``: the edit is to a statement the run cannot observe. The manifest must stay valid. Red at B1 (`e9df7ed9`) (the
  module skeleton moves); green after B2.
* ``under``: a declared boundary of the design (item 4) where B2 knowingly under-selects. The documented behaviour
  (valid after an edit that changes the run) is asserted. Red at B1 (`e9df7ed9`); green after B2; a later repair turns it red,
  and the row is then re-posed as ``recall``.
* ``over``: a declared over-selection (the statement is the unit; names are not module-qualified; effect statements
  keep module grain). Stale, though the run cannot observe the edit. Green at B1 and after B2.
* ``finding``: a recall row the design as written misses (the spec's findings section, each with the prototype's
  ``[M]`` verdict). Green at B1; red after a B2 built to the letter of the design until the design repairs it.
"""
from __future__ import annotations

import textwrap
from dataclasses import dataclass, field

#: The package every case writes under its root, beside the entry module ``fxrun``.
PACKAGE = "fx"


@dataclass(frozen=True)
class Edit:
    """Replace ``old`` by ``new`` exactly once in ``file`` (relative to the root). An edit whose ``old`` is not there
    exactly once is refused by the gate: an edit that did not apply reads as a correct precision verdict."""

    file: str
    #: ``None``: ``new`` is the file's whole new text (qa's rows).
    old: str | None
    new: str


@dataclass(frozen=True)
class Case:
    id: str
    kind: str
    files: dict[str, str]
    run: str
    edit: Edit
    #: What the edit does to the run's answer, by hand: the evidence that a ``recall`` row's edit is observable.
    effect: str
    #: The design clause or census row the case realises.
    shape: str
    #: Python-level value ``run`` returns before and after the edit, when the edit changes it (``None``: the case
    #: does not assert behaviour, e.g. a precision row, or the change is a raise).
    before: str | None = None
    after: str | None = None
    extra: dict[str, str] = field(default_factory=dict)

    def tree(self) -> dict[str, str]:
        body = textwrap.indent(textwrap.dedent(self.run).strip("\n"), "    ")
        files = {f"{PACKAGE}/__init__.py": "", **{k: textwrap.dedent(v).lstrip("\n") for k, v in self.files.items()}}
        files["fxrun.py"] = f"def run():\n{body}\n"
        return files


def _c(id: str, kind: str, shape: str, files: dict[str, str], run: str, edit: Edit, effect: str,
       before: str | None = None, after: str | None = None) -> Case:
    return Case(id, kind, files, run, edit, effect, shape, before, after)


CASES: tuple[Case, ...] = (
    # ── precision: constants A and B (the brief's memo-level control) ─────────────────────────────────────────
    _c("P1-unread-constant", "precision", "design step 2: a statement is selected only when a name it binds is read",
       {"fx/consts.py": """
            ALPHA_READ_TOL = 1e-15
            BETA_UNREAD_TOL = 2e-15
        """,
        "fx/reader.py": """
            from fx import consts


            def above_alpha(x):
                return x > consts.ALPHA_READ_TOL
        """},
       "from fx import reader\nreturn reader.above_alpha(1.0)",
       Edit("fx/consts.py", "BETA_UNREAD_TOL = 2e-15", "BETA_UNREAD_TOL = 3e-15"),
       "none: the run never reads BETA_UNREAD_TOL"),
    _c("R1-read-constant", "recall", "design step 1: a def that ran names the constant (co_names)",
       {"fx/consts.py": """
            ALPHA_READ_TOL = 1e-15
            BETA_UNREAD_TOL = 2e-15
        """,
        "fx/reader.py": """
            from fx import consts


            def above_alpha(x):
                return x > consts.ALPHA_READ_TOL
        """},
       "from fx import reader\nreturn reader.above_alpha(1.5e-15)",
       Edit("fx/consts.py", "ALPHA_READ_TOL = 1e-15", "ALPHA_READ_TOL = 2e-15"),
       "above_alpha(1.5e-15) turns False", "True", "False"),
    _c("P2-unread-class-attribute", "over", "ruling R5: a selected class header selects its body's non-def statements",
       {"fx/widgets.py": """
            class Widget:
                THETA_READ_LIMIT = 3
                THETA_UNREAD_LIMIT = 9
        """,
        "fx/factory.py": """
            from fx.widgets import Widget


            def build_widget():
                return Widget()
        """},
       "from fx import factory\nreturn factory.build_widget().THETA_READ_LIMIT",
       Edit("fx/widgets.py", "THETA_UNREAD_LIMIT = 9", "THETA_UNREAD_LIMIT = 10"),
       "none"),
    _c("P3-unread-derived-constant", "precision", "design step 2: the closure follows reads of selected statements only",
       {"fx/derived.py": """
            KAPPA_BASE_RATE = 3
            KAPPA_DOUBLED_RATE = KAPPA_BASE_RATE * 2
            KAPPA_TRIPLED_RATE = KAPPA_BASE_RATE * 3


            def doubled_rate():
                return KAPPA_DOUBLED_RATE
        """},
       "from fx import derived\nreturn derived.doubled_rate()",
       Edit("fx/derived.py", "KAPPA_TRIPLED_RATE = KAPPA_BASE_RATE * 3", "KAPPA_TRIPLED_RATE = KAPPA_BASE_RATE * 4"),
       "none"),
    _c("P4-reorder-unselected", "precision", "design step 3: the digest is of the ordered SELECTED statements",
       {"fx/order.py": """
            LAMBDA_READ_VALUE = 1
            MU_UNREAD_FIRST = 2
            MU_UNREAD_SECOND = 3


            def read_value():
                return LAMBDA_READ_VALUE
        """},
       "from fx import order\nreturn order.read_value()",
       Edit("fx/order.py", "MU_UNREAD_FIRST = 2\nMU_UNREAD_SECOND = 3", "MU_UNREAD_SECOND = 3\nMU_UNREAD_FIRST = 2"),
       "none"),
    _c("P5-unread-def-signature", "precision", "design step 2: an unread def's statement (its default) is not selected",
       {"fx/sigs.py": """
            NU_SCALE_READ = 2
            NU_SCALE_UNREAD = 5


            def scaled(x, scale=NU_SCALE_READ):
                return x * scale


            def never_called(x, scale=NU_SCALE_UNREAD):
                return x * scale
        """},
       "from fx import sigs\nreturn sigs.scaled(3)",
       Edit("fx/sigs.py", "NU_SCALE_UNREAD = 5", "NU_SCALE_UNREAD = 6"),
       "none"),

    # ── recall: one row per shape of the design's census and the brief's list ──────────────────────────────────
    _c("R2-bare-expression", "recall", "census: a bare expression statement (effect, module grain)",
       {"fx/table.py": """
            XI_TABLE = {}
        """,
        "fx/fill.py": """
            from fx.table import XI_TABLE

            XI_TABLE.update(xi_entry=3)
        """},
       "from fx import fill, table\nreturn table.XI_TABLE['xi_entry']",
       Edit("fx/fill.py", "XI_TABLE.update(xi_entry=3)", "XI_TABLE.update(xi_entry=4)"),
       "the value read turns 4", "3", "4"),
    _c("R3-subscript-store", "recall", "census: a store into a subscript (effect, module grain)",
       {"fx/table.py": """
            XI_TABLE = {}
        """,
        "fx/fill.py": """
            from fx import table

            table.XI_TABLE["xi_entry"] = 7
        """},
       "from fx import fill, table\nreturn table.XI_TABLE['xi_entry']",
       Edit("fx/fill.py", 'table.XI_TABLE["xi_entry"] = 7', 'table.XI_TABLE["xi_entry"] = 8'),
       "the value read turns 8", "7", "8"),
    _c("R4-attribute-store", "recall", "census: a store into an attribute of another module (effect, module grain)",
       {"fx/limits.py": """
            OMICRON_LIMIT = 1
        """,
        "fx/override.py": """
            from fx import limits

            limits.OMICRON_LIMIT = 9
        """},
       "from fx import override, limits\nreturn limits.OMICRON_LIMIT",
       Edit("fx/override.py", "limits.OMICRON_LIMIT = 9", "limits.OMICRON_LIMIT = 10"),
       "the value read turns 10", "9", "10"),
    _c("R5-call-decorator-order", "recall",
       "census: a first-party registering decorator, call form (effect); the edit reorders two registrations",
       {"fx/registry.py": """
            PI_HANDLERS = {}


            def register_as(key):
                def install(function):
                    PI_HANDLERS[key] = function
                    return function
                return install
        """,
        "fx/handlers.py": """
            from fx.registry import register_as


            @register_as("pi_slot")
            def pi_first():
                return 1


            @register_as("pi_slot")
            def pi_second():
                return 2
        """},
       "from fx import handlers, registry\nreturn registry.PI_HANDLERS['pi_slot']()",
       Edit("fx/handlers.py",
            '@register_as("pi_slot")\ndef pi_first():\n    return 1\n\n\n@register_as("pi_slot")\ndef pi_second():\n    return 2',
            '@register_as("pi_slot")\ndef pi_second():\n    return 2\n\n\n@register_as("pi_slot")\ndef pi_first():\n    return 1'),
       "the later registration wins: 2 turns 1 (pi_second's DefPin is unchanged)", "2", "1"),
    _c("R6-bare-decorator-order", "recall",
       "design 'decorators included': a BARE first-party decorator (the census predicate does not count it)",
       {"fx/registry.py": """
            RHO_HANDLERS = {}


            def rho_register(function):
                RHO_HANDLERS[function.__name__.split("_")[0]] = function
                return function
        """,
        "fx/handlers.py": """
            from fx.registry import rho_register


            @rho_register
            def rhoslot_first():
                return 1


            @rho_register
            def rhoslot_second():
                return 2
        """},
       "from fx import handlers, registry\nreturn registry.RHO_HANDLERS['rhoslot']()",
       Edit("fx/handlers.py",
            "@rho_register\ndef rhoslot_first():\n    return 1\n\n\n@rho_register\ndef rhoslot_second():\n    return 2",
            "@rho_register\ndef rhoslot_second():\n    return 2\n\n\n@rho_register\ndef rhoslot_first():\n    return 1"),
       "the later registration wins: 2 turns 1", "2", "1"),
    _c("R7-import-time-call-mutates-global", "recall",
       "census: a binding statement whose evaluation calls first-party code that mutates another global (effect)",
       {"fx/state.py": """
            SIGMA_LEVEL = 0


            def configure_sigma(level):
                global SIGMA_LEVEL
                SIGMA_LEVEL = level
                return level
        """,
        "fx/setup.py": """
            from fx import state

            SIGMA_HANDLE = state.configure_sigma(5)
        """},
       "from fx import setup, state\nreturn state.SIGMA_LEVEL",
       Edit("fx/setup.py", "state.configure_sigma(5)", "state.configure_sigma(6)"),
       "the level read turns 6", "5", "6"),
    _c("R8-derived-constant-chain", "recall", "design step 2: E' gains R of each selected statement",
       {"fx/derived.py": """
            KAPPA_BASE_RATE = 3
            KAPPA_DOUBLED_RATE = KAPPA_BASE_RATE * 2


            def doubled_rate():
                return KAPPA_DOUBLED_RATE
        """},
       "from fx import derived\nreturn derived.doubled_rate()",
       Edit("fx/derived.py", "KAPPA_BASE_RATE = 3", "KAPPA_BASE_RATE = 4"),
       "6 turns 8", "6", "8"),
    _c("R9-default-argument-constant", "recall", "design step 2: a def statement's defaults are in its R",
       {"fx/sigs.py": """
            NU_SCALE_READ = 2


            def scaled(x, scale=NU_SCALE_READ):
                return x * scale
        """},
       "from fx import sigs\nreturn sigs.scaled(3)",
       Edit("fx/sigs.py", "NU_SCALE_READ = 2", "NU_SCALE_READ = 3"),
       "6 turns 9 (the def's own digest names NU_SCALE_READ, not its value)", "6", "9"),
    _c("R10-reexport-renamed", "recall", "design step 2: a selected ImportFrom puts its source names in E'",
       {"fx/k.py": """
            TAU_SOURCE_LIMIT = 4
        """,
        "fx/facade.py": """
            from fx.k import TAU_SOURCE_LIMIT as TAU_PUBLIC_LIMIT
        """},
       "from fx import facade\nreturn facade.TAU_PUBLIC_LIMIT",
       Edit("fx/k.py", "TAU_SOURCE_LIMIT = 4", "TAU_SOURCE_LIMIT = 5"),
       "4 turns 5", "4", "5"),
    _c("R11-class-attribute-via-factory", "recall",
       "design step 1: an attribute read on an instance names the attribute; the run never names the class",
       {"fx/widgets.py": """
            class Widget:
                THETA_READ_LIMIT = 3
                THETA_UNREAD_LIMIT = 9
        """,
        "fx/factory.py": """
            from fx.widgets import Widget


            def build_widget():
                return Widget()
        """},
       "from fx import factory\nreturn factory.build_widget().THETA_READ_LIMIT",
       Edit("fx/widgets.py", "THETA_READ_LIMIT = 3", "THETA_READ_LIMIT = 4"),
       "3 turns 4", "3", "4"),
    _c("R12-dunder-hash", "recall", "census: a class dunder read by the interpreter (__hash__ = None)",
       {"fx/keys.py": """
            class UpsilonKey:
                __hash__ = None


            def make_upsilon():
                return UpsilonKey()
        """},
       "from fx import keys\ntry:\n    hash(keys.make_upsilon())\nexcept TypeError:\n    return 'unhashable'\nreturn 'hashable'",
       Edit("fx/keys.py", "__hash__ = None", "__hash__ = object.__hash__"),
       "'unhashable' turns 'hashable'", "'unhashable'", "'hashable'"),
    _c("R13-dunder-match-args", "recall", "census: a class dunder read by the interpreter (__match_args__)",
       {"fx/points.py": """
            class PhiPoint:
                __match_args__ = ("phi_x", "phi_y")

                def __init__(self):
                    self.phi_x = 1
                    self.phi_y = 2


            def make_phi():
                return PhiPoint()
        """},
       "from fx import points\nmatch points.make_phi():\n    case points.PhiPoint(first, second):\n        return (first, second)\nreturn None",
       Edit("fx/points.py", '__match_args__ = ("phi_x", "phi_y")', '__match_args__ = ("phi_y", "phi_x")'),
       "(1, 2) turns (2, 1)", "(1, 2)", "(2, 1)"),
    _c("R14-dunder-array-priority", "recall",
       "census: a class dunder read by a C library (__array_priority__), on an instance the run never names",
       {"fx/arrays.py": """
            class ChiOperand:
                __array_priority__ = 20.0

                def __radd__(self, other):
                    return "chi-radd"


            CHI_INSTANCE = ChiOperand()
        """},
       "import numpy as np\nfrom fx import arrays\nreturn np.ones(2) + arrays.CHI_INSTANCE",
       Edit("fx/arrays.py", "__array_priority__ = 20.0", "__array_priority__ = -20.0"),
       "numpy defers to __radd__ only while the priority is higher; the result changes type", None, None),
    _c("R15-new-shadowing-def", "recall",
       "design step 3: a new statement binding a stored name (Python binds the LAST def; the DefPin pins ordinal 0)",
       {"fx/shadow.py": """
            def omega_value():
                return 1
        """},
       "from fx import shadow\nreturn shadow.omega_value()",
       Edit("fx/shadow.py", "def omega_value():\n    return 1\n",
            "def omega_value():\n    return 1\n\n\ndef omega_value():\n    return 2\n"),
       "1 turns 2", "1", "2"),
    _c("R16-new-shadowing-rebinding", "recall", "design step 3: a new top-level rebinding of a stored name",
       {"fx/shadow.py": """
            def omega_value():
                return 1
        """},
       "from fx import shadow\nreturn shadow.omega_value()",
       Edit("fx/shadow.py", "def omega_value():\n    return 1\n",
            "def omega_value():\n    return 1\n\n\nomega_value = lambda: 2\n"),
       "1 turns 2", "1", "2"),
    _c("R17-removed-statement", "recall", "design step 3: a removed selected statement",
       {"fx/rebind.py": """
            ETA_VALUE = 1
            ETA_VALUE = 2
        """},
       "from fx import rebind\nreturn rebind.ETA_VALUE",
       Edit("fx/rebind.py", "ETA_VALUE = 1\nETA_VALUE = 2\n", "ETA_VALUE = 1\n"),
       "2 turns 1", "2", "1"),
    _c("R18-reordered-statement", "recall", "design step 3: a reordering of selected statements",
       {"fx/rebind.py": """
            ETA_VALUE = 1
            ETA_VALUE = 2
        """},
       "from fx import rebind\nreturn rebind.ETA_VALUE",
       Edit("fx/rebind.py", "ETA_VALUE = 1\nETA_VALUE = 2\n", "ETA_VALUE = 2\nETA_VALUE = 1\n"),
       "2 turns 1", "2", "1"),
    _c("R19-lambda-in-table-constant", "recall", "design step 2: a lambda body's reads are in its statement's R",
       {"fx/ops.py": """
            PSI_FACTOR = 3
            PSI_OPS = {"triple": lambda x: x * PSI_FACTOR}
        """},
       "from fx import ops\nreturn ops.PSI_OPS['triple'](2)",
       Edit("fx/ops.py", "PSI_FACTOR = 3", "PSI_FACTOR = 4"),
       "6 turns 8", "6", "8"),
    _c("R20-lambda-in-table-body", "recall", "design step 1: a module-level lambda that ran selects its statement by line",
       {"fx/ops.py": """
            PSI_FACTOR = 3
            PSI_OPS = {"triple": lambda x: x * PSI_FACTOR}
        """},
       "from fx import ops\nreturn ops.PSI_OPS['triple'](2)",
       Edit("fx/ops.py", "lambda x: x * PSI_FACTOR", "lambda x: x * PSI_FACTOR + 1"),
       "6 turns 7", "6", "7"),
    _c("R21-match-keyword-pattern", "recall", "design step 1: a match class pattern's keywords are in co_consts",
       {"fx/shapes.py": """
            class Disc:
                zeta_radius = 5


            def make_disc():
                return Disc()
        """},
       "from fx import shapes\nmatch shapes.make_disc():\n    case shapes.Disc(zeta_radius=r):\n        return r\nreturn None",
       Edit("fx/shapes.py", "zeta_radius = 5", "zeta_radius = 6"),
       "5 turns 6", "5", "6"),
    _c("R22-getattr-literal", "recall", "design step 1: getattr(x, 'literal') puts the literal in co_consts",
       {"fx/options.py": """
            class Options:
                gamma_option = 1
                gamma_unread_option = 2


            def make_options():
                return Options()
        """},
       "from fx import options\nreturn getattr(options.make_options(), 'gamma_option')",
       Edit("fx/options.py", "gamma_option = 1", "gamma_option = 3"),
       "1 turns 3", "1", "3"),

    # ── findings: recall rows the design as written misses ──────────────────────────────────────────────────────
    _c("F1-subclass-override-in-another-module", "recall",  # a finding against the first design; closed by R4
      
       "design step 3: names are restricted to identifiers in the module's source at capture; an override added to "
       "a subclass whose module never spelled the method shadows the inherited one",
       {"fx/base.py": """
            class EpsilonBase:
                def epsilon_method(self):
                    return 1
        """,
        "fx/sub.py": """
            from fx.base import EpsilonBase


            class EpsilonSub(EpsilonBase):
                pass


            def make_sub():
                return EpsilonSub()
        """},
       "from fx import sub\nreturn sub.make_sub().epsilon_method()",
       Edit("fx/sub.py", "class EpsilonSub(EpsilonBase):\n    pass\n",
            "class EpsilonSub(EpsilonBase):\n    def epsilon_method(self):\n        return 2\n"),
       "1 turns 2", "1", "2"),
    _c("F2-new-dataclass-field-iterated", "recall",  # a finding against the first design; closed by R5
      
       "design step 3: a new field read by ITERATION (dataclasses.fields, replace, __eq__, __repr__) is bound by "
       "no stored name; the decorator reads the whole class body",
       {"fx/config.py": """
            import dataclasses


            @dataclasses.dataclass(frozen=True)
            class IotaConfig:
                iota_first: int = 1


            def make_config():
                return IotaConfig()
        """},
       "import dataclasses\nfrom fx import config\nreturn len(dataclasses.fields(config.make_config()))",
       Edit("fx/config.py", "    iota_first: int = 1\n", "    iota_first: int = 1\n    iota_second: int = 2\n"),
       "1 turns 2", "1", "2"),
    _c("F3-new-enum-member-iterated", "recall",  # a finding against the first design; closed by R5
      
       "design step 3: a new Enum member read by iteration; the metaclass reads the whole body",
       {"fx/phases.py": """
            import enum


            class LambdaPhase(enum.Enum):
                LAMBDA_SOLID = 1
                LAMBDA_LIQUID = 2
        """},
       "from fx import phases\nreturn len(list(phases.LambdaPhase))",
       Edit("fx/phases.py", "    LAMBDA_LIQUID = 2\n", "    LAMBDA_LIQUID = 2\n    LAMBDA_GAS = 3\n"),
       "2 turns 3", "2", "3"),
    _c("F4-dunder-default", "recall",  # a finding against the first design; closed by R1
      
       "design step 2: a def that RAN without its name in E' (a dunder reached through an operator) does not select "
       "its own statement, so its defaults are not followed; the same shape on __init__ is rescued only because the "
       "import machinery's own code names __init__ ([M] the prototype)",
       {"fx/boxes.py": """
            MU_NEG_SCALE = 2


            class MuBox:
                def __neg__(self, scale=MU_NEG_SCALE):
                    return -scale


            def make_box():
                return MuBox()
        """},
       "from fx import boxes\nreturn -boxes.make_box()",
       Edit("fx/boxes.py", "MU_NEG_SCALE = 2", "MU_NEG_SCALE = 3"),
       "-2 turns -3", "-2", "-3"),
    _c("F5-init-subclass-registration", "recall",  # a finding against the first design; closed by R2D
      
       "census: an import-time call into first-party code with no Call node (a base's __init_subclass__ fired by a "
       "class statement's keywords; 8 sites in orpheus/, the realizer tables among them)",
       {"fx/plug.py": """
            class OmegaPluginBase:
                omega_catalog = {}

                def __init_subclass__(cls, *, omega_slot, **kwargs):
                    super().__init_subclass__(**kwargs)
                    OmegaPluginBase.omega_catalog[omega_slot] = cls
        """,
        "fx/plugins.py": """
            from fx.plug import OmegaPluginBase


            class OmegaFirst(OmegaPluginBase, omega_slot="omega"):
                pass


            class OmegaSecond(OmegaPluginBase, omega_slot="omega"):
                pass
        """},
       "from fx import plug, plugins\nreturn plug.OmegaPluginBase.omega_catalog['omega'].__name__",
       Edit("fx/plugins.py",
            'class OmegaFirst(OmegaPluginBase, omega_slot="omega"):\n    pass\n\n\nclass OmegaSecond(OmegaPluginBase, omega_slot="omega"):\n    pass\n',
            'class OmegaSecond(OmegaPluginBase, omega_slot="omega"):\n    pass\n\n\nclass OmegaFirst(OmegaPluginBase, omega_slot="omega"):\n    pass\n'),
       "'OmegaSecond' turns 'OmegaFirst'", "'OmegaSecond'", "'OmegaFirst'"),
    _c("F6-annotation-read-at-run-time", "recall",  # a finding against the first design; closed by R5
      
       "design step 1: a signature scope (a class's __annotate__, Python 3.14) that RAN after import reads every "
       "annotated statement of its scope; it is excluded from E and selects nothing",
       {"fx/records.py": """
            BetaAlias = int


            class BetaRecord:
                beta_field: BetaAlias


            def record_type():
                return BetaRecord
        """},
       "import typing\nfrom fx import records\nreturn [t.__name__ for t in typing.get_type_hints(records.record_type()).values()]",
       Edit("fx/records.py", "BetaAlias = int", "BetaAlias = float"),
       "['int'] turns ['float']", "['int']", "['float']"),

    _c("F7-def-found-by-iteration", "recall",  # carried by ruling R1 alone ([M] the ablation)
       "ruling R1: a def found by namespace iteration (as pytest finds a test) is never named, so only R1 selects its "
       "statement and follows its default",
       {"fx/runner.py": """
            OMEGA_RUN_SCALE = 2


            class OmegaSuite:
                def omega_run_first(self, scale=OMEGA_RUN_SCALE):
                    return scale
        """},
       "from fx import runner\nsuite = runner.OmegaSuite()\n"
       "return [getattr(suite, n)() for n in dir(suite) if n.startswith('omega_run_')]",
       Edit("fx/runner.py", "OMEGA_RUN_SCALE = 2", "OMEGA_RUN_SCALE = 3"),
       "[2] turns [3]", "[2]", "[3]"),
    _c("N4-frozenset-of-module-names", "recall",  # carried by ruling R6 alone
       "ruling R6: names in a frozenset constant (a set literal of strings), read from a module namespace",
       {"fx/strs.py": """
            qa_rho_alpha = 1
            qa_rho_beta = 2
        """},
       "from fx import strs\nreturn sorted(getattr(strs, n) for n in {'qa_rho_alpha', 'qa_rho_beta'})",
       Edit("fx/strs.py", "qa_rho_beta = 2", "qa_rho_beta = 5"),
       "[1, 2] turns [1, 5]", "[1, 2]", "[1, 5]"),
    _c("N5-dotted-attrgetter-to-a-module", "recall",  # carried by ruling R6 alone
       "ruling R6: a dotted string ('module.name') split into identifier tokens",
       {"fx/inner_mod.py": """
            zeta_dotted_val = 1
        """},
       "import operator\nimport fx.inner_mod\nimport fx\nreturn operator.attrgetter('inner_mod.zeta_dotted_val')(fx)",
       Edit("fx/inner_mod.py", "zeta_dotted_val = 1", "zeta_dotted_val = 2"),
       "1 turns 2", "1", "2"),

    _c("T2-map-of-a-first-party-callable", "recall",  # qa's build review T2
       "_import_time._calls_untrusted: a pure higher-order callee (map, filter, sorted, min, max) handed a first-party "
       "callable runs it at import; the statement is an effect",
       {"fx/reg.py": """
            UPSILON_REGISTRY = []


            def upsilon_register(name):
                UPSILON_REGISTRY.append(name)
        """,
        "fx/plugs.py": """
            from fx.reg import upsilon_register

            UPSILON_PLUGINS = ("a",)
            _UPSILON_DONE = list(map(upsilon_register, UPSILON_PLUGINS))
        """},
       "from fx import plugs, reg\nreturn list(reg.UPSILON_REGISTRY)",
       Edit("fx/plugs.py", 'UPSILON_PLUGINS = ("a",)', 'UPSILON_PLUGINS = ("a", "b")'),
       "['a'] turns ['a', 'b']", "['a']", "['a', 'b']"),

    # ── declared boundaries (design item 4): the documented under-selection ─────────────────────────────────────
    _c("U1-computed-getattr-name", "under",
       "design item 4: a name built from data, never present whole as a string constant in executed code (a module "
       "attribute, so no selected class header carries it: ruling R5 closes the class-attribute form)",
       {"fx/dyn.py": """
            delta_opt = 1
            DELTA_OTHER = 0
        """},
       "from fx import dyn\nname = '_'.join(['delta', 'opt'])\nreturn getattr(dyn, name)",
       Edit("fx/dyn.py", "delta_opt = 1", "delta_opt = 2"),
       "1 turns 2", "1", "2"),
    _c("U2-guarded-probe-of-a-new-name", "recall",  # boundary 2 dissolved by R4
      
       "design item 4: a new top-level name read only through a guarded probe and absent from the source at capture",
       {"fx/probe.py": """
            ZETA_PRESENT = 1
        """},
       "from fx import probe\nreturn getattr(probe, 'zeta_late_name', None)",
       Edit("fx/probe.py", "ZETA_PRESENT = 1\n", "ZETA_PRESENT = 1\nzeta_late_name = 5\n"),
       "None turns 5", "None", "5"),
    _c("U3-module-namespace-iteration", "under",
       "design item 4: vars() / __dict__ iteration over a module namespace",
       {"fx/consts2.py": """
            IOTA_A = 1
            IOTA_B = 2
        """},
       "from fx import consts2\nreturn sorted(k for k in vars(consts2) if k.startswith('IOTA_'))",
       Edit("fx/consts2.py", "IOTA_B = 2\n", "IOTA_B = 2\nIOTA_C = 3\n"),
       "two names turn three", "['IOTA_A', 'IOTA_B']", "['IOTA_A', 'IOTA_B', 'IOTA_C']"),
    _c("U4-c-code-reads-a-non-dunder-attribute", "recall",  # closed by R2D (a class with a base)
      
       "design item 4: C code reading a non-dunder class attribute (ctypes' _fields_)",
       {"fx/structs.py": """
            import ctypes


            class KappaPair(ctypes.Structure):
                _fields_ = [("kappa_a", ctypes.c_int32), ("kappa_b", ctypes.c_int32)]
        """},
       "import ctypes\nfrom fx import structs\nreturn ctypes.sizeof(structs.KappaPair)",
       Edit("fx/structs.py", '("kappa_b", ctypes.c_int32)', '("kappa_b", ctypes.c_int64)'),
       "8 turns 16", "8", "16"),

    # ── declared over-selection (stays green before and after) ─────────────────────────────────────────────────
    _c("O1-name-collision-across-modules", "over",
       "design step 2: names are not module-qualified; a statement of another pinned module binding a read name is selected",
       {"fx/one.py": """
            NU_SHARED_NAME = 1


            def read_one():
                return NU_SHARED_NAME
        """,
        "fx/two.py": """
            NU_SHARED_NAME = 100


            def two_marker():
                return 0
        """},
       "from fx import one, two\ntwo.two_marker()\nreturn one.read_one()",
       Edit("fx/two.py", "NU_SHARED_NAME = 100", "NU_SHARED_NAME = 101"),
       "none: the run reads one.NU_SHARED_NAME only"),
    _c("O2-statement-is-the-unit", "over",
       "design step 2: an ImportFrom binding two names is one statement; reading one follows both sources",
       {"fx/k.py": """
            XI_ONE = 1
            XI_TWO = 2
        """,
        "fx/facade.py": """
            from fx.k import XI_ONE, XI_TWO
        """},
       "from fx import facade\nreturn facade.XI_ONE",
       Edit("fx/k.py", "XI_TWO = 2", "XI_TWO = 3"),
       "none: the run reads XI_ONE only"),
    _c("O3-effect-statement-module-grain", "over",
       "design step 2: an effect statement is selected for every run that imported its module",
       {"fx/table.py": """
            XI_TABLE = {}
        """,
        "fx/fill.py": """
            from fx.table import XI_TABLE

            XI_TABLE.update(xi_unread_entry=3)


            def fill_marker():
                return 0
        """},
       "from fx import fill\nreturn fill.fill_marker()",
       Edit("fx/fill.py", "XI_TABLE.update(xi_unread_entry=3)", "XI_TABLE.update(xi_unread_entry=4)"),
       "none: the run never reads the table"),
)



# ── the import outcome (ruling R3 amended): an unselected statement that changes what importing the set does ────
OUTCOME_CASES: tuple[Case, ...] = (
    _c("I1-new-import-time-warning", "outcome",
       "ruling R3: a warning emitted while importing is part of the recorded outcome; the answer does not change",
       {"fx/quiet.py": """
            import warnings

            CHI_READ_VALUE = 1
            CHI_UNREAD_VALUE = 2
        """},
       "from fx import quiet\nreturn quiet.CHI_READ_VALUE",
       # A pure binding that warns at COMPILE time (an invalid escape, SyntaxWarning): it binds nothing new and calls
       # nothing, so no effect rule (R2D) selects it; only the outcome check (R3) sees it ([M] ta/ablate: the row goes
       # valid with R3 dropped). A warn() call would be an effect under R2D and would not witness R3.
       Edit("fx/quiet.py", "CHI_UNREAD_VALUE = 2", 'CHI_UNREAD_VALUE = "\\d"'),
       "none on the answer; compiling the module now emits a SyntaxWarning", "1", "1"),
    _c("I2-unselected-statement-starts-raising", "outcome",
       "ruling R3: an unselected statement that raises at import (a dead table naming a removed def)",
       {"fx/legacy.py": """
            def chi_old():
                return 1


            CHI_TOL = 1e-6
            _CHI_LEGACY = {"old": chi_old}
        """},
       "from fx import legacy\nreturn legacy.CHI_TOL",
       Edit("fx/legacy.py", "def chi_old():\n    return 1\n", ""),
       "the import raises NameError", "1e-06", "raised NameError"),
)


def _qa_case(row: dict, kind: str) -> Case:
    files = {k: v for k, v in row["files"].items()}
    run = row["run"].rstrip("\n") + "\nreturn RESULT"
    return Case(f"QA-{row['id']}", kind, files, run, Edit(row["file"], None, row["new"]),
                "measured by the gate (qa's row asserts the answer changes)", f"qa_b2_design.md row {row['id']}: {row['title']}")


def _qa_cases() -> tuple[Case, ...]:
    from ._statement_pin_qa_rows import MEMO_ROWS

    return tuple(_qa_case(row, "precision" if row["id"] == "PC2" else "recall") for row in MEMO_ROWS)


QA_CASES: tuple[Case, ...] = _qa_cases()
ALL_CASES: tuple[Case, ...] = CASES + OUTCOME_CASES + QA_CASES
BY_ID = {case.id: case for case in ALL_CASES}
