"""``tests.data`` — and the one address of its frozen artefacts."""
from pathlib import Path

#: The directory of frozen baselines this tree's gates compare against
#: (captured on a named commit). ONE spelling: a gate that builds the path
#: itself — cwd-relative, or from ``__file__`` — is a second address for the
#: same directory (#444; mirrors ``tests.numerics.NUMERICS_DATA``).
DATA_TREE_DATA = Path(__file__).resolve().parent / "data"
