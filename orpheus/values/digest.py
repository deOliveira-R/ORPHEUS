"""The one hash convention of ORPHEUS: blake2b with a 32-byte digest, as bytes and as hex.

Every digest of a byte string in ``orpheus/`` outside the closed ``derivations/`` reference code, and in ``tools/``, is
computed here. Not claimed: the closed ``orpheus/derivations/continuous/sood_registry/cache.py`` keeps its own sha256 of
a file, and the tests keep independent ``hashlib`` computations as witnesses (``tests/gates/values/test_content_identity.py``,
``tests/gates/memo/test_traced_memo_compiled_bytes.py``): a witness that called this module would agree with it
tautologically.

The module imports only the standard library: the traced memo's generating process loads it BY PATH before any
``orpheus`` import is possible (``orpheus/memo/_traced_memo_boot.py``), and registers that module as
``orpheus.values.digest`` once ``orpheus.values`` is imported, so the process holds one module object of this file."""
import hashlib
from typing import Final

#: The width, in bytes, of a digest: blake2b-256, wide enough that a collision between two contents is not a failure
#: mode a cache needs to consider.
DIGEST_SIZE: Final = 32
#: The width, in bytes, of the prefix of a digest that names a derived space (``orpheus.values.content.name_digest``).
#: A prefix of the 32-byte digest, not a blake2b of its own width: a name is a function of the content digest.
NAME_SIZE: Final = 8


def digest(data: bytes) -> bytes:
    """The digest of ``data``: exactly :data:`DIGEST_SIZE` bytes."""
    return hashlib.blake2b(data, digest_size=DIGEST_SIZE).digest()


def hexdigest(data: bytes) -> str:
    """The digest of ``data`` as lower-case hex (``2 * DIGEST_SIZE`` characters)."""
    return digest(data).hex()
