"""The harnesses this repository generates for. Adding one is a module here
implementing ``base.Harness`` and an entry in ``HARNESSES``; nothing outside
``targets/`` changes."""
from .claude_code import ClaudeCode

HARNESSES = (ClaudeCode(),)
