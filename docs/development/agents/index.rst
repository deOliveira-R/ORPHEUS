Agent definitions
=================

One page per project agent: its role (Key or Support), the workflow phases it
owns, the supports it may call, its method and its return contract. Each page
is generated, whole, into the agent's ``AGENT.md`` between ``GENERATED``
markers, so the harness's checks read every line of it; only the front matter
of an ``AGENT.md`` (tools, model, memory, the harness-specific header) is
maintained by hand.

.. toctree::
   :maxdepth: 1

   archivist
   cross-domain-attacker
   elegance-enforcer
   explorer
   literature-researcher
   method-implementer
   numerics-investigator
   qa
   test-architect
