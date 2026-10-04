"""The first code a traced memo's generating process runs (#405 P3); run by path with ``-P``.

The tracer starts before anything else is imported, so every first-party module body and def the generation
runs is recorded, the ones its own imports run included; the audit hook records every file opened for
reading. Only then are the job and :mod:`orpheus.numerics.traced_memo` loaded. The process answers on stdout
with one pickle (``{"ok": True}`` or the exception) and exits.
"""
import sys


def main() -> None:
    """Start the tracer and the audit hook, then load and run the job (everything after the first line is traced)."""
    monitoring = sys.monitoring
    tool = next(t for t in (monitoring.PROFILER_ID, 3, 4, monitoring.OPTIMIZER_ID) if monitoring.get_tool(t) is None)
    monitoring.use_tool_id(tool, "traced_memo")
    codes: set = set()

    def started(code, _offset):
        codes.add(code)
        return monitoring.DISABLE

    monitoring.register_callback(tool, monitoring.events.PY_START, started)
    monitoring.set_events(tool, monitoring.events.PY_START)
    opened: list[str] = []
    recording = [True]

    def audit(event, args):
        if recording[0] and event == "open" and args and isinstance(args[0], str):
            mode = args[1] if len(args) > 1 and isinstance(args[1], str) else "r"
            if not any(c in mode for c in "wax+"):
                opened.append(args[0])

    sys.addaudithook(audit)

    def stop_recording():
        monitoring.set_events(tool, 0)
        recording[0] = False

    import pickle  # after the tracer, so its own imports are recorded

    envelope = pickle.loads(sys.stdin.buffer.read())
    sys.path[:] = envelope["sys_path"]
    job = pickle.loads(envelope["job"])
    from orpheus.numerics import traced_memo

    result = traced_memo._generate_here(job, codes, stop_recording, opened)
    sys.stdout.buffer.write(pickle.dumps(result))
    sys.stdout.flush()


if __name__ == "__main__":  # an import (a roster walk over the package) runs nothing
    main()
