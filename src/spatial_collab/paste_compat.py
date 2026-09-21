"""Signature-only bridge for PASTE 1.4.0 and POT 0.9.6.post1.

POT >=0.9.5 passes the regularizer gradient to its line-search callback.
PASTE's five-argument callback does not use that optional cached gradient.
Keep the official objective and solver, bridge only this calling convention.
The worker is serial and the original function is restored even on failure.
"""

from contextlib import contextmanager
from functools import wraps
from importlib.metadata import version
import inspect


@contextmanager
def line_search_compatibility():
    import ot

    active = version("paste-bio") == "1.4.0" and version("POT") == "0.9.6.post1"
    original = ot.optim.cg
    if active:

        @wraps(original)
        def compatible(*args, **kwargs):
            bound = inspect.signature(original).bind(*args, **kwargs)
            callback = bound.arguments.get("line_search")
            if callback is not None:

                def adapter(
                    cost, coupling, direction, gradient, cost_value, regularizer_gradient=None, **options
                ):
                    return callback(cost, coupling, direction, gradient, cost_value, **options)

                bound.arguments["line_search"] = adapter
            return original(*bound.args, **bound.kwargs)

        ot.optim.cg = compatible
    try:
        yield {
            "active": active,
            "kind": "line_search_call_signature_only",
            "paste-bio": version("paste-bio"),
            "POT": version("POT"),
            "objective_changed": False,
        }
    finally:
        ot.optim.cg = original
