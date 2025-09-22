#
# @Author : Jean-Pascal Mercier <jean-pascal.mercier@agsis.com>
#
# @Copyright (C) 2010 Jean-Pascal Mercier
#
# All rights reserved.
#
__doc__ = """
"""

import argparse
import inspect
import os
import pickle
import sys
from typing import Any, Dict, Iterable, Sequence

import numpy as np


def _iter_required_optional(
    args: Sequence[str], defaults: Sequence[Any]
) -> tuple[Iterable[str], Iterable[tuple[str, Any]]]:
    if not defaults:
        return args, []
    split_index = len(args) - len(defaults)
    return args[:split_index], zip(args[split_index:], defaults)


def main(fct, **descr):
    """Turn a callable into a small argparse-powered CLI."""

    parser = argparse.ArgumentParser(description=fct.__doc__)
    spec = inspect.getfullargspec(fct)
    args = list(spec.args or [])
    defaults = list(spec.defaults or [])

    if "output" not in args:
        parser.add_argument("--output")
    parser.add_argument("--info")

    required_args, optional_args = _iter_required_optional(args, defaults)

    for arg in required_args:
        arg_descr: Dict[str, Any] = {}
        if arg in descr:
            arg_descr["type"] = descr[arg]
        parser.add_argument(f"--{arg}", required=True, **arg_descr)

    for arg, value in optional_args:
        arg_descr: Dict[str, Any] = {}
        if arg in descr:
            arg_descr["type"] = descr[arg]
        parser.add_argument(f"--{arg}", required=False, default=value, **arg_descr)

    ns = parser.parse_args()
    argv = [getattr(ns, a) for a in args]

    result = fct(*argv)

    info = " ".join(sys.argv)

    output_target = getattr(ns, "output", None)
    if result is not None and output_target is not None and "output" not in args:
        output_path = os.path.expanduser(output_target)
        if os.path.isdir(output_path) and isinstance(result, dict):
            for key, value in result.items():
                fname = os.path.join(output_path, key)
                with open(fname, "wb") as fdesc:
                    pickle.dump(value, fdesc, protocol=pickle.HIGHEST_PROTOCOL)
                info += f"\n     --> {fname}"
        else:
            os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
            with open(output_path, "wb") as fdesc:
                pickle.dump(result, fdesc, protocol=pickle.HIGHEST_PROTOCOL)
            info += f"\n     --> {output_path}"

    info_target = getattr(ns, "info", None)
    if info_target is not None:
        with open(info_target, "w", encoding="utf-8") as fdesc:
            fdesc.write(info)
