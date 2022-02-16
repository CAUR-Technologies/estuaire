#
# @Author : Jean-Pascal Mercier <jean-pascal.mercier@agsis.com>
#
# @Copyright (C) 2010 Jean-Pascal Mercier
#
# All rights reserved.
#
__doc__ = """
"""

import inspect
import argparse
import pickle

import os
import sys


def main(fct, **descr):
    """
    """
    parser = argparse.ArgumentParser(description = fct.__doc__)
    fct_descr = inspect.getargspec(fct)
    if "output" not in fct_descr.args:
        parser.add_argument("--output")
    parser.add_argument("--info")

    ndefaults = -len(fct_descr.defaults) if fct_descr.defaults is not None else None
    for arg in fct_descr.args[:ndefaults]:
        arg_descr = dict()
        if arg in descr:
            arg_descr['type'] = descr[arg]
        parser.add_argument("--" + arg, required = True, **arg_descr)

    if ndefaults is not None:
        for arg, value in zip(fct_descr.args[ndefaults:], fct_descr.defaults):
            arg_descr = dict()
            if arg in descr:
                arg_descr['type'] = descr[arg]
            parser.add_argument("--" + arg, required = False, default = value, **arg_descr)

    ns = parser.parse_args()
    argv = [getattr(ns, a) for a in fct_descr[0]]

    result = fct(*argv)

    info = " ".join(sys.argv)

    if (result != None) and (ns.output != None) and ("output" not in fct_descr.args):
        if os.path.isdir(ns.output):
            for k in result:
                fname = os.path.join(ns.output, k)
                pickle.dump(result[k], open(fname, 'w'), protocol = pickle.HIGHEST_PROTOCOL)
                info += '\n     --> %s' % fname
        else:
            dirname = os.path.dirname(ns.output)
            pickle.dump(result, open(ns.output, 'w'), protocol = pickle.HIGHEST_PROTOCOL)
            info += '\n     --> %s' % ns.output

    if ns.info != None:
        with open(ns.info, 'w') as f:
            f.write(info)
