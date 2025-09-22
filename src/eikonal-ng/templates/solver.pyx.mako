#
# @Author : Jean-Pascal Mercier <jean-pascal.mercier@agsis.com>
#
# @Copyright (C) 2010 Jean-Pascal Mercier
#
# All rights reserved.
#
# This file contains the python interface for the eikonal
# solver using numpy as array container.
#

__doc__ = """
Python interface to low level Eikonal Solver.

This module contains Seeded as well as Initial Value eikonal solver interface
to low level C++ Fast Marching Solver.
"""

cimport numpy as cnp
import numpy as np
from libcpp cimport bool

#
# I am not sure this is necessary since cython docs are all but precise on that
# particular matter. It really worked both ways. Decided to kept it for the
# sake of safety.
#
cdef extern from "numpy/arrayobject.h":
    void import_array() except *
import_array()


% for i in [2, 3]:
def SFMM${i}(cnp.ndarray[double, ndim = 2, mode = 'c'] seeds not None,
             cnp.ndarray[double, ndim = 1, mode = 'c'] stime not None,
             cnp.ndarray[int,    ndim = ${i}, mode = 'c'] tag not None, # VERIFY THE TYPE OF AN ENUM
             cnp.ndarray[double, ndim = ${i}, mode = 'c'] viscosity not None,
             cnp.ndarray[double, ndim = ${i}, mode = 'c'] arrival not None,
             float spacing,
             bool second_order = True):
    """
    Low level Seeded Fast 2nd Order Fast Marching Method for ${i}D Grid.

    :param seeds: an ndarray
    :param stime: an ndarra
    :param tag: an ndarray
    :param viscosity: an ndarray
    :param arrival: an ndarray
    :param spacing: the spacing of the grid
    :param second_order: a boolean
    """

    cdef ArrayDescriptor${i}_double cviscosity
    cviscosity.init(<size_t *>viscosity.shape, <double *>viscosity.data)

    cdef ArrayDescriptor${i}_tag ctag
    ctag.init(<size_t *>tag.shape, <MarchingTag *>tag.data)

    cdef OrthogonalGrid${i}_double carrival
    carrival.init(<size_t *>arrival.shape, <double *>arrival.data, spacing)

    cdef ArrayDescriptor_doublev${i} cseeds
    cseeds.init(<size_t *>seeds.shape, <doublev${i} *>seeds.data)

    with nogil:
        FMM${i}(cseeds, ctag, cviscosity, carrival, second_order)
    return arrival

def I_FMM${i}(cnp.ndarray[int,    ndim = ${i}, mode = 'c'] tag not None, # VERIFY THE TYPE OF AN ENUM
             cnp.ndarray[double, ndim = ${i}, mode = 'c'] viscosity not None,
             cnp.ndarray[double, ndim = ${i}, mode = 'c'] arrival not None,
             float spacing,
             bool second_order = True):
    """
    Low level Seeded Fast 2nd Order Fast Marching Method for ${i}D Grid.

    Parameters :    arrival - This is merely a buffer for the arrival grid.

    """

    cdef ArrayDescriptor${i}_double cviscosity
    cviscosity.init(<size_t *>viscosity.shape, <double *>viscosity.data)

    cdef ArrayDescriptor${i}_tag ctag
    ctag.init(<size_t *>tag.shape, <MarchingTag *>tag.data)

    cdef OrthogonalGrid${i}_double carrival
    carrival.init(<size_t *>arrival.shape, <double *>arrival.data, spacing)

    with nogil:
        IFMM${i}(ctag, cviscosity, carrival, second_order)
    return arrival

% endfor

def SFMM(cnp.ndarray seeds,
         cnp.ndarray time,
         cnp.ndarray tag,
         cnp.ndarray viscosity,
         cnp.ndarray arrival,
         float spacing = 1,
         **kw):
    """
    This is a redirection.
    """

    for i in range(arrival.ndim):
        if (tag.shape[i] != viscosity.shape[i]) or (tag.shape[i] != arrival.shape[i]):
            raise AttributeError("Array shape Mismatch")
    if seeds.shape[1] != arrival.ndim:
        raise AttributeError("Seeds array dimension mismatch")
% for i in [2, 3]:
    if viscosity.ndim == ${i}:
        return SFMM${i}(seeds, time, tag, viscosity, arrival, spacing, **kw)
% endfor

def I_FMM(cnp.ndarray tag,
         cnp.ndarray viscosity,
         cnp.ndarray arrival,
         float spacing = 1,
         **kw):
    """
    This is a redirection.
    """

    for i in range(arrival.ndim):
        if (tag.shape[i] != viscosity.shape[i]) or (tag.shape[i] != arrival.shape[i]):
            raise AttributeError("Array shape Mismatch")
% for i in [2, 3]:
    if viscosity.ndim == ${i}:
        return I_FMM${i}(tag, viscosity, arrival, spacing, **kw)
% endfor

# vim: filetype=pyrex
#
#
