#
# @Author : Jean-Pascal Mercier <jean-pascal.mercier@agsis.com>
#
# @Copyright (C) 2010 Jean-Pascal Mercier
#
# All rights reserved.
#
#
# TODO : Eventually, the Low level C++ primitives can handle non contiguous
#        (strided) arrays in most of the primitives. We should take a look
#        where it can cause a problem and act consequently.
#
#        Until every fct are checked, contiguous array are enforced !!!!



__doc__ = """
Python Interface to low level RayTracing and Frechet Derivative
primitives.

This module contains raytracing, traveltime calculation,
sensivity(Frechet derivatives) for from arrival and velocity orthogonal grids
or from a pre-calculated ray. Every function exists in 2 flavors, 2D and 3D
which are differentiated by the coresponding 2, 3 suffix.

The intend of this module is to interface the low level C++ ArrayDescriptor and
OrthogonalGrid directly from numpy arrays.

The only interface provided is for double precision arithmetic even though
the low level C++ function are template method and could do anything.

All the raytracing methods use a Runge-Kutta (RK4) integrator.
"""


from eikonal.solver cimport *

cimport numpy as cnp
import  numpy as np

#
# I am not sure this is necessary since cython docs are all but precise on that
# particular matter. It really worked both ways. Decided to kept it for the
# sake of safety.
#
cdef extern from "numpy/arrayobject.h":
    void import_array() except *
import_array()

%for i in [2, 3]:
def sensivity${i}(cnp.ndarray[double, ndim = ${i}, mode = 'c'] arrival not None,
                  cnp.ndarray[double, ndim = ${i}, mode = 'c'] velocity not None,
                  tuple start,
                  tuple finish,
                  cnp.ndarray[size_t, ndim = 1, mode = 'c'] indices not None,
                  cnp.ndarray[double, ndim = 1, mode = 'c'] sensivity not None,
                  float spacing,
                  double h = 1):
    """
    Low level Sensivity and Raytracing for ${i}D problems.

    This method calculates the grid sensivity for ray spawning from <start> to
    <finish>. The actual ray is never built and only indices and sensivity are
    returned.

    <indices> and <sensivity> are array buffer provided by the user. The
    Apprioriate size of the array are the maximum length of the ray (4 ** D)

    Returns :       Traveltime
                    Indices Array
                    Sensivity Array

    WARNING : There is no bound checking in both indices and sensivity. Provide
              buffer that can holds the Frechet Derivative are important or
              segfault will occur.

    WARNING : The position of the <finish> parameter should always match with
              a zero on the arrival grid.

    Parameters :    arrival  - A 2D/3D numpy array.
                    velocity - A 2D/3D numpy array the same size as arrival.
                    start    - A 2D/3D tuple
                    finish   - A 2D/3D tuple
                    indices  - A nunmpy buffer of
                    sensivity-
                    spacing  - The grid spacing
                    h        - Control the length step of the RK integration
    """


    cdef ArrayDescriptor_ulong cindices
    cindices.init(<size_t *>indices.shape, <size_t *>indices.data)

    cdef ArrayDescriptor_double csensivity
    csensivity.init(<size_t *>sensivity.shape, <double *>sensivity.data)

    cdef OrthogonalGrid${i}_double carrival
    carrival.init(<size_t *>arrival.shape, <double *>arrival.data, spacing)

    cdef ArrayDescriptor${i}_double cvelocity
    cvelocity.init(<size_t *>velocity.shape, <double *>velocity.data)

    cdef doublev${i} cstart
    for j in range(${i}):
        cstart.data[j] = start[j]

    cdef doublev${i} cfinish
    for j in range(${i}):
        cfinish.data[j] = finish[j]


    cdef double tt

    with nogil:
        tt = sensivity_RK${i}(carrival, cvelocity, cstart, cfinish,  cindices, csensivity, h)
    return tt, indices[:cindices.shape[0]].copy(), sensivity[:csensivity.shape[0]].copy()

def traveltime${i}(cnp.ndarray[double, ndim = ${i}, mode = 'c'] arrival not None,
                   cnp.ndarray[double, ndim = ${i}, mode = 'c'] velocity not None,
                   tuple start,
                   tuple finish,
                   float spacing,
                   double h = 1):
    """
    Low level TravelTime and Raytracing.


    """
    cdef OrthogonalGrid${i}_double carrival
    carrival.init(<size_t *>arrival.shape, <double *>arrival.data, spacing)

    cdef ArrayDescriptor${i}_double cvelocity
    cvelocity.init(<size_t *>velocity.shape, <double *>velocity.data)

    cdef doublev${i} cstart
    for j in range(${i}):
        cstart.data[j] = start[j]

    cdef doublev${i} cfinish
    for j in range(${i}):
        cfinish.data[j] = finish[j]

    cdef double tt

    with nogil:
        tt = traveltime_RK${i}(carrival, cvelocity, cstart, cfinish, h)
    return tt

def raytrace${i}(cnp.ndarray[double, ndim = ${i}, mode = 'c'] arrival not None,
                 cnp.ndarray[double, ndim = ${i}, mode = 'c'] velocity not None,
                 tuple start,
                 tuple finish,
                 cnp.ndarray[double,    ndim = 2, mode = 'c'] path not None,
                 float spacing,
                 double h = 1):
    """
    Low level Raytracing.
    """

    cdef OrthogonalGrid${i}_double carrival
    carrival.init(<size_t *>arrival.shape, <double *>arrival.data, spacing)

    cdef ArrayDescriptor${i}_double cvelocity
    cvelocity.init(<size_t *>velocity.shape, <double *>velocity.data)

    cdef ArrayDescriptor_doublev${i} cpath
    cpath.init(<size_t *>path.shape, <doublev${i} *>path.data)

    cdef doublev${i} cstart
    for j in range(${i}):
        cstart.data[j] = start[j]

    cdef doublev${i} cfinish
    for j in range(${i}):
        cfinish.data[j] = finish[j]

    cdef double tt

    with nogil:
        tt = raytrace_RK${i}(carrival, cvelocity, cstart, cfinish, cpath, h)

    return tt, path[:cpath.shape[0]].copy()

def ray_sensivity${i}(cnp.ndarray[double, ndim = ${i}, mode = 'c'] velocity not None,
                      cnp.ndarray[double,    ndim = 2, mode = 'c'] path not None,
                      cnp.ndarray[size_t, ndim = 1, mode = 'c'] indices not None,
                      cnp.ndarray[double, ndim = 1, mode = 'c'] sensivity not None,
                      float spacing):
    """
    Low level Sensivity (Frechet Derivative).
    """

    cdef ArrayDescriptor_ulong cindices
    cindices.init(<size_t *>indices.shape, <size_t *>indices.data)

    cdef ArrayDescriptor_double csensivity
    csensivity.init(<size_t *>sensivity.shape, <double *>sensivity.data)


    cdef OrthogonalGrid${i}_double cvelocity
    cvelocity.init(<size_t *>velocity.shape, <double *>velocity.data, spacing)

    cdef ArrayDescriptor_doublev${i} cpath
    cpath.init(<size_t *>path.shape, <doublev${i} *>path.data)

    cdef double tt
    with nogil:
        tt = c_ray_sensivity${i}(cvelocity, cpath, cindices, csensivity)
    return tt, indices[:cindices.shape[0]].copy(), sensivity[:csensivity.shape[0]].copy()

def ray_traveltime${i}(cnp.ndarray[double, ndim = ${i}, mode = 'c'] velocity not None,
                       cnp.ndarray[double,    ndim = 2, mode = 'c'] path not None,
                       float spacing):
    """
    Low level Traveltime.
    """


    cdef OrthogonalGrid${i}_double cvelocity
    cvelocity.init(<size_t *>velocity.shape, <double *>velocity.data, spacing)

    cdef ArrayDescriptor_doublev${i} cpath
    cpath.init(<size_t *>path.shape, <doublev${i} *>path.data)

    cdef double tt
    with nogil:
        tt = c_ray_traveltime${i}(cvelocity, cpath)
    return tt
%endfor


%for fname in ['traveltime', 'raytrace', 'sensivity',]:
def ${fname}(cnp.ndarray arrival, cnp.ndarray velocity, *args, **kw):
    """
    Dispatcher for ${fname} for 2D or 3D array.

    Also enforce the size matching between arrival and velocity grids.
    """
    if arrival.ndim != velocity.ndim:
        raise ValueError("Velocity and Arrival grids must have the same dimensions")
    for i in range(arrival.ndim):
        if arrival.shape[i] != velocity.shape[i]:
            raise ValueError("Velocity and Arrival Grid must have the same shape")

%for i in [2, 3]:
    if ${i} == arrival.ndim:
        return ${fname}${i}(arrival, velocity,  *args, **kw)
%endfor

%endfor

%for fname in ['ray_sensivity', 'ray_traveltime']:
def ${fname}(cnp.ndarray velocity, *args, **kw):
    """
    Dispatcher for ${fname} for 2D or 3D array.
    """
%for i in [2, 3]:
    if ${i} == velocity.ndim:
        return ${fname}${i}(velocity,  *args, **kw)
%endfor

%endfor



# vim: filetype=pyrex
#
#

