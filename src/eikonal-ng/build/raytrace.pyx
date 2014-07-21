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

def sensivity2(cnp.ndarray[double, ndim = 2, mode = 'c'] arrival not None,
                  cnp.ndarray[double, ndim = 2, mode = 'c'] velocity not None,
                  tuple start,
                  tuple finish,
                  cnp.ndarray[size_t, ndim = 1, mode = 'c'] indices not None,
                  cnp.ndarray[double, ndim = 1, mode = 'c'] sensivity not None,
                  float spacing,
                  double h = 1):
    """
    Low level Sensivity and Raytracing for 2D problems.

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

    cdef OrthogonalGrid2_double carrival
    carrival.init(<size_t *>arrival.shape, <double *>arrival.data, spacing)

    cdef ArrayDescriptor2_double cvelocity
    cvelocity.init(<size_t *>velocity.shape, <double *>velocity.data)

    cdef doublev2 cstart
    for j in xrange(2):
        cstart.data[j] = start[j]

    cdef doublev2 cfinish
    for j in xrange(2):
        cfinish.data[j] = finish[j]


    cdef double tt

    with nogil:
        tt = sensivity_RK2(carrival, cvelocity, cstart, cfinish,  cindices, csensivity, h)
    return tt, indices[:cindices.shape[0]].copy(), sensivity[:csensivity.shape[0]].copy()

def traveltime2(cnp.ndarray[double, ndim = 2, mode = 'c'] arrival not None,
                   cnp.ndarray[double, ndim = 2, mode = 'c'] velocity not None,
                   tuple start,
                   tuple finish,
                   float spacing,
                   double h = 1):
    """
    Low level TravelTime and Raytracing.


    """
    cdef OrthogonalGrid2_double carrival
    carrival.init(<size_t *>arrival.shape, <double *>arrival.data, spacing)

    cdef ArrayDescriptor2_double cvelocity
    cvelocity.init(<size_t *>velocity.shape, <double *>velocity.data)

    cdef doublev2 cstart
    for j in xrange(2):
        cstart.data[j] = start[j]

    cdef doublev2 cfinish
    for j in xrange(2):
        cfinish.data[j] = finish[j]

    cdef double tt

    with nogil:
        tt = traveltime_RK2(carrival, cvelocity, cstart, cfinish, h)
    return tt

def raytrace2(cnp.ndarray[double, ndim = 2, mode = 'c'] arrival not None,
                 cnp.ndarray[double, ndim = 2, mode = 'c'] velocity not None,
                 tuple start,
                 tuple finish,
                 cnp.ndarray[double,    ndim = 2, mode = 'c'] path not None,
                 float spacing,
                 double h = 1):
    """
    Low level Raytracing.
    """

    cdef OrthogonalGrid2_double carrival
    carrival.init(<size_t *>arrival.shape, <double *>arrival.data, spacing)

    cdef ArrayDescriptor2_double cvelocity
    cvelocity.init(<size_t *>velocity.shape, <double *>velocity.data)

    cdef ArrayDescriptor_doublev2 cpath
    cpath.init(<size_t *>path.shape, <doublev2 *>path.data)

    cdef doublev2 cstart
    for j in xrange(2):
        cstart.data[j] = start[j]

    cdef doublev2 cfinish
    for j in xrange(2):
        cfinish.data[j] = finish[j]

    cdef double tt

    with nogil:
        tt = raytrace_RK2(carrival, cvelocity, cstart, cfinish, cpath, h)

    return tt, path[:cpath.shape[0]].copy()

def ray_sensivity2(cnp.ndarray[double, ndim = 2, mode = 'c'] velocity not None,
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


    cdef OrthogonalGrid2_double cvelocity
    cvelocity.init(<size_t *>velocity.shape, <double *>velocity.data, spacing)

    cdef ArrayDescriptor_doublev2 cpath
    cpath.init(<size_t *>path.shape, <doublev2 *>path.data)

    cdef double tt
    with nogil:
        tt = c_ray_sensivity2(cvelocity, cpath, cindices, csensivity)
    return tt, indices[:cindices.shape[0]].copy(), sensivity[:csensivity.shape[0]].copy()

def ray_traveltime2(cnp.ndarray[double, ndim = 2, mode = 'c'] velocity not None,
                       cnp.ndarray[double,    ndim = 2, mode = 'c'] path not None,
                       float spacing):
    """
    Low level Traveltime.
    """


    cdef OrthogonalGrid2_double cvelocity
    cvelocity.init(<size_t *>velocity.shape, <double *>velocity.data, spacing)

    cdef ArrayDescriptor_doublev2 cpath
    cpath.init(<size_t *>path.shape, <doublev2 *>path.data)

    cdef double tt
    with nogil:
        tt = c_ray_traveltime2(cvelocity, cpath)
    return tt
def sensivity3(cnp.ndarray[double, ndim = 3, mode = 'c'] arrival not None,
                  cnp.ndarray[double, ndim = 3, mode = 'c'] velocity not None,
                  tuple start,
                  tuple finish,
                  cnp.ndarray[size_t, ndim = 1, mode = 'c'] indices not None,
                  cnp.ndarray[double, ndim = 1, mode = 'c'] sensivity not None,
                  float spacing,
                  double h = 1):
    """
    Low level Sensivity and Raytracing for 3D problems.

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

    cdef OrthogonalGrid3_double carrival
    carrival.init(<size_t *>arrival.shape, <double *>arrival.data, spacing)

    cdef ArrayDescriptor3_double cvelocity
    cvelocity.init(<size_t *>velocity.shape, <double *>velocity.data)

    cdef doublev3 cstart
    for j in xrange(3):
        cstart.data[j] = start[j]

    cdef doublev3 cfinish
    for j in xrange(3):
        cfinish.data[j] = finish[j]


    cdef double tt

    with nogil:
        tt = sensivity_RK3(carrival, cvelocity, cstart, cfinish,  cindices, csensivity, h)
    return tt, indices[:cindices.shape[0]].copy(), sensivity[:csensivity.shape[0]].copy()

def traveltime3(cnp.ndarray[double, ndim = 3, mode = 'c'] arrival not None,
                   cnp.ndarray[double, ndim = 3, mode = 'c'] velocity not None,
                   tuple start,
                   tuple finish,
                   float spacing,
                   double h = 1):
    """
    Low level TravelTime and Raytracing.


    """
    cdef OrthogonalGrid3_double carrival
    carrival.init(<size_t *>arrival.shape, <double *>arrival.data, spacing)

    cdef ArrayDescriptor3_double cvelocity
    cvelocity.init(<size_t *>velocity.shape, <double *>velocity.data)

    cdef doublev3 cstart
    for j in xrange(3):
        cstart.data[j] = start[j]

    cdef doublev3 cfinish
    for j in xrange(3):
        cfinish.data[j] = finish[j]

    cdef double tt

    with nogil:
        tt = traveltime_RK3(carrival, cvelocity, cstart, cfinish, h)
    return tt

def raytrace3(cnp.ndarray[double, ndim = 3, mode = 'c'] arrival not None,
                 cnp.ndarray[double, ndim = 3, mode = 'c'] velocity not None,
                 tuple start,
                 tuple finish,
                 cnp.ndarray[double,    ndim = 2, mode = 'c'] path not None,
                 float spacing,
                 double h = 1):
    """
    Low level Raytracing.
    """

    cdef OrthogonalGrid3_double carrival
    carrival.init(<size_t *>arrival.shape, <double *>arrival.data, spacing)

    cdef ArrayDescriptor3_double cvelocity
    cvelocity.init(<size_t *>velocity.shape, <double *>velocity.data)

    cdef ArrayDescriptor_doublev3 cpath
    cpath.init(<size_t *>path.shape, <doublev3 *>path.data)

    cdef doublev3 cstart
    for j in xrange(3):
        cstart.data[j] = start[j]

    cdef doublev3 cfinish
    for j in xrange(3):
        cfinish.data[j] = finish[j]

    cdef double tt

    with nogil:
        tt = raytrace_RK3(carrival, cvelocity, cstart, cfinish, cpath, h)

    return tt, path[:cpath.shape[0]].copy()

def ray_sensivity3(cnp.ndarray[double, ndim = 3, mode = 'c'] velocity not None,
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


    cdef OrthogonalGrid3_double cvelocity
    cvelocity.init(<size_t *>velocity.shape, <double *>velocity.data, spacing)

    cdef ArrayDescriptor_doublev3 cpath
    cpath.init(<size_t *>path.shape, <doublev3 *>path.data)

    cdef double tt
    with nogil:
        tt = c_ray_sensivity3(cvelocity, cpath, cindices, csensivity)
    return tt, indices[:cindices.shape[0]].copy(), sensivity[:csensivity.shape[0]].copy()

def ray_traveltime3(cnp.ndarray[double, ndim = 3, mode = 'c'] velocity not None,
                       cnp.ndarray[double,    ndim = 2, mode = 'c'] path not None,
                       float spacing):
    """
    Low level Traveltime.
    """


    cdef OrthogonalGrid3_double cvelocity
    cvelocity.init(<size_t *>velocity.shape, <double *>velocity.data, spacing)

    cdef ArrayDescriptor_doublev3 cpath
    cpath.init(<size_t *>path.shape, <doublev3 *>path.data)

    cdef double tt
    with nogil:
        tt = c_ray_traveltime3(cvelocity, cpath)
    return tt


def traveltime(cnp.ndarray arrival, cnp.ndarray velocity, *args, **kw):
    """
    Dispatcher for traveltime for 2D or 3D array.

    Also enforce the size matching between arrival and velocity grids.
    """
    if arrival.ndim != velocity.ndim:
        raise ValueError("Velocity and Arrival grids must have the same dimensions")
    for i in xrange(arrival.ndim):
        if arrival.shape[i] != velocity.shape[i]:
            raise ValueError("Velocity and Arrival Grid must have the same shape")

    if 2 == arrival.ndim:
        return traveltime2(arrival, velocity,  *args, **kw)
    if 3 == arrival.ndim:
        return traveltime3(arrival, velocity,  *args, **kw)

def raytrace(cnp.ndarray arrival, cnp.ndarray velocity, *args, **kw):
    """
    Dispatcher for raytrace for 2D or 3D array.

    Also enforce the size matching between arrival and velocity grids.
    """
    if arrival.ndim != velocity.ndim:
        raise ValueError("Velocity and Arrival grids must have the same dimensions")
    for i in xrange(arrival.ndim):
        if arrival.shape[i] != velocity.shape[i]:
            raise ValueError("Velocity and Arrival Grid must have the same shape")

    if 2 == arrival.ndim:
        return raytrace2(arrival, velocity,  *args, **kw)
    if 3 == arrival.ndim:
        return raytrace3(arrival, velocity,  *args, **kw)

def sensivity(cnp.ndarray arrival, cnp.ndarray velocity, *args, **kw):
    """
    Dispatcher for sensivity for 2D or 3D array.

    Also enforce the size matching between arrival and velocity grids.
    """
    if arrival.ndim != velocity.ndim:
        raise ValueError("Velocity and Arrival grids must have the same dimensions")
    for i in xrange(arrival.ndim):
        if arrival.shape[i] != velocity.shape[i]:
            raise ValueError("Velocity and Arrival Grid must have the same shape")

    if 2 == arrival.ndim:
        return sensivity2(arrival, velocity,  *args, **kw)
    if 3 == arrival.ndim:
        return sensivity3(arrival, velocity,  *args, **kw)


def ray_sensivity(cnp.ndarray velocity, *args, **kw):
    """
    Dispatcher for ray_sensivity for 2D or 3D array.
    """
    if 2 == velocity.ndim:
        return ray_sensivity2(velocity,  *args, **kw)
    if 3 == velocity.ndim:
        return ray_sensivity3(velocity,  *args, **kw)

def ray_traveltime(cnp.ndarray velocity, *args, **kw):
    """
    Dispatcher for ray_traveltime for 2D or 3D array.
    """
    if 2 == velocity.ndim:
        return ray_traveltime2(velocity,  *args, **kw)
    if 3 == velocity.ndim:
        return ray_traveltime3(velocity,  *args, **kw)




# vim: filetype=pyrex
#
#


