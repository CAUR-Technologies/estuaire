#
# @Author : Jean-Pascal Mercier <jean-pascal.mercier@agsis.com>
#
# @Copyright (C) 2010 Jean-Pascal Mercier
#
# All rights reserved.
#

from solver cimport *

cdef extern from "raytrace.hpp":
    cdef double raytrace_RK2(OrthogonalGrid2_double arrival,
                      ArrayDescriptor2_double velocity,
                      doublev2 start,
                      doublev2 finish,
                      ArrayDescriptor_doublev2 path,
                      double h) nogil

    cdef double raytrace_RK3(OrthogonalGrid3_double arrival,
                      ArrayDescriptor3_double velocity,
                      doublev3 start,
                      doublev3 finish,
                      ArrayDescriptor_doublev3 path,
                      double h) nogil

    cdef double sensivity_RK2(OrthogonalGrid2_double arrival,
                                ArrayDescriptor2_double velocity,
                                doublev2 start,
                                doublev2 finish,
                                ArrayDescriptor_ulong indices,
                                ArrayDescriptor_double sensivity,
                                double h) nogil

    cdef double sensivity_RK3(OrthogonalGrid3_double arrival,
                                ArrayDescriptor3_double velocity,
                                doublev3 start,
                                doublev3 finish,
                                ArrayDescriptor_ulong indices,
                                ArrayDescriptor_double sensivity,
                                double h) nogil

    cdef double traveltime_RK3(OrthogonalGrid3_double arrival,
                                ArrayDescriptor3_double velocity,
                                doublev3 start,
                                doublev3 finish,
                                double h) nogil

    cdef double traveltime_RK2(OrthogonalGrid2_double arrival,
                                ArrayDescriptor2_double velocity,
                                doublev2 start,
                                doublev2 finish,
                                double h) nogil

    cdef double c_ray_sensivity2(OrthogonalGrid2_double velocity,
                                 ArrayDescriptor_doublev2 path,
                                 ArrayDescriptor_ulong indices,
                                 ArrayDescriptor_double sensivity) nogil

    cdef double c_ray_sensivity3(OrthogonalGrid3_double velocity,
                                 ArrayDescriptor_doublev3 path,
                                 ArrayDescriptor_ulong indices,
                                 ArrayDescriptor_double sensivity) nogil

    cdef double c_ray_traveltime3(OrthogonalGrid3_double velocity,
                                  ArrayDescriptor_doublev3 path) nogil

    cdef double c_ray_traveltime2(OrthogonalGrid2_double velocity,
                                  ArrayDescriptor_doublev2 path) nogil




