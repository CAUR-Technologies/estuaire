#
# @Author : Jean-Pascal Mercier <jean-pascal.mercier@agsis.com>
#
# @Copyright (C) 2010 Jean-Pascal Mercier
#
# All rights reserved.
#

from libcpp cimport bool

cdef extern from "nvect.hpp":
    ctypedef struct ulongv2 "agsis::vect<std::size_t int, 2>":
        size_t data[2]

    ctypedef struct ulongv3 "agsis::vect<std::size_t, 3>":
        size_t data[3]

    ctypedef struct doublev3 "agsis::vect<double, 3>":
        double data[3]

    ctypedef struct doublev2 "agsis::vect<double, 2>":
        double data[2]

cdef extern from "narray.hpp":
    ctypedef struct ArrayDescriptor1_double "ArrayDescriptor<double, 1>":
        void init(size_t shape[1], double *data)
        size_t shape[1]

    ctypedef struct ArrayDescriptor2_double "ArrayDescriptor<double, 2>":
        void init(size_t shape[2], double *data)
        size_t shape[2]

    ctypedef struct ArrayDescriptor3_double "ArrayDescriptor<double, 3>":
        void init(size_t shape[3], double *data)
        size_t shape[3]

    ctypedef struct ArrayDescriptor3_ulong "ArrayDescriptor<size_t, 3>":
        void init(size_t shape[3], size_t *data)
        size_t shape[3]

    ctypedef struct ArrayDescriptor2_ulong "ArrayDescriptor<size_t, 2>":
        void init(size_t shape[2], size_t *data)
        size_t shape[2]


    ctypedef struct ArrayDescriptor_ulong "ArrayDescriptor<size_t, 1>":
        void init(size_t shape[1], size_t *data)
        size_t shape[1]

    ctypedef struct ArrayDescriptor_double "ArrayDescriptor<double, 1>":
        void init(size_t shape[1], double *data)
        size_t shape[1]


    ctypedef struct OrthogonalGrid3_double "OrthogonalGrid<double, 3>":
        void init(size_t shape[3], double *data, double spacing)
        size_t shape[3]

    ctypedef struct OrthogonalGrid2_double "OrthogonalGrid<double, 2>":
        void init(size_t shape[3], double *data, double spacing)
        size_t shape[3]

    ctypedef struct ArrayDescriptor_doublev2 "ArrayDescriptor<agsis::vect<double, 2>, 1>":
        void init(size_t shape[1], doublev2* data)
        size_t shape[1]

    ctypedef struct ArrayDescriptor_doublev3 "ArrayDescriptor<agsis::vect<double, 3>, 1>":
        void init(size_t shape[1], doublev3* data)
        size_t shape[1]

    ctypedef struct ArrayDescriptor2_ulongv "ArrayDescriptor<agsis::vect<std::size_t, 2>, 1>":
        void init(size_t shape[2], ulongv2* data)
        size_t shape[2]

    ctypedef struct ArrayDescriptor3_ulongv "ArrayDescriptor<agsis::vect<std::size_t, 3>, 1 >":
        void init(size_t shape[3], ulongv3* data)
        size_t shape[3]


cdef extern from "solver.hpp":
    cdef enum MarchingTag:
        FAR         = 0
        NARROW_BAND = 1
        KNOWN       = 2

    ctypedef struct ArrayDescriptor2_tag "ArrayDescriptor<MarchingTag, 2>":
        void init(size_t shape[2], MarchingTag *data)

    ctypedef struct ArrayDescriptor3_tag "ArrayDescriptor<MarchingTag, 3>":
        void init(size_t shape[3], MarchingTag *data)



    cdef void FMM3(ArrayDescriptor_doublev3,
                   ArrayDescriptor3_tag,
                   ArrayDescriptor3_double,
                   OrthogonalGrid3_double, bool) nogil

    cdef void FMM2(ArrayDescriptor_doublev2,
                   ArrayDescriptor2_tag,
                   ArrayDescriptor2_double,
                   OrthogonalGrid2_double, bool) nogil

    cdef void IFMM3(ArrayDescriptor3_tag,
                    ArrayDescriptor3_double,
                    OrthogonalGrid3_double, bool) nogil

    cdef void IFMM2(ArrayDescriptor2_tag,
                    ArrayDescriptor2_double,
                    OrthogonalGrid2_double, bool) nogil




cdef extern from "raytrace.hpp":
    cdef double sensibility_RK3(ArrayDescriptor3_double,
                                ArrayDescriptor3_double,
                                ulongv3,
                                ulongv3,
                                ArrayDescriptor3_ulong,
                                ArrayDescriprot3_ulong)

    cdef double sensibility_RK2(ArrayDescriptor2_double,
                                ArrayDescriptor2_double,
                                ulongv2,
                                ulongv2,
                                ArrayDescriptor2_ulong,
                                ArrayDescriprot2_ulong)
