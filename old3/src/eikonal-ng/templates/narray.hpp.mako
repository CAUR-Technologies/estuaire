/**
 *
 * @Author : Jean-Pascal Mercier <jean-pascal.mercier@agsis.com>
 *
 * @Copyright (C) 2010 Jean-Pascal Mercier
 *
 * All rights reserved.
 *
 *
 * This file contains the Multi-Dimentional array manipulation primitives.
 *
 */

#ifndef NARRAY_HPP
#define NARRAY_HPP

#include <stdexcept>
#include <iostream>
#include <sstream>

#include <nvect.hpp>

template <typename T>
inline std::size_t clamp(int mi, int ma, T val)
{
    return (val < ma) ? ((val > 0) ? val : 0) : (ma - 1);
}



//
//      Object Definition
//
//
template <typename T, std::size_t D, bool RC = false>
class ArrayDescriptor {
    public:
        static const std::size_t    ndim = D;
        T *                         data;

        agsis::vect<std::size_t, D>        shape;
        agsis::vect<std::size_t, D>        strides;

        ArrayDescriptor()           {}
        ArrayDescriptor(const agsis::vect<std::size_t, D> &v, T *pointer);

        T& element(const agsis::vect<std::size_t, D> &v);

        // For Cython pleasure
        void init(const agsis::vect<std::size_t, D> &v, T* pointer);

        std::size_t size() const;
        void check_boundaries(const agsis::vect<std::size_t, D> &v) const;

        T& operator[](const agsis::vect<std::size_t, D> &i);
        T& operator[](const std::size_t i);
        const T& operator[](const agsis::vect<std::size_t, D> &i) const;
        const T& operator[](const std::size_t i) const;

        const T &element(const agsis::vect<std::size_t, D> &v) const;
        std::size_t index(const agsis::vect<std::size_t, D> &v) const;

        void print_values(std::ostream &out = std::cout) const;
};

template <typename T, std::size_t D, bool RC = false>
class OrthogonalGrid : public ArrayDescriptor<T, D, RC> {
    public :
        double spacing;
        OrthogonalGrid() {}
        OrthogonalGrid(const agsis::vect<std::size_t, D> &v, T *pointer, double sp) : ArrayDescriptor<T, D, RC>(v, pointer), spacing(sp) {}
        void init(const agsis::vect<std::size_t, D> &v, T* pointer, double sp) { spacing = sp; ArrayDescriptor<T, D, RC>::init(v, pointer); }
};

template <typename T, std::size_t D, bool RC>
std::size_t ArrayDescriptor<T, D, RC>::size() const
{
    size_t size = this->shape[0];
    for (size_t i = 1; i < D; i++)
        size *= this->shape[i];
    return size;
}



inline void cubic_bspline_coefs(double x, agsis::vect<double, 4> &coefs)
{
    coefs[0] = (1 - x) * (1 - x) * (1 - x) / 6.0; // -x^3 + 3x^2 - 3x + 1
    coefs[1] = 0.5 * (x * x * x) - (x * x) + 4.0 / 6.0; // x^3 - 6x^2 + 6 
    coefs[2] = 1.0 / 6.0 + 0.5 * ((x * x) + x - (x * x * x)); // -3x^4 + 4x^2 + 3^x + 1
    coefs[3] = (x * x * x) / 6.0; // 3x^3
}


template <typename T, std::size_t D, bool RC>
inline agsis::vect<double, D> cubic_spline_gradient(const ArrayDescriptor<T, D, RC> &grid, const agsis::vect<double, D> &p)
{
    using namespace agsis;
    using namespace std;

    const vect<size_t, D> fp = static_cast<vect<size_t, D> >(p);
    const vect<double, D> dl = p - static_cast<vect<double, D> >(fp);

    size_t index = grid.index(fp);

    const T p2 = grid[index];
    double p1, p3, p4;
    double i1d, i3d, i4d;

    vect<double, D> result;

    double upper_bound, lower_bound;

    for (int i = D - 1; i >= 0; i--)
    {
        upper_bound = index + (grid.shape[i] - fp[i] - 1) * grid.strides[i];
        lower_bound = index - fp[i] * grid.strides[i];

        i1d = static_cast<double>(index) -\
                           static_cast<double>(grid.strides[i]);

        i3d = static_cast<double>(index) +\
                           static_cast<double>(grid.strides[i]);

        i4d = i3d + static_cast<double>(grid.strides[i]);

        // This section is to ensure we have a derivative a little higher than
        // 0 at the border for raytracing algorithm
        if (!(i1d > lower_bound))
        {
            p3 = grid[static_cast<size_t>(i3d)];
            p4 = grid[static_cast<size_t>(i4d)];
            p1 = p2 + abs(p3 - p2) + 1;
        }
        else
        {
            p1 = grid[static_cast<size_t>(i1d)];
            if (!(i3d < upper_bound))
            {
                p3 = p2 + abs(p2 - p1);
                p4 = p3 + abs(p2 - p1) + 1;
            }
            else
            {
                p3 = grid[static_cast<size_t>(i3d)];
                if (!(i4d < upper_bound))
                    p4 = p4 + abs(p2 - p3) + 1;
                else
                    p4 = grid[static_cast<size_t>(i4d)];
            }
        }


        result[i] = ((-0.5 * dl[i] * dl[i] + dl[i] - 0.5) * p1 \
                     +  (1.5 * dl[i] * dl[i] - 2 * dl[i]) * p2 \
                     + (-1.5 * dl[i] * dl[i] + dl[i] + 0.5) * p3 \
                     + (0.5 * dl[i] * dl[i]) * p4);

    }
    return result;
}



// This object is is for template loop unrolling
template <std::size_t N, typename T, size_t D, bool RC>
struct BSplineValue
{
    inline static double compute(const size_t index, const double weight,
                                 const ArrayDescriptor<T, D, RC> &grid,
                                 const agsis::vect<size_t, D> &fv,
                                 const agsis::vect<double, 4> c[D])
    {
        double value = 0;
        for (int i = -1; i < 3;)
        {
            value += BSplineValue<N + 1, T, D, RC>::compute(index + grid.strides[N] \
                                                           * clamp(0, grid.shape[N] - 1,
                                                           static_cast<int>(fv[N]) + i),
                                                           weight * c[N][i + 1], grid, fv, c);
            i++;

        }
        return value;
    }
};

template <typename T, size_t D, bool RC>
struct BSplineValue<D, T, D, RC>
{
    inline static double compute(const size_t index, const double weight,
                                 const ArrayDescriptor<T, D, RC> &grid,
                                 const agsis::vect<size_t, D>& fv,
                                 const agsis::vect<double, 4> c[D])
    {
        return grid.data[index] * weight;
    }
};



template <typename T, size_t D, bool RC>
inline double read_cubic_bspline(const ArrayDescriptor<T, D, RC> &grid, const agsis::vect<double, D> & v)
{
    using namespace std;
    using namespace agsis;

    const vect<size_t, D> fv = static_cast<vect<size_t, D> >(v); // clampped to grid
    const vect<double, D> dx = v - static_cast<vect<double, D> >(fv);

    vect<double, 4> c[D];
    for (size_t i = 0; i < D; i++)
    {
        cubic_bspline_coefs(dx[i], c[i]);
    }

    return BSplineValue<0, T, D, RC>::compute(0, 1.0, grid, fv, c);
}



template <typename T, std::size_t D, bool RC>
void ArrayDescriptor<T, D, RC>::print_values(std::ostream &out) const
{
    for (size_t i = 0; i < this->size(); i++)
    {
        out << this->data[i] << std::endl;
    }

}


template <typename T, std::size_t D, bool RC>
inline void ArrayDescriptor<T, D, RC>::init(const agsis::vect<std::size_t, D> &ary_shape, T *pointer)
{
    this->data = pointer;
    this->shape[D - 1] = ary_shape[D -1];
    this->strides[D - 1] = 1;
    for (int i = static_cast<long int>(D) - 2; i >= 0; i--)
    {
        this->shape[i] = ary_shape[i];
        this->strides[i] = this->strides[i + 1] * this->shape[i + 1];
    }
}



template <typename T, std::size_t D, bool RC>
ArrayDescriptor<T, D, RC>::ArrayDescriptor(const agsis::vect<std::size_t, D> &ary_shape, T *pointer)
{
    this->init(ary_shape, pointer);
}

template <typename T, std::size_t D, bool RC>
std::size_t ArrayDescriptor<T, D, RC>::index(const agsis::vect<std::size_t, D> &v) const
{
    return dot(v, this->strides);
}



template <typename T, std::size_t D, bool RC>
inline const T& ArrayDescriptor<T, D, RC>::operator[](const agsis::vect<std::size_t, D> &i) const
{
    return this->element(i);
}

template <typename T, std::size_t D, bool RC>
inline const T& ArrayDescriptor<T, D, RC>::operator[](const std::size_t i) const
{
    return this->data[i];
}



template <typename T, std::size_t D, bool RC>
inline T& ArrayDescriptor<T, D, RC>::operator[](const agsis::vect<std::size_t, D> &i)
{
    return this->element(i);
}

template <typename T, std::size_t D, bool RC>
inline T& ArrayDescriptor<T, D, RC>::operator[](const std::size_t i)
{
    return this->data[i];
}

template <typename T, std::size_t D, bool RC>
inline void ArrayDescriptor<T, D, RC>::check_boundaries(const agsis::vect<std::size_t, D> &v) const
{
    using namespace std;
    if (any(v >= this->shape))
    {
            stringstream s;
            s << "Index out of range : " << v << ", " << this->shape;
            throw out_of_range(s.str());
    }
}

template <typename T, std::size_t D, bool RC>
inline T& ArrayDescriptor<T, D, RC>::element(const agsis::vect<std::size_t, D> &v)
{
    if (RC)
        this->check_boundaries(v);

    return this->data[dot(v, this->strides)];
}

template <typename T, std::size_t D, bool RC>
inline const T& ArrayDescriptor<T, D, RC>::element(const agsis::vect<std::size_t,D> &v) const
{
    if (RC)
        this->check_boundaries(v);
    return this->data[dot(v, this->strides)];
}


template <std::size_t N, typename T, std::size_t D, bool RC>
struct zoom_loop
{
    inline static void compute(const ArrayDescriptor<T, D, RC> &input,
                               ArrayDescriptor<T, D, RC> &output,
                               agsis::vect<double, D> &pos,
                               const size_t index)
    {
        size_t loop_index = index;
        const double stride = ((double)input.shape[N]) / ((double)output.shape[N]);
        for (pos[N] = 0; pos[N] < input.shape[N]; pos[N] += stride)
        {
            loop_index += output.strides[N];
            zoom_loop<N + 1, T, D, RC>::compute(input, output, pos, loop_index);
        }
    }
};


template <typename T, std::size_t D, bool RC>
struct zoom_loop<D, T, D, RC>
{
        inline static void compute(const ArrayDescriptor<T, D, RC> &input,
                               ArrayDescriptor<T, D, RC> &output,
                               agsis::vect<double, D> &pos,
                               const size_t index)
        {
            output.data[index] = read_cubic_bspline(input, pos);
        }

};


template <typename T, std::size_t D, bool RC>
void zoom(const ArrayDescriptor<T, D, RC> input, ArrayDescriptor<T, D, RC> output)
{
    agsis::vect<double, D> pos;
    for (size_t i = 0; i < D; i++)
        pos[i] = 0;
    zoom_loop<0, T, D, RC>::compute(input, output, pos, 0);
}



// Beautiful Output
template <typename T, std::size_t D, bool RC>
inline std::ostream &operator<<( std::ostream &out, const ArrayDescriptor<T, D, RC> &ary )
{
    out << "<ArrayDescriptor D=" << D;
    out << " shape=" << ary.shape << " strides=" << ary.strides;
    out << " boundchecking=" << (RC ? "True" : "False") << ">";
    return out;
}

template <typename T, std::size_t D, bool RC>
inline std::ostream &operator<<( std::ostream &out, const OrthogonalGrid<T, D, RC> &ary )
{
    out << "<OrthogonalGrid D=" << D;
    out << " spacing=" << ary.spacing;
    out << " shape=" << ary.shape << " strides=" << ary.strides;
    out << " boundchecking=" << (RC ? "True" : "False") << ">";
    return out;
}

#endif // NARRAY_HPP

/**
 * vim: filetype=cpp
 *
 */
