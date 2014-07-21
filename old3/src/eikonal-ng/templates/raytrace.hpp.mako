/**
 * @Author : Jean-Pascal Mercier <jean-pascal.mercier@agsis.com>
 *
 * @Copyright (C) Jean-Pascal Mercier
 *
 * All rights reserved.
 *
 *
 *
 * This file contains Raytracing and Frechet Derivative (Sensitivity)
 * procedures. The Algorithms implements the Runge Kutta Method using the
 * optimal coefficient from [McLachlan and Atela, 1992]. Runge-Kutta
 * coefficients are provided for RK-4 and RK-5 methods.
 *
 */



#ifndef RAYTRACING_HPP
#define RAYTRACING_HPP

#include <narray.hpp>
#include <vector>

/*
 * This is a raytracing method using a Runge-Kutta method
 * scheme 
 */


template <std::size_t ORDER> 
struct rk
{
    static const double a[ORDER];
    static const double b[ORDER];
};


// Optimally precalculated Runge Kutta Coefficient
// See : The accuracy of symplectic integrators [McLachlan and Atela, 1992] for details
template <>
const double rk<4>::a[4] = { 0.5153528374311229364,
                            -0.085782019412973646,
                             0.4415830236164665242,
                             0.1288461583653841854 };
template <>
const double rk<4>::b[4] = { 0.1344961992774310892,
                            -0.2248198030794208058,
                             0.7563200005156682911,
                             0.3340036032863214255 };

template <>
const double rk<5>::a[5] = { 0.339839625839110000,
                            -0.088601336903027329,
                             0.5858564768259621188,
                             0.323580796554697634,
                             0.4423637942197494587 };
template <>
const double rk<5>::b[5] = { 0.1193900292875672758,
                             0.6989273703824752308,
                            -0.1713123582716007754,
                             0.010705818482359840,
                            -0.0589796254980311632};





template <std::size_t R, typename T, std::size_t  D, bool BC>
double traveltime_RK(const OrthogonalGrid<T, D, BC> &arrival,
                   const ArrayDescriptor<T, D, BC> &velocity,
                   const agsis::vect<double, D> &start,
                   const agsis::vect<double, D> &finish,
                   const double h_init = 1)

{
    ${raytrace(['tt += dl / read_cubic_bspline(velocity, last + 0.5 * dt);'])}
    return tt;
}

// This method calculate the ray 
//
//      arrival - The tim
template <std::size_t R, typename T, std::size_t  D, bool BC>
double raytrace_RK(const OrthogonalGrid<T, D, BC> &arrival,
                   const ArrayDescriptor<T, D, BC> &velocity,
                   const agsis::vect<double, D> &start,
                   const agsis::vect<double, D> &finish,
                   ArrayDescriptor<agsis::vect<double, D>, 1, BC> &path,
                   const double h_init = 1)

{
    size_t i = 0;

    ${raytrace(['tt += dl / read_cubic_bspline(velocity, last + 0.5 * dt);', 'path[i++] = last;'])}

    path[i++] = start;
    path.shape[0] = i;
    return tt;
}

template <std::size_t R, typename T, std::size_t  D, bool BC>
double sensivity_RK(const OrthogonalGrid<T, D, BC> &arrival,
                    const ArrayDescriptor<T, D, BC> &velocity,
                    const agsis::vect<double, D> &start,
                    const agsis::vect<double, D> &finish,
                    ArrayDescriptor<size_t, 1, BC> &indices,
                    ArrayDescriptor<double, 1, BC> &sensivity,
                    double h_init = 1)
{
    size_t it = 0;

    ${raytrace(['tt += integrate_segment(indices, sensivity, it, velocity, last + 0.5 * dt, dl) ;'])}

    sensivity.shape[0] = it;
    indices.shape[0] = it;
    return tt;
}

template <typename T, std::size_t D, bool BC>
double c_ray_sensivity(const OrthogonalGrid<T, D, BC> &velocity,
                       const ArrayDescriptor<agsis::vect<double, D> , 1, BC> &raypath,
                       ArrayDescriptor<size_t, 1, BC> &indices,
                       ArrayDescriptor<double, 1, BC> &sensivity)
{
    using namespace agsis;
    using namespace std;

    double dl;

    size_t it = 0;
    double tt = 0;
    for (size_t i = 1; i < raypath.size(); i++)
    {
        const vect<double, D> dt = raypath[i] - raypath[i - 1];

        dl = norm(dt) * velocity.spacing;

        tt += integrate_segment(indices, sensivity, it, velocity, raypath[i - 1] + 0.5 * dt, dl);

    }


    sensivity.shape[0] = it;
    indices.shape[0] = it;

    return tt;

}

template <typename T, std::size_t D, bool BC>
double c_ray_traveltime(const OrthogonalGrid<T, D, BC> &velocity,
                        const ArrayDescriptor<agsis::vect<double, D> , 1, BC> &raypath)
{
    using namespace agsis;
    using namespace std;

    double dl;

    double tt = 0;
    for (size_t i = 1; i < raypath.size(); i++)
    {
        agsis::vect<double, D> dt = raypath[i] - raypath[i - 1];

        dl = norm(dt) * velocity.spacing;

        tt += dl / read_cubic_bspline(velocity, (raypath[i - 1] + 0.5 * dt));
    }
    return tt;
}

template <std::size_t N, typename T, size_t D, bool RC>
struct IntegrateSegment
{
    inline static double compute(const size_t index, const double weight,
                                 const ArrayDescriptor<T, D, RC> &grid,
                                 const agsis::vect<size_t, D> &fv,
                                 const agsis::vect<double, 4> c[D],
                                 ArrayDescriptor<size_t, 1, RC> &indices,
                                 ArrayDescriptor<double, 1, RC> &sensivity,
                                 size_t &it)
    {
        double value = 0;
        size_t next_index;
        for (int i = -1; i < 3;)
        {
            next_index = index + grid.strides[N] * \
                         clamp(0, grid.shape[N] - 1, static_cast<int>(fv[N]) + i);
            value += IntegrateSegment<N + 1, T, D, RC>::compute(next_index, weight * c[N][i + 1],
                                                                grid, fv, c, indices, sensivity, it);
            i++;

        }
        return value;
    }
};

template <typename T, size_t D, bool RC>
struct IntegrateSegment<D, T, D, RC>
{
    inline static double compute(const size_t index, const double weight,
                                 const ArrayDescriptor<T, D, RC> &grid,
                                 const agsis::vect<size_t, D>& fv,
                                 const agsis::vect<double, 4> c[D],
                                 ArrayDescriptor<size_t, 1, RC> &indices,
                                 ArrayDescriptor<double, 1, RC> &sensivity,
                                 size_t &it)
    {
        indices[it] = index;
        sensivity[it++] = -weight;
        return grid.data[index] * weight;
    }
};




template <typename T, std::size_t D, bool RC>
inline double integrate_segment(ArrayDescriptor<size_t, 1, RC> &indices,
                                ArrayDescriptor<double, 1, RC> &sensivity,
                                size_t &it,
                                const ArrayDescriptor<T, D, RC> &velocity,
                                const agsis::vect<double, D> &center,
                                const double dl)
{
    using namespace agsis;
    using namespace std;

    vect<size_t, D> fv = static_cast<vect<size_t, D> >(center); // clampped to grid
    vect<double, D> dx = center - static_cast<vect<double, D> >(fv);

    vect<double, 4> c[D];
    for (size_t i = 0; i < D; i++)
    {
        cubic_bspline_coefs(dx[i], c[i]);
    }

    double w[D];
    size_t index[D];
    size_t start_it = it;

    double value = IntegrateSegment<0, T, D, RC>::compute(0, 1.0, velocity, fv, c, indices, sensivity, it);

    for (size_t i = 0; i < D * D * D * D; i++) // IMPLICITLY UNROLLED
    {
        sensivity.data[start_it + i] *=  dl / (value * value);
    }

    return dl / value;

}






// Cython Goodness
% for i in [2, 3]:
template <typename T, std::size_t  D, bool BC>
double sensivity_RK${i}(const OrthogonalGrid<T, D, BC> &arrival,
                          const ArrayDescriptor<T, D, BC> &velocity,
                          const agsis::vect<double, D> &start,
                          const agsis::vect<double, D> &finish,
                          ArrayDescriptor<size_t, 1, BC> &indices,
                          ArrayDescriptor<double, 1, BC> &sensivity,
                          const double h = 1)
{
    return sensivity_RK<4>(arrival, velocity, start, finish, indices, sensivity, h);
}

template <typename T, std::size_t  D, bool BC>
double raytrace_RK${i}(const OrthogonalGrid<T, D, BC> &arrival,
                       const ArrayDescriptor<T, D, BC> &velocity,
                       const agsis::vect<double, D> &start,
                       const agsis::vect<double, D> &finish,
                       ArrayDescriptor<agsis::vect<double, D>, 1, BC> &path,
                       const double h = 1)
{
    return raytrace_RK<4>(arrival, velocity, start, finish, path, h);
}

template <typename T, std::size_t  D, bool BC>
double traveltime_RK${i}(const OrthogonalGrid<T, D, BC> &arrival,
                       const ArrayDescriptor<T, D, BC> &velocity,
                       const agsis::vect<double, D> &start,
                       const agsis::vect<double, D> &finish,
                       const double h = 1)
{
    return traveltime_RK<4>(arrival, velocity, start, finish, h);
}

template <typename T, std::size_t D, bool BC>
double c_ray_sensivity${i}(const OrthogonalGrid<T, D, BC> &velocity,
                           const ArrayDescriptor<agsis::vect<double, D> , 1, BC> &raypath,
                           ArrayDescriptor<size_t, 1, BC> &indices,
                           ArrayDescriptor<double, 1, BC> &sensivity)
{
    return c_ray_sensivity(velocity, raypath, indices, sensivity);
}

template <typename T, std::size_t D, bool BC>
double c_ray_traveltime${i}(const OrthogonalGrid<T, D, BC> &velocity,
                             const ArrayDescriptor<agsis::vect<double, D> , 1, BC> &raypath)
{
    return c_ray_traveltime(velocity, raypath);
}
% endfor

#endif // RAYTRACING_HPP

<%def name="raytrace(extra)">
    using namespace agsis;
    using namespace std;

    vect<double, D> cpos = finish;
    vect<double, D> last = cpos;
    vect<double, D> dt;

    double tt = 0;

    double dl;

    vect<double, D> gradient;
    do
    {

        for (size_t j = 0; j < R; j++) // LOOP IS UNROLLED
        {
            gradient = cubic_spline_gradient(arrival, cpos);
            cpos -= gradient * ((h_init / norm(gradient))* rk<R>::b[j]);
        }

        dt = cpos - last;
        dl = norm(dt) * arrival.spacing;

%for e in extra:
        ${e}
%endfor
            //cerr << cpos << n << endl;

        last = cpos;


    }
    while((norm(cpos - start) > max(1.2, h_init)));
    // The value 1.2 was defined empirically and should be taken with an absolute
    // care. The curvature of the velocity field increase when approching the
    // zeros and generate a h value higher than expected. It could totally
    // messed up the cpos could lead to segfault.

    dt = start - cpos;
    dl = norm(dt) * arrival.spacing;
    cpos = start;
    //cerr << cpos << endl << endl;

%for e in extra:
    ${e}
%endfor

</%def>

/* vim: filetype=cpp
 *
 */
