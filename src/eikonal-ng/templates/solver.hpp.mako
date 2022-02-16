/**
 *
 * @Author : Jean-Pascal Mercier <jean-pascal.mercier@agsis.com>
 *
 * @Copyright (C) 2010 Jean-Pascal Mercier
 *
 * All rights reserved.
 *
 * This file contains an Implementation of the Fast Marching Method for
 * solving the Eikonal Equations. Both 1st and 2nd order  implementations
 * are provided and are functionnal. The 2nd is recommended since
 * more accurate for a minimal performance cost.
 *
 * TODO Integrate in the agsis namespace
 *
 */





#ifndef SOLVER_HPP
#define SOLVER_HPP

#include <narray.hpp>
#include <map>
#include <limits>

enum MarchingTag {
    FAR         = 0,
    NARROW_BAND = 1,
    KNOWN       = 2,
};


template <typename T, typename T2>
class emultimap : public std::multimap<T, T2>
{
    public:
        void erase_atom(T key, T2 value);
};

// This Function remove the exact occurence of an object inside the map
// This is not part of the STL standard
template <typename T, typename T2>
inline void emultimap<T, T2>::erase_atom(T key, T2 value)
{
    typename emultimap<T, T2>::iterator it = this->find(key);
    while((*it).second != value)
    {
        it++;
    }
    this->erase(it);
}


enum FiniteDifferenceSide
{
    FORWARD = 1,
    BACKWARD = -1
};


// This function simply calculate the result coefficients of the quadratic 
// equation (ax - b) ** 2
inline void quadratic_square(const double a, const double b, agsis::vect<double, 3>& p)
{
    p[0] = a * a;
    p[1] = 2 * a * b;
    p[2] = b * b;
}

// This method evaluate the First order upwind finite difference scheme. The
// DSIDE template argument determinate if the function returns the forward or
// the backward difference.
template<FiniteDifferenceSide DSIDE>
struct FirstOrderSemiFiniteDifference
{
    static const int mult = 1;
    template <typename T, std::size_t D, bool BC>
    inline static double evaluate(const ArrayDescriptor<T, D, BC> &grid,
                                  const std::size_t dim,
                                  const std::size_t index)
    {
        return -grid[index + DSIDE * grid.strides[dim]];
    }

    template <typename T, std::size_t D, bool BC>
    inline static double evaluate(const ArrayDescriptor<T, D, BC> &grid,
                                  const std::size_t dim,
                                  const agsis::vect<std::size_t, D> &v)
    {
        using namespace std;

        size_t index = grid.index(v);
        return evaluate(grid, dim, index);
    }

};


// This method evaluate the second order finite difference scheme. The
// DSIDE template argument determinate if the function returns the 
// forward of the backward difference
//
template<FiniteDifferenceSide DSIDE>
struct SecondOrderSemiFiniteDifference
{
    static constexpr double mult = 1.5;
    template <typename T, std::size_t D, bool BC>
    inline static double evaluate(const ArrayDescriptor<T, D, BC> &grid,
                                  const std::size_t dim,
                                  const std::size_t index)
    {
        using namespace std;

        const size_t stride = DSIDE * grid.strides[dim];
        return (-2 * grid[index + stride] + 0.5 * grid[index + 2 * stride]);
    }

    template <typename T, std::size_t D, bool BC>
    inline static double evaluate(const ArrayDescriptor<T, D, BC> &grid,
                                  const std::size_t dim,
                                  const agsis::vect<std::size_t, D> &v)
    {
        return evaluate(grid, dim, grid.index(v));
    }
};



// This method is only for optimization purpose if we are sure we only wants
// one root and only one. The template argument side should be set to 1 or -1
// depending if we want a positive or negative root. (default = positive)
enum Side {
    POSITIVE = 1,
    NEGATIVE = -1
};

// This method solve the quadratic equation
template <typename T>
inline agsis::vect<double, 2> quadratic_solver(const agsis::vect<T, 3> &a)
{
    agsis::vect<double, 2> result;
    double root = sqrt(a[1] * a[1] - 4 * a[0] * a[2]);
    result[0] = (-a[1] + root) / (2 * a[0]);
    result[1] = (-a[1] - root) / (2 * a[0]);
    return result;
}


// This is the solver for one half of the quadratic equation. If we only
// need the positive or the negative root we should use this function since
// it is faster.
template <Side SIDE, typename T>
inline double half_quadratic_solver(const agsis::vect<T, 3> &a)
{
    double root = sqrtf(a[1] * a[1] - 4 * a[0] * a[2]);
    return (-a[1] + SIDE * root) / (2 * a[0]);
}

// This method solves the Eikonal Equation at the given position in the
// viscosity grid using the upwind scheme. This method do not explicitly check
// for boundaries so make sure the position is at least far enough from the
// boundaries as need by the given differencing operator. There are two
// explicit template arguments.
//      DTYPE which determined the types of finite difference scheme used.
//      DSIDE which is an optimization to determine the callee last
//            updated value. We limit the calculation to this dimension
//            dimension for optimization.
template <template <FiniteDifferenceSide> class DTYPE, FiniteDifferenceSide DSIDE, typename T, std::size_t D, bool BC>
void FMM_UpwindSolve(const std::size_t index,
                     ArrayDescriptor<MarchingTag, D, BC> &tag,
                     const ArrayDescriptor<T, D, BC>& viscosity,
                     OrthogonalGrid<T, D, BC> &arrival,
                     emultimap<T, std::size_t > &band,
                     const std::size_t dim)
{
    using namespace agsis;
    using namespace std;

    const MarchingTag point_tag = tag[index];
    if (point_tag == KNOWN)
        return;

    T solution = numeric_limits<T>::infinity();
    const T point_viscosity = viscosity[index];
    const T point_arrival = arrival[index];

    //2D - 3D Solution :  PLAIN UGLY TO BE CHANGED
    vect<double, 3> p1, p2, p3;
    quadratic_square(DTYPE<DSIDE>::mult, DTYPE<DSIDE>::evaluate(arrival, dim, index), p1);
    p1[2] -= (arrival.spacing * arrival.spacing) * (point_viscosity * point_viscosity);
    for (int j = D - 1; j >= 0; j--)
    {
        if (j != static_cast<int>(dim))
        {
            quadratic_square(DTYPE<FORWARD>::mult, DTYPE<FORWARD>::evaluate(arrival, j, index), p2);
            solution = min(solution, half_quadratic_solver<POSITIVE>(p1 + p2));
            if (D >= 3)
            {
                for (int k = j - 1; k >= 0; k--)
                {
                    if (k != static_cast<int>(dim))
                    {
                        quadratic_square(DTYPE<FORWARD>::mult, DTYPE<FORWARD>::evaluate(arrival, k, index), p3);
                        solution = min(solution, half_quadratic_solver<POSITIVE>(p1 + p2 + p3));

                        quadratic_square(DTYPE<BACKWARD>::mult, DTYPE<BACKWARD>::evaluate(arrival, k, index), p3);
                        solution = min(solution, half_quadratic_solver<POSITIVE>(p1 + p2 + p3));
                    }
                }
            }


            quadratic_square(DTYPE<BACKWARD>::mult, DTYPE<BACKWARD>::evaluate(arrival, j, index), p2);
            solution = min(solution, half_quadratic_solver<POSITIVE>(p1 + p2));
            if (D >= 3)
            {
                for (int k = j - 1; k >= 0; k--)
                {
                    if (k != static_cast<int>(dim))
                    {
                        quadratic_square(DTYPE<FORWARD>::mult, DTYPE<FORWARD>::evaluate(arrival, k, index), p3);
                        solution = min(solution, half_quadratic_solver<POSITIVE>(p1 + p2 + p3));

                        quadratic_square(DTYPE<BACKWARD>::mult, DTYPE<BACKWARD>::evaluate(arrival, k, index), p3);
                        solution = min(solution, half_quadratic_solver<POSITIVE>(p1 + p2 + p3));
                    }
                }
            }
        }
    }


    // 1D Solution
    solution = min(solution, arrival[index + DSIDE * arrival.strides[dim]] + arrival.spacing * point_viscosity);

    // We verify we are in the narrow band for the modification of the
    // heap if needed.
    if (point_tag == NARROW_BAND)
    {
        // Verifying if we have a better solution than the calculated one
        if (point_arrival <= solution)
            return;
        band.erase_atom(point_arrival, index);
    }
    else
    {
        tag[index] = NARROW_BAND;
    }

    // A new solution is found. Setting the new solution and inserting
    // the node in the Narrow Band
    arrival[index] = solution;
    band.insert(pair<const T, size_t >(solution, index));
}


% for order in ["First", "Second"]:
// This method calculate the solution of the eikonal equation for the given
// velocity grid. For the sake of optimization, the arrival grid should be
// pre-filled with infinity value. The tag array should also have boundaries
// pre-setted to known.
template<typename T, std::size_t D, bool BC>
void FMM_${order}Order(const ArrayDescriptor <agsis::vect<double, D>, 1, BC> &seeds,
         ArrayDescriptor<MarchingTag, D, BC> &tag,
         const ArrayDescriptor<T, D, BC> &viscosity,
         OrthogonalGrid<T, D, BC> &arrival)
{
    using namespace agsis;
    using namespace std;

    emultimap<T, size_t>band;

    // PLAIN UGLY TO BE CHANGED
    // This is the initialization of the narrow band
    for(size_t i = 0; i < seeds.shape[0]; i++)
    {
        const vect<double, D> s = seeds[i];
        vect<double, D> fs = floor(seeds[i]);

        size_t index = arrival.index(static_cast<vect<size_t, D> >(fs));
        const T point_viscosity = read_cubic_bspline(viscosity, s);

        arrival[index] = norm(fs - s) * arrival.spacing * point_viscosity;
        tag[index] = NARROW_BAND;
        band.insert(pair<const T, size_t>(arrival[index], index));
        if (D >= 2)
        {
            if (D >= 3)
            {
                fs[D - 3] += 1;
                index += arrival.strides[D - 3];
                arrival[index] = norm(fs - s) * arrival.spacing * point_viscosity;
                tag[index] = NARROW_BAND;
                band.insert(pair<const T, size_t>(arrival[index], index));

                fs[D - 3] -= 1;
                index -= arrival.strides[D - 3];

            }

            fs[D - 2] += 1;
            index += arrival.strides[D - 2];
            arrival[index] = norm(fs - s) * arrival.spacing * point_viscosity;
            tag[index] = NARROW_BAND;
            band.insert(pair<const T, size_t>(arrival[index], index));

            if (D >= 3)
            {
                fs[D - 3] += 1;
                index += arrival.strides[D - 3];
                arrival[index] = norm(fs - s) * arrival.spacing * point_viscosity;
                tag[index] = NARROW_BAND;
                band.insert(pair<const T, size_t>(arrival[index], index));

                fs[D - 3] -= 1;
                index -= arrival.strides[D - 3];

            }

            fs[D - 2] -= 1;
            index-= arrival.strides[D - 2];
        }

        index += arrival.strides[D -1];
        fs[D - 1] += 1;
        arrival[index] = norm(fs - s) * arrival.spacing * point_viscosity;
        tag[index] = NARROW_BAND;
        band.insert(pair<const T, size_t>(arrival[index], index));
        if (D >= 2)
        {
            if (D >= 3)
            {
                fs[D - 3] += 1;
                index += arrival.strides[D - 3];
                arrival[index] = norm(fs -s) * arrival.spacing * point_viscosity;
                tag[index] = NARROW_BAND;
                band.insert(pair<const T, size_t>(arrival[index], index));

                fs[D - 3] -= 1;
                index -= arrival.strides[D - 3];

            }

            fs[D - 2] += 1;
            index += arrival.strides[D - 2];
            tag[index] = NARROW_BAND;
            arrival[index] = norm(fs -s) * arrival.spacing * point_viscosity;
            band.insert(pair<const T, size_t>(arrival[index], index));

            if (D >= 3)
            {
                fs[D - 3] += 1;
                index += arrival.strides[D - 3];
                arrival[index] = norm(fs - s) * arrival.spacing * point_viscosity;
                tag[index] = NARROW_BAND;
                band.insert(pair<const T, size_t>(arrival[index], index));
            }
        }
    }


    // The actual solver, We loop until there is no node left in the narrow
    // band
    ${eikonal_update(order = order)}
}

template<typename T, std::size_t D, bool BC>
void IFMM_${order}Order(ArrayDescriptor<MarchingTag, D, BC> &tag,
          const ArrayDescriptor<T, D, BC> &viscosity,
          OrthogonalGrid<T, D, BC> &arrival)
{
    using namespace agsis;
    using namespace std;

    emultimap<T, size_t>band;

    for (size_t index = 0; index < tag.size(); index++)
    {

        if (tag.data[index] == NARROW_BAND)
        {
            band.insert(pair<const T, size_t>(arrival.data[index], index));
        }
    }

    // The actual solver, We loop until there is no node left in the narrow
    // band
    ${eikonal_update(order = order)}

}

% endfor
// Cython Goodness
//
% for i in [2,3]:
template<typename T, std::size_t D, bool BC>
void inline FMM${i}(const ArrayDescriptor <agsis::vect<double, D>, 1, BC> &seeds,
                    ArrayDescriptor<MarchingTag, D, BC> &tag,
                    const ArrayDescriptor<T, D, BC> &viscosity,
                    OrthogonalGrid<T, D, BC> &arrival, bool second_order = true)
{
    if (second_order)
        FMM_SecondOrder(seeds, tag, viscosity, arrival);
    else
        FMM_FirstOrder(seeds, tag, viscosity, arrival);
}

template<typename T, std::size_t D, bool BC>
void inline IFMM${i}(ArrayDescriptor<MarchingTag, D, BC> &tag,
                    const ArrayDescriptor<T, D, BC> &viscosity,
                    OrthogonalGrid<T, D, BC> &arrival, bool second_order = true)
{
    if (second_order)
        IFMM_SecondOrder(tag, viscosity, arrival);
    else
        IFMM_FirstOrder(tag, viscosity, arrival);
}
% endfor


#endif // SOLVER_HPP

<%def name="eikonal_update(order = 'Second')">
    typename emultimap<T, size_t >::iterator it = band.begin();
    while(it != band.end())
    {
        const size_t index = (*it).second;

        // Setting the smallest value to KNOWN
        tag[index] = KNOWN;

        // Removing from the band
        band.erase(it);

        // Neighborhood
        for (int i = D - 1; i >= 0; i--)
        {
            FMM_UpwindSolve<${order}OrderSemiFiniteDifference, BACKWARD>(index + arrival.strides[i], tag, viscosity, arrival, band, i);
            FMM_UpwindSolve<${order}OrderSemiFiniteDifference, FORWARD>(index - arrival.strides[i], tag, viscosity, arrival, band, i);
        }

        it = band.begin();
    }
</%def>

/* vim: filetype=cpp
 *
 */
