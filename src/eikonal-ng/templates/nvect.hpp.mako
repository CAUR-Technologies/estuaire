/**
 * @Author : Jean-Pascal Mercier <jean-pascal.mercier@agsis.com>
 *
 * @Copyright (C) 2010 Jean-Pascal Mercier
 *
 * All rights reserved.
 *
 * This file contains a simple and efficient implementation of tiny vector
 * of known size.
 *
 */



#ifndef NVECT_HPP
#define NVECT_HPP

#include <iostream>
#include <cmath>

// Vector Declaration
// TODO : Adding Standard Method for generic vect
namespace agsis
{

template <typename T, std::size_t N>
union vect {
    T               data[N];
    T& operator[](const std::size_t i);
    const T& operator[](const std::size_t i) const;
};

template <typename T>
union vect<T, 1>{
    vect() : x() {};
    vect(const T a[1]) : x(a[0]) {};

    template <typename T2>
    vect(const vect<T2, 1> &a) : x(a.x) {};
    struct {
        T       x;
    };
    T               data[1];
    T& operator[](const std::size_t i);
    const T& operator[](const std::size_t i) const;
};



template <typename T>
union vect<T, 2> {
    vect() :  y(), x() {};
    vect(const T a, const T b) : y(b), x(a) {};
    vect(const T a[2]) : y(a[0]), x(a[1]) {};

    template <typename T2>
    vect(const vect<T2, 2> &a) : y(a.y), x(a.x) {};
    struct {
        T           y;
        T           x;
    };
    T data[2];
    T& operator[](const std::size_t i);
    const T& operator[](const std::size_t i) const;

};

template <typename T>
union vect<T, 3> {
    vect() :  z(), y(), x() {};
    vect(const T a, const T b, const T c) : z(c), y(b), x(a) {};
    vect(const T a[3]) : z(a[0]), y(a[1]), x(a[2]) {};

    template <typename T2>
    vect(const vect<T2, 3> &a) : z(a.z), y(a.y), x(a.x) {};
    struct {
        T           z;
        T           y;
        T           x;
    };
    T data[3];

    T& operator[](const std::size_t i);
    const T& operator[](const std::size_t i) const;
};

template <typename T, std::size_t N>
inline T& vect<T, N>::operator[](const std::size_t i)
{
    return this->data[i];
}

template <typename T, std::size_t N>
const inline T& vect<T, N>::operator[](const std::size_t i) const
{
    return this->data[i];
}

% for n in ["1", "2", "3"]:
template <typename T>
const inline T& vect<T, ${n}>::operator[](const std::size_t i) const
{
    return this->data[i];
}
template <typename T>
inline T& vect<T, ${n}>::operator[](const std::size_t i)
{
    return this->data[i];
}
% endfor

// Gener[al/ic] Vector Operations
template <typename T, std::size_t N>
inline T norm_sq(const vect<T, N> &a)
{
    return dot(a, a);
}

template <typename T, std::size_t N>
inline double norm(const vect<T, N> &a)
{
    return sqrt(norm_sq(a));
}

// Operator +, -, /, *
% for op in ["+", "-", "/", "*"]:
// With another Vector
template <typename T, std::size_t N>
vect<T, N> operator${op}(const vect<T, N> &a, const vect<T, N> &b)
{
    vect<T, N> result;
    for (std::size_t i = 0; i < N; i++)
    {
        result[i] = a[i] ${op} b[i];
    }
    return result;
}

// With a scalar
template <typename T, std::size_t N>
inline vect<T, N> operator${op}(const vect<T, N> &a, const T& b)
{
    vect<T, N> result;
    for (std::size_t i = 0; i < N; i++)
    {
        result[i] = a[i] ${op} b;
    }
    return result;
}

// With a scalar (Reversed)
template <typename T, std::size_t N>
inline vect<T, N> operator${op}(const T &a, const vect<T, N> &b)
{
    vect<T, N> result;
    for (std::size_t i = 0; i < N; i++)
    {
        result[i] = a ${op} b[i];
    }
    return result;
}

// No temporary assignment operator (Vector)
template <typename T, std::size_t N>
vect<T, N>& operator${op}=(vect<T, N>& a, const vect<T, N>& b)
{
    for (std::size_t i = 0; i < N; i++)
    {
        a[i] ${op}= b[i];
    }
    return a;
}

// No temporary assignment operator (Scalar)
template <typename T, std::size_t N>
vect<T, N>& operator${op}=(vect<T, N>& a, const T& b)
{
    for (std::size_t i = 0; i < N; i++)
    {
        a[i] ${op}= b;
    }
    return a;
}

% endfor

% for op in ["+", "-"]:
template <typename T, std::size_t N>
inline vect<T, N> operator${op}(const vect<T, N> &a)
{
    vect<T, N> result;
    for (std::size_t i = 0; i < N; i++)
    {
        result[i] = ${op} a[i];
    }
    return result;
}
% endfor

// Floor and Ceil Operator for float and double
% for dtype in ["double", "float"]:
% for func in ["floor", "ceil"]:
template <std::size_t N>
inline vect<${dtype}, N> ${func}(const vect<${dtype}, N> &a)
{
    vect<${dtype}, N> result;
    for (std::size_t i = 0; i < N; i++)
    {
        result[i] = std::${func}(a[i]);
    }
    return result;
}
% endfor
% endfor

// Comparison operators
% for op in [">", "<", ">=", "<=", "||", "&&", ">>", "<<"]:
template <typename T, std::size_t N>
inline vect<bool, N> operator${op}(const vect<T, N> &a, const vect<T, N> &b)
{
    vect<bool, N> result;
    for (std::size_t i = 0; i < N; i++)
    {
        result[i] = (a.data[i] ${op} b[i]);
    }
    return result;
}

template <typename T, std::size_t N>
inline vect<bool, N> operator${op}(const vect<T, N> &a, const T&b)
{
    vect<bool, N> result;
    for (std::size_t i = 0; i < N; i++)
    {
        result[i] = (a[i] ${op} b);
    }
    return result;
}

template <typename T, std::size_t N>
inline vect<bool, N> operator${op}(const T&a, const vect<T, N> &b)
{
    vect<bool, N> result;
    for (std::size_t i = 0; i < N; i++)
    {
        result[i] = (b ${op} a[i]);
    }
    return result;
}
% endfor
% for op in ["==", "!="]:
template <typename T, std::size_t N>
inline bool operator${op}(const vect<T, N> &a, const vect<T, N> &b)
{
    bool result;
    for (std::size_t i = 0; i < N; i++)
    {
        result = result && (a.data[i] ${op} b[i]);
    }
    return result;
}
% endfor

// Min max operator
% for op in ["max", "min"]:
template <typename T, std::size_t N>
inline T ${op}(const vect<T, N> &a)
{
    T result = a[0];
    for (std::size_t i = 1; i < N; i++)
    {
        result = std::${op}(result, a[i]);
    }
    return result;
}
% endfor

// Sum over the vector
template <typename T, std::size_t N>
inline T sum(const vect<T, N> &a)
{
    T total = a[0];
    for (std::size_t i = 1; i < N; i++)
    {
        total += a[i];
    }
    return total;
}

// Product of elements of the vector
template <typename T, std::size_t N>
inline T prod(const vect<T, N> &a)
{
    T total = a[0];
    for (std::size_t i = 1; i < N; i++)
    {
        total *= a.data[i];
    }
    return total;
}

// Dot product between 2 vectors
template <typename T, std::size_t N>
inline T dot(const vect<T, N> &a, const vect<T, N> &b)
{
    return sum(a * b);
}

// Verify if any of the element of the vector is True
template <typename T, std::size_t N>
inline bool any(const vect<T, N> &a)
{
    bool result = a[0];
    for (std::size_t i = 1; i < N; i++)
    {
        result |= a[i];
    }
    return result;
}

// Verify that all the elements of the vector are True
template <typename T, std::size_t N>
inline bool all(const vect<T, N> &a)
{
    bool result = a[0];

    for (std::size_t i = 1; i < N; i++)
    {
        result &= a[0];
    }
    return result;
}

// Clamp the vector to boundari
template <typename T, typename T2, std::size_t N>
inline vect<T, N> clamp(const vect<T, N> &a, const vect<T2, N> &b)
{
    vect<T, N> result;
    for (std::size_t i = 0; i < N; i++)
    {
        result[i] = a[i] >= b[i] ? b[i] - 1 : a[i];
    }
    return result;
}



// Generic (and pretty) Input / Output
template <typename T, std::size_t N>
inline std::ostream &operator<<( std::ostream &out, const vect<T, N> &a )
{
    out << "[" << a[0];
    for (std::size_t i = 1; i < N; i++)
    {
        out << ", " << a[i];
    }
    return out << "]";
}

} // namespace agsis


#endif // NVECT_HPP

/* vim: filetype=cpp
 *
 */

