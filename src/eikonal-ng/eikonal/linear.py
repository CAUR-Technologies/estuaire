#
# @Author : Jean-Pascal Mercier <jean-pascal.mercier@agsis.com>
#
# @Copyright (C) 2010 Jean-Pascal Mercier
#
# All rights reserved.
#
__doc__ = """
"""

import numpy as np
import scipy as sc

import scipy.sparse


stencil_D2D = np.array([[[ 1,  0, -1],
                        [ 0,  0,  0],
                        [ 1,  0, -1]]])

stencil_DD2D = np.array([[[ 0,  1,  0],
                         [ 1,  -4,  1],
                         [ 0,  1,  0]]])

stencil_D3D = np.array([[[ 1,  0,  1],
                         [ 0,  0,  0],
                         [ 1,  0,  0]],
                        [[ 0,  0,  0],
                         [ 0,  0,  0],
                         [ 0,  0,  0]],
                        [[ 0,  0,  0],
                         [ 0,  0,  0],
                         [ 0,  0,  0]],
                        [[ 0,  0, -1],
                         [ 0,  0,  0],
                         [-1,  0, -1]]])

stencil_DD3D = np.array([[[ 0,  0,  0],
                          [ 0,  1,  0],
                          [ 0,  0,  0]],
                         [[ 0,  1,  0],
                          [ 1,  -6,  1],
                          [ 0,  1,  0]],
                         [[ 0,  0,  0],
                          [ 0,  1,  0],
                          [ 0,  0,  0]]])


stencils_dict = { 2 : [stencil_D2D, stencil_DD2D],
                  3 : [stencil_D3D, stencil_DD3D] }



def linearFilter2(stencil, shape):
    return linearFilter3(np.array([stencil]), (1,) + shape)


def linearFilter3(stencil, shape):
    size = np.prod(shape)
    Xmin = (stencil.shape[2] - 1) // 2
    Ymin = (stencil.shape[1] - 1) // 2
    Zmin = (stencil.shape[0] - 1) // 2


    dias = []
    poss = []

    for k, z in enumerate(range(-Zmin, Zmin + 1)):
        for j, y in enumerate(range(-Ymin, Ymin + 1)):
            for i, x in enumerate(range(-Xmin, Xmin + 1)):
                pos = z * shape[-1] * shape[-2] + y * shape[-1] + x
                dias.append([0.0] * Xmin + [stencil[k, j, i]] * (shape[-1] - 2 * Xmin) + [0.0] * Xmin)
                dias[-1] = [0.0] * (shape[2] * Ymin) + dias[-1] * (shape[1] - 2 * Ymin) + [0.0] * (shape[2] * Ymin)
                dias[-1] = [0.0] * (shape[2] * shape[1]) * Zmin + dias[-1] * (shape[0] - 2 * Zmin) + [0.0] * (shape[2] * shape[1]) * Zmin
                if pos != 0:
                    dias[-1] = dias[-1][-pos:] + dias[-1][:-pos]
                poss.append(pos)
    return sc.sparse.dia_matrix((dias, poss), (size, size), dtype = 'float')


def linearDerivativeOP(shape):
    stencils = stencils_dict[len(shape)]
    shape = shape if len(shape) == 3 else (1,) + shape
    return [linearFilter3(stencil, shape) for stencil in stencils]


def laplaceOP(X, Y, Z):
    size = X * Y * Z
    center = ([1] + [2] * (X - 2) + [1]) * Y * Z
    right = ([0] + [-1] * (X - 1)) * Y * Z
    left = ([-1] * (X - 1) + [0]) * Y * Z

    d2x = sc.sparse.dia_matrix(([left, center, right], [-1.0, 0.0, 1.0]),
                              (size, size),
                              dtype = 'float')


    if Y == 1:
        d2y = sc.sparse.eye(size, size) * 0
    else:
        center = ([1] * X + [2] * X * (Y - 2) + [1] * X)* Z
        right = ([0] * X + [-1] * X * (Y - 1)) * Z
        left = ([-1] * X * (Y - 1) + [0] * X) * Z

        d2y = sc.sparse.dia_matrix(([left, center, right], [-X, 0.0, +X]),
                                  (size, size),
                                  dtype = 'float')

    if Z == 1:
        d2z = sc.sparse.eye(size, size) * 0
    else:
        center = ([1] * X * Y + [2] * X * Y * (Z - 2) + [1] * X * Y)
        right = ([0] * X * Y + [-1] * X * Y * (Z - 1))
        left = [-1] * X * Y * (Z - 1) + [0] * X * Y

        d2z = sc.sparse.dia_matrix(([left, center, right], [-X, 0.0, +X]),
                                  (size, size),
                                  dtype = 'float')

    return d2x + d2y + d2z

def linearXDerivativeOP(X, Y, Z):
    """
    This method calculate the linear operator for the first and second
    derivatives of a 3D grid in the X direction. (Z = 1) for the 2D one
    """
    size = X * Y * Z

    center = ([-1] * (X - 1) + [-1]) * Y * Z
    right = ([0] + [1] * (X - 1)) * Y * Z
    left = ([0] * (X - 2) + [1, 0] ) * Y * Z
    #center = ([0] * (X -1) + [0]) * Y * Z + [0] * overhead
    #right = ([0, 0] + [1] * (X - 2)) * Y * Z + [0] * overhead
    #left = ([-1] * (X - 2) + [0, 0]) * Y * Z + [0] * overhead
    d = sc.sparse.dia_matrix(([left, center, right], [-1, 0, 1]),
                             (size, size),
                             dtype ='float')

    center = ([0] + [-2] * (X - 2) + [0]) * Y * Z 
    right = ([0, 0] + [1] * (X - 2)) * Y * Z
    left = ([1] * (X - 2) + [0, 0]) * Y * Z

    d2 = sc.sparse.dia_matrix(([left, center, right], [-1.0, 0.0, 1.0]),
                              (size, size),
                              dtype = 'float')

    return d,  d2

def linearYDerivativeOP(X, Y, Z):
    """
    This method calculate the linear operator for the first and second
    partial derivatives of a 3D grid in the Y dii=rection. (Z = 1) for
    the 2D one.
    """

    size = X * Y * Z
    center = ([-1] * X * (Y - 1) + [-1] * X) * Z
    right = ([0] * X + [1] * X * (Y - 1)) * Z
    left = ([0] * X * (Y - 2) + [1] * X + [0] * X) * Z
    #center = ([0] * X * (Y -1) + [0] * X) * Z + [0] * overhead
    #right = ([0, 0] * X + [1] * X * (Y - 2)) * Z + [0] * overhead
    #left = ([-1] * X * (Y - 2) + [0, 0] * X) * Z + [0] * overhead

    d = sc.sparse.dia_matrix(([left, center, right], [-X, 0, X]),
                             (size, size),
                             dtype = 'float')

    center = ([0] * X + [-2] * X * (Y - 2) + [0] * X) * Z
    right = ([0, 0] * X  + [1] * X * (Y - 2)) * Z
    left = ([1] * X * (Y - 2) + [0, 0] * X) * Z

    dd = sc.sparse.dia_matrix(([left, center, right], [-X, 0.0, +X]),
                              (size, size),
                              dtype = 'float')
    return d, dd 

def linearZDerivativeOP(X, Y, Z):
    """
    This method calculate the linear operator for the first and second
    partial derivatives of a 3D grid in the Z direction. (Z = 1) for
    the 2D one.
    """

    size = X * Y * Z

    center = ([-1] * X * Y ) * (Z - 1) + [-1] * X * Y
    right = [0] * X * Y + ([1] * X * Y) * (Z - 1)
    left = [0] * X * Y * (Z - 2) + [1] * X * Y + [0] * X * Y

    d = sc.sparse.dia_matrix(([left, center, right], [-Y * X, 0, Y * X]),
                             (size, size),
                             dtype = 'float')

    center = [0] * X * Y + [-2] * X * Y * (Z - 2) + [0] * X * Y
    right = [0, 0] * X * Y  + [1] * X * Y * (Z - 2)
    left = [1] * X * Y * (Z - 2) + [0, 0] * X * Y

    dd = sc.sparse.dia_matrix(([left, center, right], [-X * Y, 0.0, +X * Y]),
                              (size, size),
                              dtype = 'float')

    return d, dd

def uniform_weight(residual):
    return sc.sparse.eye(residual.size, residual.size)

class GaussianWeight(object):
    def __init__(self, sigma = 0):
        self.sigma = sigma

    def __call__(self, residual):
        var = np.var(residual)
        weights = np.exp(-0.5 * residual ** 2 / var)
        threshold = np.exp(-0.5 *  self.sigma ** 2)
        tweights = np.where(weights > threshold, threshold, weights)
        return sc.sparse.dia_matrix((tweights / threshold, (0, )), shape = (residual.size, residual.size))

class BoxWeight(object):
    def __init__(self, sigma = 0):
        self.sigma = sigma

    def __call__(self, residual):
        var = np.var(residual)
        weights = np.exp(-0.5 * residual ** 2 / var)
        threshold = np.exp(-0.5 *  self.sigma ** 2)
        tweights = np.where(weights > threshold, threshold, 0)
        return sc.sparse.dia_matrix((tweights / threshold, (0, )), shape = (residual.size, residual.size))


stencils_dict = { 2 : [stencil_D2D, stencil_DD2D],
                  3 : [stencil_D3D, stencil_DD3D] }


def linearFilter2(stencil, shape):
    return linearFilter3(np.array([stencil]), (1,) + shape)


def linearFilter3(stencil, shape):
    size = np.prod(shape)
    Xmin = int((stencil.shape[2] - 1) // 2)
    Ymin = int((stencil.shape[1] - 1) // 2)
    Zmin = int((stencil.shape[0] - 1) // 2)


    dias = []
    poss = []

    for k, z in enumerate(range(-Zmin, Zmin + 1)):
        for j, y in enumerate(range(-Ymin, Ymin + 1)):
            for i, x in enumerate(range(-Xmin, Xmin + 1)):
                pos = z * shape[-1] * shape[-2] + y * shape[-1] + x
                dias.append([0.0] * Xmin + [stencil[k, j, i]] * (shape[-1] - 2 * Xmin) + [0.0] * Xmin)
                dias[-1] = [0.0] * (shape[2] * Ymin) + dias[-1] * (shape[1] - 2 * Ymin) + [0.0] * (shape[2] * Ymin)
                dias[-1] = [0.0] * (shape[2] * shape[1]) * Zmin + dias[-1] * (shape[0] - 2 * Zmin) + [0.0] * (shape[2] * shape[1]) * Zmin
                if pos != 0:
                    dias[-1] = dias[-1][-pos:] + dias[-1][:-pos]
                poss.append(pos)
    return sc.sparse.dia_matrix((dias, poss), (size, size), dtype = 'float')


def linearDerivativeOP(shape):
    stencils = stencils_dict[len(shape)]
    shape = shape if len(shape) == 3 else (1,) + shape
    dx, dxx = linearXDerivativeOP(*shape)
    dy, dyy = linearYDerivativeOP(*shape)
    dz, dzz = linearZDerivativeOP(*shape)
    return [dx + dy + dz, dxx + dyy + dzz + 2 * (dx * dz + dx * dy +  dy * dz)]
    #return [linearXDerivative(X,
    #return [linearFilter3(stencil, shape) for stencil in stencils]



class CGObj(object):
    def __init__(self,  A, residual, m0, dT, prior):
        self.A = A
        self.sT = dT + prior
        self.q = m0 * dT - (residual * A)
        self.m0 = m0

    def __call__(self, X):
        value = np.dot(X, self.fhess_p(X, X)) + np.dot(self.q, X)
        return value

    def fprime(self, X):
        j = self.fhess_p(X, X) + self.q
        return j

    def fhess_p(self, X, p):
        h = self.A.T * (self.A * p) + self.sT * p
        return h


class CGInversion(object):
    def __call__(self, m, A, residual, m0, dT, prior, maxiter = None, gtol = 1e-9, batch = 10):
        self.jnorms = jnorms = []
        q = m0 * dT - (residual * A)

        # Gradient
        sT = dT + prior
        ri = -(q)
        pi = ri

        AT = A.T.tocsr()

        if np.sum(ri) == 0:
            yield m, 0
            return

        for i in range(maxiter // batch):
            for j in range(batch):
                # This line calculate the first half of the hessian
                ridot = np.dot(ri, ri)

                Ap = sT * pi + AT * (A * pi)
                ai = ridot / np.dot(pi, Ap)
                m = m + ai * pi

                ri = ri - ai * Ap


                # Stopping criterion
                jnorms.append(np.sum(np.abs(ri)))
                #jnorm2 = np.dot(r_i1, r_i1)

                bi = np.dot(ri, ri) / ridot
                pi = ri + bi * pi

            yield m
            if jnorms[-1] < gtol:
                    return
