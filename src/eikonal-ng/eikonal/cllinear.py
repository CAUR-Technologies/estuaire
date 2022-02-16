#
# @Author : Jean-Pascal Mercier <jean-pascal.mercier@agsis.com>
#
# @Copyright (C) 2010 Jean-Pascal Mercier
#
# All rights reserved.
#
__doc__ = """
"""

# import cycl.clmath as clmath
import numpy as np

"""

float csr_spmv(__global     float   *data,
               __global     int     *indices,
               __global     int     *indptr)
{
}
__kernel_cqiter(__global    float   *A_data
                __global    int     *A_indices,
                __global    int     *A_indptr,
                __global    float   *AT_data
                __global    int     *AT_indices,
                __global    int     *AT_indptr,
                __global    float   *q,
                __global    float   *m,
                int                 size)
{
        int index = get_global_id(0);

        float pi = -q[index];
        float ri = -q[index];



}
"""


# class GPUCGInversion(object):
#     def __init__(self, ctx, A, residual, prior):
#         self.GPU_prior  = ctx.createBufferLike(prior)
#         self.A          = clmath.CLCSRMatrix(ctx, A)
#         self.AT         = clmath.CLCSRMatrix(ctx, (A.tocsc().T))
#
#         self.tmp        = ctx.createBufferLike(q)
#         self.GPU_ri     = ctx.createBufferLike(q)
#         self.m          = clmath.LinkedArray(np.zeros(q.shape, dtype ='float32'))
#
#         self.A          = A
#         self.q          = -(residual * A)
#
#
#     def batch(self, batchsize = 200):
#         """
#         Returns a batch command
#         """
#         Ap = axpy(multiply(tmp, q, pi))
#
#     def send(self):
#         """
#         Returns a send command
#         """
#         return self.GPU_A.send() + self.GPU_AT.send()


class CGInversion(object):
    def __call__(self, m, A, residual, m0, dT, prior, maxiter = None,
                 gtol = 1e-9, batch = 10):
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

        for i in range(maxiter / batch):
            for j in range(batch):
                # This line calculate the first half of the hessian
                Ap = sT * pi + AT * (A * pi)
                ai = np.dot(ri, ri) / np.dot(pi, Ap)
                m = m + ai * pi

                r_i1 = ri - ai * Ap

                # Stopping criterion
                jnorms.append(np.sum(np.abs(r_i1)))
                #jnorm2 = np.dot(r_i1, r_i1)

                bi = np.dot(r_i1, r_i1) / np.dot(ri, ri)
                ri = r_i1
                pi = ri + bi * pi

            yield m
            if jnorms[-1] < gtol:
                    return

