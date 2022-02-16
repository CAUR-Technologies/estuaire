import numpy as np
import cycl
import cycl.clmath as clmath

import scipy as sc
import scipy.sparse

import eikonal.cllinear as linear


p = cycl.getPlatforms()[0]
d = p.getDevices()[0]
c = p.createContext([d])
q = c.createCommandQueue(d)

size = 64 * 64 * 64

sp = sc.sparse.dia_matrix((([1] * size, [1] * size, [1] * size), (-1, 0, 1)),
                          shape=(size, size), dtype='float32').tocsr()
gpu_sp = clmath.CLCSRMatrix(c, sp)

a = clmath.GPULinkedArray(c, np.ones(size, dtype='float32'))
b = clmath.GPULinkedArray(c, np.ones(size, dtype='float32'))

q.enqueue(a.send() + gpu_sp.send())
q.enqueue(clmath.spmv(gpu_sp, b.GPU, a.GPU) + clmath.dot(a.GPU, b.GPU, b.GPU))
q.enqueue(b.retrieve() + a.retrieve())
q.enqueue(clmath.multiply(b.GPU, b.GPU, b.GPU))
q.enqueue(b.retrieve())

q.finish()

