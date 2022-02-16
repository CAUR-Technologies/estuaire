#
# @Author : Jean-Pascal Mercier <jean-pascal.mercier@agsis.com>
#
# @Copyright (C) 2010 Jean-Pascal Mercier
#
# All rights reserved.
#
__doc__ = """
"""

import cycl
import numpy as np
import scipy as sc

eikonal_kernel = \
"""

inline float4 quadratic_square(const float a, const float b, const float c)
{
    float4 result = {a * a, 2 * a * b, b * b + c, 0};
    return result;
}

inline float half_quadratic_solver(float4 a)
{
    float root = sqrt(a.y * a.y - 4 * a.x * a.z);
    return (-a.y + root) / (2 * a.x);
}



__kernel void eikonal3(
                       __global float * result,
                       const __global float * arrival,
                       const __global float * viscosity,
                       const int stride1)
{
    int index = 2 + get_global_id(0) + (get_global_id(1) + 2) * stride1;
    float solution = arrival[index];

    const float sqvis = viscosity[index] * viscosity[index];

    const float norm = 1.5;

    const float DX_p = -2 * arrival[index + 1] + 0.5 * arrival[index + 2];
    const float DX_m = -2 * arrival[index - 1] + 0.5 * arrival[index - 2];
    const float DY_p = -2 * arrival[index + stride1] + 0.5 * arrival[index + 2 * stride1];
    const float DY_m = -2 * arrival[index - stride1] + 0.5 * arrival[index - 2 * stride1];

    //const float norm = 1;

    //const float DX_p = - arrival[index + 1];
    //const float DX_m = - arrival[index - 1];
    //const float DY_p = - arrival[index + stride1];
    //const float DY_m = - arrival[index - stride1];
    //
    solution = min(solution, arrival[index - 1] + viscosity[index]);
    solution = min(solution, arrival[index + 1] + viscosity[index]);
    solution = min(solution, arrival[index - stride1] + viscosity[index]);
    solution = min(solution, arrival[index + stride1] + viscosity[index]);

    //solution = min(solution, half_quadratic_solver(quadratic_square(norm, DX_p, -sqvis)));
    //solution = min(solution, half_quadratic_solver(quadratic_square(norm, DX_m, -sqvis)));

    //solution = min(solution, half_quadratic_solver(quadratic_square(norm, DY_p, -sqvis)));
    //solution = min(solution, half_quadratic_solver(quadratic_square(norm, DY_m, -sqvis)));

    result[index] = solution;

}

"""

shape = (1028, 1028)

platform = cycl.getPlatforms()[0]
device = platform.getDevices()[0]

context = platform.createContext([device])

program = context.createProgramWithSource(eikonal_kernel)
try:
    program.build()
except cycl.CLError, e:
    print program.getBuildLog(device)

kernel = program.createKernel("eikonal3")
kernel.parameters = (cycl.parameter_type.MEM_TYPE,
                     cycl.parameter_type.MEM_TYPE,
                     cycl.parameter_type.MEM_TYPE,
                     cycl.parameter_type.INT_TYPE)

queue = context.createCommandQueue(device)

cpu_viscosity = np.ones(shape, dtype = 'float32')
cpu_arrival = np.empty(shape, dtype = 'float32')
cpu_arrival.fill(np.inf)
cpu_arrival[127:129, 127:129] = 0.707

gpu_viscosity = context.createBufferLike(cpu_viscosity)
gpu_arrival = context.createBufferLike(cpu_arrival)
gpu_result = context.createBufferLike(cpu_arrival)


send_arrival = cycl.CLWriteBufferNDArray(gpu_arrival, cpu_arrival)
send_result = cycl.CLWriteBufferNDArray(gpu_result, cpu_arrival)
send_viscosity = cycl.CLWriteBufferNDArray(gpu_viscosity, cpu_viscosity)
retrieve_arrival = cycl.CLReadBufferNDArray(cpu_arrival, gpu_arrival)

eiko_cmd = cycl.CLNDRangeKernel(kernel, global_work_size = (shape[0] - 4, shape[1] - 4, 1), local_work_size = (256, 1,1))

queue.enqueue(send_arrival)
queue.enqueue(send_result)
queue.enqueue(send_viscosity)


kernel.setArgs(gpu_result, gpu_arrival, gpu_viscosity, shape[1])
def eikonal():
    for i in range(1000):
        kernel.setArgs(gpu_result, gpu_arrival, gpu_viscosity, shape[1])
        queue.enqueue(eiko_cmd)
        queue.finish()
        kernel.setArgs(gpu_arrival,gpu_result, gpu_viscosity, shape[1])
        queue.enqueue(eiko_cmd)
    queue.finish()

eikonal()


queue.enqueue(retrieve_arrival)
queue.finish()

cpu_arrival = np.where(np.isinf(cpu_arrival), 0.707, cpu_arrival)
sc.misc.imshow(cpu_arrival)



