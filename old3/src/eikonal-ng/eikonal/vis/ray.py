from enthought.tvtk.tvtk import tvtk
import numpy as np


def ray_source(rays):
    data = tvtk.PolyData()
    data.points = np.concatenate(rays, axis = 0)

    length = np.cumsum([0] + [pos.shape[0] for pos in rays[:-1]])
    ranges = [np.arange(pos.shape[0]) for pos in rays]
    data.lines = [(l + r).tolist() for l, r in zip(length, ranges)]

    return data
