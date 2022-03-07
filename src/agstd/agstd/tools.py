import numpy as np


def np_load(*args, **kwargs):
    if 'allow_pickle' in kwargs.keys():
        return np.load(*args, **kwargs)
    else:
        return np.load(*args, allow_pickle=True, **kwargs)