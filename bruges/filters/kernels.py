"""
2D kernels for image processing.

:copyright: 2015 Agile Geoscience
:license: Apache 2.0
"""
import numpy as np


def gaussian(size, size_y=None):
    """
    2D Gaussian Kernel

    Args:
        size (int): the kernel size, e.g. 5 for 5x5 (in a 2D arr). Should be
            odd, rounded up if not.
        size_y (int): similar to size. If not provided, uses size as default. 

    Returns: a Gaussian kernel.
    """
    size = int(size)
    if not size_y:
        size_y = size
    else:
        size_y = int(size_y)
    x, y = np.mgrid[-size:size+1, -size_y:size_y+1]
    g = np.exp(-(x**2/float(size)+y**2/float(size_y)))
    return g / g.sum()


gaussian_kernel = gaussian
