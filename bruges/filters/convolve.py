# -*- coding: utf-8 -*-
"""
Various convolution algorithms.

:copyright: 2019 Agile Geoscience
:license: Apache 2.0
"""
import numpy as np
from util import convolve_many


def convolve(reflectivity, wavelet):
    """
    Kep trying to convolve traces to handle different shaped
    inputs. There is certainly a better way to do this.

    TODO
    - Find a better way to do this!
    - Handle the case of a 3D RC series with a wavelet bank
        to give a 5D result.
    """
    if wavelet.ndim == 1:
        try:
            # 1D reflectivity -> 1D synthetic.
            syn = np.convolve(reflectivity, wavelet, mode='same')

        except ValueError:
            try:
                # 2D reflectivity -> 2D synthetic.
                syn = convolve_many(reflectivity, wavelet)

            except ValueError:
                    # 3D reflectivity -> 3D synthetic.
                    syn = np.array([convolve_many(tr, wavelet) for tr in reflectivity])
    elif wavelet.ndim == 2:
        try:
            # 1D reflectivity, 2D wavelet bank -> 2D synthetic.
            syn = convolve_many(wavelet, reflectivity)
        except ValueError:
            # 2D reflectivity, 2D wavelet bank -> 3D synthetic.
            syn = np.array([convolve_many(reflectivity, w) for w in bank])
    else:
        raise NotImplementedError("Wavelets must be 1d or 2d arrays.")

    return syn
