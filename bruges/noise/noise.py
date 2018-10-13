# -*- coding: utf-8 -*-
"""
Noise.

:copyright: 2015 Agile Geoscience
:license: Apache 2.0
"""
import numpy as np

from bruges.util import rms


def noise_db(a, snr):
    """
    Takes an array of seismic amplitudes and SNR in dB.

     Args:
        a (array) : seismic amplitude array.
        snr (int): signal to noise ratio.

    Returns:  Noise array, the same shape as the input.

     Note: it does *not* return the input array with the noise added.

    """
    # Get the amplitude of the signal
    sigmean = rms(a)

    # Calculate the amp of the noise,
    # given the desired SNR
    noisemean = sigmean / 10.0**(snr/20.0)

    # Normal noise, centered on 0,
    # SD=sqrt(var), same shape as input
    noise = noisemean * np.random.normal(0.0, 1.0, a.shape)

    return noise
