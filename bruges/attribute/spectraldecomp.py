# -*- coding: utf-8 -*-
"""
Spectral decomposition

:copyright: 2015 Agile Geoscience
:license: Apache 2.0
"""
import numpy as np
from bruges.attribute import spectrogram


def spectraldecomp(data,
                   f=(0.1, 0.25, 0.4),
                   window_length=32,
                   dt=1,
                   window_type='hann',
                   overlap=1,
                   normalize=False):
    """
    Uses the STFT to decompose traces into normalized spectra. Only
    frequencies defined by the user will be output. Using 3
    frequencies will work for RGB color plots.

    :param data: A 1/2D array (samples, traces) of data that will
                 be decomposed.
    :keyword f: A list of frequencies to select from the spectra.
    :keyword window_length: The length of fft window to use for
                            the STFTs. Defaults to 32. Can be
                            provided in seconds if dt is specified.
    :keyword dt: The time sample interval of the traces.
    :keyword window_type: The type of window to use for the STFT. The
                          same input as scipy.signal.get_window.
    :keyword overlap: The fraction of overlap between adjacent
                      STFT windows

    :keyword normalize: Normalize the energy in each STFT window

    :returns: an array of shape (samples, traces, f)
    """

    # Do the 1D case
    if len(data.shape) == 1:
        ntraces = 1
    else:
        ntraces = data.shape[-1]

    if overlap > 1:
        overlap = 1

    zp = 4 * window_length

    # TODO We should think about removing these for loops
    for i in range(ntraces):

        if(ntraces == 1):
            spect = spectrogram(data, window_length, dt=dt,
                                window_type=window_type, overlap=overlap,
                                normalize=normalize, zero_padding=zp)
            if(i == 0):
                output = np.zeros((spect.shape[0], len(f)))
        else:
            spect = spectrogram(data[:, i], window_length, dt=dt,
                                window_type=window_type, overlap=overlap,
                                normalize=normalize, zero_padding=zp)
            if(i == 0):
                output = np.zeros((spect.shape[0], ntraces, len(f)))

        res = ((1. / dt) / 2.) / spect.shape[-1]

        # TODO again, we should think about removing this loop
        for j in range(len(f)):

            index = int(f[j] / res)

            if(ntraces == 1):
                output[:, j] = spect[:, index]
            else:
                output[:, i, j] = spect[:, index]

    return(output)
