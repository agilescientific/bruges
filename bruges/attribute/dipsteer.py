# -*- coding: utf-8 -*-
"""
A dip attribute, probably most useful for guiding other attributes.

:copyright: 2015 Agile Geoscience
:license: Apache 2.0
"""
from collections import namedtuple

import numpy as np
from bruges.attribute import energy


def dipsteer(data,
             window_length,
             stepout,
             maxlag,
             overlap=1,
             dt=1,
             return_correlation=False):
    """
    Calculates a dip field by finding the maximum correlation between
    adjacent traces.

    :param data (ndarray): A 2D seismic section (samples,traces) used to
        calculate dip.
    :param window_length (float): The length [in ms] of the window to use.
    :param stepout (int): The number of traces on either side of each point
        to average when calculating the dip.
    :param maxlag (float): The maximum amount time lag to use when correlating
        the traces.
    :keyword overlap (float): The fractional overlap for each window. A value
        of 0 uses no redudant data, a value of 1 slides the dip correlator one
        sample at a time. Defaults to 1.
    :keyword dt (float): The time sample interval in ms.
    :keyword return_correlation (bool): Whether to return the correlation
        coefficients. If you choose True, you'll get a tuple, not an ndarray.
    :returns: a dip field [samples/trace] of the same shape as the input data
        (and optionally correlation coefficients, in which case you'll get a
        tuple of ndarrays back).
    """
    maxlag = int(maxlag)
    dip = np.zeros(data.shape)
    crcf = np.zeros(data.shape)

    window_length = int(np.floor(window_length / dt))

    # Force the window length to be odd for index tracking.
    if not (window_length % 2):
        window_length += 1

    # Define time windows.
    if overlap == 1:
        stride = 1
    else:
        stride = int(window_length * (1 - overlap))
    n_windows = np.ceil((data.shape[0] - window_length) / stride) + 1

    # Normalize each trace to the same RMS energy.
    norm_factor = np.sqrt(np.abs(energy(data, window_length)))
    norm_data = data / (norm_factor + 1e-9)  # To avoid div0 error.

    # Replace the 0/0 with 0.
    norm_data = np.nan_to_num(norm_data)

    # Mid point in the data which corresponds to zero dip.
    zero_dip = (np.floor(window_length / 2.0) + maxlag)

    s = stepout + 1

    # Loop over each trace we can do a full calculation for.
    for i in np.arange(s, data.shape[-1] - s):

        i = int(i)

        # Loop over each time window.
        for j in np.arange(0, n_windows):

            start = int((j * stride) + (maxlag))
            end = start + window_length

            # Don't compute last samples if we don't have a full window.
            if (end > (norm_data.shape[0]-maxlag)):
                break

            kernel = norm_data[start: end, i]

            dips_j, crcf_j = 0, 0

            # Correlate with adjacent traces.
            for k in np.arange(1, s):

                k = int(k)

                # Do the trace on the right.
                r_trace = norm_data[start - (k*maxlag): end + (k*maxlag), i+k]

                cor_r = np.correlate(kernel, r_trace, mode='same')

                if (np.amax(cor_r) < .1):
                    dip_r = 0
                else:
                    dip_r = (np.argmax(cor_r) - zero_dip) / k

                # Do the left trace.
                l_trace = norm_data[start - (k*maxlag): end + (k*maxlag), i-k]

                cor_l = np.correlate(kernel, l_trace, mode='same')

                if (np.amax(cor_l) < .1):
                    dip_l = 0
                else:
                    dip_l = -(np.argmax(cor_l) - zero_dip) / k

                dips_j += dip_r + dip_l
                crcf_j += np.argmin(cor_l) + np.argmin(cor_r)

            # Average the result
            dips_j /= (2. * stepout)
            crcf_j /= (2. * stepout)

            # Update the output
            dip[start: start+stride, i] = dips_j
            crcf[start: start+stride, i] = crcf_j

    if return_correlation:
        DipSteer = namedtuple('DipSteer', ['dip', 'correlation_coeff'])
        return DipSteer(dip, crcf)
    else:
        return dip
