# -*- coding: utf-8 -*-
"""
Time-depth conversion.

:copyright: 2015 Agile Geoscience
:license: Apache 2.0
"""
from collections import namedtuple

from scipy.interpolate import interp1d
import numpy as np


def __convert(data, vmodel, interval, interval_new, scale, mode, return_basis=False):
    """
    Generic function for converting between scales. Use either
    time to depth or depth to time.

        Args:
        data (ndarray): The data to convert, will work with a 1 or 2D numpy
            numpy array. array(samples,traces).
        vmodel (ndarray): P-wave interval velocity model that corresponds to
            the data. Must be the same shape as data.
       interval (float): The sample interval of the input data [s] or [m].
       interval_new (float): The sample interval of the output data [m] or [s].
       mode (str): What kind of interpolation to use, defaults to 'nearest'.
        return_basis (bool): Whether to also return the new time basis.

    Returns
        ndarray: The data resampled in the depth domain.
    """
    data = np.array(data)

    if np.ndim(data) == 1:
        ntraces = 1
        nsamps = data.size
        if np.size(interval) == 1:
            basis = np.arange(nsamps)*interval
        else:
            basis = interval
        v_avg = np.cumsum(vmodel) / (np.arange(nsamps) + 1)
    else:
        ntraces = data.shape[-1]
        nsamps = data.shape[0]
        if np.size(interval) == 1:
            tr = [(np.arange(nsamps) * interval) for i in range(ntraces)]
            basis = np.transpose(np.asarray(tr))
        else:
            basis = interval
        tr = [np.arange(nsamps) + 1 for i in range(ntraces)]
        v_avg = np.cumsum(vmodel, axis=0) / np.transpose(tr)

    new_basis = basis / v_avg
    new_basis *= scale

    if np.size(interval_new) == 1:
        new_basis_lin = np.arange(np.amin(new_basis), np.amax(new_basis), interval_new)
    else:
        new_basis_lin = interval_new

    if np.ndim(data) == 1:
        inter = interp1d(new_basis, data,
                         bounds_error=False,
                         fill_value=data[-1],
                         kind=mode)
        result = inter(new_basis_lin)
    else:
        result = np.zeros((new_basis_lin.size, ntraces))
        for i in range(ntraces):
            inter = interp1d(new_basis[:, i], data[:, i],
                             bounds_error=False,
                             fill_value=data[-1, i],
                             kind=mode)
            result[:, i] += inter(new_basis_lin)

    if return_basis:
        field_names = ['data', 'basis']
        Conversion = namedtuple('Conversion', field_names)
        return Conversion(result, new_basis_lin)
    else:
        return result


def time_to_depth(data, vmodel, dt, dz, twt=True, mode="nearest", return_z=False):
    """
    Converts data from the time domain to the depth domain given a
    velocity model.

    Args:
        data (ndarray): The data to convert, will work with a 1 or 2D numpy
            numpy array. array(samples,traces).
        vmodel (ndarray): P-wave interval velocity model that corresponds to
            the data. Must be the same shape as data.
        dt (float): The sample interval of the input data [s], or an
            array of times.
        dz (float): The sample interval of the output data [m], or an
            array of depths.
        twt (bool): Use twt travel time, defaults to true.
        mode (str): What kind of interpolation to use, defaults to 'nearest'.
        return_z (bool): Whether to also return the new time basis.

    Returns
        ndarray: The data resampled in the depth domain.
    """
    if twt:
        scale = 0.5
    else:
        scale = 1.0

    # Do conversion with inverted velocity profile (slowness).
    return __convert(data,
                     vmodel=1. / vmodel,
                     interval=dt,
                     interval_new=dz,
                     scale=scale,
                     mode=mode,
                     return_basis=return_z,
                     )


def depth_to_time(data, vmodel, dz, dt, twt=True, mode="nearest", return_t=False):
    """
    Converts data from the depth domain to the time domain given a
    velocity model.

    Args:
        data (ndarray): The data to convert, will work with a 1 or 2D numpy
            numpy array. array(samples,traces).
        vmodel (ndarray): P-wave interval velocity model that corresponds to
            the data. Must be the same shape as data.
        dz (float): The sample interval of the input data [m].
        dt (float): The sample interval of the output data [s].
        twt (bool): Use twt travel time, defaults to true.
        mode (str): What kind of interpolation to use, defaults to 'nearest'.
        return_t (bool): Whether to also return the new time basis.

    Returns:
        The data resampled in the time domain.
    """
    if twt:
        scale = 2.0
    else:
        scale = 1.0

    return __convert(data,
                     vmodel=vmodel,
                     interval=dz,
                     interval_new=dt,
                     scale=scale,
                     mode=mode,
                     return_basis=return_t,
                     )
