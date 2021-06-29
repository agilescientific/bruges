import numpy as np
from scipy.ndimage import zoom
from scipy.interpolate import interp1d


def reconcile(*arrays, order=0):
    """
    Make sure 1D arrays are the same length. If not, stretch them to match
    the longest.

    Args:
        arrays (ndarray): The input arrays.
        order (int): The order of the interpolation, passed to
            scipy.ndimage.zoom. Suggestion: 0 for integers and 1 for floats.

    Returns:
        tuple of ndarrays. The reconciled arrays — all of them are now the
            same length.

    Example:
        >>> a = np.array([2, 6, 7, 7, 3])
        >>> b = np.array([3, 7, 3])
        >>> reconcile(a, b, order=0)
        (array([2, 6, 7, 7, 3]), array([3, 7, 7, 3, 3]))
    """
    maxl = max(len(arr) for arr in arrays)
    out = []
    for arr in arrays:
        if len(arr) < maxl:
            out.append(zoom(arr, zoom=maxl/len(arr), order=order))
        else:
            out.append(arr)
    return tuple(out)


def interpolate(*arrays, num=50, dists=None, kind='linear'):
    """
    Linear interpolation between 1D arrays of the same length.

    Args:
        arrays (ndarray): The 1D arrays to interpolate. All must be the same
            length. You can use the `reconcile()` function to produce them.
        num (int): The number of steps to take, so will be the width (number
            of cols) of the output array.
        dists (array-like): A list or tuple or array of the distances (any
            units) between the arrays in the real world.
        kind (str): Will be passed to scipy.interpolate.interp1d, which does
            the lateral interpolation between samples.

    Returns:
        ndarray. The result, with `num` columns. The number of rows is the
            same as the number of samples in the input arrays.

    Example:
        >>> a = np.array([2, 6, 7, 7, 3])
        >>> b = np.array([3, 7, 7, 3, 3])
        >>> interp = interpolate(a, b, num=10)
        >>> interp.shape
        (5, 10)
    """
    intervals = len(arrays) - 1
    if dists is None:
        dists = intervals * [num / intervals]
    x = np.hstack([[0], np.cumsum(dists)])
    f = interp1d(x, np.stack(arrays), axis=0, kind=kind)
    return f(np.linspace(x[0], x[-1], num=num)).T


def unreconcile(arr, sizes, dists=None, order=0):
    """
    Opposite of reconcile. Restores the various profiles (the reference arrays,
    e.g. wells) to their original lengths.

    Args:
        sizes (int): The relative lengths of the profiles in the array.
            Default returns the input array.
        dists (array-like): The relative distances between the profiles in
            the array. Sum used to calculate the output width in pixels if
            the width argument is None. If not given, the distances are
            assumed to be equal.
        order (int): The order of the spline interpolation, from 0 to 3. The
            default is 0, which gives nearest neighbour interpolation. 1 gives
            linear interpolation, etc. Use 0 for ints and 1-3 for floats.

    Returns:
        ndarray. The resulting ndarray. The array contains NaNs where there
            is no interpolated data.
    """
    if np.all(sizes[0] == np.array(sizes)):
        # Nothing to do.
        return arr

    intervals = len(sizes) - 1

    if dists is None:
        eq = arr.shape[-1] // intervals
        dists = [eq] * intervals
    assert len(dists) == intervals

    maxlen = int(np.ceil(max(sizes) * arr.shape[0]))

    dist_ = np.cumsum(dists)
    idx = arr.shape[-1] * dist_ / max(dist_)
    chunks = np.split(arr, idx[:-1].astype(int), axis=-1)

    zoomed = []
    for left, right, chunk in zip(sizes[:-1], sizes[1:], chunks):
        zooms = np.linspace(left, right, chunk.shape[-1]+1)
        for z, col in zip(zooms, chunk.T):
            new_ = zoom(col, zoom=z, order=order, mode='nearest')
            pad_width = maxlen - new_.size
            new = np.pad(new_,
                         pad_width=(0, pad_width),
                         mode='constant',
                         constant_values=np.nan,
                         )
            zoomed.append(new)

    return np.array(zoomed).T


def panel(*arrays, num=50, dists=None, order=0, kind='linear'):
    """
    Interpolate an arbitrary collection of 1D arrays.

    Args:
        num (int): The number of steps to take, so will be the width (number
            of cols) of the output array.
        dists (array-like): The relative distances between the profiles in
            the array. Sum used to calculate the output width in pixels if
            the width argument is None. If not given, the distances are
            assumed to be equal.
        order (int): The order of the interpolation, passed to
            scipy.ndimage.zoom. Suggestion: 0 for integers and 1 for floats.
        kind (str): Will be passed to scipy.interpolate.interp1d, which does
            the lateral interpolation between samples.

    Returns:
        ndarray. The interpolated panel. Contains NaNs if sizes are
            non-uniform.
    """
    sizes = np.array([len(x) for x in arrays])
    sizes = sizes / np.max(sizes)
    rec = reconcile(*arrays)
    interp = interpolate(*rec, num=num, dists=dists, kind=kind)
    panel = unreconcile(interp, sizes=sizes, dists=dists, order=order)
    return panel
