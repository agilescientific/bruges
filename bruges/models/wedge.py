from collections import namedtuple

import numpy as np
import scipy.ndimage as sn
import scipy.ndimage.morphology as morph

from bruges.util import sigmoid, root, power


def pad_func(before, after):
    """
    Padding function. Operates on vector *in place*, per the np.pad
    documentation.
    """
    def pad_with(x, pad_width, iaxis, kwargs):
        x[:pad_width[0]] = before[-pad_width[0]:]
        x[-pad_width[1]:] = after[:pad_width[1]]
        return
    return pad_with


def get_strat(strat,
              thickness,
              kind='nearest',
              position=1,
              wedge=None,
              zoom_mode='nearest'
              ):
    """
    Take a 'stratigraphy' (either an int, a tuple of ints, or a list-like of
    floats) and expand or compress it to the required thickness.

    `kind` can be 'nearest', 'linear', 'quadratic', or 'cubic'.
    """
    orders = {'nearest': 0, 'linear': 1, 'quadratic': 2, 'cubic': 3}
    order = orders.get(kind, 0)

    if isinstance(strat, int) and (order == 0):
        out = np.repeat([strat], thickness)
    elif isinstance(strat, float) and (order == 0):
        out = np.repeat([strat], thickness)
    elif isinstance(strat, tuple) and (order == 0):
        out = np.repeat(strat, int(round(thickness/len(strat))))
    else:
        if position == 0:
            wedge_zoom = wedge[1]/len(wedge[0])
            strat = strat[-int(thickness/wedge_zoom):]
        elif position == -1:
            wedge_zoom = wedge[1]/len(wedge[0])
            strat = strat[:int(thickness/wedge_zoom)]
        zoom = thickness / len(strat)
        out = sn.zoom(strat, zoom=zoom, order=order, mode=zoom_mode)

    # Guarantee correct length by adjusting bottom layer.
    missing = int(np.ceil(thickness - out.size))
    if out.size > 0 and missing > 0:
        out = np.pad(out, [0, missing], mode='edge')
    elif out.size > 0 and missing < 0:
        out = out[:missing]

    return out


def get_conforming(strat, thickness, conformance):
    """
    Function to deal with top and bottom conforming wedges.
    """
    thickness = int(np.ceil(thickness))
    if thickness == 0:
        return np.array([])
    if strat.size == thickness:
        return strat
    elif strat.size > thickness:
        if conformance == 'top':
            return strat[:thickness]
        else:
            return strat[-thickness:]
    else:
        if conformance == 'top':
            return np.pad(strat, [0, thickness-strat.size], mode='wrap')
        else:
            return np.pad(strat, [thickness-strat.size, 0], mode='wrap')
    return


def get_subwedges(target, breadth):
    """
    For a binary target (the reference trace of a wedge),
    create the range of net:gross subwedges. We do this
    with binary morphologies in the following way:

    - Erode the 'secondary' component (whatever appears
      second in the target array) one step at a time,
      until there is nothing left and the resulting
      trace contains only the primary.
    - Dilate the secondary component until there is
      nothing left of the primary and the resulting trace
      contains only the secondary.
    - Arrange the traces in order, starting with all
      primary, and ending in all secondary. The target
      trace will be somewhere in between, but not
      necessarily in the middle.

    Returns a 2D array with one target wedge trace per
    section.

    Args:
        target (array): A 1D array length N, the 'reference'
            trace for the wedge. The reference trace has
            thickness '1' in the wedge model. This trace must
            be 'binary' — i.e. it must contain exactly 2
            unique values.
        breadth (int): How many new reference traces should
            be in the output.

    Returns:
        tuple (ndarray, ndarray, int). The ndarray has shape
            N x breadth. It represents one target wedge trace
            per section in 'breadth'. The integer is the
            position of the target trace in the ndarray's
            second dimension.
    """
    try:
        components = a, b = np.unique(target)
    except ValueError:
        raise ValueError("Must be a binary (2-component) wedge.")

    out = [target]

    temp = target.copy()
    while b in temp:
        ero = morph.binary_erosion(temp == b)
        temp = components[ero.astype(int)]
        out.append(temp)

    out = out[::-1]
    ref = len(out)

    temp = target.copy()
    while a in temp:
        dil = morph.binary_dilation(temp == b)
        temp = components[dil.astype(int)]
        out.append(temp)

    arr_ = np.array(out).T
    h, w = arr_.shape
    arr = sn.zoom(arr_, zoom=(1, breadth/w), mode='nearest')

    ng = np.divide(np.sum(arr == a, axis=0), target.size)

    return arr, ng, ref * breadth / w


def wedge(depth=(30, 40, 30),
          width=(10, 80, 10),
          breadth=None,
          strat=(0, 1, 2),
          thickness=(0.0, 1.0),
          mode='linear',
          conformance='both',
          ):
    """
    Generate a wedge model.

    Args:
        depth (int or tuple): The vertical size of the model. If a 3-tuple,
            then each element corresponds to a layer. If an integer, then each
            layer of the model will be 1/3 of the thickness. Note that if the
            'right' wedge thickness is more than 1, then the total thickness
            will be greater than this value.
        width (int or tuple): The width of the model. If a 3-tuple, then each
            element corresponds to a 'zone' (left, middle, right). If an
            integer, then the zones will be 10%, 80% and 10% of the width,
            respectively.
        breadth (None or int): Not implemented. Raises an error.
        strat (tuple): Stratigraphy above, in, and below the wedge. This is the
            'reference' stratigraphy. If you give integers, you get 'solid'
            layers containing those numbers. If you give arrays, you will get
            layers of those arrays, expanded or squeezed into the layer
            thicknesses implied in `depth`.
        thickness (tuple): The wedge thickness on the left and on the right.
            Default is (0.0, 1.0) so the wedge will be thickness 0 on the left
            and the wedge thickness specified in the depth argument on the
            right. If the thickness are equal, you'll have a flat, layer-cake
            model.
        mode (str or function): What kind of interpolation to use. Default:
            'linear'. Other options are 'sigmoid', which makes a clinoform-like
            body, 'root', which makes a steep-sided bowl shape like the edge of
            a channel, and 'power', which makes a less steep-sided bowl shape.
            If you pass a function, give it an API like np.linspace: f(start,
            stop, num), where start is the left thickness, stop is the right
            thickness, and num is the width of the middle (wedge) 'zone'.
        conformance (str): 'top', 'bottom', or 'both' (the default). How you
            want the layers inside the wedge to behave. For top and bottom
            conformance, if the layer needs to be thicker than the reference.

    Returns:
        namedtuple[ndarray, ndarray, ndarray, int]: A tuple containing the
            2D wedge model, the top 'horizon', the base 'horizon', and the
            position at which the wedge has thickness 1 (i.e. is the thickness
            specfied by the middle layer depth and/or strat).
    """
    # Decide if binary (2 rocks in the wedge).
    if np.unique(strat[1]).size == 2:
        binary = True
    else:
        binary = False

    if breadth is None or breadth == 0:
        breadth = 1

    # Allow wedge to be thin-thick or thick-thin.
    left, right = thickness
    if left > right:
        left, right = right, left
        flip = True
    else:
        flip = False

    # Get all layers thicknesses.
    if isinstance(depth, int):
        L1, L2, _ = 3 * [depth//3]  # Sizes if depth is just a number.
        L3 = depth - L1 - L2
    else:
        L1, L2, L3 = map(int, depth)
    L3 += int(right * L2)  # Adjust bottom layer.

    # Get all zone widths.
    if isinstance(width, int):
        Z2 = round(0.8 * width)
        Z1 = round(max(1, width/10))
        Z3 = width - Z2 - Z1
        width = int(Z1), int(Z2), int(Z3)
    else:
        Z1, Z2, Z3 = width  # Z1 and Z3 are the bookends.

    # Deal with different interpolation patterns.
    modes = {
        'linear': np.linspace,
        'clinoform': sigmoid,
        'sigmoid': sigmoid,
        'root': root,
        'power': power,
    }
    zooms = modes.get(mode, mode)(left, right, Z2)

    # Get the reference stratigraphy in each layer.
    # The 'well log' case is tricky, because layer1 and layer3
    # need to know about the amount of zoom on the wedge layer.
    # There must be an easier way to do this.
    layer1 = get_strat(strat[0], L1, position=0, wedge=(strat[1], L2))
    layer2_ = get_strat(strat[1], L2, position=1)
    layer3 = get_strat(strat[2], L3, position=-1, wedge=(strat[1], L2))

    # Deal with width. We discard the reference breadth.
    if binary and (breadth >= 2):
        layer2s, _n2g, _ = get_subwedges(layer2_, breadth=breadth)
        layer2s = layer2s.T
    else:
        layer2s, _ = [layer2_], None

    # Make the padding function.
    padder = pad_func(layer1, layer3)

    # For everything in breadth:
    model = []
    for layer2 in layer2s:
        # Collect wedge pieces, then pad top & bottom, then stack,
        # then pad left & right.
        if conformance in ['top', 'bottom', 'base']:
            wedges = [get_conforming(layer2, z*L2, conformance) for z in zooms]
        else:
            wedges = [get_strat(layer2, thickness=z*L2) for z in zooms]
        padded = [np.pad(w, [L1, L3-w.size], mode=padder) for w in wedges]
        wedge = np.pad(np.stack(padded), [[Z1, Z3], [0, 0]], mode='edge')
        model.append(wedge.T)
    model = np.array(model)

    # Make the top and base 'horizons'.
    top = np.repeat((np.ones(np.sum(width)) * L1)[:, None], breadth, axis=-1)
    base_ = np.pad(L1 + zooms * L2, [Z1, Z3], mode='edge')
    base = np.repeat(base_[:, None], breadth, axis=-1)

    # Calculate the reference profile ('well' position).
    if left <= 1 <= right:
        ref = Z1 + np.argmin(np.abs(zooms-1))
    elif left == right == 1:
        ref = Z1 + Z2//2
    else:
        ref = -1

    # Flip if required.
    if flip:
        model = np.flip(model, axis=2)
        base = base[::-1]
        ref = sum(width) - ref

    # Move the 'breadth' dim to last.
    if model.shape[0] > 1:
        model = np.moveaxis(model, 0, -1)

    # Fix any shape discrepancy.
    # if model.shape

    # Build and return output.
    Wedge = namedtuple('Wedge', ['wedge', 'top', 'base', 'reference'])
    return Wedge(np.squeeze(model),
                 np.squeeze(top),
                 np.squeeze(base),
                 ref
                 )
