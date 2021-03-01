from collections import namedtuple

import numpy as np
import scipy.ndimage as sn

from bruges.util import sigmoid


def pad_func(before, after):
    """
    Padding function. Operates on vector *in place*,
    as per the np.pad documentation.
    """
    def pad_with(x, pad_width, iaxis, kwargs):
        x[:pad_width[0]] = before[-pad_width[0]:]
        x[-pad_width[1]:] = after[:pad_width[1]]
        return
    return pad_with


def get_strat(strat, thickness, mode='stretch', kind='nearest', position=1, wedge=None, zoom_mode='nearest'):
    """
    Take a 'stratigraphy' (either an int, a tuple of ints, or a list-like of
    floats) and expand or compress it to the required thickness.
    
    'mode' can be 'stretch' or 'crop'
    
    `kind` can be 'nearest', 'linear', 'quadratic', or 'cubic'.
    """
    orders = {'nearest': 0, 'linear': 1, 'quadratic': 2, 'cubic': 3}
    order = orders.get(kind, 0)
    
    if isinstance(strat, int) and order==0:
        out = np.repeat([strat], thickness)
    elif isinstance(strat, float) and order==0:
        out = np.repeat([strat], thickness)
    elif isinstance(strat, tuple) and order==0:
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
        return strat[:thickness] if conformance == 'top' else strat[-thickness:]
    else:
        if conformance == 'top':
            return np.pad(strat, [0, thickness-strat.size], mode='wrap')
        else:
            return np.pad(strat, [thickness-strat.size, 0], mode='wrap')
    return


def wedge(depth=(30, 40, 30),
          width=(10, 80, 10),
          breadth=None,  # Not implemented.
          strat=(0, 1, 2),
          thickness=(0.0, 1.0),
          mode='linear',
          conformance='proportional',
         ):
    """
    Generate a wedge model.
    
    Args:
        depth (int or tuple): The vertical size of the model. If a 3-tuple, then
            each element corresponds to a layer. If an integer, then each layer
            of the model will be 1/3 of the thickness. Note that if the 'right'
            wedge thickness is more than 1, then the total thickness will be
            greater than this value.
        width (int or tuple): The width of the model. If a 3-tuple, then each
            element corresponds to a 'zone' (left, middle, right). If an integer,
            then the zones will be 10%, 80% and 10% of the width, respectively.
            Note: for small models, you'll get better results with an odd number
            of pixels in the middle (wedge) zone.
        breadth (int or tuple): Not implemented. Raises an error.
        strat (tuple): Stratigraphy above, in, and below the wedge. This is the
            'reference' stratigraphy. If you give integers, you get 'solid' layers
            containing those numbers. If you give arrays, you will get layers of
            those arrays, expanded or squeezed into the layer thicknesses implied
            in `depth`. 
        thickness (tuple): The wedge thickness on the left and on the right.
            Default is (0.0, 1.0) so the wedge will be thickness 0 on the left
            and the wedge thickness specified in the depth argument on the right.
            If the thickness are equal, you'll have a flat, layer-cake model.
        mode (str): What kind of interpolation to use. Default: 'linear'. Other
            option is 'clinoform' (or synonymously 'sigmoid'), which makes a
            clinoform-like body.
        conformance (str): 'top', 'bottom', or 'proportional' (the default,
            although you can use any word other than 'top', 'bottom' or 'base').
            How you want the layers inside the wedge to behave. For top and
            bottom conformance, if the layer needs to be thicker than the
            reference position it will be wrapped (i.e. repeated) using the
            `'wrap'` mode of `np.pad()`.

    Returns:
        namedtuple[ndarray, ndarray, ndarray, int]: A tuple containing the
            2D wedge model, the top 'horizon', the base 'horizon', and the
            position at which the wedge has thickness 1 (i.e. is the thickness
            specfied by the middle layer depth and/or strat).

    TODO
    - Nearest interp ints, but linear interp floats (e.g. rock properties).
    - Breadth argument implements the third dimension.
    - If the wedge layer is a tuple of two ints, e.g. (1, 2, 1, 2, 1), then
        you are making a 'binary wedge', which has special features.
    """
    if breadth is not None:
        raise NotImplementedError("The breadth argument is not implemented yet.")

    # Allow wedge to be thin-thick or thick-thin.
    left, right = thickness
    if left > right:
        left, right = right, left
        flip = True
    else:
        flip = False

    if isinstance(depth, int):
        L1, L2, L3 = 3 * [depth//3]  # Sizes if depth is just a number.
        L3 += 1
    else:
        L1, L2, L3 = map(int, depth)
    L3 += int(right * L2)  # Adjust bottom layer.

    if isinstance(width, int):
        Z1, Z2, Z3 = width // 10, int(0.8 * zones), zones // 10
    else:
        Z1, Z2, Z3 = width  # Z1 and Z3 are the bookends.

    if mode == 'linear':
        zooms = np.linspace(left, right, Z2)
    elif mode in ['clinoform', 'sigmoid']:
        zooms = sigmoid(left, right, 1+Z2)
    else:
        raise TypeError("Mode not recognized.")

    # Get the reference stratigraphy in each layer.
    # The 'well log' case is tricky, because layer1 and layer3
    # need to know about the amount of zoom on the wedge layer.
    # There must be an easier way to do this.
    layer1 = get_strat(strat[0], L1, position=0, wedge=(strat[1], L2))
    layer2 = get_strat(strat[1], L2, position=1)
    layer3 = get_strat(strat[2], L3, position=-1, wedge=(strat[1], L2))

    padder = pad_func(layer1, layer3)

    # Collect wedge pieces, then pad top & bottom, then stack, then pad left & right.
    if conformance in ['top', 'bottom', 'base']:
        wedges = [get_conforming(layer2, z*L2, conformance) for z in zooms]
    else:
        wedges = [get_strat(layer2, thickness=z*L2) for z in zooms]
    padded = [np.pad(w, [L1, L3-w.size], mode=padder) for w in wedges] 
    wedge = np.pad(np.stack(padded), [[Z1, Z3], [0, 0]], mode='edge')

    # Make the top and base 'horizons'.
    top = np.ones(np.sum(width)) * L1
    base = np.pad(L1 + zooms * L2, [Z1, Z3], mode='edge')

    # Calculate the reference profile ('well' position).
    if left <= 1 <= right:
        ref = Z1 + np.argmin(np.abs(zooms-1))
    elif left == right == 1:
        ref = Z1 + Z2//2
    else:
        ref = -1

    if flip:
        wedge = np.flipud(wedge)
        base = base[::-1]
        ref = sum(width) - ref

    Wedge = namedtuple('Wedge', ['wedge', 'top', 'base', 'reference'])
    return Wedge(wedge.T, top, base, ref)
