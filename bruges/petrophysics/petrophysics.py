# -*- coding: utf-8 -*-
'''
===================
petrophysics.py
===================

:copyright: 2016 Agile Geoscience
:license: Apache 2.0
'''
import numpy as np


def gardner(vp, alpha=0.31, beta=0.25):
    '''
    Computes Gardner's density prediction.

    SI units only.

    Args:
        vp

    Returns:
        rhob estimate.
    '''
    return alpha * vp ** beta
