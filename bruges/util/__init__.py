#!/usr/bin/env python
# -*- coding: utf 8 -*-
"""
:copyright: 2015 Agile Geoscience
:license: Apache 2.0
"""
from .util import rms
from .util import moving_average
from .util import moving_avg_conv
from .util import moving_avg_fft
from .util import normalize
from .util import next_pow2
from .util import top_and_tail
from .util import extrapolate

__all__ = ["rms",
           "moving_average",
           "moving_avg_conv",
           "moving_avg_fft",
           "normalize",
           "next_pow2",
           "top_and_tail",
           "extrapolate",
           ]
