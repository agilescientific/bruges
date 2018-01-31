# -*- coding: utf-8 -*-
from .wavelets import sinc
from .wavelets import ricker
from .wavelets import ormsby
from .wavelets import sweep
from .wavelets import rotate_phase

from .kernels import gaussian_kernel
from .kernels import gaussian

from .filters import snn
from .filters import kuwahara
from .filters import conservative
from .filters import rms
from .filters import mean
from .filters import median
from .filters import mode

from .anisodiff import anisodiff, anisodiff3
