#!/usr/bin/env python
# -*- coding: utf 8 -*-
"""
Initialize the library.

:copyright: 2015 Agile Geoscience
:license: Apache 2.0
"""
from . import attribute
from . import filters
from . import petrophysics
from . import rockphysics
from . import reflection
from . import noise
from . import transform
from . import unit
from . import util
from .get_bruges import get_bruges

__all__ = ["attribute",
           "filters",
           "petrophysics",
           "rockphysics",
           "reflection",
           "noise",
           "transform",
           "unit",
           "util",
           "get_bruges"]

__version__ = "unknown"
try:
    from ._version import __version__
except ImportError:
    pass
