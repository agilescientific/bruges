#!/usr/bin/env python
# -*- coding: utf 8 -*-
"""
Initialize the library.

:copyright: 2015 Agile Geoscience
:license: Apache 2.0
"""
from . import filters
from . import attribute
from . import rockphysics
from . import reflection
from . import noise
from . import transform
from . import unit
from . import util

__all__ = ["filters",
           "attribute",
	   "rockphysics",
	   "reflection",
	   "noise",
	   "transform",
	   "unit",
	   "util"]

__version__ = "unknown"
try:
    from _version import __version__
except ImportError:
    pass
