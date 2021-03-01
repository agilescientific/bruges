"""
Initialize the library.

:copyright: 2015 Agile Geoscience
:license: Apache 2.0
"""
from . import attribute
from . import filters
from . import models
from . import petrophysics
from . import rockphysics
from . import reflection
from . import noise
from . import transform
from . import unit
from . import util
from .get_bruges import get_bruges
from .bruges import BrugesError

__version__ = "unknown"
try:
    from ._version import __version__
except ImportError:
    pass
