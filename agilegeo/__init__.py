"""
==================
agilegeo
==================
"""

from . import wavelet
from . import attribute
from . import avo
from . import util

__all__ = ["wavelet", "attribute", "avo", "util"]

__version__ = "unknown"
try:
    from _version import __version__
except ImportError:
    pass
