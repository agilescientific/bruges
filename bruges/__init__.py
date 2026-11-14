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

from importlib.metadata import PackageNotFoundError, version

try:
    VERSION = version(__name__)
except PackageNotFoundError:
    try:
        from ._version import version as VERSION
    except ImportError:
        # Source checkout without an install: setuptools_scm has not
        # generated _version.py yet. Do not fail the import — the
        # version is only informational here. A real install (pip
        # install .) always provides one of the two above.
        VERSION = "0.0.0+unknown"
__version__ = VERSION
