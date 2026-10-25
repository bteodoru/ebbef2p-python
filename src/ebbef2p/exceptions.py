"""Exceptions raised by EBBEF2p.

Each one also derives from the built-in exception that describes it, so
``except ValueError`` and the like keep working.
"""

import numpy as np


class EBBEF2pError(Exception):
    """Base class for all EBBEF2p errors."""


class ModelError(EBBEF2pError, ValueError):
    """The problem definition is invalid or inconsistent."""


class UnstableModelError(EBBEF2pError, np.linalg.LinAlgError):
    """The beam can move as a rigid body: the system matrix is singular."""


class ConvergenceError(EBBEF2pError, RuntimeError):
    """An iterative procedure (the gamma iteration) did not converge."""
