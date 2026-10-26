"""Exceptions raised by EBBEF2p.

Each one also derives from the built-in exception that describes it, so
``except ValueError`` and the like keep working.
"""

import numpy as np


class EBBEF2pError(Exception):
    """Base class for all EBBEF2p errors."""


class ModelError(EBBEF2pError, ValueError):
    """The problem definition is invalid or inconsistent.

    ``field`` names the offending argument when there is a single one
    (``"E"``). ``path`` locates the error in a serialized model
    (``"beams[1].E"``) when it is raised by :mod:`ebbef2p.io`. Both are None
    otherwise.
    """

    def __init__(self, message: str, *, field: str | None = None,
                 path: str | None = None) -> None:
        super().__init__(message)
        self.field = field
        self.path = path


class UnstableModelError(EBBEF2pError, np.linalg.LinAlgError):
    """The beam can move as a rigid body: the system matrix is singular."""


class ConvergenceError(EBBEF2pError, RuntimeError):
    """An iterative procedure (the gamma iteration) did not converge."""
