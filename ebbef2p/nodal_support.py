class NodalSupport():
    """Class for a node to be used.

    Args:
        constraints (:obj:`dict`): Maps node constraints. Keys
            are string identifiers to indicate the restriction 
            type and values are floats representing prescribed
            translational and rotational displacements.
        position (:obj:`float`): Location along the beam where 
            support is applied.

    Examples:
        .. code-block:: python

            NodalSupport({'uz': 0, 'ur': "NaN"}, 0) # set a roller support
            NodalSupport({'uz': 0, 'ur': 0}, 0) # set a fixed support
            NodalSupport({'uz': 1.5, 'ur': 0}, 0) # set a prescribed vertical displacement

    Todo:
        Set a constraints class
    """

    def __init__(self, constraints, position):
        self.constraints = constraints
        self.position = position

    def __str__(self):
        return f"Constraints: {self.constraints} \nPosition: {self.position}"