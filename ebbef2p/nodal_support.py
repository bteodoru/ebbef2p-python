class NodalSupport():
    """Class for a node to be used.

    Args:
        position (:obj:`float`): Location along the beam where 
            support is applied.
        uz (:obj:`float`): Translational restriction in
            vertical direction. Defaults to ``None``.
        ur (:obj:`float`): Rotational restriction. 
            Defaults to ``None``.

    Examples:
        .. code-block:: python

            from ebbef2p import NodalSupport

            # set a roller support
            NodalSupport(position=0, uz=0)

            # set a fixed support
            NodalSupport(position=0, uz=0, ur=0)

            # set a prescribed vertical displacement
            NodalSupport(position=0, uz=1.5) 
    """
    
    def __init__(self, position, uz=None, ur=None):
        self.uz = uz
        self.ur = ur
        self.position = position
    def __str__(self):
        return f"Position: {self.position} \
                \nuz: {self.uz} \nur: {self.ur}"