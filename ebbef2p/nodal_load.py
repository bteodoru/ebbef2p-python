class NodalLoad():
    """Add a :class:`~ebbef2p.nodal_load.NodalLoad` object to the 
    structure model.

    Args:
        value (:obj:`float`):  The value of the nodal load.
        position (:obj:`float`): Location along the beam where load
                is applied.
        type (:obj:`str`): The type of the nodal load. A string 
            identifier that can either be ```my``` to indicate
            a bending moment or ```fz``` to indicate a vertical 
            force. 
    """   
    
    def __init__(self, value, position, type):
        self.value = value
        self.position = position
        self.type = type