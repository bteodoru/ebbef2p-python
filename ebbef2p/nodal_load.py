class NodalLoad():
    """Add a :class:`~ebbef2p.nodal_load.NodalLoad` object to the 
    structure model.

    Args:
        value (:obj:`float`):  The value of the nodal load.
        position (:obj:`float`): Location along the beam where load
                is applied.
        type (:obj:`str`): The type of the nodal load. A string 
            identifier that can be either ```my``` to indicate
            a bending moment or ```fz``` to indicate a vertical 
            force. 
            
    Examples:
        .. code-block:: python

            # set a vertical force with 100 magnitude value acting at
            # a distance x = 3 from the left structure end
            NodalLoad(100, 3, 'fz')   

            # set a moment with 50 magnitude value acting at
            # a distance x = 1 from the left structure end
            NodalLoad(50, 1, 'my') 
    """   
    
    def __init__(self, value, position, type):
        self.value = value
        self.position = position
        self.type = type