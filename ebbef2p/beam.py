class Beam():  
    """A Beam defines a structure member.

    A beam is a slender member that is subjected only to transverse
    loadings.
    It is assumed to have homogeneous properties (e.g. E), with a 
    constant cross-section.

    Args:
        coord (:obj:`list` of :obj:`float`): List containing start
            and end coordinates of structural beam.
        E (:obj:`float`): Young's modulus of the cross-section.
        h (:obj:`float`): Height of the cross-section. 
        w (:obj:`float`): Width of the cross-section. 

    Examples:
        .. code-block:: python

            from ebbef2p import Beam

            # set a beam from 6 to 12 distance from the 
            # left structure end
            # with Young modulus E = 2.70E+07 and 
            # a rectangular cross-section (height h = 1, width w = 0.5)
            b = Beam(coord=[6, 12], E=2.70E+07, h=1, w=0.5)        
    """    

    def __init__(self, coord, E, h, w):
        """Inits the Beam class."""

        self.coord = coord
        self.E = E
        self.h = h
        self.w = w

    @property
    def length(self):
        """The length of the beam.

        Returns:
            :obj:`float`: Difference between beam coordinates.
        """

        return self.coord[1] - self.coord[0]

    @property
    def I(self):
        """Area moment of inertia of the cross-section.

        Returns:
            :obj:`float`: I
        """

        return self.w * self.h**3 / 12

    def __str__(self):
        return f"Coords: {self.coord} \nYoung's modulus: {self.E} \
                \nMoment of inertia: {self.I} \nLength: {self.length}"
