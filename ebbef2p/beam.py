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
        I (:obj:`float`): Area moment of inertia.
            cross-section. 

    Examples:
        .. code-block:: python

            from ebbef2p import Beam

            # set a beam from 6 to 12 distance from the 
            # left structure end
            # with Young modulus E = 2.70E+07 and 
            # moment of inertia I = 5.21E-03
            b = Beam([6, 12], 2.70E+07, 5.21E-03)        
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
