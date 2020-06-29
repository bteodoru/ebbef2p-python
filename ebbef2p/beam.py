class Beam():  
    """A Beam defines a structure member.

    A beam is a slender member that is subjected only to transverse
    loadings.
    It is assumed to have homogeneous properties (e.g. E), with a 
    constant cross-section.

    Args:
        coord (:obj:`list` of :obj:`float`): List containing start
            and end coordinates of structural beam.
        E (:obj:`float`): Young's modulus of the beam cross-section.
        I (:obj:`float`): Area moment of inertia of the beam 
            cross-section.           
    """    

    def __init__(self, coord, E, I):
        """Inits the Beam class."""

        self.coord = coord
        self.E = E
        self.I = I

    @property
    def length(self):
        """The length of the beam.

        Returns:
            :obj:`float`: Difference between beam coordinates.
        """

        return self.coord[1] - self.coord[0]

    def __str__(self):
        return f"Coords: {self.coord} \nYoung's modulus: {self.E} \nMoment of inertia: {self.I} \nLength: {self.length}"
