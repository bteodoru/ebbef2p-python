class ElasticFoundationSupport():

    """ Sc class

    Attributes:
        value (list): Float values 
        position (list):
        type (str): 'k' and 't' 

    """

    def __init__(self, value, position, type):
        """Inits the ElasticFoundationSupport class.

            value (list):
            position (list):
            type (str):
            
        """
        self.value = value
        self.position = position
        self.type = type
    
    def __str__(self):
        return f"Value: {self.value} \nPosition: {self.position} \nType: {self.type}"


    # def __init__(self, k, t):
    #     self.k = k
    #     self.t = t

    # def __str__(self):
    #     return f"k: {self.k} \nt: {self.t}"
