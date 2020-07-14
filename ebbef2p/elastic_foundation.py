class ElasticFoundation():
    """Elastic foundation class

    Args:
        value (:obj:`tuple` of :obj:`float`): List containing start
            and end values of the elastic foundation support.
        position (:obj:`tuple` of :obj:`float`): List containing 
            start and end coordinate values of the elastic foundation
            support.    
        type (:obj:`str`): The type of the elastic foundation support. 
            A string identifier that can be either ```k``` to indicate
            the first parameter (Winkler type) or ```t``` to indicate
            the second parameter. 
    """

    def __init__(self, value, position, type):
        self.value = value
        self.position = position
        self.type = type
    
    def __str__(self):
        return f"Value: {self.value} \nPosition: {self.position} \nType: {self.type}"