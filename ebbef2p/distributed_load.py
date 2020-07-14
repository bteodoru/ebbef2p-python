class DistributedLoad():
    """Add a :class:`~ebbef2p.distributed_load.DistributedLoad` 
    object to the structure model.

    Args:
        value (:obj:`tuple` of :obj:`float`): List containing start
            and end values of the distributed load.
        position (:obj:`tuple` of :obj:`float`): List containing 
            start and end coordinate values of the distributed load.
    
    Examples:
        .. code-block:: python

            form ebbef2p import DistributedLoad

            # set an uniformly distributed load with 10 magnitude value
            # acting from 2.5 to 3.75 distance from the 
            # left structure end
            dl = DistributedLoad((10, 10), (2.5, 3.75))
            
            # set an trianglar distributed load with 0 magnitude value 
            # at the left end and 10 magnitude value at the right end
            # acting from 2.5 to 3.75 distance from the 
            # left structure end
            dl = DistributedLoad((0, 10), (2.5, 3.75))
            
            # set an trapezoidal distributed load with 10 magnitude value 
            # at the left end and 25 magnitude value at the right end
            # acting from 2.5 to 3.75 distance from the 
            # left structure end
            dl = DistributedLoad((10, 25), (2.5, 3.75))
    """

    def __init__(self, value, position):
        self.value = value
        self.position = position

    def __str__(self):
        return f"value: {self.value} \nposition: {self.position}"