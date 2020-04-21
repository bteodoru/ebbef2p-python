class DistributedLoad():
    """
    *value* : list
        [x_start, x_end]

    *position* : list
        [x_start, x_end]

    """

    def __init__(self, value, position):
        self.value = value
        self.position = position

    def __str__(self):
        return f"value: {self.value} \nposition: {self.position}"