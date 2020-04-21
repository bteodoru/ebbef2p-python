class NodalSupport():
    def __init__(self, constraints, position):
        self.constraints = constraints
        self.position = position

    def __str__(self):
        return f"Constraints: {self.constraints} \nPosition: {self.position}"