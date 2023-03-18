class Cell:
    value = None
    back_pointer = None
    operation = None

    def __init__(self, value=None, back_pointer=None, operation=None):
        # Keep track of the cell value
        self.value = value
        # Track the back pointer (another Cell)
        self.back_pointer = back_pointer
        # operation
        self.operation = None
