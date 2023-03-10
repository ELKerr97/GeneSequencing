from Cell import Cell


# A Table to hold the cells
class Table:
    def __init__(self, num_rows, num_cols):
        self.rows = num_rows
        self.cols = num_cols
        # Populate the table with Cell objects
        self.cells = [[Cell() for j in range(num_cols)] for i in range(num_rows)]
