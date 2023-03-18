from Cell import Cell


# A Table to hold the cells
class Table:

    def __init__(self, num_rows, num_cols):
        self.num_rows = num_rows
        self.num_cols = num_cols
        # Populate the table with Cell objects
        self.cells = [[Cell() for j in range(num_cols)] for i in range(num_rows)]

    def fill_base_cases(self, indel_cost):
        insert = 0
        delete = 1
        match_sub = 2
        # Set the start cell to have value 0
        self.cells[0][0].value = 0
        # Set the start cell to have no back pointer
        self.cells[0][0].back_pointer = None
        # Set the rightmost column
        for i in range(1, self.num_rows):
            # Set value
            self.cells[i][0].value = self.cells[i-1][0].value + indel_cost
            # Set back pointer
            self.cells[i][0].back_pointer = self.cells[i-1][0]
            self.cells[i][0].operation = delete

        # Set the top row
        for i in range(1, self.num_cols):
            # Set value
            self.cells[0][i].value = self.cells[0][i-1].value + indel_cost
            # Set back pointer
            self.cells[0][i].back_pointer = self.cells[0][i-1]
            self.cells[0][i].operation = insert
