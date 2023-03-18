#!/usr/bin/python3
from Table import Table
from Cell import Cell
from which_pyqt import PYQT_VER
if PYQT_VER == 'PYQT5':
	from PyQt5.QtCore import QLineF, QPointF
elif PYQT_VER == 'PYQT4':
	from PyQt4.QtCore import QLineF, QPointF
elif PYQT_VER == 'PYQT6':
	from PyQt6.QtCore import QLineF, QPointF
else:
	raise Exception('Unsupported Version of PyQt: {}'.format(PYQT_VER))

import random

# Used to compute the bandwidth for banded version (d)
MAXINDELS = 3

# Used to implement Needleman-Wunsch scoring
MATCH = -3
INDEL = 5
SUB = 1

class GeneSequencing:

	def __init__( self ):
		self.table = None

# This is the method called by the GUI.  _seq1_ and _seq2_ are two sequences to be aligned, _banded_ is a boolean that tells
# you whether you should compute a banded alignment or full alignment, and _align_length_ tells you
# how many base pairs to use in computing the alignment

	def align( self, seq1, seq2, banded, align_length):
		self.banded = banded
		self.MaxCharactersToAlign = align_length
		self.d = 3
		self.table = None
		self.k = self.d*2 + 1

		# make the shorter one sequence 1 and the longer one sequence 2
		# this will make space O(nk) where n is the shorter of the sequences
		seq1, seq2 = self.shorter_longer_sequence(seq1, seq2)
		if banded:
			alignment1, alignment2, score = self.align_sequences_banded(seq1, seq2)
		else:
			alignment1, alignment2, score = self.align_sequences(seq1, seq2)

		# given back pointers, values, and operations from table
		return {'align_cost':score, 'seqi_first100':alignment1, 'seqj_first100':alignment2}

	def align_sequences_banded(self, seq1, seq2):
		num_rows = min(len(seq1), self.MaxCharactersToAlign)
		num_cols = min(len(seq2), self.MaxCharactersToAlign)

		# get the max and min indices as a list
		max_indices = [min(num_cols, i + self.d + 1) for i in range(num_rows)]
		min_indices = [max(0 + 1, i - self.d + 1) for i in range(num_rows)]

		# Initialize a 2D array to store the table between each pair of characters
		table = [[Cell() for j in range(num_cols + 1)] for i in range(num_rows + 1)]

		# operations to determine back pointers
		insert = "insert"
		delete = "delete"
		sub = "sub"
		match = "match"

		# Initialize the first row and column of the table array (banded)
		table[0][0].value = 0
		for i in range(1, self.d + 1):
			table[i][0].value = i*INDEL
		for j in range(1, self.d + 1):
			table[0][j].value = j*INDEL

		for i in range(1, num_rows + 1):
			for j in range(min_indices[i-1], max_indices[i-1] + 1):
				cell = table[i][j]
				# check if match
				if seq1[i - 1] == seq2[j - 1]:
					# set value
					cell.value = table[i - 1][j - 1].value + MATCH
					# set back pointer
					cell.operation = match
				else:
					top_edge_determiner = num_cols - num_rows

					# check if j at start. If so, only do delete or sub
					if i > self.d and j == min_indices[i-1]:
						cell.value = min(
							table[i - 1][j].value + INDEL,  # delete
							table[i - 1][j - 1].value + SUB  # sub
						)
						if cell.value == table[i - 1][j].value + INDEL:
							cell.operation = delete
						else:
							cell.operation = sub

					# check if j at end
					elif j == max_indices[i - 1]:

						# if j on top edge
						if top_edge_determiner <= self.d and i > num_rows - self.d + top_edge_determiner:
							# set value
							cell.value = min(
								table[i - 1][j].value + INDEL,  # delete
								table[i][j - 1].value + INDEL,  # insert
								table[i - 1][j - 1].value + SUB  # sub
							)

							# set back pointer
							if cell.value == table[i - 1][j].value + INDEL:
								cell.operation = delete
							elif cell.value == table[i][j - 1].value + INDEL:
								cell.operation = insert
							else:
								cell.operation = sub

						# if j not on top edge
						else:
							cell.value = min(
								table[i][j - 1].value + INDEL,  # insert
								table[i - 1][j - 1].value + SUB  # sub
							)
							if cell.value == table[i][j - 1].value + INDEL:
								cell.operation = insert
							else:
								cell.operation = sub

					# if j not on bottom or top edge, treat it like a normal cell
					else:
						# set value
						cell.value = min(
							table[i - 1][j].value + INDEL,  # delete
							table[i][j - 1].value + INDEL,  # insert
							table[i - 1][j - 1].value + SUB  # sub
						)

						# set back pointer
						if cell.value == table[i - 1][j].value + INDEL:
							cell.operation = delete
						elif cell.value == table[i][j - 1].value + INDEL:
							cell.operation = insert
						else:
							cell.operation = sub

		alignment1 = ""
		alignment2 = ""

		# Trace back through the table array to construct the aligned strings
		i = num_rows
		j = num_cols
		if j > i + self.d:
			score = table[i][i].value
		else:
			score = table[i][j].value
		while i > 0 or j > 0:
			cell = table[i][j]
			if cell.operation == match or cell.operation == sub:
				alignment1 = seq1[i - 1] + alignment1
				alignment2 = seq2[j - 1] + alignment2
				i -= 1
				j -= 1
			elif cell.operation == delete:
				alignment1 = seq1[i - 1] + alignment1
				alignment2 = "-" + alignment2
				i -= 1
			elif cell.operation == insert:
				alignment1 = "-" + alignment1
				alignment2 = seq2[j - 1] + alignment2
				j -= 1
			elif cell.operation == None:
				alignment1 = "-" + alignment1
				alignment2 = seq2[j - 1] + alignment2
				score += INDEL
				j -= 1

		return alignment1[:100], alignment2[:100], score

	def align_sequences(self, seq1, seq2):
		num_rows = min(len(seq1), self.MaxCharactersToAlign)
		num_cols = min(len(seq2), self.MaxCharactersToAlign)

		# Initialize a 2D array to store the table between each pair of characters
		table = [[Cell() for j in range(num_cols + 1)] for i in range(num_rows + 1)]

		insert = "insert"
		delete = "delete"
		sub = "sub"
		match = "match"

		# Initialize the first row and column of the table array
		table[0][0].value = 0
		for i in range(1, num_rows + 1):
			table[i][0].value = i*INDEL
		for j in range(1, num_cols + 1):
			table[0][j].value = j*INDEL

		# Compute the table between each pair of characters
		for j in range(1, num_cols + 1):
			for i in range(1, num_rows + 1):
				cell = table[i][j]
				# check if match
				if seq1[i - 1] == seq2[j - 1]:
					# set value
					cell.value = table[i - 1][j - 1].value + MATCH
					# set back pointer
					cell.operation = match
				else:
					# set value
					cell.value = min(table[i - 1][j].value + INDEL, table[i][j - 1].value + INDEL, table[i - 1][j - 1].value + SUB)
					# set back pointer
					if cell.value == table[i - 1][j].value + INDEL:
						cell.operation = delete
					elif cell.value == table[i][j - 1].value + INDEL:
						cell.operation = insert
					else:
						cell.operation = sub

		# Initialize variables to store the aligned strings
		alignment1 = ""
		alignment2 = ""

		# Trace back through the table array to construct the aligned strings
		i = num_rows
		j = num_cols
		score = table[i][j].value
		while i > 0 or j > 0:
			cell = table[i][j]
			if cell.operation == match or cell.operation == sub:
				alignment1 = seq1[i - 1] + alignment1
				alignment2 = seq2[j - 1] + alignment2
				i -= 1
				j -= 1
			elif cell.operation == delete:
				alignment1 = seq1[i - 1] + alignment1
				alignment2 = "-" + alignment2
				i -= 1
			else:
				alignment1 = "-" + alignment1
				alignment2 = seq2[j - 1] + alignment2
				j -= 1

		return alignment1[:100], alignment2[:100], score

	def del_cost(self, i, j, distances):
		return distances[i-1][j]

	def insert_cost(self, i, j, distances):
		return distances[i][j-1]

	def match_sub_cost(self, i, j, distances):
		return distances[i-1][j-1]

	# def get_alignments_and_score(self, seq1, seq2):
	# 	final_cell = self.table.cells[-1][-1]
	# 	alignment1 = ""
	# 	alignment2 = ""
	# 	score = final_cell.value
	# 	done = False

	# def compute_cells(self, seq1, seq2):
	# 	# Compare sequences
	# 	for i in range(1, len(seq1) + 1):
	# 		seq1letter = seq1[i - 1]
	# 		for j in range(1, len(seq2) + 1):
	# 			seq2letter = seq2[j - 1]
	# 			match = seq1letter == seq2letter
	# 			self.compute_cell(i, j, match)

	# def compute_cells_banded(self, seq1, seq2):
	# 	max_indices = [min(len(seq2) - 1, i + self.d) for i in range((len(seq1)))]
	# 	min_indices = [max(0, i - self.d) for i in range(len(seq1))]
	# 	for i in range(1, len(seq1) + 1):
	# 		seq1letter = seq1[i - 1]
	# 		for j in range(min_indices[i] + 1, max_indices[i-1] + 2):
	# 			seq2letter = seq2[j - 1]
	# 			match = seq1letter == seq2letter
	# 			lower_edge = False
	# 			upper_edge = False
	# 			# if j is at min index and i > d, only delete and match
	# 			if i > self.d and j == min_indices[0] + 1:
	# 				lower_edge = True
	# 			# if j is at max index and i <= len(seq2) - d, only insert and match
	# 			if i <= len(seq2) - self.d and j == max_indices[i-1] + 2:
	# 				upper_edge = True
	#
	# 			self.compute_cell(i, j, match, lower_edge, upper_edge)

	def compute_cell(self, i, j, match):
		insert = 0
		delete = 1
		match_sub = 2

		insert_cost = self.insert_cost(i, j)
		delete_cost = self.delete_cost(i, j)
		match_sub_cost = self.match_sub_cost(i, j, match)

		if insert_cost <= delete_cost and insert_cost <= match_sub_cost:
			self.table.cells[i][j].value = insert_cost
			self.table.cells[i][j].back_pointer = self.table.cells[i][j - 1]
			self.table.cells[i][j].operation = insert
			return insert
		elif delete_cost <= insert_cost and delete_cost <= match_sub_cost:
			self.table.cells[i][j].value = delete_cost
			self.table.cells[i][j].back_pointer = self.table.cells[i - 1][j]
			self.table.cells[i][j].operation = delete
			return delete
		elif match_sub_cost <= insert_cost and match_sub_cost <= delete_cost:
			self.table.cells[i][j].value = match_sub_cost
			self.table.cells[i][j].back_pointer = self.table.cells[i - 1][j - 1]
			self.table.cells[i][j].operation = match_sub
			return match_sub

	# def insert_cost(self, i, j):
	# 	# Same row, left one column
	# 	return self.table.cells[i][j-1].value + INDEL
	#
	# def delete_cost(self, i, j):
	# 	if self.banded and i > self.d:
	# 		# Up one row, right one column
	# 		return self.table.cells[i-1][j+1].value + INDEL
	# 	else:
	# 		# Up one row, same column
	# 		return self.table.cells[i-1][j].value + INDEL

	# def match_sub_cost(self, i, j, match):
	# 	# Check for banded cases
	# 	if self.banded and i > self.d:
	# 		if match:
	# 			# Up one row, same col
	# 			return self.table.cells[i-1][j].value + MATCH
	# 		else:
	# 			return self.table.cells[i-1][j].value + SUB
	# 	else:
	# 		if match:
	# 			# Up one row, back one col
	# 			return self.table.cells[i - 1][j - 1].value + MATCH
	# 		else:
	# 			return self.table.cells[i - 1][j - 1].value + SUB

	def shorter_longer_sequence(self, seq1, seq2):
		if len(seq1) >= len(seq2):
			return seq2, seq1
		else:
			return seq1, seq2
