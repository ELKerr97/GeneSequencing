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
		self.d = 0

# This is the method called by the GUI.  _seq1_ and _seq2_ are two sequences to be aligned, _banded_ is a boolean that tells
# you whether you should compute a banded alignment or full alignment, and _align_length_ tells you
# how many base pairs to use in computing the alignment

	def align( self, seq1, seq2, banded, align_length):
		self.banded = banded
		self.MaxCharactersToAlign = align_length
		self.d = 3
		self.table = None
		k = self.d*2 + 1

		shortest, longest = self.shorter_longer_sequence(seq1, seq2)

		# Make an empty table
		if banded:
			# Table dimensions: k x shortest
			# Time: O(kn)
			self.table = Table(len(shortest), k)
		else:
			# Time: O(nm)
			self.table = Table(len(seq1), len(seq2))

		# Fill base cases
		# Time: O(nm) or O(kn) [banded]
		self.table.fill_base_cases(banded, INDEL, self.d)

		if banded:
			alignment1, alignment2, score = "", "", 0
		else:
			self.compute_cells_non_banded(seq1, seq2)

		# todo: get alignments and score from table
		# given back pointers, values, and operations from table
		return {'align_cost':score, 'seqi_first100':alignment1, 'seqj_first100':alignment2}

	def compute_cells_non_banded(self, seq1, seq2):
		# Compare sequences
		for i in range(1, len(seq1)):
			seq1letter = seq1[i]
			for j in range(1, len(seq2)):
				seq2letter = seq2[j]
				match = seq1letter == seq2letter
				self.compute_cell(i, j, match)

	def compute_cell(self, i, j, match):
		insert = 0
		delete = 1
		match_sub = 2
		# Check if banded
		if self.banded:
			# If cell is past base cases and looking at leftmost column
			# Only delete or match
			if i > self.d + 1 and j == 0:
				delete_cost = self.delete_cost(i, j)
				match_sub_cost = self.match_sub_cost(i, j, match)
				# if delete is smaller, take delete
				if delete_cost < match_sub_cost:
					# Set the value to the new delete cost
					self.table.cells[i][j].value = delete_cost
					# Set the back pointer (up and right)
					self.table.cells[i][j].backpointer = self.table.cells[i-1][j+1]
					self.table.cells[i][j].operation = delete
					return delete
				else:
					# Set the value to the new match/sub cost
					self.table.cells[i][j].value = match_sub_cost
					# Set the back pointer (up)
					self.table.cells[i][j].backpointer = self.table.cells[i-1][j]
					self.table.cells[i][j].operation = match_sub
					return match_sub

			# if at last cell and in base case range, only insert or match
			# last index in cell is d + 1 + row
			elif i <= self.d + 1 and j == self.d + 1 + i:
				match_sub_cost = self.match_sub_cost(i, j, match)
				insert_cost = self.insert_cost(i, j)
				if insert_cost <= match_sub_cost:
					self.table.cells[i][j].value = insert_cost
					self.table.cells[i][j].backpointer = self.table.cells[i-1][j]
					self.table.cells[i][j].operation = insert
					return insert
				else:
					self.table.cells[i][j].value = match_sub_cost
					self.table.cells[i][j].backpointer = self.table.cells[i-1][j-1]
					self.table.cells[i][j].operation = match_sub
					return match_sub

			# if in sweet spot (like normal)
			else:
				insert_cost = self.insert_cost(i, j)
				delete_cost = self.delete_cost(i, j)
				match_sub_cost = self.match_sub_cost(i, j, match)

				if insert_cost <= delete_cost and insert_cost <= match_sub_cost:
					self.table.cells[i][j].value = insert_cost
					self.table.cells[i][j].backpointer = self.table.cells[i][j-1]
					self.table.cells[i][j].operation = insert
					return insert
				elif delete_cost <= insert_cost and delete_cost <= match_sub_cost:
					self.table.cells[i][j].value = delete_cost
					self.table.cells[i][j].backpointer = self.table.cells[i-1][j]
					self.table.cells[i][j].operation = delete
					return delete
				elif match_sub_cost <= insert_cost and match_sub_cost <= delete_cost:
					self.table.cells[i][j].value = match_sub_cost
					self.table.cells[i][j].backpointer = self.table.cells[i-1][j-1]
					self.table.cells[i][j].operation = match_sub
					return match_sub

		# if not banded, do normal computation
		else:
			insert_cost = self.insert_cost(i, j)
			delete_cost = self.delete_cost(i, j)
			match_sub_cost = self.match_sub_cost(i, j, match)

			if insert_cost <= delete_cost and insert_cost <= match_sub_cost:
				self.table.cells[i][j].value = insert_cost
				self.table.cells[i][j].backpointer = self.table.cells[i][j - 1]
				self.table.cells[i][j].operation = insert
				return insert
			elif delete_cost <= insert_cost and delete_cost <= match_sub_cost:
				self.table.cells[i][j].value = delete_cost
				self.table.cells[i][j].backpointer = self.table.cells[i - 1][j]
				self.table.cells[i][j].operation = delete
				return delete
			elif match_sub_cost <= insert_cost and match_sub_cost <= delete_cost:
				self.table.cells[i][j].value = match_sub_cost
				self.table.cells[i][j].backpointer = self.table.cells[i - 1][j - 1]
				self.table.cells[i][j].operation = match_sub
				return match_sub


	def insert_cost(self, i, j):
		# Same row, left one column
		return self.table.cells[i][j-1].value + INDEL

	def delete_cost(self, i, j):
		if self.banded and i > self.d + 1:
			# Up one row, right one column
			return self.table.cells[i-1][j+1].value + INDEL
		else:
			# Up one row, same column
			return self.table.cells[i-1][j].value + INDEL

	def match_sub_cost(self, i, j, match):
		# Check for banded cases
		if self.banded and i > self.d + 1:
			if match:
				# Up one row, same col
				return self.table.cells[i-1][j].value + MATCH
			else:
				return self.table.cells[i-1][j].value + SUB
		else:
			if match:
				# Up one row, back one col
				return self.table.cells[i - 1][j - 1].value + MATCH
			else:
				return self.table.cells[i - 1][j - 1].value + SUB

	def shorter_longer_sequence(self, seq1, seq2):
		if len(seq1) >= len(seq2):
			return seq2, seq1
		else:
			return seq1, seq2
