#!/usr/bin/python3

from which_pyqt import PYQT_VER
if PYQT_VER == 'PYQT5':
	from PyQt5.QtCore import QLineF, QPointF
# elif PYQT_VER == 'PYQT4':
# 	from PyQt4.QtCore import QLineF, QPointF
# elif PYQT_VER == 'PYQT6':
# 	from PyQt6.QtCore import QLineF, QPointF
else:
	raise Exception('Unsupported Version of PyQt: {}'.format(PYQT_VER))

import random

# Used to compute the bandwidth for banded version
MAXINDELS = 3

# Used to implement Needleman-Wunsch scoring
MATCH = -3
INDEL = 5
SUB = 1

# Used for checking if a square's previous pointers
TOP_INDEX = 0
LEFT_INDEX = 1
DIAGONAL_INDEX = 2


class GeneSequencing:

	def __init__( self ):
		pass

# This is the method called by the GUI.  _seq1_ and _seq2_ are two sequences to be aligned, _banded_ is a boolean
# that tells you whether you should compute a banded alignment or full alignment, and _align_length_ tells you
# how many base pairs to use in computing the alignment

	def align(self, seq1, seq2, banded, align_length):
		self.banded = banded
		self.MaxCharactersToAlign = align_length

###################################################################################################
# your code should replace these three statements and populate the three variables: score, alignment1 and alignment2
		# TODO Call my functions
		score = random.random()*100
		alignment1 = 'abc-easy  DEBUG:({} chars,align_len={}{})'.format(
			len(seq1), align_length, ',BANDED' if banded else '')
		alignment2 = 'as-123--  DEBUG:({} chars,align_len={}{})'.format(
			len(seq2), align_length, ',BANDED' if banded else '')
###################################################################################################

		return {'align_cost':score, 'seqi_first100':alignment1, 'seqj_first100':alignment2}

#   | -   A   G   T   C   G   A
# --+---------------------------
# - | 0   5   10  15  20  25  30
# A | 5  -3   2   7   12  17  22
# T | 10  2  -2  -1   4   9   14
# C | 15  7   3  -1  -4   1   6
# G | 20  12  4   4   0  -7  -2
# T | 25  17  9   1   5  -2  -6  <- Edit Distance

	def align_unrestricted(self, seq1: str, seq2: str):
		"""
		Aligns seq1 and seq2 using the edit distance algorithm with Needleman-Wunsch cost values:
		- Match: -3	(diagonal)
		- Sub:    1	(diagonal; letters don't match)
		- Indel:  5	(left or above)

		:param seq1: First sequence, as a string
		:param seq2: Second sequence, as a string
		:return: All dynamic sub-problems (edit distance comparisons), as a dictionary
			[(row, col)] : (edit-distance, [top-prev, left-prev, diagonal-prev{list of 3 booleans}])
		"""
		cache = {}

		pass

	def align_banded(self, seq1: str, seq2: str):
		"""
		Same as align_unrestricted. However, this only considers sub-alignments where corresponding
		letters are withing a bandwidth of each other.
		Bandwidth: 2d + 1 and d = 3: = 7

		Aligns seq1 and seq2 using the edit distance algorithm
		:param seq1: First sequence, as a string
		:param seq2: Second sequence, as a string
		:return:  All dynamic sub-problems (edit distance comparisons), as a dictionary
			[(row, col)] : (edit-distance, [top-prev, left-prev, diagonal-prev{list of 3 booleans}])
		"""
		cache = {}

		pass

	def build_alignments_strings(self, cache: dict, seq1_len: int, seq2_len: int):
		"""
		Computes the two alignment strings from a pre-computed cache
		Breaks ties with ordering: top, left, diagonal

		:param cache: All dynamic sub-problems (edit distance comparisons), as a dictionary
			[(row, col)] : (edit-distance, [top-prev, left-prev, diagonal-prev{list of 3 booleans}])
		:param seq1_len: Length of first sequence, as an int
		:param seq2_len: Length of second sequence, as an int
		:return: The two aligned sequences, as a tuple of two strings, or None if no alignment was possible
		"""
		pass
		# todo first check if there is a valid solution (for banded algorithm)

