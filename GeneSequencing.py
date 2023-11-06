#!/usr/bin/python3

from which_pyqt import PYQT_VER
# if PYQT_VER == 'PYQT5':
# 	from PyQt5.QtCore import QLineF, QPointF
# elif PYQT_VER == 'PYQT4':
# 	from PyQt4.QtCore import QLineF, QPointF
# elif PYQT_VER == 'PYQT6':
# 	from PyQt6.QtCore import QLineF, QPointF
# else:
# 	raise Exception('Unsupported Version of PyQt: {}'.format(PYQT_VER))

import random

# Used to compute the bandwidth for banded version
MAXINDELS = 3

# Used to implement Needleman-Wunsch scoring
MATCH = -3
INDEL = 5
SUB = 1

# Used to access edit values in a tuple in the cache
# [(row, col)] = (edit_distance, [prev_ptrs_list])
EDIT_DISTANCE_INDEX = 0
PREV_PTRS_INDEX = 1

# Used for checking if a square's previous pointers
LEFT_INDEX = 0
TOP_INDEX = 1
DIAGONAL_INDEX = 2

# Represents an indel
GAP = '-'

# Used to access row and col of a (row, col) tuple
ROW = 0
COL = 1

# Only for banded algorithm
D = 3  # num items to each side of the middle diagonal (where i == j)
BANDWIDTH = 7

NO_ALIGNMENT_POSSIBLE_MSG = "No Alignment Possible"

#   | -   A   G   T   C   G   A
# --+---------------------------
# - | 0   5   10  15  20  25  30
# A | 5  -3   2   7   12  17  22
# T | 10  2  -2  -1   4   9   14
# C | 15  7   3  -1  -4   1   6
# G | 20  12  4   4   0  -7  -2
# T | 25  17  9   1   5  -2  -6  <- Edit Distance


class GeneSequencing:

	def __init__(self):
		pass

	def align(self, seq1: str, seq2: str, banded: bool, align_length: int):
		"""
		Wrapper method to call either the unrestricted or banded alignment algorithm.

		:param seq1: The first sequence, as a string
		:param seq2: The second sequence, as a string
		:param banded: Whether the banded algorithm should be used, as a boolean
		:param align_length: The number base pairs to use in alignment computation, as an int
		:return: The alignment_cost, alignment of seq1 and seq2 with number of base pairs specified by align_length,
		stored in a dictionary
		"""
		self.banded = banded
		self.MaxCharactersToAlign = align_length

		# trim the sequences to the specified length
		seq1_trimmed = seq1[:self.MaxCharactersToAlign]
		seq2_trimmed = seq2[:self.MaxCharactersToAlign]

		# call the specified function
		if self.banded:
			cache = align_banded(seq1_trimmed, seq2_trimmed)
		else:
			cache = align_unrestricted(seq1_trimmed, seq2_trimmed)

		# get the score
		if (len(seq1_trimmed), len(seq2_trimmed)) in cache:
			score = cache[(len(seq1_trimmed), len(seq2_trimmed))][EDIT_DISTANCE_INDEX]
		else:
			score = float('inf')

		# backtrack through the result to create the alignment strings
		alignment_diff1, alignment_diff2 = build_alignments_strings(cache, seq1_trimmed, seq2_trimmed)

		# format alignment strings for easy debugging
		alignment1 = '{}  DEBUG:({} chars,align_len={}{})'.format(
			alignment_diff1, len(seq1_trimmed), align_length, ',BANDED' if banded else '')
		alignment2 = '{}  DEBUG:({} chars,align_len={}{})'.format(
			alignment_diff2, len(seq2_trimmed), align_length, ',BANDED' if banded else '')

		print("4")

		# # FIXME DELETEME
		# score = random.random()*100
		# alignment1 = 'abc-easy  DEBUG:({} chars,align_len={}{})'.format(
		# 	len(seq1), align_length, ',BANDED' if banded else '')
		# alignment2 = 'as-123--  DEBUG:({} chars,align_len={}{})'.format(
		# 	len(seq2), align_length, ',BANDED' if banded else '')

		return_val = {'align_cost': score, 'seqi_first100': alignment1, 'seqj_first100': alignment2}
		print(return_val)

		return return_val


def align_unrestricted(seq1: str, seq2: str):
	"""
	Aligns seq1 and seq2 using the edit distance algorithm with Needleman-Wunsch cost values:
	- Match: -3	(diagonal)
	- Sub:    1	(diagonal; letters don't match)
	- Indel:  5	(left or above)

	:param seq1: First sequence, as a string
	:param seq2: Second sequence, as a string
	:return: All dynamic sub-problems (edit distance comparisons), as a dictionary
		[(row, col)] : (edit-distance, [left-prev, top-prev, diagonal-prev]{list of 3 booleans})
	"""
	# base case top left
	cache = {(0, 0): (0, [False] * 3)}

	# base case rows
	for row in range(1, len(seq1) + 1):
		# rows: 0, 5, 10, 15, ...
		cache[(row, 0)] = (row * INDEL, [False, True, False])

	# base case cols
	for col in range(1, len(seq2) + 1):
		# cols: 0, 5, 10, 15, ...
		cache[(0, col)] = (col * INDEL, [True, False, False])

	# compute remaining sub-problems
	for row in range(1, len(seq1) + 1):
		for col in range(1, len(seq2) + 1):
			# find possible values for the sub-problem
			top = INDEL + cache[(row - 1, col)][EDIT_DISTANCE_INDEX]
			left = INDEL + cache[(row , col - 1)][EDIT_DISTANCE_INDEX]
			diagonal = match(row, col, seq1, seq2) + cache[(row - 1, col - 1)][EDIT_DISTANCE_INDEX]

			# result is the least cost distance
			edit_distance = min(top, left, diagonal)

			# set prev pointers to those that tied for lowest result
			prev_ptrs = [left == edit_distance, top == edit_distance, diagonal == edit_distance]

			# cache the result of the sub-problem
			cache[(row, col)] = (edit_distance, prev_ptrs)

	return cache

#   - f i r s t w o r d
# - * * * *
# s * o - - >
# d * < o - - >
# w * < - o - - >
# o   < - - o - - >
# r     < - - o - - >
# d       < - - o - - >


def align_banded(seq1: str, seq2: str):
	"""
	Same as align_unrestricted. However, this only considers sub-alignments where corresponding
	letters are withing a bandwidth of each other.
	Bandwidth: 2d + 1 and d = 3: = 7

	Aligns seq1 and seq2 using the edit distance algorithm
	:param seq1: First sequence, as a string
	:param seq2: Second sequence, as a string
	:return:  All dynamic sub-problems (edit distance comparisons), as a dictionary
		[(row, col)] : (edit-distance, [left-prev, top-prev, diagonal-prev]{list of 3 booleans})
	"""

	# base case top left
	cache = {(0, 0): (0, [False] * 3)}

	# base case rows
	for row in range(1, D + 1):
		cache[(row, 0)] = (row * INDEL, [False, True, False])

	# base case cols
	for col in range(1, D + 1):
		cache[(0, col)] = (col * INDEL, [True, False, False])

	# size of the diagonal (not including base case of null)
	diagonal_len = min(len(seq1), len(seq2))

	# compute remaining sub-problems
	for curr_middle in range(1, diagonal_len + 1):
		for col in range(curr_middle - D, curr_middle + 1 + D):
			# if within left/right bounds
			if 0 < col < len(seq2) + 1:
				# find possible values for the sub-problem
				top = float('inf') if (curr_middle - 1, col) not in cache else INDEL + cache[(curr_middle - 1, col)][EDIT_DISTANCE_INDEX]
				left = float('inf') if (curr_middle, col - 1) not in cache else INDEL + cache[(curr_middle, col - 1)][EDIT_DISTANCE_INDEX]
				diagonal = match(curr_middle, col, seq1, seq2) + cache[(curr_middle - 1, col - 1)][EDIT_DISTANCE_INDEX]

				# result is the least cost distance
				edit_distance = min(top, left, diagonal)

				# set prev pointers to those that tied for lowest result
				prev_ptrs = [left == edit_distance, top == edit_distance, diagonal == edit_distance]

				# cache the result of the sub-problem
				cache[(curr_middle, col)] = (edit_distance, prev_ptrs)

	return cache


def build_alignments_strings(cache: dict, seq1: str, seq2: str):
	"""
	Computes the two alignment strings from a pre-computed cache
	Breaks ties with ordering: left, top, diagonal

	:param cache: All dynamic sub-problems (edit distance comparisons), as a dictionary
		[(row, col)] : (edit-distance, [left-prev, top-prev, diagonal-prev]{list of 3 booleans})
	:param seq1: The first sequence, as a string
	:param seq2: The second sequence, as a string
	:return: The two aligned sequences, as a tuple of two strings, or (None, None) if no alignment was possible
	"""
	if (len(seq1), len(seq2)) not in cache:
		# no valid solution (only possible in banded)
		return NO_ALIGNMENT_POSSIBLE_MSG, NO_ALIGNMENT_POSSIBLE_MSG

	aligned1_arr = []
	aligned2_arr = []

	# iterate through the cache to build the alignment strings
	curr_row = len(seq1)
	curr_col = len(seq2)
	while curr_row != 0 or curr_col != 0:
		prev_ptrs = cache[(curr_row, curr_col)][PREV_PTRS_INDEX]

		# get letters corresponding to current row and col
		# NOTE - remember that this is offset by row 0 and col 0 being base cases
		seq1_char = seq1[curr_row - 1]
		seq2_char = seq2[curr_col - 1]

		assert prev_ptrs[LEFT_INDEX] or prev_ptrs[TOP_INDEX] or prev_ptrs[DIAGONAL_INDEX]  # FIXME DELETEME

		# add to alignment strings depending on prev square
		# NOTE - order here maintains precedence: left, dop, diagonal
		if prev_ptrs[LEFT_INDEX]:
			# came from left -> word on left (seq1) gets the gap
			aligned1_arr.append(GAP)
			aligned2_arr.append(seq2_char)

			# move left for next loop
			curr_col -= 1
		elif prev_ptrs[TOP_INDEX]:
			# came from above -> word on top (seq2) gets the gap
			aligned1_arr.append(seq1_char)
			aligned2_arr.append(GAP)

			# move up for next loop
			curr_row -= 1
		elif prev_ptrs[DIAGONAL_INDEX]:
			# came from diagonal -> both words get their letter
			aligned1_arr.append(seq1_char)
			aligned2_arr.append(seq2_char)

			# move up and left for next loop
			curr_row -= 1
			curr_col -= 1

	# reverse arrays and convert to strings
	aligned1 = "".join(aligned1_arr[::-1])
	aligned2 = "".join(aligned2_arr[::-1])

	return aligned1, aligned2


def match(row: int, col: int, seq1: str, seq2: str):
	"""
	Checks if the letters in seq1 and seq2 corresponding to the row and col in the table match.
	NOTE - The table is offset by row 0 and col 0 being occupied by null characters to support base cases

	:param row: The row with the char in seq1
	:param col: The col with the char in seq2
	:param seq1: The first sequence, as a string
	:param seq2: The second sequence, as a string
	:return: Val of MATCH if a match, else val of SUB for a substitution
	"""
	return MATCH if seq1[row - 1] == seq2[col - 1] else SUB


def tests():
	covid_a = "ATCGT"
	covid_b = "AGTCGA"
	# covid_a = "GATTACAATAGCTGCATGCA"
	# covid_b = "ATGCATGCGATCGACT"

	# # unrestricted
	# solution_cache = align_unrestricted(covid_a, covid_b)
	# print(solution_cache)
	# result1, result2 = build_alignments_strings(solution_cache, covid_a, covid_b)
	# print(f"Alignment1: {result1}")
	# print(f"Alignment2: {result2}")

	# banded
	solution_cache = align_banded(covid_a, covid_b)
	print(solution_cache)
	result1, result2 = build_alignments_strings(solution_cache, covid_a, covid_b)
	print(f"Alignment1: {result1}")
	print(f"Alignment2: {result2}")


if __name__ == "__main__":
	tests()
