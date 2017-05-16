#!/usr/bin/python3

"""
Lutz Bocher, Konrad Diedrich, Steffen Hirte, Claudius Sandmeier
GIK Aufgabe 3.1
Abgabe 18.05.2017
"""

import argparse  # to read line parameters
import collections  # to use a matrix with uninitialised elements


# read the file, parse and check the sequences
def readfile():
	parser = argparse.ArgumentParser(description=
									 'Calculate the RNA secondary structure using the Nussinov algorithm')
	parser.add_argument('filename')

	args = parser.parse_args()

	file = open(args.filename, 'r')

	with file as f:
		l = f.readlines()
		seqs = []  # sequences will be stored in this list
		for line in l:

			# exclude headers and sequence labels
			if line[0] == '#' or line[0] == '>':
				continue

			# add the line/sequence if the bases are correct
			if line[-1] == '\n':
				line = line[:-1]
			line = line.lower()
			if check_bases(line):
				seqs.append(line)

	return seqs


# exclude faulty sequences
def check_bases(line):
	allowed_bases = 'acgu'
	correct_bases = True

	try:
		for base in line:
			if base not in allowed_bases:
				raise ValueError("Wrong character '%s' detected!" % base)

	except ValueError as err:
		print(err.args)
		correct_bases = False

	return correct_bases


# energy function alpha
def alpha_function(base1, base2):
	if (base1 == 'g' and base2 == 'c') or \
			(base1 == 'c' and base2 == 'g'):
		return -3
	elif (base1 == 'a' and base2 == 'u') or \
			(base1 == 'u' and base2 == 'a'):
		return -2
	elif (base1 == 'g' and base2 == 'u') or \
			(base1 == 'u' and base2 == 'g'):
		return -1
	else:
		# return probably only needs to be larger than -1 in this case,
		# but this makes the meaning clearer
		return float("inf")

	# obtain the minimum value of cases 1, 2a, 2b and 2c


def e_min_value(i, j, e, sequence):
	equation1 = e[i + 1, j - 1] + \
				alpha_function(sequence[i - 1], sequence[j - 1])
	equation2 = e[i + 1, j]
	equation3 = e[i, j - 1]

	equation4 = float("inf")
	for k in range(i + 2, j):
		current = e[i, k - 1] + e[k, j]

		if current < equation4:
			equation4 = current

	equation_values = [equation1, equation2, equation3, equation4]

	return min(equation_values)


# evaluate all secondary structures for a sequence and
# return the minimum free energy
def calculate_mfe(sequence):
	n = len(sequence)
	l_min = 3

	# The E matrix we will use for the dynamic programming approach.
	# It allows undefined values and easy addition of additional values.
	# This is reasonably efficient since we are working at > O(n²) complexity.
	e = collections.defaultdict(float)

	for i in range(1, n + 1):
		e[i, i] = 0

	for l in range(2, n + 1):
		for i in range(1, n + 1 - l + 1):

			j = i + l - 1  # j is a substring of length l

			if (j - i - 1) < l_min:
				e[i, j] = 0
			else:
				e[i, j] = e_min_value(i, j, e, sequence)

	return e


# Trace back through the sequence using the E matrix. From this, the secondary
# structure can be derived, i. e. the loops formed.
def traceback(sequence, e, i, j):
	if i < j:
		if e[i, j] == e[i + 1, j - 1] + \
				alpha_function(sequence[i - 1], sequence[j - 1]):
			print('(%i,%i)' % (i, j))
			traceback(sequence, e, i + 1, j - 1)

		elif e[i, j] == e[i + 1, j]:
			traceback(sequence, e, i + 1, j)

		elif e[i, j] == e[i, j - 1]:
			traceback(sequence, e, i, j - 1)

		else:
			for k in range(i + 2, j):
				if e[i, j] == e[i, k - 1] + e[k, j]:
					traceback(sequence, e, i, k - 1)
					traceback(sequence, e, k, j)
					break


''' 
--------- main ----------
'''

# read the sequences, check if they are valid and preformat
sequences = readfile()

print('<fold2Dmulti>')

for sequence in sequences:
	# construct the matrix with the mfe values
	matrix_e = calculate_mfe(sequence)

	print('<fold2D>')
	print('<seq>%s</seq>' % sequence.lower())
	print('<pairs>')
	traceback(sequence, matrix_e, 1, len(sequence))  # lots of printing
	print('</pairs>')
	print('<mfe>%i</mfe>' % matrix_e[1, len(sequence)])
	print('</fold2D>')

print('</fold2Dmulti>')