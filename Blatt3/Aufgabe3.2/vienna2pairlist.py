#!/usr/bin/python3

'''
Lutz Bocher, Konrad Diedrich, Steffen Hirte, Claudius Sandmeier
GIK Aufgabe 3.2
Abgabe 18.05.2017
'''

import argparse  	# to read line parameters

# read the file, parse and check the sequences
def readfile():
	parser = argparse.ArgumentParser(description=
		'Read a list of Vienna notated RNA pairs')
	parser.add_argument('filename')

	args = parser.parse_args()

	file = open(args.filename, 'r')

	with file as f:
		viennas = f.read().splitlines() # read lines without \n

	return viennas

# the converter. create lists of the opening and closing brackets,
# insert the correct positions, pair them together by list index.
def convert_to_pairlist(vienna):

	printline = '{'

	open = []
	close = [0] * int(len(vienna) / 2)
	open_no = 0

	for i, char in enumerate(vienna):

		if char is '(':
			open.append(i)
			open_no += 1
		elif char is ')':
			close[len(open) - open_no] = i
			open_no -= 1

	for i in range(len(open)):
		printline += '({},{}),'.format(open[i], close[i])

	printline = printline[::-1].replace(',', '}', 1)[::-1]

	return printline

########### main ##########

viennas = readfile()

for line in viennas:
	print ('{}; {}'.format(len(line), convert_to_pairlist(line)))