#!/usr/bin/python3

'''
Lutz Bocher, Konrad Diedrich, Steffen Hirte, Claudius Sandmeier
GIK Aufgabe 3.2
Abgabe 18.05.2017
'''

import argparse  	# to read line parameters
import re			# regex

# read the file, parse and check the sequences
def readfile():
	parser = argparse.ArgumentParser(description=
		'Read a list of RNA pairs and convert to Vienna notation')
	parser.add_argument('filename')

	args = parser.parse_args()

	file = open(args.filename, 'r')

	with file as f:
		lines = f.read().splitlines() # read lines without \n

	return lines

# the converter. read out the length and indices for the brackets using regexp,
# replace the dots at the correct positions
def convert_to_vienna(line):

	length = re.search('^(.*)\;', line)

	if length:
		converted = ['.'] * int(length.group(1))

	open_brackets = re.findall('\((\d*),', line)
	closing_brackets = re.findall(',(\d*)\)', line)

	if open_brackets:
		for i in open_brackets:
			converted[int(i)] = '('

	if closing_brackets:
		for i in closing_brackets:
			converted[int(i)] = ')'

	return ''.join(converted)

########### main ##########

lines = readfile()

for line in lines:
	print(convert_to_vienna(line))

