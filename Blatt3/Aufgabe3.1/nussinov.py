#!/usr/bin/python3

'''
Lutz Bocher, Konrad Diedrich, Steffen Hirte, Claudius Sandmeier
GIK Aufgabe 3.1
Abgabe 18.05.2017
'''

import argparse # to read line parameters

# read the file, parse and check the sequences
def readfile():
	
	parser = argparse.ArgumentParser(description=
	'Calculate the RNA secondary structure using the Nussinov algorithm')
	parser.add_argument('filename')
    
	args = parser.parse_args()
    
	file = open(args.filename,'r')
    
	with file as f:
		l = f.readlines()
		sequences = [] # sequences will be stored in this list
		for line in l:			
			
			# exclude headers and sequence labels
			if (line[0] == '#' or line[0] == '>'):
				continue
			
			# add the line if the bases are correct
			if check_bases(line):
				sequences.append(line)

	return sequences

# exclude faulty sequences
def check_bases(line):
	
	allowed_bases = 'acgu\n'
	correct_bases = True
	
	try:		
		for base in line:
			if base.lower() not in allowed_bases:
				raise ValueError("Wrong character '{}' detected!".format(base))				
				
	except ValueError as err:
		print(err.args)
		correct_bases = False

	return correct_bases

''' 
--------- main ----------
'''

sequences = readfile()
print (sequences)