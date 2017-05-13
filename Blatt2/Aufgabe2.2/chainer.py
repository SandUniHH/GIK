#!/usr/bin/python3
# The code was written testing with python version 3.5,
# comment out the magic string if necessary.
#
# Lutz Bocher, Konrad Diedrich, Steffen Hirte, Claudius Sandmeier
# GIK Aufgabe 2.2
# Abgabe 04.05.2017
#
# Note: there will be a diff in line 59 or line 60 depending
# on whether the lines are reversed (search for lines.reverse() in this file)
# before sorting the matches. By default, the matches.txt lines are sorted by weight.
# The diff is due to same values for the first sequence (values 0 and 1 of the line),
# but different ones for the second one, because the former sequence seems to have matched
# with two sequences in the latter.

import argparse # to read line parameters

# Read length of both genomes and matches.txt
# each line has the format start1-end1-start2-end2-weight
def readfile():

    parser = argparse.ArgumentParser(description=
         'Find the optimal sequence chain using the simple chaining algorithm')
    parser.add_argument('seqlen1', type=int)
    parser.add_argument('seqlen2', type=int)
    parser.add_argument('filename')

    args = parser.parse_args()

    length_1 = args.seqlen1
    length_2 = args.seqlen2
    file = open(args.filename, 'r')

    # read the file, split the values each line and convert them to integers
    with file as f:
        l = f.readlines()
        lines = []  # the list we will be working with
        for line in l:
            lines.append([int(n) for n in line.split()])

        # lines.reverse() # uncomment to inspect the diff in line 59/60 of matches.txt and chain.txt

        # sort by end value of first sequence of the match
        lines = sorted(lines, key = lambda x: int(x[1]))

    return length_1, length_2, lines

# the main optimal score calculating function
def getOptimalScore(matches):

    overallmaxscore = 0
    bestmatch = None

    score = {} #[0] * len(matches) # initialize scores for each line with value 0
    prec = {} # create dictionary

    for j, line in enumerate(matches):

        weight = line[4]
        maxscore = weight
        prec[j] = None # default value for each dictionary entry, beginning of sequence chain

        for i in range(j):

            # fi << fj
            if (matches[i][1] < line[0] and matches[i][3] < line[2] and
                maxscore < score[i] + weight):

                maxscore = score[i] + weight
                prec[j] = i # set index of predecessor

        score[j] = maxscore

        if overallmaxscore < maxscore:
            overallmaxscore = maxscore
            bestmatch = j

    return overallmaxscore, bestmatch, prec

# calculate the values for output and construct list of matches of the chain
def calculateChain(prec, matches, bestmatch):
    chain = []
    seq1_max = 0
    seq2_max = 0
    chain_length = 0
    i = bestmatch

    while i is not None:
        seq1_max += matches[i][1] - matches[i][0] + 1
        seq2_max += matches[i][3] - matches[i][2] + 1

        # convert the integers back to a string line
        chain.append(" ".join(str(x) for x in matches[i]))

        chain_length += 1
        i = prec[i] # get the match before the current one

    return seq1_max, seq2_max, chain_length, chain

########### main ##########

# read the file and parameters
seqlen1, seqlen2, lines = readfile()

# calculate the optimal score
overallmaxscore, bestmatch, prec = getOptimalScore(lines)

# prepare the output
seq1_max, seq2_max, chain_length, chain = calculateChain(prec, lines, bestmatch)

# and finally, the output
for i in reversed(chain):
    print(i)
print('# optimal chain of length {} with score {}'.format(chain_length, overallmaxscore))
print('# {} bp covered on sequence 1 (coverage {:.2f}%)'.format(seq1_max, seq1_max / seqlen1 * 100))
print('# {} bp covered on sequence 2 (coverage {:.2f}%)'.format(seq2_max, seq2_max / seqlen2 * 100))