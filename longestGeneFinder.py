## Script Name: longestGeneFinder.py
## Purpose: Find longest sequence for each gene
##
## Created by James Pickett
## University of Connecticut
##
## Version: 1.0.0
## Last Edit 7/16/2014
## Usage: Run from terminal, filename can be included in line or a query will appear on run (Wildcard and Relative OK)

import re
import glob
import sys
import types

cmdInput = sys.argv[len(sys.argv) - 1]
target = glob.glob(cmdInput)[0]
if(cmdInput == sys.argv[0]):
	target = glob.glob(raw_input("Enter the name of the file to be parsed (Relative and wildcards OK) \n"))[0]

geneSeq = open(target,'r')  #Retrieve file contents

lastNameEnd = 0
pos = 0

name = ''
sequence = ''

seqLengthDic = {}
fastaDic = {}

notReachedEnd = True
runOne = True

geneSeqString = geneSeq.read() #Store file contents as a single string



while notReachedEnd:


	nameStart = geneSeqString.find('>', pos) #Stores the index of the first '>' character in nameStart. '>' signifies the start of a FASTA-format sequence, 'pos' indicates the position to start searching from, and is 0 at first run

	if(nameStart == -1): #Python returns a -1 if the character is not found in the string, indicating there are no more occurences and running again is unnecessary
		notReachedEnd = False
		break

	nameEnd = geneSeqString.find('\n', nameStart) #Finds the first newline character after the index stored above, which indicates the end of the FASTA sequence description and beginning of sequence data

	name = geneSeqString[nameStart:nameEnd] #Stores the entire sequence desription as name

	cNum = re.search(r"c[0-9]*", name)

	sequence = geneSeqString[nameEnd + 1:geneSeqString.find('>',nameEnd) - 1] #Stores the sequence

	if(type(cNum) != types.NoneType):
		if(str(cNum.group(0)) in seqLengthDic):
			if(seqLengthDic[str(cNum.group(0))] < len(sequence)):
				seqLengthDic[str(cNum.group(0))] = len(sequence)
				fastaDic[name] = sequence
		else:
			seqLengthDic[str(cNum.group(0))] = len(sequence)
			fastaDic[name] = sequence


	if(nameStart == -1):
		sequence = geneSeqString[geneSeqString.find('\n',pos) + 1:-1] #If there are no more sequences, the rest of the file is the current sequence

	runOne = False

	lastNameEnd = nameEnd
	pos = nameStart + 1

with open('longSequences.fasta', 'w') as storageFile:
	for key in fastaDic:
		storageFile.write(key)
		storageFile.write('\n')
		storageFile.write(fastaDic[key])
		storageFile.write('\n')