## Script Name: geneFinder.py
## Purpose: Find longest sequence for each gene if none are labelled as complete
##
## Created by James Pickett
## University of Connecticut
##
## Version: 1.0.0
## Last Edit 7/21/2014
## Usage: Run from terminal, filename can be included in line or a query will appear on run (Wildcard and Relative OK)

import re
import glob
import sys
import types

cmdInput = sys.argv[len(sys.argv) - 1]
target = glob.glob(cmdInput)[0]
if(cmdInput == sys.argv[0]):
	target = glob.glob(raw_input("Enter the name of the file to be parsed (Wildcard and Relative OK) \n"))[0]

geneSeq = open(target,'r')  #Retrieve file contents

lastNameEnd = 0
pos = 0

name = ''
sequence = ''

fastaNameDic = {}
cgDic = {}

geneSeqString = geneSeq.read() #Store file contents as a single string



while True:


	nameStart = geneSeqString.find('>', pos) #Stores the index of the first '>' character in nameStart. '>' signifies the start of a FASTA-format sequence, 'pos' indicates the position to start searching from, and is 0 at first run

	if(nameStart == -1): #Python returns a -1 if the character is not found in the string, indicating there are no more occurences and running again is unnecessary
		break

	nameEnd = geneSeqString.find('\n', nameStart) #Finds the first newline character after the index stored above, which indicates the end of the FASTA sequence description and beginning of sequence data

	name = geneSeqString[nameStart:nameEnd] #Stores the entire sequence desription as name

	cNum = str(re.search(r"c[0-9]*", name).group(0))
	gNum = str(re.search(r"g[0-9]*", name).group(0))
	geneTypeID = name.find('complete')
	cgID = cNum + gNum

	sequence = geneSeqString[nameEnd + 1:geneSeqString.find('>',nameEnd) - 1] #Stores the sequence

	if(cgID in cgDic):
		if(geneTypeID != -1):
			cgDic[cgID] = sequence
			fastaNameDic[cgID] = name
		elif(len(sequence) > len(cgDic[cgID])):
			cgDic[cgID] = sequence
			fastaNameDic[cgID] = name
	else:
		cgDic[cgID] = sequence
		fastaNameDic[cgID] = name

	lastNameEnd = nameEnd
	pos = nameStart + 1

with open('Genes.fasta', 'w') as storageFile:
	for key in cgDic:
		storageFile.write(fastaNameDic[key])
		storageFile.write('\n')
		storageFile.write(cgDic[key])
		storageFile.write('\n')
