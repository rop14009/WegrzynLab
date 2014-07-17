## Script Name: find5After.py
## Purpose: Find five largest FASTA-format sequences
##
## Created by James Pickett
## University of Connecticut
##
## Version: 1.1.0
## Last Edit 7/7/2014
##
## Usage: Run from terminal, filename can be included in line or a query will appear on run (Wildcard and Relative OK)
## Note: This script functions exactly as find5.py, but stores the resulting strings with a
## 		different filename, so it can be run on multiple sources in the same directory

import glob
import sys

cmdInput = sys.argv[len(sys.argv) - 1]
target = glob.glob(cmdInput)[0]
if(cmdInput == sys.argv[0]):
	target = glob.glob(raw_input("Enter the name of the file to be parsed (Relative and wildcards OK) \n"))[0]

geneSeq = open(target,'r')  #Retrieve file contents
storageFile = open('FiveLongestSequencesAfter.txt','w') #Change this line and save to make another "version", enabling another source from the same directory

fiveSeqs = ['','','','','']
fiveLengths = [0,0,0,0,0];

currLen = 0
lastNameEnd = 0
pos = 0
addedSeqs = 0

name = ''
sequence = ''

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
	sequence = geneSeqString[nameEnd + 1:geneSeqString.find('>',nameEnd) - 1]


	if(nameStart == -1):
		sequence = geneSeqString[geneSeqString.find('\n',pos) + 1:-1] #If there are no more sequences, the rest of the file is the current sequence

	runOne = False

	currLen = len(sequence)

	if(currLen > sorted(fiveLengths)[0]):
		cullIndex = fiveLengths.index(sorted(fiveLengths)[0])
		fiveLengths.append(currLen)

		fiveLengths.remove(sorted(fiveLengths)[0])
		fiveSeqs.pop(cullIndex)
		fiveSeqs.append(name)
		addedSeqs += 1

	lastNameEnd = nameEnd
	pos = nameStart + 1

	if(nameStart == -1):
		notReachedEnd = False

for x in range(1,6): #Writes stored sequences to previously specified file
	storageFile.write(fiveSeqs[fiveLengths.index(sorted(fiveLengths)[5 - x])])
	storageFile.write("     ")
	storageFile.write(str(sorted(fiveLengths)[5 - x]))
	storageFile.write("\n")