## Script Name: fastaSorter.py
## Purpose: Seperate FASTA input sequences into outputs above and below user-supplied threshold
##
## Created by James Pickett
## University of Connecticut
##
## Version: 1.0.0
## Last Edit 7/16/2014
## Usage: Run from terminal, filename can be included in line or a query will appear on run (Wildcard and Relative OK)

import glob
import sys
try:
	cmdInput = sys.argv[1:]
	target = glob.glob(cmdInput[0])[0]
	threshold = cmdInput[1]
except IndexError:
	target = glob.glob(raw_input("Enter the name of the file to be parsed (Relative and wildcards OK) \n"))[0]
	threshold = raw_input("Enter the length you would like to sort sequences by \n")
threshold = int(threshold)
target = str(target)
print 'Using FASTA file: ', target
print 'Sorting threshold: ', threshold


geneSeq = open(target,'r')  #Retrieve file contents

aboveList = []
belowList = []

pos = 0

name = ''
sequence = ''


geneSeqString = geneSeq.read() #Store file contents as a single string



while True:


	nameStart = geneSeqString.find('>', pos) #Stores the index of the first '>' character in nameStart. '>' signifies the start of a FASTA-format sequence, 'pos' indicates the position to start searching from, and is 0 at first run

	if(nameStart == -1): #Python returns a -1 if the character is not found in the string, indicating there are no more occurences and running again is unnecessary
		break

	nameEnd = geneSeqString.find('\n', nameStart) #Finds the first newline character after the index stored above, which indicates the end of the FASTA sequence description and beginning of sequence data

	name = geneSeqString[nameStart:nameEnd] #Stores the entire sequence desription as name
	sequence = geneSeqString[nameEnd + 1:geneSeqString.find('>',nameEnd) - 1] #Stores the sequence
	if(len(sequence) >= threshold):
		aboveList.append(name)
		aboveList.append(sequence)
	else:
		belowList.append(name)
		belowList.append(sequence)


	if(nameStart == -1):
		sequence = geneSeqString[geneSeqString.find('\n',pos) + 1:-1] #If there are no more sequences, the rest of the file is the current sequence

	pos = nameStart + 1

print belowList[:10]

with open('largeSequences.fasta','w') as storage:
	for index in range(len(aboveList)):
		storage.write("%s\n"%aboveList[index])
with open('smallSequences.fasta','w') as storage:
	for index in range(len(belowList)):
		storage.write("%s\n"%belowList[index])