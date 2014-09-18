## Script Name: headerFormatter.py
## Purpose: Edit headers in multi-fasta file to all follow the same naming convention
## Created by James Pickett
## University of Connecticut
##
## Version: 1.0.0
## Last Edit 9/17/2014
## Usage: Run from terminal, filename can be included in line or a query will appear on run (Wildcard and Relative OK)

import glob
import sys

cmdInput = sys.argv[len(sys.argv) - 1]
target = glob.glob(cmdInput)[0]
if(cmdInput == sys.argv[0]):
	target = glob.glob(raw_input("Enter the name of the file to be parsed (Relative and wildcards OK) \n"))[0]

geneSeq = open(target,'r')  #Retrieve file contents
storageName = raw_input("Enter the name of the file you would like to store output in\n")
storageFile = open(storageName,'w')     ##Commented out to avoid needless runs,  uncomment and change filename and print statements for full functionality

pos = 0

name = ''
sequence = ''

hybridHeaders = []
correctHeaders = []
numberIDs = []


geneSeqString = geneSeq.read() #Store file contents as a single string

run = 1

while True:

	nameStart = geneSeqString.find('>', pos) #Stores the index of the first '>' character in nameStart. '>' signifies the start of a FASTA-format sequence, 'pos' indicates the position to start searching from, and is 0 at first run

	if(nameStart == -1): #Python returns a -1 if the character is not found in the string, indicating there are no more occurences and running again is unnecessary
		break

	nameEnd = geneSeqString.find('\n', nameStart) #Finds the first newline character after the index stored above, which indicates the end of the FASTA sequence description and beginning of sequence data

	name = geneSeqString[nameStart:nameEnd] #Stores the entire sequence desription as name

	if name.find("hybrid") != -1:
		pos = nameStart + 1
		continue

	# correctHeaders.append(name[2:correctHeaders[0].find('_')])
	numberIDs.append(int(name[2:name.find('_')]))

	sequence = geneSeqString[nameEnd + 1:geneSeqString.find('>',nameEnd) - 1] #Stores the sequence

	storageFile.write('%s\n'%name)
	storageFile.write('%s\n'%sequence)

	if(nameStart == -1):
		sequence = geneSeqString[geneSeqString.find('\n',pos) + 1:-1] #If there are no more sequences, the rest of the file is the current sequence

	pos = nameStart + 1
	# print 'C',run
	run += 1

# print "RETEST"
# correctHeaders = sorted(correctHeaders)
# largest = correctHeaders[len(correctHeaders) - 1]
largest = sorted(numberIDs)[len(numberIDs) - 1]
# print "NORMAL: ",numberIDs
# print "SORTED: ",sorted(numberIDs)
# print largest
pos = 0

while True:

	# print "TEST"

	nameStart = geneSeqString.find('>', pos) #Stores the index of the first '>' character in nameStart. '>' signifies the start of a FASTA-format sequence, 'pos' indicates the position to start searching from, and is 0 at first run

	# print "header finding"

	if(nameStart == -1): #Python returns a -1 if the character is not found in the string, indicating there are no more occurences and running again is unnecessary
		break

	nameEnd = geneSeqString.find('\n', nameStart) #Finds the first newline character after the index stored above, which indicates the end of the FASTA sequence description and beginning of sequence data

	name = geneSeqString[nameStart:nameEnd] #Stores the entire sequence desription as name

	# print "typechecking"
	if name.find("hybrid") == -1:
		pos = nameStart + 1
		continue

	largestString = str(largest)
	name = '>c' + largestString + '_g1_i1'
	largest = int(largest) + 1

	sequence = geneSeqString[nameEnd + 1:geneSeqString.find('>',nameEnd) - 1] #Stores the sequence

	storageFile.write('%s\n'%name)
	storageFile.write('%s\n'%sequence)

	if(nameStart == -1):
		sequence = geneSeqString[geneSeqString.find('\n',pos) + 1:-1] #If there are no more sequences, the rest of the file is the current sequence

	pos = nameStart + 1
	# print "H Run"

# for x in range(1,6): #Writes stored sequences to previously specified file       ##Comemented out for same reason as storageFile line
# 	storageFile.write("     ")
# 	storageFile.write("\n")