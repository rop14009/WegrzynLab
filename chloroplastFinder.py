## Script Name: chloroplastFinder.py
## Purpose: Extract the genes from job.py output that are labelled as chloroplast by both databases
##
## Created by James Pickett
## University of Connecticut
##
## Version: 1.0.0
## Last Edit 7/16/2014
## Usage: Run from terminal, filename can be included in line or a query will appear on run (Wildcard and Relative OK)

import glob
import sys

cmdInput = sys.argv[len(sys.argv) - 1]
target = glob.glob(cmdInput)[0]
if(cmdInput == sys.argv[0]):
	target = glob.glob(raw_input("Enter the name of the file to be parsed (Relative and wildcards OK) \n"))[0]

geneSeq = open(target,'r')  #Retrieve file contents
storageFile = open('chloroplasts.txt','w')     ##Commented out to avoid needless runs,  uncomment and change filename and print statements for full functionality

lastNameEnd = 0
pos = 0

name = ''
sequence = ''

storageFile.write(geneSeq.readline())

while True:

	geneSeqString = str.lower(geneSeq.readline())
	if(geneSeqString == ''):
		break

	isChloro = geneSeqString.find('chloroplast') #Stores the index of the first 'chloroplast' substring in isChloro.
	if(isChloro == -1): #Python returns a -1 if the character is not found in the string, indicating there are no more occurences and running again is unnecessary
		continue
	isChloro = geneSeqString.find('chloroplast', isChloro + 1)
	if(isChloro == -1):
		continue

	storageFile.write('%s \n'%geneSeqString)


	if(isChloro == -1):
		sequence = geneSeqString[geneSeqString.find('\n',pos) + 1:-1] #If there are no more sequences, the rest of the file is the current sequence

	pos = isChloro + 1

# for x in range(1,6): #Writes stored sequences to previously specified file       ##Comemented out for same reason as storageFile line
# 	storageFile.write("     ")
# 	storageFile.write("\n")