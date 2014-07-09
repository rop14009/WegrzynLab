## Script Name: completeFinder.py
## Purpose: Find complete FASTA-format sequences
##
## Created by James Pickett
## University of Connecticut
##
## Version: 1.1.0
## Last Edit 7/9/2014
## Usage: Run from command prompt, and when prompted enter the name of the FASTA file exactly

import math

target = raw_input("Enter the name of the file to be parsed ") #Asks user to input a filename, and stores the name input
geneSeq = open(target,'r')  #Retrieve file contents *Only works with exact filenames* retrieval is read-only
storageFile = open('CompleteSequences.txt','w') #Declares a file to store the completed sequences in, file is write-only
statsFile = open('CompleteSequencesStats.txt', 'w') #Declares a file to record statistics on the completed sequences in, file is write-only

completeSeqs = []
sequences = []
completeLengths = []

currLen = 0
lastNameEnd = 0
numberComplete = 0
pos = 0
lengthSum = 0
mean = 0
lengthMax = 0

name = ''
sequence = ''

notReachedEnd = True

geneSeqString = geneSeq.read() #Store file contents as a single string



while notReachedEnd:

	nameStart = geneSeqString.find('>', pos) #Stores the index of the first '>' character in nameStart. '>' signifies the start of a FASTA-format sequence, 'pos' indicates the position to start searching from, and is 0 at first run

	if(nameStart == -1): #Python returns a -1 if the character is not found in the string, indicating there are no more occurences and running again is unnecessary
		notReachedEnd = False
		break

	nameEnd = geneSeqString.find('\n', nameStart) #Finds the first newline character after the index stored above, which indicates the end of the FASTA sequence description and beginning of sequence data

	name = geneSeqString[nameStart:nameEnd] #Stores the entire sequence desription as name
	currLen = name[name.find('len:') + 4:]

	if(name.lower().find("complete") != -1): #Checks if the description contains 'complete' (Not case sensitive), indicating the sequence is known to be complete
		completeSeqs.append(name) #Adds the complete sequence to a list
		numberComplete += 1
	else:
		pos = nameStart + 1
		continue

	sequence = geneSeqString[nameEnd + 1:geneSeqString.find('>', nameEnd) - 1] #Stores the sequence between the indices already declared


	if(geneSeqString.find('>', nameEnd) == -1): #Checks if there is another sequence following the current one
		sequence = geneSeqString[geneSeqString.find('\n',nameEnd) + 1:] #If there are no more sequences, the rest of the file is the current sequence

	sequences.append(sequence)

	completeLengths.append(len(sequence)) #Stores the length of the parsed sequence
	if(len(sequence) > lengthMax): #Checks if the current length is longer than the previous record, and if so, assigns the current as the new high water mark
		lengthMax = len(sequence)

	lastNameEnd = nameEnd
	pos = nameStart + 1

for x in range (0, len(completeLengths)):
	lengthSum += completeLengths[x]

mean = lengthSum / len(completeLengths)
if(len(completeLengths) % 2 != 0): #If statement to find median sequence length
	median = sorted(completeLengths)[len(completeLengths) / 2] #If the number of sequences is odd, the sequence in the middle is the median
else:
	median = (sorted(completeLengths)[len(completeLengths) / 2 - 1] + sorted(completeLengths)[len(completeLengths) / 2]) / 2 #If the number is even, the median is the average of the middle two
with open("completeSequences.fasta", "w") as fastaFile:
	for x in range (0, len(completeSeqs)): #For loop to write 
		storageFile.write(completeSeqs[x])
		storageFile.write("     ")
		storageFile.write(str(completeLengths[x]))
		storageFile.write("\n")
		storageFile.write(sequences[x])
		storageFile.write("\n")
		fastaFile.write(completeSeqs[x])
		fastaFile.write("\n")
		fastaFile.write(sequences[x])
		fastaFile.write("\n")

statsFile.write("Max: ")				#Writes to statistics file in order of: max sequence length, min length, mean and median lengths, and count of complete sequences, each on a seperate line
statsFile.write(str(lengthMax))
statsFile.write("\n")
statsFile.write("Min: ")
statsFile.write(str(sorted(completeLengths)[0]))
statsFile.write("\n")
statsFile.write("Mean: ")
statsFile.write(str(mean))
statsFile.write("\n")
statsFile.write("Median: ")
statsFile.write(str(median))
statsFile.write("\n")
statsFile.write("Count: ")
statsFile.write(str(len(completeSeqs)))