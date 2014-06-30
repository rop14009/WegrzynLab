##Find five largest FASTA-format sequences

import numpy as np
import math

target = raw_input("Enter the name of the file to be parsed ")
geneSeq = open(target,'r')  #Retrieve file contents
storageFile = open('CompleteSequences.txt','w')
statsFile = open('CompleteSequencesStats.txt', 'w')

completeSeqs = []
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

	nameStart = geneSeqString.find('>', pos)

	if(nameStart == -1):
		notReachedEnd = False
		break

	nameEnd = geneSeqString.find('\n', nameStart)

	name = geneSeqString[nameStart:nameEnd]
	currLen = name[name.find('len:') + 4:]

	if(name.lower().find("complete") != -1):
		completeSeqs.append(name)
		numberComplete += 1
	else:
		pos = nameStart + 1
		continue

	name = geneSeqString[nameStart:nameEnd]

	sequence = geneSeqString[nameEnd + 1:geneSeqString.find('>', nameEnd) - 1]


	if(geneSeqString.find('>', nameEnd) == -1):
		sequence = geneSeqString[geneSeqString.find('\n',nameEnd) + 1:]

	completeLengths.append(len(sequence))
	if(len(sequence) > lengthMax):
		lengthMax = len(sequence)

	lastNameEnd = nameEnd
	pos = nameStart + 1

	if(nameStart == -1):
		notReachedEnd = False
		break

for x in range (0, len(completeLengths)):
	lengthSum += completeLengths[x]

mean = lengthSum / len(completeLengths)
if(len(completeLengths) % 2 != 0):
	median = sorted(completeLengths)[len(completeLengths) / 2]
else:
	median = (sorted(completeLengths)[len(completeLengths) / 2 - 1] + sorted(completeLengths)[len(completeLengths) / 2]) / 2

for x in range (0, len(completeSeqs)):
	storageFile.write(completeSeqs[x])
	storageFile.write("     ")
	storageFile.write(str(completeLengths[x]))
	storageFile.write("\n")

statsFile.write("Max: ")
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