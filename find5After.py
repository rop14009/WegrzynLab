##Find five largest FASTA-format sequences

import numpy as np

target = raw_input("Enter the name of the file to be parsed ")
geneSeq = open(target,'r')  #Retrieve file contents
storageFile = open('FiveLongestSequencesAfter.txt','w')

fiveSeqs = ['','','','','']
fiveLengths = [0,0,0,0,0];

currMax = 0
currLen = 0
lastNameEnd = 0
temp2 = 0
pos = 0
addedSeqs = 0

name = ''
sequence = ''

notReachedEnd = True
runOne = True

geneSeqString = geneSeq.read() #Store file contents as a single string



while notReachedEnd:


	nameStart = geneSeqString.find('>', pos)
	if(nameStart == -1):
		notReachedEnd = False
		break

	nameEnd = str.find(geneSeqString,'\n', nameStart)

	name = geneSeqString[nameStart:nameEnd]
	sequence = geneSeqString[nameEnd + 1:geneSeqString.find('>',nameEnd) - 1]


	if(nameStart == -1):
		sequence = geneSeqString[geneSeqString.find('\n',pos) + 1:-1]

	runOne = False

	currLen = len(sequence)

	if(currLen > sorted(fiveLengths)[0]):
		cullIndex = fiveLengths.index(sorted(fiveLengths)[0])
		fiveLengths.append(currLen)

		fiveLengths.remove(sorted(fiveLengths)[0])
		fiveSeqs.pop(cullIndex)
		fiveSeqs.append(name)
		addedSeqs += 1

	currMax += 1
	lastNameEnd = nameEnd
	pos = nameStart + 1

	if(nameStart == -1):
		notReachedEnd = False

for x in range(1,6):
	storageFile.write(fiveSeqs[fiveLengths.index(sorted(fiveLengths)[5 - x])])
	storageFile.write("     ")
	storageFile.write(str(sorted(fiveLengths)[5 - x]))
	storageFile.write("\n")
print fiveSeqs
print fiveLengths