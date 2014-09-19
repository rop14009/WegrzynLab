## Script Name: interproAppender.py
## Purpose: Append info from interpro output into annotation file
## Created by James Pickett
## University of Connecticut
##
## Version: 1.2.0
## Last Edit 9/19/2014
## Usage: Run from terminal, filename can be included in line or a query will appear on run (Wildcard and Relative OK)

import glob

target = glob.glob(raw_input('Enter the name of the annotation file (Relative and wildcards OK) \n'))[0]
jobInput = glob.glob(raw_input('Enter the name of the interpro output file (Relative and wildcards OK) \n'))[0]

with open(target,'r') as inputFile:
	annoArray = inputFile.read().split('\n')
with open(jobInput,'r') as inputFile:
	interArray = inputFile.read().split('\n')


modAnnoArray = []
modInterArray = []
sequencePFams = []
targets = []
misses = []
lastContents = []
addContents = []

for fields in annoArray:
	fields = fields.split('\t')
	modAnnoArray.append(fields)
	targets.append(fields[0][:fields[0].find('|')])

for fields in interArray:
	modInterArray.append(fields.split('\t'))

for fields in modInterArray:
	fields[0] = fields[0][fields[0].find('.') + 1:]
	if len(fields) < 13:
		for i in range(13 - len(fields)):
			fields.append('')

for i in range(len(modInterArray) - 1):
	try:
		pos = targets.index(modInterArray[i][0])
	except ValueError:
		misses.append(modInterArray[i][0])
		for j in range(len(modAnnoArray[0]) - 1):
			misses[len(misses) - 1] += '\t'
		misses[len(misses) - 1] += '\t' + (modInterArray[i][len(modInterArray[i]) - 3])
		misses[len(misses) - 1] += '\t' + (modInterArray[i][len(modInterArray[i]) - 9])
		misses[len(misses) - 1] += '\t' + (modInterArray[i][len(modInterArray[i]) - 2])
		misses[len(misses) - 1] += '\t' + (modInterArray[i][len(modInterArray[i]) - 1])
		continue

	addContents.append(modInterArray[i][len(modInterArray[i]) - 3])
	addContents.append(modInterArray[i][len(modInterArray[i]) - 9])
	addContents.append(modInterArray[i][len(modInterArray[i]) - 2])
	addContents.append(modInterArray[i][len(modInterArray[i]) - 1])

	for fields in addContents:
		if fields in lastContents:
			pass
		else:
			modAnnoArray[pos].append(fields.strip('\n'))


	lastContents = addContents
	addContents = []

output = raw_input('Enter filename you would like output stored in (Relative and wildcards OK, include filetype) \n')

with open(output,'w') as outputFile:
	for fields in modAnnoArray:
	 	outputFile.write('\t'.join(fields)+'\n')

	outputFile.write('END OF ANNOTATION FILE. INTERPRO MISMATCHES LISTED BELOW.\n')

	for fields in misses:
		outputFile.write('%s\n'%fields)