##Script Name: remove_contaminants.py
##Purpose: Removes seqs annotated as contaminants by job.py from full length genes file. 
##
##Computational Genomics Lab, University of Connecticut, Storrs 
print "Getting started..."
from Bio import SeqIO
import os

#Part 0 - Get file paths
originalntseq = raw_input("Enter the path of the contaminated nucleotide FASTA file (relative is OK!): ")
contaminantsfile = raw_input("Enter the path of the contaminants job.py out file (relative is OK!): ")

#Part 1 - Simplifies FASTA header lines to just ID number. 
outputfile=(open('centroids.fasta.transdecoder.cds.mod','w'))
with open(originalntseq,'r') as pepfile:
		print "Loading original FASTA protein file...."
		for line in pepfile:
			if line.startswith('>'):
				x = line.split('|', 1)[0]
				#print x
				outputfile.write(x)
				outputfile.write('\n')
			else:
				#print line
				outputfile.write(line)
		print "FASTA modification complete!"
outputfile.close()

#Part 2 - Make list of contaminant ID numbers. 
nameslist=open('ntcontaminants.list','w')
with open(contaminantsfile,'r') as contaminants_table:
	print "Loading contaminants table..."
	for line in contaminants_table:
		firstcol = line.split('\t',1)[0]
		idnum = firstcol.split('|', 1)[0]
		nameslist.write(idnum)
		nameslist.write('\n')
print "Contaminant names list generated!"
nameslist.close()

#Import list of header lines.
print "Removing contaminants..."
headerlines=[]
with open('ntcontaminants.list', 'r') as namefile:
	headerlines = namefile.read().splitlines()

#Go through each record in the fasta file and if a record's header line matches, add it to a list.
matches=[]
for record in SeqIO.parse('centroids.fasta.transdecoder.cds.mod', "fasta"):
	if record.id not in headerlines: #Records not in headers of contaminants must be clean hits/
		matches.append(record)
#Print to an output fasta file
print "Writing uncontaminated FASTA file..."
output=open("uncontaminatedseqs.fasta", "w")
SeqIO.write(matches, output, 'fasta')
output.close()
print "Cleaning up..."	
os.remove("centroids.fasta.transdecoder.cds.mod")
os.remove("ntcontaminants.list")
print "Done!"