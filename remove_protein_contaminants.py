##Script Name: remove_protein_contaminants.py
##Purpose: Removes protein seqs annotated as contaminants by job.py from the original FASTA protein file (*.transdecoder.pep)
##
##Author: Ethan Baker
## Computational Genomics Lab, University of Connecticut
##
##Last Update 7/10/14
 
from Bio import SeqIO
print "Getting started..."
originalprotseq = raw_input("Enter the path of the contaminated protein FASTA file (relative is OK!): ")
contaminantsfile = raw_input("Enter the path of the contaminants job.py out file (relative is OK!): ")

#Part 1 - Simplifies FASTA header lines to just ID number. 
outputfile=(open('centroids.fasta.transdecoder.pep.mod','w'))
with open(originalprotseq,'r') as pepfile:
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
		print "FASTA modification complete."
outputfile.close()

#Part 2 - Make list of contaminant ID numbers. 
nameslist=open('proteincontaminants.list','w')
with open(contaminantsfile,'r') as contaminants_table:
	print "Loading contaminants table..."
	for line in contaminants_table:
		firstcol = line.split('\t',1)[0]
		idnum = firstcol.split('|', 1)[0]
		nameslist.write('cds.')
		nameslist.write(idnum)
		nameslist.write('\n')
print "Contaminant names list generated!"
nameslist.close()

#Part 3 - Load contaminant IDs to list.
headerlines=[]
with open('proteincontaminants.list', 'r') as namefile:
	headerlines = namefile.read().splitlines()

#Part 4 -  Remove contaminants from protein seq file
uncontaminated=[]
print 'Removing contaminants...'
for record in SeqIO.parse('centroids.fasta.transdecoder.pep.mod', "fasta"):
	if record.id not in headerlines:
		uncontaminated.append(record)
print "Generating clean sequence file..."
finaloutput=open("uncontaminatedseqs.pep","w")
SeqIO.write(uncontaminated, finaloutput, "fasta")
finaloutput.close()	
print "Done!"		