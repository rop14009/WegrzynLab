#This script is unfinished.
#Start with list of uncontaminated header lines, could pipe through w/ bash file
from Bio import SeqIO

#Import list of header lines.
headerlines=[]
with open('clean_hits.names', 'r') as namefile:
	headerlines = namefile.read().splitlines()
	#Checked to see if above code works, it does. 
	#print headerlines

#Go through each record in the fasta file and if a record's header line matches, add it to a list.
matches=[]
for record in SeqIO.parse("centroids.fasta.transdecoder.cds", "fasta"):
	if record.id in headerlines:
		matches.append(record)
#Print to an output fasta file
output=open("uncontaminatedseqs.fasta", "w")
SeqIO.write(matches, output, 'fasta')
output.close()
		