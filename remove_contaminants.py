#Script relies on a list of uncontaminated names. Must be 
#Start with list of uncontaminated header lines, could pipe through w/ bash file
from Bio import SeqIO

#Import list of header lines.
headerlines=[]
with open('contaminants.names', 'r') as namefile:
	headerlines = namefile.read().splitlines()

#Go through each record in the fasta file and if a record's header line matches, add it to a list.
matches=[]
for record in SeqIO.parse("centroids.fasta.transdecoder.cds", "fasta"):
	if record.id not in headerlines: #Records not in headers of contaminants must be clean hits/
		matches.append(record)
#Print to an output fasta file
output=open("uncontaminatedseqs.fasta", "w")
SeqIO.write(matches, output, 'fasta')
output.close()
		