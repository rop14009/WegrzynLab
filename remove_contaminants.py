#Purpose: Removes seqs annotated as contaminants from full length genes file. 
#Script relies on a list of uncontaminated names. See remove_contaminants_pipe.sh
#Computational Genomics Lab, University of Connecticut, Storrs 
from Bio import SeqIO

#Import list of header lines.
headerlines=[]
with open('contaminants.names', 'r') as namefile:
	headerlines = namefile.read().splitlines()

#Go through each record in the fasta file and if a record's header line matches, add it to a list.
matches=[]
for record in SeqIO.parse(raw_input("Enter FASTA file path (relative OK): "), "fasta"):
	if record.id not in headerlines: #Records not in headers of contaminants must be clean hits/
		matches.append(record)
#Print to an output fasta file
output=open("uncontaminatedseqs.fasta", "w")
SeqIO.write(matches, output, 'fasta')
output.close()
		