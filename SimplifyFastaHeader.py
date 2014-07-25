from Bio import SeqIO
originalprotseq = raw_input("Enter the path of the original protein FASTA file (relative is OK!): ")
outputfile=(open('Protein_Genes_simple.fasta','w'))
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