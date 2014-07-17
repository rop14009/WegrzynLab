#!/bin/bash -l
#SBATCH -J trinotate
#SBATCH -o trinotate-%j.output
#SBATCH -e trinotate-%j.err
#SBATCH -n 28
#SBATCH -N 1
#SBATCH -p bigmem
#SBATCH --mail-type=ALL
#SBATCH --mail-user=eagbaker@ucdavis.edu

#PURPOSE: Runs entire trinotate pipeline. Requires Trinity.fasta and transdecoder.pep (with removed contaminants file) to be in same dir from where script will be run. 
#Based on the Trinotate documentation at http://http://trinotate.sourceforge.net/ 
#Ethan Baker, Computational Genomics Lab, University of Connecticut

##################
###Load modules###
##################
module load gmap
module load blast
module load hmmer/2.3.2 
module load trinity
module load trinotate
module load ncbi-blast+
module load express
module load signalp
module load rnammer
module load tmhmm
module load SQLite


blastx -query Trinity.fasta -db /share/databases/uniprot_sprot.fasta -num_threads 14 -max_target_seqs 1 -outfmt 6 > blastx.outfmt6 &
blastp -query transdecoder.pep -db /share/databases/uniprot_sprot.fasta -num_threads 14 -max_target_seqs 1 -outfmt 6 > blastp.outfmt6
wait && echo "Blast complete."

echo "Starting signalp"
signalp -f short -n signalp.out transdecoder.pep && echo "signalp complete."

echo"starting tmhmm"
tmhmm --short < transdecoder.pep > tmhmm.out && echo "tmhmm complete."
wait

echo "Starting rnamer."
perl /share/apps/trinotate-20131110/util/rnammer_support/RnammerTranscriptome.pl --transcriptome Trinity.fasta --path_to_rnammer /share/apps/rnammer_v1.2/rnammer
wait && echo "rnamer complete"

echo "Downloading latest databases..."
wget "http://sourceforge.net/projects/trinotate/files/TRINOTATE_RESOURCES/20140708/Trinotate.20140708.swissprot.sqlite.gz/download" -O Trinotate.sqlite.gz
wait && echo "Download complete."
gunzip Trinotate.sqlite.gz
wait && echo "Sucessfully decompressed."

perl /share/apps/trinityrnaseq-r20140413p1/util/support_scripts/get_Trinity_gene_to_trans_map.pl Trinity.fasta >  Trinity.fasta.gene_trans_map
echo "Gene map created."

echo "Loading gene map into sql..."
Trinotate Trinotate.sqlite init --gene_trans_map Trinity.fasta.gene_trans_map --transcript_fasta Trinity.fasta --transdecoder_pep transdecoder.pep
echo "Loading blast..."
# load protein hits
Trinotate Trinotate.sqlite LOAD_swissprot_blastp blastp.outfmt6
# load transcript hits
Trinotate Trinotate.sqlite LOAD_swissprot_blastx blastx.outfmt6
echo "Loading signalp, rnnmer, tmhmm data..."
Trinotate Trinotate.sqlite LOAD_pfam TrinotatePFAM.out
Trinotate Trinotate.sqlite LOAD_tmhmm tmhmm.out
Trinotate Trinotate.sqlite LOAD_signalp signalp.out
echo "Done."

echo "Generating Annotation Report"
Trinotate Trinotate.sqlite report -E 0.0001 > trinotate_annotation_report.xls
echo "DONE!" 