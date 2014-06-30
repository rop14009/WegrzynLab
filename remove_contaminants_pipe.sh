#This code is unfinished. 
module load python
module load Biopython

cut -f1 AnnotationReport-* > clean_hits.names

python remove_contaminants.py

