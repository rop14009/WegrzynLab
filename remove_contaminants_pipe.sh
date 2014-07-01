#This code is unfinished. 
module load python
module load Biopython

cut -f1 AnnotationReport-* | cut -f1 -d " " > clean_hits.names

python remove_contaminants.py

