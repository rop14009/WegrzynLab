#!/bin/bash -l
#SBATCH -J rmcontam
#SBATCH -o rmcontam-%j.output
#SBATCH -e rmcontam-%j.err
#SBATCH -n 4
#SBATCH -N 1
#SBATCH -p bigmem

module load python
module load Biopython

cut -f1 contaminated-*.txt | cut -f1 -d " " > contaminants.names

python remove_contaminants.py

