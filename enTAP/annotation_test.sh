#!/bin/bash -l
# NOTE the -l flag!
######################
##SLURM CONFIGURATION#
######################

#SBATCH --job-name=annotation
#SBATCH --output=shell05-%j.output
#SBATCH --error=shell05-%j.error
#SBATCH --partition=bigmemh
#SBATCH --mail-type=ALL # other options are ALL, NONE, BEGIN, FAIL
#SBATCH --mail-user=samuel.ginzburg@uconn.edu

#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=30000
#SBATCH --time=2-12:00:00
#

####################
##END CONFIGURATION$
####################

##Log Start time
date
hostname

##Load your modules here
module load benchmarks

##Your bash commands here

python -u report_generator.py > test_output.txt

##log ending time
date
