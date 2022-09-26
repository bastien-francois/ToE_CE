#!/bin/bash
# Specify batch queue
#PBS -q	h12 
# Specify execution shell
#PBS -S /bin/bash
# Request that stdin and stdout are merged in the same output file
#PBS -j eo
#PBS -l mem=50gb
#PBS -l vmem=55gb

module load R/4.0.5 hdf5/1.10.7

cd /home/bfrancois/Compound_Event/TCE6/Code

R CMD BATCH 1_2_Download_GCMs_0.5deg_1850_2100.R out_1_2_Download.out
