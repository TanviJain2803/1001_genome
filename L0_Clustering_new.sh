#!/bin/bash
#$ -N tanvi_test
#$ -cwd
#$ -pe smp 12
#$ -l h_vmem=12G
#$ -V
#$ -e error_log_file
#$ -o output_log_file
#$ -q gale.q

cd /gale/ddn/ddn_neomorph/tjain

Rscript L0_Clustering.R


