#!/bin/bash
#$ -N tanvi_test
#$ -cwd
#$ -pe smp 4
#$ -l h_vmem=12G
#$ -V
#$ -e error_log_file_clus
#$ -o output_log_file_clus
#$ -q gale.q

cd /gale/ddn/ddn_neomorph/tjain/1001_genome

Rscript clustree_findMarkers.R "$1"
