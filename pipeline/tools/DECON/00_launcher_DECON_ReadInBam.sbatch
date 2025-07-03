#!/bin/bash

#SBATCH -c 32
#SBATCH --mem=40G
#SBATCH -t 6-23:59:00
#SBATCH -o log_DECON_ReadInBams.out
#SBATCH -e log_DECON_ReadInBams.err
#SBATCH -J DECON_ReadInBams

apptainer exec PATH_TO_CONTAINER_DIR/r-exomedepth_1.1.16--r44h9f7ea62_4.sif Rscript ReadInBams.R \
--bams "$1" \
--bed "$2" \
--fasta "$3" \
--out "$4"