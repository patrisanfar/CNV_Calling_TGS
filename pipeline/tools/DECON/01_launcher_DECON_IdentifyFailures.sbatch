#!/bin/bash

#SBATCH -c 32
#SBATCH --mem=40G
#SBATCH -t 6-23:59:00
#SBATCH -o log_DECON_IdentifyFailures.out
#SBATCH -e log_DECON_IdentifyFailures.err
#SBATCH -J DECON_IdentifyFailures

apptainer exec PATH_TO_CONTAINER_DIR/r-exomedepth_1.1.16--r44h9f7ea62_4.sif Rscript IdentifyFailures.R \
--RData "$1" \
--mincorr 0.98 \
--mincov 100 \
--out "$2"