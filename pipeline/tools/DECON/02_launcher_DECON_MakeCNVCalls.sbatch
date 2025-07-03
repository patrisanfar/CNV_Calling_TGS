#!/bin/bash

#SBATCH -c 32
#SBATCH --mem=40G
#SBATCH -t 6-23:59:00
#SBATCH -o log_DECON_makeCNVCalls.out
#SBATCH -e log_DECON_makeCNVCalls.err
#SBATCH -J DECON_makeCNVCalls

apptainer exec PATH_TO_CONTAINER_DIR/r-exomedepth_1.1.16--r44h9f7ea62_4.sif Rscript makeCNVcalls.R \
--RData "$1" \
--out "$2"

