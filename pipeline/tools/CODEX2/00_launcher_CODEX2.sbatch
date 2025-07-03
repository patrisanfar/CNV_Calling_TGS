#!/bin/bash

#SBATCH -c 32
#SBATCH --mem=40G
#SBATCH -t 6-23:59:00
#SBATCH -o log_CODEX2.out
#SBATCH -e log_CODEX2.err
#SBATCH -J CODEX2

apptainer exec PATH_TO_CONTAINERS_DIR/fjd_variantcalling.sif Rscript CODEX2.R \
-d "$1" \
-n "$2" \
-b "$3" \
-o "$4" 