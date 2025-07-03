#!/bin/bash

#SBATCH -c 32
#SBATCH --mem=40G
#SBATCH -t 6-23:59:00
#SBATCH -o log_CLINCNV.out
#SBATCH -e log_CLINCNV.err
#SBATCH -J CLINCNV

apptainer exec PATH_TO_CONTAINER_DIR/ClinCNV_v1.18.3.sif Rscript clinCNV.R \
--bed "$1" \
--normal "$2" \
--normalOfftarget "$3" \
--bedOfftarget "$4" \
--out "$5"