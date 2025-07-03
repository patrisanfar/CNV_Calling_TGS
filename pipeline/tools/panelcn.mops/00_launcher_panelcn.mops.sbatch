#!/bin/bash

#SBATCH -c 32
#SBATCH --mem=40G
#SBATCH -t 6-23:59:00
#SBATCH -o log_panelcn.mops.out
#SBATCH -e log_panelcn.mops.err
#SBATCH -J panelcn.mop

apptainer exec PATH_TO_CONTAINERS_DIR/bioconductor-panelcn.mops_1.28.0--r44hdfd78af_0.sif Rscript panelcnMOPS.R \
--dir "$1" \
--name "$2" \
--bed "$3" \
--outputdir "$4"