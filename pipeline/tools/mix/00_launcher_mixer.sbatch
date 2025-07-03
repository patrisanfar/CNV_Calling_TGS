#!/bin/bash

#SBATCH -c 32
#SBATCH --mem=40G
#SBATCH -t 6-23:59:00
#SBATCH -o log_Mixer.out
#SBATCH -e log_Mixer.err
#SBATCH -J Mixer

apptainer exec PATH_TO_CONTAINERS_DIR/fjd_variantcalling.sif Rscript CNV_result_mixer.R \
--name "$1" \
--bed "$2" \
--output "$3"