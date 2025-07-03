#!/bin/bash

#SBATCH -c 32
#SBATCH --mem=40G
#SBATCH -t 6-23:59:00
#SBATCH -o log_finalMerge.out
#SBATCH -e log_finalMerge.err
#SBATCH -J finalMerge

apptainer exec PATH_TO_CONTAINERS_DIR/verse_4.4.2.sif RScript CNVfinalMerge_v2.R \
-d "$1" \
-a "$2"