#!/bin/bash

#SBATCH -c 32
#SBATCH --mem=40G
#SBATCH -t 6-23:59:00
#SBATCH -o log_convading.out
#SBATCH -e log_convading.err
#SBATCH -J convading

apptainer exec PATH_TO_CONTAINER_DIR/fjd_variantcalling.sif python CoNVading_pipeline.py "$1" "$2" "$3" "$4" "$5"
