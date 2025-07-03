#!/bin/bash

#SBATCH -c 32
#SBATCH --mem=40G
#SBATCH -t 99:00:00
#SBATCH -o log_TGS_getBams.out
#SBATCH -e log_TGS_getBams.err
#SBATCH -J Sarek
 
# TMP directories
export NXF_APPTAINER_TMPDIR=$PWD/tmp/$SLURM_JOB_ID
export TMPDIR=$PWD/tmp/$SLURM_JOB_ID

# Nextflow settings
SCRATCH=false
STAGEINMODE="symlink"

# Environment
shopt -s expand_aliases
source ~/.bash_aliases
conda_activate
conda activate env_nf

# Load config variables
source config_sarek.sh

nextflow run nf-core/sarek -r 3.5.1 \
 -profile apptainer \
 -process.scratch=${SCRATCH} \
 -process.stageInMode=${STAGEINMODE} \
 -process.maxForks=200 \
 -process.errorStrategy='ignore' \
 -process.time='2h' \
 -process.cache='lenient' \
 -c /path/to/sarek.config \
 --outdir "$SAREK_OUTPUT_DIR" \
 --input "$SAMPLES_CSV" \
 --wes \
 --intervals "$PANEL" \
 --igenomes_base "$IGENOMES_PATH" \
 --trim_fastq \
 --aligner bwa-mem2 \
 --save_mapped  \
 --save_output_as_bam  \
 --genome GATK.GRCh38