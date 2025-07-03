#!/bin/bash
## config.sh
# Configuration file for CNV detection pipeline

## Input parameters
# TO EDIT by user
# BAM directory
DIR_BAM="/path/to/BAMs"
# Project name
PROJECT_NAME="PatriTFM"
# Panel
PANEL="path/to/panel.bed"
# Reference genome (fasta)
FASTA="path/to/genome.fasta"

## Panelcn.mops
# Panel
# Example of format-> 1	11023093	11023393	TARDBP.E6.chr1.11023093.11023393
PANELCNMOPS_bed="bedPath"
# Output directory for panelcn.mops results
PANELCNMOPS_output="outputPath"

# CODEX2
# Output directory for CODEX2 results
CODEX2_output="outputPath"

# CLINCNV
# BED file with GC content (from preprocessing step)
CLINCNV_bed="path/to/clinCNV_gc.bed"
# Coverage file for on-target regions
CLINCNV_normal="path/to/coverage_ontarget.txt"
# Coverage file for off-target regions
CLINCNV_normalOfftarget="path/to/coverage_offtarget.txt"
# BED file for off-target regions
CLINCNV_bedOfftarget="path/to/bed_offtarget.bed"
# Output directory for CLINCNV results
CLINCNV_output="outputPath"

# CONVADING
# Index file for the reference genome
CONVADING_fai="path/to/genome.fasta.fai"
# Output directory for CONVADING results
CONVADING_output="outputPath"

# DECON
# Output directory for DECON results
DECON_output="outputPath"

# MIXER
# Output directory for MIXER results
MIXER_output="outputPath"

# ANNOTATION
# Directory for annotation step
ANNOTATION_dir="path"
# CNV table to annotate
ANNOTATION_file=".combined.txt"
# Annotated CNV table (final output)
ANNOTATION_annotatedFile=".combinedAnnotated.txt"