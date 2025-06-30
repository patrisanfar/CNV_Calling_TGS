# Pipeline for Copy Number Variant Detection in Targeted Sequencing
### Development of a pipeline for Copy Number Variant detection in targeted sequencing

This repository contains the pipeline developed as part of a master’s thesis focused on detecting copy number variants (CNVs) from targeted NGS data. The workflow integrates five open-source tools (DECoN, CODEX2, CoNVaDING, panelcn.MOPS and ClinCNV) into a unified framework for CNV calling and annotation.

Once implemented, the pipeline was applied to a Parkinson’s disease cohort to assess its performance and explore its potential for genetic characterization.

For a better understanding of the tools employed and the structure of the pipeline, please consult the workflow section below.

Each script includes information about its purpose, input requirements and filtering strategies.

## Author
Patricia Sánchez Fariña

## 1. First steps

Pre-processing and alignment are performed using the nf-core/sarek framework, a workflow developed by the nf-core community in Nextflow for germline variant detection from whole genome, exome, or targeted sequencing data.

This step is executed by launching the `sarek_launch.sbatch` script, using FASTQ files as input.
