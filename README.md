# Pipeline for Copy Number Variant Detection in Targeted Sequencing
### Development of a pipeline for Copy Number Variant detection in targeted sequencing

This repository contains the pipeline developed as part of a master’s thesis focused on detecting copy number variants (CNVs) from targeted NGS data. The workflow integrates five open-source tools (DECoN, CODEX2, CoNVaDING, panelcn.MOPS and ClinCNV) into a unified framework for CNV calling and annotation.

Once implemented, the pipeline was applied to a Parkinson’s disease cohort to assess its performance and explore its potential for genetic characterization.

For a better understanding of the tools employed and the structure of the pipeline, please consult the workflow section below.

Each script includes information about its purpose, input requirements and filtering strategies.

## Author
Patricia Sánchez Fariña

## 1. First steps

Pre-processing and alignment were performed using the [nf-core/sarek](https://github.com/nf-core/sarek) framework, a workflow written in Nextflow by the nf-core community for germline variant detection from whole genome, exome, or targeted panel sequencing data.

This step was executed by launching the `sarek_launch.sbatch` script, which takes raw FASTQ files as input and produces preprocessed BAM files aligned to the `GRCh38` reference genome.

Once preprocessed BAMs are available, the full CNV detection and annotation workflow is launched using the script `run_CNV_pipeline.sh`. This main script manages all subsequent steps (CNV calling, merging and annotation) through SLURM job submission, setting the order of execution via job dependencies.

## 2. CNV calling

Copy number variation detection was performed using five open-source tools suitable for use with targeted sequencing data:

- [DECoN](https://github.com/RahmanTeam/DECoN)
- [CODEX2](https://github.com/yuchaojiang/CODEX2)
- [CoNVaDING](https://github.com/molgenis/CoNVaDING)
- [panelcn.MOPS](https://github.com/bioinf-jku/panelcn.mops)
- [ClinCNV](https://github.com/imgag/ClinCNV)

Each method was executed separately using tool-specific launcher scripts (e.g., 00_launcher_panelcn.mops.sbatch) with a shared BED file defining the panel regions. All jobs were submitted by run_CNV_pipeline.sh and executed in a high-performance computing environment using SLURM. Containerized environments were handled via Apptainer to ensure full reproducibility.

The outputs of each tool consist of CNV calls per sample, including genomic coordinates, CNV type (deletion/duplication), and quality metrics.

## 3. Mixing

