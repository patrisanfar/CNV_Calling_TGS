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

All paths and parameters used throughout the pipeline are defined in the file `config.sh`.

## 2. CNV calling

Copy number variation detection was performed using five open-source tools suitable for use with targeted sequencing data:

- [DECoN](https://github.com/RahmanTeam/DECoN)
- [CODEX2](https://github.com/yuchaojiang/CODEX2)
- [CoNVaDING](https://github.com/molgenis/CoNVaDING)
- [panelcn.MOPS](https://github.com/bioinf-jku/panelcn.mops)
- [ClinCNV](https://github.com/imgag/ClinCNV)

Each tool was run independently using tool-specific launcher scripts (e.g., `00_launcher_panelcn.mops.sbatch`). All receive their parameters from the central configuration file.

These tools were executed automatically within the `run_CNV_pipeline.sh` script. Each tool was run independently using tool-specific launcher scripts (e.g., `00_launcher_panelcn.mops.sbatch`) , with input parameters loaded from the shared `config.sh` configuration file.

The output of each tool includes per-sample CNV calls, specifying genomic coordinates, variant type (deletion/duplication), and quality metrics. These results are later combined in the mixing step to identify consensus CNVs.

## 3. Mixing

After CNV detection has finished, the pipeline proceeds to merge the results from the five tools. This step is handled by the script `00_launcher_mixer.sh`.

The merging strategy relies on exon midpoint overlap to group CNVs affecting the same gene, sample, and variant type. For each unified event, the script records how many tools detected it and includes relevant quality metrics. The resulting CNV table includes the variant type, coordinates, the number of tools that detected each event and tool-specific quality metrics.

This step is automatically launched within the main pipeline script `run_CNV_pipeline.sh` after the CNV calling stage.

## 4. Annotation 

Once the CNV merging step is completed, the pipeline performs functional and clinical annotation of the resulting variants. This is done in two sequential steps:

- The script `00_launcher_annotsv.sh` runs AnnotSV to annotate the unified CNV list using the `GRCh38` reference genome.

- Then, `00_launcher_ResultAnnotationMix.sh` combines the annotation output with the original CNV table, producing a final annotated file that unifies both structural and functional information.

These steps are also submitted by `run_CNV_pipeline.sh`. The resulting annotated CNVs are suitable for downstream analysis and clinical interpretation.

## 5. Workflow

The complete pipeline structure is summarized in the diagram below:

![Workflow](https://github.com/user-attachments/assets/e39902c7-7aa5-4b33-be83-10f47cd7caf9)

