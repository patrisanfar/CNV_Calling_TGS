# Pipeline for Copy Number Variant Detection in Targeted Sequencing
### Development of a pipeline for Copy Number Variant detection in targeted sequencing

This repository contains the CNV detection pipeline developed as part of a master’s thesis project focused on detecting copy number variants (CNVs) from targeted sequencing data. The workflow integrates five open-source tools (DECoN, CODEX2, CoNVaDING, panelcn.MOPS and ClinCNV) into a unified pipeline for CNV detection, merging and annotation.

Once implemented, the pipeline is applied to a Parkinson’s disease cohort to evaluate its performance and detect clinically relevant variants.

Each script includes information about its purpose and inputs required. For a better understanding of the tools employed and the structure of the pipeline, please consult the workflow section below.

## Author
Patricia Sánchez Fariña

## 1. First steps

Pre-processing and alignment of raw FASTQ files are performed using the [nf-core/sarek](https://github.com/nf-core/sarek), a Nextflow-based workflow developed by the nf-core community for germline variant detection in whole-genome, exome or panel sequencing data.

This step is executed by launching the `sarek_launch.sbatch` script, which takes FASTQ files as input and generates BAM files aligned to the `GRCh38` reference genome. It has an associated configuration file called `config_sarek.sh`, which defines input paths and parameters used.

Once the BAMs are generated, the main CNV analysis is launched using the script `run_CNV_pipeline.sh`. This central script manages all subsequent steps (CNV calling, merging and annotation) through SLURM job submission, setting the order of execution via job dependencies.

All paths and parameters used throughout the pipeline are defined in the file `config.sh`.

## 2. CNV calling

The pipeline includes five tools commonly used for CNV detection in targeted NGS data:

- [DECoN](https://github.com/RahmanTeam/DECoN)
- [CODEX2](https://github.com/yuchaojiang/CODEX2)
- [CoNVaDING](https://github.com/molgenis/CoNVaDING)
- [panelcn.MOPS](https://github.com/bioinf-jku/panelcn.mops)
- [ClinCNV](https://github.com/imgag/ClinCNV)

Each tool is run independently using tool-specific launcher scripts (e.g., `00_launcher_panelcn.mops.sbatch`). All of them use input variables loaded from the shared `config.sh` file. Each launcher internally runs the original scripts provided by the official repositories of the tools, which are linked above.

The outputs include per-sample CNV calls, specifying genomic coordinates, CNV type (deletion/duplication) and quality metrics. These results are later combined in the mixing step to identify CNVs detected by multiple tools.

## 3. CNV merging

After CNV calling, the pipeline merges the results from all tools using the `00_launcher_mixer.sh` script.

This step identifies overlapping CNVs based on exon midpoints, grouping together events from different tools that affect the same gene, sample and CNV type. 

The output is a unified CNV table summarizing all detected events across methods.

This step is automatically launched within the main pipeline script `run_CNV_pipeline.sh` after the CNV calling stage.

## 4. Annotation 

Once the CNV merging step is completed, CNVs are annotated in two steps:

- `00_launcher_annotsv.sh` runs AnnotSV to add functional and clinical information based on the `GRCh38` reference genome.

- `00_launcher_ResultAnnotationMix.sh` combines the annotation output with the original CNV table of the mixing step to produce a final annotated file.

Both annotation steps are automatically launched by `run_CNV_pipeline.sh`, producing a final CNV file ready for downstream analysis and potential clinical use.

## 5. Workflow

The structure of the complete pipeline is shown in the diagram below:

![Workflow](https://github.com/user-attachments/assets/e39902c7-7aa5-4b33-be83-10f47cd7caf9)

