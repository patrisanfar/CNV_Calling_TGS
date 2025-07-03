#!/bin/bash
# Author: Patricia Sánchez Fariña
# Project: TFM MADOBIS 2025
# Main launcher script for CNV detection and annotation pipeline

# Load configuration variables
source config.sh

echo "[!] Launching CNV detection methods..."

# Method 1: panelcn.mops
jid1=$(sbatch --parsable "$SCRIPT_DIR/panelcn.mops/00_launcher_panelcn.mops.sbatch" \
  "$DIR_BAM" "$PROJECT_NAME" "$PANELCNMOPS_bed" "$PANELCNMOPS_output")

# Method 2: CODEX2
jid2=$(sbatch --parsable "$SCRIPT_DIR/CODEX2/00_launcher_CODEX2.sbatch" \
  "$DIR_BAM" "$PROJECT_NAME" "$PANEL" "$CODEX2_output")

# Method 3: CLINCNV
jid3=$(sbatch --parsable "$SCRIPT_DIR/CLINCNV/00_launcher_CLINCNV.sbatch" \
  "$CLINCNV_bed" "$CLINCNV_normal" "$CLINCNV_normalOfftarget" \
  "$CLINCNV_bedOfftarget" "$CLINCNV_output")

# Method 4: CONVADING
jid4=$(sbatch --parsable "$SCRIPT_DIR/CONVADING/00_launcher_convading.sbatch" \
  "$DIR_BAM" "$PANEL" "$CONVADING_output" "$PROJECT_NAME" "$CONVADING_fai")

# Method 5: DECoN (consists of 3 sequential internal steps)
jid5=$(bash "$SCRIPT_DIR/DECON/00_launcher_DECON_pipeline.sh")

# Collect all job IDs to build SLURM dependency chain
all_jobs="${jid1}:${jid2}:${jid3}:${jid4}:${jid5}"

# Merge CNV results after all methods have finished
jid_integrate=$(sbatch --parsable --dependency=afterok:${all_jobs} "$SCRIPT_DIR/mixer/00_launcher_mixer.sh" \
"$PROJECT_NAME" "$PANEL" "$MIXER_output")

# Annotate merge CNV results
jid_annotation=$(sbatch --dependency=afterok:${jid_integrate} "$SCRIPT_DIR/annotation/00_launcher_annotsv.sh" \
"$ANNOTATION_dir" "$FASTA" "$ANNOTATION_file")

# Combine structural and functional data into a final annotated table
sbatch --dependency=afterok:${jid_annotation} "$SCRIPT_DIR/annotation/00_launcher_ResultAnnotationMix.sh" \
"$ANNOTATION_file" "$ANNOTATION_annotatedFile"